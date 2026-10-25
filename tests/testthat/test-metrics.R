# `riskRegression::Score()` evaluates its formula in the calling frame, so the
# survival namespace has to be attached for the comparison to be meaningful.
suppressPackageStartupMessages(library(survival))

cox_fixture <- function(n = 200L, seed = 11L) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  eta <- 0.8 * x1 - 0.4 * x2
  event <- stats::rexp(n, rate = exp(eta) / 12)
  censoring <- stats::rexp(n, rate = 1 / 18)
  time <- pmin(event, censoring)
  status <- as.integer(event <= censoring)
  data <- data.frame(time = time, status = status, x1 = x1, x2 = x2)
  fit <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2,
    data = data,
    x = TRUE,
    y = TRUE
  )
  list(data = data, fit = fit, eta = eta)
}

# Predicted survival from a coxph fit at the requested time points, one row per
# patient, matching what `survfit()` and therefore `riskRegression` use.
cox_survival <- function(fit, data, times) {
  summary_fit <- summary(
    survival::survfit(fit, newdata = data),
    times = times,
    extend = TRUE
  )
  t(summary_fit$surv)
}

test_that("surv_cindex() equals survival::concordance()", {
  fixture <- cox_fixture()
  data <- fixture$data
  mine <- surv_cindex(data$time, data$status, fixture$eta)
  reference <- survival::concordance(
    survival::Surv(data$time, data$status) ~ fixture$eta,
    reverse = TRUE
  )$concordance
  expect_equal(mine, as.numeric(reference), tolerance = 1e-10)

  # A risk score that ranks perfectly gives 1; the reversed score gives 0.
  perfect <- -data$time
  expect_equal(surv_cindex(data$time, data$status, perfect), 1, tolerance = 1e-8)
  expect_equal(surv_cindex(data$time, data$status, -perfect), 0, tolerance = 1e-8)

  # No comparable pair is reported as missing, not as a number.
  expect_true(is.na(surv_cindex(c(1, 2), c(0L, 0L), c(0.1, 0.2))))
})

test_that("the IPCW Brier score matches riskRegression::Score()", {
  skip_if_not_installed("riskRegression")
  skip_if_not_installed("prodlim")
  fixture <- cox_fixture()
  data <- fixture$data
  times <- c(6, 12, 24)
  survival <- cox_survival(fixture$fit, data, times)

  mine <- vapply(seq_along(times), function(index) {
    surv_brier(data$time, data$status, survival[, index], times[[index]])
  }, numeric(1))

  score <- riskRegression::Score(
    object = list(cox = fixture$fit),
    formula = Surv(time, status) ~ 1,
    data = data,
    times = times,
    metrics = c("brier", "auc"),
    cens.model = "marginal",
    conf.int = FALSE,
    verbose = FALSE
  )
  reference <- subset(score$Brier$score, model == "cox")$Brier
  expect_equal(mine, as.numeric(reference), tolerance = 1e-10)
})

test_that("the time-dependent AUC matches riskRegression::Score()", {
  skip_if_not_installed("riskRegression")
  skip_if_not_installed("prodlim")
  fixture <- cox_fixture()
  data <- fixture$data
  times <- c(6, 12, 24)
  survival <- cox_survival(fixture$fit, data, times)

  mine <- vapply(seq_along(times), function(index) {
    surv_auc(data$time, data$status, survival[, index], times[[index]])
  }, numeric(1))

  score <- riskRegression::Score(
    object = list(cox = fixture$fit),
    formula = Surv(time, status) ~ 1,
    data = data,
    times = times,
    metrics = c("brier", "auc"),
    cens.model = "marginal",
    conf.int = FALSE,
    verbose = FALSE
  )
  reference <- subset(score$AUC$score, model == "cox")$AUC
  expect_equal(mine, as.numeric(reference), tolerance = 1e-10)
})

test_that("the AUC is the area under the ROC curve the package draws", {
  fixture <- cox_fixture()
  data <- fixture$data
  times <- c(6, 12, 24)
  survival <- cox_survival(fixture$fit, data, times)

  for (index in seq_along(times)) {
    curve <- surv_roc_curve(data$time, data$status, survival[, index], times[[index]])
    area <- sum(
      diff(curve$fpr) *
        (utils::head(curve$tpr, -1L) + utils::tail(curve$tpr, -1L)) / 2
    )
    expect_equal(area, surv_auc(data$time, data$status, survival[, index], times[[index]]),
                 tolerance = 1e-10)
    expect_true(all(curve$fpr >= 0 & curve$fpr <= 1))
    expect_true(all(curve$tpr >= 0 & curve$tpr <= 1))
    expect_false(is.unsorted(curve$fpr))
  }
})

test_that("without censoring the IPCW weights are one", {
  set.seed(4)
  n <- 120L
  time <- stats::rexp(n, rate = 1 / 10)
  status <- rep(1L, n)
  survival <- stats::runif(n, 0.1, 0.9)
  time_point <- 8

  expected <- mean((as.numeric(time <= time_point) - (1 - survival))^2)
  expect_equal(
    surv_brier(time, status, survival, time_point, censoring = "none"),
    expected,
    tolerance = 1e-12
  )
  expect_equal(
    surv_brier(time, status, survival, time_point, censoring = "marginal"),
    expected,
    tolerance = 1e-10
  )
})

test_that("a perfect prediction has a Brier score of zero", {
  set.seed(5)
  n <- 100L
  time <- stats::rexp(n, rate = 1 / 10)
  status <- rep(1L, n)
  time_point <- 12
  # Survival is either 1 (event after t) or 0 (event before t): the true value.
  survival <- as.numeric(time > time_point)
  expect_equal(
    surv_brier(time, status, survival, time_point, censoring = "none"),
    0,
    tolerance = 1e-12
  )
})

test_that("surv_ibs() integrates the Brier score over the time points", {
  fixture <- cox_fixture()
  data <- fixture$data
  times <- c(6, 12, 18, 24)
  survival <- cox_survival(fixture$fit, data, times)

  scores <- vapply(seq_along(times), function(index) {
    surv_brier(data$time, data$status, survival[, index], times[[index]])
  }, numeric(1))
  expected <- sum(
    diff(times) * (utils::head(scores, -1L) + utils::tail(scores, -1L)) / 2
  ) / diff(range(times))

  expect_equal(
    surv_ibs(data$time, data$status, survival, times),
    expected,
    tolerance = 1e-12
  )
  # One time point is just the Brier score.
  expect_equal(
    surv_ibs(data$time, data$status, survival[, 1, drop = FALSE], times[1]),
    scores[[1L]],
    tolerance = 1e-12
  )
})

test_that("calc_surv_metrics() reports the same numbers as its parts", {
  fixture <- cox_fixture()
  data <- fixture$data
  times <- c(6, 12, 24)
  survival <- cox_survival(fixture$fit, data, times)

  metrics <- calc_surv_metrics(
    data$time, data$status, fixture$eta, survival, times
  )
  expect_equal(metrics$overall$cindex, surv_cindex(data$time, data$status, fixture$eta))
  expect_equal(metrics$overall$n, nrow(data))
  expect_equal(metrics$overall$n_events, sum(data$status))
  expect_equal(metrics$by_time$brier[2], surv_brier(data$time, data$status, survival[, 2], 12))
  expect_equal(metrics$by_time$auc[2], surv_auc(data$time, data$status, survival[, 2], 12))
  expect_identical(
    names(metrics$row),
    c("n", "n_events", "cindex", "ibs", "brier_6", "auc_6",
      "brier_12", "auc_12", "brier_24", "auc_24")
  )

  # A learner with no survival distribution still gets a C-index.
  risk_only <- calc_surv_metrics(
    data$time, data$status, fixture$eta, NULL, times
  )
  expect_equal(risk_only$overall$cindex, metrics$overall$cindex)
  expect_true(is.na(risk_only$overall$ibs))
  expect_true(all(is.na(risk_only$by_time$brier)))
})

test_that("metric inputs are validated", {
  expect_error(surv_cindex(c(1, 2), c(0, 1), c(0.1)), "same length")
  expect_error(surv_cindex(c(1, 2), c(0, 2), c(0.1, 0.2)), "0/1")
  expect_error(surv_brier(c(1, 2), c(0, 1), c(0.5, 1.5), 1), "between")
  expect_error(
    surv_auc(c(1, 2, 3), c(0, 1, 0), c(0.1, 0.2, 0.3), -1),
    "time_point"
  )
  expect_error(surv_ibs(c(1, 2), c(0, 1), matrix(0.5, nrow = 2L), c(1, 2)),
               "one column per time point")
})
