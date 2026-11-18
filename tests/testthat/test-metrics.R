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
  expect_equal(
    metrics$by_time$calib_slope[2],
    surv_calibration_slope(data$time, data$status, survival[, 2], 12)$slope
  )
  expect_equal(
    metrics$by_time$calib_intercept[2],
    surv_calibration_slope(data$time, data$status, survival[, 2], 12)$intercept
  )
  expect_identical(
    names(metrics$row),
    c("n", "n_events", "cindex", "ibs", "uno_c", "iauc", "iauc_uno",
      "brier_6", "auc_6", "auc_uno_6", "citl_6", "e50_6", "e90_6",
      "calib_intercept_6", "calib_slope_6",
      "brier_12", "auc_12", "auc_uno_12", "citl_12", "e50_12", "e90_12",
      "calib_intercept_12", "calib_slope_12",
      "brier_24", "auc_24", "auc_uno_24", "citl_24", "e50_24", "e90_24",
      "calib_intercept_24", "calib_slope_24")
  )
  # Without training data, Uno's estimators stay missing by design.
  expect_true(is.na(metrics$overall$uno_c))
  expect_true(all(is.na(metrics$by_time$auc_uno)))
  expect_false(is.na(metrics$overall$iauc))

  # A learner with no survival distribution still gets a C-index.
  risk_only <- calc_surv_metrics(
    data$time, data$status, fixture$eta, NULL, times
  )
  expect_equal(risk_only$overall$cindex, metrics$overall$cindex)
  expect_true(is.na(risk_only$overall$ibs))
  expect_true(all(is.na(risk_only$by_time$brier)))
})

test_that("Uno's C and AUC match survAUC exactly", {
  skip_if_not_installed("survAUC")
  library(survAUC)
  set.seed(11)
  data <- sim_surv_data(400, seed = 2)
  index <- sample.int(nrow(data), 200)
  truth <- sim_surv_truth(data, c(3, 6, 12))
  train <- data[-index, ]
  test <- data[index, ]
  tau <- 12

  ref_c <- survAUC::UnoC(
    survival::Surv(train$time, train$status),
    survival::Surv(test$time, test$status),
    truth$risk[index],
    time = tau
  )
  expect_equal(
    surv_cindex_uno(train$time, train$status, test$time, test$status,
                    truth$risk[index], tau = tau),
    ref_c,
    tolerance = 1e-10
  )
  expect_equal(
    surv_cindex_uno(train$time, train$status, test$time, test$status,
                    truth$risk[index], tau = tau),
    ref_c,
    tolerance = 1e-10
  )
  # With no `tau`, both implementations truncate at the largest observed test
  # time.
  ref_c_default <- survAUC::UnoC(
    survival::Surv(train$time, train$status),
    survival::Surv(test$time, test$status),
    truth$risk[index]
  )
  expect_equal(
    surv_cindex_uno(train$time, train$status, test$time, test$status,
                    truth$risk[index]),
    ref_c_default,
    tolerance = 1e-10
  )

  for (t0 in c(3, 6, 12)) {
    ref_a <- survAUC::AUC.uno(
      survival::Surv(train$time, train$status),
      survival::Surv(test$time, test$status),
      truth$risk[index],
      times = t0
    )$auc
    expect_equal(
      surv_auc_uno(train$time, train$status, test$time, test$status,
                   truth$survival[index, which(c(3, 6, 12) == t0)], t0),
      ref_a,
      tolerance = 1e-10
    )
  }

  # Risk-only evaluation still yields Uno's C.
  risk_only <- calc_surv_metrics(
    test$time, test$status, truth$risk[index], NULL, c(3, 6),
    time_train = train$time, status_train = train$status, tau = tau
  )
  expect_equal(risk_only$overall$uno_c, ref_c, tolerance = 1e-10)
  expect_true(is.na(risk_only$overall$iauc))
})

test_that("calibration-in-the-large, E50 and E90 are reproducible", {
  data <- sim_surv_data(300, seed = 3)
  truth <- sim_surv_truth(data, c(6, 12))
  calibration <- surv_calibration(data$time, data$status,
                                  truth$survival[, 2], time_point = 12)
  # Manual recomputation from the group table.
  expect_equal(
    calibration$citl,
    sum(calibration$table$n * calibration$table$difference) /
      sum(calibration$table$n)
  )
  # ICI stays the patient-weighted absolute difference.
  expect_equal(
    calibration$ici,
    sum(calibration$table$n * abs(calibration$table$difference)) /
      sum(calibration$table$n)
  )
  expect_true(calibration$e50 <= calibration$e90)
  expect_true(calibration$e90 >= calibration$ici)

  # The curves summary carries the same columns.
  curves <- surv_calibration_curves(data$time, data$status, truth$survival,
                                    c(6, 12))
  expect_true(all(c("ici", "citl", "e50", "e90") %in% names(curves$summary)))
})

test_that("risk groups and the stratified KM plot", {
  data <- sim_surv_data(300, seed = 4)
  truth <- sim_surv_truth(data, 12)
  groups <- surv_risk_groups(truth$risk, k = 3)
  expect_s3_class(groups, "factor")
  expect_equal(length(groups), nrow(data))
  expect_equal(levels(groups), c("Q1", "Q2", "Q3"))
  # Equal-sized groups, ordered from low to high risk.
  expect_equal(as.integer(table(groups)), rep(100L, 3L))
  expect_lt(mean(truth$risk[groups == "Q1"]), mean(truth$risk[groups == "Q3"]))
  expect_error(surv_risk_groups(c(1, 2, 2, 2), k = 3), "vary enough")
  expect_error(surv_risk_groups(c(NA, 1, 2), k = 2), "missing")

  spec <- surv_spec(time_points = 12)
  frame <- data
  frame$grp <- groups
  plot <- plot_surv_km(frame, spec, group = "grp", logrank = TRUE)
  expect_s3_class(plot, "ggplot")
  caption <- plot$labels$caption
  expect_true(grepl("Log-rank", caption))
  expect_true(grepl("p = ", caption))
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
