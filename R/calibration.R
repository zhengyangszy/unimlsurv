# L6 calibration layer -------------------------------------------------------

#' Time-point calibration
#'
#' Compares predicted survival probabilities with observed survival at one time
#' point, in groups of patients with similar predictions. The observed value is
#' the Kaplan-Meier estimate inside the group, so censoring is respected rather
#' than ignored.
#'
#' The summary is the integrated calibration index (ICI): the number of patients
#' in a group times the absolute difference between observed and predicted
#' survival, summed over groups and divided by the total. A value of 0 means the
#' predictions sit on the Kaplan-Meier curve.
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param survival Predicted survival probability at `time_point`.
#' @param time_point The time point being calibrated.
#' @param groups Number of prediction groups.
#' @param censoring Censoring model, kept for symmetry with the other metrics.
#'
#' @return A list with `table` (one row per group), `ici`, `n`, `events`, and a
#'   `message` when the calibration is not estimable.
#' @references
#' Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW (2019).
#' Calibration: the Achilles heel of predictive analytics. *BMC Medicine*, 17,
#' 230. doi:10.1186/s12916-019-1466-7
#'
#' Austin PC, Steyerberg EW (2019). The integrated calibration index (ICI) and
#' related metrics for assessing the calibration of survival models with binary
#' outcomes. *Statistics in Medicine*, 38(28), 5476-5494. doi:10.1002/sim.8356
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, 36)
#' calibration <- surv_calibration(data$time, data$status, truth$survival[, 1],
#'                                 time_point = 36)
#' calibration$ici
surv_calibration <- function(time, status, survival, time_point, groups = 10L,
                             censoring = "marginal") {
  .check_survival_inputs(time, status, survival = survival)
  .assert_scalar_number(time_point, "time_point", min = 0)
  .assert_scalar_number(groups, "groups", integer = TRUE, min = 2)

  n <- length(time)
  empty <- data.frame(
    group = integer(), n = integer(), predicted = numeric(),
    observed = numeric(), difference = numeric(),
    stringsAsFactors = FALSE
  )
  if (!any(time <= time_point) || !any(time > time_point)) {
    return(list(
      table = empty,
      ici = NA_real_,
      n = n,
      events = sum(status == 1L),
      message = "The time point does not split the cohort into cases and controls."
    ))
  }

  breaks <- unique(stats::quantile(
    survival,
    probs = seq(0, 1, length.out = groups + 1L),
    na.rm = TRUE,
    names = FALSE
  ))
  if (length(breaks) < 3L) {
    return(list(
      table = empty,
      ici = NA_real_,
      n = n,
      events = sum(status == 1L),
      message = "The predicted survival does not vary enough to form groups."
    ))
  }
  bin <- cut(survival, breaks = breaks, include.lowest = TRUE)
  rows <- lapply(split(seq_len(n), bin), function(index) {
    if (length(index) == 0L) {
      return(NULL)
    }
    km <- .km_survival(time[index], status[index])(time_point)
    data.frame(
      group = as.integer(bin[index][[1L]]),
      n = length(index),
      predicted = mean(survival[index]),
      observed = as.numeric(km),
      difference = as.numeric(km) - mean(survival[index]),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(table) || nrow(table) == 0L) {
    return(list(
      table = empty,
      ici = NA_real_,
      n = n,
      events = sum(status == 1L),
      message = "No calibration group could be formed."
    ))
  }
  rownames(table) <- NULL
  ici <- sum(table$n * abs(table$difference)) / sum(table$n)
  list(
    table = table,
    ici = ici,
    n = n,
    events = sum(status == 1L),
    message = NULL
  )
}

# The complementary log-log link turns a survival probability into the scale on
# which a proportional hazards model is additive: log(-log(S(t|x))) is the linear
# predictor plus the baseline term. A predicted survival of exactly 1 would give
# -Inf, so it is truncated the same way `.safe_logit()` truncates probabilities in
# `uniml`.
.safe_cloglog <- function(survival) {
  eps <- 1e-12
  safe <- pmin(pmax(survival, eps), 1 - eps)
  pmin(pmax(log(-log(safe)), -50), 50)
}

#' Time-point recalibration coefficients
#'
#' Fits the survival analogue of a logistic recalibration model at one time point
#' and returns its slope and intercept. The outcome "had the event by
#' `time_point`" is regressed on the complementary log-log of the predicted
#' survival probability, with a complementary log-log link:
#'
#' `cloglog(P(T <= t | x)) = intercept + slope * cloglog(1 - S(t | x))`
#'
#' `1 - S(t | x)` is the predicted event probability, so a perfectly calibrated
#' prediction satisfies `intercept = 0` and `slope = 1`, exactly as in the binary
#' case. A slope below 1 means the predictions are too extreme (over-dispersed), a
#' slope above 1 means they are too timid, and a non-zero intercept means they are
#' systematically too high or too low.
#'
#' Patients censored before `time_point` carry no information about the outcome
#' and are dropped; the patients who remain are weighted by the inverse of the
#' Kaplan-Meier censoring survival, which is the same inverse-probability-of-
#' censoring weighting that `surv_brier()` uses. Without it, censoring that
#' depends on the covariates would bias the coefficients.
#'
#' With fewer than 20 patients or fewer than 5 events or 5 non-events at the time
#' point the coefficients are returned as `NA` with a `message`, because a
#' two-parameter fit on a handful of patients produces a number nobody should
#' report. The same applies when the fit separates the outcome or fails to
#' converge. When the predicted survival is constant the slope is not identified,
#' so only the intercept is returned, which in that case is calibration-in-the-
#' large.
#'
#' Passing that floor is not the same as having enough data. A slope from a few
#' dozen patients can be far from 1 while the model is fine, because the spread of
#' the predicted survival on the cloglog scale is small and the fit mostly
#' amplifies noise. Treat a slope as worth reading at roughly 100 patients with a
#' few dozen events; `n` and `events` are returned precisely so a reader can
#' judge, and a report that gives the slope without them is not reproducible.
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param survival Predicted survival probability at `time_point`.
#' @param time_point The time point being calibrated.
#' @param censoring Censoring model: `"marginal"` or `"none"`.
#'
#' @return A list with `slope`, `intercept`, `n`, `events`, and a `message` when
#'   the coefficients are not estimable.
#' @references
#' van Houwelingen JC, le Cessie S (1990). Predictive value of statistical models.
#' *Statistics in Medicine*, 9(11), 1303-1325. doi:10.1002/sim.4780091109
#'
#' Royston P, Altman DG (2013). External validation of a Cox prognostic model:
#' principles and methods. *BMC Medical Research Methodology*, 13, 33.
#' doi:10.1186/1471-2288-13-33
#'
#' Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW (2019).
#' Calibration: the Achilles heel of predictive analytics. *BMC Medicine*, 17,
#' 230. doi:10.1186/s12916-019-1466-7
#' @export
#' @examples
#' data <- sim_surv_data(n = 300, seed = 1)
#' truth <- sim_surv_truth(data, 36)
#' # The simulated truth is calibrated by construction: slope near 1, intercept
#' # near 0.
#' surv_calibration_slope(data$time, data$status, truth$survival[, 1],
#'                        time_point = 36)
surv_calibration_slope <- function(time, status, survival, time_point,
                                   censoring = "marginal") {
  .check_survival_inputs(time, status, survival = survival)
  .assert_scalar_number(time_point, "time_point", min = 0)
  event_at <- status == 1L & time <= time_point
  survived <- time > time_point
  out <- list(
    slope = NA_real_,
    intercept = NA_real_,
    n = length(time),
    events = sum(event_at),
    message = NULL
  )
  if (!any(event_at) || !any(survived)) {
    out$message <- paste(
      "The time point does not split the cohort into cases and controls."
    )
    return(out)
  }
  if (sum(event_at) < 5L || sum(survived) < 5L || length(time) < 20L) {
    out$message <- paste(
      "Too few cases, controls, or patients at this time point for a",
      "two-parameter recalibration fit."
    )
    return(out)
  }
  G <- .censoring_survival(time, status, censoring)
  weight <- numeric(length(time))
  weight[event_at] <- 1 / G(time[event_at])
  weight[survived] <- 1 / G(time_point)
  keep <- event_at | survived
  outcome <- as.integer(event_at[keep])
  linear <- .safe_cloglog(survival[keep])
  weights <- weight[keep]
  if (length(unique(linear)) < 2L) {
    # A constant prediction has no slope to estimate, but the offset-only fit
    # still says whether the level is right; that is calibration-in-the-large.
    shifted <- tryCatch(
      suppressWarnings(stats::glm(
        outcome ~ 1,
        offset = linear,
        family = stats::binomial(link = "cloglog"),
        weights = weights
      )),
      error = function(e) NULL
    )
    if (!is.null(shifted) && isTRUE(shifted$converged)) {
      out$intercept <- unname(stats::coef(shifted)[[1L]])
    }
    out$message <- paste(
      "The predicted survival does not vary at this time point, so only the",
      "intercept is estimable."
    )
    return(out)
  }
  fit <- tryCatch(
    suppressWarnings(stats::glm(
      outcome ~ linear,
      family = stats::binomial(link = "cloglog"),
      weights = weights
    )),
    error = function(e) NULL
  )
  if (is.null(fit) || !isTRUE(fit$converged)) {
    out$message <- "The recalibration fit did not converge."
    return(out)
  }
  coefficients <- stats::coef(fit)
  if (length(coefficients) < 2L || anyNA(coefficients)) {
    out$message <- "The recalibration fit separated the outcome or is singular."
    return(out)
  }
  out$intercept <- unname(coefficients[[1L]])
  out$slope <- unname(coefficients[[2L]])
  out
}

#' Calibration at every requested time point
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param survival A matrix of predicted survival probabilities.
#' @param time_points The time points.
#' @param groups Number of prediction groups.
#' @param censoring Censoring model.
#'
#' @return A list with `table` (a long table over time points), `summary` (one
#'   row per time point, including the `surv_calibration_slope()` coefficients),
#'   and `curves` (the per-group data used by the plot).
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, c(12, 36, 60))
#' calibration <- surv_calibration_curves(
#'   data$time, data$status, truth$survival, c(12, 36, 60)
#' )
#' calibration$summary
surv_calibration_curves <- function(time, status, survival, time_points,
                                    groups = 10L, censoring = "marginal") {
  survival <- as.matrix(survival)
  if (ncol(survival) != length(time_points)) {
    .abort("`survival` must have one column per time point.")
  }
  tables <- list()
  summaries <- list()
  for (index in seq_along(time_points)) {
    result <- surv_calibration(
      time, status, survival[, index], time_points[[index]],
      groups = groups, censoring = censoring
    )
    table <- result$table
    if (nrow(table) > 0L) {
      table$time <- time_points[[index]]
      tables[[length(tables) + 1L]] <- table
    }
    slope <- surv_calibration_slope(
      time, status, survival[, index], time_points[[index]],
      censoring = censoring
    )
    summaries[[length(summaries) + 1L]] <- data.frame(
      time = time_points[[index]],
      ici = result$ici,
      groups = nrow(table),
      slope = slope$slope,
      intercept = slope$intercept,
      n = result$n,
      events = result$events,
      message = result$message %||% slope$message %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }
  combined <- if (length(tables) == 0L) {
    data.frame(
      group = integer(), n = integer(), predicted = numeric(),
      observed = numeric(), difference = numeric(), time = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, tables)
    rownames(out) <- NULL
    out[, c("time", "group", "n", "predicted", "observed", "difference")]
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  list(table = combined, summary = summary)
}
