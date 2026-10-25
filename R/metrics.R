# L5 metrics layer -----------------------------------------------------------
#
# The conventions are fixed here once, so that a table, a figure, and a paper
# cannot disagree:
#
# * `cindex` is Harrell's concordance index: the probability that a patient who
#   had the event earlier has the higher risk score, computed over comparable
#   pairs only. Higher is better.
# * `brier(t)` is Graf's inverse-probability-of-censoring weighted Brier score at
#   time `t`. Lower is better.
# * `auc(t)` is the cumulative/dynamic time-dependent AUC: a case is someone with
#   an event at or before `t`, a control is someone still event-free after `t`,
#   and both are weighted by the inverse censoring survival.
# * `ibs` is the Brier score integrated over the requested time points with the
#   trapezoidal rule.
#
# Every one of them is pinned by a test against `riskRegression::Score()` and
# `survival::concordance()`, computed on the same predictions.

.km_survival <- function(time, event, floor = 0) {
  if (length(time) == 0L) {
    return(function(u) rep(1, length(u)))
  }
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  function(u) {
    values <- summary(fit, times = u, extend = TRUE)$surv
    values[is.na(values)] <- 1
    pmax(values, floor)
  }
}

# The censoring distribution is estimated by Kaplan-Meier with the event
# indicator flipped, which is the "marginal" model of Graf (1999) and of
# `pec::pec(cens.model = "marginal")`. The floor keeps a zero-probability tail
# from producing infinite weights.
.censoring_survival <- function(time, status, censoring = "marginal") {
  .assert_choice(censoring, "censoring", c("marginal", "none"))
  if (identical(censoring, "none")) {
    return(function(u) rep(1, length(u)))
  }
  .km_survival(time, 1 - status, floor = 1e-9)
}

.check_survival_inputs <- function(time, status, survival = NULL, risk = NULL) {
  .assert_numeric_vector(time, "time", min = 0)
  if (!is.numeric(status) || !all(status %in% c(0, 1))) {
    .abort("`status` must be a 0/1 numeric vector.")
  }
  if (length(time) != length(status)) {
    .abort("`time` and `status` must have the same length.")
  }
  if (!is.null(survival)) {
    .assert_numeric_vector(survival, "survival", min = 0, max = 1)
    if (length(survival) != length(time)) {
      .abort("`survival` and `time` must have the same length.")
    }
  }
  if (!is.null(risk)) {
    .assert_numeric_vector(risk, "risk")
    if (length(risk) != length(time)) {
      .abort("`risk` and `time` must have the same length.")
    }
  }
  invisible(TRUE)
}

#' Harrell's concordance index
#'
#' @param time Follow-up times.
#' @param status Event indicator, `1` for an event and `0` for censoring.
#' @param risk Risk score; higher means shorter survival.
#'
#' @return A single numeric value, or `NA` when no pair is comparable.
#' @references
#' Harrell FE, Lee KL, Mark DB (1996). Multivariable prognostic models: issues in
#' developing models, evaluating assumptions and adequacy, and measuring and
#' reducing errors. *Statistics in Medicine*, 15(4), 361-387.
#'
#' Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011). On the C-statistics
#' for evaluating overall adequacy of risk prediction procedures with censored
#' survival data. *Statistics in Medicine*, 30(10), 1105-1117.
#' doi:10.1002/sim.4154
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, 36)
#' with(data, surv_cindex(time, status, truth$risk))
surv_cindex <- function(time, status, risk) {
  .check_survival_inputs(time, status, risk = risk)
  frame <- data.frame(time = time, status = status, risk = risk)
  estimate <- tryCatch(
    survival::concordance(
      survival::Surv(time, status) ~ risk,
      data = frame,
      reverse = TRUE
    ),
    error = function(e) e
  )
  if (inherits(estimate, "error")) {
    return(NA_real_)
  }
  value <- as.numeric(estimate$concordance)
  if (length(value) != 1L || !is.finite(value)) {
    return(NA_real_)
  }
  value
}

#' Time-dependent Brier score
#'
#' The inverse-probability-of-censoring weighted Brier score of Graf and
#' colleagues at one time point.
#'
#' @param time Follow-up times.
#' @param status Event indicator, `1` for an event and `0` for censoring.
#' @param survival Predicted survival probability at `time_point`, one value per
#'   patient.
#' @param time_point The time point at which the score is evaluated.
#' @param censoring Censoring model: `"marginal"` (Kaplan-Meier) or `"none"`.
#'
#' @return A single numeric value, or `NA` when it is not estimable.
#' @references
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). Assessment and comparison
#' of prognostic classification schemes for survival data. *Statistics in
#' Medicine*, 18(17-18), 2529-2545.
#'
#' Gerds TA, Schumacher M (2006). Consistent estimation of the expected Brier
#' score in general survival models with right-censored event times. *Biometrical
#' Journal*, 48(6), 1029-1040. doi:10.1002/bimj.200610301
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, 36)
#' surv_brier(data$time, data$status, truth$survival[, 1], time_point = 36)
surv_brier <- function(time, status, survival, time_point, censoring = "marginal") {
  .check_survival_inputs(time, status, survival = survival)
  .assert_scalar_number(time_point, "time_point", min = 0)
  if (time_point <= min(time) && !any(status == 1L & time <= time_point)) {
    return(NA_real_)
  }
  G <- .censoring_survival(time, status, censoring)
  event_at <- status == 1L & time <= time_point
  survived <- time > time_point
  weight <- ifelse(event_at, 1 / G(time), 0) + ifelse(survived, 1 / G(time_point), 0)
  observed <- ifelse(event_at, 1, 0)
  mean(weight * (observed - (1 - survival))^2)
}

#' Time-dependent area under the ROC curve
#'
#' The cumulative/dynamic definition: a case is a patient with an event at or
#' before `time_point`, a control is a patient still event-free after it, and the
#' pairs are weighted by the inverse censoring survival.
#'
#' @param time Follow-up times.
#' @param status Event indicator, `1` for an event and `0` for censoring.
#' @param survival Predicted survival probability at `time_point`.
#' @param time_point The time point at which the AUC is evaluated.
#' @param censoring Censoring model: `"marginal"` or `"none"`.
#'
#' @return A single numeric value, or `NA` when there is no case or no control.
#' @references
#' Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC curves.
#' *Biometrics*, 61(1), 92-105. doi:10.1111/j.0006-341X.2005.030814.x
#'
#' Uno H, Cai T, Tian L, Wei LJ (2007). Evaluating prediction rules for t-year
#' survivors with censored regression models. *Journal of the American
#' Statistical Association*, 102(478), 527-537. doi:10.1198/016214507000000149
#'
#' Blanche P, Dartigues JF, Jacqmin-Gadda H (2013). Estimating and comparing
#' time-dependent areas under receiver operating characteristic curves for
#' censored event times with competing risks. *Statistics in Medicine*, 32(30),
#' 5381-5397. doi:10.1002/sim.5958
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, 36)
#' surv_auc(data$time, data$status, truth$survival[, 1], time_point = 36)
surv_auc <- function(time, status, survival, time_point, censoring = "marginal") {
  .check_survival_inputs(time, status, survival = survival)
  .assert_scalar_number(time_point, "time_point", min = 0)
  case <- status == 1L & time <= time_point
  control <- time > time_point
  if (!any(case) || !any(control)) {
    return(NA_real_)
  }
  G <- .censoring_survival(time, status, censoring)
  weight_case <- 1 / G(time[case])
  weight_control <- 1 / G(rep(time_point, sum(control)))
  score_case <- survival[case]
  score_control <- survival[control]

  # Weighted Mann-Whitney statistic with ties counted as half, computed from the
  # sorted control scores so the cost stays O(n log n).
  ordering <- order(score_control)
  sorted <- score_control[ordering]
  cumulative <- cumsum(weight_control[ordering])
  total <- cumulative[length(cumulative)]
  le <- findInterval(score_case, sorted)
  lt <- findInterval(score_case, sorted, left.open = TRUE)
  weight_le <- ifelse(le > 0L, cumulative[pmax(le, 1L)], 0)
  weight_lt <- ifelse(lt > 0L, cumulative[pmax(lt, 1L)], 0)
  credit <- (total - weight_le) + 0.5 * (weight_le - weight_lt)
  denominator <- sum(weight_case) * total
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  as.numeric(sum(weight_case * credit) / denominator)
}

#' Integrated Brier score
#'
#' Integrates the time-dependent Brier score over the requested time points with
#' the trapezoidal rule, then divides by the width of the interval, so the result
#' is on the same scale as a single Brier score.
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param survival A matrix of predicted survival probabilities, one row per
#'   patient and one column per time point.
#' @param time_points The time points, matching the columns of `survival`.
#' @param censoring Censoring model.
#'
#' @return A single numeric value.
#' @references
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). Assessment and comparison
#' of prognostic classification schemes for survival data. *Statistics in
#' Medicine*, 18(17-18), 2529-2545.
#'
#' Gerds TA, Schumacher M (2006). Consistent estimation of the expected Brier
#' score in general survival models with right-censored event times. *Biometrical
#' Journal*, 48(6), 1029-1040. doi:10.1002/bimj.200610301
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, c(12, 36, 60))
#' surv_ibs(data$time, data$status, truth$survival, c(12, 36, 60))
surv_ibs <- function(time, status, survival, time_points, censoring = "marginal") {
  survival <- as.matrix(survival)
  if (ncol(survival) != length(time_points)) {
    .abort("`survival` must have one column per time point.")
  }
  scores <- vapply(seq_along(time_points), function(index) {
    surv_brier(
      time, status, survival[, index], time_points[[index]],
      censoring = censoring
    )
  }, numeric(1))
  if (length(scores) == 1L) {
    return(scores[[1L]])
  }
  keep <- is.finite(scores)
  if (sum(keep) < 2L) {
    return(NA_real_)
  }
  times <- time_points[keep]
  values <- scores[keep]
  width <- diff(range(times))
  if (width <= 0) {
    return(mean(values))
  }
  area <- sum(
    diff(times) * (utils::head(values, -1L) + utils::tail(values, -1L)) / 2
  )
  area / width
}

#' Calculate survival metrics
#'
#' The one function every table and figure goes through. It reports the
#' concordance index, the integrated Brier score, and, at each requested time
#' point, the Brier score and the cumulative/dynamic AUC.
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param risk Risk score, used for the concordance index.
#' @param survival A matrix of predicted survival probabilities, or a vector when
#'   `time_points` has length one.
#' @param time_points The time points.
#' @param censoring Censoring model.
#'
#' @return A list with `overall` (a one-row data frame), `by_time` (one row per
#'   time point), `row` (a wide one-row data frame with `cindex`, `ibs`, and one
#'   `brier_<t>`/`auc_<t>` column per time point), and `n`/`n_events`.
#' @references
#' Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
#' Development, Validation, and Updating*, 2nd edition. Springer.
#' doi:10.1007/978-3-030-16399-0
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, c(12, 36, 60))
#' metrics <- calc_surv_metrics(
#'   data$time, data$status, truth$risk, truth$survival, c(12, 36, 60)
#' )
#' metrics$overall
calc_surv_metrics <- function(time, status, risk = NULL, survival = NULL,
                              time_points, censoring = "marginal") {
  .assert_numeric_vector(time_points, "time_points", min = 0)
  if (!is.null(survival)) {
    survival <- as.matrix(survival)
    if (ncol(survival) != length(time_points)) {
      if (length(survival) == length(time) && length(time_points) == 1L) {
        survival <- matrix(survival, ncol = 1L)
      } else {
        .abort("`survival` must have one column per time point.")
      }
    }
  }
  if (is.null(survival)) {
    # A learner that predicts only a risk score still has a concordance index;
    # everything that needs survival probabilities is reported as unavailable.
    by_time <- data.frame(
      time = as.numeric(time_points),
      brier = NA_real_,
      auc = NA_real_,
      n_case = vapply(time_points, function(t0) sum(status == 1L & time <= t0),
                      integer(1)),
      n_control = vapply(time_points, function(t0) sum(time > t0), integer(1)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    overall <- data.frame(
      n = length(time),
      n_events = sum(status == 1L),
      cindex = if (is.null(risk)) NA_real_ else surv_cindex(time, status, risk),
      ibs = NA_real_,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    row <- overall
    for (index in seq_len(nrow(by_time))) {
      label <- format(by_time$time[[index]], trim = TRUE)
      row[[paste0("brier_", label)]] <- NA_real_
      row[[paste0("auc_", label)]] <- NA_real_
    }
    return(list(overall = overall, by_time = by_time, row = row))
  }
  by_time <- data.frame(
    time = as.numeric(time_points),
    brier = vapply(seq_along(time_points), function(index) {
      surv_brier(time, status, survival[, index], time_points[[index]],
                 censoring = censoring)
    }, numeric(1)),
    auc = vapply(seq_along(time_points), function(index) {
      surv_auc(time, status, survival[, index], time_points[[index]],
               censoring = censoring)
    }, numeric(1)),
    n_case = vapply(time_points, function(t0) {
      sum(status == 1L & time <= t0)
    }, integer(1)),
    n_control = vapply(time_points, function(t0) {
      sum(time > t0)
    }, integer(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  cindex <- if (is.null(risk)) NA_real_ else surv_cindex(time, status, risk)
  ibs <- surv_ibs(time, status, survival, time_points, censoring = censoring)
  overall <- data.frame(
    n = length(time),
    n_events = sum(status == 1L),
    cindex = cindex,
    ibs = ibs,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  row <- overall
  for (index in seq_len(nrow(by_time))) {
    label <- format(by_time$time[[index]], trim = TRUE)
    row[[paste0("brier_", label)]] <- by_time$brier[[index]]
    row[[paste0("auc_", label)]] <- by_time$auc[[index]]
  }
  list(overall = overall, by_time = by_time, row = row)
}
