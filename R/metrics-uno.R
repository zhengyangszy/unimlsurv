# L5 metrics layer · Uno estimators -------------------------------------------
#
# Trapezoidal integral over the time grid, divided by the grid width, so the
# result stays on the same scale as a single per-time value (the same
# convention as `surv_ibs()`).
.trapz_scaled <- function(time_points, values) {
  if (length(values) == 1L) {
    return(values[[1L]])
  }
  keep <- is.finite(values)
  if (sum(keep) < 2L) {
    return(NA_real_)
  }
  times <- time_points[keep]
  vals <- values[keep]
  width <- diff(range(times))
  if (width <= 0) {
    return(mean(vals))
  }
  area <- sum(
    diff(times) * (utils::head(vals, -1L) + utils::tail(vals, -1L)) / 2
  )
  area / width
}

# Uno's IPCW concordance index and time-dependent AUC -------------------------
#
# Harrell's C is computed over comparable pairs of the evaluation set, so its
# null behaviour depends on the censoring distribution of that set. Uno and
# colleagues reweight the same pairs by the inverse of the censoring survival
# estimated on the **training** data, which makes the estimator robust to a
# censoring distribution that differs between development and validation.
# The censoring distribution must therefore come from the training set: that is
# part of the definition, not an implementation detail.
#
# Both estimators are implemented here from the definitions in the papers and
# pinned by tests against `survAUC::UnoC()` and `survAUC::AUC.uno()`, which use
# independent C implementations of the same estimators.

#' Uno's concordance index
#'
#' The inverse-probability-of-censoring weighted concordance index of Uno and
#' colleagues at a truncation horizon `tau`: among patients of the evaluation
#' set who had an event by `tau`, pairs are formed with every patient followed
#' at least that long, and each pair is weighted by the inverse censoring
#' survival of the earlier patient, estimated on the **training** data.
#' Higher is better, 0.5 is a coin flip.
#'
#' @param time_train Follow-up times of the training data.
#' @param status_train Event indicators of the training data (`1` = event).
#' @param time_test Follow-up times of the evaluation data.
#' @param status_test Event indicators of the evaluation data.
#' @param risk_test Risk scores of the evaluation data; higher means shorter
#'   survival.
#' @param tau Truncation horizon. Defaults to the largest observed time in the
#'   evaluation data, as in `survAUC::UnoC()`. The censoring survival is
#'   floored at `1e-9`, so a `tau` beyond the training follow-up still runs but
#'   the weights become extreme; keep `tau` within the training follow-up.
#'
#' @return A single numeric value, or `NA` when no comparable pair exists.
#' @references
#' Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011). On the C-statistics
#' for evaluating overall adequacy of risk prediction procedures with censored
#' survival data. *Statistics in Medicine*, 30(10), 1105-1117.
#' doi:10.1002/sim.4154
#'
#' Potapov S, Adler W, Schmid M, Bertrand F (2025). survAUC: Estimators of
#' Prediction Accuracy for Time-to-Event Data. R package version 1.4-0.
#' doi:10.32614/CRAN.package.survAUC
#' @export
#' @examples
#' data <- sim_surv_data(n = 400, seed = 1)
#' index <- sample.int(nrow(data), 200)
#' truth <- sim_surv_truth(data, 36)
#' surv_cindex_uno(
#'   time_train = data$time[-index], status_train = data$status[-index],
#'   time_test = data$time[index], status_test = data$status[index],
#'   risk_test = truth$risk[index]
#' )
surv_cindex_uno <- function(time_train, status_train, time_test, status_test,
                            risk_test, tau = NULL) {
  .check_survival_inputs(time_train, status_train)
  .check_survival_inputs(time_test, status_test, risk = risk_test)
  if (is.null(tau)) {
    tau <- max(time_test)
  } else {
    .assert_scalar_number(tau, "tau", min = 0)
  }
  # The training censoring survival, evaluated as a step function so that
  # times beyond the training follow-up keep the last value of the curve
  # (the way `survAUC` evaluates it). Cases whose weight is not finite
  # (G = 0, which happens when `tau` runs past the training follow-up) are
  # excluded, matching the `G[i] > 0` condition of the reference
  # implementation.
  censor_fit <- survival::survfit(survival::Surv(time_train, 1 - status_train) ~ 1)
  censor_times <- censor_fit$time
  censor_surv <- censor_fit$surv
  G_at <- function(u) {
    index <- findInterval(u, censor_times)
    ifelse(index == 0L, 1, censor_surv[pmax(index, 1L)])
  }
  # Uno's weights are the inverse *squared* censoring survival at each event
  # patient's own observed time, risk ties are discarded (not counted as
  # half), and the pair pool is every evaluation patient observed strictly
  # later: this is the convention `survAUC::UnoC()` implements, and the tests
  # pin it exactly.
  event_index <- which(status_test == 1L & time_test < tau)
  if (length(event_index) == 0L) {
    return(NA_real_)
  }
  censoring_at <- G_at(time_test[event_index])
  keep <- censoring_at > 0
  event_index <- event_index[keep]
  censoring_at <- censoring_at[keep]
  if (length(event_index) == 0L) {
    return(NA_real_)
  }
  weight_event <- 1 / censoring_at^2
  numerator <- 0
  denominator <- 0
  for (i in seq_along(event_index)) {
    pool <- which(time_test > time_test[event_index[i]])
    if (length(pool) == 0L) {
      next
    }
    risk_pool <- risk_test[pool]
    risk_i <- risk_test[event_index[i]]
    concordant <- sum(risk_pool < risk_i)
    weight <- weight_event[[i]]
    numerator <- numerator + weight * concordant
    denominator <- denominator + weight * length(pool)
  }
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  numerator / denominator
}

#' Uno's time-dependent area under the ROC curve
#'
#' The inverse-probability-of-censoring weighted cumulative/dynamic AUC of Uno
#' and colleagues at one time point: cases are patients of the evaluation set
#' with an event by `time_point`, controls are patients still event-free after
#' it, and the case weights are the inverse censoring survival at `time_point`,
#' estimated on the **training** data. This answers the same question as
#' [surv_auc()] with a censoring estimate that does not depend on the
#' evaluation set.
#'
#' @param time_train Follow-up times of the training data.
#' @param status_train Event indicators of the training data.
#' @param time_test Follow-up times of the evaluation data.
#' @param status_test Event indicators of the evaluation data.
#' @param survival_test Predicted survival probability at `time_point`, one
#'   value per evaluation patient.
#' @param time_point The time point at which the AUC is evaluated.
#'
#' @return A single numeric value, or `NA` when there is no case or no control.
#' @references
#' Uno H, Cai T, Tian L, Wei LJ (2007). Evaluating prediction rules for t-year
#' survivors with censored regression models. *Journal of the American
#' Statistical Association*, 102(478), 527-537.
#' doi:10.1198/016214507000000149
#'
#' Potapov S, Adler W, Schmid M, Bertrand F (2025). survAUC: Estimators of
#' Prediction Accuracy for Time-to-Event Data. R package version 1.4-0.
#' doi:10.32614/CRAN.package.survAUC
#' @export
#' @examples
#' data <- sim_surv_data(n = 400, seed = 1)
#' index <- sample.int(nrow(data), 200)
#' truth <- sim_surv_truth(data, 36)
#' surv_auc_uno(
#'   time_train = data$time[-index], status_train = data$status[-index],
#'   time_test = data$time[index], status_test = data$status[index],
#'   survival_test = truth$survival[index, 1], time_point = 36
#' )
surv_auc_uno <- function(time_train, status_train, time_test, status_test,
                         survival_test, time_point) {
  .check_survival_inputs(time_train, status_train)
  .check_survival_inputs(time_test, status_test, survival = survival_test)
  .assert_scalar_number(time_point, "time_point", min = 0)
  case <- status_test == 1L & time_test <= time_point
  control <- time_test > time_point
  if (!any(case) || !any(control)) {
    return(NA_real_)
  }
  G <- .censoring_survival(time_train, status_train, "marginal")
  # Each case is weighted by the inverse training censoring survival at its
  # own observed time (not at the horizon); controls are unweighted. This is
  # the estimator `survAUC::AUC.uno()` implements via its sensitivity curve.
  weight_case <- 1 / G(time_test[case])
  score_case <- survival_test[case]
  score_control <- survival_test[control]

  # Weighted Mann-Whitney statistic with ties counted as half; controls are
  # unweighted, cases carry 1/G(t) from the training censoring model.
  ordering <- order(score_control)
  sorted <- score_control[ordering]
  cumulative <- cumsum(rep(1, length(sorted)))
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

#' Integrated Uno's time-dependent AUC
#'
#' Integrates [surv_auc_uno()] over the requested time points with the
#' trapezoidal rule and divides by the width of the interval, on the same scale
#' as a single AUC. Like [surv_ibs()], it therefore depends on the grid the
#' user chose, and the documentation says so.
#'
#' @param time_train Follow-up times of the training data.
#' @param status_train Event indicators of the training data.
#' @param time_test Follow-up times of the evaluation data.
#' @param status_test Event indicators of the evaluation data.
#' @param survival A matrix of predicted survival probabilities (evaluation
#'   patients × time points).
#' @param time_points The time points.
#'
#' @return A single numeric value.
#' @references
#' Uno H, Cai T, Tian L, Wei LJ (2007). Evaluating prediction rules for t-year
#' survivors with censored regression models. *Journal of the American
#' Statistical Association*, 102(478), 527-537.
#' doi:10.1198/016214507000000149
#' @export
#' @examples
#' data <- sim_surv_data(n = 400, seed = 1)
#' index <- sample.int(nrow(data), 200)
#' truth <- sim_surv_truth(data, c(12, 36, 60))
#' surv_iauc_uno(
#'   time_train = data$time[-index], status_train = data$status[-index],
#'   time_test = data$time[index], status_test = data$status[index],
#'   survival = truth$survival[index, ], time_points = c(12, 36, 60)
#' )
surv_iauc_uno <- function(time_train, status_train, time_test, status_test,
                          survival, time_points) {
  survival <- as.matrix(survival)
  if (ncol(survival) != length(time_points)) {
    .abort("`survival` must have one column per time point.")
  }
  scores <- vapply(seq_along(time_points), function(index) {
    surv_auc_uno(
      time_train, status_train, time_test, status_test,
      survival[, index], time_points[[index]]
    )
  }, numeric(1))
  .trapz_scaled(time_points, scores)
}
