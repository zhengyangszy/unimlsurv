# L6 decision-curve layer ----------------------------------------------------

#' Time-point decision-curve analysis
#'
#' Net benefit of using a survival model to decide who gets an intervention, at
#' one time point and over a range of threshold probabilities. Events and
#' non-events are weighted by the inverse censoring survival, so a patient who was
#' censored before the time point does not count as a certain non-event.
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param survival Predicted survival probability at `time_point`.
#' @param time_point The time point.
#' @param thresholds Threshold probabilities defining "treat".
#' @param censoring Censoring model.
#'
#' @return A data frame with `threshold`, `strategy`, and `net_benefit`.
#' @references
#' Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel method for
#' evaluating prediction models. *Medical Decision Making*, 26(6), 565-574.
#' doi:10.1177/0272989X06295361
#'
#' Vickers AJ, Cronin AM, Elkin EB, Gonen M (2008). Extensions to decision curve
#' analysis, a novel method for evaluating diagnostic tests, prediction models and
#' molecular markers. *BMC Medical Informatics and Decision Making*, 8, 53.
#' doi:10.1186/1472-6947-8-53
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, 36)
#' curve <- surv_dca(data$time, data$status, truth$survival[, 1], time_point = 36)
#' head(curve)
surv_dca <- function(time, status, survival, time_point,
                     thresholds = seq(0.01, 0.99, by = 0.01),
                     censoring = "marginal") {
  .check_survival_inputs(time, status, survival = survival)
  .assert_scalar_number(time_point, "time_point", min = 0)
  .assert_numeric_vector(thresholds, "thresholds", min = 0, max = 1)
  thresholds <- sort(unique(as.numeric(thresholds)))

  event_at <- status == 1L & time <= time_point
  survived <- time > time_point
  G <- .censoring_survival(time, status, censoring)
  weight_event <- ifelse(event_at, 1 / G(time), 0)
  weight_survived <- ifelse(survived, 1 / G(time_point), 0)
  n <- length(time)
  risk <- 1 - survival

  rows <- lapply(thresholds, function(threshold) {
    if (threshold >= 1) {
      return(NULL)
    }
    odds <- threshold / (1 - threshold)
    flagged <- risk >= threshold
    true_positive <- sum(weight_event[flagged]) / n
    false_positive <- sum(weight_survived[flagged]) / n
    model <- true_positive - false_positive * odds
    all_positive <- sum(weight_event) / n - sum(weight_survived) / n * odds
    data.frame(
      threshold = c(threshold, threshold, threshold),
      strategy = c("model", "treat_all", "treat_none"),
      net_benefit = c(model, all_positive, 0),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(data.frame(
      threshold = numeric(), strategy = character(), net_benefit = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Decision curves at every requested time point
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param survival A matrix of predicted survival probabilities.
#' @param time_points The time points.
#' @param thresholds Threshold probabilities.
#' @param censoring Censoring model.
#'
#' @return A data frame with `time`, `threshold`, `strategy`, and `net_benefit`.
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, c(12, 36))
#' curves <- surv_dca_curves(data$time, data$status, truth$survival, c(12, 36))
#' head(curves)
surv_dca_curves <- function(time, status, survival, time_points,
                            thresholds = seq(0.01, 0.99, by = 0.01),
                            censoring = "marginal") {
  survival <- as.matrix(survival)
  if (ncol(survival) != length(time_points)) {
    .abort("`survival` must have one column per time point.")
  }
  rows <- lapply(seq_along(time_points), function(index) {
    curve <- surv_dca(
      time, status, survival[, index], time_points[[index]],
      thresholds = thresholds, censoring = censoring
    )
    if (nrow(curve) == 0L) {
      return(NULL)
    }
    curve$time <- time_points[[index]]
    curve[, c("time", "threshold", "strategy", "net_benefit")]
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(data.frame(
      time = numeric(), threshold = numeric(), strategy = character(),
      net_benefit = numeric(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
