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
#'   row per time point), and `curves` (the per-group data used by the plot).
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
    summaries[[length(summaries) + 1L]] <- data.frame(
      time = time_points[[index]],
      ici = result$ici,
      groups = nrow(table),
      message = result$message %||% NA_character_,
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
