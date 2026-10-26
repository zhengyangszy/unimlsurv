# Bootstrap intervals for the survival metrics -------------------------------

#' Bootstrap confidence intervals for survival metrics
#'
#' Percentile bootstrap intervals for the concordance index, the integrated Brier
#' score, and every per-time Brier score and time-dependent AUC. Patients are
#' resampled with replacement, which is the unit of analysis; resamples that lose
#' all events, all censoring, or the cases or controls a time point needs
#' contribute a missing value instead of an error, and `n_effective` counts the
#' resamples that produced a number.
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param risk Optional risk score, used for the concordance index.
#' @param survival A matrix of predicted survival probabilities, or `NULL` for a
#'   risk-only model.
#' @param time_points The time points.
#' @param censoring Censoring model.
#' @param R Number of bootstrap resamples.
#' @param level Confidence level.
#' @param seed Random seed.
#'
#' @return A data frame with `metric`, `mean`, `se`, `lower`, `upper`,
#'   `n_effective`, and `formatted`.
#' @references
#' Efron B, Tibshirani RJ (1993). *An Introduction to the Bootstrap*. Chapman &
#' Hall. doi:10.1201/9780429246593
#'
#' Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
#' Development, Validation, and Updating*, 2nd edition. Springer.
#' doi:10.1007/978-3-030-16399-0
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, c(12, 36))
#' surv_ci(data$time, data$status, truth$risk, truth$survival, c(12, 36), R = 40)
surv_ci <- function(time, status, risk = NULL, survival = NULL, time_points,
                    censoring = "marginal", R = 500L, level = 0.95,
                    seed = 2026L) {
  .check_survival_inputs(time, status, survival = survival, risk = risk)
  .assert_numeric_vector(time_points, "time_points", min = 0)
  .assert_scalar_number(R, "R", integer = TRUE, min = 10)
  .assert_scalar_number(level, "level", min = 0.5, max = 0.999)
  .assert_scalar_number(seed, "seed", integer = TRUE, min = 1)
  if (!is.null(survival)) {
    survival <- as.matrix(survival)
    if (ncol(survival) != length(time_points)) {
      .abort("`survival` must have one column per time point.")
    }
  }
  point <- calc_surv_metrics(
    time, status, risk, survival, time_points, censoring = censoring
  )$row
  names_used <- names(point)
  values <- matrix(
    NA_real_,
    nrow = as.integer(R),
    ncol = length(names_used),
    dimnames = list(NULL, names_used)
  )
  n <- length(time)
  set.seed(as.integer(seed))
  for (index in seq_len(as.integer(R))) {
    sample_index <- sample.int(n, n, replace = TRUE)
    resample <- tryCatch(
      calc_surv_metrics(
        time[sample_index],
        status[sample_index],
        if (is.null(risk)) NULL else risk[sample_index],
        if (is.null(survival)) NULL else survival[sample_index, , drop = FALSE],
        time_points,
        censoring = censoring
      )$row,
      error = function(e) NULL
    )
    if (is.null(resample)) {
      next
    }
    values[index, ] <- vapply(
      names_used,
      function(name) as.numeric(resample[[name]]),
      numeric(1)
    )
  }
  .surv_interval_table(values, point, level = level)
}

# Shared percentile summary: point estimate, bootstrap mean and standard error,
# interval, and the number of resamples that produced a usable number.
.surv_interval_table <- function(values, point, level = 0.95,
                                 ignore = c("n", "n_events", "n_case",
                                            "n_control")) {
  alpha <- (1 - level) / 2
  columns <- setdiff(colnames(values), ignore)
  rows <- lapply(columns, function(name) {
    column <- values[, name]
    usable <- column[is.finite(column)]
    estimate <- as.numeric(point[[name]])
    if (length(usable) < 2L) {
      return(data.frame(
        metric = name,
        mean = estimate,
        se = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        n_effective = as.integer(length(usable)),
        formatted = if (is.na(estimate)) "NA" else formatC(estimate, format = "f", digits = 3),
        stringsAsFactors = FALSE
      ))
    }
    lower <- as.numeric(stats::quantile(usable, probs = alpha, names = FALSE))
    upper <- as.numeric(stats::quantile(usable, probs = 1 - alpha, names = FALSE))
    data.frame(
      metric = name,
      mean = estimate,
      se = stats::sd(usable),
      lower = lower,
      upper = upper,
      n_effective = as.integer(length(usable)),
      formatted = sprintf("%.3f (%.3f-%.3f)", estimate, lower, upper),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
