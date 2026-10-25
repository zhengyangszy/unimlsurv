# Simulated cohort -----------------------------------------------------------

#' Simulate a survival cohort
#'
#' Generates a cohort with a known linear predictor, so every statistic in this
#' package can be checked against the truth rather than against another
#' implementation of the same assumption. Event times follow a Weibull model with
#' proportional hazards, censoring is exponential, and administrative censoring
#' can be added so that the cohort looks like a trial with a fixed follow-up.
#'
#' The true linear predictor is stored as the attribute `"eta"`, and
#' [sim_surv_truth()] turns it into the true risk score and the true survival
#' probabilities at any time points.
#'
#' @param n Number of patients.
#' @param beta Named vector of log hazard ratios for the continuous features.
#' @param shape Weibull shape; `1` gives exponential event times.
#' @param scale Weibull scale; larger values give longer survival.
#' @param censor_rate Rate of the exponential censoring distribution. Set to `0`
#'   to simulate without censoring.
#' @param admin_time Optional administrative censoring time.
#' @param missing Proportion of feature values to set to `NA`, for demonstrating
#'   imputation.
#' @param seed Random seed.
#'
#' @return A data frame with `time`, `status`, and the feature columns
#'   `x1`, `x2`, ..., `x<k>`, and `sex`.
#' @references
#' Bender R, Augustin T, Blettner M (2005). Generating survival times to simulate
#' Cox proportional hazards models. *Statistics in Medicine*, 24(11), 1713-1723.
#' doi:10.1002/sim.2059
#'
#' Burton A, Altman DG, Royston P, Holder RL (2006). The design of simulation
#' studies in medical statistics. *Statistics in Medicine*, 25(24), 4279-4292.
#' doi:10.1002/sim.2328
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' head(data)
#' attr(data, "eta")[1:3]
sim_surv_data <- function(n = 500L, beta = c(x1 = 0.8, x2 = -0.4),
                          shape = 1, scale = 12, censor_rate = 1 / 18,
                          admin_time = NULL, missing = 0, seed = 2026L) {
  .assert_scalar_number(n, "n", integer = TRUE, min = 10)
  .assert_numeric_vector(beta, "beta")
  if (is.null(names(beta))) {
    names(beta) <- paste0("x", seq_along(beta))
  }
  .assert_scalar_number(shape, "shape", min = 1e-08)
  .assert_scalar_number(scale, "scale", min = 1e-08)
  .assert_scalar_number(censor_rate, "censor_rate", min = 0)
  .assert_scalar_number(missing, "missing", min = 0, max = 0.8)
  .assert_scalar_number(seed, "seed", integer = TRUE, min = 1)
  if (!is.null(admin_time)) {
    .assert_scalar_number(admin_time, "admin_time", min = 1e-08)
  }

  set.seed(as.integer(seed))
  k <- length(beta)
  features <- matrix(
    stats::rnorm(n * k),
    nrow = n,
    ncol = k,
    dimnames = list(NULL, names(beta))
  )
  sex <- factor(sample(c("female", "male"), n, replace = TRUE))
  eta <- as.numeric(features %*% beta) + ifelse(sex == "male", 0.3, 0)

  event_time <- scale * (-log(stats::runif(n)) / exp(eta))^(1 / shape)
  if (censor_rate > 0) {
    censor_time <- stats::rexp(n, rate = censor_rate)
  } else {
    censor_time <- rep(Inf, n)
  }
  follow_up <- if (is.null(admin_time)) Inf else admin_time
  observed <- pmin(event_time, censor_time, follow_up)
  status <- as.integer(event_time <= pmin(censor_time, follow_up))

  data <- data.frame(
    time = as.numeric(observed),
    status = status,
    features,
    sex = sex,
    stringsAsFactors = FALSE
  )
  if (missing > 0 && k >= 1L) {
    target <- c("x1", "x2")[seq_len(min(2L, k))]
    for (column in target) {
      index <- sample.int(n, size = max(1L, floor(missing * n)))
      data[[column]][index] <- NA_real_
    }
  }
  attr(data, "eta") <- eta
  attr(data, "shape") <- shape
  attr(data, "scale") <- scale
  attr(data, "beta") <- beta
  data
}

#' True risk and true survival of a simulated cohort
#'
#' @param data A data frame created by [sim_surv_data()].
#' @param time_points Time points at which the true survival is evaluated.
#'
#' @return A list with `risk` (the true linear predictor, higher meaning shorter
#'   survival) and `survival` (a matrix with one row per patient and one column
#'   per time point).
#' @export
#' @examples
#' data <- sim_surv_data(n = 100, seed = 1)
#' truth <- sim_surv_truth(data, c(12, 36))
#' dim(truth$survival)
sim_surv_truth <- function(data, time_points) {
  eta <- attr(data, "eta")
  shape <- attr(data, "shape")
  scale <- attr(data, "scale")
  if (is.null(eta) || is.null(shape) || is.null(scale)) {
    .abort("`data` must come from `sim_surv_data()`.")
  }
  .assert_numeric_vector(time_points, "time_points", min = 0)
  survival <- outer(
    eta,
    time_points,
    function(e, t) exp(-(t / scale)^shape * exp(e))
  )
  colnames(survival) <- paste0("t", format(time_points, trim = TRUE))
  list(risk = eta, survival = survival)
}
