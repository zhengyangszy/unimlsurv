test_that("recalibration recovers slope 1 and intercept 0 on the true model", {
  data <- sim_surv_data(n = 4000, seed = 21, censor_rate = 0.02, admin_time = 120)
  time_point <- 20
  truth <- sim_surv_truth(data, time_point)
  coefficients <- surv_calibration_slope(
    data$time, data$status, truth$survival[, 1], time_point
  )
  expect_null(coefficients$message)
  expect_true(is.finite(coefficients$slope))
  expect_equal(coefficients$slope, 1, tolerance = 0.1)
  expect_equal(coefficients$intercept, 0, tolerance = 0.1)
  expect_identical(coefficients$n, nrow(data))
  expect_identical(coefficients$events, sum(data$status == 1L & data$time <= time_point))
})

test_that("the recalibration slope follows the direction of the miscalibration", {
  data <- sim_surv_data(n = 4000, seed = 22, censor_rate = 0.02, admin_time = 120)
  time_point <- 20
  truth <- sim_surv_truth(data, time_point)
  event <- 1 - truth$survival[, 1]
  stretch <- function(scale) {
    1 - stats::plogis(scale * stats::qlogis(pmin(pmax(event, 1e-8), 1 - 1e-8)))
  }
  over <- surv_calibration_slope(data$time, data$status, stretch(1.6), time_point)
  under <- surv_calibration_slope(data$time, data$status, stretch(0.5), time_point)
  # Over-dispersion means the predictions are too extreme: slope below 1.
  expect_lt(over$slope, 0.95)
  # Under-dispersion means they are too timid: slope above 1.
  expect_gt(under$slope, 1.05)
  # Both are fixed by raising or lowering the level, so the sign is diagnostic.
  expect_lt(over$intercept, 0)
  expect_gt(under$intercept, 0)
})

test_that("a constant prediction reports only the level", {
  data <- sim_surv_data(n = 600, seed = 23)
  coefficients <- surv_calibration_slope(
    data$time, data$status, rep(0.7, nrow(data)), 24
  )
  expect_true(is.na(coefficients$slope))
  expect_true(is.finite(coefficients$intercept))
  expect_match(coefficients$message, "does not vary")
})

test_that("recalibration refuses to fit on too little information", {
  data <- sim_surv_data(n = 300, seed = 24)
  # Every patient has an event before the time point: there is no control.
  only_events <- surv_calibration_slope(
    data$time, rep(1L, nrow(data)), rep(0.5, nrow(data)), max(data$time) + 1
  )
  expect_true(is.na(only_events$slope))
  expect_match(only_events$message, "cases and controls")

  small <- surv_calibration_slope(
    data$time[1:10], data$status[1:10], rep(0.5, 10),
    mean(range(data$time[1:10]))
  )
  expect_true(is.na(small$slope))
  expect_match(small$message, "Too few")
})

test_that("calibration curves carry the coefficients next to the ICI", {
  data <- sim_surv_data(n = 800, seed = 25, censor_rate = 0.02, admin_time = 120)
  time_points <- c(12, 20)
  truth <- sim_surv_truth(data, time_points)
  curves <- surv_calibration_curves(
    data$time, data$status, truth$survival, time_points
  )
  expect_true(all(c("ici", "slope", "intercept", "n", "events") %in%
                    names(curves$summary)))
  expect_identical(nrow(curves$summary), 2L)
  expect_true(all(is.finite(curves$summary$slope)))
  expect_true(all(abs(curves$summary$slope - 1) < 0.2))
})

test_that("calc_surv_metrics() and surv_ci() report the coefficients", {
  data <- sim_surv_data(n = 600, seed = 26, censor_rate = 0.02, admin_time = 120)
  time_points <- c(12, 20)
  truth <- sim_surv_truth(data, time_points)
  metrics <- calc_surv_metrics(
    data$time, data$status, truth$risk, truth$survival, time_points
  )
  expect_true(all(c("calib_slope", "calib_intercept") %in% names(metrics$by_time)))
  expect_true(all(c("calib_slope_12", "calib_intercept_12",
                    "calib_slope_20", "calib_intercept_20") %in%
                    names(metrics$row)))
  expect_equal(
    metrics$row$calib_slope_12,
    metrics$by_time$calib_slope[[1L]],
    tolerance = 1e-12
  )

  # A risk-only model cannot be calibrated: the columns stay, the values do not.
  risk_only <- calc_surv_metrics(
    data$time, data$status, truth$risk, NULL, time_points
  )
  expect_true(all(is.na(risk_only$by_time$calib_slope)))
  expect_true(all(is.na(risk_only$row$calib_slope_12)))

  intervals <- surv_ci(
    data$time, data$status, truth$risk, truth$survival, time_points, R = 30
  )
  expect_true(all(c("calib_slope_12", "calib_intercept_20") %in% intervals$metric))
  slope_row <- intervals[intervals$metric == "calib_slope_12", , drop = FALSE]
  expect_true(is.finite(slope_row$lower))
  expect_true(is.finite(slope_row$upper))
  expect_lt(slope_row$lower, slope_row$upper)
})

test_that("the calibration plot can annotate the coefficients", {
  data <- sim_surv_data(n = 400, seed = 27)
  time_points <- c(12, 24)
  truth <- sim_surv_truth(data, time_points)
  curves <- surv_calibration_curves(
    data$time, data$status, truth$survival, time_points
  )
  annotated <- plot_surv_calibration(curves, annotate = TRUE)
  plain <- plot_surv_calibration(curves, annotate = FALSE)
  expect_s3_class(annotated, "ggplot")
  expect_match(annotated$labels$caption, "ICI")
  expect_match(annotated$labels$caption, "slope")
  expect_null(plain$labels$caption)
  # A bare table has nothing to annotate from, and that is not an error.
  expect_null(plot_surv_calibration(curves$table, annotate = TRUE)$labels$caption)
})
