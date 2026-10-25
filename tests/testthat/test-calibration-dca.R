test_that("calibration separates a good predictor from a miscalibrated one", {
  data <- sim_surv_data(n = 300, seed = 6)
  truth <- sim_surv_truth(data, 36)$survival[, 1]

  good <- surv_calibration(data$time, data$status, truth, 36, groups = 5L)
  inverted <- surv_calibration(data$time, data$status, 1 - truth, 36, groups = 5L)

  expect_true(is.finite(good$ici))
  expect_true(is.finite(inverted$ici))
  expect_gt(inverted$ici, good$ici)
  expect_equal(sum(good$table$n), nrow(data))
  expect_true(all(good$table$n > 0L))
})

test_that("calibration reports observed Kaplan-Meier survival per group", {
  data <- sim_surv_data(n = 300, seed = 7)
  truth <- sim_surv_truth(data, 36)
  calibration <- surv_calibration(
    data$time, data$status, truth$survival[, 1], 36, groups = 5L
  )
  groups <- calibration$table
  expect_identical(nrow(groups), 5L)
  expect_true(all(groups$observed >= 0 & groups$observed <= 1))
  expect_true(all(groups$predicted >= 0 & groups$predicted <= 1))
  # Ordered by predicted survival, the groups must be ordered too.
  expect_false(is.unsorted(groups$predicted))
  expect_equal(calibration$ici, sum(groups$n * abs(groups$difference)) / sum(groups$n))
})

test_that("calibration refuses to invent groups it cannot form", {
  time <- c(5, 6, 7, 8)
  status <- c(1L, 1L, 1L, 1L)
  result <- surv_calibration(time, status, rep(0.5, 4), 6)
  expect_true(is.na(result$ici))
  expect_true(!is.null(result$message))
  expect_identical(nrow(result$table), 0L)
})

test_that("surv_calibration_curves() stacks the time points", {
  data <- sim_surv_data(n = 200, seed = 8)
  truth <- sim_surv_truth(data, c(12, 36))
  curves <- surv_calibration_curves(data$time, data$status, truth$survival, c(12, 36))
  expect_setequal(unique(curves$table$time), c(12, 36))
  expect_identical(nrow(curves$summary), 2L)
  expect_true(all(c("ici", "groups", "message") %in% names(curves$summary)))
})

test_that("decision curves obey their definitions", {
  data <- sim_surv_data(n = 300, seed = 9)
  truth <- sim_surv_truth(data, 36)
  time_point <- 36
  curve <- surv_dca(data$time, data$status, truth$survival[, 1], time_point)

  model <- curve[curve$strategy == "model", , drop = FALSE]
  treat_all <- curve[curve$strategy == "treat_all", , drop = FALSE]
  treat_none <- curve[curve$strategy == "treat_none", , drop = FALSE]

  expect_identical(nrow(model), nrow(treat_all))
  expect_true(all(treat_none$net_benefit == 0))
  # Treating everybody equals the model at a threshold of zero.
  expect_equal(model$net_benefit[[1L]], treat_all$net_benefit[[1L]], tolerance = 1e-10)
  # At very low thresholds almost everybody is flagged.
  expect_gt(model$net_benefit[[1L]], 0)
  # Net benefit falls as the threshold rises.
  expect_true(model$net_benefit[[1L]] >= utils::tail(model$net_benefit, 1L))

  # The model curve equals the treat-all curve when every risk is above the cut.
  flat <- rep(0.01, nrow(data))
  flat_curve <- surv_dca(data$time, data$status, flat, time_point,
                         thresholds = c(0.5, 0.9))
  expect_equal(
    flat_curve$net_benefit[flat_curve$strategy == "model"],
    flat_curve$net_benefit[flat_curve$strategy == "treat_all"],
    tolerance = 1e-10
  )
})

test_that("surv_dca_curves() stacks the time points and keeps labels", {
  data <- sim_surv_data(n = 200, seed = 10)
  truth <- sim_surv_truth(data, c(12, 36))
  curves <- surv_dca_curves(data$time, data$status, truth$survival, c(12, 36))
  expect_setequal(unique(curves$time), c(12, 36))
  expect_setequal(unique(curves$strategy), c("model", "treat_all", "treat_none"))
  expect_true(all(c("time", "threshold", "strategy", "net_benefit") %in% names(curves)))
})

test_that("the simulated cohort is reproducible and has a known truth", {
  first <- sim_surv_data(n = 150, seed = 11)
  second <- sim_surv_data(n = 150, seed = 11)
  expect_identical(first$time, second$time)
  expect_identical(first$status, second$status)
  expect_false(identical(first$time, sim_surv_data(n = 150, seed = 12)$time))

  truth <- sim_surv_truth(first, c(12, 36, 60))
  expect_identical(dim(truth$survival), c(150L, 3L))
  expect_true(all(truth$survival >= 0 & truth$survival <= 1))
  # Survival is non-increasing in time.
  expect_true(all(truth$survival[, 3] <= truth$survival[, 1]))
  # The true risk is a real predictor.
  expect_gt(surv_cindex(first$time, first$status, truth$risk), 0.6)

  expect_error(sim_surv_truth(data.frame(time = 1, status = 1), 12), "sim_surv_data")
  expect_error(sim_surv_data(n = 5), "n")
  expect_error(sim_surv_data(missing = 2), "missing")
})
