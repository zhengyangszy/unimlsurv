test_that("surv_ci() brackets the point estimate for every metric", {
  data <- sim_surv_data(n = 200, seed = 41)
  truth <- sim_surv_truth(data, c(12, 36))
  ci <- surv_ci(data$time, data$status, truth$risk, truth$survival, c(12, 36),
                R = 60, seed = 5)

  expect_true(all(c("metric", "mean", "se", "lower", "upper", "n_effective",
                    "formatted") %in% names(ci)))
  expect_true(all(c("cindex", "ibs", "brier_12", "auc_12", "brier_36", "auc_36")
                  %in% ci$metric))
  usable <- ci[is.finite(ci$lower), , drop = FALSE]
  expect_true(all(usable$lower <= usable$upper))
  expect_true(all(usable$n_effective > 0L))
  # The reported mean is the point estimate, not the bootstrap average.
  point <- calc_surv_metrics(data$time, data$status, truth$risk, truth$survival,
                             c(12, 36))$row
  expect_equal(ci$mean[ci$metric == "cindex"], point$cindex, tolerance = 1e-12)
  expect_equal(ci$mean[ci$metric == "ibs"], point$ibs, tolerance = 1e-12)
})

test_that("surv_ci() validates its inputs", {
  data <- sim_surv_data(n = 80, seed = 42)
  truth <- sim_surv_truth(data, 12)
  expect_error(
    surv_ci(data$time, data$status, truth$risk, truth$survival, c(12, 36), R = 20),
    "one column per time point"
  )
  expect_error(
    surv_ci(data$time, data$status, truth$risk, truth$survival, 12, R = 5),
    "R"
  )
  expect_error(
    surv_ci(data$time, data$status, truth$risk, truth$survival, 12, level = 1.2),
    "level"
  )
})

test_that("a workflow reports intervals for the comparison set", {
  set.seed(43)
  data <- sim_surv_data(n = 240, seed = 43)
  index <- sample.int(nrow(data))
  spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
  spec <- surv_set_learners(spec, "coxph", tune = FALSE)
  spec <- surv_set_tuning(spec, folds = 2L, evals = 1L)

  result <- run_surv_workflow(
    spec,
    train = data[index[1:180], ],
    validation = data[index[181:240], ],
    ci_R = 30L
  )
  expect_s3_class(result$ci, "data.frame")
  expect_identical(unique(result$ci$set), "validation")
  expect_identical(unique(result$ci$model), "coxph")
  expect_true("cindex" %in% result$ci$metric)

  # `ci_R = 0` skips them.
  skipped <- run_surv_workflow(
    spec,
    train = data[index[1:180], ],
    validation = data[index[181:240], ],
    ci_R = 0L
  )
  expect_null(skipped$ci)
})
