test_that("plot_surv_km() draws one curve per group", {
  data <- sim_surv_data(n = 200, seed = 31)
  spec <- surv_spec(time = "time", status = "status")

  overall <- plot_surv_km(data, spec)
  expect_s3_class(overall, "ggplot")
  grouped <- plot_surv_km(data, spec, group = "sex", conf_int = TRUE)
  expect_s3_class(grouped, "ggplot")
  expect_identical(length(unique(grouped$data$strata)), 2L)

  expect_error(plot_surv_km(data, spec, group = "nope"), "grouping column")
  expect_error(
    plot_surv_km(data, surv_spec(time = "time", status = "nothing")),
    "does not contain"
  )
})

test_that("the metric plots accept the tables they document", {
  data <- sim_surv_data(n = 200, seed = 32)
  truth <- sim_surv_truth(data, c(12, 36))
  curves <- surv_calibration_curves(data$time, data$status, truth$survival, c(12, 36))
  calibration_plot <- plot_surv_calibration(curves$table)
  expect_s3_class(calibration_plot, "ggplot")

  dca <- surv_dca_curves(data$time, data$status, truth$survival, c(12, 36))
  dca_plot <- plot_surv_dca(dca)
  expect_s3_class(dca_plot, "ggplot")

  expect_error(plot_surv_calibration(data.frame(x = 1)), "predicted")
  expect_error(plot_surv_dca(data.frame(x = 1)), "threshold")
})

test_that("the workflow plots use the stored predictions", {
  set.seed(33)
  data <- sim_surv_data(n = 260, seed = 33)
  index <- sample.int(nrow(data))
  spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
  spec <- surv_set_learners(spec, "coxph", tune = FALSE)
  spec <- surv_set_tuning(spec, folds = 2L, evals = 1L)
  result <- run_surv_workflow(
    spec,
    train = data[index[1:180], ],
    validation = data[index[181:260], ]
  )

  expect_s3_class(plot_surv_roc(result), "ggplot")
  expect_s3_class(plot_surv_time_metrics(result, set = "validation"), "ggplot")
  expect_s3_class(plot_surv_calibration(result), "ggplot")
  expect_s3_class(plot_surv_dca(result), "ggplot")
  expect_error(plot_surv_time_metrics(result, set = "nope"), "No per-time metrics")
  expect_error(plot_surv_roc(result, set = "nope"), "No predictions")
})
