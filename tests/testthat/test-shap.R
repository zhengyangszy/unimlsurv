shap_fixture <- function() {
  set.seed(51)
  data <- sim_surv_data(n = 200, seed = 51)
  index <- sample.int(nrow(data))
  spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
  spec <- surv_set_learners(spec, "coxph", tune = FALSE)
  spec <- surv_set_tuning(spec, folds = 2L, evals = 1L)
  result <- run_surv_workflow(
    spec,
    train = data[index[1:150], ],
    validation = data[index[151:200], ],
    ci_R = 0L
  )
  list(data = data, result = result, validation = data[index[151:200], ])
}

test_that("surv_shap() explains the risk score by default", {
  skip_if_not_installed("kernelshap")
  skip_if_not_installed("shapviz")
  fixture <- shap_fixture()

  explained <- surv_shap(
    fixture$result, fixture$validation, model = "coxph",
    n_explain = 20L, bg_n = 10L, seed = 1
  )
  expect_s3_class(explained, "unimlsurv_shap")
  expect_identical(explained$model, "coxph")
  expect_identical(explained$quantity, "risk score")
  expect_null(explained$time_point)
  expect_identical(nrow(explained$data), 20L)
  expect_true(all(c("x1", "x2") %in% explained$feature_names))

  values <- shapviz::get_shap_values(explained$shap)
  expect_identical(nrow(values), 20L)
  # SHAP values add up to the prediction minus the baseline.
  expect_equal(
    rowSums(values) + as.numeric(shapviz::get_baseline(explained$shap)),
    explained$values,
    tolerance = 1e-6
  )
  expect_identical(unclass(print(explained)), unclass(explained))
})

test_that("surv_shap() can explain survival at a time point", {
  skip_if_not_installed("kernelshap")
  skip_if_not_installed("shapviz")
  fixture <- shap_fixture()

  explained <- surv_shap(
    fixture$result, fixture$validation, model = "coxph", time_point = 36,
    n_explain = 15L, bg_n = 10L, seed = 2
  )
  expect_identical(explained$time_point, 36)
  expect_match(explained$quantity, "survival at t = 36")
  # Explaining S(t) is not explaining the risk score.
  risk <- surv_shap(fixture$result, fixture$validation, model = "coxph",
                    n_explain = 15L, bg_n = 10L, seed = 2)
  expect_false(isTRUE(all.equal(explained$values, risk$values)))
  expect_true(all(explained$values >= 0 & explained$values <= 1))
})

test_that("the survival SHAP plots all draw", {
  skip_if_not_installed("kernelshap")
  skip_if_not_installed("shapviz")
  fixture <- shap_fixture()
  explained <- surv_shap(fixture$result, fixture$validation, model = "coxph",
                         n_explain = 20L, bg_n = 10L, seed = 3)

  expect_s3_class(plot_surv_shap_importance(explained), "ggplot")
  expect_s3_class(plot_surv_shap_beeswarm(explained), "ggplot")
  expect_s3_class(plot_surv_shap_dependence(explained, "x1"), "ggplot")
  expect_s3_class(plot_surv_shap_waterfall(explained, row = 2), "ggplot")

  # A bare shapviz object is accepted too.
  expect_s3_class(plot_surv_shap_importance(explained$shap), "ggplot")
  expect_error(plot_surv_shap_importance("nope"), "shapviz")
})

test_that("surv_shap() validates what it is asked to explain", {
  skip_if_not_installed("kernelshap")
  skip_if_not_installed("shapviz")
  fixture <- shap_fixture()

  expect_error(surv_shap(fixture$result, fixture$validation, model = "nope"),
               "Unknown model")
  expect_error(surv_shap(fixture$result, "nope"), "data frame")
  expect_error(
    surv_shap(fixture$result, fixture$validation[, c("x1", "x2")]),
    "must contain"
  )
  expect_error(
    surv_shap(fixture$result, fixture$validation, n_explain = 1L),
    "n_explain"
  )
})

test_that("a risk-only model can only be explained on the risk score", {
  skip_if_not_installed("kernelshap")
  skip_if_not_installed("shapviz")
  skip_if_not_installed("rpart")
  set.seed(52)
  data <- sim_surv_data(n = 160, seed = 52)
  index <- sample.int(nrow(data))
  spec <- surv_spec(time = "time", status = "status", time_points = 12)
  spec <- surv_set_learners(spec, "rpart", tune = FALSE)
  spec <- surv_set_tuning(spec, folds = 2L, evals = 1L)
  result <- run_surv_workflow(
    spec,
    train = data[index[1:120], ],
    validation = data[index[121:160], ],
    ci_R = 0L
  )
  expect_error(
    surv_shap(result, data[index[121:160], ], model = "rpart", time_point = 12),
    "risk score"
  )
  explained <- surv_shap(result, data[index[121:160], ], model = "rpart",
                         n_explain = 10L, bg_n = 10L, seed = 4)
  expect_identical(explained$quantity, "risk score")
})
