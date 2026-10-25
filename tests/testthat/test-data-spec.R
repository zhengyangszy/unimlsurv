test_that("surv_spec() records the outcome columns and time points", {
  spec <- surv_spec(time = "time", status = "status", time_points = c(36, 12, 36))
  expect_s3_class(spec, "unimlsurv_spec")
  expect_s3_class(spec, "uniml_spec")
  expect_identical(spec$task_type, "survival")
  expect_identical(spec$outcome, "status")
  expect_null(spec$positive)
  expect_identical(spec$survival$time, "time")
  expect_identical(spec$survival$status, "status")
  # Time points are de-duplicated and sorted, so a table and a figure agree.
  expect_identical(spec$survival$time_points, c(12, 36))
  expect_true(validate_surv_spec(spec))
  expect_identical(unclass(print(spec)), unclass(spec))

  expect_error(surv_spec(time = "y", status = "y"), "must be different")
  expect_error(surv_spec(time_points = numeric()), "time_points")
  expect_error(surv_spec(censoring_model = "nope"), "censoring_model")
  expect_error(surv_spec(seed = 0), "seed")
})

test_that("survival setters keep the survival block", {
  spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
  updated <- surv_set_time_points(spec, c(3, 6))
  expect_identical(updated$survival$time_points, c(3, 6))
  expect_identical(surv_set_censoring(spec, "none")$survival$censoring_model, "none")

  trained <- surv_set_learners(spec, "coxph", tune = FALSE)
  expect_true("coxph" %in% names(trained$survival$learners))
  expect_false(trained$survival$learners$coxph$tune)

  # Preprocessing and tuning are the uniml setters, so the specification keeps
  # one shape across the two packages.
  scaled <- uniml::set_preprocess(spec, scale = "none", impute = "none")
  expect_identical(scaled$preprocess$scale, "none")
  expect_identical(scaled$survival$time_points, c(12, 36))

  tuned <- surv_set_tuning(spec, folds = 3L, evals = 5L)
  expect_identical(tuned$tuning$folds, 3L)
  expect_identical(tuned$survival$tuning_metric, "cindex")
  expect_error(surv_set_tuning(spec, metric = "brier"), "concordance index")
})

test_that("as_surv_data() checks the outcome and keeps the features", {
  data <- sim_surv_data(n = 120, seed = 3)
  spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
  prepared <- as_surv_data(data, spec)

  expect_s3_class(prepared, "unimlsurv_data")
  expect_identical(prepared$meta$n, 120L)
  expect_identical(prepared$features, c("x1", "x2", "sex"))
  expect_identical(prepared$meta$n_events, sum(data$status))
  expect_true(validate_surv_data(prepared))
  expect_identical(unclass(print(prepared)), unclass(prepared))

  # Missing outcomes are dropped, not imputed.
  broken <- data
  broken$time[1] <- NA
  broken$status[2] <- NA
  expect_identical(as_surv_data(broken, spec)$meta$n, 118L)

  # A non-positive time is not a follow-up time.
  broken <- data
  broken$time[3] <- 0
  expect_error(as_surv_data(broken, spec), "strictly positive")

  # Character features are refused rather than silently encoded.
  broken <- data
  broken$label <- as.character(seq_len(nrow(data)))
  expect_error(as_surv_data(broken, spec), "Unsupported feature type")

  # Logical features become factors.
  broken <- data
  broken$flag <- broken$sex == "male"
  expect_identical(
    class(as_surv_data(broken, spec)$data$flag),
    "factor"
  )

  # A validation set with one status value is still a valid set to predict.
  only_censored <- data[data$status == 0L, , drop = FALSE]
  expect_error(as_surv_data(only_censored, spec), "single status value")
  expect_s3_class(
    as_surv_data(only_censored, spec, "validation", require_both = FALSE),
    "unimlsurv_data"
  )

  # The event value can be something other than 1.
  recoded <- data
  recoded$status <- ifelse(data$status == 1L, "dead", "alive")
  expect_error(as_surv_data(recoded, spec), "logical or numeric")
})

test_that("make_surv_task() builds an mlr3 survival task", {
  data <- sim_surv_data(n = 100, seed = 4)
  spec <- surv_spec(time = "time", status = "status")
  task <- make_surv_task(data, spec)

  expect_s3_class(task, "TaskSurv")
  expect_identical(task$nrow, 100L)
  expect_identical(task$target_names, c("time", "status"))
  expect_setequal(task$feature_names, c("x1", "x2", "sex"))
  expect_error(
    make_surv_task(data, surv_spec(time = "nope", status = "status")),
    "does not contain"
  )

  # A task round-trips.
  expect_s3_class(make_surv_task(task, spec), "TaskSurv")
})

test_that("survival_matrix() and survival_risk() expose the prediction", {
  data <- sim_surv_data(n = 100, seed = 5)
  spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
  task <- make_surv_task(data, spec)
  prediction <- mlr3::lrn("surv.coxph")$train(task)$predict(task)

  survival <- survival_matrix(prediction, c(12, 36))
  expect_identical(dim(survival), c(100L, 2L))
  expect_identical(colnames(survival), c("t12", "t36"))
  expect_true(all(survival >= 0 & survival <= 1))
  # Survival cannot increase with time.
  expect_true(all(survival[, 2] <= survival[, 1]))

  risk <- survival_risk(prediction)
  expect_identical(length(risk), 100L)
  expect_true(all(is.finite(risk)))

  expect_error(survival_matrix(prediction, "time"), "time_points")
  expect_error(survival_matrix(prediction, numeric()), "time_points")
})
