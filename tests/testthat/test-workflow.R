workflow_data <- function(seed = 21L, n = 300L) {
  data <- sim_surv_data(n = n, seed = seed)
  set.seed(seed)
  index <- sample.int(nrow(data))
  list(
    data = data,
    train = data[index[1:200], , drop = FALSE],
    validation = data[index[201:250], , drop = FALSE],
    test = data[index[251:300], , drop = FALSE]
  )
}

workflow_spec <- function() {
  spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
  spec <- surv_set_learners(spec, "coxph", tune = FALSE)
  spec <- surv_set_tuning(spec, folds = 2L, evals = 1L)
  spec
}

test_that("run_surv_workflow() scores every model on every set", {
  split <- workflow_data()
  result <- run_surv_workflow(workflow_spec(), split$train, split$validation, split$test)

  expect_s3_class(result, "uniml_result")
  expect_identical(result$data_info$task_type, "survival")
  expect_setequal(unique(result$metrics$set), c("train", "validation", "test"))
  expect_identical(unique(result$metrics$model), "coxph")
  expect_true(all(c("cindex", "ibs", "brier_12", "auc_36") %in% names(result$metrics)))

  # Selection uses the validation set; the test set only reports.
  expect_identical(result$meta$select_on, "validation")
  expect_identical(unique(result$ranking$set), "validation")
  expect_identical(result$ranking$rank, 1L)
  expect_equal(
    result$ranking$score,
    result$metrics$cindex[result$metrics$set == "validation"]
  )

  # Calibration and decision curves are computed on the comparison set only.
  expect_identical(unique(result$calibration$metrics$set), "validation")
  expect_true(all(result$calibration$metrics$model == "coxph"))
  expect_true(all(c("model", "time", "threshold", "strategy", "net_benefit") %in%
                    names(result$dca)))

  expect_true(is.list(result$models))
  expect_true("coxph" %in% names(result$models))
  expect_s3_class(result$models$coxph$learner, "Learner")
  expect_true(all(is.na(result$meta$failures$stage)) == FALSE || nrow(result$meta$failures) == 0L)
})

test_that("the test and external sets never influence the fit", {
  split <- workflow_data()
  spec <- workflow_spec()
  base <- run_surv_workflow(spec, split$train, split$validation)

  # Scrambling the validation outcome must not change the model, only the
  # numbers reported on that set.
  scrambled_validation <- split$validation
  scrambled_validation$status <- 1L - scrambled_validation$status
  changed <- run_surv_workflow(spec, split$train, scrambled_validation)

  train_base <- base$metrics$cindex[base$metrics$set == "train"]
  train_changed <- changed$metrics$cindex[changed$metrics$set == "train"]
  expect_equal(train_base, train_changed, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    base$metrics$cindex[base$metrics$set == "validation"],
    changed$metrics$cindex[changed$metrics$set == "validation"]
  )))

  external <- split$test
  with_external <- run_surv_workflow(spec, split$train, split$validation,
                                    external = external)
  expect_true("external" %in% with_external$metrics$set)
  expect_equal(
    with_external$metrics$cindex[with_external$metrics$set == "train"],
    base$metrics$cindex[base$metrics$set == "train"],
    tolerance = 1e-12
  )
})

test_that("predict_surv() reproduces the stored validation predictions", {
  split <- workflow_data()
  result <- run_surv_workflow(workflow_spec(), split$train, split$validation)

  predicted <- predict_surv(result, split$validation)
  expect_s3_class(predicted, "data.frame")
  expect_true(all(c("risk", "S_12", "S_36") %in% names(predicted)))
  expect_identical(nrow(predicted), nrow(split$validation))
  expect_true(all(predicted$S_36 <= predicted$S_12))

  stored <- result$meta$predictions[["coxph:validation"]]
  expect_equal(predicted$risk, stored$risk, tolerance = 1e-10)
  expect_equal(predicted$S_36, stored$survival[, 2], tolerance = 1e-10)

  expect_error(predict_surv(result, split$validation, model = "nope"), "Unknown model")
  expect_error(
    predict_surv(result, split$validation[, c("x1", "x2", "sex")]),
    "must contain"
  )
})

test_that("a discrimination-only learner is reported, not silently dropped", {
  skip_if_not_installed("rpart")
  split <- workflow_data()
  spec <- workflow_spec()
  spec <- surv_set_learners(spec, c("coxph", "rpart"), tune = FALSE)
  result <- run_surv_workflow(spec, split$train, split$validation)

  expect_setequal(unique(result$metrics$model), c("coxph", "rpart"))
  rpart_rows <- result$metrics[result$metrics$model == "rpart", , drop = FALSE]
  expect_true(all(is.finite(rpart_rows$cindex)))
  expect_true(all(is.na(rpart_rows$ibs)))
  expect_true(any(grepl("risk score only", result$meta$notes)))

  # The ranking still orders both models by their concordance index.
  expect_identical(nrow(result$ranking), 2L)
  expect_error(
    predict_surv(result, split$validation, model = "rpart"),
    "risk score only"
  )
})

test_that("cross validation pools the out-of-fold predictions", {
  split <- workflow_data()
  spec <- workflow_spec()
  prepared <- as_surv_data(split$train, spec)
  learners <- resolve_surv_learners("coxph", tune = FALSE)

  cv <- surv_cross_validate(prepared, learners, spec, folds = 3L)
  expect_identical(nrow(cv), 1L)
  expect_true(all(c("model", "cindex", "ibs", "brier_12") %in% names(cv)))
  expect_true(is.finite(cv$cindex))
})

test_that("the learner registry is honest about capabilities", {
  registry <- surv_learner_registry()
  expect_true(all(c("key", "label", "package", "available", "tunable") %in%
                    names(registry)))
  expect_true("coxph" %in% registry$key)
  expect_true(all(surv_learner_registry(available_only = TRUE)$available))

  learners <- resolve_surv_learners(c("coxph", "rpart"), tune = FALSE)
  expect_identical(names(learners), c("coxph", "rpart"))
  expect_true(learners$coxph$distr)
  expect_false(learners$rpart$distr)
  expect_error(resolve_surv_learners("nope"), "Unknown survival learner")
  expect_error(resolve_surv_learners(list(a = 1)), "must be an mlr3 learner")

  # A tuning search space exists only where it means something.
  tuned <- resolve_surv_learners("rpart", tune = TRUE)
  expect_true(tuned$rpart$tune)
  expect_false(resolve_surv_learners("kaplan", tune = TRUE)$kaplan$tune)
})

test_that("workflow inputs are validated", {
  split <- workflow_data()
  spec <- workflow_spec()
  expect_error(
    run_surv_workflow(uniml::uniml_spec("binary", "y", "yes"), split$train, split$validation),
    "task_type"
  )
  expect_error(
    run_surv_workflow(spec, split$train[, c("x1", "x2")], split$validation),
    "does not contain"
  )
})
