test_that("a survival workflow never runs a feature selector", {
  set.seed(61)
  data <- sim_surv_data(n = 200, seed = 61)
  index <- sample.int(nrow(data))
  spec <- surv_spec(time = "time", status = "status", time_points = 12)
  spec <- surv_set_learners(spec, "coxph", tune = FALSE)
  spec <- surv_set_tuning(spec, folds = 2L, evals = 1L)

  # `surv_spec()` asks for no selection, because uniml's selectors screen the
  # outcome column, which for a survival model is the event indicator. How that
  # is spelled depends on the installed uniml, so the test asks the package.
  expect_identical(spec$selection$methods, .surv_no_selection())

  # Even a specification that asks for one cannot smuggle it into the workflow:
  # the run overrides it, so the survival path is deterministic.
  spec$selection$methods <- "none"
  result <- run_surv_workflow(
    spec,
    train = data[index[1:150], ],
    validation = data[index[151:200], ],
    ci_R = 0L
  )
  expect_true(is.finite(result$metrics$cindex[result$metrics$set == "validation"]))
  expect_identical(result$meta$spec$selection$methods, .surv_no_selection())
})
