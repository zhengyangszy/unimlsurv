# L3 task layer --------------------------------------------------------------
#
# From 0.5.0 the task and prediction objects are native (no mlr3proba). The
# public names and the value shapes are preserved: `make_surv_task()` returns a
# `unimlsurv_task`, `train_surv()` returns a fitted learner, `survival_matrix()`
# returns an n x time_points matrix named `t<point>`, and `survival_risk()`
# returns the score the concordance index ranks (higher = shorter survival).
# Foreign mlr3 learners keep working through an escape hatch that requires
# `mlr3proba` (Suggests).

#' Build a native survival task
#'
#' A native task is a lightweight list with class `unimlsurv_task`: the data
#' frame (features plus the time and status columns), the column names, and the
#' row/event counts. The workflow engine treats tasks opaquely and hands them
#' back to the adapter's train and predict hooks, so nothing outside this
#' package needs to know the internals.
#'
#' @param data A data frame, a `unimlsurv_data` object, or an mlr3 `TaskSurv`.
#' @param spec A specification created by [surv_spec()].
#'
#' @return A native task with class `unimlsurv_task`.
#' @references
#' Sonabend R, Kirchner P, Bender A, et al. (2021). mlr3proba: an R package for
#' machine learning in survival analysis. *Bioinformatics*, 37(17), 2789-2791.
#' doi:10.1093/bioinformatics/btab039
#' @export
#' @examples
#' data <- sim_surv_data(n = 60, seed = 1)
#' spec <- surv_spec(time = "time", status = "status")
#' task <- make_surv_task(data, spec)
#' task$n
make_surv_task <- function(data, spec) {
  validate_surv_spec(spec)
  block <- spec$survival
  if (inherits(data, "unimlsurv_task")) {
    if (!identical(data$time_col, block$time) ||
        !identical(data$status_col, block$status)) {
      .abort("`spec` does not match the time and status columns of `task`.")
    }
    return(data)
  }
  if (inherits(data, "TaskSurv")) {
    # An mlr3 task is still accepted: convert it to the native frame.
    frame <- as.data.frame(data)
    if (!all(c(block$time, block$status) %in% names(frame))) {
      .abort("`spec` does not match the time and status columns of `data`.")
    }
    data <- frame
  }
  prepared <- if (inherits(data, "unimlsurv_data")) {
    validate_surv_data(data)
    data
  } else {
    as_surv_data(data, spec)
  }
  frame <- prepared$data[, c(prepared$features, block$time, block$status),
                         drop = FALSE]
  status_values <- as.numeric(frame[[block$status]])
  structure(
    list(
      data = frame,
      features = prepared$features,
      time_col = block$time,
      status_col = block$status,
      n = nrow(frame),
      n_events = sum(status_values == 1L)
    ),
    class = "unimlsurv_task"
  )
}

#' Print a native survival task
#'
#' @param x A `unimlsurv_task`.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.unimlsurv_task <- function(x, ...) {
  cat(sprintf(
    paste0(
      "<unimlsurv_task> rows: %d; features: %d; events: %d; time: %s; ",
      "status: %s\n"
    ),
    x$n, length(x$features), x$n_events, x$time_col, x$status_col
  ))
  invisible(x)
}

# A native task can still be handed to foreign mlr3 learners through the
# escape hatch; that path requires mlr3proba (Suggests).
.task_surv_to_mlr3 <- function(task) {
  if (!requireNamespace("mlr3proba", quietly = TRUE)) {
    .abort(paste0(
      "Foreign mlr3 survival learners need the mlr3proba package, which is ",
      "not on CRAN. Install it from GitHub (mlr-org/mlr3proba), or use the ",
      "native learners (coxph, kaplan, rpart)."
    ))
  }
  mlr3proba::as_task_surv(
    task$data,
    time = task$time_col,
    event = task$status_col,
    id = "unimlsurv_task"
  )
}

#' Fit and predict a survival learner
#'
#' Native learners fit `survival::coxph`, a Kaplan-Meier baseline, or
#' `rpart(method = "exp")` directly. Foreign mlr3 learners (registered through
#' [surv_register_learner()]) are trained through the mlr3 stack, which needs
#' `mlr3proba`.
#'
#' @param task A native task from [make_surv_task()].
#' @param learner A native learner from the registry, or a foreign mlr3 learner.
#' @param params Optional named list of hyperparameters applied before training.
#'
#' @return A fitted native learner (class `unimlsurv_fitted`) or the trained
#'   mlr3 learner.
#' @references
#' Bender A, Groll A, Scheipl F (2018). A general machine learning framework for
#' survival analysis. arXiv:1806.06108.
#' @export
#' @examples
#' data <- sim_surv_data(n = 80, seed = 3)
#' spec <- surv_spec(time = "time", status = "status")
#' task <- make_surv_task(data, spec)
#' trained <- train_surv(task, resolve_surv_learners("coxph", tune = FALSE)$coxph$learner)
#' trained$engine
train_surv <- function(task, learner, params = NULL) {
  if (inherits(learner, "unimlsurv_native_learner")) {
    if (!inherits(task, "unimlsurv_task")) {
      .abort("`task` must be a native task from make_surv_task().")
    }
    return(.train_native(task, learner, params))
  }
  if (inherits(learner, "Learner")) {
    # Foreign mlr3 learners keep the mlr3 path; mlr3proba is the price of the
    # escape hatch and lives in Suggests.
    if (!inherits(task, "TaskSurv")) {
      task <- .task_surv_to_mlr3(task)
    }
    if (!is.null(params) && length(params) > 0L) {
      if (!is.list(params) || is.null(names(params))) {
        .abort("`params` must be a named list.")
      }
      learner$param_set$values <- utils::modifyList(
        learner$param_set$values,
        params
      )
    }
    learner$train(task)
    return(learner)
  }
  .abort(paste(
    "`learner` must be a native learner from surv_learner_registry() or a",
    "foreign mlr3 survival learner."
  ))
}

#' Extract survival probabilities from a prediction
#'
#' @param prediction A native prediction from [predict_surv_learner()], or an
#'   mlr3 `PredictionSurv` (foreign learners).
#' @param time_points Numeric vector of time points.
#'
#' @return A matrix with `length(time_points)` columns, named after the time
#'   points.
#' @references
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). Assessment and comparison
#' of prognostic classification schemes for survival data. *Statistics in
#' Medicine*, 18(17-18), 2529-2545.
#' @export
#' @examples
#' data <- sim_surv_data(n = 60, seed = 1)
#' spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
#' task <- make_surv_task(data, spec)
#' fitted <- train_surv(task, resolve_surv_learners("coxph", tune = FALSE)$coxph$learner)
#' head(survival_matrix(predict_surv_learner(fitted, data, c(12, 36)), c(12, 36)))
survival_matrix <- function(prediction, time_points) {
  .assert_numeric_vector(time_points, "time_points", min = 0)
  if (inherits(prediction, "unimlsurv_prediction")) {
    return(prediction$survival)
  }
  if (inherits(prediction, "PredictionSurv")) {
    if (!requireNamespace("mlr3proba", quietly = TRUE)) {
      .abort(paste(
        "mlr3 predictions need the mlr3proba package, which is not on CRAN.",
        "Use the native learners, or install mlr3proba from GitHub."
      ))
    }
    if (is.null(prediction$distr)) {
      .abort(paste(
        "The prediction carries no survival distribution.",
        "Set `predict_type = \"distr\"` on the learner."
      ))
    }
    values <- prediction$distr$survival(time_points)
    matrix_values <- t(as.matrix(values))
    colnames(matrix_values) <- paste0("t", format(time_points, trim = TRUE))
    rownames(matrix_values) <- NULL
    return(matrix_values)
  }
  .abort("`prediction` must be a native or mlr3 survival prediction.")
}

#' Extract the risk score from a prediction
#'
#' The risk score is the quantity a concordance index ranks: a higher value means
#' a shorter expected survival.
#'
#' @param prediction A native prediction from [predict_surv_learner()], or an
#'   mlr3 `PredictionSurv` (foreign learners).
#'
#' @return A numeric vector.
#' @export
#' @examples
#' data <- sim_surv_data(n = 60, seed = 1)
#' spec <- surv_spec(time = "time", status = "status")
#' task <- make_surv_task(data, spec)
#' fitted <- train_surv(task, resolve_surv_learners("coxph", tune = FALSE)$coxph$learner)
#' head(survival_risk(predict_surv_learner(fitted, data, c(12, 36, 60))))
survival_risk <- function(prediction) {
  if (inherits(prediction, "unimlsurv_prediction")) {
    return(as.numeric(prediction$risk))
  }
  if (inherits(prediction, "PredictionSurv")) {
    if (!is.null(prediction$crank)) {
      return(as.numeric(prediction$crank))
    }
    if (!is.null(prediction$lp)) {
      return(as.numeric(prediction$lp))
    }
    .abort(paste(
      "The prediction has neither `crank` nor `lp`;",
      "a survival learner must provide a risk score."
    ))
  }
  .abort("`prediction` must be a native or mlr3 survival prediction.")
}
