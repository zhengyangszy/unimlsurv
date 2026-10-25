# L3 task layer --------------------------------------------------------------

#' Build an mlr3 survival task
#'
#' @param data A data frame, a `unimlsurv_data` object, or an mlr3 `TaskSurv`.
#' @param spec A specification created by [surv_spec()].
#'
#' @return An `mlr3proba::TaskSurv` with the id `"unimlsurv_task"`.
#' @references
#' Sonabend R, Kirchner P, Bender A, et al. (2021). mlr3proba: an R package for
#' machine learning in survival analysis. *Bioinformatics*, 37(17), 2789-2791.
#' doi:10.1093/bioinformatics/btab039
#' @export
#' @examples
#' data <- sim_surv_data(n = 60, seed = 1)
#' spec <- surv_spec(time = "time", status = "status")
#' task <- make_surv_task(data, spec)
#' task$nrow
make_surv_task <- function(data, spec) {
  validate_surv_spec(spec)
  block <- spec$survival
  if (inherits(data, "TaskSurv")) {
    task <- data$clone(deep = TRUE)
    if (!identical(task$target_names, c(block$time, block$status))) {
      .abort("`spec` does not match the time and status columns of `data`.")
    }
    return(task)
  }
  prepared <- if (inherits(data, "unimlsurv_data")) {
    validate_surv_data(data)
    data
  } else {
    as_surv_data(data, spec)
  }
  frame <- prepared$data[, c(prepared$features, block$time, block$status),
                         drop = FALSE]
  task <- mlr3proba::as_task_surv(
    frame,
    time = block$time,
    event = block$status,
    id = "unimlsurv_task"
  )
  task
}

#' Fit and predict a survival learner
#'
#' @param task An `mlr3proba::TaskSurv`.
#' @param learner An `mlr3` survival learner.
#' @param params Optional named list of hyperparameters applied before training.
#'
#' @return The trained `Learner`.
#' @references
#' Bender A, Groll A, Scheipl F (2018). A general machine learning framework for
#' survival analysis. arXiv:1806.06108.
#' @export
#' @examples
#' data <- sim_surv_data(n = 80, seed = 3)
#' spec <- surv_spec(time = "time", status = "status")
#' task <- make_surv_task(data, spec)
#' trained <- train_surv(task, mlr3::lrn("surv.coxph"))
#' trained$predict(task)$predict_types
train_surv <- function(task, learner, params = NULL) {
  if (!inherits(task, "TaskSurv")) {
    .abort("`task` must be an mlr3 survival task.")
  }
  if (!inherits(learner, "Learner")) {
    .abort("`learner` must be an mlr3 learner.")
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
  learner
}

#' Extract survival probabilities from a prediction
#'
#' `mlr3proba` returns a distribution object. Everything downstream of this
#' function works on a plain numeric matrix with one row per patient and one
#' column per time point, which is also what gets exported to disk.
#'
#' @param prediction An `mlr3proba::PredictionSurv`.
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
#' prediction <- mlr3::lrn("surv.coxph")$train(task)$predict(task)
#' head(survival_matrix(prediction, c(12, 36)))
survival_matrix <- function(prediction, time_points) {
  if (!inherits(prediction, "PredictionSurv")) {
    .abort("`prediction` must be an mlr3 survival prediction.")
  }
  .assert_numeric_vector(time_points, "time_points", min = 0)
  if (is.null(prediction$distr)) {
    .abort(paste(
      "The prediction carries no survival distribution.",
      "Set `predict_type = \"distr\"` on the learner."
    ))
  }
  values <- prediction$distr$survival(time_points)
  # `distr$survival()` returns time points in rows; everything downstream of
  # here uses one row per patient.
  matrix_values <- t(as.matrix(values))
  colnames(matrix_values) <- paste0("t", format(time_points, trim = TRUE))
  rownames(matrix_values) <- NULL
  matrix_values
}

#' Extract the risk score from a prediction
#'
#' The risk score is the quantity a concordance index ranks: a higher value means
#' a shorter expected survival. `mlr3proba` exposes it as `crank`, and falls back
#' to the linear predictor when `crank` is absent.
#'
#' @param prediction An `mlr3proba::PredictionSurv`.
#'
#' @return A numeric vector.
#' @export
#' @examples
#' data <- sim_surv_data(n = 60, seed = 1)
#' spec <- surv_spec(time = "time", status = "status")
#' task <- make_surv_task(data, spec)
#' head(survival_risk(mlr3::lrn("surv.coxph")$train(task)$predict(task)))
survival_risk <- function(prediction) {
  if (!inherits(prediction, "PredictionSurv")) {
    .abort("`prediction` must be an mlr3 survival prediction.")
  }
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
