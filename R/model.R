# L7 tuning and cross-validation ---------------------------------------------

.surv_measure <- function(metric = "cindex") {
  if (!identical(metric, "cindex")) {
    .abort("Only `\"cindex\"` can be used as a survival tuning measure.")
  }
  mlr3::msr("surv.cindex")
}

#' Tune a survival learner
#'
#' Random search (or any `mlr3tuning` tuner) over the learner's search space,
#' scored by the concordance index. Learners without a search space are returned
#' untuned rather than being given one that would not mean anything.
#'
#' @param task An `mlr3proba::TaskSurv`.
#' @param learner_spec An entry from [resolve_surv_learners()].
#' @param spec A survival specification.
#'
#' @return A list with `learner`, `best_params`, `best_metric`, `status`, and
#'   `measure`.
#' @references
#' Probst P, Boulesteix AL, Bischl B (2019). Tunability: importance of
#' hyperparameters. *Journal of Machine Learning Research*, 20(53), 1-32.
#'
#' Bischl B, Binder M, Lang M, et al. (2021). Hyperparameter optimization: a
#' brief introduction. arXiv:2107.05847.
#' @export
#' @examples
#' \dontrun{
#' data <- sim_surv_data(n = 150, seed = 1)
#' spec <- surv_spec(time = "time", status = "status")
#' task <- make_surv_task(data, spec)
#' learners <- resolve_surv_learners("rpart", tune = TRUE)
#' tuned <- tune_surv(task, learners$rpart, spec)
#' tuned$status
#' }
tune_surv <- function(task, learner_spec, spec) {
  if (!inherits(task, "TaskSurv")) {
    .abort("`task` must be an mlr3 survival task.")
  }
  validate_surv_spec(spec)
  learner <- learner_spec$learner
  space <- learner_spec$search_space
  untuned <- list(
    learner = learner$clone(deep = TRUE),
    best_params = list(),
    best_metric = NA_real_,
    status = "untuned",
    measure = NULL
  )
  if (!isTRUE(learner_spec$tune) || is.null(space)) {
    return(untuned)
  }
  measure <- .surv_measure(spec$survival$tuning_metric %||% "cindex")
  folds <- max(2L, as.integer(spec$tuning$folds))
  repeats <- max(1L, as.integer(spec$tuning$repeats))
  resampling <- if (repeats > 1L) {
    mlr3::rsmp("repeated_cv", folds = folds, repeats = repeats)
  } else {
    mlr3::rsmp("cv", folds = folds)
  }
  tuner <- tryCatch(
    mlr3tuning::tnr(spec$tuning$tuner),
    error = function(e) NULL
  )
  if (is.null(tuner)) {
    .abort(sprintf("Unknown tuning strategy \"%s\".", spec$tuning$tuner))
  }
  terminator <- if (identical(spec$tuning$terminator, "seconds")) {
    mlr3tuning::trm("run_time", secs = as.numeric(spec$tuning$seconds))
  } else {
    mlr3tuning::trm("evals", n_evals = as.integer(spec$tuning$evals))
  }
  tunable <- mlr3tuning::auto_tuner(
    tuner = tuner,
    learner = learner$clone(deep = TRUE),
    resampling = resampling,
    measure = measure,
    terminator = terminator,
    search_space = space
  )
  set.seed(as.integer(spec$seed))
  tunable$train(task)
  instance <- tunable$tuning_instance
  result <- instance$result
  if (is.null(result) || nrow(result) == 0L) {
    return(list(
      learner = learner$clone(deep = TRUE),
      best_params = list(),
      best_metric = NA_real_,
      status = "failed",
      measure = measure
    ))
  }
  best_params <- result$learner_param_vals[[1L]]
  best_metric <- as.numeric(result[[measure$id]][[1L]])
  list(
    learner = learner$clone(deep = TRUE),
    best_params = best_params,
    best_metric = best_metric,
    status = "tuned",
    measure = measure
  )
}

#' Cross-validate a survival model
#'
#' Repeated cross-validation scored on the out-of-fold predictions, which is the
#' only honest way to report a C-index from resampling: averaging fold-wise
#' C-indices weights folds by their number of comparable pairs, which nobody
#' wants.
#'
#' @param data A `unimlsurv_data` object.
#' @param learners A resolved learner list from [resolve_surv_learners()].
#' @param spec A survival specification.
#' @param folds Number of folds.
#' @param repeats Number of repetitions.
#'
#' @return A data frame with one row per model, repetition, and fold, plus the
#'   `cindex`, `ibs`, and per-time `brier_<t>`/`auc_<t>` columns.
#' @references
#' Varma S, Simon R (2006). Bias in error estimation when using cross-validation
#' for model selection. *BMC Bioinformatics*, 7, 91.
#' doi:10.1186/1471-2105-7-91
#'
#' Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
#' Development, Validation, and Updating*, 2nd edition. Springer.
#' doi:10.1007/978-3-030-16399-0
#' @export
#' @examples
#' \dontrun{
#' data <- as_surv_data(sim_surv_data(n = 200, seed = 1),
#'                      surv_spec(time = "time", status = "status"))
#' spec <- surv_spec(time = "time", status = "status")
#' learners <- resolve_surv_learners("coxph", tune = FALSE)
#' surv_cross_validate(data, learners, spec, folds = 3)
#' }
surv_cross_validate <- function(data, learners, spec, folds = 5L, repeats = 1L) {
  validate_surv_data(data)
  validate_surv_spec(spec)
  .assert_scalar_number(folds, "folds", integer = TRUE, min = 2)
  .assert_scalar_number(repeats, "repeats", integer = TRUE, min = 1)
  task <- make_surv_task(data, spec)
  time_points <- spec$survival$time_points
  resampling <- if (repeats > 1L) {
    mlr3::rsmp("repeated_cv", folds = as.integer(folds), repeats = as.integer(repeats))
  } else {
    mlr3::rsmp("cv", folds = as.integer(folds))
  }

  rows <- list()
  for (key in names(learners)) {
    entry <- learners[[key]]
    learner <- entry$learner$clone(deep = TRUE)
    if (!is.null(entry$search_space) && isTRUE(entry$tune)) {
      tuned <- tryCatch(
        tune_surv(task, entry, spec),
        error = function(e) e
      )
      if (inherits(tuned, "error")) {
        next
      }
      if (length(tuned$best_params) > 0L) {
        learner$param_set$values <- utils::modifyList(
          learner$param_set$values,
          tuned$best_params
        )
      }
    }
    set.seed(as.integer(spec$seed))
    result <- tryCatch(
      mlr3::resample(task, learner, resampling, store_models = FALSE),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      next
    }
    prediction <- result$prediction()
    survival <- survival_matrix(prediction, time_points)
    truth <- task$truth()
    metrics <- calc_surv_metrics(
      time = as.numeric(truth[, 1L]),
      status = as.integer(truth[, 2L]),
      risk = survival_risk(prediction),
      survival = survival,
      time_points = time_points,
      censoring = spec$survival$censoring_model
    )
    row <- data.frame(
      model = key,
      repetition = 1L,
      fold = NA_integer_,
      metrics$row,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- row
  }
  if (length(rows) == 0L) {
    return(data.frame(
      model = character(), repetition = integer(), fold = integer(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
