# SHAP for survival predictions ----------------------------------------------
#
# A survival learner predicts a risk score and, when it can, a whole survival
# curve. There is no single "the prediction" to attribute, so the quantity being
# explained has to be named:
#
# * `time_point = NULL` explains the risk score (`crank`), which is what the
#   concordance index ranks. This is the quantity the published survival SHAP
#   examples explain.
# * `time_point = t` explains the predicted survival probability at that time,
#   S(t). The same patient can have a feature that raises the risk score and
#   lowers S(365), so the two are not interchangeable.
#
# Kernel SHAP is used because survival learners are not tree ensembles in
# general; a TreeSHAP shortcut would only cover some of them and would silently
# change the estimator for the rest.

#' SHAP values for a survival model
#'
#' Explains one model of a survival workflow with Kernel SHAP, either on the risk
#' score or on the predicted survival probability at one time point.
#'
#' `newdata` must contain the feature columns the model was trained on. The two
#' outcome columns are used only to build the mlr3 task that scoring needs; they
#' do not enter the explanation.
#'
#' @param result A `uniml_result` returned by [run_surv_workflow()].
#' @param newdata Data frame of patients to explain.
#' @param model Model key, or `"best"` for the top-ranked model.
#' @param time_point Time point at which to explain S(t), or `NULL` to explain the
#'   risk score.
#' @param n_explain Number of rows to explain; larger values are subsampled.
#' @param bg_n Number of background rows used as the reference distribution.
#' @param seed Random seed.
#' @param ... Further arguments passed to [kernelshap::kernelshap()].
#'
#' @return An object of class `unimlsurv_shap` with the `shapviz` object, the
#'   explained quantity, the time point, the model key, and the explained values.
#' @references
#' Lundberg SM, Lee SI (2017). A unified approach to interpreting model
#' predictions. *Advances in Neural Information Processing Systems*, 30.
#'
#' Lundberg SM, Erion G, Chen H, et al. (2020). From local explanations to global
#' understanding with explainable AI for trees. *Nature Machine Intelligence*,
#' 2(1), 56-67. doi:10.1038/s42256-019-0138-9
#'
#' Covert I, Lee SI (2021). Improving KernelSHAP: practical Shapley value
#' estimation using linear regression. *Proceedings of AISTATS*, PMLR 130,
#' 3457-3465.
#'
#' Sonabend R, Kirchner P, Bender A, et al. (2021). mlr3proba: an R package for
#' machine learning in survival analysis. *Bioinformatics*, 37(17), 2789-2791.
#' doi:10.1093/bioinformatics/btab039
#' @export
#' @examples
#' \dontrun{
#' result <- run_surv_workflow(spec, train, validation)
#' explained <- surv_shap(result, newdata = validation, model = "coxph")
#' plot_surv_shap_beeswarm(explained)
#' explained_sixty <- surv_shap(result, validation, model = "coxph", time_point = 60)
#' plot_surv_shap_dependence(explained_sixty, "x1")
#' }
surv_shap <- function(result, newdata, model = "best", time_point = NULL,
                      n_explain = 100L, bg_n = 50L, seed = 2026L, ...) {
  if (!requireNamespace("kernelshap", quietly = TRUE) ||
      !requireNamespace("shapviz", quietly = TRUE)) {
    .abort("Packages \"kernelshap\" and \"shapviz\" are required for SHAP.")
  }
  if (!inherits(result, "uniml_result")) {
    .abort("`result` must be a `uniml_result` from `run_surv_workflow()`.")
  }
  spec <- result$meta$spec
  validate_surv_spec(spec)
  block <- spec$survival
  models <- result$models
  if (is.null(models) || length(models) == 0L) {
    .abort("The result contains no fitted models.")
  }
  if (identical(model, "best")) {
    ranking <- result$ranking
    if (is.null(ranking) || nrow(ranking) == 0L) {
      .abort("No ranking is available; supply `model` explicitly.")
    }
    model <- as.character(ranking$model[[1L]])
  }
  .assert_scalar_character(model, "model")
  if (!model %in% names(models)) {
    .abort(sprintf(
      "Unknown model \"%s\"; available: %s.",
      model, paste(names(models), collapse = ", ")
    ))
  }
  entry <- models[[model]]
  if (!"distr" %in% entry$learner$predict_types && !is.null(time_point)) {
    .abort(sprintf(
      "Model \"%s\" cannot predict survival probabilities, so only the risk score can be explained; leave `time_point` NULL.",
      model
    ))
  }
  .assert_scalar_number(n_explain, "n_explain", integer = TRUE, min = 2)
  .assert_scalar_number(bg_n, "bg_n", integer = TRUE, min = 2)
  .assert_scalar_number(seed, "seed", integer = TRUE, min = 1)
  if (!is.null(time_point)) {
    .assert_scalar_number(time_point, "time_point", min = 0)
    if (!time_point %in% block$time_points) {
      warning(
        sprintf(
          "`time_point` = %s is not one of the reporting time points (%s).",
          format(time_point),
          paste(format(block$time_points), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  if (!is.data.frame(newdata)) {
    .abort("`newdata` must be a data frame.")
  }
  if (!all(c(block$time, block$status) %in% names(newdata))) {
    .abort(sprintf(
      "`newdata` must contain \"%s\" and \"%s\".",
      block$time, block$status
    ))
  }

  prepared <- .surv_shap_features(newdata, result, entry, spec)
  frame <- prepared$frame
  time_values <- prepared$time
  status_values <- prepared$status
  set.seed(as.integer(seed))
  if (nrow(frame) > n_explain) {
    explain_index <- sort(sample.int(nrow(frame), as.integer(n_explain)))
  } else {
    explain_index <- seq_len(nrow(frame))
  }
  X <- frame[explain_index, , drop = FALSE]
  explain_time <- time_values[explain_index]
  explain_status <- status_values[explain_index]

  prediction_function <- function(object, X) {
    features <- as.data.frame(X, check.names = FALSE)
    # Kernel SHAP calls the prediction function on the background rows too, so the
    # outcome columns cannot be the explained patients' own. A survival learner's
    # prediction is a function of the features only, so a placeholder of the right
    # length is enough; the outcome never enters the explanation.
    features[[block$time]] <- rep(stats::median(time_values), nrow(features))
    features[[block$status]] <- rep(1L, nrow(features))
    task <- mlr3proba::as_task_surv(
      features,
      time = block$time,
      event = block$status,
      id = "unimlsurv_shap"
    )
    prediction <- object$predict(task)
    if (is.null(time_point)) {
      return(as.numeric(survival_risk(prediction)))
    }
    as.numeric(survival_matrix(prediction, time_point)[, 1L])
  }

  bg_index <- if (nrow(frame) > bg_n) {
    sort(sample.int(nrow(frame), as.integer(bg_n)))
  } else {
    seq_len(nrow(frame))
  }
  settings <- list(
    object = entry$learner,
    X = X,
    bg_X = frame[bg_index, , drop = FALSE],
    pred_fun = prediction_function,
    bg_n = as.integer(bg_n),
    verbose = FALSE,
    seed = as.integer(seed)
  )
  overrides <- list(...)
  if (length(overrides) > 0L) {
    settings[names(overrides)] <- overrides
  }
  explanation <- do.call(kernelshap::kernelshap, settings)
  values <- as.numeric(prediction_function(entry$learner, X))
  structure(
    list(
      shap = shapviz::shapviz(explanation),
      quantity = if (is.null(time_point)) "risk score" else
        sprintf("survival at t = %s", format(time_point, trim = TRUE)),
      time_point = time_point,
      model = model,
      feature_names = colnames(X),
      values = values,
      data = X
    ),
    class = c("unimlsurv_shap", "uniml_shap")
  )
}

# The learner was trained on the preprocessed feature columns, so the explanation
# has to happen in exactly that space: preprocess with the stored preprocessor,
# then keep the trained feature names in the trained order.
.surv_shap_features <- function(newdata, result, entry, spec) {
  block <- spec$survival
  if (is.null(entry$preprocessor)) {
    .abort("The fitted model carries no preprocessor, so new data cannot be prepared.")
  }
  transformed <- uniml::apply_preprocessor(
    entry$preprocessor,
    newdata[, unique(c(
      setdiff(names(newdata), c(block$time, block$status, spec$id, spec$exclude)),
      block$time, block$status
    )), drop = FALSE]
  )
  required <- entry$task$feature_names
  missing_columns <- setdiff(required, names(transformed))
  if (length(missing_columns) > 0L) {
    .abort(sprintf(
      "The preprocessor did not produce the trained column(s): %s.",
      paste(missing_columns, collapse = ", ")
    ))
  }
  list(
    frame = as.data.frame(transformed[, required, drop = FALSE]),
    time = newdata[[block$time]],
    status = as.integer(newdata[[block$status]])
  )
}

#' Print a survival SHAP explanation
#'
#' @param x A `unimlsurv_shap` object.
#' @param ... Additional arguments, unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.unimlsurv_shap <- function(x, ...) {
  cat("<unimlsurv_shap>\n")
  cat(sprintf("model: %s\n", x$model))
  cat(sprintf("explained: %s\n", x$quantity))
  cat(sprintf("rows: %d; features: %d\n", nrow(x$data), length(x$feature_names)))
  cat(sprintf("mean |SHAP|: %.4f\n", mean(abs(shapviz::get_shap_values(x$shap)))))
  invisible(x)
}

.surv_shap_object <- function(x) {
  if (inherits(x, "unimlsurv_shap")) {
    return(x$shap)
  }
  if (inherits(x, c("shapviz", "mshapviz"))) {
    return(x)
  }
  .abort("`x` must be a `unimlsurv_shap` or a `shapviz` object.")
}

#' Plot survival SHAP importance
#'
#' @param x A `unimlsurv_shap` object.
#' @param ... Arguments passed to [uniml::plot_shap_importance()].
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' plot_surv_shap_importance(explained)
#' }
plot_surv_shap_importance <- function(x, ...) {
  uniml::plot_shap_importance(.surv_shap_object(x), ...)
}

#' @rdname plot_surv_shap_importance
#' @export
plot_surv_shap_beeswarm <- function(x, ...) {
  uniml::plot_shap_beeswarm(.surv_shap_object(x), ...)
}

#' @rdname plot_surv_shap_importance
#' @param feature Feature name for the dependence plot.
#' @export
plot_surv_shap_dependence <- function(x, feature, ...) {
  uniml::plot_shap_dependence(.surv_shap_object(x), feature, ...)
}

#' @rdname plot_surv_shap_importance
#' @param row Row of the explained data to show for the waterfall and force plots.
#' @export
plot_surv_shap_waterfall <- function(x, row = 1L, ...) {
  uniml::plot_shap_waterfall(.surv_shap_object(x), row = row, ...)
}

#' @rdname plot_surv_shap_importance
#' @export
plot_surv_shap_force <- function(x, row = 1L, ...) {
  uniml::plot_shap_force(.surv_shap_object(x), row = row, ...)
}
