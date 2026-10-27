#' Build a survival specification
#'
#' A survival specification is the `uniml_spec` object you already know, with
#' `task_type = "survival"` plus a `survival` block that records the time column,
#' the status column, and the time points at which survival probabilities are
#' reported. Registrations, preprocessing, feature selection, tuning, and
#' evaluation are configured with the `uniml` setters, so the two packages share
#' one mental model.
#'
#' The status column is stored as `spec$outcome`, which keeps the object
#' compatible with [uniml::assemble_result()] and [uniml::export_results()]. The
#' time column and the time points live in `spec$survival`.
#'
#' @param time Name of the follow-up time column.
#' @param status Name of the event indicator column; `1` means the event
#'   happened, `0` means the observation was censored.
#' @param time_points Time points, in the units of `time`, at which survival
#'   probabilities, the Brier score, the time-dependent AUC, calibration, and
#'   decision curves are reported. They are always chosen by the user; the
#'   package never invents them.
#' @param id Optional row identifier column.
#' @param exclude Optional columns to keep out of the model.
#' @param event Value of the status column that means "event". Defaults to `1`;
#'   pass `TRUE` for logical status columns.
#' @param censoring_model How the censoring distribution is estimated for the
#'   inverse-probability-of-censoring weights: `"marginal"` (Kaplan-Meier on the
#'   whole cohort, the default) or `"none"` (no weighting, only valid when
#'   censoring is absent or independent of everything).
#' @param seed Random seed recorded for reproducible tuning.
#'
#' @return A `uniml_spec` object with class `c("unimlsurv_spec", "uniml_spec")`.
#' @references
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). Assessment and comparison
#' of prognostic classification schemes for survival data. *Statistics in
#' Medicine*, 18(17-18), 2529-2545.
#'
#' Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
#' Development, Validation, and Updating*, 2nd edition. Springer.
#' doi:10.1007/978-3-030-16399-0
#' @export
#' @examples
#' spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
#' spec$survival$time_points
#' spec$task_type
surv_spec <- function(time = "time", status = "status",
                      time_points = c(12, 36, 60), id = NULL,
                      exclude = NULL, event = 1, censoring_model = "marginal",
                      seed = 2026L) {
  .assert_scalar_character(time, "time")
  .assert_scalar_character(status, "status")
  if (identical(time, status)) {
    .abort("`time` and `status` must be different columns.")
  }
  .assert_numeric_vector(time_points, "time_points", min = 0)
  .assert_choice(censoring_model, "censoring_model", c("marginal", "none"))
  .assert_scalar_number(seed, "seed", integer = TRUE, min = 1)

  spec <- uniml::uniml_spec(
    task_type = "survival",
    outcome = status,
    id = id,
    exclude = exclude,
    seed = seed
  )
  spec$survival <- list(
    time = time,
    status = status,
    event = event,
    time_points = sort(unique(as.numeric(time_points))),
    censoring_model = censoring_model
  )
  # `uniml`'s feature selection screens on the outcome column, which for a
  # survival model is the event indicator: it would ignore the time entirely and
  # rank features by a binary outcome. Survival selection is therefore off by
  # default; every feature is used unless the user selects them.
  spec$selection$methods <- "none"
  class(spec) <- c("unimlsurv_spec", class(spec))
  validate_surv_spec(spec)
  spec
}

#' Validate a survival specification
#'
#' @param x A specification created by [surv_spec()].
#'
#' @return Invisibly returns `TRUE` when valid.
#' @export
#' @examples
#' validate_surv_spec(surv_spec(time = "time", status = "dead"))
validate_surv_spec <- function(x) {
  if (!inherits(x, "uniml_spec")) {
    .abort("`x` must be a `uniml_spec` created by `surv_spec()`.")
  }
  uniml::validate_uniml_spec(x)
  if (!identical(x$task_type, "survival")) {
    .abort("A survival specification must have `task_type = \"survival\"`.")
  }
  block <- x$survival
  if (!is.list(block)) {
    .abort("The specification has no `survival` block; create it with `surv_spec()`.")
  }
  required <- c("time", "status", "event", "time_points", "censoring_model")
  missing_fields <- setdiff(required, names(block))
  if (length(missing_fields) > 0L) {
    .abort(sprintf(
      "The survival block is missing: %s.",
      paste(missing_fields, collapse = ", ")
    ))
  }
  .assert_scalar_character(block$time, "survival$time")
  .assert_scalar_character(block$status, "survival$status")
  if (identical(block$time, block$status)) {
    .abort("`survival$time` and `survival$status` must be different columns.")
  }
  .assert_numeric_vector(block$time_points, "survival$time_points", min = 0)
  .assert_choice(block$censoring_model, "survival$censoring_model",
                 c("marginal", "none"))
  invisible(TRUE)
}

.survival_block <- function(spec) {
  validate_surv_spec(spec)
  spec$survival
}

#' Print a survival specification
#'
#' @param x A specification created by [surv_spec()].
#' @param ... Additional arguments, unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.unimlsurv_spec <- function(x, ...) {
  block <- x$survival
  cat("<unimlsurv_spec>\n")
  cat(sprintf("time: %s; status: %s (event = %s)\n",
              block$time, block$status, paste(block$event, collapse = ", ")))
  cat(sprintf("time points: %s\n",
              paste(format(block$time_points), collapse = ", ")))
  cat(sprintf("censoring model: %s\n", block$censoring_model))
  cat(sprintf("impute: %s; scale: %s; encode: %s\n",
              x$preprocess$impute, x$preprocess$scale, x$preprocess$encode))
  cat(sprintf("tuning folds: %s; evals: %s; metric: %s\n",
              x$tuning$folds, x$tuning$evals, x$tuning$metric))
  invisible(x)
}

#' Set the reporting time points
#'
#' @param spec A survival specification.
#' @param time_points Numeric vector of time points, in the units of the follow-up
#'   time.
#'
#' @return The updated specification.
#' @export
#' @examples
#' surv_set_time_points(surv_spec(time = "time", status = "status"), c(6, 24))
surv_set_time_points <- function(spec, time_points) {
  validate_surv_spec(spec)
  .assert_numeric_vector(time_points, "time_points", min = 0)
  spec$survival$time_points <- sort(unique(as.numeric(time_points)))
  spec
}

#' Set the censoring model used for inverse-probability weights
#'
#' @param spec A survival specification.
#' @param censoring_model `"marginal"` or `"none"`.
#'
#' @return The updated specification.
#' @export
#' @examples
#' spec <- surv_spec(time = "time", status = "status")
#' surv_set_censoring(spec, "none")$survival$censoring_model
surv_set_censoring <- function(spec, censoring_model = "marginal") {
  validate_surv_spec(spec)
  .assert_choice(censoring_model, "censoring_model", c("marginal", "none"))
  spec$survival$censoring_model <- censoring_model
  spec
}

#' Set the survival learners
#'
#' @param spec A survival specification.
#' @param learners Learner keys from [surv_learner_registry()], or a named list of
#'   `mlr3` survival learners.
#' @param tune Whether to tune the learners that have a search space.
#'
#' @return The updated specification.
#' @export
#' @examples
#' spec <- surv_spec(time = "time", status = "status")
#' surv_set_learners(spec, "coxph")$survival$learners$coxph$key
surv_set_learners <- function(spec, learners = c("coxph", "rpart"), tune = TRUE) {
  validate_surv_spec(spec)
  .assert_scalar_logical(tune, "tune")
  spec$survival$learners <- resolve_surv_learners(learners, tune = tune)
  spec
}

#' Set the tuning control
#'
#' A thin wrapper around [uniml::set_tuning()]. Tuning uses the concordance index,
#' which is the only one of the reported statistics that does not need the
#' training set to be recomputed: the inverse-probability weights of the Brier
#' score depend on the censoring distribution of the data the model was fitted on,
#' so the Brier score is reported, not tuned.
#'
#' @param spec A survival specification.
#' @param ... Arguments passed to [uniml::set_tuning()].
#' @param metric Measure used for tuning; only `"cindex"` is implemented.
#'
#' @return The updated specification.
#' @export
#' @examples
#' spec <- surv_spec(time = "time", status = "status")
#' surv_set_tuning(spec, folds = 3, evals = 5)$tuning$folds
surv_set_tuning <- function(spec, ..., metric = "cindex") {
  validate_surv_spec(spec)
  if (!identical(metric, "cindex")) {
    .abort(sprintf(
      paste(
        "Survival tuning uses the concordance index, not \"%s\".",
        "The Brier score is a reporting metric: its inverse-probability weights",
        "depend on the censoring distribution of the training data, so tuning on",
        "it would compare models fitted on different data."
      ),
      metric
    ))
  }
  updated <- uniml::set_tuning(spec, ...)
  updated$survival$tuning_metric <- metric
  validate_surv_spec(updated)
  updated
}
