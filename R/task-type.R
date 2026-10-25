#' The survival task type
#'
#' `unimlsurv` adds `"survival"` to the task-type registry of `uniml`, so that
#' the shared specification, ranking, and export layers know what a survival
#' outcome is. The registration happens when the package is loaded; calling
#' [register_survival_task_type()] again is harmless.
#'
#' @return Invisibly returns the registered `uniml_task_type` object.
#' @references
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). Assessment and comparison
#' of prognostic classification schemes for survival data. *Statistics in
#' Medicine*, 18(17-18), 2529-2545.
#' doi:10.1002/(SICI)1097-0258(19990915/30)18:17/18<2529::AID-SIM274>3.0.CO;2-5
#'
#' Harrell FE, Lee KL, Mark DB (1996). Multivariable prognostic models: issues in
#' developing models, evaluating assumptions and adequacy, and measuring and
#' reducing errors. *Statistics in Medicine*, 15(4), 361-387.
#'
#' Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
#' Development, Validation, and Updating*, 2nd edition. Springer.
#' doi:10.1007/978-3-030-16399-0
#' @export
#' @examples
#' register_survival_task_type()
#' uniml::get_task_type("survival")$mlr3_task
register_survival_task_type <- function() {
  spec <- uniml::task_type_spec(
    id = "survival",
    mlr3_task = "surv",
    outcome_class = "time/status",
    positive = FALSE,
    metrics = c("cindex", "brier", "auc", "ibs"),
    composite_metrics = c("cindex"),
    default_learners = c("coxph", "rpart", "kaplan"),
    supports_shap = FALSE,
    supports_dca = TRUE,
    supports_calibration = TRUE,
    implemented = TRUE
  )
  uniml::register_task_type(spec, overwrite = TRUE)
  invisible(uniml::get_task_type("survival"))
}
