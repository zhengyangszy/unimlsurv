#' unimlsurv: survival analysis with mlr3
#'
#' A time-to-event counterpart to `uniml`, built on the same layering, the same
#' specification object, and the same `mlr3` machinery. See
#' `vignette("unimlsurv")` for a worked example and
#' `system.file("REFERENCES.md", package = "unimlsurv")` for the literature the
#' statistics follow.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # `mlr3proba` adds the `surv.*` learners to mlr3's dictionary from its own
  # `.onLoad()`. Importing the package only makes it available, so load the
  # namespace explicitly; otherwise `mlr3::lrn("surv.coxph")` fails on a session
  # that has not touched `mlr3proba` yet.
  loadNamespace("mlr3proba")
  register_survival_task_type()
  invisible()
}
