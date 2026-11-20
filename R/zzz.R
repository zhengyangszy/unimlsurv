#' unimlsurv: survival analysis in the shape of uniml
#'
#' A time-to-event counterpart to `uniml`, built on the same layering, the same
#' specification object, and the same workflow engine. Native engines
#' (`survival`, `rpart`) do the fitting; see `vignette("unimlsurv")` for a
#' worked example and `system.file("REFERENCES.md", package = "unimlsurv")`
#' for the literature the statistics follow.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  register_survival_task_type()
  invisible()
}
