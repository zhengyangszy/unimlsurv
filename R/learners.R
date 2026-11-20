# L4 learner layer -----------------------------------------------------------

# From 0.5.0 the registry lists native engines (implemented in learners-native.R
# against `survival`). Foreign mlr3 survival learners stay available through
# `surv_register_learner()`, which is an advanced path: their training needs
# `mlr3proba`, which is not on CRAN and therefore lives in Suggests.

# Curated entries: a short alias, a readable label, and a parameter space that
# is worth searching, in the package's own plain-R format (no paradox):
# list(param = list(type = "categorical", values = ...)) or
# list(param = list(type = "numeric"|"integer", lower =, upper =, logscale =)).
.surv_curated_entries <- function() {
  list(
    coxph = list(
      key = "coxph",
      label = "Cox proportional hazards",
      engine = "coxph",
      package = "survival",
      params = list(),
      param_space = list(
        ties = list(type = "categorical",
                    values = c("efron", "breslow", "exact"))
      )
    ),
    kaplan = list(
      key = "kaplan",
      label = "Kaplan-Meier (no covariates)",
      engine = "kaplan",
      package = "survival",
      params = list(),
      param_space = NULL
    )
  )
}

#' Available survival learners
#'
#' The native engines this package ships, with a short alias you can pass to
#' [surv_set_learners()]. Each one has a readable label and, where meaningful, a
#' tuning parameter space in the package's own format. Users of the `mlr3`
#' ecosystem can add their own survival learners with [surv_register_learner()];
#' those require the `mlr3proba` package, which is not on CRAN.
#'
#' `distr` tells you whether the learner can produce survival probabilities. A
#' learner with `distr = FALSE` can only give a risk score, so it reports a
#' concordance index and nothing else.
#'
#' @param available_only Whether to keep only the learners that can be built in
#'   this session. Native engines always can; this exists for call compatibility.
#' @param curated_only Whether to keep only the curated entries.
#'
#' @return A data frame with `key`, `label`, `engine`, `package`, `available`,
#'   `distr`, and `tunable`.
#' @references
#' Sonabend R, Kirchner P, Bender A, et al. (2021). mlr3proba: an R package for
#' machine learning in survival analysis. *Bioinformatics*, 37(17), 2789-2791.
#' doi:10.1093/bioinformatics/btab039
#' @export
#' @examples
#' registry <- surv_learner_registry()
#' registry[, c("key", "label", "engine", "tunable")]
surv_learner_registry <- function(available_only = FALSE, curated_only = FALSE) {
  .assert_scalar_logical(available_only, "available_only")
  .assert_scalar_logical(curated_only, "curated_only")
  entries <- .surv_curated_entries()
  out <- data.frame(
    key = vapply(entries, `[[`, character(1), "key"),
    label = vapply(entries, `[[`, character(1), "label"),
    engine = vapply(entries, `[[`, character(1), "engine"),
    package = vapply(entries, `[[`, character(1), "package"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # The native engines only need `survival`, which is a hard dependency, so they
  # are always available and always carry a survival distribution.
  out$available <- TRUE
  out$distr <- TRUE
  out$tunable <- vapply(entries, function(entry) {
    !is.null(entry$param_space) && length(entry$param_space) > 0L
  }, logical(1))
  if (available_only) {
    out <- out[out$available, , drop = FALSE]
    row.names(out) <- NULL
  }
  out
}

#' Resolve survival learners into a workflow-ready list
#'
#' @param learners Learner keys from [surv_learner_registry()], or a named list
#'   of learners: native entries built by the registry, or foreign `mlr3`
#'   survival learners (the advanced path that needs `mlr3proba`).
#' @param tune Whether to attach a tuning parameter space.
#'
#' @return A named list of learner specifications, each with `key`, `learner`,
#'   `search_space`, `tune`, and `distr`.
#' @export
#' @examples
#' names(resolve_surv_learners(c("coxph", "kaplan"), tune = FALSE))
resolve_surv_learners <- function(learners = c("coxph", "kaplan"), tune = TRUE) {
  .assert_scalar_logical(tune, "tune")
  entries <- .surv_curated_entries()
  if (is.character(learners)) {
    .assert_character_vector(learners, "learners", allow_empty = FALSE)
    unknown <- setdiff(learners, names(entries))
    if (length(unknown) > 0L) {
      .abort(sprintf(
        paste(
          "Unknown survival learner(s): %s. The native engines are %s;",
          "mlr3 learners can be added with surv_register_learner()."
        ),
        paste(unknown, collapse = ", "),
        paste(sprintf("\"%s\"", names(entries)), collapse = ", ")
      ))
    }
    learners <- stats::setNames(learners, learners)
  } else if (is.list(learners)) {
    if (is.null(names(learners)) || any(!nzchar(names(learners)))) {
      .abort("A learner list must be named.")
    }
    supported <- vapply(learners, function(value) {
      inherits(value, "unimlsurv_native_learner") || inherits(value, "Learner")
    }, logical(1))
    if (any(!supported)) {
      .abort(sprintf(
        paste(
          "Every element of `learners` must be a native learner from the",
          "registry or an mlr3 survival learner; not: %s."
        ),
        paste(names(learners)[!supported], collapse = ", ")
      ))
    }
  } else {
    .abort("`learners` must be a character vector or a named list of learners.")
  }

  resolved <- lapply(seq_along(learners), function(index) {
    key <- names(learners)[index]
    value <- learners[[index]]
    if (is.character(value)) {
      entry <- entries[[value]]
      learner <- .new_native_learner(
        entry$key, entry$engine, entry$label, entry$params
      )
      space <- if (tune && !is.null(entry$param_space)) {
        entry$param_space
      } else {
        NULL
      }
      list(
        key = key,
        learner = learner,
        search_space = space,
        tune = tune && !is.null(space),
        distr = TRUE
      )
    } else if (inherits(value, "unimlsurv_native_learner")) {
      list(
        key = key,
        learner = value,
        search_space = NULL,
        tune = FALSE,
        distr = TRUE
      )
    } else {
      if (!identical(value$task_type, "surv") &&
          !"surv" %in% value$properties) {
        .abort(sprintf(
          "Learner \"%s\" is not a survival learner; its task type is \"%s\".",
          key, value$task_type
        ))
      }
      available <- value$predict_types
      value$predict_type <- if ("distr" %in% available) {
        "distr"
      } else {
        available[[1L]]
      }
      # Foreign mlr3 learners keep their paradox search space; tuning them runs
      # through mlr3tuning inside tune_surv(), which needs mlr3proba.
      space <- if (tune && length(value$param_set$ids()) > 0L) {
        value$param_set$search_space()
      } else {
        NULL
      }
      list(
        key = key,
        learner = value,
        search_space = space,
        tune = tune && !is.null(space) && length(space$ids()) > 0L,
        distr = "distr" %in% value$predict_types
      )
    }
  })
  names(resolved) <- vapply(resolved, `[[`, character(1), "key")
  if (anyDuplicated(names(resolved)) > 0L) {
    .abort("Learner keys must be unique.")
  }
  resolved
}

#' Register a survival learner under a key
#'
#' Use this to add a learner the registry does not ship: a native engine with
#' custom defaults, or a foreign `mlr3` survival learner. The foreign path is
#' advanced: training it needs the `mlr3proba` package, which is not on CRAN.
#'
#' @param spec A survival specification.
#' @param key A short learner key.
#' @param learner A native learner from [surv_learner_registry()], or an `mlr3`
#'   survival learner.
#' @param search_space Optional parameter space. For foreign mlr3 learners this
#'   is a [paradox::ParamSet]; native learners take the package's own plain-list
#'   format (see [surv_learner_registry()]).
#'
#' @return The updated specification.
#' @export
#' @examples
#' spec <- surv_spec(time = "time", status = "status")
#' spec <- surv_register_learner(spec, "mycox", resolve_surv_learners("coxph", tune = FALSE)$coxph$learner)
#' names(spec$survival$learners)
surv_register_learner <- function(spec, key, learner, search_space = NULL) {
  validate_surv_spec(spec)
  .assert_scalar_character(key, "key")
  native <- inherits(learner, "unimlsurv_native_learner")
  if (!native && !inherits(learner, "Learner")) {
    .abort(paste(
      "`learner` must be a native learner from surv_learner_registry() or an",
      "mlr3 survival learner."
    ))
  }
  current <- spec$survival$learners
  if (is.null(current)) {
    current <- list()
  }
  if (native) {
    current[[key]] <- list(
      key = key,
      learner = learner,
      search_space = search_space,
      tune = !is.null(search_space) && length(search_space) > 0L,
      distr = TRUE
    )
  } else {
    current[[key]] <- list(
      key = key,
      learner = learner,
      search_space = search_space,
      tune = !is.null(search_space),
      distr = "distr" %in% learner$predict_types
    )
  }
  spec$survival$learners <- current
  spec
}
