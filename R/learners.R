# L4 learner layer -----------------------------------------------------------

.surv_learner_entries <- function() {
  list(
    coxph = list(
      key = "coxph",
      label = "Cox proportional hazards",
      package = "survival",
      build = function() {
        mlr3::lrn("surv.coxph")
      },
      search_space = function(learner) {
        paradox::ps(ties = paradox::p_fct(c("efron", "breslow", "exact")))
      }
    ),
    rpart = list(
      key = "rpart",
      label = "Recursive partitioning (discrimination only)",
      package = "rpart",
      build = function() {
        mlr3::lrn("surv.rpart")
      },
      search_space = function(learner) {
        paradox::ps(
          minsplit = paradox::p_int(10L, 60L),
          minbucket = paradox::p_int(5L, 30L),
          cp = paradox::p_dbl(1e-04, 0.1, logscale = TRUE),
          maxdepth = paradox::p_int(2L, 10L)
        )
      }
    ),
    kaplan = list(
      key = "kaplan",
      label = "Kaplan-Meier (no covariates)",
      package = "survival",
      build = function() {
        mlr3::lrn("surv.kaplan")
      },
      search_space = function(learner) NULL
    ),
    ranger = list(
      key = "ranger",
      label = "Random survival forest",
      package = "mlr3extralearners",
      build = function() {
        mlr3::lrn("surv.ranger")
      },
      search_space = function(learner) {
        paradox::ps(
          mtry.ratio = paradox::p_dbl(0.1, 0.9),
          min.node.size = paradox::p_int(3L, 20L),
          sample.fraction = paradox::p_dbl(0.5, 1)
        )
      }
    ),
    xgboost = list(
      key = "xgboost",
      label = "Gradient boosting",
      package = "mlr3extralearners",
      build = function() {
        mlr3::lrn("surv.xgboost")
      },
      search_space = function(learner) {
        paradox::ps(
          eta = paradox::p_dbl(0.01, 0.3, logscale = TRUE),
          max_depth = paradox::p_int(2L, 8L),
          subsample = paradox::p_dbl(0.5, 1),
          nrounds = paradox::p_int(50L, 500L)
        )
      }
    ),
    glmnet = list(
      key = "glmnet",
      label = "Penalised Cox model",
      package = "mlr3extralearners",
      build = function() {
        mlr3::lrn("surv.glmnet")
      },
      search_space = function(learner) {
        paradox::ps(alpha = paradox::p_dbl(0, 1), lambda = paradox::p_dbl(-4, 0))
      }
    )
  )
}

#' Available survival learners
#'
#' Lists the learner keys that [surv_set_learners()] understands, together with
#' the package each one needs and whether it has a tuning search space.
#'
#' @param available_only Whether to keep only the learners whose package is
#'   installed.
#'
#' @return A data frame with `key`, `label`, `package`, `available`, and
#'   `tunable`.
#' @export
#' @examples
#' head(surv_learner_registry())
surv_learner_registry <- function(available_only = FALSE) {
  .assert_scalar_logical(available_only, "available_only")
  entries <- .surv_learner_entries()
  out <- data.frame(
    key = vapply(entries, `[[`, character(1), "key"),
    label = vapply(entries, `[[`, character(1), "label"),
    package = vapply(entries, `[[`, character(1), "package"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  built <- lapply(entries, .try_build_surv_learner)
  out$available <- !vapply(built, is.null, logical(1))
  out$distr <- vapply(seq_along(entries), function(index) {
    if (!out$available[[index]]) {
      return(NA)
    }
    .surv_has_distr(built[[index]])
  }, logical(1))
  out$tunable <- vapply(seq_along(entries), function(index) {
    if (!out$available[[index]]) {
      return(NA)
    }
    !is.null(entries[[index]]$search_space(built[[index]]))
  }, logical(1))
  if (available_only) {
    out <- out[out$available, , drop = FALSE]
    row.names(out) <- NULL
  }
  out
}

# A learner may be registered in mlr3's dictionary only after its package has
# been loaded, so the namespace is loaded first and any failure is reported as
# "not available" rather than as an error the user cannot act on.
.try_build_surv_learner <- function(entry) {
  tryCatch({
    if (!identical(entry$package, "survival") &&
        !identical(entry$package, "rpart")) {
      loadNamespace(entry$package)
    }
    entry$build()
  }, error = function(e) NULL)
}

.surv_has_distr <- function(learner) {
  "distr" %in% learner$predict_types
}

.surv_set_predict_type <- function(learner) {
  available <- learner$predict_types
  learner$predict_type <- if ("distr" %in% available) "distr" else available[[1L]]
  learner
}

#' Resolve survival learners into a workflow-ready list
#'
#' @param learners Learner keys from [surv_learner_registry()], or a named list
#'   of `mlr3` survival learners.
#' @param tune Whether to build a tuning search space.
#'
#' @return A named list of learner specifications, each with `key`, `learner`,
#'   `search_space`, and `tune`.
#' @export
#' @examples
#' names(resolve_surv_learners("coxph", tune = FALSE))
resolve_surv_learners <- function(learners = c("coxph", "rpart"), tune = TRUE) {
  .assert_scalar_logical(tune, "tune")
  entries <- .surv_learner_entries()
  if (is.character(learners)) {
    .assert_character_vector(learners, "learners", allow_empty = FALSE)
    unknown <- setdiff(learners, names(entries))
    if (length(unknown) > 0L) {
      .abort(sprintf(
        "Unknown survival learner(s): %s. Available: %s.",
        paste(unknown, collapse = ", "),
        paste(names(entries), collapse = ", ")
      ))
    }
    learners <- stats::setNames(learners, learners)
  } else if (is.list(learners)) {
    if (is.null(names(learners)) || any(!nzchar(names(learners)))) {
      .abort("A learner list must be named.")
    }
    bad <- !vapply(learners, inherits, logical(1), "Learner")
    if (any(bad)) {
      .abort(sprintf(
        "Every element of `learners` must be an mlr3 learner; not: %s.",
        paste(names(learners)[bad], collapse = ", ")
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
      raw <- .try_build_surv_learner(entry)
      if (is.null(raw)) {
        .abort(sprintf(
          paste(
            "The learner \"%s\" could not be built.",
            "It needs the package \"%s\"; install it, or check",
            "`surv_learner_registry()` for what is available."
          ),
          value, entry$package
        ))
      }
      built <- .surv_set_predict_type(raw)
      space <- if (tune) entry$search_space(built) else NULL
      list(
        key = key,
        learner = built,
        search_space = space,
        tune = tune && !is.null(space),
        distr = .surv_has_distr(built)
      )
    } else {
      if (!identical(value$task_type, "surv") &&
          !"surv" %in% value$properties) {
        .abort(sprintf(
          "Learner \"%s\" is not a survival learner; its task type is \"%s\".",
          key, value$task_type
        ))
      }
      value <- .surv_set_predict_type(value)
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
        distr = .surv_has_distr(value)
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
#' @param spec A survival specification.
#' @param key A short learner key.
#' @param learner An `mlr3` survival learner.
#' @param search_space Optional [paradox::ParamSet] used for tuning.
#'
#' @return The updated specification.
#' @export
#' @examples
#' spec <- surv_spec(time = "time", status = "status")
#' spec <- surv_register_learner(spec, "mycox", mlr3::lrn("surv.coxph"))
#' names(spec$survival$learners)
surv_register_learner <- function(spec, key, learner, search_space = NULL) {
  validate_surv_spec(spec)
  .assert_scalar_character(key, "key")
  if (!inherits(learner, "Learner")) {
    .abort("`learner` must be an mlr3 learner.")
  }
  current <- spec$survival$learners
  if (is.null(current)) {
    current <- list()
  }
  current[[key]] <- list(
    key = key,
    learner = learner,
    search_space = search_space,
    tune = !is.null(search_space)
  )
  spec$survival$learners <- current
  spec
}
