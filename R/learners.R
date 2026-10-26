# L4 learner layer -----------------------------------------------------------

# Curated entries: a short alias, a readable label, and a search space that is
# worth searching. Everything else that `mlr3` exposes as `surv.*` is discovered
# at runtime and usable, just without a curated space.
.surv_curated_entries <- function() {
  list(
    coxph = list(
      key = "coxph",
      label = "Cox proportional hazards",
      package = "survival",
      mlr3_key = "surv.coxph",
      search_space = function(learner) {
        paradox::ps(ties = paradox::p_fct(c("efron", "breslow", "exact")))
      }
    ),
    rpart = list(
      key = "rpart",
      label = "Recursive partitioning (discrimination only)",
      package = "rpart",
      mlr3_key = "surv.rpart",
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
      mlr3_key = "surv.kaplan",
      search_space = function(learner) NULL
    ),
    ranger = list(
      key = "ranger",
      label = "Random survival forest",
      package = "mlr3extralearners",
      mlr3_key = "surv.ranger",
      search_space = function(learner) {
        paradox::ps(
          mtry.ratio = paradox::p_dbl(0.1, 0.9),
          min.node.size = paradox::p_int(3L, 20L),
          sample.fraction = paradox::p_dbl(0.5, 1)
        )
      }
    ),
    cb = list(
      key = "cb",
      label = "Component-wise likelihood boosting",
      package = "mlr3extralearners",
      mlr3_key = "surv.coxboost",
      search_space = function(learner) {
        paradox::ps(
          stepno = paradox::p_int(50L, 500L),
          penalty = paradox::p_dbl(10, 1000, logscale = TRUE)
        )
      }
    ),
    glmnet = list(
      key = "glmnet",
      label = "Penalised Cox model (cross-validated lambda)",
      package = "mlr3extralearners",
      mlr3_key = "surv.cv_glmnet",
      search_space = function(learner) {
        paradox::ps(alpha = paradox::p_dbl(0, 1))
      }
    ),
    xgboost = list(
      key = "xgboost",
      label = "Gradient boosting, Cox objective",
      package = "mlr3extralearners",
      mlr3_key = "surv.xgboost.cox",
      search_space = function(learner) {
        paradox::ps(
          eta = paradox::p_dbl(0.01, 0.3, logscale = TRUE),
          max_depth = paradox::p_int(2L, 8L),
          subsample = paradox::p_dbl(0.5, 1),
          nrounds = paradox::p_int(50L, 500L)
        )
      }
    ),
    rfsrc = list(
      key = "rfsrc",
      label = "Random survival forests (randomForestSRC)",
      package = "mlr3extralearners",
      mlr3_key = "surv.rfsrc",
      search_space = function(learner) {
        paradox::ps(
          mtry = paradox::p_int(1L, 20L),
          nodesize = paradox::p_int(3L, 20L),
          ntree = paradox::p_int(200L, 1000L)
        )
      }
    ),
    aorsf = list(
      key = "aorsf",
      label = "Oblique random survival forest",
      package = "mlr3extralearners",
      mlr3_key = "surv.aorsf",
      search_space = function(learner) {
        paradox::ps(
          mtry = paradox::p_int(1L, 20L),
          n_tree = paradox::p_int(100L, 500L)
        )
      }
    ),
    cforest = list(
      key = "cforest",
      label = "Conditional inference forest",
      package = "mlr3extralearners",
      mlr3_key = "surv.cforest",
      search_space = function(learner) {
        paradox::ps(mtry = paradox::p_int(1L, 20L))
      }
    ),
    gamboost = list(
      key = "gamboost",
      label = "Likelihood-based boosting with splines",
      package = "mlr3extralearners",
      mlr3_key = "surv.gamboost",
      search_space = function(learner) {
        paradox::ps(mstop = paradox::p_int(50L, 500L))
      }
    ),
    survdnn = list(
      key = "survdnn",
      label = "Deep neural network",
      package = "mlr3extralearners",
      mlr3_key = "surv.survdnn",
      search_space = function(learner) {
        paradox::ps(
          epochs = paradox::p_int(20L, 200L),
          batch_size = paradox::p_int(16L, 128L),
          learning_rate = paradox::p_dbl(1e-04, 0.01, logscale = TRUE)
        )
      }
    )
  )
}

.surv_learning_keys <- function() {
  # `mlr3extralearners` registers its learners in its own `.onLoad()`, so the
  # dictionary only lists them once that namespace has been loaded.
  loadNamespace("mlr3proba")
  if (requireNamespace("mlr3extralearners", quietly = TRUE)) {
    loadNamespace("mlr3extralearners")
  }
  keys <- tryCatch(mlr3::mlr_learners$keys(), error = function(e) character())
  sort(grep("^surv\\.", keys, value = TRUE))
}

# A learner object can be constructed even when the package it delegates to is
# missing; `$packages` is what decides whether it can actually be trained.
.surv_missing_packages <- function(learner) {
  packages <- learner$packages
  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
}

# Every `surv.*` learner mlr3 knows about, minus the ones a curated alias already
# covers. The alias is the mlr3 key without its prefix, which is what users see in
# the mlr3 documentation.
.surv_discovered_entries <- function() {
  curated <- .surv_curated_entries()
  covered <- vapply(curated, `[[`, character(1), "mlr3_key")
  keys <- setdiff(.surv_learning_keys(), covered)
  if (length(keys) == 0L) {
    return(list())
  }
  entries <- lapply(keys, function(key) {
    alias <- sub("^surv\\.", "", key)
    list(
      key = alias,
      label = key,
      package = NA_character_,
      mlr3_key = key,
      search_space = function(learner) NULL
    )
  })
  stats::setNames(entries, vapply(entries, `[[`, character(1), "key"))
}

.surv_learner_entries <- function() {
  entries <- c(.surv_curated_entries(), .surv_discovered_entries())
  # A curated alias always wins if a discovered one would shadow it.
  entries[!duplicated(names(entries))]
}

.surv_build_entry <- function(entry) {
  if (!is.null(entry$mlr3_key)) {
    if (!is.na(entry$package) &&
        !entry$package %in% c("survival", "rpart")) {
      loadNamespace(entry$package)
    }
    return(mlr3::lrn(entry$mlr3_key))
  }
  entry$build()
}

# A learner may be registered in mlr3's dictionary only after its package has
# been loaded, so the namespace is loaded first and any failure is reported as
# "not available" rather than as an error the user cannot act on.
.try_build_surv_learner <- function(entry) {
  # Building a learner whose optional package is missing emits a warning and then
  # fails; both are expected here, because the registry reports availability.
  suppressWarnings(tryCatch(.surv_build_entry(entry), error = function(e) NULL))
}

.surv_has_distr <- function(learner) {
  "distr" %in% learner$predict_types
}

.surv_set_predict_type <- function(learner) {
  available <- learner$predict_types
  learner$predict_type <- if ("distr" %in% available) "distr" else available[[1L]]
  learner
}

#' Available survival learners
#'
#' Everything `mlr3` exposes as a `surv.*` learner, with a short alias you can pass
#' to [surv_set_learners()]. The first block is curated: it has a readable label
#' and a tuning search space. The rest is discovered at run time and is fitted at
#' its defaults unless you register a search space with
#' [surv_register_learner()].
#'
#' `distr` tells you whether the learner can produce survival probabilities. A
#' learner with `distr = FALSE` can only give a risk score, so it reports a
#' concordance index and nothing else.
#'
#' @param available_only Whether to keep only the learners that can be built in
#'   this session.
#' @param curated_only Whether to keep only the curated entries.
#'
#' @return A data frame with `key`, `label`, `mlr3_key`, `package`, `available`,
#'   `distr`, and `tunable`.
#' @references
#' Sonabend R, Kirchner P, Bender A, et al. (2021). mlr3proba: an R package for
#' machine learning in survival analysis. *Bioinformatics*, 37(17), 2789-2791.
#' doi:10.1093/bioinformatics/btab039
#' @export
#' @examples
#' registry <- surv_learner_registry()
#' sum(registry$available)
#' subset(registry, distr, select = c(key, mlr3_key))[1:5, ]
surv_learner_registry <- function(available_only = FALSE, curated_only = FALSE) {
  .assert_scalar_logical(available_only, "available_only")
  .assert_scalar_logical(curated_only, "curated_only")
  entries <- if (curated_only) {
    .surv_curated_entries()
  } else {
    .surv_learner_entries()
  }
  built <- lapply(entries, .try_build_surv_learner)
  out <- data.frame(
    key = vapply(entries, `[[`, character(1), "key"),
    label = vapply(entries, `[[`, character(1), "label"),
    mlr3_key = vapply(entries, `[[`, character(1), "mlr3_key"),
    package = vapply(entries, function(entry) {
      if (is.na(entry$package)) NA_character_ else entry$package
    }, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out$available <- vapply(seq_along(entries), function(index) {
    !is.null(built[[index]]) && length(.surv_missing_packages(built[[index]])) == 0L
  }, logical(1))
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

#' Resolve survival learners into a workflow-ready list
#'
#' @param learners Learner keys from [surv_learner_registry()], or a named list
#'   of `mlr3` survival learners.
#' @param tune Whether to build a tuning search space.
#'
#' @return A named list of learner specifications, each with `key`, `learner`,
#'   `search_space`, `tune`, and `distr`.
#' @export
#' @examples
#' names(resolve_surv_learners(c("coxph", "kaplan"), tune = FALSE))
resolve_surv_learners <- function(learners = c("coxph", "rpart"), tune = TRUE) {
  .assert_scalar_logical(tune, "tune")
  entries <- .surv_learner_entries()
  if (is.character(learners)) {
    .assert_character_vector(learners, "learners", allow_empty = FALSE)
    unknown <- setdiff(learners, names(entries))
    if (length(unknown) > 0L) {
      .abort(sprintf(
        "Unknown survival learner(s): %s. See `surv_learner_registry()` for the %d available keys.",
        paste(unknown, collapse = ", "),
        length(entries)
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
      if (!is.null(raw)) {
        missing_packages <- .surv_missing_packages(raw)
        if (length(missing_packages) > 0L) {
          .abort(sprintf(
            "The learner \"%s\" needs the package(s) %s, which are not installed.",
            value, paste(sprintf("\"%s\"", missing_packages), collapse = ", ")
          ))
        }
      }
      if (is.null(raw)) {
        .abort(sprintf(
          paste(
            "The learner \"%s\" could not be built.",
            "It needs the package \"%s\"; install it, or run",
            "`surv_learner_registry(available_only = TRUE)` for what works here."
          ),
          value, entry$package %||% entry$mlr3_key
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
#' Use this when you want a learner that is not in the registry, or when you want
#' a search space for one of the discovered learners that ships without one.
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
    tune = !is.null(search_space),
    distr = .surv_has_distr(learner)
  )
  spec$survival$learners <- current
  spec
}
