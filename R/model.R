# L7 tuning and cross-validation ---------------------------------------------
#
# Tuning is a repeated cross-validated search implemented here: candidate
# settings are drawn from the registry's plain-R parameter spaces, every
# candidate is scored by the package's own concordance index on pooled
# out-of-fold risk predictions (the Varma and Simon honest layout), and the
# winner is returned for the workflow's train hook to fit. Foreign mlr3
# learners keep the mlr3tuning path, which needs mlr3proba (Suggests).

# 折切分：按事件指示分层，事件与非事件各自均匀铺到 folds 里 ----
.surv_cv_pairs <- function(task, folds, repeats, seed) {
  status <- as.numeric(task$data[[task$status_col]])
  n <- length(status)
  pairs <- list()
  for (rep in seq_len(repeats)) {
    set.seed(as.integer(seed) + rep - 1L)
    fold_of <- integer(n)
    for (group in list(which(status == 1L), which(status != 1L))) {
      if (length(group) == 0L) {
        next
      }
      group <- group[sample.int(length(group))]
      fold_of[group] <- rep(seq_len(folds), length.out = length(group))
    }
    for (f in seq_len(folds)) {
      test_idx <- which(fold_of == f)
      pairs[[length(pairs) + 1L]] <- list(
        train = setdiff(seq_len(n), test_idx),
        test = test_idx,
        repetition = as.integer(rep),
        fold = as.integer(f)
      )
    }
  }
  pairs
}

.surv_task_subset <- function(task, index) {
  structure(
    list(
      data = task$data[index, , drop = FALSE],
      features = task$features,
      time_col = task$time_col,
      status_col = task$status_col,
      n = length(index),
      n_events = sum(as.numeric(task$data[[task$status_col]][index]) == 1L)
    ),
    class = "unimlsurv_task"
  )
}

# 候选生成：全分类空间铺网格，含数值维度时随机采样 evals 个 ----
.surv_tune_candidates <- function(space, evals, seed) {
  if (is.null(space) || length(space) == 0L) {
    return(list(list()))
  }
  numeric_dims <- vapply(space, function(d) {
    !identical(d$type, "categorical")
  }, logical(1))
  if (!any(numeric_dims)) {
    combos <- expand.grid(
      lapply(space, function(d) d$values),
      stringsAsFactors = FALSE
    )
    set.seed(as.integer(seed))
    if (nrow(combos) > evals) {
      combos <- combos[sample.int(nrow(combos), evals), , drop = FALSE]
    }
    candidates <- lapply(seq_len(nrow(combos)), function(i) {
      as.list(combos[i, , drop = FALSE])
    })
    return(candidates)
  }
  set.seed(as.integer(seed))
  lapply(seq_len(evals), function(i) {
    lapply(space, function(d) {
      if (identical(d$type, "categorical")) {
        return(sample(d$values, 1L))
      }
      value <- if (isTRUE(d$logscale)) {
        stats::runif(1L, log(d$lower), log(d$upper))
      } else {
        stats::runif(1L, d$lower, d$upper)
      }
      if (isTRUE(d$logscale)) {
        value <- exp(value)
      }
      if (identical(d$type, "integer")) {
        value <- max(d$lower, min(d$upper, round(value)))
      }
      value
    })
  })
}

# The old mlr3tuning path, kept for foreign mlr3 learners only. It needs
# mlr3proba (the surv measures live in its dictionary) and mlr3tuning; both are
# GitHub-only or optional, so they load through the namespace, not ::.
.tune_surv_mlr3 <- function(task, learner_spec, spec, measure_id = "cindex") {
  for (pkg in c("mlr3proba", "mlr3tuning", "mlr3")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      .abort(sprintf(paste0(
        "Tuning a foreign mlr3 survival learner needs %s. The native ",
        "engines tune without it."
      ), pkg))
    }
  }
  ns_mlr3 <- asNamespace("mlr3")
  ns_tuning <- asNamespace("mlr3tuning")
  learner <- learner_spec$learner
  rsmp <- getExportedValue(ns_mlr3, "rsmp")
  msr <- getExportedValue(ns_mlr3, "msr")
  tnr <- getExportedValue(ns_tuning, "tnr")
  trm <- getExportedValue(ns_tuning, "trm")
  auto_tuner <- getExportedValue(ns_tuning, "auto_tuner")
  resampling <- rsmp("cv", folds = max(2L, as.integer(spec$tuning$folds)))
  tuner <- tnr(spec$tuning$tuner)
  terminator <- if (identical(spec$tuning$terminator, "seconds")) {
    trm("run_time", secs = as.numeric(spec$tuning$seconds))
  } else {
    trm("evals", n_evals = as.integer(spec$tuning$evals))
  }
  tunable <- auto_tuner(
    tuner = tuner,
    learner = learner$clone(deep = TRUE),
    resampling = resampling,
    measure = msr(sprintf("surv.%s", measure_id)),
    terminator = terminator,
    search_space = learner_spec$search_space
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
      measure = measure_id
    ))
  }
  list(
    learner = learner$clone(deep = TRUE),
    best_params = result$learner_param_vals[[1L]],
    best_metric = as.numeric(result[[sprintf("surv.%s", measure_id)]][[1L]]),
    status = "tuned",
    measure = measure_id
  )
}

#' Tune a survival learner
#'
#' A repeated cross-validated search over the learner's parameter space, scored
#' by the package's own concordance index on pooled out-of-fold risk
#' predictions. Learners without a space are returned untuned rather than being
#' given one that would not mean anything. Foreign mlr3 learners tune through
#' `mlr3tuning`, which requires `mlr3proba`.
#'
#' @param task A native task from [make_surv_task()].
#' @param learner_spec An entry from [resolve_surv_learners()].
#' @param spec A survival specification.
#'
#' @return A list with `learner`, `best_params`, `best_metric`, `status`, and
#'   `measure`. The learner is returned unfitted; the workflow's train hook
#'   fits it with `best_params`.
#' @references
#' Probst P, Boulesteix AL, Bischl B (2019). Tunability: importance of
#' hyperparameters. *Journal of Machine Learning Research*, 20(53), 1-32.
#'
#' Varma S, Simon R (2006). Bias in error estimation when using cross-validation
#' for model selection. *BMC Bioinformatics*, 7, 91. doi:10.1186/1471-2105-7-91
#' @export
#' @examples
#' \dontrun{
#' data <- sim_surv_data(n = 150, seed = 1)
#' spec <- surv_spec(time = "time", status = "status")
#' task <- make_surv_task(data, spec)
#' learners <- resolve_surv_learners("coxph", tune = TRUE)
#' tuned <- tune_surv(task, learners$coxph, spec)
#' tuned$best_params
#' }
tune_surv <- function(task, learner_spec, spec) {
  learner <- learner_spec$learner
  untuned <- list(
    learner = learner,
    best_params = list(),
    best_metric = NA_real_,
    status = "untuned",
    measure = "cindex"
  )
  if (!isTRUE(learner_spec$tune) || is.null(learner_spec$search_space) ||
      length(learner_spec$search_space) == 0L) {
    return(untuned)
  }
  if (inherits(learner, "Learner")) {
    # Foreign learners need a TaskSurv and the mlr3 measure dictionary.
    if (!inherits(task, "TaskSurv")) {
      task <- .task_surv_to_mlr3(task)
    }
    return(.tune_surv_mlr3(task, learner_spec, spec))
  }
  if (!inherits(task, "unimlsurv_task")) {
    .abort("`task` must be a native task from make_surv_task().")
  }
  folds <- max(2L, as.integer(spec$tuning$folds))
  repeats <- max(1L, as.integer(spec$tuning$repeats %||% 1L))
  evals <- max(1L, as.integer(spec$tuning$evals))
  pairs <- .surv_cv_pairs(task, folds, repeats, spec$seed)
  candidates <- .surv_tune_candidates(learner_spec$search_space, evals,
                                      spec$seed)
  time <- as.numeric(task$data[[task$time_col]])
  status <- as.integer(as.numeric(task$data[[task$status_col]]))
  scores <- vapply(candidates, function(candidate) {
    oof_risk <- rep(NA_real_, task$n)
    for (pair in pairs) {
      fitted <- tryCatch(
        .train_native(.surv_task_subset(task, pair$train), learner, candidate),
        error = function(e) e
      )
      if (inherits(fitted, "error")) {
        return(NA_real_)
      }
      predicted <- tryCatch(
        predict_surv_learner(
          fitted, task$data[pair$test, , drop = FALSE],
          spec$survival$time_points
        ),
        error = function(e) e
      )
      if (inherits(predicted, "error")) {
        return(NA_real_)
      }
      oof_risk[pair$test] <- predicted$risk
    }
    if (anyNA(oof_risk)) {
      return(NA_real_)
    }
    surv_cindex(time, status, oof_risk)
  }, numeric(1))
  if (all(is.na(scores))) {
    return(list(
      learner = learner,
      best_params = list(),
      best_metric = NA_real_,
      status = "failed",
      measure = "cindex"
    ))
  }
  best <- which.max(scores)
  list(
    learner = learner,
    best_params = candidates[[best]],
    best_metric = as.numeric(scores[[best]]),
    status = "tuned",
    measure = "cindex"
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
  pairs <- .surv_cv_pairs(task, as.integer(folds), as.integer(repeats),
                          spec$seed)

  rows <- list()
  for (key in names(learners)) {
    entry <- learners[[key]]
    learner <- entry$learner
    if (!is.null(entry$search_space) && isTRUE(entry$tune)) {
      tuned <- tryCatch(
        tune_surv(task, entry, spec),
        error = function(e) e
      )
      if (inherits(tuned, "error")) {
        next
      }
      if (length(tuned$best_params) > 0L) {
        learner <- if (inherits(learner, "unimlsurv_native_learner")) {
          .new_native_learner(
            learner$key, learner$engine, learner$label,
            utils::modifyList(learner$params, tuned$best_params)
          )
        } else {
          learner
        }
      }
    }
    oof_risk <- rep(NA_real_, task$n)
    oof_survival <- matrix(NA_real_, nrow = task$n, ncol = length(time_points))
    keep <- TRUE
    for (pair in pairs) {
      fitted <- tryCatch(
        train_surv(.surv_task_subset(task, pair$train), learner),
        error = function(e) e
      )
      if (inherits(fitted, "error")) {
        keep <- FALSE
        break
      }
      predicted <- tryCatch(
        predict_surv_learner(fitted, task$data[pair$test, , drop = FALSE],
                             time_points),
        error = function(e) e
      )
      if (inherits(predicted, "error")) {
        keep <- FALSE
        break
      }
      oof_risk[pair$test] <- predicted$risk
      oof_survival[pair$test, ] <- predicted$survival
    }
    if (!keep || anyNA(oof_risk)) {
      next
    }
    truth <- task$data[, c(task$time_col, task$status_col)]
    metrics <- calc_surv_metrics(
      time = as.numeric(truth[[task$time_col]]),
      status = as.integer(as.numeric(truth[[task$status_col]])),
      risk = oof_risk,
      survival = oof_survival,
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
