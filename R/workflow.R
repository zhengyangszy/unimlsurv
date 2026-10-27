# L8 workflow ----------------------------------------------------------------

.surv_adapter <- function() {
  list(
    prepare = function(data, spec, name = "data", levels = NULL,
                        require_two = identical(name, "train")) {
      prepared <- as_surv_data(
        data, spec, name,
        require_both = identical(name, "train")
      )
      list(
        data = prepared$data,
        outcome = spec$survival$status,
        levels = NULL,
        survival = prepared
      )
    },
    make_task = function(data, spec, features) {
      block <- spec$survival
      frame <- data[, c(features, block$time, block$status), drop = FALSE]
      mlr3proba::as_task_surv(
        frame,
        time = block$time,
        event = block$status,
        id = "unimlsurv_task"
      )
    },
    predict = function(learner, task, spec) {
      time_points <- spec$survival$time_points
      prediction <- learner$predict(task)
      survival <- if (is.null(prediction$distr)) {
        NULL
      } else {
        survival_matrix(prediction, time_points)
      }
      list(
        prediction = prediction,
        value = list(
          survival = survival,
          risk = survival_risk(prediction)
        )
      )
    },
    pool = function(values) {
      if (length(values) == 1L) {
        return(values[[1L]])
      }
      risk <- Reduce(`+`, lapply(values, `[[`, "risk")) / length(values)
      survival <- NULL
      if (!any(vapply(values, function(value) is.null(value$survival), logical(1)))) {
        survival <- Reduce(`+`, lapply(values, `[[`, "survival")) / length(values)
      }
      list(survival = survival, risk = risk)
    },
    features = function(data, spec) {
      block <- spec$survival
      setdiff(names(data), c(block$time, block$status))
    },
    metric_names = function(spec) {
      c("cindex", "ibs", "brier", "auc")
    },
    metrics = function(value, truth, spec, metrics, context = list()) {
      time <- context$time
      status <- context$status
      if (is.null(time) || is.null(status)) {
        .abort("The survival metric adapter needs `context$time` and `context$status`.")
      }
      time_points <- context$time_points %||% spec$survival$time_points
      result <- calc_surv_metrics(
        time = as.numeric(time),
        status = as.integer(status),
        risk = value$risk,
        survival = value$survival,
        time_points = time_points,
        censoring = spec$survival$censoring_model
      )
      as.list(result$row)
    },
    train = function(task, learner, params, spec) {
      train_surv(task, learner, params)
    },
    tune = function(task, learner_spec, spec) {
      tune_surv(task, learner_spec, spec)
    }
  )
}

#' Run a survival workflow
#'
#' The survival counterpart of [uniml::run_workflow()]. It runs the shared
#' [uniml::workflow_engine()], so preprocessing, imputation, and tuning are fitted
#' on the training data only and reused untouched on the validation, test, and
#' external sets, then scores every model with the survival metrics this package
#' defines.
#'
#' The validation set selects the best model; the test and external sets report.
#' Nothing from them is used for fitting, tuning, or model selection.
#'
#' @param spec A specification created by [surv_spec()].
#' @param train Training data frame.
#' @param validation Validation data frame; used for model selection.
#' @param test Optional test data frame; reporting-only.
#' @param external Optional external data frame; reporting-only.
#' @param ci_R Bootstrap resamples for the confidence intervals on the comparison
#'   set. Set to `0` to skip them.
#'
#' @return A `uniml_result` with `metrics` (one row per model and set), `dca`
#'   (time-point decision curves on the comparison set), `ranking` (models ordered
#'   by the concordance index), `calibration` (time-point calibration), `models`
#'   (fitted learners by name), and `meta` (the specification, predictions, and
#'   time points).
#' @references
#' Collins GS, Moons KGM, Dhiman P, et al. (2024). TRIPOD+AI statement: updated
#' guidance for reporting clinical prediction models that use regression or
#' machine learning methods. *BMJ*, 385, e078378. doi:10.1136/bmj-2023-078378
#'
#' Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
#' Development, Validation, and Updating*, 2nd edition. Springer.
#' doi:10.1007/978-3-030-16399-0
#' @export
#' @examples
#' \dontrun{
#' data <- sim_surv_data(n = 300, seed = 1)
#' index <- sample.int(nrow(data))
#' spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36, 60))
#' spec <- surv_set_learners(spec, "coxph", tune = FALSE)
#' result <- run_surv_workflow(
#'   spec,
#'   train = data[index[1:200], ],
#'   validation = data[index[201:300], ]
#' )
#' result$ranking
#' }
run_surv_workflow <- function(spec, train, validation, test = NULL,
                              external = NULL, ci_R = 200L) {
  validate_surv_spec(spec)
  .assert_scalar_number(ci_R, "ci_R", integer = TRUE, min = 0)
  block <- spec$survival
  learners <- block$learners
  if (is.null(learners)) {
    learners <- resolve_surv_learners(c("coxph", "rpart"), tune = FALSE)
  }
  time_points <- block$time_points

  # The preprocessor must leave the two outcome columns alone. `vars` restricts
  # it to the feature columns, which are known only once the training frame has
  # been read.
  prepared_train <- as_surv_data(train, spec, "train")
  spec$preprocess$vars <- prepared_train$features
  # And there is no survival feature selector: `uniml`'s selectors screen the
  # outcome column, which here would be the event indicator, so they would ignore
  # the follow-up time entirely. Set it here as well as in `surv_spec()`, so a
  # specification assembled by hand or by another setter cannot reintroduce one.
  spec$selection$methods <- .surv_no_selection()

  engine <- uniml::workflow_engine(
    spec,
    train,
    validation,
    test,
    external,
    adapter = .surv_adapter(),
    learners = learners
  )
  if (length(engine$states) == 0L) {
    failures <- .surv_failure_table(engine$failures)
    detail <- if (nrow(failures) == 0L) {
      "No failure was logged, so the problem is upstream of the model loop."
    } else {
      paste0(
        "\n  ",
        paste(
          sprintf("%s [%s]: %s", failures$model, failures$stage, failures$message),
          collapse = "\n  "
        )
      )
    }
    .abort(paste0("No model could be fitted.", detail))
  }

  sets <- engine$prepared
  metric_rows <- list()
  by_time_rows <- list()
  predictions <- list()
  notes <- character()
  for (key in names(engine$states)) {
    run_states <- engine$states[[key]]
    available <- unique(unlist(lapply(run_states, function(state) {
      names(state$values)
    })))
    for (set in names(sets)) {
      if (is.null(sets[[set]]) || !set %in% available) {
        next
      }
      set_states <- run_states[vapply(
        run_states,
        function(state) set %in% names(state$values),
        logical(1)
      )]
      pooled <- engine$adapter$pool(
        lapply(set_states, function(state) state$values[[set]]$value)
      )
      if (is.null(pooled$survival)) {
        notes <- c(notes, sprintf(
          paste(
            "Model \"%s\" predicts a risk score only (predict type \"%s\"),",
            "so Brier, AUC, calibration, and decision curves are not reported",
            "for it; the concordance index is."
          ),
          key,
          paste(run_states[[1L]]$learner$predict_types, collapse = "/")
        ))
      }
      truth <- sets[[set]]$data
      time <- truth[[block$time]]
      status <- as.integer(truth[[block$status]])
      scored <- calc_surv_metrics(
        time = as.numeric(time),
        status = status,
        risk = pooled$risk,
        survival = pooled$survival,
        time_points = time_points,
        censoring = block$censoring_model
      )
      row <- data.frame(
        model = key,
        set = set,
        scored$row,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      metric_rows[[length(metric_rows) + 1L]] <- row
      by_time <- scored$by_time
      by_time$model <- key
      by_time$set <- set
      by_time_rows[[length(by_time_rows) + 1L]] <- by_time[, c(
        "model", "set", "time", "brier", "auc", "calib_intercept",
        "calib_slope", "n_case", "n_control"
      )]
      predictions[[paste0(key, ":", set)]] <- list(
        model = key,
        set = set,
        time = as.numeric(time),
        status = status,
        risk = pooled$risk,
        survival = pooled$survival
      )
    }
  }
  metrics <- if (length(metric_rows) == 0L) {
    data.frame(model = character(), set = character(), stringsAsFactors = FALSE)
  } else {
    out <- do.call(rbind, metric_rows)
    rownames(out) <- NULL
    out
  }

  select_on <- if ("validation" %in% metrics$set) {
    "validation"
  } else if ("test" %in% metrics$set) {
    "test"
  } else {
    NA_character_
  }
  ranking <- if (is.na(select_on)) {
    NULL
  } else {
    uniml::rank_models(
      metrics,
      select_on = select_on,
      best_model_metric = "composite",
      task_type = "survival"
    )
  }

  calibration <- .surv_workflow_calibration(
    predictions, spec, select_on, names(engine$states)
  )
  dca <- .surv_workflow_dca(predictions, spec, select_on)
  ci <- if (ci_R > 0L) {
    .surv_workflow_ci(predictions, spec, select_on, as.integer(ci_R))
  } else {
    NULL
  }

  uniml::assemble_result(
    data_info = list(
      task_type = "survival",
      outcome = block$status,
      time = block$time,
      time_points = time_points,
      positive = NULL,
      n_train = nrow(sets$train$data),
      n_validation = nrow(sets$validation$data),
      n_test = if (is.null(sets$test)) 0L else nrow(sets$test$data),
      n_external = if (is.null(sets$external)) 0L else nrow(sets$external$data)
    ),
    preprocess = engine$preprocessor,
    selection = engine$selection,
    tuning = engine$tuning,
    metrics = metrics,
    ci = ci,
    calibration = calibration,
    dca = dca,
    ranking = ranking,
    models = .surv_models(engine),
    meta = list(
      spec = spec,
      predictions = predictions,
      time_points = time_points,
      by_time = if (length(by_time_rows) == 0L) {
        NULL
      } else {
        do.call(rbind, by_time_rows)
      },
      failures = .surv_failure_table(engine$failures),
      notes = unique(notes),
      select_on = select_on
    ),
    select_on = if (is.na(select_on)) "auto" else select_on,
    best_model_metric = "composite"
  )
}

.surv_workflow_ci <- function(predictions, spec, select_on, R) {
  if (is.na(select_on)) {
    return(NULL)
  }
  rows <- list()
  for (entry in predictions) {
    if (!identical(entry$set, select_on)) {
      next
    }
    intervals <- tryCatch(
      surv_ci(
        entry$time,
        entry$status,
        entry$risk,
        entry$survival,
        spec$survival$time_points,
        censoring = spec$survival$censoring_model,
        R = R,
        seed = spec$seed
      ),
      error = function(e) NULL
    )
    if (is.null(intervals)) {
      next
    }
    intervals$model <- entry$model
    intervals$set <- entry$set
    rows[[length(rows) + 1L]] <- intervals[, c(
      "model", "set", "metric", "mean", "se", "lower", "upper",
      "n_effective", "formatted"
    )]
  }
  if (length(rows) == 0L) {
    return(NULL)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.surv_failure_table <- function(log) {
  empty <- data.frame(
    model = character(), stage = character(), message = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(log) || length(log$rows) == 0L) {
    return(empty)
  }
  out <- do.call(rbind, log$rows)
  out <- unique(out)
  rownames(out) <- NULL
  out
}

.surv_models <- function(engine) {
  out <- list()
  for (key in names(engine$states)) {
    run_states <- engine$states[[key]]
    entry <- list(
      learner = run_states[[1L]]$learner,
      tuned = run_states[[1L]]$tuned,
      preprocessor = run_states[[1L]]$preprocessor,
      selection = run_states[[1L]]$selection,
      task = run_states[[1L]]$task
    )
    if (length(run_states) > 1L) {
      entry$learners <- lapply(run_states, function(state) state$learner)
      entry$imputations <- vapply(
        run_states, function(state) as.integer(state$imputation), integer(1)
      )
    }
    out[[key]] <- entry
  }
  out
}

.surv_workflow_calibration <- function(predictions, spec, select_on, models) {
  if (is.na(select_on)) {
    return(NULL)
  }
  time_points <- spec$survival$time_points
  rows <- list()
  summaries <- list()
  for (entry in predictions) {
    if (!identical(entry$set, select_on) || is.null(entry$survival)) {
      next
    }
    curves <- surv_calibration_curves(
      entry$time, entry$status, entry$survival, time_points,
      censoring = spec$survival$censoring_model
    )
    tables <- curves$table
    if (nrow(tables) > 0L) {
      tables$model <- entry$model
      tables$set <- entry$set
      rows[[length(rows) + 1L]] <- tables[, c(
        "model", "set", "time", "group", "n", "predicted", "observed",
        "difference"
      )]
    }
    summary <- curves$summary
    summary$model <- entry$model
    summary$set <- entry$set
    summaries[[length(summaries) + 1L]] <- summary[, c(
      "model", "set", "time", "ici", "slope", "intercept", "groups", "n",
      "events", "message"
    )]
  }
  if (length(rows) == 0L && length(summaries) == 0L) {
    return(NULL)
  }
  out <- if (length(rows) == 0L) {
    data.frame(
      model = character(), set = character(), time = numeric(),
      group = integer(), n = integer(), predicted = numeric(),
      observed = numeric(), difference = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    combined <- do.call(rbind, rows)
    rownames(combined) <- NULL
    combined
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  list(
    method = "kaplan_meier_groups",
    # Survival calibration has no per-model fitted parameters; the empty entries
    # keep the shape that `uniml::export_results()` expects.
    parameters = stats::setNames(rep(list(list()), length(models)), models),
    models = models,
    set = select_on,
    metrics = out,
    summary = summary
  )
}

.surv_workflow_dca <- function(predictions, spec, select_on) {
  if (is.na(select_on)) {
    return(NULL)
  }
  rows <- list()
  for (entry in predictions) {
    if (!identical(entry$set, select_on) || is.null(entry$survival)) {
      next
    }
    curve <- surv_dca_curves(
      entry$time, entry$status, entry$survival, spec$survival$time_points,
      censoring = spec$survival$censoring_model
    )
    if (nrow(curve) == 0L) {
      next
    }
    curve$model <- entry$model
    rows[[length(rows) + 1L]] <- curve[, c(
      "model", "time", "threshold", "strategy", "net_benefit"
    )]
  }
  if (length(rows) == 0L) {
    return(NULL)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Predict survival probabilities with a fitted workflow
#'
#' @param result A `uniml_result` returned by [run_surv_workflow()].
#' @param newdata A data frame with the same feature columns as the training data.
#' @param model Model key, or `"best"` to use the top-ranked model.
#' @param time_points Time points; defaults to the ones recorded in the result.
#'
#' @return A data frame with `time`, `risk`, and one `S_<t>` column per time point.
#' @export
#' @examples
#' \dontrun{
#' result <- run_surv_workflow(spec, train, validation)
#' predict_surv(result, test)
#' }
predict_surv <- function(result, newdata, model = "best",
                         time_points = NULL) {
  if (!inherits(result, "uniml_result")) {
    .abort("`result` must be a `uniml_result` from `run_surv_workflow()`.")
  }
  spec <- result$meta$spec
  validate_surv_spec(spec)
  if (is.null(time_points)) {
    time_points <- result$meta$time_points
  }
  .assert_numeric_vector(time_points, "time_points", min = 0)
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
  block <- spec$survival
  if (!is.data.frame(newdata)) {
    .abort("`newdata` must be a data frame.")
  }
  if (!"distr" %in% entry$learner$predict_types) {
    eligible <- names(Filter(
      function(candidate) "distr" %in% candidate$learner$predict_types,
      models
    ))
    .abort(sprintf(
      paste(
        "Model \"%s\" predicts a risk score only, so it cannot produce survival",
        "probabilities. Use `survival_risk()` on a prediction instead, or pick a",
        "model that carries a distribution%s."
      ),
      model,
      if (length(eligible) == 0L) {
        " (none of the fitted models does)"
      } else {
        sprintf(": %s", paste(sprintf("\"%s\"", eligible), collapse = ", "))
      }
    ))
  }
  missing_columns <- setdiff(c(block$time, block$status), names(newdata))
  if (length(missing_columns) > 0L) {
    .abort(sprintf(
      "`newdata` must contain: %s.",
      paste(missing_columns, collapse = ", ")
    ))
  }
  features <- setdiff(
    names(newdata),
    c(block$time, block$status, spec$id, spec$exclude)
  )
  frame <- uniml::apply_preprocessor(
    entry$preprocessor,
    newdata[, c(features, block$time, block$status), drop = FALSE]
  )
  # The learner was trained on exactly these columns; anything else is either a
  # new feature it has never seen or a feature it expects and cannot find.
  required <- entry$task$feature_names
  missing_training <- setdiff(required, names(frame))
  if (length(missing_training) > 0L) {
    .abort(sprintf(
      paste(
        "The preprocessor did not produce the column(s) the model was trained on:",
        "%s.",
        "The new data probably lack a level or a column seen during training."
      ),
      paste(missing_training, collapse = ", ")
    ))
  }
  frame <- frame[, c(required, block$time, block$status), drop = FALSE]
  task <- mlr3proba::as_task_surv(
    frame,
    time = block$time,
    event = block$status,
    id = "unimlsurv_newdata"
  )
  prediction <- entry$learner$predict(task)
  survival <- survival_matrix(prediction, time_points)
  out <- data.frame(risk = survival_risk(prediction), stringsAsFactors = FALSE)
  for (index in seq_along(time_points)) {
    out[[paste0("S_", format(time_points[[index]], trim = TRUE))]] <-
      survival[, index]
  }
  out
}
