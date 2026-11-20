# Native survival engines -----------------------------------------------------
#
# From 0.5.0 the survival stack is native: a lightweight task class (task.R)
# and thin wrappers around the original engines on CRAN (`survival`, `rpart`).
# The tutorial literature the package follows uses exactly these engines, and
# the statistics the package reports are computed in its own metrics layer
# either way, so replacing the mlr3 plumbing changes nothing downstream:
# predictions keep the same value shape (a risk vector and an n x time_points
# survival matrix with `t<point>` column names).
#
# `mlr3proba` was LGPL-3 and not on CRAN, so its code is not reused; the
# engines below are written directly against `survival` and `rpart`.

.surv_native_engines <- function() {
  c("coxph", "kaplan")
}

# The registry hands out these objects; `engine` selects the fitting function
# below and `params` carries the hyperparameters.
.new_native_learner <- function(key, engine, label, params = list()) {
  structure(
    list(key = key, engine = engine, label = label, params = params),
    class = "unimlsurv_native_learner"
  )
}

.train_native <- function(task, learner, params) {
  merged <- utils::modifyList(learner$params, params %||% list())
  frame <- task$data
  survival_formula <- stats::as.formula(paste0(
    "survival::Surv(", task$time_col, ", ", task$status_col, ") ~ ."
  ))
  fit <- switch(learner$engine,
    coxph = survival::coxph(
      survival_formula,
      data = frame,
      ties = merged$ties %||% "efron",
      x = FALSE, y = TRUE, model = TRUE
    ),
    kaplan = survival::survfit(
      stats::as.formula(paste0(
        "survival::Surv(", task$time_col, ", ", task$status_col, ") ~ 1"
      )),
      data = frame,
      se.fit = FALSE
    ),

    .abort(sprintf("Unknown native engine \"%s\".", learner$engine))
  )
  # rpart reads its per-leaf Kaplan-Meier curves after fitting (it needs
  # fit$where); the other engines need nothing at this stage.
  leaf_fits <- if (identical(learner$engine, "rpart")) {
    .surv_leaf_fits(fit, frame, task$time_col, task$status_col)
  } else {
    NULL
  }
  structure(
    list(
      key = learner$key,
      engine = learner$engine,
      label = learner$label,
      params = merged,
      fit = fit,
      features = task$features,
      time_col = task$time_col,
      status_col = task$status_col,
      leaf_fits = leaf_fits
    ),
    class = "unimlsurv_fitted"
  )
}

#' Print a fitted native learner
#'
#' @param x A `unimlsurv_fitted`.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.unimlsurv_fitted <- function(x, ...) {
  cat(sprintf(
    "<unimlsurv_fitted> engine: %s; params: %s\n",
    x$engine,
    if (length(x$params)) {
      paste(names(x$params), x$params, sep = "=", collapse = ", ")
    } else {
      "(none)"
    }
  ))
  invisible(x)
}

# A helper that forces the survival summary into (n rows x m times). `summary()`
# lays the matrix out time-major (times in rows, curves in columns), the same
# layout the old distr6 path had, so a matching shape gets transposed.
.surv_matrix_shape <- function(values, n, m) {
  if (is.null(dim(values))) {
    values <- as.numeric(values)
    if (n == 1L && m == 1L) {
      values <- matrix(values, nrow = 1L, ncol = 1L)
    } else if (n == 1L) {
      values <- matrix(values, nrow = 1L, ncol = m)
    } else if (m == 1L) {
      values <- matrix(values, nrow = n, ncol = 1L)
    } else {
      # A plain vector with n > 1 and m > 1 is time-major, like summary().
      values <- matrix(values, nrow = m, ncol = n)
    }
  }
  out <- as.matrix(values)
  if (nrow(out) == m && ncol(out) == n) {
    out <- t(out)
  }
  if (nrow(out) != n || ncol(out) != m) {
    .abort(sprintf(
      "Survival prediction shape mismatch: got %d x %d, expected %d x %d.",
      nrow(out), ncol(out), n, m
    ))
  }
  out
}

.surv_named <- function(values, time_points) {
  colnames(values) <- paste0("t", format(time_points, trim = TRUE))
  rownames(values) <- NULL
  values
}

#' Predict survival probabilities and a risk score from a fitted learner
#'
#' The native counterpart of the old mlr3 prediction path. `survival` holds one
#' row per patient and one column per time point; `risk` is the score the
#' concordance index ranks, with a higher value meaning a shorter expected
#' survival. To predict from a full workflow result, use the existing
#' [predict_surv()].
#'
#' @param fitted A fitted native learner from [train_surv()].
#' @param newdata A data frame with the training feature columns.
#' @param time_points Numeric vector of time points.
#'
#' @return A `unimlsurv_prediction`: a list with `survival` (matrix) and `risk`
#'   (numeric vector).
#' @references
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). Assessment and comparison
#' of prognostic classification schemes for survival data. *Statistics in
#' Medicine*, 18(17-18), 2529-2545.
#' @export
#' @examples
#' data <- sim_surv_data(n = 60, seed = 1)
#' spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
#' task <- make_surv_task(data, spec)
#' fitted <- train_surv(task, surv_learner_registry()$coxph)
#' head(predict_surv_learner(fitted, data, c(12, 36))$survival)
predict_surv_learner <- function(fitted, newdata, time_points) {
  if (!inherits(fitted, "unimlsurv_fitted")) {
    .abort("`fitted` must be a fitted native learner from train_surv().")
  }
  .assert_numeric_vector(time_points, "time_points", min = 0)
  frame <- newdata[, fitted$features, drop = FALSE]
  n <- nrow(frame)
  m <- length(time_points)
  out <- switch(fitted$engine,
    coxph = {
      lp <- as.numeric(stats::predict(fitted$fit, newdata = frame, type = "lp"))
      sf <- survival::survfit(fitted$fit, newdata = frame)
      values <- summary(sf, times = time_points)$surv
      list(survival = values, risk = lp)
    },
    kaplan = {
      values <- summary(fitted$fit, times = time_points)$surv
      values <- matrix(
        rep(as.numeric(values), times = n),
        nrow = n,
        ncol = m,
        byrow = TRUE
      )
      list(survival = values, risk = rep(0, n))
    },

    .abort(sprintf("Unknown native engine \"%s\".", fitted$engine))
  )
  structure(
    list(
      survival = .surv_named(.surv_matrix_shape(out$survival, n, m),
                             time_points),
      risk = as.numeric(out$risk),
      engine = fitted$engine
    ),
    class = "unimlsurv_prediction"
  )
}
