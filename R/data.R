# L2 data layer --------------------------------------------------------------

#' Prepare a survival data set
#'
#' Splits a data frame into the feature columns and the two outcome columns,
#' checks that the follow-up time is positive and that the event indicator is
#' binary, and records which columns must never be treated as predictors.
#'
#' Rows whose time or status is missing are dropped rather than imputed: an
#' imputed outcome is not an outcome. Features may still be missing; imputation
#' of features is handled by [uniml::fit_preprocessor()].
#'
#' @param data A data frame.
#' @param spec A specification created by [surv_spec()].
#' @param name Data-set label used in messages.
#' @param require_both Whether events and censoring are both required. A training
#'   set needs them; a validation set that happens to contain only censored
#'   patients is still a valid set to predict.
#'
#' @return An object of class `unimlsurv_data`.
#' @references
#' Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
#' Development, Validation, and Updating*, 2nd edition. Springer.
#' doi:10.1007/978-3-030-16399-0
#' @export
#' @examples
#' data <- sim_surv_data(n = 60, seed = 1)
#' spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36))
#' prepared <- as_surv_data(data, spec)
#' prepared$meta$n
as_surv_data <- function(data, spec, name = "data", require_both = TRUE) {
  validate_surv_spec(spec)
  block <- spec$survival
  if (!is.data.frame(data)) {
    .abort(sprintf("`%s` must be a data frame.", name))
  }
  for (column in c(block$time, block$status)) {
    if (!column %in% names(data)) {
      .abort(sprintf("`%s` does not contain the column \"%s\".", name, column))
    }
  }
  time <- data[[block$time]]
  if (!is.numeric(time)) {
    .abort(sprintf("The time column \"%s\" must be numeric.", block$time))
  }
  status <- data[[block$status]]

  keep <- !is.na(time) & !is.na(status)
  dropped <- sum(!keep)
  if (dropped > 0L) {
    data <- data[keep, , drop = FALSE]
    time <- time[keep]
    status <- status[keep]
  }
  if (nrow(data) == 0L) {
    .abort(sprintf("`%s` has no rows with both a time and a status.", name))
  }
  if (any(time <= 0)) {
    .abort(sprintf("The time column \"%s\" must be strictly positive.", block$time))
  }

  event <- block$event
  status_values <- unique(status)
  if (is.logical(status)) {
    event_flag <- as.integer(status)
  } else if (is.numeric(status)) {
    if (!all(status_values %in% c(0, 1)) && !any(status_values %in% event)) {
      .abort(sprintf(
        "The status column \"%s\" must be 0/1 or contain the `event` value.",
        block$status
      ))
    }
    event_flag <- if (all(status_values %in% c(0, 1))) {
      as.integer(status)
    } else {
      as.integer(status %in% event)
    }
  } else {
    .abort(sprintf(
      "The status column \"%s\" must be logical or numeric.",
      block$status
    ))
  }
  if (require_both && length(unique(event_flag)) < 2L) {
    .abort(sprintf(
      "`%s` has a single status value; a survival outcome needs both events and censoring.",
      name
    ))
  }

  reserved <- unique(c(block$time, block$status, spec$id, spec$exclude))
  features <- setdiff(names(data), reserved)
  if (length(features) == 0L) {
    .abort("At least one feature column is required.")
  }
  types <- vapply(data[features], function(column) {
    if (is.numeric(column)) {
      "numeric"
    } else if (is.factor(column)) {
      "categorical"
    } else if (is.logical(column)) {
      "logical"
    } else if (is.character(column)) {
      "character"
    } else {
      "other"
    }
  }, character(1))
  unsupported <- names(types)[!types %in% c("numeric", "categorical", "logical")]
  if (length(unsupported) > 0L) {
    .abort(sprintf(
      paste(
        "Unsupported feature type(s) for: %s.",
        "Convert character columns to factors before modelling; the survival",
        "pipeline does not guess whether they are ordered."
      ),
      paste(unsupported, collapse = ", ")
    ))
  }
  if (any(types == "logical")) {
    convert <- names(types)[types == "logical"]
    data[convert] <- lapply(data[convert], function(column) {
      factor(as.integer(column), levels = c(0, 1))
    })
  }

  out <- structure(
    list(
      data = data,
      features = features,
      time = as.numeric(time),
      status = event_flag,
      time_name = block$time,
      status_name = block$status,
      time_points = block$time_points,
      spec = spec,
      meta = list(
        n = nrow(data),
        n_features = length(features),
        n_dropped = dropped,
        n_events = sum(event_flag == 1L),
        censoring_rate = mean(event_flag == 0L),
        types = types
      )
    ),
    class = c("unimlsurv_data", "list")
  )
  validate_surv_data(out)
  out
}

#' Validate a prepared survival data set
#'
#' @param x An object created by [as_surv_data()].
#'
#' @return Invisibly returns `TRUE` when valid.
#' @export
#' @examples
#' data <- sim_surv_data(n = 40, seed = 2)
#' spec <- surv_spec(time = "time", status = "status")
#' validate_surv_data(as_surv_data(data, spec))
validate_surv_data <- function(x) {
  if (!inherits(x, "unimlsurv_data")) {
    .abort("`x` must be a `unimlsurv_data` object created by `as_surv_data()`.")
  }
  if (length(x$time) != length(x$status)) {
    .abort("`x$time` and `x$status` must have the same length.")
  }
  if (!is.numeric(x$time) || anyNA(x$time) || any(x$time <= 0)) {
    .abort("`x$time` must be a positive numeric vector without missing values.")
  }
  if (!is.numeric(x$status) || !all(x$status %in% c(0, 1))) {
    .abort("`x$status` must be a 0/1 numeric vector.")
  }
  if (!is.character(x$features) || length(x$features) == 0L) {
    .abort("`x$features` must be a non-empty character vector.")
  }
  invisible(TRUE)
}

#' Print a prepared survival data set
#'
#' @param x An object created by [as_surv_data()].
#' @param ... Additional arguments, unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.unimlsurv_data <- function(x, ...) {
  cat("<unimlsurv_data>\n")
  cat(sprintf("rows: %d; features: %d\n", x$meta$n, x$meta$n_features))
  cat(sprintf(
    "time: %s (median %.1f); status: %s (%d events, %.1f%% censored)\n",
    x$time_name,
    stats::median(x$time),
    x$status_name,
    x$meta$n_events,
    100 * x$meta$censoring_rate
  ))
  cat(sprintf("time points: %s\n",
              paste(format(x$time_points), collapse = ", ")))
  if (x$meta$n_dropped > 0L) {
    cat(sprintf("dropped for missing time/status: %d\n", x$meta$n_dropped))
  }
  invisible(x)
}
