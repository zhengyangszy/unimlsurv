# L8 Nomogram for survival models ---------------------------------------------
#
# The same clinical figure as `uniml::plot_nomogram()`, for the survival
# counterpart: one scale per variable, a total linear-predictor scale, and -
# when the fitted model carries its baseline cumulative hazard and the user
# chose a horizon - the survival probability that scale maps to. The mapping
# from coefficients back to original units is the same linear
# re-parametrisation: numeric variables are unscaled through the fitted
# scaling step, categorical variables are shown relative to their reference
# level.

.surv_nomogram_model <- function(result, model) {
  key <- if (identical(model, "best")) {
    result$ranking$model[[which.min(result$ranking$rank)]]
  } else {
    model
  }
  entry <- result$models[[key]]
  if (is.null(entry)) {
    .abort(sprintf("Unknown model \"%s\".", key))
  }
  fitted <- if (inherits(entry$learner, "unimlsurv_fitted")) {
    entry$learner$fit
  } else {
    entry$learner$model
  }
  if (is.null(fitted) || !inherits(fitted, "coxph")) {
    .abort(paste0(
      "A nomogram reads coefficients, so it is implemented for `coxph` ",
      "models only; model \"", key, "\" is not one."
    ))
  }
  list(key = key, coxph = fitted)
}

.surv_nomogram_coefficients <- function(coxph_fit) {
  coefficients <- stats::coef(coxph_fit)
  if (is.null(coefficients) || length(coefficients) == 0L) {
    .abort("The Cox model has no readable coefficients.")
  }
  coefficients
}

#' Nomogram table for a fitted Cox model
#'
#' The survival counterpart of [uniml::nomogram_table()]: every variable with
#' its Cox coefficient mapped back to the original units, and the total linear
#' predictor scale. When `time_point` is given and the fitted object carries
#' the data needed for the baseline cumulative hazard, the table also reports
#' the survival probability that each end of the scale maps to, through
#' `S(t | x) = exp(-H0(t) * exp(lp))`.
#'
#' @param result A `uniml_result` from [run_surv_workflow()].
#' @param model Model key or `"best"`; only `coxph` models are supported.
#' @param time_point Optional horizon for the survival row. The reporting
#'   horizons are the ones in the specification.
#'
#' @return A list with `variables`, `intercept` (always 0 for Cox; the scale
#'   is anchored at the sample mean), `total`, `probability` (when
#'   `time_point` is given), and `model`.
#' @references
#' Harrell FE Jr (2015). *Regression Modeling Strategies*, 2nd edition.
#' Springer. doi:10.1007/978-3-319-19425-7
#'
#' Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
#' Development, Validation, and Updating*, 2nd edition. Springer.
#' doi:10.1007/978-3-030-16399-0
#' @export
#' @examples
#' \dontrun{
#' result <- run_surv_workflow(spec, train, validation)
#' table <- surv_nomogram_table(result, model = "coxph", time_point = 36)
#' table$variables
#' }
surv_nomogram_table <- function(result, model = "best", time_point = NULL) {
  chosen <- .surv_nomogram_model(result, model)
  coefficients <- .surv_nomogram_coefficients(chosen$coxph)
  encoded <- names(coefficients)
  features <- coefficients

  scale_params <- result$preprocess$steps$scale$params$differences
  encoder_levels <- result$preprocess$steps$encode$params$levels
  if (is.null(encoder_levels)) {
    encoder_levels <- list()
  }
  column_to_variable <- list()
  for (variable in names(encoder_levels)) {
    for (level in encoder_levels[[variable]]) {
      column_to_variable[[paste0(variable, "_", level)]] <- variable
    }
  }
  variable_of <- vapply(encoded, function(column) {
    mapped <- column_to_variable[[column]]
    if (is.null(mapped)) column else mapped
  }, character(1))

  observed_train <- result$preprocess$train

  raw_range <- function(column) {
    difference <- scale_params[[column]]
    if (!is.null(difference)) {
      return(c(difference$center - 2 * difference$scale,
               difference$center + 2 * difference$scale))
    }
    values <- observed_train[[column]]
    if (is.null(values)) {
      return(c(0, 1))
    }
    stats::quantile(values, c(0.025, 0.975), names = FALSE)
  }

  rows <- list()
  for (variable in unique(variable_of)) {
    columns <- encoded[variable_of == variable]
    # A Cox model returns NA coefficients for columns aliased away (with
    # full-level encoding the reference dummy is exactly that); they are not
    # estimable, so they are skipped rather than drawn.
    keep <- vapply(columns, function(column) {
      coefficient <- features[[column]]
      is.finite(coefficient) && coefficient != 0
    }, logical(1))
    columns <- columns[keep]
    if (length(columns) == 0L) {
      next
    }
    is_categorical <- variable %in% names(encoder_levels)
    if (!is_categorical) {
      column <- columns[[1L]]
      difference <- scale_params[[column]]
      per_raw <- if (!is.null(difference)) {
        features[[column]] / difference$scale
      } else {
        features[[column]]
      }
      range <- raw_range(column)
      rows[[length(rows) + 1L]] <- data.frame(
        variable = variable,
        level = NA_character_,
        value_low = range[1L],
        value_high = range[2L],
        points_low = per_raw * range[1L],
        points_high = per_raw * range[2L],
        coefficient = features[[column]],
        stringsAsFactors = FALSE
      )
    } else {
      reference_column <- paste0(variable, "_", encoder_levels[[variable]][[1L]])
      reference <- if (reference_column %in% encoded) {
        features[[reference_column]]
      } else {
        0
      }
      for (column in columns) {
        level <- sub(paste0("^", variable, "_"), "", column)
        contribution <- features[[column]] - reference
        rows[[length(rows) + 1L]] <- data.frame(
          variable = variable,
          level = level,
          value_low = 0,
          value_high = 1,
          points_low = 0,
          points_high = contribution,
          coefficient = features[[column]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  variables <- do.call(rbind, rows)
  if (nrow(variables) == 0L) {
    .abort("Every coefficient was shrunk to zero; there is nothing to draw.")
  }
  total_low <- sum(pmin(variables$points_low, 0))
  total_high <- sum(pmax(variables$points_high, 0))

  probability <- NULL
  if (!is.null(time_point)) {
    .assert_scalar_number(time_point, "time_point", min = 0)
    hazard <- tryCatch(
      survival::basehaz(chosen$coxph, centered = FALSE),
      error = function(e) NULL
    )
    if (!is.null(hazard) && nrow(hazard) > 0L) {
      at_horizon <- hazard[hazard$time <= time_point, , drop = FALSE]
      if (nrow(at_horizon) > 0L) {
        h0 <- at_horizon$hazard[[nrow(at_horizon)]]
        probability <- c(
          low = exp(-h0 * exp(total_low)),
          high = exp(-h0 * exp(total_high))
        )
      }
    }
  }
  list(
    model = chosen$key,
    variables = variables,
    intercept = 0,
    total = c(low = total_low, high = total_high),
    probability = probability,
    time_point = time_point
  )
}

#' Nomogram figure for a fitted Cox model
#'
#' The survival counterpart of [uniml::plot_nomogram()]: one scale per
#' variable, the total linear-predictor scale, and - when `time_point` is
#' given and the baseline hazard is available - a survival probability row
#' through `S(t | x) = exp(-H0(t) * exp(lp))`. Read the figure as a *ranking*
#' tool: the probability row is the Cox model's own calibration, and the
#' recalibration coefficients in the workflow result say how much to trust it.
#'
#' @param result A `uniml_result` from [run_surv_workflow()].
#' @param model Model key or `"best"`; `coxph` models only.
#' @param time_point Optional horizon for the survival row; one of the
#'   specification's reporting horizons.
#' @param palette Colour for the scales.
#' @param theme A `ggplot2` theme.
#' @param title Optional title.
#'
#' @return A `ggplot` object.
#' @references
#' Harrell FE Jr (2015). *Regression Modeling Strategies*, 2nd edition.
#' Springer. doi:10.1007/978-3-319-19425-7
#' @export
#' @examples
#' \dontrun{
#' result <- run_surv_workflow(spec, train, validation)
#' plot_surv_nomogram(result, model = "coxph", time_point = 36)
#' }
plot_surv_nomogram <- function(result, model = "best", time_point = NULL,
                               palette = "bmj",
                               theme = uniml::theme_uniml(),
                               title = NULL) {
  table <- surv_nomogram_table(result, model, time_point)
  variables <- table$variables
  span <- variables$points_high - variables$points_low
  zero_span <- span <= 0
  span[zero_span] <- 1
  variables$axis_high <- ifelse(
    zero_span, 0,
    100 * (variables$points_high - variables$points_low) / max(span)
  )
  # Every scale starts at 0 on the shared points axis; its extent is axis_high.
  variables$axis_low <- 0
  variable_rows <- seq_len(nrow(variables))
  variables$y <- max(variable_rows) + 1L - variable_rows
  total_y <- 0.25
  probability_y <- if (is.null(table$probability)) -0.75 else -0.95
  logit_low <- table$total[["low"]]
  logit_span <- table$total[["high"]] - table$total[["low"]]

  colours <- .surv_palette(1)
  axis_points <- seq(0, 100, by = 10)
  total_labels <- round(logit_low + axis_points / 100 * logit_span, 2)
  tick_labels <- data.frame(
    y = total_y,
    x = axis_points,
    label = as.character(total_labels)
  )
  rows_frame <- data.frame(
    y = c(variables$y, total_y, probability_y),
    x_low = c(variables$axis_low, rep(0, 2)),
    x_high = c(variables$axis_high, rep(100, 2))
  )
  graph <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = rows_frame,
      ggplot2::aes(x = .data$x_low, xend = .data$x_high,
                   y = .data$y, yend = .data$y),
      colour = colours,
      linewidth = 0.9,
      lineend = "butt"
    ) +
    ggplot2::geom_segment(
      data = variables,
      ggplot2::aes(x = .data$axis_low, xend = .data$axis_high,
                   y = .data$y, yend = .data$y),
      colour = "grey55",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = variables,
      ggplot2::aes(x = .data$axis_low, y = .data$y),
      colour = colours,
      size = 1.6
    ) +
    ggplot2::geom_point(
      data = variables,
      ggplot2::aes(x = .data$axis_high, y = .data$y),
      colour = colours,
      size = 1.6
    ) +
    ggplot2::geom_text(
      data = tick_labels,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.6,
      vjust = 1.4,
      colour = "grey25"
    ) +
    ggplot2::geom_text(
      data = variables,
      ggplot2::aes(
        x = .data$axis_low, y = .data$y,
        label = ifelse(is.na(.data$level), sprintf("%.0f", .data$value_low),
                       .data$level)
      ),
      size = 2.6,
      vjust = -0.8,
      colour = "grey25"
    ) +
    ggplot2::geom_text(
      data = variables,
      ggplot2::aes(
        x = .data$axis_high, y = .data$y,
        label = ifelse(
          is.na(.data$level), sprintf("%.0f", .data$value_high),
          ifelse(abs(.data$points_high) < 1e-9, "",
                 sprintf("%+.2f", .data$points_high))
        )
      ),
      size = 2.6,
      vjust = -0.8,
      colour = "grey25"
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(variables$y, total_y, probability_y),
      labels = c(
        ifelse(is.na(variables$level), variables$variable,
               paste0(variables$variable, ": ", variables$level)),
        "Total points (linear predictor)",
        if (is.null(table$probability)) NULL else
          sprintf("Survival at t = %s",
                  format(table$time_point, trim = TRUE))
      ),
      expand = ggplot2::expansion(mult = 0.08)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 100, by = 10),
      labels = seq(0, 100, by = 10),
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::labs(
      x = "Points",
      y = NULL,
      title = if (is.null(title)) {
        sprintf("Nomogram - model \"%s\"", table$model)
      } else {
        title
      }
    ) +
    uniml::theme_uniml()
  graph
}
