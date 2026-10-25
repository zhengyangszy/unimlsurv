# L6 plots -------------------------------------------------------------------

.surv_palette <- function(n) {
  colours <- c(
    "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
    "#8C564B", "#17BECF", "#E377C2"
  )
  if (n <= length(colours)) {
    colours[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, "Dark 3")
  }
}

#' Kaplan-Meier curves
#'
#' @param data A data frame, or a `unimlsurv_data` object.
#' @param spec A specification created by [surv_spec()].
#' @param group Optional column name to compare groups.
#' @param conf_int Whether to draw the pointwise confidence band.
#' @param palette Optional vector of colours.
#'
#' @return A `ggplot` object.
#' @references
#' Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete
#' observations. *Journal of the American Statistical Association*, 53(282),
#' 457-481. doi:10.1080/01621459.1958.10501452
#'
#' Pocock SJ, Clayton TC, Altman DG (2002). Survival plots of time-to-event
#' outcomes in clinical trials: good practice and pitfalls. *The Lancet*,
#' 359(9318), 1686-1689. doi:10.1016/S0140-6736(02)08594-X
#' @export
#' @examples
#' data <- sim_surv_data(n = 150, seed = 1)
#' spec <- surv_spec(time = "time", status = "status")
#' plot_surv_km(data, spec, group = "sex")
plot_surv_km <- function(data, spec, group = NULL, conf_int = FALSE,
                         palette = NULL) {
  validate_surv_spec(spec)
  block <- spec$survival
  frame <- if (inherits(data, "unimlsurv_data")) data$data else data
  if (!is.data.frame(frame)) {
    .abort("`data` must be a data frame or a `unimlsurv_data` object.")
  }
  for (column in c(block$time, block$status)) {
    if (!column %in% names(frame)) {
      .abort(sprintf("`data` does not contain the column \"%s\".", column))
    }
  }
  frame <- data.frame(
    time = as.numeric(frame[[block$time]]),
    status = as.integer(frame[[block$status]])
  )
  if (is.null(group)) {
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = frame)
    curves <- .surv_curve_frame(fit)
  } else {
    .assert_scalar_character(group, "group")
    if (!group %in% names(data)) {
      .abort(sprintf("`data` does not contain the grouping column \"%s\".", group))
    }
    frame$group <- factor(data[[group]])
    fit <- survival::survfit(survival::Surv(time, status) ~ group, data = frame)
    curves <- .surv_curve_frame(fit)
  }
  colours <- palette %||% .surv_palette(length(unique(curves$strata)))
  graph <- ggplot2::ggplot(
    curves,
    ggplot2::aes(x = .data$time, y = .data$surv, colour = .data$strata)
  )
  if (conf_int && all(c("lower", "upper") %in% names(curves))) {
    graph <- graph + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$strata),
      alpha = 0.15,
      colour = NA
    )
  }
  graph +
    ggplot2::geom_step(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Time",
      y = "Survival probability",
      colour = if (is.null(group)) NULL else group,
      fill = if (is.null(group)) NULL else group
    ) +
    uniml::theme_uniml()
}

.surv_curve_frame <- function(fit) {
  summary_fit <- summary(fit, censored = TRUE)
  labels <- if (is.null(summary_fit$strata)) {
    rep("all", length(summary_fit$time))
  } else {
    # With strata, `summary.survfit()` returns one label per output row.
    sub("^[^=]*=", "", as.character(summary_fit$strata))
  }
  frame <- data.frame(
    time = as.numeric(summary_fit$time),
    surv = as.numeric(summary_fit$surv),
    strata = factor(labels, levels = unique(labels)),
    stringsAsFactors = FALSE
  )
  if (!is.null(summary_fit$lower)) {
    frame$lower <- as.numeric(summary_fit$lower)
    frame$upper <- as.numeric(summary_fit$upper)
  }
  frame
}

#' Time-dependent ROC curve
#'
#' Cumulative/dynamic ROC curve with inverse-probability-of-censoring weights.
#' The area under this curve is the value [surv_auc()] returns, which is pinned by
#' a test.
#'
#' @param time Follow-up times.
#' @param status Event indicator.
#' @param survival Predicted survival probability at `time_point`.
#' @param time_point The time point.
#' @param censoring Censoring model.
#'
#' @return A data frame with `threshold`, `fpr`, and `tpr`, ordered by threshold.
#' @references
#' Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC curves.
#' *Biometrics*, 61(1), 92-105. doi:10.1111/j.0006-341X.2005.030814.x
#'
#' Fawcett T (2006). An introduction to ROC analysis. *Pattern Recognition
#' Letters*, 27(8), 861-874. doi:10.1016/j.patrec.2005.10.010
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, 36)
#' curve <- surv_roc_curve(data$time, data$status, truth$survival[, 1], 36)
#' head(curve)
surv_roc_curve <- function(time, status, survival, time_point,
                           censoring = "marginal") {
  .check_survival_inputs(time, status, survival = survival)
  .assert_scalar_number(time_point, "time_point", min = 0)
  case <- status == 1L & time <= time_point
  control <- time > time_point
  if (!any(case) || !any(control)) {
    return(data.frame(
      threshold = numeric(), fpr = numeric(), tpr = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  G <- .censoring_survival(time, status, censoring)
  weight_case <- 1 / G(time[case])
  weight_control <- 1 / G(rep(time_point, sum(control)))
  risk <- 1 - survival
  cuts <- sort(unique(risk), decreasing = TRUE)
  if (length(cuts) == 0L) {
    return(data.frame(
      threshold = numeric(), fpr = numeric(), tpr = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(cuts, function(cut) {
    flagged <- risk >= cut
    data.frame(
      threshold = 1 - cut,
      tpr = sum(weight_case * flagged[case]) / sum(weight_case),
      fpr = sum(weight_control * flagged[control]) / sum(weight_control),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$fpr, out$tpr), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Plot the time-dependent ROC curves of a workflow
#'
#' @param result A `uniml_result` from [run_surv_workflow()].
#' @param set Which data set to plot.
#' @param models Optional model keys; defaults to all of them.
#' @param palette Optional vector of colours.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' result <- run_surv_workflow(spec, train, validation)
#' plot_surv_roc(result)
#' }
plot_surv_roc <- function(result, set = NULL, models = NULL, palette = NULL) {
  entries <- .surv_predictions(result, set = set, models = models)
  entries <- Filter(function(entry) !is.null(entry$survival), entries)
  if (length(entries) == 0L) {
    .abort(paste(
      "No predictions with survival probabilities are stored for the requested",
      "set and models. A learner that predicts only a risk score cannot be plotted."
    ))
  }
  spec <- result$meta$spec
  time_points <- result$meta$time_points
  rows <- list()
  for (entry in entries) {
    for (index in seq_along(time_points)) {
      curve <- surv_roc_curve(
        entry$time, entry$status, entry$survival[, index], time_points[[index]],
        censoring = spec$survival$censoring_model
      )
      if (nrow(curve) == 0L) {
        next
      }
      curve$model <- entry$model
      curve$time_point <- time_points[[index]]
      rows[[length(rows) + 1L]] <- curve
    }
  }
  if (length(rows) == 0L) {
    .abort("No ROC curve could be computed for the requested time points.")
  }
  frame <- do.call(rbind, rows)
  frame$panel <- sprintf("%s (t = %s)", frame$model,
                         format(frame$time_point, trim = TRUE))
  colours <- palette %||% .surv_palette(length(unique(frame$panel)))
  ggplot2::ggplot(
    frame,
    ggplot2::aes(x = .data$fpr, y = .data$tpr, colour = .data$panel)
  ) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey70") +
    ggplot2::geom_path(linewidth = 0.7) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
    uniml::theme_uniml()
}

#' Plot time-point calibration
#'
#' @param calibration Either `result$calibration$metrics`, a calibration table
#'   with `predicted` and `observed` columns, or a `uniml_result`.
#' @param set Which data set to plot when a result is supplied.
#' @param models Optional model keys.
#' @param palette Optional vector of colours.
#'
#' @return A `ggplot` object.
#' @references
#' Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW (2019).
#' Calibration: the Achilles heel of predictive analytics. *BMC Medicine*, 17,
#' 230. doi:10.1186/s12916-019-1466-7
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, c(12, 36))
#' curves <- surv_calibration_curves(data$time, data$status, truth$survival,
#'                                   c(12, 36))
#' plot_surv_calibration(curves$table)
plot_surv_calibration <- function(calibration, set = NULL, models = NULL,
                                  palette = NULL) {
  frame <- .surv_calibration_table(calibration, set = set, models = models)
  if (is.null(frame) || nrow(frame) == 0L) {
    .abort("No calibration table is available to plot.")
  }
  frame$panel <- sprintf("%s (t = %s)", frame$model,
                         format(frame$time, trim = TRUE))
  colours <- palette %||% .surv_palette(length(unique(frame$panel)))
  ggplot2::ggplot(
    frame,
    ggplot2::aes(x = .data$predicted, y = .data$observed, colour = .data$panel)
  ) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey70") +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(ggplot2::aes(group = .data$panel), linewidth = 0.6) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(
      x = "Predicted survival",
      y = "Observed survival (Kaplan-Meier)",
      colour = NULL
    ) +
    uniml::theme_uniml()
}

#' Plot time-point decision curves
#'
#' @param dca Either `result$dca`, a decision-curve table, or a `uniml_result`.
#' @param set Which data set to plot when a result is supplied.
#' @param models Optional model keys.
#' @param palette Optional vector of colours.
#'
#' @return A `ggplot` object.
#' @references
#' Vickers AJ, van Calster B, Steyerberg EW (2016). Net benefit approaches to the
#' evaluation of prediction models, molecular markers, and diagnostic tests.
#' *BMJ*, 352, i6. doi:10.1136/bmj.i6
#' @export
#' @examples
#' data <- sim_surv_data(n = 200, seed = 1)
#' truth <- sim_surv_truth(data, 36)
#' curve <- surv_dca(data$time, data$status, truth$survival[, 1], 36)
#' plot_surv_dca(curve)
plot_surv_dca <- function(dca, set = NULL, models = NULL, palette = NULL) {
  frame <- .surv_dca_table(dca, set = set, models = models)
  if (is.null(frame) || nrow(frame) == 0L) {
    .abort("No decision curve is available to plot.")
  }
  frame$panel <- sprintf("%s (t = %s)", frame$model,
                         format(frame$time, trim = TRUE))
  frame$strategy <- factor(
    frame$strategy,
    levels = c("model", "treat_all", "treat_none"),
    labels = c("Model", "Treat all", "Treat none")
  )
  colours <- palette %||% .surv_palette(length(unique(frame$panel)))
  ggplot2::ggplot(
    frame,
    ggplot2::aes(
      x = .data$threshold, y = .data$net_benefit,
      colour = .data$panel, linetype = .data$strategy
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(
      x = "Threshold probability",
      y = "Net benefit",
      colour = NULL,
      linetype = NULL
    ) +
    uniml::theme_uniml()
}

#' Plot the Brier score and AUC over time
#'
#' @param result A `uniml_result` from [run_surv_workflow()].
#' @param set Which data set to plot.
#' @param models Optional model keys.
#' @param palette Optional vector of colours.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' result <- run_surv_workflow(spec, train, validation)
#' plot_surv_time_metrics(result, set = "validation")
#' }
plot_surv_time_metrics <- function(result, set = "validation", models = NULL,
                                   palette = NULL) {
  frame <- result$meta$by_time
  if (is.null(frame) || nrow(frame) == 0L) {
    .abort("The result carries no per-time metrics.")
  }
  if (!is.null(set)) {
    frame <- frame[frame$set == set, , drop = FALSE]
  }
  if (!is.null(models)) {
    frame <- frame[frame$model %in% models, , drop = FALSE]
  }
  if (nrow(frame) == 0L) {
    .abort("No per-time metrics are available for the requested set and models.")
  }
  long <- rbind(
    data.frame(frame[, c("model", "time")],
               metric = "Brier score", value = frame$brier,
               stringsAsFactors = FALSE),
    data.frame(frame[, c("model", "time")],
               metric = "Time-dependent AUC", value = frame$auc,
               stringsAsFactors = FALSE)
  )
  colours <- palette %||% .surv_palette(length(unique(long$model)))
  ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$time, y = .data$value, colour = .data$model)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(x = "Time", y = NULL, colour = NULL) +
    uniml::theme_uniml()
}

.surv_predictions <- function(result, set = NULL, models = NULL) {
  entries <- result$meta$predictions
  if (is.null(entries)) {
    .abort("The result carries no predictions.")
  }
  entries <- Filter(function(entry) {
    (is.null(set) || identical(entry$set, set)) &&
      (is.null(models) || entry$model %in% models)
  }, entries)
  if (is.null(set)) {
    sets <- unique(vapply(entries, `[[`, character(1), "set"))
    preferred <- c("validation", "test", "external", "train")
    chosen <- preferred[preferred %in% sets][[1L]]
    entries <- Filter(function(entry) identical(entry$set, chosen), entries)
  }
  entries
}

.surv_calibration_table <- function(calibration, set = NULL, models = NULL) {
  frame <- if (inherits(calibration, "uniml_result")) {
    calibration$calibration$metrics
  } else if (is.list(calibration) && !is.null(calibration$metrics)) {
    calibration$metrics
  } else if (is.data.frame(calibration)) {
    calibration
  } else {
    NULL
  }
  if (is.null(frame)) {
    return(NULL)
  }
  if (!all(c("predicted", "observed") %in% names(frame))) {
    .abort("A calibration table needs `predicted` and `observed` columns.")
  }
  if (is.null(frame$model)) {
    frame$model <- "model"
  }
  if (is.null(frame$time)) {
    frame$time <- NA_real_
  }
  if (!is.null(set) && !is.null(frame$set)) {
    frame <- frame[frame$set == set, , drop = FALSE]
  }
  if (!is.null(models)) {
    frame <- frame[frame$model %in% models, , drop = FALSE]
  }
  frame
}

.surv_dca_table <- function(dca, set = NULL, models = NULL) {
  frame <- if (inherits(dca, "uniml_result")) {
    dca$dca
  } else if (is.data.frame(dca)) {
    dca
  } else {
    NULL
  }
  if (is.null(frame)) {
    return(NULL)
  }
  required <- c("threshold", "strategy", "net_benefit")
  if (!all(required %in% names(frame))) {
    .abort(sprintf(
      "A decision-curve table needs: %s.", paste(required, collapse = ", ")
    ))
  }
  if (is.null(frame$time)) {
    frame$time <- NA_real_
  }
  if (is.null(frame$model)) {
    frame$model <- "model"
  }
  if (!is.null(set) && !is.null(frame$set)) {
    frame <- frame[frame$set == set, , drop = FALSE]
  }
  if (!is.null(models)) {
    frame <- frame[frame$model %in% models, , drop = FALSE]
  }
  frame
}
