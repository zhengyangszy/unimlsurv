# unimlsurv 0.2.2

Time-point recalibration coefficients, the survival counterpart of the binary
calibration slope and intercept that `uniml` reports.

* New `surv_calibration_slope()` fits
  `cloglog(P(T <= t | x)) = intercept + slope * cloglog(1 - S(t | x))` at one time
  point: a complementary log-log link, because that is the scale on which a
  proportional hazards model is additive (van Houwelingen & le Cessie 1990;
  Royston & Altman 2013). Perfectly calibrated predictions give `slope = 1` and
  `intercept = 0`. Patients censored before the time point carry no information
  and are dropped; the rest are weighted by the inverse Kaplan-Meier censoring
  survival, the same weighting `surv_brier()` uses.
* A slope below 1 means the predictions are too extreme, above 1 too timid, and a
  non-zero intercept means they are systematically too high or too low.
* Fewer than 20 patients, fewer than 5 events or 5 non-events, separation, or
  non-convergence return `NA` with a `message` instead of a number nobody should
  report. A constant prediction has no slope; only the intercept (calibration-in-
  the-large) is returned in that case.
* `calc_surv_metrics()` gains `calib_intercept` and `calib_slope` in `by_time` and
  `calib_intercept_<t>`/`calib_slope_<t>` in `row`, so `surv_ci()` bootstraps them
  and the workflow reports them next to the Brier score and the AUC.
* `surv_calibration_curves()$summary` now carries `ici`, `slope`, `intercept`,
  `n`, `events`, and `message` side by side, one row per time point.
* `run_surv_workflow()` keeps them in `result$calibration$summary`. That table is
  exported as `calibration/summary.csv`, which needs `uniml (>= 0.4.0)`; the
  floor is raised so the numbers cannot be silently dropped on export.
* `plot_surv_calibration()` gains `annotate` (default `TRUE`), which writes the
  ICI, the slope and the intercept of every curve into the caption. It also
  accepts the list returned by `surv_calibration_curves()` directly.

# unimlsurv 0.2.1

* `surv_spec()` and `run_surv_workflow()` ask uniml for `methods = "none"` rather
  than an empty method vector. An empty vector used to fall through uniml's
  selector machinery and abort with "`selections` must be a non-empty list" when
  the vignettes were built in a fresh session. The workflow sets it again at the
  point the engine is called, so a specification assembled elsewhere cannot
  reintroduce a selector. Requires uniml (>= 0.3.4).

# unimlsurv 0.2.0

Three things a survival model report needs and did not have: intervals, model
interpretation, and the full learner list.

## Confidence intervals

* New `surv_ci()` bootstraps the concordance index, the integrated Brier score,
  and every per-time Brier score and time-dependent AUC, resampling patients.
  Resamples that cannot estimate a time point return `NA` and are counted in
  `n_effective` rather than silently dropping out.
* `run_surv_workflow()` gains `ci_R` (default 200) and reports the intervals for
  the comparison set in `result$ci`. `ci_R = 0` skips them.

## SHAP

* New `surv_shap()` explains a fitted survival model with Kernel SHAP. The
  quantity being explained has to be named, because a survival prediction is not
  a single number: leave `time_point = NULL` to explain the risk score the
  concordance index ranks, or pass a time point to explain S(t).
* `plot_surv_shap_importance()`, `plot_surv_shap_beeswarm()`,
  `plot_surv_shap_dependence()`, `plot_surv_shap_waterfall()`, and
  `plot_surv_shap_force()` draw them on `uniml`'s SHAP figures.
* Kernel SHAP rather than a tree shortcut, because survival learners are not all
  tree ensembles and the estimator should not change with the learner.

## Every survival learner mlr3 has

* `surv_learner_registry()` now lists all 29 `surv.*` learners that `mlr3`,
  `mlr3proba`, and `mlr3extralearners` expose in the session, not just a curated
  handful. Twelve of them keep a short alias and a tuning search space; the rest
  are discovered at run time and fitted at their defaults unless a search space is
  registered with `surv_register_learner()`.
* `available` now reflects whether the learner's own package is installed, rather
  than whether an object could be constructed, so a missing optional package is
  reported when the learner is chosen instead of when it is trained.
* `distr` reports whether the learner can produce survival probabilities; a
  risk-only learner reports a concordance index and nothing else.

## Feature selection is off by default

* `surv_spec()` no longer inherits `uniml`'s `lasso` selection. That selector
  screens features against the outcome column, which for a survival model is the
  event indicator: it would ignore the time entirely. Every feature is used unless
  the user selects them.

# unimlsurv 0.1.0

First release. The survival stack is implemented end to end, on top of the
`uniml` layers.

## Task type

* `surv_spec()` builds a `uniml_spec` with `task_type = "survival"` and a
  `survival` block holding the time column, the status column, the event value,
  the reporting time points, and the censoring model.
* Loading the package registers `"survival"` in the `uniml` task-type registry,
  so `uniml::get_task_type("survival")` and `uniml::rank_models()` know what a
  survival outcome is.
* `as_surv_data()` splits features from the two outcome columns, drops rows with
  a missing time or status instead of imputing them, refuses non-positive follow-up
  times, and refuses character predictors rather than guessing whether they are
  ordered.
* `make_surv_task()` builds an `mlr3proba::TaskSurv`.

## Metrics

* `surv_cindex()` is Harrell's concordance index, pinned against
  `survival::concordance(reverse = TRUE)`.
* `surv_brier()` is Graf's inverse-probability-of-censoring weighted Brier score
  with a marginal Kaplan-Meier censoring model, pinned against
  `riskRegression::Score()` to machine precision.
* `surv_auc()` is the cumulative/dynamic time-dependent AUC with IPCW, pinned
  against `riskRegression::Score()` and against the area under the ROC curve that
  `surv_roc_curve()` draws.
* `surv_ibs()` integrates the Brier score over the requested time points.
* `calc_surv_metrics()` reports all of them and returns both a wide row and a
  per-time table.
* A learner that predicts only a risk score still gets a concordance index; the
  metrics that need survival probabilities are reported as unavailable, never
  filled in.

## Calibration and decisions

* `surv_calibration()` compares predicted survival with Kaplan-Meier survival
  inside prediction groups, and summarises with the integrated calibration index.
* `surv_calibration_curves()` does that at every requested time point.
* `surv_dca()` and `surv_dca_curves()` compute net benefit at a time point with
  IPCW weights, against the treat-all and treat-none strategies.

## Workflow

* `run_surv_workflow()` runs `uniml::workflow_engine()` with a survival adapter,
  so imputation, scaling, encoding, and tuning are fitted on the training data
  only and reused untouched everywhere else. The validation set selects; the test
  and external sets report.
* `surv_cross_validate()` scores out-of-fold predictions rather than averaging
  fold-wise C-indices.
* `predict_surv()` predicts new patients with the fitted preprocessor and learner.
* `surv_learner_registry()`, `resolve_surv_learners()`, and
  `surv_register_learner()` cover Cox, recursive partitioning, Kaplan-Meier, and
  the `mlr3extralearners` survival learners that are installed.

## Figures

* `plot_surv_km()`, `plot_surv_roc()`, `plot_surv_calibration()`,
  `plot_surv_dca()`, and `plot_surv_time_metrics()`, all on `uniml::theme_uniml()`.

## Simulation

* `sim_surv_data()` simulates a Weibull proportional-hazards cohort with
  exponential censoring, administrative censoring, and optional missingness, and
  stores the true linear predictor.
* `sim_surv_truth()` turns that into the true risk score and true survival
  probabilities at any time points, so every metric can be checked against the
  truth.

## Scope

Competing risks, Uno's IPCW C-index, SHAP for survival predictions, and tuning on
the Brier score are deliberately absent. Each needs a decision, not just code;
`inst/REFERENCES.md` and the vignettes say which decision and why.
