# unimlsurv

**English** | [简体中文](README.zh.md)

[![R-CMD-check](https://github.com/zhengyangszy/unimlsurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhengyangszy/unimlsurv/actions/workflows/R-CMD-check.yaml)

Survival analysis with `mlr3`, in the same shape as [`uniml`](https://github.com/zhengyangszy/uniml).

`unimlsurv` is the time-to-event sibling of `uniml`. It keeps the layering —
specification, data, task, learners, tuning, metrics, plots, workflow — and the
same objects, so a user who knows one knows the other. What changes is the
outcome, and therefore what the numbers mean.

```r
library(unimlsurv)
library(mlr3)

data <- sim_surv_data(n = 400, seed = 2026)

spec <- surv_spec(time = "time", status = "status", time_points = c(3, 6, 12))
spec <- surv_set_learners(spec, c("coxph", "rpart"), tune = FALSE)

set.seed(1)
index <- sample.int(nrow(data))
result <- run_surv_workflow(
  spec,
  train = data[index[1:280], ],
  validation = data[index[281:340], ],
  test = data[index[341:400], ]
)

result$metrics[, c("model", "set", "cindex", "uno_c", "ibs", "iauc")]
result$ranking
```

## What it reports

| Statistic | Function | Definition | Checked against |
| --- | --- | --- | --- |
| Concordance index | `surv_cindex()` | Harrell's C over comparable pairs | `survival::concordance()` |
| Uno's concordance index | `surv_cindex_uno()` | IPCW-weighted C at a horizon, weights from the training data | `survAUC::UnoC()` |
| Brier score at `t` | `surv_brier()` | Graf's IPCW Brier score | `riskRegression::Score()` |
| Integrated Brier score | `surv_ibs()` | Brier integrated over your time points | hand trapezoid |
| Time-dependent AUC | `surv_auc()` | cumulative/dynamic AUC, IPCW | `riskRegression::Score()`; the area under `surv_roc_curve()` |
| Uno's AUC at `t` | `surv_auc_uno()` | IPCW AUC with training-based case weights | `survAUC::AUC.uno()` |
| Integrated AUCs | `surv_iauc_uno()` and the `iauc` column | trapezoid over your time points | hand trapezoid |
| Calibration at `t` | `surv_calibration()` | predicted vs Kaplan-Meier, by prediction group; ICI, `citl`, `e50`, `e90` | Austin et al. 2019/2020 definitions |
| Recalibration at `t` | `surv_calibration_slope()` | cloglog slope and intercept, IPCW | a true model gives slope 1, intercept 0 |
| Stratified Kaplan-Meier | `surv_risk_groups()` + `plot_surv_km(logrank = TRUE)` | risk terciles with the log-rank test | `survival::survdiff()` |
| Net benefit at `t` | `surv_dca()` | IPCW net benefit vs treat-all / treat-none | definition of net benefit |
| Confidence intervals | `surv_ci()` | percentile bootstrap over patients | the point estimate it brackets |

Every statistic is pinned by a test against an independent implementation, and
every one is documented with the paper it follows (DOIs verified against
Crossref). See `system.file("REFERENCES.md", package = "unimlsurv")`.

## Explaining a survival model

A survival prediction is not one number, so `surv_shap()` makes you say which
quantity you are explaining: the risk score the C-index ranks
(`time_point = NULL`), or the predicted survival probability at a horizon
(`time_point = 6`). The two are not interchangeable.

```r
explained <- surv_shap(result, validation, model = "coxph")
plot_surv_shap_beeswarm(explained)
plot_surv_shap_dependence(explained, "x1")

explained_6 <- surv_shap(result, validation, model = "coxph", time_point = 6)
plot_surv_shap_importance(explained_6)
```

Kernel SHAP is used because survival learners are not all tree ensembles; a tree
shortcut would silently change the estimator depending on the model. Since
0.3.0, explanations are reported per **variable** — the one-hot columns of each
factor are summed back, the same convention published SHAP figures use.

## Handing the model to someone else

Subgroup analyses are deliberately out of scope: which subgroups matter is a
clinical decision. `uniml::save_result()` writes the fitted learners, the fitted
preprocessor, the selected features, and the specification to one file, so whoever
opens it can score any subset and use the same metric functions.

```r
path <- uniml::save_result(result, "model.rds")
restored <- uniml::read_result(path)
scored <- predict_surv(restored, older_patients)
```

## Every survival learner mlr3 has

`surv_learner_registry()` lists all 29 `surv.*` learners the installed `mlr3`
packages expose, with a short alias, whether it can produce survival
probabilities, and whether it carries a curated tuning space. Twelve have a
curated alias and search space; the rest are usable at their defaults.


## What it deliberately does not do

* **Competing risks.** A death from another cause is an event of a different
  type, not a censoring. That changes the estimand, so it is out of scope.
* **Net reclassification improvement.** NRI is reserved by the Steyerberg
  framework for the "does a new marker add value" question, not for reporting
  one model.
* **D-calibration and IPCW log-loss.** Rarely reported in prediction papers;
  they can be borrowed from `mlr3proba`'s measure set if a reviewer asks.
* **Distribution composition.** `mlr3proba`'s `PipeOpDistrCompositor` can bolt a
  full survival distribution onto a risk-only model; that changes what "the
  model's Brier score" means, so it stays an opt-in design decision.
* **Subgroup analyses.** Which subgroups matter is a clinical decision, not a
  statistical default. You get a saved model and a scorer instead of invented
  groups.
* **Tuning on the Brier score.** Its inverse-probability weights depend on the
  censoring distribution of the training data, so tuning on it would compare
  numbers computed on different data. Tuning uses the concordance index.
* **`surv.rpart` and other risk-only learners** are supported, but they report a
  C-index and Uno's C and nothing else. The package says so in
  `result$meta$notes` instead of filling in a number it cannot compute.

## Installation

`unimlsurv` is on GitHub only, because `mlr3proba` is not on CRAN.

```r
# install.packages("remotes")
remotes::install_github("zhengyangszy/unimlsurv")
```

## Documentation

* `vignette("unimlsurv")` — a full walk-through on a simulated cohort.
* `vignette("unimlsurv-zh")` — the same in Chinese.
* `vignette("unimlsurv-metrics")` — every metric explained: what question it
  answers, how to read it, when not to read it.
* `vignette("unimlsurv-metrics-zh")` — the same in Chinese.
* `?surv_spec`, `?run_surv_workflow`, and the help page of every metric.

## Relationship to uniml

`unimlsurv` imports `uniml` and registers `"survival"` in its task-type
registry. Preprocessing, imputation, scaling, encoding, tuning, and export are
the `uniml` implementations; the workflow runs `uniml::workflow_engine()` with a
survival adapter. What `unimlsurv` adds is the survival statistics and figures,
and the rule that the preprocessor must leave the two outcome columns alone.

## Citation

```r
citation("unimlsurv")
```
