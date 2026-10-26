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

spec <- surv_spec(time = "time", status = "status", time_points = c(12, 36, 60))
spec <- surv_set_learners(spec, c("coxph", "rpart"), tune = FALSE)

set.seed(1)
index <- sample.int(nrow(data))
result <- run_surv_workflow(
  spec,
  train = data[index[1:280], ],
  validation = data[index[281:340], ],
  test = data[index[341:400], ]
)

result$metrics[, c("model", "set", "cindex", "ibs")]
result$ranking
```

## What it reports

| Statistic | Function | Definition | Checked against |
| --- | --- | --- | --- |
| Concordance index | `surv_cindex()` | Harrell's C over comparable pairs | `survival::concordance()` |
| Brier score at `t` | `surv_brier()` | Graf's IPCW Brier score | `riskRegression::Score()` |
| Integrated Brier score | `surv_ibs()` | Brier integrated over your time points | hand trapezoid |
| Time-dependent AUC | `surv_auc()` | cumulative/dynamic AUC, IPCW | `riskRegression::Score()`; the area under `surv_roc_curve()` |
| Calibration at `t` | `surv_calibration()` | predicted vs Kaplan-Meier, by prediction group | definition of ICI |
| Net benefit at `t` | `surv_dca()` | IPCW net benefit vs treat-all / treat-none | definition of net benefit |
| Confidence intervals | `surv_ci()` | percentile bootstrap over patients | the point estimate it brackets |

Every statistic is pinned by a test against an independent implementation, and
every one is documented with the paper it follows. See
`system.file("REFERENCES.md", package = "unimlsurv")`.

## Explaining a survival model

A survival prediction is not one number, so `surv_shap()` makes you say which
quantity you are explaining: the risk score the C-index ranks
(`time_point = NULL`), or the predicted survival probability at a horizon
(`time_point = 36`). The two are not interchangeable.

```r
explained <- surv_shap(result, validation, model = "coxph")
plot_surv_shap_beeswarm(explained)
plot_surv_shap_dependence(explained, "x1")

explained_36 <- surv_shap(result, validation, model = "coxph", time_point = 36)
plot_surv_shap_importance(explained_36)
```

Kernel SHAP is used because survival learners are not all tree ensembles; a tree
shortcut would silently change the estimator depending on the model.

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
* **Uno's IPCW C-index.** A different estimand from Harrell's C; reporting both
  under one name is how the literature gets confused. Only Harrell's C is here.
* **Subgroup analyses.** Which subgroups matter is a clinical decision, not a
  statistical default. You get a saved model and a scorer instead of invented
  groups.
* **Tuning on the Brier score.** Its inverse-probability weights depend on the
  censoring distribution of the training data, so tuning on it would compare
  numbers computed on different data. Tuning uses the concordance index.
* **`surv.rpart` and other risk-only learners** are supported, but they report a
  C-index and nothing else. The package says so in `result$meta$notes` instead of
  filling in a number it cannot compute.

## Installation

`unimlsurv` is on GitHub only, because `mlr3proba` is not on CRAN.

```r
# install.packages("remotes")
remotes::install_github("zhengyangszy/unimlsurv")
```

## Documentation

* `vignette("unimlsurv")` — a full walk-through on a simulated cohort.
* `vignette("unimlsurv-zh")` — the same in Chinese.
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
