# unimlsurv

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

Every statistic is pinned by a test against an independent implementation, and
every one is documented with the paper it follows. See
`system.file("REFERENCES.md", package = "unimlsurv")`.

## What it deliberately does not do

* **Competing risks.** A death from another cause is an event of a different
  type, not a censoring. That changes the estimand, so it is out of scope.
* **Uno's IPCW C-index.** A different estimand from Harrell's C; reporting both
  under one name is how the literature gets confused. Only Harrell's C is here.
* **SHAP for survival predictions.** Shapley values for a function-valued
  prediction are a research topic.
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
