# unimlsurv

[English](README.md) | **简体中文**

[![R-CMD-check](https://github.com/zhengyangszy/unimlsurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhengyangszy/unimlsurv/actions/workflows/R-CMD-check.yaml)

用 `mlr3` 做生存分析，结构和 [`uniml`](https://github.com/zhengyangszy/uniml) 完全一致。

`unimlsurv` 是 `uniml` 的生存分析版本。分层（规格、数据、任务、学习器、调参、指标、
绘图、工作流）和对象都一样，会用一个就会用另一个。变的是结局，以及因此变化的
"数字代表什么"。

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

## 报告哪些指标

| 指标 | 函数 | 定义 | 对照验证 |
| --- | --- | --- | --- |
| C 指数 | `surv_cindex()` | Harrell C，只算可比对 | `survival::concordance()` |
| 时点 Brier | `surv_brier()` | Graf 的 IPCW Brier | `riskRegression::Score()` |
| 积分 Brier | `surv_ibs()` | 在你给的时间点上积分 | 手算梯形积分 |
| 时依 AUC | `surv_auc()` | 累积/动态 AUC，IPCW | `riskRegression::Score()`；以及 `surv_roc_curve()` 曲线下面积 |
| 时点校准 | `surv_calibration()` | 预测值 vs 分组内 Kaplan-Meier | ICI 定义 |
| 时点重新校准 | `surv_calibration_slope()` | cloglog 斜率与截距，IPCW 加权 | 真模型给出斜率 1、截距 0 |
| 时点净获益 | `surv_dca()` | IPCW 净获益，对照全干预/不干预 | 净获益定义 |
| 置信区间 | `surv_ci()` | 按病人重抽样的百分位 bootstrap | 它所包围的点估计 |

每个指标都有测试钉在独立实现上，每个都写明所依据的文献，见
`system.file("REFERENCES.md", package = "unimlsurv")`。

## 解释生存模型

生存预测不是一个数，所以 `surv_shap()` 要求你说清楚解释的是哪个量：C 指数排序的
**风险评分**（`time_point = NULL`），还是某个时点的**预测生存概率**
（`time_point = 36`）。两者不能互换。

```r
explained <- surv_shap(result, validation, model = "coxph")
plot_surv_shap_beeswarm(explained)
plot_surv_shap_dependence(explained, "x1")

explained_36 <- surv_shap(result, validation, model = "coxph", time_point = 36)
plot_surv_shap_importance(explained_36)
```

用 Kernel SHAP 而不是树模型捷径：生存学习器不都是树集成，捷径会让估计量随模型变化。

## 把模型交给别人

亚组分析刻意不在范围内：哪些亚组有意义是临床决定。`uniml::save_result()` 把拟合好的
学习器、预处理器、选中的特征和规格写进一个文件，拿到它的人可以对任何子集打分，并用
同一套指标函数。

```r
path <- uniml::save_result(result, "model.rds")
restored <- uniml::read_result(path)
scored <- predict_surv(restored, older_patients)
```

## mlr3 能支持的生存学习器都在

`surv_learner_registry()` 列出当前安装的 `mlr3` 系列包暴露的**全部 29 个** `surv.*`
学习器，连同短别名、能否给生存概率、有没有预置调参空间。其中 12 个有别名和调参空间，
其余用默认参数即可使用。


## 刻意没做的部分

* **竞争风险。** 死于其他原因属于另一类事件，不是删失。这改变估计量本身，不在
  本包范围。
* **Uno 的 IPCW C 指数。** 和 Harrell C 是两个不同的估计量，用同一个名字混着报，
  正是文献里最容易出错的地方。本包只报 Harrell C。
* **亚组分析。** 哪些亚组有意义是临床决定，不是统计默认值。包给你一个存好的模型和
  一个打分函数，而不是替你编组。
* **用 Brier 调参。** 逆概率权重依赖训练数据的删失分布，用它调参等于拿不同数据
  上的数互相比。调参用 C 指数。
* **`surv.rpart` 这类只给风险分数的学习器**可以跑，但只报 C 指数，其他指标留空，
  并在 `result$meta$notes` 里说明原因，而不是填一个算不出来的数。

## 安装

因为 `mlr3proba` 不在 CRAN 上，`unimlsurv` 只通过 GitHub 发布。

```r
# install.packages("remotes")
remotes::install_github("zhengyangszy/unimlsurv")
```

## 文档

* `vignette("unimlsurv")` —— 英文完整示例。
* `vignette("unimlsurv-zh")` —— 中文完整示例。
* `?surv_spec`、`?run_surv_workflow`，以及每个指标的帮助页。

## 和 uniml 的关系

`unimlsurv` 依赖 `uniml`，并在它的任务类型注册表里注册 `"survival"`。预处理、
插补、标准化、编码、调参、导出用的都是 `uniml` 的实现；工作流跑的是
`uniml::workflow_engine()` 加一个生存适配器。`unimlsurv` 新增的是生存统计量、
生存图，以及"预处理器必须放过两列结局"这条规则。

## 引用

```r
citation("unimlsurv")
```
