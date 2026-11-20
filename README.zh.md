# unimlsurv

[English](README.md) | **简体中文**

[![R-CMD-check](https://github.com/zhengyangszy/unimlsurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhengyangszy/unimlsurv/actions/workflows/R-CMD-check.yaml)

生存分析，结构与 [`uniml`](https://github.com/zhengyangszy/uniml) 完全一致，原生引擎（`survival`、`rpart`）、依赖全部来自 CRAN。

## 安装

先从 CRAN 装 `uniml`，再从 GitHub 装 `unimlsurv`：

```r
install.packages("uniml")                             # CRAN
# install.packages("remotes")
remotes::install_github("zhengyangszy/unimlsurv")     # GitHub
```

unimlsurv 目前发 GitHub（CRAN 提交审核中）；Imports 里的依赖全部可从 CRAN 安装。

`unimlsurv` 是 `uniml` 的生存分析版本。分层（规格、数据、任务、学习器、调参、指标、
绘图、工作流）和对象都一样，会用一个就会用另一个。变的是结局，以及因此变化的
数字含义。

```r
library(unimlsurv)


data <- sim_surv_data(n = 400, seed = 2026)

spec <- surv_spec(time = "time", status = "status", time_points = c(3, 6, 12))
spec <- surv_set_learners(spec, c("coxph", "kaplan"), tune = FALSE)

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

## 它报告什么

| 统计量 | 函数 | 定义 | 验证对照 |
| --- | --- | --- | --- |
| 一致性指数 | `surv_cindex()` | 可比较配对上的 Harrell C | `survival::concordance()` |
| Uno 一致性指数 | `surv_cindex_uno()` | 截断点处的 IPCW 加权 C，权重取自训练集 | `survAUC::UnoC()` |
| `t` 时点 Brier | `surv_brier()` | Graf 反删失概率加权 Brier | `riskRegression::Score()` |
| 积分 Brier | `surv_ibs()` | Brier 在你选时点上的积分 | 手工梯形 |
| 时点 AUC | `surv_auc()` | cumulative/dynamic AUC，IPCW | `riskRegression::Score()`；`surv_roc_curve()` 曲线下面积 |
| Uno 时点 AUC | `surv_auc_uno()` | case 权重取自训练集删失曲线 | `survAUC::AUC.uno()` |
| 积分 AUC | `surv_iauc_uno()` 及 `iauc` 列 | 你选时点上的梯形积分 | 手工梯形 |
| `t` 时点校准 | `surv_calibration()` | 预测 vs 组内 Kaplan-Meier；ICI、`citl`、`e50`、`e90` | Austin 等 2019/2020 的定义 |
| `t` 时点再校准 | `surv_calibration_slope()` | cloglog 斜率与截距，IPCW | 真值模型给 slope 1、intercept 0 |
| 分层 Kaplan-Meier | `surv_risk_groups()` + `plot_surv_km(logrank = TRUE)` | 风险三分位 + log-rank 检验 | `survival::survdiff()` |
| `t` 时点净获益 | `surv_dca()` | IPCW 净获益 vs 全治 / 全不治 | 净获益定义 |
| 置信区间 | `surv_ci()` | 按病人重抽样的百分位 bootstrap | 它所夹住的点估计 |

每个统计量都被测试钉在独立实现上，每个都有其所依循的论文（DOI 逐条过
Crossref 核实）。见 `system.file("REFERENCES.md", package = "unimlsurv")`。

## 解释生存模型

生存预测不止一个数，`surv_shap()` 强制你说清解释的是哪个量：C 指数排序的
风险分（`time_point = NULL`），或某时点的预测生存概率（`time_point = 6`）。
两者不可互换。

```r
explained <- surv_shap(result, validation, model = "coxph")
plot_surv_shap_beeswarm(explained)
plot_surv_shap_dependence(explained, "x1")

explained_6 <- surv_shap(result, validation, model = "coxph", time_point = 6)
plot_surv_shap_importance(explained_6)
```

用 Kernel SHAP 是因为生存学习器不全是树模型；树的捷径会随模型悄悄更换估计量。
从 0.3.0 起，解释按**原变量**给出——每个因子的 one-hot 列加总回原变量，与论文
里 SHAP 图的惯例一致。

## 把模型交出去

亚组分析刻意不在范围内：哪些亚组重要是临床判断。`uniml::save_result()` 把拟合
好的学习器、预处理器、入选特征与规格写进一个文件，拿到的人可以对任何子集打分，
并使用同一批指标函数。

```r
path <- uniml::save_result(result, "model.rds")
restored <- uniml::read_result(path)
scored <- predict_surv(restored, older_patients)
```

## 原生生存引擎

`surv_learner_registry()` 列出原生引擎——`coxph`（经 `survival::coxph` 的 Cox
比例风险）与 `kaplan`（无协变量 Kaplan-Meier 基线）——含各自的调参空间。高级
用户可用 `surv_register_learner()` 注册自己的 mlr3 生存学习器；该路径需要
GitHub 上的 `mlr3proba` 包。

## 刻意不做的

* **竞争风险。**死于其他原因是另一类型的事件，不是删失。这改变的是估计量本身，
  所以不在范围内。
* **NRI（净重分类改善）**。Steyerberg 框架把它留给"新标志物是否增值"的问题，
  不是单个模型的常规报告。
* **D-calibration 与 IPCW log-loss**。预测模型论文里罕见；审稿人若要求，可以借
  直接算。
* **分布组合**。给纯风险分模型借另一个学习器的基线"装"上完整生存分布；这会
  改变"该模型的 Brier"的含义，所以保持为可另行选择的设计决定。
  决定。
* **亚组分析。**哪些亚组重要是临床判断，不是统计默认。给你的是存档模型和打分器，
  不是替你发明的分组。
* **用 Brier 调参。**其反删失权重取决于训练数据的删失分布，用它调参会拿不同口径
  的数互相比。调参用 C 指数。
  C 指数。
* **注册的纯风险分学习器**只报 C 指数和 Uno C。包把这个写进
  `result$meta$notes`，而不是填一个算不出来的数。
## 安装

## 安装

```r
# install.packages("remotes")
remotes::install_github("zhengyangszy/unimlsurv")
```

## 文档

* `vignette("unimlsurv")` —— 模拟队列上的完整教程。
* `vignette("unimlsurv-zh")` —— 中文版。
* `vignette("unimlsurv-metrics")` —— 每个指标逐个讲透：回答什么问题、怎么读、
  什么时候别读。
* `vignette("unimlsurv-metrics-zh")` —— 中文版。
* `?surv_spec`、`?run_surv_workflow`，以及每个指标的帮助页。

## 与 uniml 的关系

`unimlsurv` 依赖 `uniml` 并在其任务类型注册表里注册 `"survival"`。预处理、插补、
标准化、编码、调参与导出都是 `uniml` 的实现；工作流跑的是
`uniml::workflow_engine()` 加一个生存适配器。`unimlsurv` 添加的是生存统计量与
图形，以及"预处理器必须放过两个结局列"的规则。

## 引用

```r
citation("unimlsurv")
```
