# References

Every statistic in `unimlsurv` follows a named definition. This file records
which one, so a reported number can be traced back to a paper.

## Concordance index

Harrell FE, Lee KL, Mark DB (1996). Multivariable prognostic models: issues in
developing models, evaluating assumptions and adequacy, and measuring and
reducing errors. *Statistics in Medicine*, 15(4), 361-387.

Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011). On the C-statistics for
evaluating overall adequacy of risk prediction procedures with censored survival
data. *Statistics in Medicine*, 30(10), 1105-1117. doi:10.1002/sim.4154

*What the package does*: `surv_cindex()` is Harrell's C, computed over comparable
pairs only, with a higher risk score meaning a shorter expected survival. It is
pinned by a test against `survival::concordance(reverse = TRUE)`.
Uno's IPCW C is not implemented; it reweights the same pairs by the inverse
censoring survival and answers a different question.

## Brier score and integrated Brier score

Brier GW (1950). Verification of forecasts expressed in terms of probability.
*Monthly Weather Review*, 78(1), 1-3.

Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). Assessment and comparison of
prognostic classification schemes for survival data. *Statistics in Medicine*,
18(17-18), 2529-2545.

Gerds TA, Schumacher M (2006). Consistent estimation of the expected Brier score
in general survival models with right-censored event times. *Biometrical
Journal*, 48(6), 1029-1040. doi:10.1002/bimj.200610301

*What the package does*: `surv_brier()` is Graf's inverse-probability-of-censoring
weighted Brier score with a marginal Kaplan-Meier censoring model, which is
`pec::pec(cens.model = "marginal")` and `riskRegression::Score()`'s default. The
test compares the two to machine precision. `surv_ibs()` integrates the Brier
score over the time points the user asked for, with the trapezoidal rule, and is
documented as such because it therefore depends on that grid.

## Time-dependent AUC

Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC curves.
*Biometrics*, 61(1), 92-105. doi:10.1111/j.0006-341X.2005.030814.x

Uno H, Cai T, Tian L, Wei LJ (2007). Evaluating prediction rules for t-year
survivors with censored regression models. *Journal of the American Statistical
Association*, 102(478), 527-537. doi:10.1198/016214507000000149

Blanche P, Dartigues JF, Jacqmin-Gadda H (2013). Estimating and comparing
time-dependent areas under receiver operating characteristic curves for censored
event times with competing risks. *Statistics in Medicine*, 32(30), 5381-5397.
doi:10.1002/sim.5958

*What the package does*: `surv_auc()` is the cumulative/dynamic AUC with inverse
censoring weights: a case is an event at or before `t`, a control is a patient
still event-free after `t`. It is pinned twice — against
`riskRegression::Score()` and against the area under the ROC curve that
`surv_roc_curve()` draws, so the number and the figure cannot drift apart.

## Calibration

Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW (2019).
Calibration: the Achilles heel of predictive analytics. *BMC Medicine*, 17, 230.
doi:10.1186/s12916-019-1466-7

Austin PC, Steyerberg EW (2019). The integrated calibration index (ICI) and
related metrics for assessing the calibration of survival models with binary
outcomes. *Statistics in Medicine*, 38(28), 5476-5494. doi:10.1002/sim.8356

*What the package does*: `surv_calibration()` groups patients by predicted
survival and compares the group mean prediction with the group Kaplan-Meier
estimate at the time point. The summary is the ICI, the patient-weighted mean
absolute difference. Groups are formed from the predictions, never from the
outcome, so calibration cannot peek at what it is being judged on.

## Decision curves

Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel method for
evaluating prediction models. *Medical Decision Making*, 26(6), 565-574.
doi:10.1177/0272989X06295361

Vickers AJ, Cronin AM, Elkin EB, Gonen M (2008). Extensions to decision curve
analysis, a novel method for evaluating diagnostic tests, prediction models and
molecular markers. *BMC Medical Informatics and Decision Making*, 8, 53.
doi:10.1186/1472-6947-8-53

Vickers AJ, van Calster B, Steyerberg EW (2016). Net benefit approaches to the
evaluation of prediction models, molecular markers, and diagnostic tests. *BMJ*,
352, i6. doi:10.1136/bmj.i6

*What the package does*: `surv_dca()` computes net benefit at one time point with
inverse-probability-of-censoring weights, against the treat-all and treat-none
strategies. A patient censored before the time point is not counted as a certain
non-event.

## Survival estimation and simulation

Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete
observations. *Journal of the American Statistical Association*, 53(282),
457-481. doi:10.1080/01621459.1958.10501452

Cox DR (1972). Regression models and life-tables. *Journal of the Royal
Statistical Society: Series B*, 34(2), 187-220.

Bender R, Augustin T, Blettner M (2005). Generating survival times to simulate
Cox proportional hazards models. *Statistics in Medicine*, 24(11), 1713-1723.
doi:10.1002/sim.2059

Burton A, Altman DG, Royston P, Holder RL (2006). The design of simulation
studies in medical statistics. *Statistics in Medicine*, 25(24), 4279-4292.
doi:10.1002/sim.2328

*What the package does*: `plot_surv_km()` draws Kaplan-Meier curves.
`sim_surv_data()` simulates a Weibull proportional-hazards cohort with
exponential censoring and stores the true linear predictor, so every metric can
be checked against the truth rather than against another implementation of the
same assumption.

## Confidence intervals

Efron B, Tibshirani RJ (1993). *An Introduction to the Bootstrap*. Chapman &
Hall. doi:10.1201/9780429246593

*What the package does*: `surv_ci()` resamples patients with replacement and
reports percentile intervals for the concordance index, the integrated Brier
score, and every per-time Brier score and AUC. A resample that cannot estimate a
time point (no cases, no controls, or no events) contributes `NA` and is counted
in `n_effective`, so an interval is never quietly built from a handful of
resamples.

## Model explanation

Lundberg SM, Lee SI (2017). A unified approach to interpreting model predictions.
*Advances in Neural Information Processing Systems*, 30.

Covert I, Lee SI (2021). Improving KernelSHAP: practical Shapley value estimation
using linear regression. *Proceedings of AISTATS*, PMLR 130, 3457-3465.

*What the package does*: `surv_shap()` computes Kernel SHAP values for a fitted
survival model. Because a survival prediction is a function, the quantity being
explained has to be named: the risk score when `time_point` is `NULL`, and the
predicted survival probability at that horizon otherwise. Kernel SHAP is used
rather than a tree-specific shortcut so the estimator does not change with the
learner. Explanations are properties of the model, not causal effects; the
vignette says so in the same words a paper should.

## Validation and reporting

Steyerberg EW (2019). *Clinical Prediction Models: A Practical Approach to
Development, Validation, and Updating*, 2nd edition. Springer.
doi:10.1007/978-3-030-16399-0

Collins GS, Moons KGM, Dhiman P, et al. (2024). TRIPOD+AI statement: updated
guidance for reporting clinical prediction models that use regression or machine
learning methods. *BMJ*, 385, e078378. doi:10.1136/bmj-2023-078378

Varma S, Simon R (2006). Bias in error estimation when using cross-validation for
model selection. *BMC Bioinformatics*, 7, 91. doi:10.1186/1471-2105-7-91

*What the package does*: the validation set selects, the test and external sets
report. Cross-validation is scored on out-of-fold predictions rather than by
averaging fold-wise C-indices.

## Software

Sonabend R, Kirchner P, Bender A, et al. (2021). mlr3proba: an R package for
machine learning in survival analysis. *Bioinformatics*, 37(17), 2789-2791.
doi:10.1093/bioinformatics/btab039

Bender A, Groll A, Scheipl F (2018). A general machine learning framework for
survival analysis. arXiv:1806.06108.

Lang M, Binder M, Richter J, et al. (2019). mlr3: a modern object-oriented machine
learning framework in R. *Journal of Open Source Software*, 4(44), 1903.
doi:10.21105/joss.01903
