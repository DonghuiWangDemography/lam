# Latent Attitude Method for Stata

Estimate a common latent attitude trend from survey series that use different ordinal response scales.

This repository contains the original replication materials for Wang, Xie, and Huang's 2023 paper and a general Stata implementation of its constant-cutpoint model. Version 0.1.0 replaces the manually specified survey likelihood blocks with one command.

```stata
lam response, survey(source) time(year)
lam_plot
```

[中文说明](README.zh-CN.md) · [Full package documentation](stata-lam/README.md) · [Testing guide](stata-lam/docs/testing-en.md)

## Install

Run the following in Stata.

```stata
net install lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
help lam
```

The package uses ado and Mata source and has no third-party dependencies. It declares Stata 16 compatibility and has been tested on StataNow 19.5 for macOS. Earlier Stata versions and other operating systems have not been runtime-tested.

## Try the package

After installation, fetch and run the standalone synthetic example from a working directory where you want its output saved. Save any unsaved data first because the example loads a new dataset.

```stata
net get lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
do synthetic.do
```

The example creates three overlapping survey series, fits the model, verifies category probabilities, checks frequency-count equivalence, and saves a PNG figure and result tables under `results`.

To reproduce the China application, download or clone this repository and substitute its local path below.

```stata
cd "/path/to/lam/stata-lam"
do examples/china.do
matrix list e(trend), format(%9.4f)
```

The example reads `lam-china.dta` at the repository root. It estimates 40 period means from 13 survey series and 108,036 observations. The log likelihood is approximately -155574.5183. Outputs include `results/china_trend.png`, `results/china_trend.csv`, `results/china_cutpoints.csv`, and `results/china_lam.ster`.

## Use your own data

```stata
lam response, survey(source) time(year) nolog
matrix list e(trend)
matrix list e(cutpoints)
predict double attitude, xb
predict double probability, pr outcome(2)
lam_plot
graph export "trend.png", replace width(1800)
```

For integer counts, use `lam response [fw=count], survey(source) time(year)`. The outcome must be numeric and ordered, the survey identifier may be numeric or string, and the period must be numeric. Nonconsecutive response codes and irregular periods are supported. Larger response values must have the same substantive direction across surveys.

Each survey has its own ordered cutpoints, held constant over time. Surveys in the same period share a latent population mean. The latent residual variance is fixed at one, and the earliest observed period is normalized to mean zero. `baseperiod(#)` selects a different reference period.

Survey series must be connected through overlapping periods, directly or through other surveys. They must measure the same construct and represent the same underlying population in those periods. Software checks cannot establish those substantive assumptions. The command checks graph connectivity and rejects separated response patterns without a finite maximum-likelihood estimate.

`predict, xb` returns a period's population mean on the latent scale, not a respondent-specific latent score. Default confidence intervals use the observed information matrix. Bootstrap inference is separate, and resampling aggregated count rows does not resample respondents. See the [full documentation](stata-lam/README.md) for inference and prediction details.

Version 0.1.0 implements the main constant-cutpoint specification. It does not implement covariates, sampling weights, complex-survey or cluster-robust inference, logistic links, time-varying cutpoints, or interpolation of missing periods.

## Verification and computation

The original and new estimators were compared using all 89 free parameters and their full covariance matrices. Maximum absolute differences were 0.000008832 for coefficients and 0.0000003143 for covariance entries. Five test suites cover original-model replication, agreement with official `oprobit`, analytical derivatives, identification and separation, frequency counts, prediction, saved estimates, plotting, and bootstrap compatibility.

The estimator aggregates the 108,036 China records into 483 occupied survey-period-response cells without changing their likelihood contribution. It also supplies analytical gradients and Hessians. These changes reduce repeated computation while retaining the same model.

Run the tests from the repository root.

```stata
cd "/path/to/lam"
do stata-lam/tests/test_lam.do
do stata-lam/tests/test_china.do
do stata-lam/tests/test_derivatives.do
do stata-lam/tests/test_separation.do
do stata-lam/tests/test_postestimation.do
```

Tests deliberately trigger some errors and check their return codes. Each suite must reach its final `PASSED` marker. See the [English testing guide](stata-lam/docs/testing-en.md), [中文测试教程](stata-lam/docs/testing-zh.md), and [validation details](stata-lam/tests/README.md).

## Repository contents

| Location | Purpose |
| --- | --- |
| `stata-lam/` | General Stata package, documentation, examples, and tests |
| `lam-ml.do` | Original China application code |
| `lam-sim.do` | Original simulation example |
| `lam-china.dta` | Original China application data |
| `2020-05-28-simulation-random.csv` | Original simulation input |

The original files remain as replication materials. The legacy do-files contain machine-specific paths. Use `stata-lam/examples` for the current runnable workflow. The installable package does not include the original empirical dataset.

## Reference and contact

Wang, Donghui, Yu Xie, and Junming Huang. 2023. *Trend Analysis with Pooled Data from Different Survey Series. The Latent Attitude Method.* Sociological Methodology. [DOI 10.1177/00811750231193641](https://doi.org/10.1177/00811750231193641).

Report problems through [GitHub Issues](https://github.com/DonghuiWangDemography/lam/issues). Include the Stata version, command, error message, and a small reproducible example when possible. Maintainer profile is [DonghuiWangDemography](https://github.com/DonghuiWangDemography).

See the [release notes](stata-lam/docs/release-notes.md) and [method audit](stata-lam/docs/method-audit.md). No explicit software license is included in this version.
