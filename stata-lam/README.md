# LAM for Stata

`lam` estimates a shared latent attitude trend from ordinal responses collected by different survey series. It implements the constant-cutpoint model in Wang, Xie, and Huang’s 2023 paper, *Trend Analysis with Pooled Data from Different Survey Series. The Latent Attitude Method*, Sociological Methodology, [doi.org/10.1177/00811750231193641](https://doi.org/10.1177/00811750231193641).

One command replaces the survey-specific likelihood blocks in the original do-files.

```stata
lam response, survey(source) time(year)
lam_plot
```

Version 0.1.0 uses Stata ado and Mata source. It requires no third-party packages or compiled plugin. The source targets Stata 16 or later. Runtime validation uses StataNow 19.5 for macOS. Earlier versions and other operating systems have not been runtime-tested.

## Start here

Install directly from GitHub.

```stata
net install lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
help lam
```

To run a standalone example, change to a working directory where you want to save files and fetch the ancillary files. Save any unsaved data first because the example starts with `clear all`.

```stata
net get lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
do synthetic.do
```

The example generates data, fits LAM, checks category probabilities and frequency-count equivalence, saves estimates and result tables, and exports `results/synthetic_trend.png`. `net get` places ancillary files in the current directory using their leaf filenames.

For the original China example, download or clone the [LAM repository](https://github.com/DonghuiWangDemography/lam), then run the following from its `stata-lam` folder.

```stata
adopath ++ "."
do examples/china.do
```

This example first looks for `../lam-main/lam-china.dta`, then `../lam-china.dta`. It supports both data layouts, writes outputs under `results`, and leaves the source data unchanged. The empirical dataset is in the repository but is not included in `net install` or `net get`.

An existing local copy can also be installed. Replace the example path with your repository directory.

```stata
net install lam, from("/path/to/lam/stata-lam") replace
```

中文入门见 [quickstart-zh.md](docs/quickstart-zh.md)。Step-by-step testing instructions are available in [English](docs/testing-en.md) and [中文](docs/testing-zh.md). Paths under `examples/` and `docs/` refer to the source repository layout.

## Model and identification

For survey `q`, period `t`, and ordered response category `k`, LAM uses

```text
Pr(Y = k | q, t) = Phi(tau[q,k] - mu[t]) - Phi(tau[q,k-1] - mu[t])
```

The latent residual variance is fixed at one. Each survey has its own ordered cutpoints, held constant over time. All surveys observed in a period share that period’s latent mean. By default, the earliest observed period has mean zero. `baseperiod(#)` changes this reference period.

Surveys must form a connected overlap network. An indirect path through other surveys is sufficient. A common year shared by every survey is unnecessary. Disconnected groups cannot share a single identified latent scale, so the command rejects them.

These data checks do not establish sample representativeness or measurement consistency. Researchers must justify that the series measure the same construct and represent the same population in overlapping periods. Larger response codes must indicate the same attitude direction across all series. The software does not infer substantive direction from value labels.

## Inputs

```stata
lam response [if] [in] [fweight], survey(source) time(year) ///
    [baseperiod(#) level(#) iterate(#) nolog]
```

| Input | Meaning |
| --- | --- |
| `response` | Numeric ordinal outcome, ordered by its numeric values within each survey |
| `survey()` | Numeric or string survey-series identifier |
| `time()` | Numeric period variable, with no requirement for consecutive periods |
| `baseperiod()` | Observed period whose latent mean is fixed at zero |
| `level()` | Confidence level, defaulting to Stata’s current confidence level |
| `iterate()` | Maximum optimizer iterations, default 300 |
| `nolog` | Suppress iteration output |

Missing outcomes, survey identifiers, periods, and weights are excluded. Zero-frequency rows do not contribute. Each survey needs at least two observed categories in the estimation sample. Nonconsecutive numeric response codes are allowed. If a nominal category is absent from an entire survey in the estimation sample, the command cannot estimate its separate boundaries and uses the observed support.

Frequency counts can replace expanded respondent rows without changing the model likelihood or its observed-information standard errors.

```stata
lam category [fw=count], survey(source) time(year)
```

Frequency weights must be nonnegative integer counts. They are not sampling weights. This version does not implement `pweight`, `svy`, cluster-robust inference, covariates, or time-varying cutpoints. It also does not interpolate unobserved periods.

## Results and prediction

```stata
matrix list e(trend)
matrix list e(cutpoints)
predict double attitude, xb
predict double attitude_se, stdp
predict double probability, pr outcome(2)
lam_plot, title("Latent attitude trend")
graph export "trend.png", replace width(1800)
```

`e(trend)` contains period, estimate, standard error, and confidence limits. The reference-period row is zero by construction, including its standard error and confidence limits. These zeros describe a normalization, not a known population attitude.

`e(cutpoints)` identifies each cutpoint by its internal survey index and the original response code immediately below it. `e(survey1)`, `e(survey2)`, and so forth map internal indices to survey values. `e(categories)` records category codes and their ranks. `e(periods)` maps sorted period indices to original period values.

`e(b)` and `e(V)` contain free period means and cutpoints on their natural scale. Mean coefficients are named `mu:t1`, `mu:t2`, and so forth using sorted period indices, with the reference period omitted. Cutpoints are named `cut1:c1`, `cut1:c2`, and so forth. These stored results support standard commands such as `estimates store` and `lincom`.

`predict, xb` gives the fitted population mean for an observed period. It is not a respondent-specific posterior attitude estimate. `predict, stdp` gives the standard error of that mean. `predict, pr outcome(#)` returns the probability of the requested original response code for each record’s survey and period. An unknown period, survey, or category returns missing where required for that prediction. The response variable need not be present for prediction.

By default, confidence intervals and `e(V)` use the observed information matrix. If bootstrapping individual survey records, retain the survey-period sampling structure and verify that response support and coefficient dimensions remain constant across resamples. Resampling aggregated count rows is not equivalent to resampling respondents.

After the standard `bootstrap` prefix or a manual covariance replacement, `predict, stdp` and `lam_plot` use the current `e(V)`. The matrices `e(trend)` and `e(cutpoints)` are snapshots from the original fit and are not refreshed by an external prefix. `lam_plot` returns the current plot table in `r(trend)`. Under bootstrap, its intervals use a normal approximation with bootstrap standard errors and are not percentile or bias-corrected intervals.

`e(surveyvalues)` preserves numeric survey identifiers at full precision. String identifiers are stored literally in `e(survey#)` and used directly by prediction.

The estimator also rejects separated response patterns for which no finite maximum-likelihood estimate exists. This check is distinct from requiring a connected overlap network.

For the supplied China data, internal grouping reduces 108,036 record contributions to 483 positive survey-period-category cells while preserving the likelihood exactly. The optimizer uses an analytical gradient and Hessian through Mata's `d2` interface, reducing repeated evaluations of the same probabilities.

Validation details and release boundaries are in [release-notes.md](docs/release-notes.md). On macOS or Linux, run `sh tests/run_tests.sh` from any directory using the script's path. Set `STATA_BIN` if Stata is not found automatically. The scripts can also be run directly inside Stata as described in [tests/README.md](tests/README.md).

## Distribution files

| File | Role |
| --- | --- |
| `lam.ado` | Estimation command and Stata interface |
| `lam_mata.mata` | General likelihood and numerical implementation |
| `lam_p.ado` | Postestimation predictions |
| `lam_plot.ado` | Trend plot with confidence interval |
| `lam.sthlp` | Estimation and postestimation help |
| `lam_plot.sthlp` | Plotting help |
| `lam.pkg` and `stata.toc` | Stata installation metadata |
| `examples` | Reproducible empirical and synthetic examples |
| `tests` | Validation scripts and comparison with the original likelihood |
| `docs` | Method audit and Chinese quick start |

Both supplied examples have been run in StataNow 19.5. The synthetic example
passes its category-probability and frequency-input equivalence checks. A
separate installation check loads the command from an isolated ado directory,
fits a model, predicts probabilities, and exports a PNG graph.

The package implements the paper’s main constant-cutpoint specification. The original exploratory alternatives are not silently translated into the new estimator. In particular, unrestricted time-varying thresholds can be confounded with period means.

Maintainer [Donghui Wang](https://github.com/DonghuiWangDemography). Report reproducible problems through the [repository issues](https://github.com/DonghuiWangDemography/lam/issues). Version 0.1.0 does not include a software license.
