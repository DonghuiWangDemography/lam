# Version 0.1.0

23 October 2026. LAM is distributed from the [DonghuiWangDemography/lam repository](https://github.com/DonghuiWangDemography/lam) in its `stata-lam` directory.

```stata
net install lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
```

The package implements the paper's constant-cutpoint ordered-probit model as a general Stata estimation command. It accepts arbitrary survey identifiers, irregular numeric periods, ordered numeric response codes, and nonnegative integer frequency counts. Estimation uses exact grouped likelihood contributions, analytical derivatives, ordered cutpoint transformations, and model-based observed-information uncertainty.

The command checks both the survey-overlap network and existence of a finite maximum-likelihood estimate. It rejects disconnected or separated data before optimization. It preserves the input dataset and provides standard stored estimates, predictions, a trend plot, examples, and installation metadata.

The source declares Stata 16 compatibility and has been executed on StataNow 19.5 for macOS. Earlier Stata versions and other operating systems have not been runtime-tested. No compiled binaries or third-party Stata dependencies are required.

## Verification

- The original China example is reproduced with 108,036 observations and 89 free parameters. The original log likelihood is -155574.5182879295 and the new log likelihood is -155574.5182877433. Maximum absolute differences are 0.000008832 for coefficients and 0.0000003143 for covariance entries.
- The single-survey special case agrees with Stata's official ordered probit estimator, including its covariance matrix.
- Central finite differences verify the analytical gradient and Hessian away from the optimum. Stable log probabilities agree with an independent R calculation in the Gaussian tail.
- Frequency counts agree with expanded individual records. Tests also cover baseline changes, indirect overlap, single-period models, irregular category and period codes, missing values, separation, and unchanged input data.
- Predictions agree with ordered probit probabilities and sum to one on each survey's observed support. Stored-estimate reloads, unknown prediction inputs, exact survey identifiers, plotting, and a stratified bootstrap smoke test pass.
- Installation in an isolated ado directory and both example scripts have run successfully. The generated PNG figures have been visually inspected.

Detailed commands, tolerances, and output are recorded in the test scripts and logs. Successful execution checks explicit final markers, because the Stata batch process exit code alone can miss do-file errors.

## Distribution boundary

The source archive contains software, documentation, and test scripts. The original article and empirical data remain in their original folders and are not redistributed inside the installable package. The China replication test requires the supplied original data alongside the source tree.

The examples and China replication test support both the local workspace's `lam-main/lam-china.dta` and the original GitHub repository's root-level `lam-china.dta`. Both layouts have been tested and produce identical trend tables. Separate Chinese and English testing tutorials explain how to run the examples and verify all success markers.

Maintainer [Donghui Wang](https://github.com/DonghuiWangDemography). Version 0.1.0 does not include a software license. Author credits describe the statistical method and do not imply that every paper coauthor has reviewed this software implementation.

The initial implementation does not include covariates, sampling weights, complex-survey or cluster-robust inference, logistic links, formal overlap-equivalence tests, or the original exploratory time-varying-threshold specifications. Those extensions require separate statistical design and validation.
