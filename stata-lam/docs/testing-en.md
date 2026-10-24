# Testing LAM in Stata

Start with the supplied China example, then run the complete test suite. Download or clone the [LAM repository](https://github.com/DonghuiWangDemography/lam) and replace `/path/to/lam` with its local root directory, which contains `stata-lam` and the original data. Examples and test scripts execute `clear all`, so save any unsaved data in Stata before running them.

Open Stata's Do-file Editor, paste one complete code block, and click Do. Alternatively, enter the commands individually in the Command window. Do not copy the Markdown code fences.

## Run the original-data example

```stata
cd "/path/to/lam/stata-lam"
do examples/china.do
```

The example loads the original data, estimates LAM, predicts, plots the trend, and saves results. It does not modify the source data. It supports data stored at `lam-main/lam-china.dta` under the project directory or directly at `lam-china.dta` in the GitHub repository root.

A graph titled Americans' attitudes toward China should appear. The following files are written under `stata-lam/results`.

| File | Contents |
| --- | --- |
| `china_trend.png` | Latent trend and 95% confidence interval |
| `china_trend.csv` | Periods, estimates, standard errors, and confidence limits |
| `china_trend.dta` | The same trend table as a Stata dataset |
| `china_cutpoints.csv` | Survey cutpoints and standard errors |
| `china_lam.ster` | Saved Stata estimates |

Run these checks immediately after the example. Successful `assert` commands are silent. The last line prints only after all these checks pass.

```stata
assert e(converged) == 1
assert e(N) == 108036
assert e(n_surveys) == 13
assert e(n_periods) == 40
assert e(k) == 89
assert abs(e(ll) - (-155574.5182877433)) < 0.0001
display %18.4f e(ll)
display "CHINA QUICK CHECK PASSED"
```

The log likelihood should be approximately `-155574.5183`. The estimation header may use scientific notation. The explicit display command gives additional digits. Small numerical differences across platforms are allowed.

```stata
matrix list e(trend), format(%9.4f)
matrix list e(cutpoints), format(%9.4f)
ereturn list
```

The five trend columns are period, estimate, standard error, lower limit, and upper limit. The 1974 mean is fixed at zero to identify the scale. Its zero standard error and interval follow from that normalization. Other means measure latent attitude differences relative to 1974.

`e(b)` contains free coefficients and `e(V)` contains their covariance matrix. This example uses observed-information confidence intervals, not bootstrap intervals.

## Call the estimator directly

The following is a complete session supporting either original-data layout.

```stata
cd "/path/to/lam/stata-lam"
adopath ++ "."
local china "../lam-main/lam-china.dta"
capture confirm file "`china'"
if _rc local china "../lam-china.dta"
use "`china'", clear
lam k, survey(q) time(s) nolog
predict double fitted_mean, xb
predict double fitted_mean_se, stdp
predict double probability_category2, pr outcome(2)
lam_plot
```

The inputs are ordinal response `k`, survey identifier `q`, and period `s`. Replace those names for a new dataset. `fitted_mean` is the population latent mean for a period, not a respondent-specific latent score. `probability_category2` is the fitted probability of response code 2 for the record's survey and period.

Reload the data and rerun the entire block to repeat it. Repeating only `predict` with an existing output variable name gives a variable-already-defined error. That is not an estimation failure.

## Run a data-free synthetic example

```stata
cd "/path/to/lam/stata-lam"
do examples/synthetic.do
```

This generates 14,000 records across three connected survey series and ten periods. The log likelihood is approximately `-15409.43571`. The script checks that fitted category probabilities sum to one and that frequency-count input reproduces the individual-record likelihood.

This example has no separate overall pass marker. It should run through `end of do-file` without stopping at an `assertion is false` error. Main outputs are `results/synthetic_trend.png` and `results/synthetic_trend.csv`. The CSV includes known means in `truth` and estimation differences in `error`. This example generates a fresh reproducible dataset rather than replicating the original simulation's exact random sample.

## Run all five test suites

Change the working directory to the repository root, one level above `stata-lam`. This differs from the examples' working directory. Paste this entire block into one do-file and run it.

```stata
cd "/path/to/lam"
do stata-lam/tests/test_lam.do
do stata-lam/tests/test_china.do
do stata-lam/tests/test_derivatives.do
do stata-lam/tests/test_separation.do
do stata-lam/tests/test_postestimation.do
display "ALL FIVE LAM TEST SUITES COMPLETED"
```

Each suite must reach its final success marker.

| Test script | Final success marker |
| --- | --- |
| `test_lam.do` | `ALL CORE TESTS PASSED` |
| `test_china.do` | `CHINA REPLICATION PASSED` |
| `test_derivatives.do` | `DERIVATIVE AND TAIL TESTS PASSED` |
| `test_separation.do` | `SEPARATION AND PERIOD RECODING TESTS PASSED` |
| `test_postestimation.do` | `ALL POSTESTIMATION TESTS PASSED` |

Tests deliberately try disconnected surveys, invalid weights, separated responses, and invalid graph options. Error messages during these checks can be expected. The scripts capture the errors and verify their return codes. Success means all five final markers appear. If execution stops before a marker, inspect the last lines of that suite's log.

Logs are written beside the scripts under `stata-lam/tests`, using the same filename with a `.log` suffix. The China test independently re-estimates the original likelihood and compares coefficients, the full covariance matrix, and the maximized likelihood. Other suites test agreement with official `oprobit`, derivatives, count equivalence, overlap, identification, prediction, and plotting.

On macOS or Linux, the supplied terminal runner is another option.

```sh
cd "/path/to/lam"
sh stata-lam/tests/run_tests.sh
```

The final message should be `All available LAM checks passed.`. If the runner reports that China replication was skipped because its data are absent, that run does not establish that China replication passed. The runner checks explicit success markers because Stata batch process exit codes alone do not reliably establish that every assertion ran.

## Install for routine use

After the examples work, install from the local package directory.

```stata
net install lam, from("/path/to/lam/stata-lam") replace
which lam
help lam
```

`which lam` shows the file Stata is using. If an earlier version remains loaded after editing source, save your data and restart Stata before repeating the test.

Version 0.1.0 can also be installed directly from GitHub. The package is in the `stata-lam` directory, alongside the original paper's code and data at the repository root.

```stata
net install lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
```

This installs the command and help files. To fetch optional documentation and the standalone synthetic example, change to the directory where you want to save them and run the following commands.

```stata
net get lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
do synthetic.do
```

`net get` places ancillary files directly in the current directory, so the downloaded example is run with `do synthetic.do`. The original China dataset is not part of the installable package. Run full replication from a project or repository copy that includes those data.
