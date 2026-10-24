# Validation

Run these do-files from the workspace or repository root that contains `stata-lam`. The China test looks for `lam-main/lam-china.dta` first, then `lam-china.dta` directly in that root.

```stata
do stata-lam/tests/test_lam.do
do stata-lam/tests/test_china.do
do stata-lam/tests/test_derivatives.do
do stata-lam/tests/test_separation.do
do stata-lam/tests/test_postestimation.do
```

An absolute workspace directory can also be supplied as the first argument to either do-file. The original files under `lam-main` are read but never modified.

`test_lam.do` checks agreement with Stata's ordered probit command, including its observed-information covariance matrix. It also checks synthetic parameter recovery, indirect overlap, reference-period invariance, arbitrary ordered numeric category codes, string survey identifiers, sample restrictions, unchanged input data, frequency-count equivalence, and rejection of unidentified or invalid input.

`test_china.do` re-estimates the original likelihood with its original natural-threshold parameterization and original separate ordered-probit initial values. The original evaluator is reproduced in `lam_reference.ado` with its program renamed so that it cannot shadow the package command. Exact frequency counts reduce the computation without altering the likelihood. All free coefficients, the full observed-information covariance matrix, and the maximized log likelihood are compared against the new implementation.

The original model fits 108,036 observations across 40 periods and 13 survey series. There are 39 free period means and 50 thresholds. Its independently re-estimated log likelihood is -155574.5182879295, which agrees with the paper's rounded -155574.5.

`test_derivatives.do` evaluates analytical first and second derivatives away from the optimum using central finite differences. It also verifies a far-tail probability against an independent R calculation with `pnorm` and log probabilities.

`test_separation.do` rejects unbounded period means and joint mean-threshold drift, accepts a nearby identified model, and verifies invariance to period recoding that includes negative and zero values.

`test_postestimation.do` compares predicted means, standard errors, and category probabilities with ordered probit. It checks estimates saved and restored after clearing Mata, predictions on new data, unknown categories and identifiers, empty prediction samples, double-precision numeric identifiers, and string identifiers containing spaces, quotes, and literal macro punctuation. It also verifies that prediction and plotting use the current covariance matrix after a manual repost or a stratified 20-replicate bootstrap. Plotting must preserve the data and estimates even when a graph option is invalid.

The initial validation completed on StataNow 19.5. The new implementation gave a log likelihood of -155574.5182877433. Maximum absolute differences from the original evaluator were 8.8317240e-06 for coefficients and 3.1426512e-07 for covariance entries. Maximum finite-difference discrepancies were 1.0456034e-08 for the gradient and 1.3113095e-08 for the Hessian. All five test suites passed.

The do-files stop at the first failed assertion. Successful runs end with `ALL CORE TESTS PASSED`, `CHINA REPLICATION PASSED`, `DERIVATIVE AND TAIL TESTS PASSED`, `SEPARATION AND PERIOD RECODING TESTS PASSED`, and `ALL POSTESTIMATION TESTS PASSED`. Inspect these markers in the text logs because Stata batch exit codes alone do not reliably indicate whether every assertion ran.

Numerical tolerances allow small differences between Stata optimizers and platforms. Agreement with the original likelihood validates implementation of the model. It does not establish sample representativeness, time-invariant measurement, or other substantive assumptions for a new dataset.
