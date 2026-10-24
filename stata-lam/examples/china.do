version 16.0
clear all
set more off

* Run from the stata-lam directory in the workspace or GitHub checkout.
adopath ++ "."
capture mkdir results
local datafile "../lam-main/lam-china.dta"
capture confirm file "`datafile'"
if _rc local datafile "../lam-china.dta"
capture confirm file "`datafile'"
if _rc {
    display as error "Run from stata-lam with lam-china.dta in the parent folder or ../lam-main."
    exit 601
}

use "`datafile'", clear
lam k, survey(q) time(s) nolog
estimates save "results/china_lam.ster", replace
matrix list e(trend)
matrix list e(cutpoints)

predict double latent_mean, xb
predict double latent_mean_se, stdp
lam_plot, title("Americans' attitudes toward China") ///
    xtitle("Year") ytitle("Latent mean relative to 1974")
graph export "results/china_trend.png", replace width(1800)

matrix trend = e(trend)
matrix cutpoints = e(cutpoints)
preserve
clear
svmat double trend, names(col)
export delimited using "results/china_trend.csv", replace
save "results/china_trend.dta", replace
restore
preserve
clear
svmat double cutpoints, names(col)
export delimited using "results/china_cutpoints.csv", replace
restore

* These confidence intervals use the observed information matrix.
* Individual-record bootstrap is a separate analysis and is not run here.
* Resampling aggregated count rows would not resample respondents.
