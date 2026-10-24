version 16.0
clear all
set more off
set seed 20261023

* Run from the source directory, or run the net-get copy after net install.
adopath ++ "."
capture mkdir results

* Three connected survey series with different ordinal response scales.
* Each retained survey-period contains 1,000 independent observations.
set obs 30000
generate byte source = ceil(_n / 10000)
generate int year = 2000 + floor(mod(_n - 1, 10000) / 1000)
generate byte period = year - 1999
keep if (source == 1 & inrange(period, 1, 5)) | ///
    (source == 2 & inrange(period, 4, 8)) | ///
    (source == 3 & inrange(period, 7, 10))
generate double true_mean = 0.5 * sin((period - 1) * 0.7) + 0.04 * (period - 1)
generate double latent = true_mean + rnormal()
generate byte response = .
replace response = 1 + (latent > -1.5) + (latent > 0) + (latent > 1.6) if source == 1
replace response = 1 + (latent > -2.2) + (latent > -0.8) + ///
    (latent > 0.9) + (latent > 1.4) if source == 2
replace response = 1 + (latent > -1.4) + (latent > -0.05) + (latent > 1.5) if source == 3

lam response, survey(source) time(year) nolog
assert e(converged) == 1
estimates save "results/synthetic_lam.ster", replace
predict double estimated_mean, xb

* Category probabilities sum to one on each survey's fitted support.
generate double probability_sum = 0
forvalues k = 1/5 {
    predict double p`k', pr outcome(`k')
    replace probability_sum = probability_sum + p`k' if !missing(p`k')
}
assert abs(probability_sum - 1) < 1e-10

lam_plot, title("Synthetic latent attitude trend") ///
    xtitle("Year") ytitle("Latent mean relative to 2000")
graph export "results/synthetic_trend.png", replace width(1800)

matrix trend = e(trend)
preserve
clear
svmat double trend, names(col)
generate double truth = 0.5 * sin((time - 2000) * 0.7) + 0.04 * (time - 2000)
generate double error = estimate - truth
export delimited using "results/synthetic_trend.csv", replace
save "results/synthetic_trend.dta", replace
restore

* Counts produce the same likelihood without expanding individuals.
scalar individual_ll = e(ll)
preserve
contract source year response, freq(count)
lam response [fw=count], survey(source) time(year) nolog
assert abs(e(ll) - individual_ll) < 1e-6
restore
