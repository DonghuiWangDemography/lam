version 16.0
clear all
set more off
set seed 417031
args root
if `"`root'"' == "" local root `"`c(pwd)'"'
adopath ++ `"`root'/stata-lam"'
capture log close
log using `"`root'/stata-lam/tests/test_lam.log"', text replace

* A single survey must reproduce official ordered probit
set obs 8000
generate int year = 2000 + mod(_n,4)
generate byte survey = 1
generate double latent = .25*(year-2000) + rnormal()
generate byte response = 1 + (latent > -.8) + (latent > .2) + (latent > 1.3)
quietly oprobit response ib2000.year
scalar ll_oprobit = e(ll)
matrix b_oprobit = e(b)
matrix V_oprobit = e(V)
quietly lam response, survey(survey) time(year) nolog
assert e(converged)==1
assert abs(e(ll)-ll_oprobit)<1e-6
matrix b_lam = e(b)
matrix V_lam = e(V)
mata st_numscalar("b_difference",max(abs(st_matrix("b_lam")-st_matrix("b_oprobit")[,2..cols(st_matrix("b_oprobit"))])))
mata st_numscalar("V_difference",max(abs(st_matrix("V_lam")-st_matrix("V_oprobit")[2..rows(st_matrix("V_oprobit")),2..cols(st_matrix("V_oprobit"))])))
assert b_difference < 2e-5
assert V_difference < 2e-5
display as result "PASS single survey agrees with official oprobit"

* With one period the cutpoints equal empirical cumulative normal quantiles
preserve
replace year=2000
quietly lam response, survey(survey) time(year) nolog
matrix one_period_cuts=e(cutpoints)
assert e(n_periods)==1
assert colsof(e(b))==3
forvalues j=1/3 {
    quietly count if response<=`j'
    scalar expected_cut=invnormal(r(N)/_N)
    assert abs(one_period_cuts[`j',3]-expected_cut)<1e-6
}
restore
display as result "PASS one-period cutpoint identification"

* Indirect overlap must be sufficient for a common scale
clear
set obs 15000
generate byte survey = ceil(_n/5000)
generate int year = 2000 + 2*(survey-1) + mod(_n,3)
generate double latent = -.15*(year-2000) + rnormal()
generate byte response = 1 + (latent > -.8) + (latent > .4) if survey==1
replace response = 1 + (latent > -1.2) + (latent > -.3) + (latent > .9) if survey==2
replace response = 1 + (latent > -.6) if survey==3
generate long rowid = _n
quietly datasignature
local original_signature `"`r(datasignature)'"'
quietly lam response, survey(survey) time(year) nolog
assert e(converged)==1
assert e(N)==15000
assert e(sample)==1
matrix b_raw = e(b)
matrix V_raw = e(V)
matrix trend_raw = e(trend)
matrix cuts_raw = e(cutpoints)
scalar ll_raw = e(ll)
quietly datasignature
assert `"`r(datasignature)'"' == `"`original_signature'"'
assert rowid==_n
assert trend_raw[1,2]==0
assert trend_raw[1,3]==0
matrix expected_cuts = (-.8\.4\-1.2\-.3\.9\-.6)
mata st_numscalar("recovery_mu_difference",max(abs(st_matrix("trend_raw")[,2]:+.15*(st_matrix("trend_raw")[,1]:-2000))))
mata st_numscalar("recovery_cut_difference",max(abs(st_matrix("cuts_raw")[,3]-st_matrix("expected_cuts"))))
assert recovery_mu_difference<.12
assert recovery_cut_difference<.12
display as result "PASS indirect overlap and data preservation"
display as result "PASS synthetic latent mean and threshold recovery"

* Changing the reference period must preserve fit and translate the latent scale
quietly lam response, survey(survey) time(year) baseperiod(2003) nolog
matrix trend_rebased = e(trend)
assert abs(e(ll)-ll_raw)<1e-6
mata st_numscalar("rebase_difference",max(abs(st_matrix("trend_rebased")[,2]-(st_matrix("trend_raw")[,2]:-st_matrix("trend_raw")[4,2]))))
assert rebase_difference<2e-5
assert trend_rebased[4,2]==0
assert abs(trend_rebased[1,3]-sqrt(V_raw[3,3]))<2e-5
forvalues j=2/7 {
    local p=`j'-1
    scalar expected_se = sqrt(V_raw[`p',`p']+V_raw[3,3]-2*V_raw[`p',3])
    assert abs(trend_rebased[`j',3]-expected_se)<2e-5
}
display as result "PASS reference period invariance"

* String survey identifiers and arbitrary ordered response codes
generate str8 source = "group" + string(survey)
generate double response_codes = response*10-30
quietly lam response_codes, survey(source) time(year) nolog
matrix b_recoded = e(b)
mata st_numscalar("recode_difference",max(abs(st_matrix("b_recoded")-st_matrix("b_raw"))))
assert abs(e(ll)-ll_raw)<1e-6
assert recode_difference<2e-5
display as result "PASS string sources and ordered response recoding"

* Missing observations and if restrictions must produce a reliable e(sample)
preserve
replace response = . in 1
replace year = . in 2
replace survey = . in 3
quietly lam response if rowid<=14999, survey(survey) time(year) nolog
assert e(sample)==0 in 1/3
assert e(sample)==0 in 15000
assert e(N)==14996
restore
display as result "PASS missingness and sample restrictions"

* Frequency counts must equal expanded observations
generate long freq=1
collapse (sum) freq, by(survey year response)
quietly lam response [fw=freq], survey(survey) time(year) nolog
matrix b_freq=e(b)
matrix V_freq=e(V)
mata st_numscalar("frequency_b_difference",max(abs(st_matrix("b_freq")-st_matrix("b_raw"))))
mata st_numscalar("frequency_V_difference",max(abs(st_matrix("V_freq")-st_matrix("V_raw"))))
assert abs(e(ll)-ll_raw)<1e-6
assert e(N)==15000
assert frequency_b_difference<2e-5
assert frequency_V_difference<2e-5
display as result "PASS frequency-count equivalence"

* Invalid estimation requests must not silently fit a different model
capture noisily lam response [fw=freq], survey(survey) time(year) baseperiod(1900) nolog
assert _rc!=0
preserve
replace year=year+100 if survey==3
capture noisily lam response [fw=freq], survey(survey) time(year) nolog
assert _rc!=0
restore
preserve
replace response=1 if survey==2
capture noisily lam response [fw=freq], survey(survey) time(year) nolog
assert _rc!=0
restore
preserve
replace freq=1.5 in 1
capture noisily lam response [fw=freq], survey(survey) time(year) nolog
assert _rc!=0
restore
capture noisily lam response [pw=freq], survey(survey) time(year) nolog
assert _rc!=0
capture noisily lam response [aw=freq], survey(survey) time(year) nolog
assert _rc!=0
display as result "PASS invalid inputs are rejected"
display as result "ALL CORE TESTS PASSED"
log close
