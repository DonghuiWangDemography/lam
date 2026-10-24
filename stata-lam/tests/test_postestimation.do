version 16.0
clear all
set more off
set seed 733105
args root
if `"`root'"' == "" local root `"`c(pwd)'"'
adopath ++ `"`root'/stata-lam"'
capture log close
log using `"`root'/stata-lam/tests/test_postestimation.log"', text replace

* Compare every single-survey prediction against official ordered probit
set obs 5000
generate int year=2000+mod(_n,4)
generate byte survey=1
generate double latent=.2*(year-2000)+rnormal()
generate byte response=1+(latent>-.8)+(latent>.3)+(latent>1.1)
quietly oprobit response ib2000.year
predict double official_xb, xb
predict double official_se, stdp
forvalues k=1/4 {
    predict double official_p`k', pr outcome(`k')
}
quietly lam response, survey(survey) time(year) nolog
predict double package_xb, xb
predict double package_se, stdp
assert abs(package_xb-official_xb)<2e-5
assert abs(package_se-official_se)<2e-5
forvalues k=1/4 {
    predict double package_p`k', pr outcome(`k')
    assert abs(package_p`k'-official_p`k')<2e-5
}
assert abs(package_p1+package_p2+package_p3+package_p4-1)<1e-12
predict double empty_selection if 0, pr outcome(1)
assert missing(empty_selection)
display as result "PASS predictions agree with official oprobit"

* Persistence requires no live Mata state or original estimation data
tempfile saved_model
estimates save `"`saved_model'"', replace
mata: mata clear
estimates clear
ereturn clear
estimates use `"`saved_model'"'
predict double restored_xb, xb
predict double restored_p, pr outcome(2)
assert restored_xb==package_xb
assert restored_p==package_p2
preserve
clear
input year survey
2000 1
2001 1
2002 99
9999 1
. 1
2003 .
end
predict double new_xb, xb
predict double new_se, stdp
predict double new_p, pr outcome(2)
predict double absent_outcome, pr outcome(99)
assert !missing(new_xb) in 1/3
assert !missing(new_xb) in 6
assert missing(new_xb) in 4/5
assert !missing(new_p) in 1/2
assert missing(new_p) in 3/6
assert missing(absent_outcome)
assert new_xb==0 in 1
assert new_se==0 in 1
restore
display as result "PASS saved estimates and new-data prediction behavior"

* Predictions and graphs must use the currently posted covariance matrix
matrix scaled_V=4*e(V)
program define lam_test_repost, eclass
    ereturn repost V=scaled_V
end
lam_test_repost
predict double doubled_se, stdp
assert abs(doubled_se-2*package_se)<1e-12
mata st_matrix("current_trend",lam_current_trend(95))
matrix original_trend=e(trend)
forvalues j=1/4 {
    assert abs(current_trend[`j',3]-2*original_trend[`j',3])<1e-12
}
quietly datasignature
local data_signature `"`r(datasignature)'"'
matrix before_plot_b=e(b)
lam_plot, name(lam_validation_plot, replace)
quietly datasignature
assert `"`r(datasignature)'"'==`"`data_signature'"'
matrix after_plot_b=e(b)
mata st_numscalar("plot_b_difference",max(abs(st_matrix("before_plot_b")-st_matrix("after_plot_b"))))
assert plot_b_difference==0
capture noisily lam_plot, invalidgraphoption
assert _rc!=0
quietly datasignature
assert `"`r(datasignature)'"'==`"`data_signature'"'
display as result "PASS current covariance and plotting preserve data and estimates"

* Different category supports and string survey IDs containing spaces and quotes
clear
set obs 6000
generate byte q=1+mod(_n,2)
generate int year=2000+mod(floor((_n-1)/2),3)
generate double latent=.15*(year-2000)+rnormal()
generate double response=1+(latent>-.7)+(latent>.6) if q==1
replace response=1+(latent>.2) if q==2
replace response=response*10-30
generate str40 source="alpha " + char(34) + "quoted" + char(34) + " survey" if q==1
replace source="beta survey" if q==2
quietly lam response, survey(source) time(year) nolog
predict double p_low, pr outcome(-20)
predict double p_mid, pr outcome(-10)
predict double p_high, pr outcome(0)
assert abs(p_low+p_mid+p_high-1)<1e-12 if q==1
assert abs(p_low+p_mid-1)<1e-12 if q==2
assert missing(p_high) if q==2
assert !missing(p_low)
generate double precise_source=100+q/13
quietly lam response, survey(precise_source) time(year) nolog
predict double precise_low, pr outcome(-20)
assert abs(precise_low-p_low)<1e-12
replace source="tick "+char(96)+"x"+char(39)+" "+char(36)+"z" if q==1
quietly lam response, survey(source) time(year) nolog
predict double literal_low, pr outcome(-20)
assert abs(literal_low-p_low)<1e-12
display as result "PASS mixed category supports and survey identifiers"

* Standard bootstrap prefix must retain metadata for predictions and plotting
quietly bootstrap _b, reps(20) seed(19035) strata(q year) nodots nowarn: ///
    lam response, survey(source) time(year) nolog
assert e(N_reps)==20
assert "`e(vce)'"=="bootstrap"
predict double bootstrap_se, stdp
predict double bootstrap_p, pr outcome(-20)
assert !missing(bootstrap_se)
assert !missing(bootstrap_p)
assert bootstrap_se==0 if year==2000
matrix bootstrap_V=e(V)
assert abs(bootstrap_se-sqrt(bootstrap_V[1,1]))<1e-12 if year==2001
assert abs(bootstrap_se-sqrt(bootstrap_V[2,2]))<1e-12 if year==2002
mata st_matrix("bootstrap_trend",lam_current_trend(95))
assert abs(bootstrap_trend[2,3]-sqrt(bootstrap_V[1,1]))<1e-12
assert abs(bootstrap_trend[3,3]-sqrt(bootstrap_V[2,2]))<1e-12
quietly datasignature
local bootstrap_signature `"`r(datasignature)'"'
lam_plot, name(lam_validation_bootstrap, replace)
quietly datasignature
assert `"`r(datasignature)'"'==`"`bootstrap_signature'"'
display as result "PASS bootstrap prefix predictions and plotting"
display as result "ALL POSTESTIMATION TESTS PASSED"
log close
