version 16.0
clear all
set more off
args root
if `"`root'"' == "" local root `"`c(pwd)'"'
adopath ++ `"`root'/stata-lam"'
capture log close
log using `"`root'/stata-lam/tests/test_china.log"', text replace
do `"`root'/stata-lam/tests/reference_china.do"' `"`root'"'
matrix reference_free_b = reference_b[1,2...]
matrix reference_free_V = reference_V[2...,2...]
local datafile `"`root'/lam-main/lam-china.dta"'
capture confirm file `"`datafile'"'
if _rc local datafile `"`root'/lam-china.dta"'
use `"`datafile'"', clear
quietly lam k, survey(q) time(s) nolog
assert e(converged)==1
assert e(N)==reference_N
assert e(rank)==89
matrix package_b=e(b)
matrix package_V=e(V)
scalar ll_difference=abs(e(ll)-reference_ll)
mata st_numscalar("maximum_b_difference",max(abs(st_matrix("package_b")-st_matrix("reference_free_b"))))
mata st_numscalar("maximum_V_difference",max(abs(st_matrix("package_V")-st_matrix("reference_free_V"))))
assert ll_difference<1e-5
assert maximum_b_difference<2e-4
assert maximum_V_difference<2e-4
display as result "PASS complete China likelihood, parameter and covariance replication"
display as result "Original log likelihood = " %20.10f reference_ll
display as result "Package log likelihood = " %20.10f e(ll)
display as result "Maximum coefficient difference = " %14.10g maximum_b_difference
display as result "Maximum covariance difference = " %14.10g maximum_V_difference
display as result "CHINA REPLICATION PASSED"
log close
