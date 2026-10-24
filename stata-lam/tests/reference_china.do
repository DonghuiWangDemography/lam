version 16.0
set more off
args root
if `"`root'"' == "" local root `"`c(pwd)'"'
adopath ++ `"`root'/stata-lam/tests"'
local datafile `"`root'/lam-main/lam-china.dta"'
capture confirm file `"`datafile'"'
if _rc local datafile `"`root'/lam-china.dta"'
use `"`datafile'"', clear
count
scalar reference_N = r(N)
generate long freq = 1
collapse (sum) freq, by(k q s)
quietly oprobit k i.s [fw=freq]
matrix ref_b0 = e(b)
local nt = colsof(ref_b0)-9
matrix ref_mu = ref_b0[1,1..`nt']
levelsof q, local(surveys)
local equations
foreach j of local surveys {
    quietly oprobit k i.s [fw=freq] if q==`j'
    quietly levelsof k if q==`j', local(categories)
    local nc = wordcount("`categories'")-1
    matrix this_b = e(b)
    local start = colsof(this_b)-`nc'+1
    matrix ref_tau`j' = this_b[1,`start'...]
    local eqnames
    local cnames
    forvalues h=1/`nc' {
        local equations `equations' (tau`j'_`h':)
        local eqnames `eqnames' tau`j'_`h'
        local cnames `cnames' _cons
    }
    matrix colnames ref_tau`j' = `cnames'
    matrix coleq ref_tau`j' = `eqnames'
}
matrix ref_start = ref_mu
foreach j of local surveys {
    matrix ref_start = ref_start, ref_tau`j'
}
ml model lf lam_reference (k:k=i.s, noconstant) `equations' [fw=freq]
ml init ref_start
ml maximize, difficult nolog
assert e(converged)==1
scalar reference_ll = e(ll)
matrix reference_b = e(b)
matrix reference_V = e(V)
display as result "Original likelihood reference fitted"
display as result "N = " %12.0f reference_N
display as result "Log likelihood = " %20.10f reference_ll
display as result "Free parameters = " e(rank)
estimates save `"`root'/stata-lam/tests/reference_china.ster"', replace
