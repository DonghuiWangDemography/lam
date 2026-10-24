version 16.0
clear all
set more off
args root
if `"`root'"' == "" local root `"`c(pwd)'"'
adopath ++ `"`root'/stata-lam"'
capture log close
log using `"`root'/stata-lam/tests/test_separation.log"', text replace

* A period with only the highest category has an unbounded latent mean
input byte y byte t byte q int n
0 1 1 50
1 1 1 50
1 2 1 100
end
capture noisily lam y [fw=n], survey(q) time(t) nolog
assert _rc==430

* The analogous lower-tail separation must also be rejected
replace y=0 if t==2
capture noisily lam y [fw=n], survey(q) time(t) nolog
assert _rc==430

* Joint drift of a period mean and a survey threshold is also unbounded
clear
input byte y byte t byte q int n
0 1 1 50
1 1 1 50
1 2 1 100
0 1 2 100
0 2 2 50
1 2 2 50
end
capture noisily lam y [fw=n], survey(q) time(t) nolog
assert _rc==430

* Adding the missing opposing response closes the separation direction
local last=_N+1
set obs `last'
replace y=0 in `last'
replace t=2 in `last'
replace q=1 in `last'
replace n=1 in `last'
quietly lam y [fw=n], survey(q) time(t) nolog
assert e(converged)==1

* Negative and zero period codes must be mapped without rank collisions
clear
input byte y byte t byte q int n
0 -1 1 50
1 -1 1 50
0 0 1 40
1 0 1 60
0 1 1 30
1 1 1 70
end
quietly lam y [fw=n], survey(q) time(t) nolog
matrix before_recode=e(b)
matrix periods_before=e(periods)
assert periods_before[1,1]==-1
assert periods_before[1,2]==0
assert periods_before[1,3]==1
replace t=t+11
quietly lam y [fw=n], survey(q) time(t) nolog
matrix after_recode=e(b)
mata st_numscalar("period_recode_difference",max(abs(st_matrix("before_recode")-st_matrix("after_recode"))))
assert period_recode_difference<1e-10
* A connected finite likelihood need not have N greater than parameter count
clear
input byte y byte t byte q
1 1 1
3 1 1
2 2 1
end
quietly lam y, survey(q) time(t) nolog
assert e(N)==e(k)
assert e(converged)==1
assert abs(e(ll)-3*ln(1/3))<1e-8

display as result "SEPARATION AND PERIOD RECODING TESTS PASSED"
log close
