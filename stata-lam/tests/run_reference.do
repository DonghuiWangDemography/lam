version 16.0
clear all
set more off
capture log close
log using "stata-lam/tests/reference_china.log", text replace
do "stata-lam/tests/reference_china.do"
display as result "REFERENCE COMPLETE"
log close
