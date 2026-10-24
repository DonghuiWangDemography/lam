*! version 0.1.0 23oct2026
program define lam_p
    version 16.0
    if "`e(cmd)'" != "lam" & "`e(cmdname)'" != "lam" error 301
    syntax newvarname [if] [in] [, XB STDP PR OUTcome(numlist max=1)]
    local nopts = ("`xb'"!="") + ("`stdp'"!="") + ("`pr'"!="")
    if `nopts' > 1 {
        di as error "specify only one of xb, stdp, or pr"
        exit 198
    }
    if "`pr'" != "" & "`outcome'" == "" {
        di as error "pr requires outcome() with an original numeric response code"
        exit 198
    }
    if "`outcome'" != "" & "`pr'" == "" {
        di as error "outcome() requires pr"
        exit 198
    }
    if "`outcome'" == "" local outcome .
    local stat "xb"
    if "`stdp'" != "" local stat "stdp"
    if "`pr'" != "" local stat "pr"
    confirm numeric variable `e(timevar)'
    if "`stat'" == "pr" {
        confirm `e(surveytype)' variable `e(surveyvar)'
    }
    marksample touse, novarlist
    capture mata: lam_version()
    if _rc {
        quietly findfile lam_mata.mata
        quietly do `"`r(fn)'"'
    }
    quietly generate `typlist' `varlist' = .
    mata: lam_predict("`varlist'", "`touse'", "`stat'", `outcome')
    if "`stat'" == "xb" label variable `varlist' "LAM population mean for the period"
    if "`stat'" == "stdp" label variable `varlist' "Standard error of LAM period mean"
    if "`stat'" == "pr" label variable `varlist' "LAM probability of response `outcome'"
end
