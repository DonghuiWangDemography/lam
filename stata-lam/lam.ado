*! version 0.1.0 23oct2026
program define lam, eclass sortpreserve
    version 16.0
    if replay() {
        if "`e(cmd)'" != "lam" error 301
        syntax [, Level(cilevel)]
        if "`level'" == "" local level = e(level)
        _lam_display, level(`level')
        exit
    }
    local cmdline `"lam `0'"'
    syntax varlist(min=1 max=1 numeric) [if] [in] [fweight], ///
        SURvey(varname) TIME(varname) ///
        [BASEPeriod(numlist max=1) Level(cilevel) ITERate(integer 300) noLOG]
    confirm numeric variable `time'
    if "`baseperiod'" == "" local baseperiod .
    if `iterate' < 1 {
        di as error "iterate() must be a positive integer"
        exit 198
    }
    marksample touse
    markout `touse' `time'
    quietly replace `touse' = 0 if missing(`survey')
    tempvar w qid
    quietly {
        if "`weight'" == "" generate double `w' = 1
        else generate double `w' `exp'
        count if `touse' & !missing(`w') & (`w' < 0 | `w' != floor(`w'))
        if r(N) {
            noisily di as error "frequency weights must be nonnegative integers"
            exit 401
        }
        replace `touse' = 0 if missing(`w') | `w' == 0
        count if `touse'
        if r(N) == 0 error 2000
        local N_rows = r(N)
        egen long `qid' = group(`survey') if `touse'
        summarize `qid', meanonly
        local Q = r(max)
    }
    capture mata: lam_version()
    if _rc {
        quietly findfile lam_mata.mata
        quietly do `"`r(fn)'"'
    }
    tempname b V trend cuts cats periods stats
    mata: lam_fit("`varlist'", "`qid'", "`time'", "`w'", "`touse'", ///
        `baseperiod', `level', `iterate', "`log'" == "nolog", ///
        "`b'", "`V'", "`trend'", "`cuts'", "`cats'", "`periods'", "`stats'")
    ereturn post `b' `V', esample(`touse') depname(`varlist') obs(`=el(`stats',1,1)')
    ereturn scalar N_rows = `N_rows'
    ereturn scalar N_cells = `stats'[1,2]
    ereturn scalar n_surveys = `Q'
    ereturn scalar n_periods = colsof(`periods')
    ereturn scalar k = colsof(e(b))
    ereturn scalar rank = colsof(e(b))
    ereturn scalar ll = `stats'[1,3]
    ereturn scalar converged = `stats'[1,4]
    ereturn scalar iterations = `stats'[1,5]
    ereturn scalar baseperiod = `stats'[1,6]
    ereturn scalar level = `level'
    ereturn scalar aic = -2*e(ll) + 2*e(k)
    ereturn scalar bic = -2*e(ll) + ln(e(N))*e(k)
    ereturn matrix trend = `trend'
    ereturn matrix cutpoints = `cuts'
    ereturn matrix categories = `cats'
    ereturn matrix periods = `periods'
    capture confirm string variable `survey'
    if !_rc ereturn local surveytype "string"
    else ereturn local surveytype "numeric"
    tempvar map_sample
    quietly generate byte `map_sample' = e(sample)
    mata: lam_surveymap("`survey'", "`qid'", "`map_sample'")
    ereturn local surveyvar "`survey'"
    ereturn local timevar "`time'"
    ereturn local wtype "`weight'"
    ereturn local wexp "`exp'"
    ereturn local vce "oim"
    ereturn local vcetype "OIM"
    ereturn local predict "lam_p"
    ereturn local properties "b V"
    ereturn local title "Latent attitude method"
    ereturn local cmdline `"`cmdline'"'
    ereturn local cmd "lam"
    _lam_display, level(`level')
end

program define _lam_display
    version 16.0
    syntax [, Level(cilevel)]
    di as text _newline "Latent attitude method"
    di as text "Observations    " as result %12.0fc e(N) ///
        as text "     Survey series " as result %6.0f e(n_surveys)
    di as text "Log likelihood  " as result %12.5f e(ll) ///
        as text "     Time periods  " as result %6.0f e(n_periods)
    di as text "Baseline period " as result %12.0g e(baseperiod)
    di as text "Constant survey-specific cutpoints and unit latent variance"
    ereturn display, level(`level')
    di as text "Period estimates including the baseline are in e(trend)"
    di as text "Survey and response mappings are in e(survey#) and e(categories)"
end
