*! version 0.1.0 23oct2026
program define lam_plot, rclass
    version 16.0
    if "`e(cmd)'" != "lam" & "`e(cmdname)'" != "lam" error 301
    syntax [, Level(numlist max=1 >0 <100) *]
    if "`level'" == "" local level = e(level)
    capture mata: lam_version()
    if _rc {
        quietly findfile lam_mata.mata
        quietly do `"`r(fn)'"'
    }
    tempname trend
    mata: st_matrix("`trend'",lam_current_trend(`level'))
    matrix colnames `trend' = time estimate se ll ul
    preserve
    quietly {
        clear
        svmat double `trend'
        rename (`trend'1 `trend'2 `trend'3 `trend'4 `trend'5) ///
            (period estimate se lower upper)
    }
    capture noisily twoway ///
        (rarea lower upper period, color(navy%15) lwidth(none)) ///
        (connected estimate period, lcolor(navy) mcolor(navy) msymbol(O) msize(small)), ///
        xtitle("`e(timevar)'") ytitle("Latent attitude relative to `e(baseperiod)'") ///
        legend(order(2 "Estimate" 1 "`level'% confidence interval")) ///
        graphregion(color(white)) `options'
    local rc = _rc
    restore
    if !`rc' return matrix trend = `trend'
    exit `rc'
end
