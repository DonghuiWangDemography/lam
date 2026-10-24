version 16.0
clear all
set more off
args root
if `"`root'"' == "" local root `"`c(pwd)'"'
capture log close
log using `"`root'/stata-lam/tests/test_derivatives.log"', text replace
quietly do `"`root'/stata-lam/lam_mata.mata"'

* Check gradients and Hessians away from the optimum with central differences
mata
void lam_test_derivatives()
{
    real rowvector p, pp, pm, g, gp, gm, gn
    real matrix ranges, C, H, Hp, Hm, Hn
    real scalar v, vp, vm, h, j
    p=(-.3,.45,-1.1,ln(.9),ln(1.3),-.7,ln(1.8))
    ranges=(3,5\6,7)
    C=(0,0,3,12\0,3,4,20\0,4,5,30\0,5,0,10\
       1,0,3,15\1,3,4,25\1,4,5,25\1,5,0,8\
       1,0,6,15\1,6,7,35\1,7,0,18\
       2,0,6,12\2,6,7,28\2,7,0,20)
    lam_eval(2,p,C,ranges,v,g,H)
    gn=J(1,cols(p),.)
    Hn=J(cols(p),cols(p),.)
    h=1e-5
    for (j=1; j<=cols(p); j++) {
        pp=p; pm=p
        pp[j]=pp[j]+h
        pm[j]=pm[j]-h
        lam_eval(1,pp,C,ranges,vp,gp,Hp)
        lam_eval(1,pm,C,ranges,vm,gm,Hm)
        gn[j]=(vp-vm)/(2*h)
        Hn[,j]=(gp-gm)'/(2*h)
    }
    st_numscalar("gradient_difference",max(abs(g-gn)))
    st_numscalar("hessian_difference",max(abs(H-Hn)))
    st_numscalar("hessian_asymmetry",max(abs(H-H')))
    st_numscalar("tail_logprob",lam_loginterval(30,31))
    st_numscalar("reflected_tail_logprob",lam_loginterval(-31,-30))
}
lam_test_derivatives()
end
assert gradient_difference<1e-6
assert hessian_difference<1e-6
assert hessian_asymmetry<1e-10
assert abs(tail_logprob-(-454.3212439563433))<1e-10
assert abs(tail_logprob-reflected_tail_logprob)<1e-10
display as result "Maximum finite-difference gradient discrepancy = " %14.10g gradient_difference
display as result "Maximum finite-difference Hessian discrepancy = " %14.10g hessian_difference
display as result "DERIVATIVE AND TAIL TESTS PASSED"
log close
