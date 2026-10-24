*! version 0.1.0 23oct2026
version 16.0
mata:

// Log probability between finite standardized cutpoints, using the nearer tail.
real scalar lam_loginterval(real scalar lo, real scalar hi)
{
    real scalar a, d
    if (hi <= lo) return(.)
    if (lo > 0) {
        a = lnnormal(-lo)
        d = lnnormal(-hi) - a
    }
    else {
        a = lnnormal(hi)
        d = lnnormal(lo) - a
    }
    if (d >= 0) return(.)
    if (d > -ln(2)) return(a + ln(-expm1(d)))
    return(a + ln1p(-exp(d)))
}

// Transform unconstrained log gaps to ordered natural cutpoints.
void lam_transform(real rowvector p, real matrix ranges,
    real rowvector b, real matrix J)
{
    real scalar q, a, z, j, k
    b = p
    J = I(cols(p))
    for (q=1; q<=rows(ranges); q++) {
        a = ranges[q,1]
        z = ranges[q,2]
        for (j=a+1; j<=z; j++) {
            b[j] = b[j-1] + exp(p[j])
            for (k=a; k<=j; k++) {
                if (k==a) J[j,k] = 1
                else J[j,k] = exp(p[k])
            }
        }
    }
}

// C columns are mean index (zero for baseline), lower and upper cut indexes,
// and frequency. Analytical derivatives are first computed in natural units.
void lam_eval(real scalar todo, real rowvector p, real matrix C,
    real matrix ranges, real scalar value, real rowvector g, real matrix H)
{
    real rowvector b, gn, ids, s
    real matrix J, Hn, A
    real scalar i, m, l, u, mu, lo, hi, lp, rl, ru, n, j, k, q, a, z
    lam_transform(p, ranges, b, J)
    value = 0
    if (hasmissing(b)) {
        value=.
        return
    }
    gn = J(1,cols(p),0)
    Hn = J(cols(p),cols(p),0)
    for (i=1; i<=rows(C); i++) {
        m=C[i,1]; l=C[i,2]; u=C[i,3]; n=C[i,4]
        mu=0
        if (m) mu=b[m]
        lo=0; hi=0; rl=0; ru=0
        if (l) lo=b[l]-mu
        if (u) hi=b[u]-mu
        if (!l) lp=lnnormal(hi)
        else if (!u) lp=lnnormal(-lo)
        else lp=lam_loginterval(lo,hi)
        if (missing(lp)) {
            value=.
            return
        }
        value=value+n*lp
        if (!todo) continue
        if (l) rl=exp(lnnormalden(lo)-lp)
        if (u) ru=exp(lnnormalden(hi)-lp)
        ids=(m,l,u)
        s=(rl-ru,-rl,ru)
        for (j=1; j<=3; j++) {
            if (ids[j]) gn[ids[j]]=gn[ids[j]]+n*s[j]
        }
        if (todo<2) continue
        A=(lo*rl-hi*ru, -lo*rl, hi*ru \
           -lo*rl, lo*rl, 0 \
           hi*ru, 0, -hi*ru) - s'*s
        for (j=1; j<=3; j++) {
            if (!ids[j]) continue
            for (k=1; k<=3; k++) {
                if (ids[k]) Hn[ids[j],ids[k]]=Hn[ids[j],ids[k]]+n*A[j,k]
            }
        }
    }
    if (!todo) return
    g=gn*J
    if (todo<2) return
    H=J'*Hn*J
    for (q=1; q<=rows(ranges); q++) {
        a=ranges[q,1]; z=ranges[q,2]
        for (j=a+1; j<=z; j++) H[j,j]=H[j,j]+exp(p[j])*sum(gn[|j\z|])
    }
}

void lam_fit(string scalar yvar, string scalar qvar, string scalar tvar,
    string scalar wvar, string scalar sample, real scalar baseline,
    real scalar level, real scalar maxiter, real scalar nolog,
    string scalar bname, string scalar vname, string scalar trname,
    string scalar cname, string scalar catname, string scalar pname,
    string scalar sname)
{
    real matrix D, cells, info, cats, ranges, C, Vp, V, Jac, trend, cuts, adj
    real matrix H, chol
    real colvector times, ylevels, ix, seen, next, counts, mind, K, se
    real rowvector p, b, g
    string matrix stripe
    transmorphic S
    real scalar Q, T, base, i, q, j, a, z, k, r, N, nc, csum, pnum, crit, val

    D=st_data(.,(yvar,qvar,tvar,wvar),sample)
    N=sum(D[,4]); Q=max(D[,2]); times=uniqrows(sort(D[,3],1)); T=rows(times)
    if (missing(baseline)) baseline=times[1]
    ix=selectindex(times:==baseline)
    if (!rows(ix)) {
        errprintf("baseperiod() must be observed in the estimation sample\n")
        _error(198)
    }
    base=ix[1]
    counts=D[,3]
    for (i=1; i<=T; i++) {
        ix=selectindex(counts:==times[i]); D[ix,3]=J(rows(ix),1,i)
    }
    // Retain exact original ordinal codes, including nonconsecutive codes.
    cats=J(0,3,.)
    K=J(Q,1,.)
    for (q=1; q<=Q; q++) {
        ix=selectindex(D[,2]:==q)
        ylevels=uniqrows(sort(D[ix,1],1)); K[q]=rows(ylevels)
        if (K[q]<2) {
            errprintf("each survey must have at least two observed response categories\n")
            _error(459)
        }
        cats=cats \ (J(K[q],1,q),ylevels,(1::K[q]))
        // Map from an unmodified copy so original codes cannot collide with ranks.
        counts=D[ix,1]
        for (j=1; j<=K[q]; j++) D[ix[selectindex(counts:==ylevels[j])],1]=J(sum(counts:==ylevels[j]),1,j)
    }
    D=sort(D,(2,3,1))
    // panelsetup() uses a scalar grouping variable; mark the tuple boundaries.
    counts=J(rows(D),1,1)
    for (i=2; i<=rows(D); i++) counts[i]=any(D[i,1..3]:!=D[i-1,1..3])
    counts=runningsum(counts)
    info=panelsetup(counts,1)
    cells=J(rows(info),4,.)
    for (i=1; i<=rows(info); i++) {
        cells[i,1..3]=D[info[i,1],1..3]
        cells[i,4]=sum(D[|info[i,1],4 \ info[i,2],4|])
    }
    nc=rows(cells)
    // Connectivity of the bipartite survey-period graph is required.
    adj=J(Q,T,0)
    for (i=1; i<=nc; i++) adj[cells[i,2],cells[i,3]]=1
    seen=J(Q,1,0); seen[1]=1
    do {
        next=((adj*((adj'*seen):>0)):>0)
        if (all(next:==seen)) break
        seen=next
    } while (1)
    if (!all(seen)) {
        errprintf("survey overlap graph is disconnected\n")
        errprintf("fit each connected component separately or supply bridging observations\n")
        _error(459)
    }
    mind=J(T,1,0); k=0
    for (i=1; i<=T; i++) if (i!=base) {
        k++
        mind[i]=k
    }
    pnum=T-1+sum(K:-1)
    ranges=J(Q,2,.)
    p=J(1,pnum,0); stripe=J(pnum,2,"")
    for (i=1; i<=T; i++) if (mind[i]) stripe[mind[i],.]=("mu","t"+strofreal(i))
    cuts=J(sum(K:-1),4,.)
    r=0
    for (q=1; q<=Q; q++) {
        a=k+1; z=k+K[q]-1; ranges[q,.]=(a,z)
        ix=selectindex(cells[,2]:==q)
        counts=J(K[q],1,0)
        for (i=1; i<=rows(ix); i++) counts[cells[ix[i],1]]=counts[cells[ix[i],1]]+cells[ix[i],4]
        csum=0
        ylevels=select(cats[,2],cats[,1]:==q)
        for (j=1; j<K[q]; j++) {
            csum=csum+counts[j]
            p[a+j-1]=invnormal(csum/sum(counts))
            stripe[a+j-1,.]=("cut"+strofreal(q),"c"+strofreal(j))
            r++; cuts[r,1..2]=(q,ylevels[j])
        }
        // Backwards conversion preserves the natural starting thresholds.
        for (j=z; j>a; j--) p[j]=ln(p[j]-p[j-1])
        k=z
    }
    C=J(nc,4,0)
    for (i=1; i<=nc; i++) {
        q=cells[i,2]; j=cells[i,1]
        C[i,1]=mind[cells[i,3]]
        if (j>1) C[i,2]=ranges[q,1]+j-2
        if (j<K[q]) C[i,3]=ranges[q,1]+j-1
        C[i,4]=cells[i,4]
    }
    // A finite MLE requires no improving recession direction. Each occupied
    // response interval imposes d_lower <= d_mean <= d_upper. The directed
    // inequality graph must be strongly connected after fixing the baseline.
    // All observed interior categories already imply adjacent cut ordering.
    adj=J(pnum+1,pnum+1,0)
    for (i=1; i<=nc; i++) {
        j=cells[i,3]
        if (C[i,2]) adj[C[i,2]+1,j]=1
        if (C[i,3]) adj[j,C[i,3]+1]=1
    }
    for (r=1; r<=2; r++) {
        seen=J(pnum+1,1,0); seen[base]=1
        do {
            next=((seen+adj*seen):>0)
            if (all(next:==seen)) break
            seen=next
        } while (1)
        if (!all(seen)) {
            errprintf("response patterns are separated; a finite maximum likelihood estimate does not exist\n")
            errprintf("inspect extreme-only periods and survey-period response support\n")
            _error(430)
        }
        adj=adj'
    }
    S=optimize_init()
    optimize_init_evaluator(S,&lam_eval())
    optimize_init_evaluatortype(S,"d2")
    optimize_init_argument(S,1,C)
    optimize_init_argument(S,2,ranges)
    optimize_init_params(S,p)
    optimize_init_technique(S,"nr")
    optimize_init_conv_maxiter(S,maxiter)
    optimize_init_conv_ptol(S,1e-9)
    optimize_init_conv_vtol(S,1e-11)
    optimize_init_conv_nrtol(S,1e-7)
    optimize_init_valueid(S,"Log likelihood")
    if (nolog) optimize_init_tracelevel(S,"none")
    p=optimize(S)
    if (!optimize_result_converged(S)) {
        errprintf("LAM did not converge; inspect sparse cells or increase iterate()\n")
        _error(430)
    }
    lam_eval(2,p,C,ranges,val,g,H)
    chol=cholesky(-H)
    if (hasmissing(chol)) {
        errprintf("observed information is not positive definite; finite estimates are not identified\n")
        _error(430)
    }
    Vp=invsym(-H)
    if (hasmissing(Vp) | any(diagonal(Vp):<=0)) {
        errprintf("observed information is singular; check sparse or separated response patterns\n")
        _error(430)
    }
    lam_transform(p,ranges,b,Jac)
    V=Jac*Vp*Jac'; V=(V+V')/2
    se=sqrt(diagonal(V))
    crit=invnormal(1-(1-level/100)/2)
    trend=J(T,5,0); trend[,1]=times
    for (i=1; i<=T; i++) if (mind[i]) {
        j=mind[i]; trend[i,2..5]=(b[j],se[j],b[j]-crit*se[j],b[j]+crit*se[j])
    }
    cuts[,3]=b[|T\pnum|]'; cuts[,4]=se[|T\pnum|]
    st_matrix(bname,b); st_matrixcolstripe(bname,stripe)
    st_matrix(vname,V); st_matrixcolstripe(vname,stripe); st_matrixrowstripe(vname,stripe)
    st_matrix(trname,trend)
    st_matrixcolstripe(trname,(J(5,1,""),("time"\"estimate"\"se"\"ll"\"ul")))
    st_matrix(cname,cuts)
    st_matrixcolstripe(cname,(J(4,1,""),("survey_id"\"category"\"threshold"\"se")))
    st_matrix(catname,cats)
    st_matrixcolstripe(catname,(J(3,1,""),("survey_id"\"category"\"rank")))
    st_matrix(pname,times')
    st_matrix(sname,(N,nc,val,optimize_result_converged(S),optimize_result_iterations(S),baseline))
}
// Keep exact identifiers out of Stata macro expansion and display formatting.
void lam_surveymap(string scalar survey, string scalar idvar, string scalar sample)
{
    real colvector ids, vals, qdata, ix
    string colvector strings
    real scalar q, Q
    ids=st_data(.,idvar,sample); Q=max(ids)
    vals=J(Q,1,.)
    if (st_global("e(surveytype)")=="string") {
        strings=st_sdata(.,survey,sample)
        for (q=1; q<=Q; q++) {
            ix=selectindex(ids:==q)
            st_global("e(survey"+strofreal(q)+")",strings[ix[1]])
        }
    }
    else {
        qdata=st_data(.,survey,sample)
        for (q=1; q<=Q; q++) {
            ix=selectindex(ids:==q); vals[q]=qdata[ix[1]]
            st_global("e(survey"+strofreal(q)+")",strtrim(strofreal(vals[q],"%24.17g")))
        }
    }
    st_matrix("e(surveyvalues)",vals)
}

// Reconstruct from e(b) and e(V) so bootstrap/estimates restore use current VCE.
real matrix lam_current_trend(real scalar level)
{
    real matrix tr, V
    real rowvector b
    real colvector times, ix
    string matrix stripe
    real scalar i, j, z
    times=st_matrix("e(periods)")'
    b=st_matrix("e(b)"); V=st_matrix("e(V)")
    stripe=st_matrixcolstripe("e(b)")
    tr=J(rows(times),5,0); tr[,1]=times
    z=invnormal(1-(1-level/100)/2)
    for (i=1; i<=rows(times); i++) {
        if (times[i]==st_numscalar("e(baseperiod)")) continue
        ix=selectindex((stripe[,1]:=="mu") :& (stripe[,2]:=="t"+strofreal(i)))
        if (rows(ix)!=1) _error(498)
        j=ix[1]
        tr[i,2]=b[j]; tr[i,3]=sqrt(V[j,j])
        tr[i,4]=tr[i,2]-z*tr[i,3]; tr[i,5]=tr[i,2]+z*tr[i,3]
    }
    return(tr)
}

void lam_predict(string scalar newvar, string scalar sample,
    string scalar stat, real scalar outcome)
{
    real matrix tr, cats, sv
    real rowvector b
    string matrix stripe
    real colvector t, fit, qid, ix, qvalues, rank, cutidx
    string colvector qstrings
    real scalar i, q, r, j, lo, hi, K, lp
    tr=lam_current_trend(95)
    cats=st_matrix("e(categories)"); b=st_matrix("e(b)")
    stripe=st_matrixcolstripe("e(b)")
    t=st_data(.,st_global("e(timevar)"),sample)
    fit=J(rows(t),1,.)
    for (i=1; i<=rows(tr); i++) {
        ix=selectindex(t:==tr[i,1])
        if (stat=="stdp") fit[ix]=J(rows(ix),1,tr[i,3])
        else fit[ix]=J(rows(ix),1,tr[i,2])
    }
    if (stat!="pr") {
        st_store(.,newvar,sample,fit)
        return
    }
    qid=J(rows(t),1,.)
    if (st_global("e(surveytype)")=="string") {
        qstrings=st_sdata(.,st_global("e(surveyvar)"),sample)
        for (q=1; q<=st_numscalar("e(n_surveys)"); q++) {
            ix=selectindex(qstrings:==st_global("e(survey"+strofreal(q)+")"))
            qid[ix]=J(rows(ix),1,q)
        }
    }
    else {
        sv=st_matrix("e(surveyvalues)")
        qvalues=st_data(.,st_global("e(surveyvar)"),sample)
        for (q=1; q<=st_numscalar("e(n_surveys)"); q++) {
            ix=selectindex(qvalues:==sv[q,1])
            qid[ix]=J(rows(ix),1,q)
        }
    }
    for (i=1; i<=rows(t); i++) {
        q=qid[i]
        if (missing(q) | missing(fit[i])) {
            fit[i]=.
            continue
        }
        ix=selectindex((cats[,1]:==q) :& (cats[,2]:==outcome))
        if (!rows(ix)) {
            fit[i]=.
            continue
        }
        r=cats[ix[1],3]
        K=sum(cats[,1]:==q)
        cutidx=selectindex(stripe[,1]:=="cut"+strofreal(q))
        if (r==1) lp=lnnormal(b[cutidx[1]]-fit[i])
        else if (r==K) lp=lnnormal(fit[i]-b[cutidx[K-1]])
        else lp=lam_loginterval(b[cutidx[r-1]]-fit[i],b[cutidx[r]]-fit[i])
        fit[i]=exp(lp)
    }
    st_store(.,newvar,sample,fit)
}

// Define the load sentinel last so an incomplete compile is never accepted.
real scalar lam_version() return(1)
end
