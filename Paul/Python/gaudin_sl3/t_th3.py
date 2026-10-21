from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
T,rep=gk.load_rep(lam)
prs=gk.pairs_of(lam); n=len(prs)
evs=ev_cov(L)
for vals in ({a_s:sp.Rational(3,7),ap_s:sp.Rational(-5,4),t_s:sp.Rational(2,9)},{a_s:sp.Rational(-2,3),ap_s:sp.Rational(7,5),t_s:sp.Rational(-3,11)}):
    tv=vals[t_s]
    Yt=sum((coef[k]*tv**k for k in range(len(coef))),sp.zeros(*coef[0].shape))
    ev=evs.subs(vals)
    th=Yt.T.LUsolve(ev.T)   # th Y = ev  => Y^T th^T = ev^T
    bad=0
    for i in range(n):
        for j in range(n):
            ref=sum(sp.Rational(c.numerator,c.denominator)*tv**nn*vals[a_s]**p*vals[ap_s]**q for (nn,p,q),c in T[(i,j)].items())
            if th[i*n+j]!=ref: bad+=1
    print('L',L,'mismatches',bad,'of',n*n)
