"""Q^{(lam)} from the universal normal symbol exp(S): apply to F_a^lam G_a'^lam, numeric a,a'. Compare with exact data."""
import sys, pickle, itertools, time
import sympy as sp
from sympy import Rational as R, ff
from mk2_expfit import *
import master_kernel2 as m
def brk(f,g): return sp.expand(sum(sp.diff(f,p)*sp.diff(g,q) for p,q in zip(PV,V)))
def Dlist(nmax, D3=None, D4=None):
    D={1:D1}
    D[2]=brk(D1,D1)
    if nmax>=3:
        D[3]=sp.expand(sp.Rational(1,2)*(brk(D[2],D[1])+brk(D[1],D[2])))
    import os
    for n in range(4,nmax+1):
        f='master_kernel2_D%d_gen.pkl'%n
        if os.path.exists(f):
            D[n]=sym_from_terms(pickle.load(open(f,'rb')))
        else:
            break
    return D
def exp_series(D,nmax):
    """coefficient of t^n in exp(sum t^n D_n), as polynomials"""
    t=sp.Symbol('t')
    S=sum(D[n]*t**n for n in D)
    ser=sp.Integer(1); term=sp.Integer(1)
    out={0:sp.Integer(1)}
    # exp series truncated
    tot=0; pw=sp.Integer(1)
    for k in range(0,nmax+1):
        tot+=pw/sp.factorial(k)
        pw=sp.expand(pw*S)
        pw=sum(pw.coeff(t,j)*t**j for j in range(nmax+1))
    tot=sp.expand(tot)
    return {n:sp.expand(tot.coeff(t,n)) for n in range(nmax+1)}
def apply_to_power(symb, lam, a, ap):
    """symb (poly in V,PV) applied to F_a^lam G_a'^lam; returns (P) with result = P * F^(lam-?)...: returns Q_n * Pi^{-(lam-n)} style: we return Q = result expression (full polynomial)"""
    F=1-a*E1+a**2*E2; G=1-ap*H1+ap**2*H2
    P=sp.Poly(sp.expand(symb),*PV)
    out=0
    for (k1,k2,l1,l2),co in P.terms():
        i=k1+k2; j=l1+l2
        fac=ff(lam,i)*ff(lam,j)*(-a)**k1*(a**2)**k2*(-ap)**l1*(ap**2)**l2
        if lam-i<0 or lam-j<0: continue
        out+=co*fac*F**(lam-i)*G**(lam-j)
    return sp.expand(out)
if __name__=='__main__':
    nmax=int(sys.argv[1]); Lmax=int(sys.argv[2])
    a=R(1,3); ap=R(1,5)
    D=Dlist(nmax)
    Phi=exp_series(D,nmax)
    for L in range(int(sys.argv[3]) if len(sys.argv)>3 else 1,Lmax+1):
        Q=m.symbol_orders_exact(L,nmax=nmax)
        for n in range(1,nmax+1):
            pred=apply_to_power(Phi[n],L,a,ap)
            exact=sp.expand(Q.get(n,0).subs({m.a_:a,m.ap_:ap})) if n in Q else 0
            print('lam',L,'n',n,'operator Q_n == exact Q_n :',sp.expand(pred-exact)==0,flush=True)
