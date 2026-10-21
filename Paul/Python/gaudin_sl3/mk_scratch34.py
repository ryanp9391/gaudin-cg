import sys; sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as Fr
import general_kernel as G
from mk_lib import wedge2
from mk_scratch33 import Fmat, Tnum
a0,ap0,t0,s0=Fr(2,5),Fr(3,7),Fr(1,3),Fr(1,5)
for l1 in (1,2,3):
    n=l1+1; lam=(l1,0,0)
    T,rep=G.load_rep(lam); prs=G.pairs_of(lam)
    F2=wedge2(Fmat(n),prs)
    for sgn in (1,-1):
        s=sgn*s0
        U=(sp.Rational(s.numerator,s.denominator)*F2).exp()
        lhs0=Tnum(T,prs,a0,ap0,t0)
        for tname,tfun in (('t',lambda: t0),):
            a1=a0/(1-a0*s); ap1=ap0/(1-ap0*s)
            fac=sp.Rational(((1-a0*s)*(1-ap0*s))**(2*l1)).limit_denominator(10**12) if False else sp.Rational(((1-a0*s)*(1-ap0*s)).numerator,((1-a0*s)*(1-ap0*s)).denominator)**(2*l1)
            l=Tnum(T,prs,a1,ap1,t0)*fac
            r=U*lhs0*U.T
            d=l-r
            print(l1,'sgn',sgn,'ZERO' if d==sp.zeros(*d.shape) else 'nonzero',flush=True)
