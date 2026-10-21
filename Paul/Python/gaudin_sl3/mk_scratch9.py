import sys
sys.path.insert(0,'.')
import sympy as sp
from math import comb
from fractions import Fraction as Fr
import general_kernel as G
for l1 in range(1,6):
    n=l1+1
    lam=(l1,0,0)
    c=lambda k,l: 0 if k==l else (sp.Rational(G.ckl_closed(lam,k,l)) if k<l else -sp.Rational(G.ckl_closed(lam,l,k)))
    nu=sp.symbols('nu1:%d'%(n+1))
    eqs=[]
    for k in range(n+1):
        for l in range(k+1,n+1):
            lhs=(k+l-1)*c(k,l)
            rhs=0
            if k>=1: rhs+=nu[k-1]*c(k-1,l)
            if l>=1 and l-1>=0: rhs+=nu[l-1]*c(k,l-1)
            eqs.append(lhs-rhs)
    sol=sp.solve(eqs,nu,dict=True)
    print(l1,sol)
