import sys, itertools, pickle
import sympy as sp
from sympy import Rational as R
from mk2_k7 import *
import mk2_expfit as X
lam=int(sys.argv[1]); a=R(1,3); ap=R(1,5)
ms,idx=space(lam)
Qn=build(lam,a,ap)
E1o,E2o,Om=casimir(lam)
n=Om.shape[0]
vs={k:vec(d,lam) for k,d in Qn.items()}
# N1 as matrix: apply symbol x to E-monomials
def Nsym_matrix(symb):
    M=sp.zeros(n,n)
    for (u,v),c in idx.items():
        f=X.E1**u[1]*X.E2**u[2]*X.H1**v[1]*X.H2**v[2]
        g=sp.Poly(X.apply_symbol(symb,f),X.E1,X.E2,X.H1,X.H2)
        for (i,j,k,l),co in g.terms():
            u2=(lam-i-j,i,j); v2=(lam-k-l,k,l)
            M[idx[(u2,v2)],c]+=co
    return M
A=Nsym_matrix(X.D1)
v0=vs[0]
z=a*ap/(ap-a)
print('A v0 == Q1 / ? :', sp.simplify((A*v0-vs[1]/1)==sp.zeros(n,1)) if False else (A*v0==vs[1]))
print('[A,Om]==0 ?',(A*Om-Om*A)==sp.zeros(n,n))
P=projectors(lam,Om)
WQ=sp.Matrix.hstack(*[Om**k*vs[kk] for kk in sorted(vs) for k in range(lam+2)])
rQ={mu:sp.Matrix.hstack(*[p*Om**k*vs[kk] for kk in sorted(vs) for k in range(1)]).rank() for mu,p in P.items()}
print('dim W_Q',WQ.rank(),{mu:sp.Matrix.hstack(*[p*vs[kk] for kk in sorted(vs)]).rank() for mu,p in P.items()})
# Krylov with words in (Om,A) applied to v0, increasing A-count
def span_basis(vecs):
    M=sp.Matrix.hstack(*vecs)
    return M.columnspace()
cur=[v0]
def close_under_Om(vecs):
    B=span_basis(vecs)
    allv=list(B)
    while True:
        new=allv+[Om*b for b in allv]
        Bn=span_basis(new)
        if len(Bn)==len(allv): return Bn
        allv=list(Bn)
W=close_under_Om(cur)
for j in range(1,9):
    inQ=sp.Matrix.hstack(*W,*[vs[k] for k in vs]).rank()==len(W)
    print('A-count',j-1,'dim',len(W),'W_Q contained:',inQ, ' per irrep dims',{mu:sp.Matrix.hstack(*[p*w for w in W]).rank() for mu,p in P.items()},flush=True)
    W=close_under_Om(list(W)+[A*w for w in W])
