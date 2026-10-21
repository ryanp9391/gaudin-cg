import sys
sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as Fr
import general_kernel as G
from master_kernel import *
def Nmat(n):
    M=sp.zeros(n+1,n+1)
    for k in range(n): M[k+1,k]=-(n-k)
    return M
def wedge2(M,prs):
    idx={p:i for i,p in enumerate(prs)}
    out=sp.zeros(len(prs),len(prs))
    for (k,l),j in idx.items():
        for r in range(M.shape[0]):
            if M[r,k]!=0 and r!=l:
                s=1; a,b=r,l
                if a>b: a,b=b,a; s=-1
                out[idx[(a,b)],j]+=s*M[r,k]
            if M[r,l]!=0 and r!=k:
                s=1; a,b=k,r
                if a>b: a,b=b,a; s=-1
                out[idx[(a,b)],j]+=s*M[r,l]
    return out
def Tnum(T,prs,a,ap,t):
    return sp.Matrix(len(prs),len(prs),lambda i,j: sp.Rational(G.eval_T(T[(i,j)],Fr(a),Fr(ap),Fr(t)).numerator, G.eval_T(T[(i,j)],Fr(a),Fr(ap),Fr(t)).denominator))
if __name__!="__main__": raise SystemExit
a0,ap0,t0,s0=Fr(2,5),Fr(3,7),Fr(1,3),Fr(1,5)
for l1 in (1,2,3):
    n=l1+1; lam=(l1,0,0)
    T,rep=G.load_rep(lam); prs=G.pairs_of(lam)
    N2=wedge2(Nmat(n),prs)
    U=(sp.Rational(s0.numerator,s0.denominator)*N2).exp()
    lhs0=Tnum(T,prs,a0,ap0,t0)
    for name,tp in (('t',t0),('t aa/((a+s)(a+s))',t0*a0*ap0/((a0+s0)*(ap0+s0)))):
        l=Tnum(T,prs,a0+s0,ap0+s0,tp)
        r=U*lhs0*U.T
        print(l1,name,'zero' if (l-r)==sp.zeros(*l.shape) else 'NONZERO',flush=True)
