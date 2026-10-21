"""cumulant fit with numeric a,a': U_n(lam) = Pi^n [t^n] log( sum Q_m t^m / Pi^lam ), expect polynomial in lam of degree n+1 with no constant/linear term"""
import sys, pickle, itertools, time
import sympy as sp
from sympy import Rational as R
import master_kernel2 as m
a=R(sys.argv[1]); ap=R(sys.argv[2]); Lmax=int(sys.argv[3]); NMAX=int(sys.argv[4])
gens=(m.E1,m.E2,m.H1,m.H2)
Pi=sp.Poly(sp.expand(m.Pi_.subs({m.a_:a,m.ap_:ap})),*gens)
Pipow={0:sp.Poly(1,*gens)}
def pipow(k):
    if k not in Pipow: Pipow[k]=pipow(k-1)*Pi
    return Pipow[k]
def Pn_of(L,n):
    Q=m.symbol_orders_exact(L,nmax=NMAX)
    out={}
    for k in range(0,NMAX+1):
        q=sp.Poly(sp.expand(Q.get(k,0).subs({m.a_:a,m.ap_:ap})) if k in Q else 0,*gens)
        if k<=L:
            qq,r=sp.div(q,pipow(L-k)); assert r.is_zero,(L,k)
            out[k]=qq
        else:
            out[k]=q*pipow(k-L)
    return out
def cumulants(P,N):
    # P[m] polynomials (P[0]=1), t-series S = sum P_m s^m (s=t/Pi), log S
    # use recursion: n W_n = n P_n - sum_{k=1}^{n-1} k W_k P_{n-k}
    W={}
    for n in range(1,N+1):
        acc=n*P[n]
        for k in range(1,n):
            acc=acc-k*W[k]*P[n-k]
        W[n]=acc*sp.Rational(1,n)
    return W
t0=time.time()
data={}
for L in range(int(sys.argv[5]) if len(sys.argv)>5 else 1,Lmax+1):
    P=Pn_of(L,NMAX); P[0]=sp.Poly(1,*gens)
    data[L]=cumulants(P,NMAX)
    print('lam',L,'done',round(time.time()-t0,1),flush=True)
pickle.dump({L:{n:w.as_expr() for n,w in d.items()} for L,d in data.items()},open('master_kernel2_cum_%s_%s.pkl'%(sys.argv[1].replace('/','o'),sys.argv[2].replace('/','o')),'wb'))
lam=sp.Symbol('lam')
for n in range(1,NMAX+1):
    Ls=list(range(1,Lmax+1))
    # unknowns: coefficients of lam^2..lam^(n+1) ; also test lam^1 (should vanish) by including it
    ks=list(range(1,n+2))
    V=sp.Matrix([[L**k for k in ks] for L in Ls])
    # solve least squares exactly using first len(ks) points, check on rest
    use=Ls[:len(ks)]; chk=Ls[len(ks):]
    Vu=sp.Matrix([[L**k for k in ks] for L in use]).inv() if len(use)==len(ks) else None
    if Vu is None:
        print('n',n,'not enough data'); continue
    coeffs=[]
    for ki in range(len(ks)):
        coeffs.append(sum((Vu[ki,i]*data[use[i]][n] for i in range(len(use))),sp.Poly(0,*gens)))
    ok=all((sum((coeffs[ki]*L**ks[ki] for ki in range(len(ks))),sp.Poly(0,*gens))-data[L][n]).is_zero for L in chk)
    print('n',n,'fit degree',n+1,'(incl lam^1) checks on',chk,'all ok:',ok,' lam^1 coeff zero:',coeffs[0].is_zero,' term counts',[len(c.terms()) for c in coeffs])
