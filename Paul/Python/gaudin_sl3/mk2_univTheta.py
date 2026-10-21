import pickle
from fractions import Fraction as Fr
from math import comb
import sympy as sp
import general_kernel as G
import master_kernel2 as m
import mk2_Nsym as S
Nbig=S.Nmat(5); Lbig=5
prsb=G.pairs_of((Lbig,0,0)); nb=len(prsb)
index={p:i for i,p in enumerate(prsb)}
def ckl(lam,k,l): return Fr((-1)**(k+l-1)*comb(lam,l-1)*comb(lam+1,k)*(l-k),l)
def theta_from_N(lam):
    prs=G.pairs_of((lam,0,0)); n=len(prs)
    lev=lambda p:p[0]+p[1]-1
    T={}
    for ip,p in enumerate(prs):
        for iq,q in enumerate(prs):
            col=index[p]*nb+index[q]
            d={}
            for r in prs:
                for s in prs:
                    row=index[r]*nb+index[s]
                    v=Nbig[row,col]
                    if v==0: continue
                    nn=lev(r)+lev(s)-lev(p)-lev(q)
                    key=(nn,lev(r),lev(s))
                    val=Fr(int(v.p),int(v.q))*ckl(lam,*r)*ckl(lam,*s)
                    d[key]=d.get(key,0)+val
            T[(ip,iq)]={k:v for k,v in d.items() if v!=0}
    return T
for lam in (2,3,4,5):
    T=theta_from_N(lam)
    T0=m.get_T(lam)
    bad=0
    for k in T0:
        a={kk:Fr(v) for kk,v in T0[k].items() if v!=0}
        if a!=T[k]: bad+=1
    print('lam',lam,'Theta from universal N (L=5) and closed-form c_kl(lam): mismatching entries',bad,'of',len(T0))
