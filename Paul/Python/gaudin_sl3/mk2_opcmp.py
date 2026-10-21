import sys, itertools, pickle, sympy as sp
import mk2_opfit as F
from mk2_opfit import *
c={(1,1):E1-H1,(1,2):E2-H2,(2,1):E2-H2,(2,2):E2*H1-E1*H2}
dE={1:(1,0),2:(0,1)}
def N1sq_normal():
    """normal-ordered square (coefficients left): dict (e,d)->coeff"""
    out={}
    for (i,j),ci in c.items():
        for (k,l),ck in c.items():
            poly=sp.Poly(sp.expand(ci*ck),*V)
            d=(dE[i][0]+dE[k][0],dE[i][1]+dE[k][1],dE[j][0]+dE[l][0],dE[j][1]+dE[l][1])
            for e,co in poly.terms():
                out[(e,d)]=out.get((e,d),0)+co
    return out
if __name__=='__main__':
    L=int(sys.argv[1]);n=2
    terms,cs,sol,shape=F.fit(L,n,2)
    s=list(sol)[0]
    N2={t:v for t,v in zip(terms,s) if v!=0}
    sq=N1sq_normal()
    diff={}
    for k in set(N2)|set(sq):
        v=N2.get(k,0)-sp.Rational(1,2)*sq.get(k,0)
        if v!=0: diff[k]=v
    print('N2 terms',len(N2),'normal-ordered N1^2/2 terms',len(sq),'difference terms',len(diff))
    # group difference by derivative order
    byd={}
    for (e,d),v in diff.items(): byd.setdefault(d,[]).append((e,v))
    for d,l in sorted(byd.items()): print(d,l)
    pickle.dump((N2,diff),open('master_kernel2_N2op.pkl','wb'))
