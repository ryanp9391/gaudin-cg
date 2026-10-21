import pickle, time
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
import general_kernel as G
import master_kernel2 as m
import mk2_Nsym as S
def Y0_Y1_from_pts(L):
    Y1d=pickle.load(open('master_kernel2_Y1pts%d.pkl'%L,'rb'))[0]
    D=pickle.load(open('master_kernel2_Y1inv%d.pkl'%L,'rb'))
    return Y1d,D['kappa'],D['rho']
def N_from_pts(L):
    Y1d,kappa,rho=Y0_Y1_from_pts(L)
    n=len(kappa)
    r1={};r0={}
    for (r,c),v in Y1d.items():
        q=QQ(v.numerator,v.denominator)
        r1.setdefault(r,{})[c]=q
        if rho[r]==kappa[c]: r0.setdefault(r,{})[c]=q
    d1=DomainMatrix(r1,(n,n),QQ); d0=DomainMatrix(r0,(n,n),QQ)
    M=d1*d0.inv()
    return M.inv().to_Matrix(), M
if __name__=='__main__':
    t=time.time()
    N6,M6=N_from_pts(6); print('N6',time.time()-t,flush=True)
    pickle.dump({k:v for k,v in N6.todok().items()},open('master_kernel2_N6.pkl','wb'))
    # universality check vs L=5
    N5=S.Nmat(5)
    prs6=G.pairs_of((6,0,0)); n6=len(prs6)
    keep=[i for i,p in enumerate(prs6) if p[1]<=6]
    sel=[i*n6+j for i in keep for j in keep]
    print('N^(6)|_5 == N^(5):', N6[sel,sel]==N5)
