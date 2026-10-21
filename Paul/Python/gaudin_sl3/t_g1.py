from jet_basis_sym import *
import sys
L=int(sys.argv[1])
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
Y0i=coef[0].inv(); n=coef[0].shape[0]
T=[Y0i*c for c in coef]   # T_k = Y0^{-1} Y_k ; T_0=1
K=len(T)-1
# series of T^{-1}
def mul(P,Q,deg):
    R=[sp.zeros(n,n) for _ in range(deg+1)]
    for i,p in enumerate(P):
        for j,q in enumerate(Q):
            if i+j<=deg: R[i+j]+=p*q
    return R
deg=K+2
Tp=T+[sp.zeros(n,n)]*(deg+1-len(T))
Ti=[sp.eye(n)]+[sp.zeros(n,n) for _ in range(deg)]
# Ti = sum (-X)^k
X=[sp.zeros(n,n)]+Tp[1:]
pw=[sp.eye(n)]+[sp.zeros(n,n) for _ in range(deg)]
Ti=[sp.eye(n)]+[sp.zeros(n,n) for _ in range(deg)]
for k in range(1,deg+1):
    pw=mul(pw,X,deg)
    for i in range(deg+1): Ti[i]+= (-1)**k*pw[i]
dT=[ (k+1)*Tp[k+1] for k in range(deg)]
Gr=mul(Ti,dT+[sp.zeros(n,n)],deg-1)   # T^{-1} T'
Gl=mul(dT+[sp.zeros(n,n)],Ti,deg-1)   # T' T^{-1}
for k in range(deg):
    print('k',k,'right-log-deriv nnz',sum(1 for x in Gr[k] if x!=0),'left nnz',sum(1 for x in Gl[k] if x!=0), 'A1-sector nnz', sum(1 for x in A if x!=0) if k==0 else '')
pickle.dump((Gr,Gl),open('jet_basis_sym_G%d.pkl'%L,'wb'))
