from jet_basis_sym import *
import sys
L=int(sys.argv[1])
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
Y0i=coef[0].inv()
Ms=[coef[k]*Y0i for k in range(len(coef))]
M1=Ms[1]
print('M2==M1^2/2', Ms[2]==M1*M1/2)
# try general: is M(t) = exp(t M1 + ...)?  compute log series
t=sp.symbols('t')
n=M1.shape[0]
# log(1+X) with X = sum t^k M_k : coefficients up to t^L
X=[sp.zeros(n,n)]+Ms[1:]
def mul(P,Q,deg):
    R=[sp.zeros(n,n) for _ in range(deg+1)]
    for i,p in enumerate(P):
        for j,q in enumerate(Q):
            if i+j<=deg: R[i+j]+=p*q
    return R
deg=L
logs=[sp.zeros(n,n) for _ in range(deg+1)]
pw=X
for k in range(1,deg+1):
    for i in range(deg+1): logs[i]+= sp.Rational((-1)**(k+1),k)*pw[i]
    pw=mul(pw,X,deg)
for k in range(1,deg+1):
    print('log coeff t^%d nnz'%k, sum(1 for x in logs[k] if x!=0), ' M_k nnz',sum(1 for x in Ms[k] if x!=0))
# M_k vs M1^k/k! ratio structure
for k in range(2,deg+1):
    Pk=M1**k/sp.factorial(k)
    print('k',k,'M_k - M1^k/k! nnz',sum(1 for x in (Ms[k]-Pk) if x!=0))
