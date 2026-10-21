from jet_basis_sym import *
import sys
L=int(sys.argv[1])
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
Y0i=coef[0].inv(); n=int(round(coef[0].shape[0]**0.5))
Ms=[coef[k]*Y0i for k in range(len(coef))]
def schmidt(M):
    # M[(i j),(k l)] -> R[(i k),(j l)]
    R=sp.zeros(n*n,n*n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    R[i*n+k,j*n+l]=M[i*n+j,k*n+l]
    return R
Pm=sp.zeros(n*n,n*n)
for i in range(n):
    for j in range(n): Pm[i*n+j,j*n+i]=1
for k in range(1,len(Ms)):
    print('M%d operator-Schmidt rank'%k, schmidt(Ms[k]).rank(), ' with swap', schmidt(Ms[k]*Pm).rank())
