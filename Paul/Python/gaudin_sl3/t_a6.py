from jet_basis_sym import *
import sys, itertools
from gaudin_chain import kron
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
R=gk.one(lam); d=R.d; I=SM.eye(d)
alpha,s,hw,Eg=dictionary(L)
def dense(M):
    X=sp.zeros(d*d,d*d)
    for i,row in M.r.items():
        for c,v in row.items(): X[i,c]=sp.Rational(v.numerator,v.denominator)
    return X
E={k:dense_ for k,dense_ in R.E.items()}
A1=A
A2=coef[0].inv()*coef[2]
B=A2-A1*A1/2
w3=[tuple(alpha[p][k]+alpha[q][k] for k in range(3)) for p in range(d) for q in range(d)]
shifts=lambda X: set(tuple(w3[r][k]-w3[c][k] for k in range(3)) for r in range(d*d) for c in range(d*d) if X[r,c]!=0)
print('A2 shifts',shifts(A2)); print('B shifts',shifts(B))
print('B nnz',sum(1 for x in B if x!=0),'A2 nnz',sum(1 for x in A2 if x!=0))
