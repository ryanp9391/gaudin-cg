from jet_basis_sym import *
import sys
from gaudin_chain import kron
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
R=gk.one(lam); d=R.d; I=SM.eye(d)
def dense(M):
    X=sp.zeros(d*d,d*d)
    for i,row in M.r.items():
        for c,v in row.items(): X[i,c]=sp.Rational(v.numerator,v.denominator)
    return X
def C(i,j,swap=False):
    tot=SM(d*d)
    for k in (1,2,3):
        tot=tot+(kron(R.E[(k,j)],R.E[(i,k)]) if swap else kron(R.E[(i,k)],R.E[(k,j)]))
    return dense(tot)
ops={}
for (i,j) in ((2,1),(3,2),(3,1)):
    ops['C%d%d'%(i,j)]=C(i,j); ops["C'%d%d"%(i,j)]=C(i,j,True)
    ops['E1_%d%d'%(i,j)]=dense(kron(R.E[(i,j)],I)); ops['E2_%d%d'%(i,j)]=dense(kron(I,R.E[(i,j)]))
names=list(ops); 
Mx=sp.Matrix([[ops[n][r,c] for n in names] for r in range(d*d) for c in range(d*d)])
b=sp.Matrix([A[r,c] for r in range(d*d) for c in range(d*d)])
sol=sp.linsolve((Mx,b),*sp.symbols('c0:%d'%len(names)))
print(names); print(sol)
