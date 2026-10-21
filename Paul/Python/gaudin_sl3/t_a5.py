from jet_basis_sym import *
import sys, itertools
from gaudin_chain import kron
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
R=gk.one(lam); d=R.d; I=SM.eye(d)
def dense(M):
    X=sp.zeros(d*d,d*d)
    for i,row in M.r.items():
        for c,v in row.items(): X[i,c]=sp.Rational(v.numerator,v.denominator)
    return X
E={k:v for k,v in R.E.items()}
def wsh(i,j):
    w=[0,0,0]; w[i-1]+=1; w[j-1]-=1; return w
targets=[(-1,1,0),(0,-1,1)]
for tgt in targets:
    ops={}
    for (i,j) in itertools.product((1,2,3),repeat=2):
        for (k,l) in itertools.product((1,2,3),repeat=2):
            w=tuple(a+b for a,b in zip(wsh(i,j),wsh(k,l)))
            if w==tgt: ops[('B',i,j,k,l)]=dense(kron(E[(i,j)],E[(k,l)]))
    for (i,j) in itertools.product((1,2,3),repeat=2):
        if tuple(wsh(i,j))==tgt:
            ops[('S1',i,j)]=dense(kron(E[(i,j)],I)); ops[('S2',i,j)]=dense(kron(I,E[(i,j)]))
    names=list(ops)
    print(tgt,len(names))
    # restrict A to that weight shift
    alpha,s,hw,Eg=dictionary(L)
    w3=[tuple(alpha[p][k]+alpha[q][k] for k in range(3)) for p in range(d) for q in range(d)]
    Asub=sp.zeros(d*d,d*d)
    for r in range(d*d):
        for c in range(d*d):
            if tuple(w3[r][k]-w3[c][k] for k in range(3))==tgt: Asub[r,c]=A[r,c]
    Mx=sp.Matrix([[ops[n][r,c] for n in names] for r in range(d*d) for c in range(d*d)])
    b=sp.Matrix([Asub[r,c] for r in range(d*d) for c in range(d*d)])
    sol=sp.linsolve((Mx,b),*sp.symbols('c0:%d'%len(names)))
    print(names)
    print(sol)
