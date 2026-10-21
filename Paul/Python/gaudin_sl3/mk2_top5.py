import pickle, sympy as sp
from mk2_opQ import *
D=Dlist(4)
D5=sym_from_terms(pickle.load(open('master_kernel2_D5_gen_L5.pkl','rb')))
def top(expr,deg):
    P=sp.Poly(sp.expand(expr),*PV)
    return sum(co*sp.prod([p**k for p,k in zip(PV,mon)]) for mon,co in P.terms() if sum(mon)==deg)
D4t=top(D[4],5)
rec5=sp.expand(sp.Rational(1,4)*(brk(D[1],D4t)+brk(D[2],D[3])+brk(D[3],D[2])+brk(D4t,D[1])))
D5top=top(D5,6)
print('D5 top == (1/4) sum <Di,Dj>_1 (tree):',sp.expand(D5top-rec5)==0)
D5low=sp.expand(D5-D5top)
print('D5 low-degree terms count',len(sp.Poly(D5low,*(V+PV)).terms()))
def br2(f,g):
    return sp.expand(sp.Rational(1,2)*sum(sp.diff(f,p1,p2)*sp.diff(g,q1,q2) for p1,q1 in zip(PV,V) for p2,q2 in zip(PV,V)))
D4low=D[4]-D4t
# candidate loop-1 sources: <D1,D4>_2,<D4,D1>_2,<D2,D3>_2,<D3,D2>_2 ; and <D_i,D_j>_1 with one loop part (D4low with D1)
cands={'b2_14':br2(D[1],D4t),'b2_41':br2(D4t,D[1]),'b2_23':br2(D[2],D[3]),'b2_32':br2(D[3],D[2]),
       'b1_1L4':brk(D[1],D4low),'b1_L41':brk(D4low,D[1])}
names=list(cands); cs=sp.symbols('k0:%d'%len(names))
expr=sp.expand(D5low-sum(c*cands[n] for c,n in zip(cs,names)))
sol=sp.solve(sp.Poly(expr,*(V+PV)).coeffs(),cs,dict=True)
print(names); print(sol)
