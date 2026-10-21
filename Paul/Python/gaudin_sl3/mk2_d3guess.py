import sympy as sp, pickle
from mk2_expfit import *
D2=getD2(); D3=sym_from_terms(pickle.load(open('master_kernel2_D3.pkl','rb')))
x=D1
def br(f,g): return sum(sp.diff(f,p)*sp.diff(g,q) for p,q in zip(PV,V))
t1=br(D2,x); t2=br(x,D2)
t3=sum(sp.diff(x,p1,p2)*sp.diff(x,q1)*sp.diff(x,q2) for p1,q1 in zip(PV,V) for p2,q2 in zip(PV,V))
t4=sum(sp.diff(x,p1)*sp.diff(x,p2)*sp.diff(x,q1,q2) for p1,q1 in zip(PV,V) for p2,q2 in zip(PV,V))
al,be,ga,de=sp.symbols('al be ga de')
expr=sp.expand(D3-(al*t1+be*t2+ga*t3+de*t4))
eqs=sp.Poly(expr,*(V+PV)).coeffs()
print(sp.solve(eqs,[al,be,ga,de],dict=True))
print('t1-t2-t3+t4 ==0 ?', sp.expand(t1-t2-t3+t4)==0)
print('D3 == (t1+t2)/2 ?', sp.expand(D3-(t1+t2)/2)==0)
print('D2 == br(x,x)', sp.expand(D2-br(x,x))==0)
