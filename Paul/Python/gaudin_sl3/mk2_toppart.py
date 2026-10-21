import pickle, sympy as sp
from mk2_opQ import *
D=Dlist(4)
def top(expr,deg):
    P=sp.Poly(sp.expand(expr),*PV)
    return sum(co*sp.prod([p**k for p,k in zip(PV,mon)]) for mon,co in P.terms() if sum(mon)==deg)
D4top=top(D[4],5); D4low=D[4]-D4top
rec4=sp.expand(sp.Rational(1,3)*(brk(D[1],D[3])+brk(D[2],D[2])+brk(D[3],D[1])))
print('D4 top-degree part == recursion (1/3)sum <Di,Dj>:',sp.expand(D4top-rec4)==0)
# try general combination of the first-order brackets
al,be,ga=sp.symbols('al be ga')
expr=sp.expand(D4top-(al*brk(D[1],D[3])+be*brk(D[2],D[2])+ga*brk(D[3],D[1])))
eqs=sp.Poly(expr,*(V+PV)).coeffs()
print(sp.solve(eqs,[al,be,ga],dict=True))
def br2(f,g):
    return sp.expand(sp.Rational(1,2)*sum(sp.diff(f,p1,p2)*sp.diff(g,q1,q2) for p1,q1 in zip(PV,V) for p2,q2 in zip(PV,V)))
b13,b22,b31=sp.symbols('b13 b22 b31')
expr=sp.expand(D4low-(b13*br2(D[1],D[3])+b22*br2(D[2],D[2])+b31*br2(D[3],D[1])))
print('second order fit:',sp.solve(sp.Poly(expr,*(V+PV)).coeffs(),[b13,b22,b31],dict=True))
print('br2(D1,D1) alone', br2(D[1],D[1]))
print('D4low',sp.factor(D4low))
print('br2(D1,D3)',sp.factor(br2(D[1],D[3])))
print('br2(D2,D2)',sp.factor(br2(D[2],D[2])))
