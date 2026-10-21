import sys, pickle
import sympy as sp
l1=2
Q=pickle.load(open('mk_symbol_%d.pkl'%l1,'rb'))
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
f=lambda c:(1-c*x1)*(1-c*x2); g=lambda c:(1-c*y1)*(1-c*y2)
basis=[f(a)**2*g(ap)**2, f(a)*f(ap)*g(a)*g(ap), f(ap)**2*g(a)**2]
c=sp.symbols('c0:3')
expr=sp.expand(Q-sum(ci*b for ci,b in zip(c,basis)))
eqs=sp.Poly(expr,x1,x2,y1,y2).coeffs()
sol=sp.solve(eqs,c,dict=True)
print(sol)
