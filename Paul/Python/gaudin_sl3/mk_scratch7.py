import sys, pickle
import sympy as sp
l1=int(sys.argv[1])
Q=pickle.load(open('mk_symbol_%d.pkl'%l1,'rb'))
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
P=sp.Poly(Q,t)
d=P.degree()
print('deg',d)
for n in range(d,d-2,-1):
    c=sp.factor(P.coeff_monomial(t**n))
    s=str(c)
    print(n, s if len(s)<1500 else s[:1500])
