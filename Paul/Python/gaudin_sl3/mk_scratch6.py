import sys, pickle
import sympy as sp
l1=int(sys.argv[1])
Q=pickle.load(open('mk_symbol_%d.pkl'%l1,'rb'))
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
sub={x1:x1/(1+a*x1),x2:x2/(1+a*x2),y1:y1/(1+ap*y1),y2:y2/(1+ap*y2)}
G=Q.subs(sub,simultaneous=True)
G=G*((1+a*x1)*(1+a*x2))**l1*((1+ap*y1)*(1+ap*y2))**l1
G=sp.factor(sp.cancel(sp.together(G)))
Gp=sp.Poly(sp.expand(G),t)
for n,c in sorted(((m[0],c) for m,c in Gp.terms())):
    print('t^%d:'%n, sp.factor(c))
