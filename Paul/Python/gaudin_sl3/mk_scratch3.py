import sys
sys.path.insert(0,'.')
from master_kernel import *
import sympy as sp
l1=int(sys.argv[1])
lam=(l1,0,0)
T,rep=G.load_rep(lam)
prs=G.pairs_of(lam)
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
S=0
for i,p in enumerate(prs):
    for j,q in enumerate(prs):
        e=sum(sp.Rational(c.numerator,c.denominator)*t**n*a**pp*ap**qq for (n,pp,qq),c in T[(i,j)].items())
        S+=e*(x1**p[0]*x2**p[1]-x1**p[1]*x2**p[0])*(y1**q[0]*y2**q[1]-y1**q[1]*y2**q[0])
S=sp.expand(S)
Q=sp.cancel(S/((x2-x1)*(y2-y1)))
Q=sp.expand(Q)
import pickle
pickle.dump(Q,open('mk_symbol_%d.pkl'%l1,'wb'))
print(sp.factor(Q.subs(t,0)))
Q1=sp.expand(Q.coeff(t,1)); print('t1:',sp.factor(Q1))
print('t2:',sp.factor(sp.expand(Q.coeff(t,2))))
