import sys, pickle
import sympy as sp
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
f=lambda c:(1-c*x1)*(1-c*x2); g=lambda c:(1-c*y1)*(1-c*y2)
Pi=f(a)*g(ap)
Q={l:pickle.load(open("mk_symbol_%d.pkl"%l,"rb")) for l in (1,2,3)}
for n in (1,2,3):
    for l in (1,2,3):
        c=sp.Poly(Q[l],t).coeff_monomial(t**n)
        k=l-n
        if k>=0:
            q,r=sp.div(sp.Poly(c,x1,x2,y1,y2),sp.Poly(sp.expand(Pi**k),x1,x2,y1,y2))
            print('n',n,'l',l,'div by Pi^%d'%k,'rem zero' if r.is_zero else 'rem nonzero', flush=True)
        else:
            q,r=sp.div(sp.Poly(c,x1,x2,y1,y2),sp.Poly(1,x1,x2,y1,y2)); print('n',n,'l',l,'skip')
        if l>=1 and k>=0 and r.is_zero:
            pickle.dump(q.as_expr(),open('mk_Pn_%d_%d.pkl'%(n,l),'wb'))
