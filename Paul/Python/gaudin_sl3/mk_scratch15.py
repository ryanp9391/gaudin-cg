import sys, pickle, time
import sympy as sp
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
vars_=(x1,x2,y1,y2,a,ap)
l=4
Q=pickle.load(open('mk_symbol_4.pkl','rb'))
Pq=sp.Poly(Q,t)
lin=[sp.Poly(1-a*x1,*vars_),sp.Poly(1-a*x2,*vars_),sp.Poly(1-ap*y1,*vars_),sp.Poly(1-ap*y2,*vars_)]
for n in (1,2,3):
    c=sp.Poly(Pq.coeff_monomial(t**n),*vars_)
    t0=time.time()
    for rep in range(l-n):
        for f in lin:
            q,r=sp.div(c,f)
            assert r.is_zero,(n,rep)
            c=q
    pickle.dump(c.as_expr(),open('mk_Pn_%d_%d.pkl'%(n,l),'wb'))
    print(n,'ok',time.time()-t0,flush=True)
