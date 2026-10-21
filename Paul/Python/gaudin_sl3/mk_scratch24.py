import pickle, sympy as sp, sys
from mk_scratch20 import to_e
a,ap,t,e1,e2,h1,h2=sp.symbols('a ap t e1 e2 h1 h2')
u,up,v,vp=sp.symbols('u up v vp')
s1=sp.solve([1-a*e1+a**2*e2-u,1-ap*e1+ap**2*e2-up],[e1,e2],dict=True)[0]
s2=sp.solve([1-ap*h1+ap**2*h2-v,1-a*h1+a**2*h2-vp],[h1,h2],dict=True)[0]
sub={**s1,**s2}
for l in [int(c) for c in sys.argv[1:]]:
    Q=pickle.load(open('mk_symbol_%d.pkl'%l,'rb'))
    Qe=to_e(Q)
    Qu=sp.together(Qe.subs(sub,simultaneous=True))
    num,den=sp.fraction(Qu)
    print(l,'den',sp.factor(den))
    Qp=sp.Poly(sp.expand(num),u,up,v,vp)
    res={}
    for m,c in Qp.terms():
        res[m]=sp.factor(c/den)
    pickle.dump(res,open('mk_Qu_%d.pkl'%l,'wb'))
    print(l,'nterms',len(res), 'monomials(u,up,v,vp):',sorted(res)[:40])
