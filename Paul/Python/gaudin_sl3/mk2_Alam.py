import sys, pickle, sympy as sp
from sympy import Rational as R
import master_kernel2 as m
x1,x2,y1,y2=sp.symbols('x1 x2 y1 y2')
a=R(1,3); ap=R(1,5); zz=a*ap/(ap-a)
sE=sp.solve([1-a*m.E1+a**2*m.E2-x1,1-ap*m.E1+ap**2*m.E2-x2],[m.E1,m.E2],dict=True)[0]
sH=sp.solve([1-a*m.H1+a**2*m.H2-y1,1-ap*m.H1+ap**2*m.H2-y2],[m.H1,m.H2],dict=True)[0]
sub={**sE,**sH}
out={}
for L in range(1,int(sys.argv[1])+1):
    Q=m.symbol_orders_exact(L)
    co=[]
    for n in sorted(Q):
        q=sp.expand(Q[n].subs({m.a_:a,m.ap_:ap}).subs(sub,simultaneous=True))
        p=sp.Poly(q,x1,x2,y1,y2)
        c=p.coeff_monomial(x1**L*y2**L)
        co.append(sp.nsimplify(c/zz**n))
    out[L]=co
    print(L,co,flush=True)
pickle.dump(out,open('master_kernel2_Alam.pkl','wb'))
