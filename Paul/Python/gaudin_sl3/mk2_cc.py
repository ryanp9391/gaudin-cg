import sympy as sp, master_kernel2 as m
from mk2_coord import *
K=m.K_; Pi=m.Pi_
e1,e2,h1,h2=E1,E2,H1,H2
C=sp.cancel((-a*ap*(a-ap)*(e1-h1)*Pi + K*((a-ap)*(a*ap*(e1*h1-e2-h2)-a*e1-ap*h1+2))))
Cc=tocoord(C)
Dm=ap-a
p=sp.Poly(sp.expand(Cc*Dm**2/(a*ap)**0),x1,x2,y1,y2)
for mon,co in sorted(p.terms()):
    print(mon,sp.factor(co))
print('K',sp.factor(tocoord(K)))
