import pickle, sympy as sp, master_kernel2 as m
from mk2_coord import *
W3=pickle.load(open('master_kernel2_W3.pkl','rb'))   # = W3 * (ap-a)^3 as polynomial in lam, E,H,a,ap
p=sp.Poly(W3,lam)
cs={k:p.coeff_monomial(lam**k) for k in (2,3,4)}
# rewrite lam^2 d2 + lam^3 d3 + lam^4 d4 in basis lam^2(lam-1)(lam-2), lam^2(lam-1), lam^2:
# lam^2(lam-1)(lam-2)=lam^4-3lam^3+2lam^2 ; lam^2(lam-1)=lam^3-lam^2
d4=cs[4]; d3=cs[3]+3*d4; d2=cs[2]-2*d4+ (d3)  # check below
# solve generally: W = A*lam^2(lam-1)(lam-2)+B*lam^2(lam-1)+C*lam^2
A_=cs[4]; B_=cs[3]+3*cs[4]; C_=cs[2]+B_-2*A_
chk=sp.expand(A_*(lam**4-3*lam**3+2*lam**2)+B_*(lam**3-lam**2)+C_*lam**2-sum(cs[k]*lam**k for k in cs))
print('basis check',chk==0)
Dm=ap-a
res={}
for nm,e in (('A',A_),('B',B_),('C',C_)):
    c=tocoord(e)
    pp=sp.Poly(sp.expand(c),x1,x2,y1,y2)
    res[nm]=pp
    print(nm,len(pp.terms()),'terms')
pickle.dump({k:v.as_expr() for k,v in res.items()},open('master_kernel2_W3coord.pkl','wb'))
for nm in 'ABC':
    print('=====',nm)
    for mon,co in sorted(res[nm].terms()):
        print(mon,sp.factor(co))
