import pickle, sympy as sp, master_kernel2 as m
from mk2_coord import *
c=pickle.load(open('master_kernel2_P3coef.pkl','rb'))
Dm=ap-a
nonmono=0; tot=0
out=[]
for k in range(4):
    e=tocoord(c[k]*Dm**3)
    p=sp.Poly(e,x1,x2,y1,y2)
    for mon,co in p.terms():
        tot+=1
        f=sp.factor(co)
        if len(sp.Add.make_args(sp.expand(co)))>1: nonmono+=1
    out.append(p)
    print('lam^%d coefficient: %d terms'%(k,len(p.terms())))
print('non-monomial coefficient count',nonmono,'of',tot)
pickle.dump(out,open('master_kernel2_P3coord.pkl','wb'))
p=out[3]
for mon,co in p.terms(): print(mon,sp.factor(co))
