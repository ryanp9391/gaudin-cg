import pickle, sympy as sp, sys
a,ap,t,z=sp.symbols('a ap t z')
l=int(sys.argv[1])
res=pickle.load(open('mk_Qu_%d.pkl'%l,'rb'))
for m in sorted(res):
    c=res[m]
    # substitute t = z (ap-a)/(a ap)
    cz=sp.factor(sp.simplify(c.subs(t,z*(ap-a)/(a*ap))))
    print(m,cz)
