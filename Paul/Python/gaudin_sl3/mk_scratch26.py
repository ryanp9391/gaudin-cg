import pickle, sympy as sp, sys
a,ap,t,z=sp.symbols('a ap t z')
w,wp,nu=sp.symbols('w wp nu')
l=int(sys.argv[1])
res=pickle.load(open('mk_Qu_%d.pkl'%l,'rb'))
out={}
for m in sorted(res):
    c=res[m]
    # u=U/w^2 -> u^m0 ... factor w^{-2 m0} wp^{-2 m1}... u=U/w^2, up=U'/wp^2, v=V/wp^2, vp=V'/w^2
    fac=(w**2)**(-m[0])*(wp**2)**(-m[1])*(wp**2)**(-m[2])*(w**2)**(-m[3])
    cc=c.subs({a:1/w,ap:1/wp,t:z*(w-wp)},simultaneous=True)  # t = z (w-wp)
    cc=sp.factor(sp.simplify(cc*fac*(w*wp)**(2*l)))
    out[m]=cc
    print(m,cc)
pickle.dump(out,open('mk_Qw_%d.pkl'%l,'wb'))
