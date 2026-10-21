import pickle, sympy as sp, master_kernel2 as m
lam=sp.Symbol('lam')
a,ap=m.a_,m.ap_
x1,x2,y1,y2=sp.symbols('x1 x2 y1 y2')
E1,E2,H1,H2=m.E1,m.E2,m.H1,m.H2
sE=sp.solve([1-a*E1+a**2*E2-x1,1-ap*E1+ap**2*E2-x2],[E1,E2],dict=True)[0]
sH=sp.solve([1-ap*H1+ap**2*H2-y2,1-a*H1+a**2*H2-y1],[H1,H2],dict=True)[0]
# note: x1=F_a(xi), x2=F_a'(xi), y1=G_a(eta), y2=G_a'(eta)
sub={**sE,**{H1:sH[H1],H2:sH[H2]}}
def tocoord(expr):
    return sp.expand(sp.cancel(sp.together(expr.subs(sub,simultaneous=True))))
if __name__=='__main__':
    z=sp.Symbol('z')
    c=pickle.load(open('master_kernel2_P3coef.pkl','rb'))
    K=m.K_
    for nm,e in (('K',K),):
        print(nm,tocoord(e))
    Pi=m.Pi_
    e1,e2,h1,h2=E1,E2,H1,H2
    A=sp.cancel(-a*ap*(a-ap)*(e1-h1)*Pi + K*(K/2 + (a-ap)*(a*ap*(e1*h1-e2-h2)-a*e1-ap*h1+2)))
    Ac=tocoord(A); print('A',sp.factor(Ac))
    B=tocoord(K**2/2)
    # test dependence: substitute t-scaling: z=a*ap/(ap-a); P_2/z^2
    P2c=tocoord(A+lam*K**2/2)
    zz=a*ap/(ap-a)
    print('A/z^2',sp.factor(sp.cancel(Ac/zz**2)))
