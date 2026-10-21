import sys, itertools, pickle
import sympy as sp, mpmath as mp
mp.mp.dps=60
l1=2
Q=pickle.load(open('mk_symbol_%d.pkl'%l1,'rb'))
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
vals={a:sp.Rational(2,5),ap:sp.Rational(3,7),t:sp.Rational(1,3),x1:sp.Rational(1,5),x2:sp.Rational(-2,3),y1:sp.Rational(3,4),y2:sp.Rational(2,9)}
exact=mp.mpf(sp.N(Q.subs(vals),50))
def cand(eps, rule):
    av=[mp.mpf(2)/5+ k*eps for k in (0,1)]; apv=[mp.mpf(3)/7+k*eps*mp.mpf('1.37') for k in (0,1)]
    tt=mp.mpf(1)/3
    X=[mp.mpf(1)/5,mp.mpf(-2)/3]; Y=[mp.mpf(3)/4,mp.mpf(2)/9]
    pts=[('s',i,1/av[i]) for i in range(2)]+[('p',i,1/apv[i]) for i in range(2)]
    tot=0
    for S in itertools.combinations(range(4),2):
        w=1
        for u in range(2):
            for v in range(2,4):
                z=tt/(1/av[u]-1/apv[v-2])
                inS_u=u in S; inS_v=v in S
                if inS_u and not inS_v: w*=(1+z)
                elif (not inS_u) and inS_v: w*=(-z)
        F=1
        for u in S:
            c=av[u] if u<2 else apv[u-2]
            for xi in X: F*=(1-c*xi)
        G=1
        for u in range(4):
            if u not in S:
                c=av[u] if u<2 else apv[u-2]
                for yy in Y: G*=(1-c*yy)
        tot+=w*F*G
    return tot
for eps in (mp.mpf('1e-6'),mp.mpf('1e-8')):
    print(eps, cand(eps,0), exact)
