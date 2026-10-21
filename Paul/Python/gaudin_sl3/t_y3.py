from jet_basis_sym import *
res=pickle.load(open('jet_basis_sym_Y200.pkl','rb'))
# zeta independence
for k,(th,Y) in res.items():
    ref=res.get(('P',k[1]))
    if ref and k[0]!='P': print(k, Y==ref[1])
# t-fit
pts=[(1/(th[0]-th[1]),Y) for (tag,key),(th,Y) in res.items() if tag=='P']
d={}
for tt,Y in pts:
    if tt in d: print('same t, equal Y:', d[tt]==Y)
    d[tt]=Y
pts=list(d.items())
print(len(pts), 'distinct t', len(set(p[0] for p in pts)))
n=len(pts[0][1]); m=len(pts[0][1][0])
import sympy as sp
t=sp.symbols('t')
use=pts[:7]; chk=pts[7:]
Ypoly=[[None]*m for _ in range(n)]
maxdeg=0
for i in range(n):
    for j in range(m):
        xs=[sp.Rational(p[0].numerator,p[0].denominator) for p in use]
        ys=[sp.Rational(p[1][i][j].numerator,p[1][i][j].denominator) for p in use]
        poly=sp.interpolate(list(zip(xs,ys)),t)
        poly=sp.expand(poly)
        for p in chk:
            assert poly.subs(t,sp.Rational(p[0].numerator,p[0].denominator))==sp.Rational(p[1][i][j].numerator,p[1][i][j].denominator),(i,j)
        Ypoly[i][j]=poly
        if poly!=0: maxdeg=max(maxdeg,sp.degree(poly,t))
print('maxdeg',maxdeg)
pickle.dump(Ypoly,open('jet_basis_sym_Y200_poly.pkl','wb'))
