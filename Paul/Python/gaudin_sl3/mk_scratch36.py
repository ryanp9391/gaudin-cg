import pickle, sympy as sp, sys
exec(open('mk_scratch35.py').read().split('# span test')[0])
Om=[(l-j)**2-j for j in range(l+1)]
def proj(f,j):
    # f homogeneous bidegree (l,l) polynomial in X,Y ; apply prod_{k!=j}(Omega-Om_k)/(Om_j-Om_k)
    g=f
    for k in range(l+1):
        if k==j: continue
        g=sp.expand((Omega(g)-Om[k]*g)/(Om[j]-Om[k]))
    return g
# homogenise Q_n: bidegree (l,l) in (X,Y): each monomial e1^a e2^b -> x0^(l-a-b) ...
def homog(f):
    P=sp.Poly(f,e1,e2,h1,h2); out=0
    for (a1,a2,b1,b2),c in P.terms():
        out+=c*x0**(l-a1-a2)*e1**a1*e2**a2*y0**(l-b1-b2)*h1**b1*h2**b2
    return out
Qn=[homog(q) for q in Qn]
mons=None
for j in range(l+1):
    vs=[proj(q,j) for q in Qn]
    allm=sorted({m for b in vs for m in sp.Poly(b,x0,e1,e2,y0,h1,h2).monoms()})
    M=sp.Matrix([[sp.Poly(b,x0,e1,e2,y0,h1,h2).coeff_monomial(m) for m in allm] for b in vs])
    print('component j=',j,'Omega eig',Om[j],'rank of span over orders:',M.rank(),' (orders',len(vs),')')
