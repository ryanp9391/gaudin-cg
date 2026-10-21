import pickle, sympy as sp
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
vars_=(x1,x2,y1,y2,a,ap)
P=lambda n,l: sp.Poly(pickle.load(open('mk_Pn_%d_%d.pkl'%(n,l),'rb')),*vars_)
for n,ls in ((1,(1,2,3,4)),(2,(2,3,4)),(3,(3,4))):
    ps=[P(n,l) for l in ls]
    mons=sorted({m for p in ps for m in p.monoms()})
    M=sp.Matrix([[p.coeff_monomial(m) for m in mons] for p in ps])
    print(n,ls,'rank',M.rank(),'nmons',len(mons))
    print(' total degrees in (x,y):',sorted({sum(m[:4]) for m in mons}),' a-degs',sorted({m[4] for m in mons}),'ap',sorted({m[5] for m in mons}))
