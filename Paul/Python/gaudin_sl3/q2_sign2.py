import sys,glob,re; sys.path.insert(0,'.')
import gaudin_norm_gen as g, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps=80
tot=0;bad=0
for lam in [(1,0,0),(2,0,0),(3,0,0),(1,1,0),(2,2,0),(3,3,0),(2,1,0),(3,1,0),(3,2,0),(4,1,0),(4,2,0),(4,3,0),(4,4,0),(4,0,0),(5,1,0)]:
    tags=['P','Pp','X1'] if lam[0]<=3 else ['P']
    nb=0;n=0
    for tag in tags:
        recs=g.table(lam,tag)
        for r in recs:
            N,J=g.unpack(r); dj=mp.det(J); n+=1
            if (dj>0)!=(r['sector'][1]%2==0): nb+=1
    tot+=n;bad+=nb
    print(lam,'states tested',n,'violations of sign(detJ6)=(-1)^mu2:',nb,flush=True)
print('total',tot,'bad',bad)
