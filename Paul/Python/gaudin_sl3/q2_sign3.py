import sys; sys.path.insert(0,'.')
import gaudin_norm_gen as g, mpmath as mp
mp.mp.dps=80
sg=lambda x:1 if x>0 else -1
for lam in [(2,1,0),(3,1,0),(3,2,0),(2,2,0),(3,0,0)]:
  for tag in ['P','Pp','X1','X2','X3','q2a','q2b','q2c','q2d']:
    recs=g.table(lam,tag); z=g.ZET[tag]
    s12,s13,s23=sg(z[0]-z[1]),sg(z[0]-z[2]),sg(z[1]-z[2])
    pat=set(); pat2=set()
    for r in recs:
        N,J=g.unpack(r); m=r['sector']; dj=mp.det(J)
        base=sg(dj)*(-1)**m[1]
        pat.add(base)                                  # relative to (-1)^mu2
        pat2.add(base*s12**m[2]*s13**m[1]*s23**m[0])   # with sign(zeta differences)^mu
    print(lam,tag,'order signs (12,13,23)=',(s12,s13,s23),'sign(detJ)*(-1)^mu2 takes',sorted(pat),' with zeta-diff signs ^mu:',sorted(pat2))
