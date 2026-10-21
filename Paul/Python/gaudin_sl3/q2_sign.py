import sys; sys.path.insert(0,'.')
import gaudin_norm_gen as g, mpmath as mp
from collections import defaultdict
mp.mp.dps=80
for lam in [(2,1,0),(3,1,0),(3,2,0)]:
  for tag in ['P','Pp','X1']:
    recs=g.table(lam,tag)
    agree=0; rat=[]; pos=0; byd=defaultdict(list); bysec=defaultdict(list)
    for i,r in enumerate(recs):
        N,J=g.unpack(r)
        dj=mp.det(J)
        B1=mp.det(J[0:3,0:3]); B2=mp.det(J[3:6,3:6])
        agree+= (dj*B1*B2>0); rat.append(dj/(B1*B2))
        pos+= dj>0
        bysec[r['sector']].append(int(dj>0))
        # sign of B1*B2, of off-diagonal blocks
    # sector purity
    pure=sum(1 for s,v in bysec.items() if len(set(v))==1)
    print(lam,tag,'n=%d  #detJ>0=%d  sign(detJ)=sign(B1 B2) in %d/%d ; detJ/(B1B2) range [%s, %s] ; sectors with constant sign: %d/%d'%(len(recs),pos,agree,len(recs),mp.nstr(min(rat),5),mp.nstr(max(rat),5),pure,len(bysec)))
    if tag=='P': print('   sign by sector (mu: +/- counts):',{s:(sum(v),len(v)-sum(v)) for s,v in sorted(bysec.items())})
