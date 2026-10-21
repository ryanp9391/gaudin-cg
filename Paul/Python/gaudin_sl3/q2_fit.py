import sys,time; sys.path.insert(0,'.')
import gaudin_norm_gen as g, mpmath as mp, numpy as np
from fractions import Fraction as Fr
mp.mp.dps=80
TAGS=['P','Pp','X1','X2','X3','q2a','q2b','q2c','q2d','q2e','q2f','q2g','q2h']
def c0data(lam,tag,th=g.TH0, exps=None):
    recs=g.table(lam,tag,th)
    z=[g.mpf(x) for x in g.ZET[tag]]
    d12,d13,d23=abs(z[0]-z[1]),abs(z[0]-z[2]),abs(z[1]-z[2])
    vals=[];rows=[]
    for r in recs:
        N,J=g.unpack(r); m=r['sector']; dj=mp.det(J)
        vals.append(N/abs(dj)/(d23**m[0]*d13**m[1]*d12**m[2]))
        rows.append((m,mp.log(N/abs(dj)),dj))
    return z,(d12,d13,d23),vals,rows
if __name__=='__main__':
    lams=[(3,1,0),(3,2,0)]
    for lam in lams:
        print('=====',lam)
        Cs=[];Ds=[]
        for tag in TAGS:
            try: z,d,vals,rows=c0data(lam,tag)
            except Exception as e: print(tag,'missing',e); continue
            spread=max(abs(v/vals[0]-1) for v in vals)
            D=d[0]*d[1]*d[2]
            Cs.append(vals[0]); Ds.append(D)
            # regression of log(N/|J|) on mu
            A=np.array([[1.0,m[0],m[1],m[2]] for m,l,dj in rows]); y=np.array([float(l) for m,l,dj in rows])
            c,*_=np.linalg.lstsq(A,y,rcond=None)
            L=[float(mp.log(x)) for x in (d[2],d[1],d[0])]
            print(tag,'n=%d c0=%s spread=%s'%(len(rows),mp.nstr(vals[0],15),mp.nstr(spread,3)),'regress coef diff (mu1-mu3,mu2-mu3):',c[1]-c[3],c[2]-c[3],'vs',L[0]-L[2],L[1]-L[2], '#detJ>0',sum(1 for m,l,dj in rows if dj>0))
        # fit c0 = C * D^(-c)
        lx=np.array([float(mp.log(D)) for D in Ds]); ly=np.array([float(mp.log(C)) for C in Cs])
        A=np.vstack([np.ones_like(lx),lx]).T
        co,*_=np.linalg.lstsq(A,ly,rcond=None)
        print('fit log c0 = logC + s*log D: C=%.10f, s=%.8f'%(np.exp(co[0]),co[1]))
