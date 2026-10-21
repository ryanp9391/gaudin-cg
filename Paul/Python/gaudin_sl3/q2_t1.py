import sys,time; sys.path.insert(0,'.')
import gaudin_norm_gen as g, mpmath as mp
for lam in [(2,1,0),(3,1,0),(3,2,0)]:
    R=g.rep(lam)
    w4=R.w4['gmet1']
    gm=[R.gd[R.w2m[i]] for i in range(R.d)]
    print(lam,'metric matches W4 GMET1:', [int(x) for x in gm]==[int(x) for x in w4], flush=True)
for lam in [(2,1,0)]:
    recs=g.table(lam,'P')
    zeta=g.ZET['P']; z=[g.mpf(x) for x in zeta]
    d12,d13,d23=z[0]-z[1],z[0]-z[2],z[1]-z[2]
    vals=[]
    for r in recs:
        N,J=g.unpack(r); m=r['sector']
        vals.append(N/abs(mp.det(J))/(abs(d23)**m[0]*abs(d13)**m[1]*abs(d12)**m[2]))
    print(vals[0], max(abs(v/vals[0]-1) for v in vals), vals[0]*abs(d12*d13*d23)**6)
    print(max(abs(mp.mpf(q)) for r in recs for q in r['Fres']))
