import sys,time; sys.path.insert(0,'.')
import gaudin_norm_gen as g, mpmath as mp
for lam in [(3,3,0),(4,2,0),(4,3,0),(4,4,0),(5,1,0),(2,0,0),(3,0,0),(4,0,0),(1,0,0)]:
    t=time.time()
    try:
        recs=g.table(lam,'P')
    except Exception as e:
        print(lam,'FAIL',repr(e)[:100],flush=True); continue
    z=[g.mpf(x) for x in g.ZET['P']]
    d12,d13,d23=abs(z[0]-z[1]),abs(z[0]-z[2]),abs(z[1]-z[2]); D=d12*d13*d23
    vals=[]
    for r in recs:
        N,J=g.unpack(r); m=r['sector']
        vals.append(N/abs(mp.det(J))/(d23**m[0]*d13**m[1]*d12**m[2]))
    c=2*(lam[0]+lam[1])
    print(lam,len(recs),'c0 spread',mp.nstr(max(abs(v/vals[0]-1) for v in vals),3),'C(if c=2|L|)=',mp.nstr(vals[0]*D**c,20),'%.0fs'%(time.time()-t),flush=True)
