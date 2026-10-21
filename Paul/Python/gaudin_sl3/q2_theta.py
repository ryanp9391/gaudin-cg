import sys,time,glob,re; sys.path.insert(0,'.')
import gaudin_norm_gen as g, mpmath as mp
from fractions import Fraction as Fr
from math import factorial as f
mp.mp.dps=80
def Cpred(l1,l2):
    return Fr(f(l2)*f(l1-l2)*f(l1+1),(l2+1)*(l1-l2+1))**2
def run(lam,tag,th=g.TH0):
    recs=g.table(lam,tag,th)
    z=[g.mpf(x) for x in g.ZET[tag]]
    d12,d13,d23=abs(z[0]-z[1]),abs(z[0]-z[2]),abs(z[1]-z[2]); D=d12*d13*d23
    c=2*(lam[0]+lam[1])
    vals=[]
    for r in recs:
        N,J=g.unpack(r); m=r['sector']
        vals.append(N/abs(mp.det(J))/(d23**m[0]*d13**m[1]*d12**m[2]))
    C=mp.mpf(Cpred(lam[0],lam[1]).numerator)/Cpred(lam[0],lam[1]).denominator
    # full formula residual: N/(C |J| prod d^(mu-c)) - 1
    res=max(abs(v*D**c/C-1) for v in vals)
    sp=max(abs(v/vals[0]-1) for v in vals)
    return len(recs),sp,res
print('## (C) theta dependence at tag P')
for lam in [(3,2,0)]:
    files=sorted(glob.glob('general_kernel_pt_%d%d%d_P_*_80.pkl'%lam))
    ths=[]
    for fn in files:
        m=re.match(r'.*_P_(-?\d+)\.(\d+)_(-?\d+)\.(\d+)_80\.pkl',fn)
        ths.append((Fr(int(m.group(1)),int(m.group(2))),Fr(int(m.group(3)),int(m.group(4)))))
    for th in ths:
        n,a,b=run(lam,'P',th)
        print(lam,'theta=(%s,%s)'%th,'spread %s residual %s'%(mp.nstr(a,2),mp.nstr(b,2)),flush=True)
