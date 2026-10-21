import sys,time,pickle; sys.path.insert(0,'.')
import gaudin_norm_gen as g, nu_general as NG, mpmath as mp
from fractions import Fraction as Fr
lam=tuple(int(c) for c in sys.argv[1]); which=sys.argv[2] if len(sys.argv)>2 else ''   # '' = phi, '2' = phi'
phi=NG.PHI if which=='' else NG.PHI2
NPTS=24
mp.mp.dps=200
R=NG.rep_of(lam)
pts=NG.load_pts(lam,which,phi,NPTS)
res=pickle.load(open('nu_general_res_%d%d%d%s.pkl'%(lam+(which,)),'rb'))
l1,l2=lam[0],lam[1]
F=NG.fitter(NPTS)
c=2*(l1+l2)
Cc=mp.mpf(g.Cpred(l1,l2).numerator)/g.Cpred(l1,l2).denominator if hasattr(g,'Cpred') else None
from math import factorial as f
Cf=Fr(f(l2)*f(l1-l2)*f(l1+1),(l2+1)*(l1-l2+1))**2
Cc=mp.mpf(Cf.numerator)/Cf.denominator
out=[]
for bi,rec in enumerate(res['br']):
    t0=time.time()
    kB=rec['kB']; XB=rec['XB']; mu=rec['mu']
    p=2*(2*l2+kB)
    Nf=[];Nd=[]
    for pt in pts:
        b=pt['branches'][bi]
        zeta=pt['zeta']
        lau=list(NG.laurent_of(b['lam2'],b['lam3'],0))+list(NG.laurent_of(b['lam2'],b['lam3'],1))
        Fv,J=g.jac_at(lam,lau,zeta,NG.THETA)
        z=[g.mpf(x) for x in zeta]
        d12,d13,d23=abs(z[0]-z[1]),abs(z[0]-z[2]),abs(z[1]-z[2])
        Nf.append(Cc*abs(mp.det(J))*d23**(mu[0]-c)*d13**(mu[1]-c)*d12**(mu[2]-c))
        v=R.full(b)
        n2=mp.fsum(R.gdm[i//R.d]*R.gdm[i%R.d]*v[i]**2 for i in range(R.d**2) if v[i]!=0)
        xx=R.xcov2(zeta)
        ox=mp.fsum(xx[i]*v[i] for i in range(R.d**2) if v[i]!=0)
        Nd.append(n2/ox**2)
    dev=max(abs(a/b-1) for a,b in zip(Nf,Nd))
    eps=[pt['eps'] for pt in pts]
    ser=F.fit([mp.mpf(e.numerator)/e.denominator**1*0+ (mp.mpf(e.numerator)/e.denominator)**p*n for e,n in zip(eps,Nf)])
    ser_d=F.fit([(mp.mpf(e.numerator)/e.denominator)**p*n for e,n in zip(eps,Nd)])
    lim=ser[0]; limd=ser_d[0]
    # check lower powers vanish: eps^(p-1)*N would diverge; check p-minimal: fit of eps^p N has finite value
    print(mu,rec['copy'],'k_B',kB,'p=2(2l2+kB)=',p,'formula vs direct N max rel dev (24 pts) %s'%mp.nstr(dev,3),'| lim eps^p N (formula) %s ; direct %s ; 1/X_B^2 %s ; rel diff %s'%(mp.nstr(lim,12),mp.nstr(limd,12),mp.nstr(1/XB**2,12),mp.nstr(lim*XB**2-1,3)),'%.0fs'%(time.time()-t0),flush=True)
