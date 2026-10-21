import sys; sys.path.insert(0,'.')
from master_kernel import *
lam=(4,0,0)
th=G.TH0
pt=load_pt(lam,'P',th,80)
w1,w1p=Fr(5,2),Fr(7,4); a,ap=1/w1,1/w1p; t=1/(th[0]-th[1])
t0=time.time()
Th=theta_from_pkl(4,a,ap,t,None)
r,F,y=residuals(lam,pt,w1,w1p,Th)
print('full pkl: max rel residual',mp.nstr(max(r),5),time.time()-t0)
for up in (1,2):
    Th=theta_from_pkl(4,a,ap,t,up)
    r,_,_=residuals(lam,pt,w1,w1p,Th,F,y); print('trunc',up,mp.nstr(max(r),5),mp.nstr(min(r),5))
