import sys; sys.path.insert(0,'.')
from master_kernel import *
lam=(5,0,0)
w1,w1p=Fr(5,2),Fr(7,4); a,ap=1/w1,1/w1p
for th12 in (Fr(100),Fr(1000)):
    th=G.theta_pair(th12); t=1/th12
    mp.mp.dps=40
    pt=load_pt(lam,'P',th,40)
    t0=time.time()
    F,y=state_data(lam,pt,w1,w1p)
    print('state data',time.time()-t0,flush=True)
    for up in (0,1,2):
        Th=theta_from_symbol_orders(5,a,ap,t,up)
        r,_,_=residuals(lam,pt,w1,w1p,Th,F,y)
        print('theta12',th12,'orders<=',up,'max rel resid',mp.nstr(max(r),5),'median',mp.nstr(sorted(r)[len(r)//2],5),flush=True)
