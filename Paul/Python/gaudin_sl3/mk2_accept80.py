import sys; sys.path.insert(0,'.')
from master_kernel import *
lam=(5,0,0)
T=pickle.load(open('master_kernel2_rep_500.pkl','rb'))
prs=G.pairs_of(lam)
def theta_at(a,ap,t):
    return {(p,q): sum(c*t**n*a**pp*ap**qq for (n,pp,qq),c in T[(i,j)].items()) for i,p in enumerate(prs) for j,q in enumerate(prs)}
w1,w1p=Fr(5,2),Fr(7,4); a,ap=1/w1,1/w1p
th12=Fr(100); th=G.theta_pair(th12); t=1/th12
mp.mp.dps=80
pt=load_pt(lam,'P',th,80)
F,y=state_data(lam,pt,w1,w1p)
r,_,_=residuals(lam,pt,w1,w1p,theta_at(a,ap,t),F,y)
ys=[abs(v) for v in y]; ymax=max(ys)
rs=sorted(r)
print('80 digits P theta12=100: n',len(r),'max rel',mp.nstr(max(r),4),'median',mp.nstr(rs[len(rs)//2],4),'n>1e-40',sum(1 for v in r if v>1e-40),'n>1e-70',sum(1 for v in r if v>1e-70))
print('max abs resid/ymax',mp.nstr(max(r[k]*ys[k] for k in range(len(r)))/ymax,4))
# truncated comparison for reference
