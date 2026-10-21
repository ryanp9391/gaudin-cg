import sys; sys.path.insert(0,'.')
from master_kernel import *
import master_kernel2 as m2
lam=(5,0,0)
T=pickle.load(open('master_kernel2_rep_500.pkl','rb'))
prs=G.pairs_of(lam)
def theta_at(a,ap,t,maxn=None):
    return {(p,q): sum(c*t**n*a**pp*ap**qq for (n,pp,qq),c in T[(i,j)].items() if maxn is None or n<=maxn) for i,p in enumerate(prs) for j,q in enumerate(prs)}
w1,w1p=Fr(5,2),Fr(7,4); a,ap=1/w1,1/w1p
for tag,th12 in (('P',Fr(100)),('P',Fr(1000)),('P',Fr(31,35)),('Pp',Fr(31,35))):
    th=G.theta_pair(th12); t=1/th12
    mp.mp.dps=40
    try: pt=load_pt(lam,tag,th,40)
    except Exception as e: print('missing',tag,th12,e,flush=True); continue
    F,y=state_data(lam,pt,w1,w1p)
    r,_,_=residuals(lam,pt,w1,w1p,theta_at(a,ap,t),F,y)
    print(tag,'theta12',th12,'nstates',len(r),'max rel resid',mp.nstr(max(r),5),'median',mp.nstr(sorted(r)[len(r)//2],5),flush=True)
print('--- diagnostics')
for tag,th12 in (('P',Fr(100)),('P',Fr(1000)),('P',Fr(31,35)),('Pp',Fr(31,35))):
    th=G.theta_pair(th12); t=1/th12
    pt=load_pt(lam,tag,th,40)
    F,y=state_data(lam,pt,w1,w1p)
    r,_,_=residuals(lam,pt,w1,w1p,theta_at(a,ap,t),F,y)
    ys=[abs(v) for v in y]; ymax=max(ys)
    bad=[(k) for k,v in enumerate(r) if v>1e-18]
    print(tag,th12,'n>1e-18:',len(bad),'n>1e-24',sum(1 for v in r if v>1e-24),'ok flags bad',sum(1 for s in pt['states'] if not s['ok']),
      'abs resid/ymax max',mp.nstr(max(r[k]*ys[k] for k in range(len(r)))/ymax,4),
      'min|y|/ymax',mp.nstr(min(ys)/ymax,4),'|y| of worst',[mp.nstr(ys[k]/ymax,3) for k in bad[:5]])
