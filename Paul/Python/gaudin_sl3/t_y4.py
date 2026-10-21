from jet_basis_sym import *
import time, sys, glob, re
L=int(sys.argv[1]); lam=(L,0,0); nmax=int(sys.argv[2])
fs=sorted(glob.glob('general_kernel_pt_%d00_P_*_80.pkl'%L))
out={}
cache='jet_basis_sym_Ypts%d.pkl'%L
if os.path.exists(cache): out=pickle.load(open(cache,'rb'))
for f in fs:
    key=re.match(r'general_kernel_pt_%d00_P_(.+)_80.pkl'%L,f).group(1)
    th=th_from_key(key); t=1/(th[0]-th[1])
    if t in out or len(out)>=nmax: continue
    pt=gk.get_point2(lam,'P',gk.ZETAS2['P'],th)
    t0=time.time(); Y,Psi=jet_Y(lam,pt); Yr,bad=recon(Y)
    print(key,'bad',bad,'%.0fs'%(time.time()-t0),flush=True)
    out[t]=Yr
    pickle.dump(out,open(cache,'wb'))
