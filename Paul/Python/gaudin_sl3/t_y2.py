from jet_basis_sym import *
import time, glob, re
lam=(2,0,0)
res={}
for f in sorted(glob.glob('general_kernel_pt_200_*_80.pkl')):
    m=re.match(r'general_kernel_pt_200_(\w+?)_(.+)_80.pkl',f); tag,key=m.group(1),m.group(2)
    th=th_from_key(key)
    pt=gk.get_point2(lam,tag,gk.ZETAS2[tag],th)
    t0=time.time(); Y,Psi=jet_Y(lam,pt); Yr,bad=recon(Y)
    res[(tag,key)]=(th,Yr)
    print(tag,key,'bad',bad,'%.1fs'%(time.time()-t0),'cond',mp.nstr(mp.mpf(1),2),flush=True)
pickle.dump(res,open('jet_basis_sym_Y200.pkl','wb'))
