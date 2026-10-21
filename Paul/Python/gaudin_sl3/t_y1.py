from jet_basis_sym import *
import time
lam=(1,0,0)
for tag,key in (('P','2.7_-3.5'),('Pp','2.7_-3.5'),('X1','2.7_-3.5'),('P','2.7_-23.42')):
    th=th_from_key(key)
    pt=gk.get_point2(lam,tag,gk.ZETAS2[tag],th)
    t0=time.time(); Y,Psi=jet_Y(lam,pt)
    Yr,bad=recon(Y)
    print(tag,key,th,'bad',bad,time.time()-t0)
    if tag=='P' and key=='2.7_-3.5': Y0=Yr
    for i,row in enumerate(Yr): print(i,[str(x) for x in row])
