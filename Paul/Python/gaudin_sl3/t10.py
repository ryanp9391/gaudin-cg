from untwist_limit import *
a = analyse_nu(PHI, ''); b = analyse_nu(PHI2, '2')
print('phi vs phi2 recovered slice max diff per branch:', [mp.nstr(max(abs(x-y) for x,y in zip(p['slice_q'], q['slice_q'])),3) for p,q in zip(a,b)])
print('X_B:', [(mp.nstr(p['XB'],10), mp.nstr(q['XB'],10)) for p,q in zip(a,b)])
