from jet_basis_sym import *
import sys
L=int(sys.argv[1])
nu,cond=nu_formulas(L)
print('no-log condition:',sp.factor(cond))
for kl,v in nu.items(): print(kl, sp.factor(v))
