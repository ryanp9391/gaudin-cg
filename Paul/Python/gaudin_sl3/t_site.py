from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
alpha,s,hw,Eg=dictionary(L); d=len(alpha); prs=gk.pairs_of(lam); n=len(prs)
Y0=coef[0]
S=sp.diag(*[s[p]*s[q] for p in range(d) for q in range(d)])
Ym=Y0*S.inv()   # acts on monomial coefficients c[alpha;beta]
# one-site: row (i,j) = y_i (x) y_j ; take y_i = row (i,i0) restricted ... use entries
# find a one-site matrix y (n x d) with Ym[(i,j),(p,q)] = y[i,p] y[j,q]
import itertools
# pick (j0,q0) with nonzero y[j0,q0]
y=sp.zeros(n,d)
# determine column-wise: Ym[(i,j),(p,q)]
for j0 in range(n):
    for q0 in range(d):
        blk=sp.Matrix(n,d,lambda i,p: Ym[i*n+j0,p*d+q0])
        if any(x!=0 for x in blk): break
    else: continue
    break
c=Ym[0*n+j0,0]  # placeholder
# normalise: use first nonzero element
import math
# y[:, :] = blk / y[j0,q0]; fix scale by requiring y[0,hw]=1
y=blk
sc=y[0,hw]
y=y/sc
# check factorisation
ok=all(sp.simplify(Ym[i*n+j,p*d+q]-y[i,p]*y[j,q])==0 for i in range(n) for j in range(n) for p in range(d) for q in range(d))
print('factorises',ok)
for i,kl in enumerate(prs):
    terms=['%s*c%s'%(y[i,p],''.join(str(k+1)*alpha[p][k] for k in range(3))) for p in range(d) if y[i,p]!=0]
    print(kl,' + '.join(terms))
