from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0); prs=gk.pairs_of(lam); n=len(prs)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
Y0i=coef[0].inv()
cs=[sp.Rational(gk.ckl_closed(lam,k,l).numerator,gk.ckl_closed(lam,k,l).denominator) for (k,l) in prs]
print('c_kl',dict(zip(prs,cs)))
# tilde basis: jt = c * j.  <ij|_t = sum M[(ij),(kl)] <kl|_0, <ij|_0 = j_i j_j = jt_i jt_j/(c_i c_j)
# define Mt[(ij),(kl)] = M * c_i c_j /(c_k c_l)  for  <ij~|_t := c_i c_j <ij|_t = sum Mt jt_k jt_l
kk=int(sys.argv[2])
M=coef[kk]*Y0i
for r in range(n*n):
    terms=[]
    for c in range(n*n):
        if M[r,c]!=0:
            v=M[r,c]*cs[r//n]*cs[r%n]/(cs[c//n]*cs[c%n])
            terms.append('%s [%s%s|%s%s]'%(v,*prs[c//n],*prs[c%n]))
    if terms: print('[%s%s|%s%s] <- '%(*prs[r//n],*prs[r%n]),' + '.join(terms))
