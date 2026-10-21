from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0); prs=gk.pairs_of(lam); n=len(prs)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
Y0i=coef[0].inv()
for k in range(1,len(coef)):
    M=coef[k]*Y0i
    print('=== M_%d  (<ij|_t = sum_n t^n sum M[(ij),(kl)] <kl|_0)'%k)
    for r in range(n*n):
        terms=[]
        for c in range(n*n):
            if M[r,c]!=0:
                terms.append('%s [%s%s|%s%s]'%(M[r,c],*prs[c//n],*prs[c%n]))
        if terms: print('[%s%s|%s%s] <- '%(*prs[r//n],*prs[r%n]),' + '.join(terms))
