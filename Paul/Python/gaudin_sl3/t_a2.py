from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
alpha,s,hw,Eg=dictionary(L); d=len(alpha)
S=sp.diag(*[s[p]*s[q] for p in range(d) for q in range(d)])
Am=S*A*S.inv()
names=[(alpha[p],alpha[q]) for p in range(d) for q in range(d)]
def nm(m): return 'x'+''.join(str(i+1)*m[0][i] for i in range(3)) if False else ''
for r in range(d*d):
    for c in range(d*d):
        if Am[r,c]!=0: print('c[%s;%s] <- %s * c[%s;%s]'%(names[r][0],names[r][1],Am[r,c],names[c][0],names[c][1]))
