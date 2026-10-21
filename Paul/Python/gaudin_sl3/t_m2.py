from jet_basis_sym import *
import sys
Ms={}
for L in (1,2,3):
    lam=(L,0,0); prs=gk.pairs_of(lam); n=len(prs)
    coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
    print('L',L,'t-degree',len(coef)-1)
    Y0i=coef[0].inv()
    Ms[L]=[ (coef[k]*Y0i, prs) for k in range(1,len(coef))]
def entries(L,k):
    M,prs=Ms[L][k-1]; n=len(prs); out={}
    for r in range(n*n):
        for c in range(n*n):
            if M[r,c]!=0: out[(prs[r//n],prs[r%n],prs[c//n],prs[c%n])]=M[r,c]
    return out
for k in (1,2):
    e1=entries(1,k) if k<=len(Ms[1]) else {}
    e2=entries(2,k); e3=entries(3,k)
    # restrict e3 to labels with all indices<=3 (L=2 labels)
    lab2=set(gk.pairs_of((2,0,0)))
    r3={key:v for key,v in e3.items() if all(x in lab2 for x in key)}
    print('k',k,'L=2 vs L=3 restricted equal:',r3==e2, len(r3),len(e2))
    lab1=set(gk.pairs_of((1,0,0)))
    r2={key:v for key,v in e2.items() if all(x in lab1 for x in key)}
    print('   L=1 vs L=2 restricted equal:',r2==e1,len(r2),len(e1))
if len(Ms[3])>2:
    print('M3 nnz L=3',len(entries(3,3)))
