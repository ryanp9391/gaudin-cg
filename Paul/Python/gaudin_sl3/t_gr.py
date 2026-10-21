from jet_basis_sym import *
for L in (1,2,3):
    lam=(L,0,0); prs=gk.pairs_of(lam); n=len(prs)
    coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
    alpha,s,hw,Eg=dictionary(L); d=len(alpha)
    wj=[k+l-1 for k,l in prs]
    rw=[wj[i]+wj[j] for i in range(n) for j in range(n)]
    cw=[alpha[p][1]+2*alpha[p][2]+alpha[q][1]+2*alpha[q][2] for p in range(d) for q in range(d)]
    bad=0;tot=0;diffs=set()
    for k,C in enumerate(coef):
        for r in range(n*n):
            for c in range(n*n):
                if C[r,c]!=0:
                    tot+=1; diffs.add((k,rw[r]-cw[c]))
    print(L,'(t power, rowlevel - colweight) pairs:',sorted(diffs))
