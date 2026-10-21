from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
alpha,s,hw,Eg=dictionary(L); d=len(alpha)
w3=[tuple(alpha[p][k]+alpha[q][k] for k in range(3)) for p in range(d) for q in range(d)]
def sect(X,sh):
    Z=sp.zeros(*X.shape)
    for r in range(X.shape[0]):
        for c in range(X.shape[1]):
            if X[r,c]!=0 and tuple(w3[r][k]-w3[c][k] for k in range(3))==sh: Z[r,c]=X[r,c]
    return Z
A21=sect(A,(-1,1,0)); A32=sect(A,(0,-1,1))
A2=coef[0].inv()*coef[2]
for sh,pr in (((-2,2,0),A21*A21),((0,-2,2),A32*A32),((-1,0,1),A21*A32+A32*A21)):
    S=sect(A2,sh)
    rat=set()
    for r in range(S.shape[0]):
        for c in range(S.shape[1]):
            if pr[r,c]!=0: rat.add(S[r,c]/pr[r,c])
            elif S[r,c]!=0: rat.add('extra')
    print(sh,'ratios A2/prod',rat)
    for nm,pr2 in (('A21A32',A21*A32),('A32A21',A32*A21)):
        if sh==(-1,0,1):
            rat=set()
            for r in range(S.shape[0]):
                for c in range(S.shape[1]):
                    if pr2[r,c]!=0: rat.add(S[r,c]/pr2[r,c])
                    elif S[r,c]!=0: rat.add('extra')
            print('   ',nm,rat)
