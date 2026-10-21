from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
Y0=coef[0]; 
A2=Y0.inv()*coef[2] if len(coef)>2 else None
print('A^2 nnz',sum(1 for x in A*A if x!=0),'A^3 zero?',(A*A*A).is_zero_matrix, 'A^2 zero?',(A*A).is_zero_matrix)
if A2 is not None:
    print('A2 nnz',sum(1 for x in A2 if x!=0))
    ratios=set()
    for i in range(A2.shape[0]):
        for j in range(A2.shape[1]):
            if (A*A)[i,j]!=0: ratios.add(A2[i,j]/(A*A)[i,j])
            elif A2[i,j]!=0: ratios.add('extra')
    print('A2/A^2 ratios',ratios)
# weight-grading
alpha,s,hw,Eg=dictionary(L); d=len(alpha)
# principal weight w=alpha2+2alpha3 of two-site basis element
wgt=[alpha[p][1]+2*alpha[p][2]+alpha[q][1]+2*alpha[q][2] for p in range(d) for q in range(d)]
degs=set(wgt[r]-wgt[c] for r in range(d*d) for c in range(d*d) if A[r,c]!=0)
print('A changes principal weight by',degs)
if A2 is not None: print('A2 changes by',set(wgt[r]-wgt[c] for r in range(d*d) for c in range(d*d) if A2[r,c]!=0))
# gl3 weights
w3=[tuple(alpha[p][k]+alpha[q][k] for k in range(3)) for p in range(d) for q in range(d)]
print('A gl3 weight shifts',set(tuple(w3[r][k]-w3[c][k] for k in range(3)) for r in range(d*d) for c in range(d*d) if A[r,c]!=0))
