import sys; sys.path.insert(0,'.')
from master_kernel import *
a,ap,t=Fr(2,5),Fr(3,7),Fr(1,3)
for l1 in (2,3,4):
    X=theta_from_symbol_orders(l1,a,ap,t,2); Y=theta_from_pkl(l1,a,ap,t,2)
    bad=sum(1 for k in X if X[k]!=Y[k]); print(l1,'mismatch orders<=2:',bad,'of',len(X))
    X=theta_from_symbol_orders(l1,a,ap,t,1); Y=theta_from_pkl(l1,a,ap,t,1)
    print(l1,'order<=1 mismatch',sum(1 for k in X if X[k]!=Y[k]))
