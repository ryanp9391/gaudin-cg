from jet_basis_sym import *
L=1; lam=(1,0,0)
coef,A=pickle.load(open('jet_basis_sym_A1.pkl','rb'))
alpha,s,hw,Eg=dictionary(1); d=3
E=ev_cov(1)
Pm=sp.zeros(9,9)
for p in range(3):
    for q in range(3): Pm[p*3+q,q*3+p]=1
lhs=(E*A).applyfunc(sp.simplify)
rhs=(a_s*ap_s/(a_s-ap_s)*E*(sp.eye(9)-Pm)).applyfunc(sp.simplify)
print('E A == aa\'/(a-a\') E (1-P):', (lhs-rhs).applyfunc(sp.simplify).is_zero_matrix)
Y0=coef[0]
th0=(E*Y0.inv()).applyfunc(sp.factor)   # theta (x) theta'
print('theta x theta\' =',th0)
z=t_s*a_s*ap_s/(ap_s-a_s)
# swapped
th0s=th0.subs({a_s:ap_s,ap_s:a_s},simultaneous=True)
Th=pickle.load(open('jet_basis_sym_Theta1.pkl','rb'))
# swapped jets: theta'(x)theta = theta'_i theta_j  -> entry (i,j) = th0[(j,i)] with a<->a' ... build explicitly
thv=[th0[0,i*3+0].subs(ap_s,1)  for i in range(3)]  # not robust; construct from c_kl
prs=gk.pairs_of(lam)
tha=[sp.Rational(gk.ckl_closed(lam,k,l).numerator,gk.ckl_closed(lam,k,l).denominator)*a_s**(k+l-1) for (k,l) in prs]
thb=[x.subs(a_s,ap_s) for x in tha]
Rm=sp.Matrix([[ (1+z)*tha[i]*thb[j]-z*thb[i]*tha[j] for i in range(3) for j in range(3)]])
print('Theta == (1+z) th(x)th\' - z th\'(x)th with z=t a a\'/(a\'-a):', (Th-Rm).applyfunc(sp.simplify).is_zero_matrix)
print('c_kl',[(kl,c) for kl,c in zip(prs,[x.subs(a_s,1) for x in tha])])
