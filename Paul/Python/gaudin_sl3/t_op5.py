from jet_basis_sym import *
import sys
L=int(sys.argv[1])
res=pickle.load(open('jet_basis_sym_Yops%d.pkl'%L,'rb'))
coef=poly_fit_in_t(res,L)
print('t-degree',len(coef)-1)
# the operator-derived Y is in the GT companion basis with the same ordering as jets_Y -> compare to state-derived coef
try:
    c0,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
    print('equal to state-derived Y(t):',all(c0[k]==coef[k] for k in range(len(c0))) and len(c0)==len(coef))
except Exception as e: print('no state-derived file',e)
vals=[{a_s:sp.Rational(3,7),ap_s:sp.Rational(-5,4),t_s:sp.Rational(2,9)},{a_s:sp.Rational(-2,3),ap_s:sp.Rational(7,5),t_s:sp.Rational(-3,11)}]
print('Theta(operator-derived Y) vs cached kernel mismatches:',theta_check(L,coef,vals))
pickle.dump(coef,open('jet_basis_sym_Ycoef_ops%d.pkl'%L,'wb'))
