import sys; sys.path.insert(0,'.')
import nu_general as NG, sympy as sp
for lam in [(2,1),(3,1),(3,2),(4,1),(4,2)]:
    nl=NG.frobenius(lam)['nolog']
    print(lam)
    for nm,f in zip('A B1 B2'.split(),nl):
        p=sp.Poly(f,NG.Hs,NG.As,NG.Bs)
        tot=p.total_degree()
        top=[(m,sp.factor(c)) for m,c in p.terms() if sum(m)==tot]
        # coefficient of b alone
        cb=p.coeff_monomial(NG.Bs)
        print(' ',nm,'deg',tot,'top:',top[:6],' coef of b:',sp.factor(cb))
