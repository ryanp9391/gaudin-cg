"""Test: normal symbol of the universal operator N(t) is exp(sum_n t^n D_n), D_n of total p-degree n+1.  Fit D_3 (and D_4) assuming the hypothesis."""
import sys, itertools, pickle, time
import sympy as sp
import mk2_Nsym as S
import mk2_opfit as F
E1,E2,H1,H2=S.E1,S.E2,S.H1,S.H2
pe1,pe2,ph1,ph2=sp.symbols('pe1 pe2 ph1 ph2')
V=[E1,E2,H1,H2]; PV=[pe1,pe2,ph1,ph2]
def sym_from_terms(d):
    return sum(sp.Rational(v)*sp.prod([x**k for x,k in zip(V,e)])*sp.prod([p**k for p,k in zip(PV,dd)]) for (e,dd),v in d.items())
def apply_symbol(symb, f):
    """normal-ordered action: coefficients (in V) multiply after derivatives (p -> d/dV)"""
    P=sp.Poly(sp.expand(symb),*(PV))
    out=0
    for dd,co in P.terms():
        g=f
        for v,k in zip(V,dd):
            if k: g=sp.diff(g,v,k)
        out+=co*g
    return sp.expand(out)
D1=(E1-H1)*pe1*ph1+(E2-H2)*(pe1*ph2+pe2*ph1)+(E2*H1-E1*H2)*pe2*ph2
def getD2(L=4):
    N2,diff=pickle.load(open('master_kernel2_N2op.pkl','rb'))
    return sym_from_terms(diff)
def fit_D(n, known_rest, L, order_terms):
    """N_n data = [exp-series coefficient] ; known_rest = symbol of everything except D_n at order n. unknown D_n over order_terms"""
    N=S.Nmat(L)
    out,prs=S.sym_op(L,N,n)
    cs=sp.symbols('c0:%d'%len(order_terms))
    Dn=sum(cs[i]*sp.prod([x**k for x,k in zip(V,e)])*sp.prod([p**k for p,k in zip(PV,d)]) for i,(e,d) in enumerate(order_terms))
    eqs=[]
    for (r,s),res in out.items():
        f=S.iota(*prs[r],E1,E2)*S.iota(*prs[s],H1,H2)
        lhs=apply_symbol(known_rest+Dn,f)
        diff=sp.Poly(sp.expand(lhs-res),*V)
        eqs+=list(diff.coeffs())
    A,b=sp.linear_eq_to_matrix(eqs,cs)
    sol=sp.linsolve((A,b),cs)
    return cs,sol,A.shape
def terms_total(n):
    """(e,d) with |d_E|>=1,|d_H|>=1, |d|=n+1, coefficient weight = wt(d)-n"""
    terms=[]
    tot=n+1
    for kE in range(1,tot):
        kH=tot-kE
        for dE in [(i,kE-i) for i in range(kE+1)]:
            for dH in [(i,kH-i) for i in range(kH+1)]:
                d=(dE[0],dE[1],dH[0],dH[1])
                wd=d[0]+2*d[1]+d[2]+2*d[3]
                w=wd-n
                if w<0: continue
                for e in F.weight_monos(w): terms.append((e,d))
    return terms
if __name__=='__main__':
    n=int(sys.argv[1]); L=int(sys.argv[2])
    D2=getD2()
    # check D2 assumption: N2 = D1^2/2 + D2 (by construction)
    if n==3:
        rest=D1**3/6+D1*D2
    if n==4:
        D3=sym_from_terms(pickle.load(open('master_kernel2_D3.pkl','rb')))
        rest=D1**4/24+D1**2*D2/2+D2**2/2+D1*D3
    t=time.time()
    terms=terms_total(n)
    print('unknowns',len(terms),flush=True)
    cs,sol,shape=fit_D(n,rest,L,terms)
    print('system',shape,'time',time.time()-t)
    if sol==sp.EmptySet: print('NO SOLUTION (hypothesis exp(sum t^n D_n) with D_n of degree n+1 fails)')
    else:
        s=list(sol)[0]
        free=set().union(*[x.free_symbols for x in s])
        print('free',len(free))
        D3={t_:v for t_,v in zip(terms,s) if v!=0}
        pickle.dump(D3,open('master_kernel2_D%d.pkl'%n,'wb'))
        print('D%d nonzero terms'%n,len(D3))
