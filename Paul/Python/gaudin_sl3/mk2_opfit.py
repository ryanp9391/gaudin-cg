import sys, itertools, pickle, time
import sympy as sp
import mk2_Nsym as S
E1,E2,H1,H2=S.E1,S.E2,S.H1,S.H2
V=[E1,E2,H1,H2]; WT=[1,2,1,2]
def weight_monos(w):
    out=[]
    for e in itertools.product(range(w+1),range(w//2+1),range(w+1),range(w//2+1)):
        if e[0]+2*e[1]+e[2]+2*e[3]==w: out.append(e)
    return out
def build_ansatz(n, order):
    """terms (coef-monomial exponents e, derivative orders d) with weight(e) - weight(d) = -n, |d_E|<=order, |d_H|<=order"""
    terms=[]
    for dE in itertools.product(range(order+1),repeat=2):
        if sum(dE)>order: continue
        for dH in itertools.product(range(order+1),repeat=2):
            if sum(dH)>order: continue
            d=(dE[0],dE[1],dH[0],dH[1])
            wd=d[0]+2*d[1]+d[2]+2*d[3]
            w=wd-n
            if w<0: continue
            for e in weight_monos(w):
                terms.append((e,d))
    return terms
def apply_term(term, f):
    e,d=term
    g=f
    for v,k in zip(V,d):
        if k: g=sp.diff(g,v,k)
    return sp.expand(g*sp.prod([v**k for v,k in zip(V,e)]))
def fit(L, n, order, N=None):
    N=N if N is not None else S.Nmat(L)
    out,prs=S.sym_op(L,N,n)
    terms=build_ansatz(n,order)
    cs=sp.symbols('c0:%d'%len(terms))
    eqs=[]
    for (r,s),res in out.items():
        f=S.iota(*prs[r],E1,E2)*S.iota(*prs[s],H1,H2)
        lhs=sum((cs[i]*apply_term(t,f) for i,t in enumerate(terms)),0)
        diff=sp.Poly(sp.expand(lhs-res),*V)
        eqs+= [c for c in diff.coeffs()]
    A,b=sp.linear_eq_to_matrix(eqs,cs)
    sol=sp.linsolve((A,b),cs)
    return terms,cs,sol,A.shape
if __name__=='__main__':
    L=int(sys.argv[1]); n=int(sys.argv[2]); order=int(sys.argv[3])
    t=time.time()
    terms,cs,sol,shape=fit(L,n,order)
    print('system',shape,'unknowns',len(terms),'time',time.time()-t)
    if sol==sp.EmptySet: print('NO SOLUTION')
    else:
        s=list(sol)[0]
        free=set().union(*[x.free_symbols for x in s])
        print('free params',len(free))
        for t_,v in zip(terms,s):
            if v!=0: print(t_,v)
