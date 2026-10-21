"""Fast exact fit of the grade-n part of the universal operator N(t) (symbol Phi_n) with dict polynomials + sparse QQ solve.
   Phi_n = known_rest(n) + D_n with D_n general (p-degrees dmin..dmax)."""
import sys, itertools, pickle, time
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
import mk2_Nsym as S
import general_kernel as G
import master_kernel2 as m
import mk2_opfit as F
import mk2_expfit as X

def padd(p,q,c=1):
    for k,v in q.items():
        p[k]=p.get(k,0)+c*v
        if p[k]==0: del p[k]
def pmul(p,q):
    out={}
    for k1,v1 in p.items():
        for k2,v2 in q.items():
            k=tuple(a+b for a,b in zip(k1,k2)); out[k]=out.get(k,0)+v1*v2
    return {k:v for k,v in out.items() if v!=0}
def todict(expr, gens):
    P=sp.Poly(sp.expand(expr),*gens)
    return {mon:Fr(int(co.p),int(co.q)) for mon,co in P.terms()}
def iota_dict(k,l,site):
    e1,e2=(S.E1,S.E2) if site==1 else (S.H1,S.H2)
    return todict(S.iota(k,l,e1,e2), (e1,e2))   # keys (i,j)
def inputs(L):
    prs=G.pairs_of((L,0,0)); npr=len(prs)
    io1=[iota_dict(k,l,1) for k,l in prs]; io2=[iota_dict(k,l,2) for k,l in prs]
    out={}
    for r in range(npr):
        for s in range(npr):
            out[(r,s)]={(a[0],a[1],b[0],b[1]):va*vb for a,va in io1[r].items() for b,vb in io2[s].items()}
    return out,prs
def apply_term(term, f):
    """term=(e,d) coefficient monomial e (4-tuple), derivative orders d (4-tuple) applied to dict poly f"""
    e,d=term; out={}
    for mon,c in f.items():
        cc=c; ok=True
        new=[]
        for x,dd,ee in zip(mon,d,e):
            if dd>x: ok=False;break
            cc*= (factorial(x)//factorial(x-dd))
            new.append(x-dd+ee)
        if ok: out[tuple(new)]=out.get(tuple(new),0)+cc
    return out
def apply_symbol_dict(symb_terms, f):
    """symb_terms dict {(e,d):coeff}"""
    out={}
    for (e,d),c in symb_terms.items():
        r=apply_term((e,d),f)
        for k,v in r.items():
            out[k]=out.get(k,0)+c*v
    return {k:v for k,v in out.items() if v!=0}
def sym_to_terms(symb):
    P=sp.Poly(sp.expand(symb),*(X.V+X.PV))
    return {(mon[:4],mon[4:]):Fr(int(co.p),int(co.q)) for mon,co in P.terms()}
def data_output(L,N,n):
    """exact grade-n output polys for each input (r,s), dict polys"""
    prs=G.pairs_of((L,0,0)); npr=len(prs)
    lev=[p[0]+p[1]-1 for p in prs]
    io1=[iota_dict(k,l,1) for k,l in prs]; io2=[iota_dict(k,l,2) for k,l in prs]
    outp={}
    # N as dict rows
    for r in range(npr):
        for s in range(npr):
            row=r*npr+s
            acc={}
            for c in range(npr*npr):
                v=N[row,c]
                if v==0: continue
                p,q=divmod(c,npr)
                if lev[r]+lev[s]-lev[p]-lev[q]!=n: continue
                vv=Fr(int(v.p),int(v.q))
                for a,va in io1[p].items():
                    for b,vb in io2[q].items():
                        k=(a[0],a[1],b[0],b[1]); acc[k]=acc.get(k,0)+vv*va*vb
            outp[(r,s)]={k:v for k,v in acc.items() if v!=0}
    return outp
def terms_general(n,dmin,dmax,minside=1):
    terms=[]
    for tot in range(dmin,dmax+1):
        for kE in range(minside,tot-minside+1):
            kH=tot-kE
            for dE in [(i,kE-i) for i in range(kE+1)]:
                for dH in [(i,kH-i) for i in range(kH+1)]:
                    d=(dE[0],dE[1],dH[0],dH[1])
                    w=d[0]+2*d[1]+d[2]+2*d[3]-n
                    if w<0: continue
                    for e in F.weight_monos(w): terms.append((e,d))
    return terms
def fit_general(n,L,rest_terms,terms,N=None):
    N=N if N is not None else S.Nmat(L)
    outp=data_output(L,N,n)
    inp,prs=inputs(L)
    # right-hand side: data - rest applied
    rows={}   # (input, outmono) -> row index
    cols=len(terms)
    entries=[]  # (rowidx, col, val)
    rhs={}
    def rid(key):
        if key not in rows: rows[key]=len(rows)
        return rows[key]
    for key,f in inp.items():
        b=dict(outp[key]); 
        rr=apply_symbol_dict(rest_terms,f) if rest_terms else {}
        padd(b,rr,-1)
        for k,v in b.items(): rhs[rid((key,k))]=v
        for ci,t in enumerate(terms):
            r=apply_term(t,f)
            for k,v in r.items():
                entries.append((rid((key,k)),ci,v))
    nr=len(rows)
    print('system',nr,'x',cols,'nnz',len(entries),flush=True)
    # augmented sparse matrix
    rowsd={}
    for i,j,v in entries:
        rowsd.setdefault(i,{})
        rowsd[i][j]=rowsd[i].get(j,0)+v
    for i,v in rhs.items():
        rowsd.setdefault(i,{})
        rowsd[i][cols]=rowsd[i].get(cols,0)+v
    Rr={i:{j:QQ(v.numerator,v.denominator) for j,v in d.items() if v!=0} for i,d in rowsd.items()}
    for i in range(nr):
        Rr.setdefault(i,{})
    dm=DomainMatrix(Rr,(nr,cols+1),QQ)
    t=time.time()
    rref,piv=dm.rref()
    print('rref time',time.time()-t,'rank',len(piv),flush=True)
    if cols in piv: return None,None
    sdm=rref.rep.to_sdm()
    sol={}
    freecols=[j for j in range(cols) if j not in piv]
    for ri,rowd in sdm.items():
        pj=min(rowd) 
        if pj>=cols: continue
        rhsv=rowd.get(cols,0)
        dependents=[j for j in rowd if j!=pj and j<cols]
        sol[pj]=(rhsv,dependents)
    return sol,freecols,terms
