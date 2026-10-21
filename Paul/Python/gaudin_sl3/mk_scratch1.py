import sys, itertools, time
sys.path.insert(0,'.')
from master_kernel import *
from sympy import Rational, Matrix, QQ
from sympy.polys.matrices import DomainMatrix
tab = get_tab()
def monos(d, nv=4):
    out=[]
    for e in itertools.product(range(d+1), repeat=nv):
        if sum(e)<=d: out.append(e)
    return out
def fitdeg(l1, m, mp_, dmax=8):
    R = tab[l1]; n=l1+1
    prs=[(k,l) for k in range(n+1) for l in range(k+1,n+1)]
    rows=[]; rhs=[]
    pts=[]
    for p in prs:
        for q in prs:
            pts.append((p,q,R.get((p[0],p[1],q[0],q[1],m,mp_),Fr(0))))
    for d in range(0,dmax+1):
        ms=monos(d)
        if len(ms)>len(pts): return None
        A=[[Rational(p[0]**e[0]*p[1]**e[1]*q[0]**e[2]*q[1]**e[3]) for e in ms]+[Rational(v.numerator,v.denominator)] for p,q,v in pts]
        M=DomainMatrix.from_Matrix(Matrix(A)).convert_to(QQ)
        r1=M.rank(); 
        M0=M[:, :len(ms)]
        r0=M0.rank()
        if r1==r0: return d, len(ms)
    return None
if __name__=='__main__':
    l1=int(sys.argv[1])
    for (m,mp_) in sorted({(k[4],k[5]) for k in tab[l1]}, key=lambda x:(x[0]+x[1],x)):
        t=time.time(); print(l1,(m,mp_),fitdeg(l1,m,mp_),round(time.time()-t,1),flush=True)
