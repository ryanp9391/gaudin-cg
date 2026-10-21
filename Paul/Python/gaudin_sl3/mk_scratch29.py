import pickle, sympy as sp, sys, itertools
z,w,wp=sp.symbols('z w wp')
l=int(sys.argv[1])
out=pickle.load(open('mk_Qw_%d.pkl'%l,'rb'))
def parts(l): return [(l-b-c,b,c) for b in range(l+1) for c in range(l+1-b)]
rows=parts(l)   # (a0,a1,a2): U_inf^a0 U^a1 U'^a2  ; U = P(w), U' = P(wp)
cols=parts(l)   # (c0,c1,c2): V_inf^c0 V^c1 V'^c2 ; V=Q(wp), V'=Q(w)
M=sp.zeros(len(rows),len(cols))
for (a,b,c,d),co in out.items():
    r=rows.index((l-a-b,a,b)); cidx=cols.index((l-c-d,c,d))
    M[r,cidx]=sp.factor(co.subs(w,wp+1))
sp.pprint(M.applyfunc(sp.factor),wrap_line=False) if l==2 else None
print('rows/cols order',rows)
print('det',sp.factor(M.det()))
pickle.dump((rows,M),open('mk_Cmat_%d.pkl'%l,'wb'))
