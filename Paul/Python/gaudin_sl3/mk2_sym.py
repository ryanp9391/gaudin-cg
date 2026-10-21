import sys, time, pickle
import sympy as sp
import master_kernel2 as m
gens = (m.E1, m.E2, m.H1, m.H2, m.a_, m.ap_)
res = {}
for L in range(1, int(sys.argv[1]) + 1):
    t0 = time.time()
    Q = m.symbol_orders_exact(L)
    for n, q in sorted(Q.items()):
        q = sp.Poly(sp.expand(q), *gens)
        if n <= L:
            Pn = sp.Poly(sp.expand(m.Pi_ ** (L - n)), *gens)
            qq, r = sp.div(q, Pn)
            res[(L, n)] = (qq.as_expr(), r.is_zero)
            print(L, n, 'divisible', r.is_zero, 'terms', len(qq.terms()), round(time.time() - t0, 1), flush=True)
        else:
            res[(L, n)] = (q.as_expr(), None)
            print(L, n, 'raw (n>L)', len(q.terms()), flush=True)
pickle.dump(res, open('master_kernel2_Pn.pkl', 'wb'))
