"""Frobenius at theta_s of the dual oper P''' - S1 P'' + S2 P' - S3 P = 0: nu = (1, kappa, pi, kappa pi, eps) from Laurent data only (symbolic)."""
import sympy as sp
x = sp.symbols('x')
s1_, e2, e3, Hs, Ht, As, At, Bs, Bt, dl = sp.symbols('sigma1 e2 e3 H_s H_t a_s a_t b_s b_t delta')


def taylor_s(NT=6):
    """Taylor coefficients (in x = u - theta_s) of s_i = x^2 S_i,  delta = theta_s - theta_t  (x_t = x + delta)"""
    S1 = s1_ + 3 / x + 3 / (x + dl)
    S2 = e2 + Hs / x + 3 / x ** 2 + Ht / (x + dl) + 3 / (x + dl) ** 2
    S3 = e3 + Bs / x + As / x ** 2 + Bt / (x + dl) + At / (x + dl) ** 2
    out = []
    for S in (S1, S2, S3):
        ex = sp.series(sp.together(x ** 2 * S), x, 0, NT).removeO()
        out.append([sp.simplify(ex.coeff(x, k)) for k in range(NT)])
    return out


def frob(s1, s2, s3, f, start, free, N):
    """solve recursion for series f (dict m->value, preset), coefficient of x^N eq gives f_{N+1}"""
    pass


def solve_nu(NT=6):
    s1, s2, s3 = taylor_s(NT)
    def coef(Nn, fvals):
        # equation at order x^Nn:  f_{N+1}(N+1)N(N-1) - sum_k s1_k f_{N+2-k}(N+2-k)(N+1-k) + sum_k s2_k f_{N+1-k}(N+1-k) - sum_k s3_k f_{N-k} = 0
        tot = 0
        F = lambda m: fvals.get(m, 0)
        tot += F(Nn + 1) * (Nn + 1) * Nn * (Nn - 1)
        for k in range(NT):
            m = Nn + 2 - k
            tot -= s1[k] * F(m) * m * (m - 1)
            m = Nn + 1 - k
            tot += s2[k] * F(m) * m
            m = Nn - k
            tot -= s3[k] * F(m)
        return tot
    # f0-type: f0 = 1, f2 = 0 (free), unknown f1, f3
    k1, c3, c2 = sp.symbols('kappa c3 c2')
    fv = {0: 1, 1: k1, 2: 0, 3: c3}
    e0 = sp.simplify(coef(0, fv))      # fixes kappa
    sol_k = sp.solve(e0, k1)[0]
    fv[1] = sol_k
    e2_ = sp.simplify(coef(1, fv))     # no-log condition at m=2 (f2 coefficient drops out)
    e3_ = sp.simplify(coef(2, fv))     # fixes c3
    sol_c3 = sp.solve(e3_, c3)[0]
    fv[3] = sol_c3
    # f2-type: f2 = 1, f3 = d3
    d3 = sp.symbols('d3')
    gv = {2: 1, 3: d3}
    e_d = coef(2, gv)
    sol_d3 = sp.solve(e_d, d3)[0]
    kap = sp.simplify(sol_k)
    pi_ = sp.simplify(3 * sol_d3)
    c3v = sp.simplify(sol_c3)
    eps_ = sp.simplify(-6 * c3v + 6 * 0 * sol_d3)        # c2 = 0 choice
    return dict(kappa=kap, pi=pi_, eps=eps_, nolog1=sp.simplify(e2_), c3=c3v, d3=sp.simplify(sol_d3))


if __name__ == '__main__':
    r = solve_nu()
    for k, v in r.items():
        print(k, '=', sp.factor(v))
