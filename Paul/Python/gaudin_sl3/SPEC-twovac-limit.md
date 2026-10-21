# Spec D2: the Gaudin two-vacua formula as the formal h → 0 limit of the XXX ten-product bracket formula (2026-10-16)

Delegated to a general-purpose subagent (sonnet), Python only, exact arithmetic (sympy / Fraction) until the final numerical test. Read
first: `Paul/Mathematica/Gaudin/Gaudin-sl3-round2-plan.md` §0, §1 (binding: the ladder → jet dictionary and the warnings), §5, §8; the
round-1 `RESULTS-xxx-limit.md` (P2: the XXX port at z = e^{hζ}, series class, x₀ limit, Baxter → oper) and `RESULTS-Qside.md` §0, §2
(Gaudin Q conventions, Wronskians). **The XXX formula:** `Paul/Mathematica/Gaudin/ref/xxx-210-infinf-formula-cells.txt` (extract of the
XXX notebook `L2-210-two-vacua-brackets.wb`: setup cells with the ladder/bracket/Bk2 definitions and the K tensor, the ten products, the
closed forms in (s, t) [cell 62: kernel laws κ(s,t) and the wedge vectors f, g, f′, g′ per product with the gauge and (e, e′) per block],
the Q-side dictionary [cells 49–57], the summary [cell 70]); the full dump `ref-L2-210-two-vacua-brackets-all.txt` and
`ref-L2-two-vacua-untwist-all.txt` (frame tie `kappaTie`, `Pphys`, `x0Of`, `ΩD`, `nOf`) are in the same directory. The prose version:
`Paul/Mathematica/new-conventions/Two-site-rotation-brackets-notes.md` §1 and §5.

**Code to reuse (read-only, import; never edit):** `xxx_chain.py` (`Ser` formal series in h, the XXX chain, `gauss_r`, `covec_Omega`),
`xxx_limit.py` (`Paux_exact`, `Gc_exact`, `gprime`, `mp_x0_onesite`, `truncated_exp`), `gaudin_chain.py` and `qside.py` (Gaudin states,
Q solvers, vacuum) for the numerical test in L3. **Files you own:** `twovac_limit.py` (runnable top-to-bottom, `python -X utf8
twovac_limit.py`), `RESULTS-twovac-limit.md`. Environment: Python 3.13, sympy 1.14, mpmath 1.3; no `uv`. Parse Mathematica InputForm with
`sympy.parsing.mathematica.parse_mathematica` (symbols `sS` = s, `tt` = t; `aS` = a, `bS` = a′ in the raw cells) and VERIFY every parse
(L0) before using it.

## 0. The XXX identity and the scaling (what is being limited)

XXX, companion frame: Φ_c(a,a′) = gaussCov{a,0,a}⊗gaussCov{a′,0,a′} = Ω_c[(1/a,∞)]⊗Ω_c[(1/a′,∞)], and for every eigenstate Ψ_c (companion
frame, Ψ_c = (P⁻¹⊗P⁻¹)ψ_d, P = P_phys(z)):

    ⟨Φ_c|Ψ_c⟩/⟨x₀^c|Ψ_c⟩ · M¹_vac M²_vac = Σ_{c,c′∈{0,1}} Q³(θ₁−ch) Q³(θ₂−c′h) Σ_{k=1}^{10} Bk2^{(1)}[F_k,G_k] Bk2^{(2)}[F′_k,G′_k],

Bk2^{(s)}[F,G] = ⟨Q³F⟩_s⟨Q²G⟩_s − ⟨Q²F⟩_s⟨Q³G⟩_s, ⟨Q^jF⟩_s = Σ_{m=−3}^{0} F_m Q^j(θ_s+mh), F_m = a^m f_m(s,t), G_m = a^{m+e} g_m(s,t) on
ladder 1, F′_n = a′^n f′_n, G′_n = a′^{n+e′} g′_n on ladder 2, (e,e′) = (7,5),(7,4),(6,5),(6,4) for (c,c′) = (0,0),(0,1),(1,0),(1,1), s = a′/a,
t = θ₁₂/h, M^{(s)}_vac = minor_s(−3,−2)Q³(θ_s−h), minor_s(m,m′) = Q²(θ_s+m′h)Q³(θ_s+mh) − Q³(θ_s+m′h)Q²(θ_s+mh) (check the sign/order
convention in the setup cells), Q^j the dual Q's of the state (Q^j = z_j^{u/h}q^j, monic q^j of degree μ_j), x₀^c = e₁⊗e₁ in the
companion frame. Dictionary to the diagonal frame (cells 49–57): ⟨Φ_c|Ψ_c⟩/⟨x₀^c|Ψ_c⟩ = ⟨Φ_d|ψ_d⟩/(κκ′) with Φ_d = Ω_D[z,(w₁,∞)]⊗Ω_D[z,(w₁′,∞)],
κ = kappaTie(z,w), κ′ = kappaTie(z,w′), and ⟨x₀^c|Ψ_c⟩ = ⟨x₀^d|ψ_d⟩ with x₀^d = x0Of(z). Gaudin scaling: z = e^{hζ}, w₁ = e^{hω₁} (so
a = 1/w₁ = e^{−hω₁}, a′ = e^{−hω₁′}, s = e^{−h(ω₁′−ω₁)}), t = θ₁₂/h; Q^j(u;h) = e^{ζ_j u}[q^j_0(u) + h q^j_1(u) + h² q^j_2(u) + …] with
q^j_0 the Gaudin dual polynomial (P^j = e^{ζ_j u}q^j_0 solves the Gaudin oper). **Warnings (plan §1):** the companion frame degenerates at
h = 0 (P_aux(z) → rank 1), so work with the diagonal-frame form of the identity; gaussCov{a,0,a} → e₁ᵀρ(e^{E₁₂}e^{E₂₃}) as h → 0, so the
ω-dependence of the family is O(h) and the leading nonvanishing order of the right-hand side can be high with cancellations between the
products; the Wronskian of the Gaudin dual Q's vanishes at θ_s (so minor_s(−3,−2) = (5/2)h²W′(θ_s) + O(h³), not O(h)).

## 1. Tasks

**L0 — parse and verify the closed forms.** From cell 62: the 40 kernel laws κ^{(cc′)}_{(mm′),(nn′)}(s,t) and the ten products' wedge
vectors (f, g, f′, g′; the gauge and (e,e′)). Reconstruct each block A^{(cc′)}_{(mm′),(nn′)} = Σ_k (f_k∧g_k)_{mm′}(f′_k∧g′_k)_{nn′} (s-frame,
exponentials stripped) and compare with the printed κ laws: must agree identically in (s,t) (sympy `simplify`/`cancel` to 0). Also
reproduce one numeric check from cell 31/32 if the data are in the extract (optional). Record the exact conventions of minor, Bk2, ladder
point order and the wedge (F∧G)_{mm′} = F_mG_{m′} − F_{m′}G_m from the setup cells.

**L1 — every h-dependent ingredient, exactly.** (a) kappaTie(z,w) := ∏_k Det[(n(z,w)·P₃(z))_{1..k,1..k}]^{λ_k−λ_{k+1}} (cell 252 of the
untwist dump; identify n = nOf and P₃ = Pphys3 there) at z = e^{hζ}, w = (e^{hω₁}, ∞, w₃) — w₂ = ∞ limit first (symbolically), then the
h-series (exact: `Ser` or sympy series) — its leading power and coefficient; same for κ′. (b) ⟨x₀^d|ψ_d⟩: P2 found ‖x₀‖ ~ h per site and
x₀ → x₀^G projectively; get the exact leading coefficient c₀h²·x₀^G (one-site: the leading term of e₁ᵀP_phys(z)⁻¹ as a series). (c) The
ladder/vacuum pieces as formal series: with jets J^j_{n,k} := (d/du)^k q^j_n(θ_s) as independent symbols (n = 0,1,2 correction order,
k = 0..3, j = 1,2,3) and the exponential e^{ζ_jθ_s} factored, write Q^j(θ_s+mh) = e^{ζ_j(θ_s+mh)}Σ_{n,k} h^{n+k} m^k J^j_{n,k}/k!; expand
f_m(e^{−hδ}, θ₁₂/h) (δ = ω₁′−ω₁) and a^m = e^{−hω₁m}, a^{m+e} etc. as Laurent series in h (sympy; note t = θ₁₂/h makes the kernel coefficients
Laurent in h — keep enough orders: the (s−1) factors are O(h), the 1/(s−1) pole of product 7 is O(1/h)). (d) The direct-side dictionary
is not needed, but record the Gaudin facts to be used for simplification: W₀ := J²_{0,0}J³_{0,1} − J³_{0,0}J²_{0,1} = 0 (Wronskian zero),
J^j_{0,1} = κ_s J^j_{0,0} (Frobenius; D1 derives κ_s — use it as a symbol).

**L2 — the formal limit.** Assemble RHS·κκ′·⟨x₀^d|ψ_d⟩/(M¹M²) as a Laurent series in h with coefficients polynomial in the jets and rational in
(ω₁, ω₁′, θ₁₂, ζ_j, e^{ζθ}) — product by product (10) and block by block (4). Determine the leading order h^N of the full RHS (after all
cancellations; use the Wronskian zero W₀ = 0 as a substitution rule to expose cancellations) and compare with the LHS ⟨Φ_d|ψ_d⟩ = O(h⁰)
(the Gaudin diagonal-frame overlap is finite and ω-dependent). Verdict: (i) N = 0 and the coefficient involves only n = 0 jets → the
Gaudin formula; print it as a bilinear form in {J³_{0,c}(θ₁)·p_{kl}(θ₁)} ⊗ {J³_{0,c′}(θ₂)·p_{k′l′}(θ₂)} (Plücker p_{kl} = J²_kJ³_l − J³_kJ²_l)
with kernel entries as explicit functions of (ω₁, ω₁′, θ₁₂) (and ζ? there should be no ζ in the kernel — if ζ appears, say where), simplified
with W₀ = 0 and J_{0,1} = κ_sJ_{0,0}; also in the direct-side basis {Q₁^{(k)}(θ_s)J³_{0,c}(θ_s)} via W′ = (ζ₃−ζ₂)e^{Σζθ_s}Π′(θ_s)Q₁(θ_s) and its
derivatives. (ii) N ≠ 0 or n ≥ 1 jets enter at the leading order → report exactly which q₁ (q₂) jets enter and with which coefficients;
derive the ODE for q^j_1 from the order-h⁴ term of the XXX dual Baxter expansion (extend P2's X4 symbolic expansion by one order, with
t₁, t₂ expanded to h⁴ as unknown scalar functions A₄, B₄ and the X2 dictionary for A₂, A₃, B₂, B₃ — the A₄, B₄ terms must drop out of
the q₁ equation or be expressible through S_k and their derivatives; say what happens), and STOP there with a clear statement of what is
needed to finish. Keep the computation exact; if sympy is too slow on the full product, work block by block and truncate the jet symbols
to the orders that can contribute (count the h-powers first: a bookkeeping table of the minimal h-order of every factor).

**L3 — numerical test (only in case (i)).** Gaudin states at ζ = P (P1's eigenstates/Q's via import), Φ∞ = Ω^∞[ζ,ω₁]⊗Ω^∞[ζ,ω₁′] with the
plan's r_k^∞ at (ω₁,ω₁′) = (5/2,7/4), x₀^G (P2's construction, basis order mapped to gl3_reps order), LHS Qnum(B) = ⟨Φ∞|Ψ_B⟩/⟨x₀^G|Ψ_B⟩ vs
the limit formula on all 64 states (60 digits): max relative deviation; the state-independent constant between them (the κκ′ and x₀
normalisation constants) — report it and whether it is the one predicted by L1(a,b). If `RESULTS-twovac-fit.md` (agent D1) exists, compare
your kernel with its fitted kernel at the same point (same monomial basis; modulo the null space it reports).

## 2. Reporting

`RESULTS-twovac-limit.md`: L0 (parse verification), L1 (the h-orders table of all ingredients, κκ′ and x₀ leading coefficients), L2 (the
leading order N, the verdict, the Gaudin formula in both bases or the q₁ obstruction and its equation), L3 (numbers). Everything exact
where it can be; numbers with precision. No notebooks, no git, no Wolfram, no edits outside your two files. If an expectation is refuted
and two attempts do not fix it, stop there, keep what works, report what you saw.
