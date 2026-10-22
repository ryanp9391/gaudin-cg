# Spec G5: the oper side of the gl(4) missing labels — the dictionary (V₃, V₄ᵃ, V₄ᵇ) ↔ the 4th-order Fuchsian oper, its polynomial solutions, the independent Q-system solve, and the July-2026 gl(4) Gaudin ODE (2026-10-22)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram kernel, no wolfbook tools, no git, no notebooks — the one Mathematica notebook you need is read
as TEXT). Work package WP-G5 of `Paul/Mathematica/Gaudin/gl4-slice-engine-plan.md`. Environment as in SPEC-G1/G2 (Python 3.13, numpy/scipy/sympy/mpmath, no gmpy2;
`python -X utf8`; run from `Paul/Python/gaudin_gl4/`; `perf_counter`; ≤ 2 heavy processes). **Never edit existing files; import G1's and G2's modules read-only.**
Files you own: `oper_gl4.py`, `test_oper_gl4.py`, `RESULTS-G5-oper-gl4.md`, `oper_gl4_output.txt`, `../../Mathematica/Gaudin/ref/gl4-oper-dictionary.txt`.
Another agent (G3, `SPEC-G3-slice-ops-gl4.md`) works concurrently on the slice reduction; do not touch its files (`slice_ops_gl4.py`, `RESULTS-G3-*`, `ref/gl4-cg-copies-*`).

Read first: this file; the plan §0 (row "integrability cross-check") and WP-G5; `RESULTS-G2-gaudin-gl4-ops.md` §0, §2, §4, §5, §6 (esp. 6c, the table "the triple generates
everything"), §7, §10 (hand-off notes); the gl(3) conventions and oper facts in `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.md` §1.1–1.2 and §2.1–2.4 (oper, Q-functions,
local exponents, no-log conditions, Wronskians) and its "untwisted facts" (search the file for "untwisted": Fuchsian at ∞, three polynomial solutions of degrees μ + (2,1,0),
a₁ = −a₂ = (μ₁+2)(μ₂+1)μ₃/(2θ₁₂)); the gl(3) Q-system solver report `Paul/Mathematica/Gaudin/Experiments/Gaudin-Qsystem-scale.md` §0–§2 (the pencil method and ITS BUG: the
GCD of two of three consistency polynomials over-counts; all of them are needed) and `Paul/Python/gaudin_sl3/RESULTS-jacobi-cg.md` §3 (which Baxter system alone has the right
degree); the July-2026 gl(4) block of `HANDOFF.md` ("State as of 2026-07-22 (gl(4) L=2 untwisted τ1,τ2,τ3 ALL DONE …)", around line 1990) and the notebook
`Paul/Mathematica/old-conventions/old/Gaudin/Clean/GaudinSolver.wb` (JSON text: grep for `gaudinODE`, `degQ`, `muhat`, `Cas2`, `Cas4`, `wsub`, `irreps`).

## 0. What is established and what this package adds

G2 built the Laurent coefficients S^{(1)}_{k,j} (k = 2,3,4; j = 1..k) at θ₁ of the column-ordered cdet(∂_u − L(u)) = ∂⁴ − S₁∂³ + S₂∂² − S₃∂ + S₄ of the untwisted two-site gl(4)
Gaudin Lax matrix L(u) = E⁽¹⁾/(u−θ₁) + E⁽²⁾/(u−θ₂) (t = θ₁₂; x = u − θ₁), exactly, as term lists (`ref/gl4-gaudin-ops.txt`), proved (exactly, on 14 products) that they commute,
and that on every HW space M_μ of V_λ⊗V_λ′ each of them is an INTEGER combination of the t-free triple **V₃ = t S_{3,2}, V₄ᵃ = t S_{4,3}, V₄ᵇ = t² S_{4,2}** plus a scalar
(G2 §6c; the θ₂ coefficients by site swap and t → −t are included in that table), that the triple has simple joint spectrum on every sector, and it tabulated the joint
eigenvalue triples (18 digits, `gaudin_gl4_ops_output.txt`; exact certificates via `analyse_sector`). Facts you may use: S_k(u) has poles of order ≤ k at θ₁, θ₂ only and
S_k = O(u^{−k}) at ∞, so **S_k(u) = Σ_{s=1,2} Σ_{j=1}^{k} S^{(s)}_{k,j}/(u−θ_s)^j exactly**; S₁ = n₁/(u−θ₁) + n₂/(u−θ₂) with n₁ = |λ|, n₂ = |λ′|; S^{(1)}_{k,j} = t^{j−k}·Ŝ_{k,j};
the top poles are central scalars (S^{(1)}_{k,k} = (−1)^k a_{N−k} with Π_c(θ − y_c) = Σ_k a_k θ(θ−1)…(θ−k+1), y_c = λ_c + N − c; e.g. S_{2,2} = e₂(λ) + Σ_k (k−1)λ_k);
t·S^{(1)}_{2,1} = n₁n₂ − Ω with Ω = Σ e¹_ij e²_ji = (C₂(μ) − C₂(λ) − C₂(λ′))/2 on weight μ, C₂(λ) = Σ_i λ_i(λ_i + N + 1 − 2i).
On a joint eigenvector with eigenvalues (v₃, v₄ᵃ, v₄ᵇ) the cdet becomes a SCALAR 4th-order Fuchsian differential operator (the oper) with singular points θ₁, θ₂, ∞ and
three accessory parameters — the three missing labels. For gl(3) the analogous dictionary was V = a₁ − η₁/2 (a₁ = θ₁₂ × double-pole coefficient of S₃), the untwisted oper
had polynomial solutions of degrees μ + (2,1,0), and the Q-system consistency polynomial equalled the characteristic polynomial of V (`qsSolve`, Jacobi chain). **This package
does the gl(4) version:** (1) the complete dictionary (all 18 Laurent coefficients, scalars included, as functions of (v₃,v₄ᵃ,v₄ᵇ; λ, λ′, μ)); (2) exponents and polynomial
solutions for every copy; (3) the INDEPENDENT solve: the set of (v₃,v₄ᵃ,v₄ᵇ) for which the oper has polynomial solutions of the right degrees must be exactly the m joint
eigenvalue triples; (4) the comparison with the July-2026 gl(4) Gaudin ODE of `GaudinSolver.wb`.

## 1. Tasks and acceptance (every number goes into RESULTS-G5)

**G5.1 The dictionary.** For the products [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [2,2,1,0]⊗[2,1,1,0], [1,1,0,0]⊗[2,1,0,0], [4,2,1,0]² (m ≤ 5) and every irrep μ (m ≥ 1): restrict all
18 Laurent coefficients (θ₁ and θ₂; S₂, S₃, S₄) to M_μ (G2's `HWExact`, `restrict_exact`; or your own exact null space), subtract the integer combination of the restricted triple
from G2's §6c table, and verify the remainder is a scalar matrix — record the scalar **c^{(s)}_{k,j}(λ,λ′,μ)** exactly (Fractions). Fit closed forms: in terms of the Casimir
eigenvalues C₂, C₃, C₄ of λ, λ′, μ (the trace forms G2 uses: C_k = Σ e_{i₂i₁}e_{i₃i₂}⋯e_{i₁i_k}; compute their eigenvalues on V_λ on the HW vector with `GLNRep`), n₁, n₂, Ω and t
(powers fixed by homogeneity): state which scalars you pinned to closed forms (verify each on ALL sectors, ≥ 60 sectors) and which remain tabulated. Write the dictionary to
`ref/gl4-oper-dictionary.txt` (documented format: per operator the integer coefficients of (V₃,V₄ᵃ,V₄ᵇ), the power of t, the scalar formula or the per-sector table).

**G5.2 Exponents and polynomial solutions.** With the oper D(v) := ∂⁴ − S₁∂³ + S₂∂² − S₃∂ + S₄ built from the dictionary for a sector and a triple v = (v₃,v₄ᵃ,v₄ᵇ):
(i) the local exponents at θ₁ are {λ_c + 4 − c}_{c=1..4} and at θ₂ {λ′_c + 4 − c} (indicial polynomial from the most singular terms; verify symbolically in v — they must not
depend on v); (ii) at ∞ the indicial polynomial is I(n) = Π_c (n − (μ_c + 4 − c)), so polynomial solutions can only have degrees d_j = μ_j + 4 − j (verify); (iii) for EVERY copy
of every sector listed (joint triple from G2, exact if rational else 40 digits): D(v) has a polynomial solution of degree exactly d_j for each j = 1..4 — the null space of D(v)
on polynomials of degree ≤ d_j has dimension #{k : d_k ≤ d_j} (exact rank where v is rational; mpmath SVD at 40 digits with a gap ≥ 10²⁰ between kept and dropped singular values
otherwise); (iv) the Wronskian of the four polynomial solutions is ∝ (u−θ₁)^{|λ|}(u−θ₂)^{|λ′|} (the gl(3) analogue is W ∝ Π^{|Λ|}); (v) negative tests: a random triple, and each
known triple perturbed by 10⁻³ in one component, admit NO polynomial solution of degree d₁ (smallest singular value ≥ 10⁻⁶ relative). **Acceptance: (iii) holds for all copies
(≥ 150), (v) rejects all perturbations.**

**G5.3 The independent solve (the point of the package; the gl(3) analogue is Q0's `qsSolve`).** Treat v as unknown. D(v) is AFFINE in v. After multiplying by
(u−θ₁)⁴(u−θ₂)⁴ the condition "a polynomial solution of degree ≤ d_j exists" is a rank condition on an affine matrix pencil M_j(v) of size (d_j + 5) × (d_j + 1) (explain the
sizes). Find the set Σ = {v : all four degree conditions hold, with the null-space dimensions of (iii)}; **expected: Σ = exactly the m joint triples of G2** (completeness of the
Bethe/oper description; the two-site case is affinely generic). Method, exact where feasible: the Frobenius recursion at ∞ (p_{d_j} = 1; the coefficient of u^{d_j+4−k}
determines p_{d_j−k} unless d_j − k is another exponent d_i — a resonance: one consistency polynomial in (v, free parameters) and one free parameter; plus the bottom
equations), eliminate the free parameters (they enter linearly), intersect over j by Gröbner bases / resultants (sympy; time-box 10 min per sector), and compare with G2's
triples (≥ 20 digits). Do this exactly on all sectors of [2,1,0,0]², [2,1,1,0]², [2,2,1,0]⊗[2,1,1,0], [1,1,0,0]⊗[2,1,0,0] and the m ≤ 3 sectors of [3,2,1,0]². For the larger
sectors ((5,4,2,1) of [3,2,1,0]² with m = 4 — where V₃ is degenerate, so this is the interesting case — and the m = 4, 5 sectors of [4,2,1,0]²): if the exact elimination is out
of reach, certify at least that each of G2's triples is an isolated, nondegenerate solution (Jacobian of the condition system has rank 3 there) and search for further
solutions numerically from ≥ 200 random real starts in the box |v_i| ≤ 3·max|spectrum| (Newton/Levenberg–Marquardt on the smallest singular values or on the consistency
polynomials); report what was certified and what was not. **Also report the gl(4) version of the pencil-bug analysis:** for each small sector, how many solutions each SINGLE
degree condition j has (the condition for P₁ alone, P₂ alone, …) — i.e. which subsets of the four conditions already cut out exactly the spectrum. **Acceptance:** #Σ = m with
Σ = G2's triples on every exactly-solved sector; the certification table for the rest.

**G5.4 The July-2026 gl(4) Gaudin ODE.** `GaudinSolver.wb` (JSON; read the cell sources as text, e.g. with Python's `json` module or grep) holds `gaudinODE`: the h → 0 limit of
the XXX Baxter equation for the multiplicity-free products [S,S,0,0]⊗[P,P,0,0], a 4th-order ODE in z for the DIRECT Q-functions Q₀ with polynomial solutions of degrees
`degQ[S,P,μ]_k = (S+P+3) − muhat_k`, muhat_k = μ_k + 4 − k (`degQ` is in the notebook), coefficients in `C2 = Cas2[μ]`, `C4 = Cas4[μ]`, `S`, `P`, the sites `w1`, `w2` (the
substitution `wsub`), leading coefficient 8(z−w₁)²(z−w₂)² (HANDOFF). Transcribe it to sympy (document the transcription; the notebook's own definitions of `Cas2`, `Cas4`,
`muhat` fix the Casimir conventions — transcribe those too) and (a) reproduce the notebook's claim: polynomial solutions of degrees `degQ` exist for every irrep of
[S,S,0,0]⊗[P,P,0,0] with (S,P) ∈ {(1,1), (2,1), (2,2), (3,2), (3,3)} and only for those degrees; (b) relate it to OUR oper for the same sectors (all m = 1, so our (v₃,v₄ᵃ,v₄ᵇ) are
the Casimir scalars of G2's table): the expected relation is **formal adjoint + gauge** (gl(3): direct oper = Π^{p} ∘ (formal adjoint) ∘ Π^{−p}, p = λ₁ − λ₃; here try
(z−w₁)^{p₁}(z−w₂)^{p₂} with p_s ∈ {S or P, S−1, …} — scan small integers) and an overall factor: determine p₁, p₂ and the dictionary between their (Cas2, Cas4) and our scalars,
and verify the two operators coincide as differential operators (coefficient by coefficient, exact) on all tested sectors; the degrees must map as polynomial solutions of the
adjoint (our degrees μ_j + 4 − j vs their (S+P) − μ_k + k − 1). If the transcription cannot be made reliable, do the degree-level comparison only and say so.

**G5.5 `RESULTS-G5-oper-gl4.md`** (style of G1/G2): files, the dictionary (table + closed forms found), the exponent/degree verification (counts), the independent-solve table
(sector, m, #solutions, agreement digits, Jacobian ranks, single-condition counts), the July-2026 comparison (p₁, p₂, the Casimir dictionary, agreement), surprises, not done.
`test_oper_gl4.py`: a few-second regression test (one sector dictionary, one copy's polynomial solutions, one exact solve). Output of all runs in `oper_gl4_output.txt`. CLI
stages: `dict | exps | solve | wb | all`. Report failures with numbers; never paper over a failing test. Time-box the whole package to about 2 hours of wall clock; if G5.3 or
G5.4 is only partly done, deliver the rest in full and say what is missing.

## 2. Pitfalls

* θ₂ data: site swap AND t → −t (G2 §2); S₂'s simple poles are opposite: S^{(2)}_{2,1} = −S^{(1)}_{2,1}; all residue relations (S_k = O(u^{−k})) must hold for your oper — assert them.
* The oper's coefficients are operators on the HW space; it becomes scalar only on a joint eigenvector — never build it from a non-eigenvector's expectation values.
* The exponents at ∞ are integers and all resonate: logs are the generic case, polynomial solutions the exceptional one — that is the content of the Bethe/no-log conditions.
* Gauge freedom: a polynomial solution of degree d_j is defined modulo the lower-degree ones (the notebook's "gauge fixing"); count null-space dimensions, not vectors.
* Precision: use exact arithmetic wherever the triple is rational (many m = 1 and some m = 2 sectors); otherwise mpmath ≥ 40 digits and explicit gap criteria.
* Tooling trap (G1, G2 both hit it): inside Bash heredocs `\t`, `\n` in Python source get un-escaped — write files with the file tools.
