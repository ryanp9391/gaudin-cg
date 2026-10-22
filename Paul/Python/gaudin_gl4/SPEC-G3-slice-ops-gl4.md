# Spec G3: reduction of the gl(4) Gaudin missing-label operators to the PRV slice, and the canonical copies (2026-10-22; final — written after RESULTS-G1 and RESULTS-G2)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram, no wolfbook tools, no git, no notebooks). Work package WP-G3 of
`Paul/Mathematica/Gaudin/gl4-slice-engine-plan.md`. Environment as in SPEC-G1/G2 (Python 3.13, numpy/scipy/sympy/mpmath, no gmpy2; `python -X utf8`; run from
`Paul/Python/gaudin_gl4/`; `perf_counter`; ≤ 2 heavy processes). **Never edit existing files; import G1's and G2's modules read-only.**
Files you own: `slice_ops_gl4.py`, `test_slice_ops_gl4.py`, `RESULTS-G3-slice-ops-gl4.md`, `slice_ops_gl4_output.txt`, `../../Mathematica/Gaudin/ref/gl4-cg-copies-*.txt`.
Another agent (G5, `SPEC-G5-oper-gl4.md`) works concurrently on the oper/Q-system side; it does not touch your files, do not touch its (`oper_gl4.py`, `RESULTS-G5-*`).

Read first: this file; the plan §0–§2 (row "reduction to the slice"); `RESULTS-G1-slice-engine-gl4.md` §0, §1, §12 (usage) and `RESULTS-G2-gaudin-gl4-ops.md` §0, §2, §5, §6, §7, §10
(the facts and the APIs you build on); the gl(3) template `Paul/Python/gaudin_sl3/slice_engine.py` §4 (`vred_apply`, `V_red`, `copies`) and `RESULTS-slice-engine.md` §0;
the term-list file `Paul/Mathematica/Gaudin/ref/gl4-gaudin-ops.txt` (its `#` header documents the format).

## 0. Why this exists, and what G1/G2 established

G1 (`glN_reps.py`, `slice_engine_gl4.py`) gives, for every irrep μ of V_λ⊗V_λ′, the exact slice space M_μ ⊂ V_λ[μ−λ′] (dim m = multiplicity; `prv_space`) and the full
HW vector of any slice (`hw_from_slice`, exact twin `hw_from_slice_exact`). G2 (`gaudin_gl4_ops.py`) gives the three commuting two-site operators
**V₃ = t·S_{3,2}, V₄ᵃ = t·S_{4,3}, V₄ᵇ = t²·S_{4,2}** (t = θ₁₂; S_{k,j} = Laurent coefficients at θ₁ of the column-ordered cdet(∂−L)), as exact term lists in e¹_ij, e²_kl.
Established by G2 (exact, 14 products, 180 sectors with m ≥ 2, m ≤ 12): the triple is **exactly t-independent**; the three operators commute and are self-adjoint for the
invariant form G⊗G on every HW space; **the joint spectrum is simple on every sector**; **V₃ alone is degenerate on some λ = λ′ sectors** (site-exchange parity forces zero
modes: e.g. [3,2,1,0]², μ = (4,4,2,2), m = 2, V₃ = 81·1; μ = (5,4,2,1), m = 4, eigenvalue 57 twice), while {V₃, V₄ᵇ} already separates everything tested; the eigenvalues are
real. Term counts: V₃ 88, V₄ᵃ 152, V₄ᵇ 324; **V₃ and V₄ᵃ have at most ONE site-2 factor per term, V₄ᵇ at most TWO**; the other Laurent coefficients have more (S_{3,1}: 156
terms, S_{4,1}: 412 terms, up to three site-2 factors) and are integer combinations of the triple plus scalars on every HW space (G2 §6c: t²·S_{3,1} = −2V₃ + scalar,
t³·S_{4,1} = −3(V₄ᵃ + V₄ᵇ) + scalar, t²·S_{4,2}^{(θ₂)} = 3V₄ᵃ + 2V₄ᵇ + scalar, …).
What is missing is the step that made gl(3) n_W-free: **the action of each operator on the slice alone**, V_red,k : M_μ → M_μ, from the term list, without any two-site
vector. For gl(3) every term had one site-2 factor (one-line rule); here the rule needs the normal ordering below.

APIs (take the final names from the modules' docstrings): G1 — `gr.GLNRep.of(lam)` (`.E[(i,j)]` columns of {row: Fraction}, `.apply((i,j), {idx: Fraction})`, `.wsp`, `.top`,
`.g`, `.dense((i,j))`, `.float_ops()`), `se.prv_space(lam, lamp, mu)` → `.W .basis .free .m .w .vec(k) .coords(v) .float_basis()`, `se.hw_from_slice(lam, lamp, mu, vfloat)` → HW
(`.rows .cols .vals` (nW, k), plain GT basis), `se.hw_from_slice_exact(lam, lamp, mu, basis_vec)` → ({(ν,q): Fraction}, #conditions, violations), `se.hw_from_slice_rounded`,
`se.hw_residual(hw, extra=True)`, `se.parse_ref`, `se.Recursion.of`, `se.ExactRecursion.of`. G2 — `load_terms(path)` → {name: {(w₁, w₂): {power: Fraction}}} with w₁ the
site-1 word and w₂ the site-2 word (tuples of (i,j), product order, leftmost factor acts last; e¹ and e² commute); `HWExact(R1, R2, mu)` (exact HW basis of a sector),
`triple_ops(R1, R2)` / the Kronecker-sum builder (apply a term list to a two-site vector), `restrict_exact`, `analyse_sector` (exact certificates, joint eigenvalue triples);
its symbolic stages (`cdet_S`, `laurent`, …) are general in N and can produce the N = 3 term list of t·S_{3,2} for the convention test.

## 1. The reduction theorem (state it in RESULTS in your own words, then implement it)

Setting (G1 conventions): Ψ = Σ_q ψ_q ⊗ e_q a HW vector of weight μ in V_λ⊗V_λ′, slice ψ_top = ψ_{q_top} ∈ M_μ ⊂ V_λ[μ−λ′]; ⟨top| the coefficient functional of e_top on
site 2 (slice(Φ) = (id ⊗ ⟨top|)Φ). Facts: (i) ⟨top| e²_{ab} = 0 for a > b (a lowering operator never produces the top state); (ii) ⟨top| e²_{aa} = λ′_a ⟨top|; (iii) for a < b,
e²_{ab} Ψ = −e¹_{ab} Ψ (the HW condition E^{tot}_{ab}Ψ = 0), and e²_{ab}(W¹Ψ) = W¹e²_{ab}Ψ = −W¹e¹_{ab}Ψ for any site-1 word W¹ (sites commute).
Hence for a term c·W¹·W² (site-1 factors collected to the left in their original relative order, site-2 factors to the right in theirs): normal-order the site-2 word W² in
U(gl(N)) as Σ c_w L·H·R (L lowering, H Cartan, R raising; fixed order inside each class) with [E_ij,E_kl] = δ_jk E_il − δ_li E_kj; then
**slice(c W¹ L H R Ψ) = 0 if L ≠ ∅; otherwise = c · Π_{E_aa∈H} λ′_a · (−1)^{|R|} · W¹ e¹_{r_|R|} ⋯ e¹_{r_1} ψ_top for R = e²_{r_1}⋯e²_{r_|R|}**
(transferred factors in REVERSED order, to the right of W¹: the rightmost site-2 factor acts on Ψ first and each transfer passes through the site-1 operators already
produced). So V_red,k = Σ_terms (site-1 word, scalar coefficient in ℚ[λ′]) — a universal polynomial in site-1 generators, independent of λ, μ except through its action on
V_λ[μ−λ′]; it maps M_μ into M_μ (the operators preserve HW vectors) and the three V_red,k commute on M_μ. Products of site-1 words act right to left on vectors.

## 2. Tasks and acceptance (every number goes into RESULTS-G3)

**G3.1 Normal-ordering engine** `normal_order(word)` for words of gl(N) generators (tuples of (i,j); class order lowering < Cartan < raising; PBW reordering by commutators,
memoised; returns {word: coefficient}). Test: ≥ 200 random words of length ≤ 4 in gl(4): the reordered sum equals the original word as an exact matrix on [2,1,0,0] and
[1,1,0,0] (`GLNRep`, exact); plus the length-3 site-2 words occurring in S_{4,1}.

**G3.2 Reduction** `reduce_terms(terms, lamp)` → site-1 term list (exact Fractions), `vred_apply(rep1, red_terms, v)` (exact, dict vectors), `V_red(lam, lamp, mu, name)` →
exact m×m matrix in G1's PRV basis (column j = V_red b_j in PRV coordinates via `PRV.coords`, which raises if the image leaves M_μ). **Convention test on N = 3:** the N = 3
term list of t·S_{3,2} (from G2's symbolic engine) reduced by your engine must equal `slice_engine.vred_apply` + (η₁/2)·1 on every PRV basis vector of every irrep of
[2,1,0]², [3,1,0]², [4,2,0]² — exactly (η₁ = |Λ|² + cas₂(Λ) − cas₂(μ)/2; this is θ₁₂a₁ = V + η₁/2). Also reduce the N = 3 S_{3,1}, which has two-factor site-2 words, and check
t²·S_{3,1,red} = −2·V_red(t S_{3,2}) + scalar on M_μ (G2 §6c pattern holds for N = 3 too — verify the integer by the direct route if unsure).

**G3.3 Direct-route verification on gl(4):** for every irrep of [2,1,0,0]², [2,1,1,0]², [2,2,1,0]⊗[2,1,1,0], [1,1,0,0]⊗[2,1,0,0], [3,2,0,0]⊗[2,1,1,0], and all 16 irreps of
[3,2,1,0]² (its 8 sectors with m ≥ 2 included): reconstruct the full HW vectors of the PRV basis (exact twin where dim·dim′ ≤ 4096, else float), apply the full two-site
operator (G2's builder, term list from the file), read the slice, compare with V_red b_j: **exactly equal / ≤ 1e-12 relative**, for V₃, V₄ᵃ, V₄ᵇ AND for S_{3,1}, S_{4,1} (the
multi-factor ones, at t = 1: they test the general rule). Then on M_μ: V_red,k preserves M_μ; the three commute (exact); their joint spectrum = G2's joint eigenvalue triples
(exact characteristic polynomials / 30-digit eigenvalues); t²·S_{3,1,red} = −2V_{3,red} + c·1 and t³·S_{4,1,red} = −3(V_{4a,red} + V_{4b,red}) + c′·1 with scalars c, c′ (report them).

**G3.4 Canonical copies** `copies(lam, lamp, mu, dps=50)`: joint eigenvectors of (V_red,3, V_red,4a, V_red,4b) on M_μ by eigen-decomposition (mpmath) of a generic rational
combination (V₃ alone is NOT enough — see §0), then verify every V_red,k is diagonal in that basis to 10^{−(dps−10)}; order copies lexicographically by (ev₃, ev₄ᵃ, ev₄ᵇ)
ascending; slice normalised to first nonzero component (ascending GT index) = 1; full vectors by G1's recursion. Checks: real eigenvalues; each V_red,k self-adjoint for the
restricted invariant form on M_μ (G-Gram of the PRV basis; report); G⊗G-orthogonality of the copies' full vectors (Gram off-diagonal ≤ 1e-12 relative; distinct joint
eigenvalues + self-adjointness make this a theorem — confirm numerically); agreement with G2's joint eigenbasis up to scale (slices of G2's exact HW eigenvectors vs yours) on
all sectors of [3,2,1,0]² and [2,1,0,0]². For λ = λ′ tabulate the exchange structure (V₃ − c/2 spectrum symmetric; which copies are paired).

**G3.5 Timing:** t_a = n_W-free part (PRV spaces + the three reduced operators on every M_μ + joint eigen, all irreps; count the cost of `reduce_terms` once per λ′ separately
from the per-sector application) and t_b = + full reconstruction of all copies, for: the scaling family [2a,a,0,0]², a = 1..6 (m = a+1; a = 6: dim 6860, 322 irreps, Σm = 728,
max n_W 29 028); the stress family [3a,2a,a,0]², a = 1, 2 (a = 2: dim 729, 76 irreps, Σm = 297, max m = 12, max w = 15, max n_W 4743) and a = 3 slices only if t_a < 2 min
(dim 4096, 213 irreps, Σm = 1592, max m = 27, w ≤ 32); and [S,S,0,0]², S = 4, 8 (m = 1 everywhere: the operators are not needed, but report that V_red,3 on each M_μ is then a
scalar and what it is). Table vs a: #irreps, Σm, max m, max w, max n_W, t_a (with breakdown), t_b. State honestly where the reduction cost sits relative to the gl(3) 24-term V_red.

**G3.6 Reference files** `ref/gl4-cg-copies-<l1l2l3l4>[-<l1'..>].txt` for [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [2,2,1,0]⊗[2,1,1,0], [4,2,0,0]², [6,4,2,0]² (slices only for the last):
`#` header (conventions, the reduction rule, the operator definitions with G2's normalisation, ordering/normalisation of copies, script, date); per sector (descending μ):
`SECTOR mu m w nW`, `PRVBASIS` × m (as G1), `VRED3`, `VRED4A`, `VRED4B` (m×m exact in the PRV basis), `CHARPOLY3` (exact char poly of VRED3, ascending powers) and `CHARPOLYJ`
(of the generic combination used, state it), `COPY mu k ev3 ev4a ev4b {{nu,val},...}` (30 digits), `HWRESID`, `GNORM2`, and `VEC` where dim·dim′ ≤ 2000. `verify_ref`: re-parse
from the text and check with G1's and G2's builders independently of your reduction code (string conditions, the full two-site operators on the reconstructed copies, HW
residuals, GNORM2): 0 mismatches at 1e-12.

**G3.7 `RESULTS-G3-slice-ops-gl4.md`** (style of G1/G2): files, the theorem as implemented (sign/order rule), the N = 3 convention test, the gl(4) direct-route table (max
deviation per product and operator), commutativity/simplicity on M_μ, the scalars c, c′, the self-adjointness/G-orthogonality findings, the exchange structure, the timing
table, surprises, not done. `test_slice_ops_gl4.py`: a few-second exact regression test (normal ordering, N = 3 convention, one gl(4) sector direct route, one copies check).
CLI stages: `order | reduce | gl3 | direct | copies | bench | ref | verify | all`. Report failures with numbers; never paper over a failing test.

## 3. Pitfalls

* Site-2 Cartan factors evaluate to λ′ (site 2's HW); the weight space and the transferred operators live on site 1 (V_λ) — do not mix λ and λ′ (for λ ≠ λ′ G1 puts the factor
  with the smaller weight spaces at site 1; the file's conventions fix which is which).
* Terms whose site-2 word has an odd number of transferred factors pick up a sign; a lowering factor anywhere in the normal-ordered word kills the term, but NOT before normal
  ordering (e²_{21}e²_{12} = e²_{12}e²_{21} − e²_{11} + e²_{22} contributes through its Cartan part). Test the whole rule on N = 3 (G3.2) before touching N = 4.
* Keep everything exact until the eigen-decomposition; compare float results relative to the vector norm in the orthonormal basis.
* `PRV.coords` is the M_μ-membership test — use it, and count how often `reduce_terms` output words act trivially on a weight space (sparsity report).
* Tooling trap (G1, G2 both hit it): inside Bash heredocs `\t`, `\n` in Python source get un-escaped — write files with the file tools.
