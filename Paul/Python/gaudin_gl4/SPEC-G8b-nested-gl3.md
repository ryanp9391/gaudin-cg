# Spec G8b: nesting gl(3) ⊂ gl(4) — the paths of the weight space carry the gl(3) Jacobi chains; the nested basis; the depth-1 sectors in closed form (2026-10-24)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram, no wolfbook tools, no git, no notebooks). Work package G8b of
`Paul/Mathematica/Gaudin/gl4-jacobi-plan.md` (read its §0 and §1 H-B first). Environment: Python 3.13, numpy, scipy, sympy, mpmath; run everything from
`Paul/Python/gaudin_gl4/` as `python -X utf8 ...`; at most two heavy processes at any time; `time.perf_counter()`. **Never edit existing files; import G1–G5 and the gl(3)
code read-only.** Files you own: `nested_gl3_gl4.py`, `test_nested_gl3_gl4.py`, `RESULTS-G8b-nested-gl3.md`, `nested_gl3_gl4_output.txt`, `nested_gl3_gl4_tables.txt`.
Two other agents (G8a, G8c) work concurrently in the same directory on their own files — do not touch `lattice_recursion_gl4*`, `relabelled_gl4*`. Time-box ≈ 2.5 hours;
deliver what is finished and say precisely what is not.

Tooling trap (G1–G4 all hit it): Python source written through Bash heredocs loses `\n` escapes. Write code with the file tools and check `python -c "import ast; ast.parse(open('nested_gl3_gl4.py').read())"`.

## 0. Setting (binding conventions; see RESULTS-G4 §1, §3.3, §5.1 and the docstring of `glN_reps.py`)

GT pattern of gl(4) = tuple of 4 rows, rows[0] = (m₁₁,), rows[1] = (m₁₂, m₂₂), rows[2] = (m₁₃, m₂₃, m₃₃) =: κ, rows[3] = λ. `R = gr.GLNRep.of(lam)` (`gr = glN_reps`):
`.pats`, `.index`, `.wsp[weight]` (global indices), `.E[(i,j)]` (list over columns of {row: Fraction}), `.apply(ij, vec)`. The weight space of a sector W = V_λ[ν],
ν = μ − λ′ (`R.wsp[nu]`), decomposes as W = ⊕_κ P_κ, κ running over the row-3 patterns that occur, P_κ = the patterns of W with that row 3, ordered by m₁₂ (a path).
Fact to verify first (B0): the gl(3) GT formulas for E₁₂, E₂₁, E₂₃, E₃₂ involve rows 1–3 only, so **the κ-block of these generators equals the gl(3) GLNRep of κ**:
`R3 = gr.GLNRep.of(kappa)` (N = 3) has patterns (rows[0], rows[1], κ) and the map (rows[0], rows[1]) ↔ (rows[0], rows[1], κ, λ) identifies P_κ with the gl(3) weight space
V_κ[ν₁₂₃] (ν₁₂₃ = (ν₁,ν₂,ν₃)) and the matrix elements EXACTLY (check all four generators on 3 reps). The gl(4) weight ν has ν₄ = |λ| − |κ| =: n (the depth).

The reduced gl(4) triple on W: `ops = {n: so.reduced_op(n, lamp) for n in so.TRIPLE}`, `mats = so.weight_space_mats(R, ops, W)` (exact w×w, M[row][col] = coefficient
of e_row in V e_col; `so = slice_ops_gl4`, TRIPLE = ('V3','V4a','V4b')). The gl(3) reduced vertex operator for ANY λ′₃ (not only λ = λ′): `so.reduced_op('V', lamp3)` with
len(lamp3) = 3 (name 'V' = t S_{3,2} of G2's symbolic engine for N = 3; it equals `slice_engine.vred_apply` + η₁/2, G3 §1, and `structure_gl4.FAM[3]` uses it), acting on
`gr.GLNRep.of(kappa)` through `so.weight_space_mats(R3, {'V': op3}, W3)`. The gl(3) closed-form chain for λ = λ′ (reference only): `Paul/Python/gaudin_sl3/jacobi_cg.py`
(`chain(lam, mu).alpha/.u/.d`, `slice_from_V`); for the structural comparison use the general-λ′ operator above, not the λ = λ′ chain.
PRV: `se.prv_space(lam, lamp, mu)` (`.W .basis .free .m .w`), strings E_i^{k_i} with k_i = λ′_i − λ′_{i+1} + 1 (`se.power_images(R, (i,i+1), k, W)`, `se._rref_nullspace`);
copies `so.copies(lam, lamp, mu, dps=50)` (`ev`, `slice` {global index: mpf}); whole-space spectra `structure_gl4.ident_sector`. gl(3) multiplicities for the κ-sectors:
`gr.racah_speiser(kappa, lamp3)` (dict μ₃ → m) — the gl(3) sector is (κ, λ′₁₂₃, μ₁₂₃) with μ₁₂₃ = (μ₁, μ₂, μ₃) (its weight is ν₁₂₃ = μ₁₂₃ − λ′₁₂₃ ✓).

Products: P_small = [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², [3,2,1,0]⊗[4,2,1,0] (both orders), [4,2,2,0]², [4,3,1,0]²; P_big = [6,4,2,0]² (76 sectors,
w ≤ 15, m ≤ 12), [5,3,1,0]², [5,3,2,1]². Also three generic rational non-dominant λ′ (`structure_gl4.GENERIC_LAMPS`) for the structural identities of B1 (the reduction is
polynomial in λ′ and the triple commutes on every weight space for every λ′ — G3 §7b).

## 1. Questions

**B0. The identification** (above) — exact, 3 reps, all four generators, every κ-block. Also: the within-path tridiagonality of all three operators on every P_κ (should be
automatic from the stencil; count the nonzero entries with |Δm₁₂| ≥ 2 — must be 0).

**B1. Is the within-path block the gl(3) chain?** For every weight space of P_small ∪ P_big, every κ, every operator k ∈ TRIPLE: J^{(k)}_κ := the κ-diagonal block of M_k
(rows and columns in P_κ), versus V^{(3)}_κ := the matrix of `so.reduced_op('V', lamp[:3])` on V_κ[ν₁₂₃] in the same (identified) basis. Test, exactly:
  (i) the 1-parameter ansatz J = a·V^{(3)} + b·1: the ratios of corresponding off-diagonal entries are all equal (= a(κ)) and the diagonal difference J − a V^{(3)} is
      scalar (= b(κ)); (ii) if (i) fails, the ansatz J = a·V^{(3)} + b·1 + c·D where D = diag of the gl(2)₁₂ Casimir ξ = x(x+2), x = m₁₂ − m₂₂ (the only natural diagonal
      besides 1 on a path), and (iii) the weakest: the off-diagonal PRODUCTS J[k+1,k]J[k,k+1] are a(κ)²·(the products of V^{(3)}) (diagonal-similarity invariant).
  Report, per operator: in how many (weight space, κ) pairs (with |P_κ| ≥ 2 resp. ≥ 3 so that the test is non-trivial) each ansatz holds, exactly. Where (i) holds for V₃:
  fit closed forms for a(κ), b(κ) as rational functions of (κ, λ, λ′, ν) (sympy interpolation on lines, or guess from G4's factor structure: Q₆ of the (±1,0,0) class
  is a polynomial in (κ, λ) only, RESULTS-G4 §3.3; test the guess a(κ) ∝ Q₆-ratios), and VERIFY on held-out reps ([5,3,2,1]², [5,3,1,0]², the generic λ′). Where it fails,
  show the smallest failing example in full (the two tridiagonal matrices). For V₄ᵃ, V₄ᵇ the same; given G4's M^{V₄ᵃ} = (y_k+3)M^{V₃} on the hopping classes, also test
  J^{(4a)}_κ = (Σ_k (y_k + 3) − something)·J^{(3)}_κ + const, i.e. whether J^{(4a)} ∈ span(J^{(3)}, 1) and J^{(4b)} ∈ span(J^{(3)}, 1, D).

**B2. The nested subspace and the nested basis.** M̃ := ker E₁₂^{k₁} ∩ ker E₂₃^{k₂} ⊂ W (void strings skipped as in `prv_space`). Verify exactly: (a) M̃ = ⊕_κ M^{(3)}_κ where
M^{(3)}_κ ⊂ P_κ is the gl(3) PRV space of (κ, λ′₁₂₃, μ₁₂₃) (compute it with the N = 3 machinery on `GLNRep.of(kappa)` and map it into W); dim M̃ = Σ_κ m₃(κ) with
m₃ = `gr.racah_speiser(kappa, lamp3)[mu123]` (0 if absent); (b) the triple preserves M̃ (BGG) and M_μ ⊂ M̃ is spanned by the m physical eigenvectors (the joint spectrum
of the triple on M̃ contains Spec(M_μ); the other dim M̃ − m eigenvectors violate E₃₄^{k₃}); (c) the nested basis: for each κ the gl(3) canonical copies of
(κ, λ′₁₂₃, μ₁₂₃) — eigenvectors of V^{(3)}_κ on M^{(3)}_κ, 50 digits, normalised first nonzero = 1 (build them with `so.copies` on the gl(3) sector if that works for N = 3,
else by eigen-decomposition of the exact m₃×m₃ PRV matrix as `so.joint_diag` does) — mapped into W: |κ, s⟩, s = 1..m₃(κ). Then compute in that basis (B = matrix of the
nested vectors, T_k = B⁺ M_k B restricted to M̃, 50 digits):
  * the κ-diagonal blocks of T_k: diagonal iff B1(i) holds for that κ (report the off-diagonal norm);
  * the hopping blocks T_k[κ′, κ] for κ′ = κ ± root: rank (vs min(m₃(κ), m₃(κ′))), the fraction of entries above 1e-30, and the rank-one test (is T[κ′,κ] =
    f(κ′,s′)·g(κ,s)? compare with the products of the END components of the gl(3) copies: the last component (largest m₁₂) and the first, since hopping between paths
    enters through the ends of the chains in the lattice basis — test T[κ′,κ]_{s′s} ∝ p^{(κ′)}_{end}(s′)·p^{(κ)}_{end}(s) for each of the four end-pairs);
  * the physical copies in the nested basis: c(κ, s) (50 digits), tabulated for all sectors of [3,2,1,0]², [4,2,1,0]², and the sectors of [6,4,2,0]² with m ≥ 6; the rank of
    the matrix c[κ][s] (rows κ, columns s with missing → 0); any visible structure (zeros, equal magnitudes, factorisation).
  * the relation between the gl(4) eigenvalue v₃ of a physical copy and the gl(3) eigenvalues {v^{(3)}(κ,s)} it is built from: e.g. is v₃ − Σ_κ,s |c|²-weighted average a
    simple function? (Only if something obvious shows up; do not fit blindly.)

**B3. Depth structure and the depth-1 sectors.** Tabulate, per product, for every sector: n = ν₄ = μ₄ − λ′₄, the number of paths #κ, w, m, dim M̃. Confirm n = 0 ⇒ one
path and the gl(4) problem = the gl(3) problem (the triple restricted equals (a V^{(3)} + b, …) with G5's μ₄ = 0 facts: V₄ᵃ = V₄ᵇ = 0 when μ₄ = 0; note n = 0 means
μ₄ = λ′₄, not μ₄ = 0, so say what V₄ᵃ, V₄ᵇ are on n = 0 sectors with λ′₄ > 0). For n = 1 (κ ∈ {λ₁₂₃ − e_i}, ≤ 3 paths): write out the structure completely for all depth-1
sectors of [3,2,1,0]², [4,2,1,0]², [6,4,2,0]²: the paths, the within-path chains (B1), the hopping blocks between the three paths in the LATTICE basis (they are given by
G4's universal product formula for the Δm₁₂ = ±1 classes and the cubic-Q class for Δm₁₂ = 0 — evaluate `structure_gl4.universal_F('V3', pt, d)` on these entries and confirm
against the matrices), and a solution strategy: with B1 true, the depth-1 eigenvector equations are three coupled gl(3) chains; say whether the physical copies at depth 1
are determined by (the gl(3) copies of the three κ) + a small linear problem whose size you state, and implement it as `depth1_copies(lam, lamp, mu)` if it is that simple.
Validate against `so.copies`. If B1 is false, do the same with the actual within-path blocks J_κ (their eigenvectors are then "intermediate labels" instead of the gl(3)
copies) — still a nested construction; say so.

## 2. Deliverables and acceptance

* `nested_gl3_gl4.py` with a CLI `python -X utf8 nested_gl3_gl4.py ident | within | nested | depth | all` (and `run <stage,...>`), tables to `nested_gl3_gl4_tables.txt`,
  printed output of the final run to `nested_gl3_gl4_output.txt`.
* `test_nested_gl3_gl4.py`: ≥ 6 exact regression tests (< 60 s): B0 on one rep; within-path tridiagonality; the B1 verdict (whichever way it goes) on two (space, κ) pairs;
  M̃ = ⊕ M^{(3)}_κ on two sectors; the triple preserves M̃; one depth-1 sector's hopping entries = the universal formula.
* `RESULTS-G8b-nested-gl3.md`: §0 headline (B1 verdict per operator with counts; closed forms if any, with held-out counts; the nested-basis structure of the hopping
  blocks in one sentence; the depth-1 statement), §1 conventions, §2–§5 per question with tables, §6 surprises and corrections (to this brief or to G4), §7 not done /
  limits, §8 usage. Every claim exact or "50 digits". If B1 fails, the headline must say what the within-path block IS (its closest description), not only that it is
  not the gl(3) chain.
