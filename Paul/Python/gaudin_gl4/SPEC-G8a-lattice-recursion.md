# Spec G8a: block-Jacobi structure of the gl(4) triple with respect to lattice gradings — the forward recursion that builds a slice from its eigenvalue triple (2026-10-24)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram, no wolfbook tools, no git, no notebooks). Work package G8a of
`Paul/Mathematica/Gaudin/gl4-jacobi-plan.md` (read its §0 and §1 H-A first). Environment: Python 3.13, numpy, scipy, sympy, mpmath; run everything from
`Paul/Python/gaudin_gl4/` as `python -X utf8 ...`; at most two heavy processes at any time; `time.perf_counter()`. **Never edit existing files; import G1–G5 read-only.**
Files you own: `lattice_recursion_gl4.py`, `test_lattice_recursion_gl4.py`, `RESULTS-G8a-lattice-recursion.md`, `lattice_recursion_gl4_output.txt`,
`lattice_recursion_gl4_tables.txt`. Two other agents (G8b, G8c) work concurrently in the same directory on their own files — do not touch `nested_gl3_gl4*`,
`relabelled_gl4*`. Time-box ≈ 2.5 hours; deliver what is finished and say precisely what is not.

Tooling trap (G1–G4 all hit it): Python source written through Bash heredocs loses `\n` escapes. Write code with the file tools and check `python -c "import ast; ast.parse(open('lattice_recursion_gl4.py').read())"`.

## 0. Setting (binding conventions; see RESULTS-G4 §1 and the docstring of `glN_reps.py`)

Site 1 = V_λ in the plain rational GT basis of `glN_reps.GLNRep.of(lam)`: pattern = tuple of 4 rows, rows[0] = (m₁₁,), rows[1] = (m₁₂, m₂₂), rows[2] = (m₁₃, m₂₃, m₃₃),
rows[3] = λ; `R.pats[i]`, `R.index[pattern]`, `R.wsp[weight]` = numpy array of global indices of the weight space. A sector (λ, λ′, μ) has the weight space
W = V_λ[ν], ν = μ − λ′, `W = R.wsp[nu]`; the row sums are fixed by ν (m₁₁ = ν₁, m₁₂ + m₂₂ = ν₁ + ν₂, m₁₃ + m₂₃ + m₃₃ = ν₁ + ν₂ + ν₃), free coordinates
x = (m₁₂, m₁₃, m₂₃) (`structure_gl4.coords(pat)` if it does that — check — else write your own). The three reduced operators on W:
`ops = {n: so.reduced_op(n, lamp) for n in so.TRIPLE}` (`so = slice_ops_gl4`, TRIPLE = ('V3','V4a','V4b')), `mats = so.weight_space_mats(R, ops, W)` →
exact w×w lists of Fractions with **M[row][col] = coefficient of e_row in V e_col** (columns are images). The E₁₂ cut: the string E₁₂^{k₁}, k₁ = λ′₁ − λ′₂ + 1,
kills a GT vector iff m₁₂ ≤ μ₁ − λ′₂ (E₁₂ acts inside the gl(2)₁₂ irrep of spin (m₁₂ − m₂₂)/2 with m₁₁ fixed) — so W₁ := span{x ∈ W : m₁₂ ≤ μ₁ − λ′₂} is a coordinate
subspace, invariant under the triple (BGG intertwining, G3 §7a); verify both on data (`se.power_images(R, (1,2), k1, W)`; `se = slice_engine_gl4`). M_μ ⊂ W₁
(`se.prv_space(lam, lamp, mu)`: `.W .basis .free .m .w`), the canonical copies `so.copies(lam, lamp, mu, dps=50)` (each: `ev` = (v₃,v₄ᵃ,v₄ᵇ) mpf, `slice` =
{global GT index: mpf}, first nonzero component = 1), and the whole-weight-space joint spectrum `structure_gl4.ident_sector(lam, lamp, mu)` (`js['trip']`,
`js['U']` columns = eigenvectors in the GT basis, `phys` flags, `S` violated strings, `active` strings). M̃ := ker E₁₂^{k₁} ∩ ker E₂₃^{k₂} ⊂ W (build it like
`se.prv_space` with i ∈ {1,2} only, using `se.power_images` and `se._rref_nullspace`) — G8b studies it; here it is just a third invariant subspace to grade.

Products: P_small = [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², [3,2,1,0]⊗[4,2,1,0] and [4,2,1,0]⊗[3,2,1,0], [4,2,2,0]², [4,3,1,0]²;
P_big = [6,4,2,0]² (all 76 sectors, w ≤ 15, m ≤ 12), [5,3,1,0]², [5,3,2,1]². Sectors = all dominant μ with non-empty W (`gr.candidates(lam, lamp)` or
`structure_gl4.sector_list(lam, lamp, 'all')`); keep the m = 0 ones for the grading statistics but mark them.

## 1. Questions (answer each with data; exact unless "50 digits" is said)

**A1. Gradings.** For g = a·m₁₂ + b·m₁₃ + c·m₂₃ with (a,b,c) ∈ {−2,…,2}³ \ {0}, one of each ±pair (62 gradings; m₃₃-gradings are included as −m₁₃−m₂₃), and for every
weight space of P_small ∪ P_big with w ≥ 2, on each of the three spaces S ∈ {W, W₁, M̃} (restrict the matrices to the coordinate subspace W₁; for M̃, which is not a
coordinate subspace, grade W and treat M̃ only in A3):
  * band β = max |g(target) − g(source)| over the nonzero entries of the three matrices (structural bound: |a| + max(|b|,|c|,|b−c|); report the data value);
  * the levels L_d (distinct g-values ascending), their sizes n_d, the bottom size n₀ and the top size;
  * **forward determinacy** for band-1 gradings: at every level d, rank of the stacked (3n_d × n_{d+1}) matrix F_d = [M_k[L_d, L_{d+1}]]_{k ∈ TRIPLE} equals n_{d+1}
    (exact rref over ℚ); for band-2 gradings the greedy version: process levels ascending, at level d use its 3n_d equations to determine all not-yet-determined
    unknowns in levels ≤ d + 2 (full column rank of the stacked block on those unknowns); record at which level (if any) it fails and by how much (rank deficit);
  * the same in the descending direction (g → −g is in the list, so this is automatic).
  Tables: per space S and grading: number of weight spaces where the recursion is fully determined, the distribution of n₀ (how many have n₀ = 1), the worst rank
  deficit. Then per weight space: the best grading(s) (fully determined, smallest n₀), and the relation to the active strings A (from `ident_sector`): is there a rule
  "A → grading"? Specifically test the natural candidates: +m₁₂ and −m₁₂ on W₁ (bottom = far from the cut / at the cut), ±m₁₃, ±m₂₃, ±m₃₃, and the diagonal ones
  ±(m₁₂ + m₁₃), ±(m₁₂ − m₁₃), ±(m₁₃ + m₂₃ + m₁₂) etc.
  Side check (the plan's Vandermonde fact): for 50 random (source x, sector) pairs where all three Δm₁₂ = +1 root-move targets exist, the 3×3 matrix of
  (V₃, V₄ᵃ, V₄ᵇ)-entries to those targets has determinant = (Π_α c_α)·(y₂−y₁)(y₂−y₃) with y_k = m_{k3} − k and c_α the V₃ entries (exact) — confirm or correct.

**A2. The recursion as a formula.** Implement `slice_from_triple(lam, lamp, mu, trip, grading, space='W1', dps=50)`: given the eigenvalue triple of a copy (from
`copies`), run the forward recursion in mpmath: ψ on the bottom level = the initial data (if n₀ = 1: the number 1; if n₀ > 1: see A4), then level by level solve the
stacked least-squares/exact system F_d ψ_{d+1} = [(v_k − M_k[L_d,L_d])ψ_d − M_k[L_d,L_{d−1}]ψ_{d−1}]_k (use the exact left inverse of F_d computed once per sector over ℚ
as `se.exact_left_inverse` does, or mpmath QR), and return the slice in G3's normalisation (first nonzero component in ascending GT index = 1). Validate: every copy of
every sector of P_small and of [6,4,2,0]² with a fully determined grading and n₀ = 1 is reproduced to ≤ 1e-40 (50 digits). Report the residual of the unused equations
(should vanish at a true eigenvalue triple: these are the characteristic equations) and its size at a perturbed triple. Timing: the recursion per copy vs the eigen-decomposition.

**A3. Polynomials and the characteristic equations.** For the sectors with n₀ = 1 on W₁ (or W) and a fully determined band-1 grading, run the recursion symbolically
(sympy, exact rationals, symbols v3, v4a, v4b): ψ_x = P_x(v₃,v₄ᵃ,v₄ᵇ) for every lattice point x. Report: the total degree of P_x vs the level of x and vs G4's Krylov
degrees (`structure_gl4.left_krylov`; G4 found the generic counting degree d_gen(m)); the unused equations C_j(v) (the rows of each level's system beyond the rank, plus
the top level's 3 n_top equations) as polynomials; check that their common zeros on the joint spectrum: evaluate C_j at every joint triple of W₁ (physical and spurious,
`ident_sector`) — all vanish; and at the physical triples only; compare the ideal they generate (sympy `groebner`, lex or grevlex) with the vanishing ideal of the joint
spectrum of W₁ (dimension of ℚ[v]/J should be dim W₁ if they generate it; report when it is larger). Do this for ≥ 10 sectors across P_small incl. m ≥ 3; for [6,4,2,0]²
pick 3 sectors with m ≥ 6. For M̃ ⊂ W₁ the physical triples are a subset: say which of the C_j are needed. This is the gl(4) analogue of "continuant = Q-system polynomial"
(RESULTS-jacobi-cg §4); state it in those terms.

**A4. When n₀ > 1.** For band-1 gradings with full determinacy but n₀ > 1 (the common case, presumably): the recursion maps the initial data ψ₀ ∈ ℚ^{n₀} linearly to the
whole vector, ψ = T(v)ψ₀, T polynomial in v; the characteristic equations become C(v)ψ₀ = 0 with C a polynomial matrix; at a joint eigenvalue rank C(v) = n₀ − 1 and
ker C(v) ∝ the copy's bottom-level component. Verify this on data (50 digits) for ≥ 20 sectors with n₀ ∈ {2, 3}, and report the minors of C(v) that vanish (the n₀ × n₀
minors generate a multiple of the ideal). Is there a grading for which n₀ = 1 in EVERY sector of a product (for a natural choice of space S)? If not, what is the minimum
of n₀ over gradings per sector, as a function of (w, m, A)? (A table.)

**A5. Separability diagnostic (cheap, 15 min).** For every physical copy of P_small: the matrix Ψ[κ][m₁₂] of its slice (rows κ = row 3, columns m₁₂; entries missing where
the pattern does not exist), its rank over the κ's with the full m₁₂-range and over all (missing → 0), 50 digits with a 1e-30 threshold. Rank 1 would mean ψ(m₁₂, κ) =
f(m₁₂)g(κ). Report the distribution; also the same for the matrix Ψ[m₁₃][(m₂₃, m₁₂)] and Ψ[m₂₃][(m₁₃, m₁₂)].

## 2. Deliverables and acceptance

* `lattice_recursion_gl4.py` with a CLI `python -X utf8 lattice_recursion_gl4.py gradings | recursion | charpoly | n0 | separability | all` (and `run <stage,...>`), writing
  the per-stage tables to `lattice_recursion_gl4_tables.txt` and the printed output of the final run to `lattice_recursion_gl4_output.txt`; `slice_from_triple` as a
  public function with a docstring giving the recursion in formulas.
* `test_lattice_recursion_gl4.py`: ≥ 6 exact regression tests (< 60 s): the E₁₂ cut is ker E₁₂^{k₁} on two sectors; W₁ invariant; the Vandermonde determinant on one
  source; forward determinacy on two sectors with the best grading; `slice_from_triple` reproduces two copies to 1e-40; one characteristic-equation check.
* `RESULTS-G8a-lattice-recursion.md`: §0 headline (numbers: for each space S the best grading and the fraction of weight spaces where it is fully determined, n₀ = 1
  fraction; the recursion validation count and precision; the characteristic-equation statement), §1 conventions, §2–§6 per question with tables, §7 surprises and
  corrections (to this brief or to G4), §8 not done / limits, §9 usage. Every claim exact or "50 digits".
* A HANDOFF block is NOT yours; the session will write it from your RESULTS §0.
