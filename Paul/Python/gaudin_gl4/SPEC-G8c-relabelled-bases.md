# Spec G8c: the 24 relabelled Gelfand–Tsetlin chains of gl(4) — bandwidth and coordinate cuts of the triple in each, and the smallest non-chain sectors by hand (2026-10-24)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram, no wolfbook tools, no git, no notebooks). Work package G8c of
`Paul/Mathematica/Gaudin/gl4-jacobi-plan.md` (read its §0 and §1 H-C first). Environment: Python 3.13, numpy, scipy, sympy, mpmath; run everything from
`Paul/Python/gaudin_gl4/` as `python -X utf8 ...`; at most two heavy processes at any time; `time.perf_counter()`. **Never edit existing files; import G1–G5 read-only.**
Files you own: `relabelled_gl4.py`, `test_relabelled_gl4.py`, `RESULTS-G8c-relabelled-bases.md`, `relabelled_gl4_output.txt`, `relabelled_gl4_tables.txt`. Two other
agents (G8a, G8b) work concurrently in the same directory on their own files — do not touch `lattice_recursion_gl4*`, `nested_gl3_gl4*`. Time-box ≈ 2.5 hours; deliver
what is finished and say precisely what is not.

Tooling trap (G1–G4 all hit it): Python source written through Bash heredocs loses `\n` escapes. Write code with the file tools and check `python -c "import ast; ast.parse(open('relabelled_gl4.py').read())"`.

## 0. Setting (binding conventions; see RESULTS-G4 §1, §5.1 and the docstring of `glN_reps.py`)

`R = gr.GLNRep.of(lam)` (`gr = glN_reps`): plain rational GT basis adapted to the chain gl(4) ⊃ gl(3)_{123} ⊃ gl(2)_{12} ⊃ gl(1)_{1}; pattern = tuple of 4 rows
(rows[0] = (m₁₁,), rows[1] = (m₁₂,m₂₂), rows[2] = (m₁₃,m₂₃,m₃₃), rows[3] = λ); `.pats .index .wsp[weight] .E[(i,j)] .apply`. The reduced operators are universal
polynomials in the site-1 generators: `red = so.reduced_op(name, lamp).red` is a dict {word: Fraction} with word = tuple of generators (i, j) (check the exact layout in
`slice_ops_gl4.ReducedOp` / `word_cols`; `so = slice_ops_gl4`, names 'V3', 'V4a', 'V4b', `so.TRIPLE`). `so.weight_space_mats(R, ops, W)` gives the exact w×w matrices
(M[row][col] = coefficient of e_row in V e_col) of ReducedOps on a weight space `W = R.wsp[nu]`, ν = μ − λ′. Strings: E_{i,i+1}^{k_i}, k_i = λ′_i − λ′_{i+1} + 1;
`se.prv_space`, `se.power_images(R, ij, k, W)`, `se._rref_nullspace` (`se = slice_engine_gl4`); copies `so.copies(lam, lamp, mu, dps=50)`; whole-space spectra and
graph statistics `structure_gl4.ident_sector(lam, lamp, mu)`, `structure_gl4.graph_stats(M)` (path test, RCM bandwidth), `structure_gl4.jspec`.

**Relabelled chains.** For σ ∈ S₄ (a permutation of {1,2,3,4}) the σ-chain is gl(4) ⊃ gl({σ₁,σ₂,σ₃}) ⊃ gl({σ₁,σ₂}) ⊃ gl({σ₁}), and the σ-GT basis of V_λ is the GT basis
adapted to it. Implementation without building new representations: the matrix of an operator X = polynomial in the E_ab in the σ-basis equals the matrix of the
relabelled polynomial φ_σ(X), E_ab ↦ E_{σ(a)σ(b)} (or σ⁻¹ — fix the convention by the check below), in the STANDARD basis of the same `GLNRep.of(lam)`, acting on the
weight space whose weight is the correspondingly permuted ν (the pattern's eigenvalue under E_{σ(i)σ(i)} must be ν_i). Reason: the relabelled representation
π∘φ_σ is isomorphic to π (same highest weight) through a Weyl group element w_σ ∈ GL(V_λ), and the standard GT basis of π∘φ_σ is the σ-GT basis of π. The gl(3)
"chain 2" of RESULTS-jacobi-cg §0 is σ = w₀ (E_ij → E_{4−i,4−j} in gl(3)). **Convention check (do first):** the joint spectrum of the relabelled triple on the relabelled
PRV space M^σ = ∩_i ker E_{σ(i)σ(i+1)}^{k_i} ∩ W^σ must equal Spec(M_μ) (50 digits, every sector of [3,2,1,0]²); the wrong convention (σ vs σ⁻¹, or permuting ν wrongly)
fails this. Note E_{σ(i)σ(i+1)} may be a LOWERING operator (σ(i) > σ(i+1)): that is fine — E_ab^k with a > b is still the string; implement `power_images` for any (a,b).

Products: P_small = [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², [3,2,1,0]⊗[4,2,1,0] (both orders), [4,2,2,0]², [4,3,1,0]²; P_big = [6,4,2,0]² (76 sectors,
w ≤ 15, m ≤ 12), [5,3,1,0]², [5,3,2,1]². Every statement exact (Fractions / exact rref) unless "50 digits" is said.

## 1. Questions

**C0. Implementation and convention check** (above). Also verify for σ = id that everything equals G3/G4, and for σ = w₀ = (4,3,2,1) that the E₃₄ string becomes a
coordinate cut (it maps to E₂₁ — the gl(2)₁₂ lowering operator — and E₂₁^k e_pattern = 0 iff m₁₁ − k < m₂₂ … derive and verify the cut inequality on data).

**C1. Bandwidth and chain structure per σ.** For every weight space (w ≥ 2) of P_small ∪ P_big and all 24 σ: the stencil of the triple in the σ-basis (the set of
displacements of the free coordinates (m₁₂, m₁₃, m₂₃) of the σ-patterns with a nonzero entry; is it still the 21-point set, or does the relabelled reduction change it?),
`graph_stats` of V₃ and of the union graph of the triple on W^σ: connected? path (tridiagonal in some order)? RCM bandwidth; the same on the coordinate subspace cut
out by the strings that are coordinate cuts in that σ (see C2), and on M^σ when it is a coordinate subspace. Tables: per σ the number of weight spaces where V₃ is a
path (G4: 29/44 for m = w and 19/106 for m < w at σ = id); per weight space the best σ (min bandwidth); the fraction of weight spaces where SOME σ gives a path, split by
m = w / m < w and by the active-string class A. Is there a rule A → σ (e.g. A = {3} → w₀, A = {2} → a chain with gl(2)_{23} innermost)? Test every candidate rule you can
state on the data and report hits/misses.

**C2. Coordinate cuts.** For every σ and every string i: is ker E_{σ(i)σ(i+1)}^{k_i} ∩ W^σ a coordinate subspace (spanned by σ-GT basis vectors)? Expectation: iff
{σ(i), σ(i+1)} = {1,2} (the innermost gl(2)); verify exactly on all sectors, and give the cut inequality in the σ-coordinates for the cases it is. For sectors with TWO
active strings, is there a σ in which both are coordinate cuts? (Expectation: never, as in gl(3) — but check; a counter-example would matter.) For sectors with one
active string, M^σ is a coordinate sub-lattice in the right σ: tabulate its shape (the inequalities) and whether the triple on it is a path.

**C3. The smallest non-chain sectors by hand.** (a) (5,4,2,1) of [3,2,1,0]² (m = w = 4, V₃-graph = K₄, G4 §5.1): list the 4 patterns and their coordinates, the three
exact 4×4 matrices, the joint eigenvalues (exact: factor the characteristic polynomial of V₃ over ℚ — G3 §8 says V₃ has a double eigenvalue 57 here; give the
eigenvectors of the triple exactly if the field is small, else 50 digits), the site-exchange parity P on this sector (λ = λ′: how P acts on the 4 patterns / on the
copies), the matrices in all 24 σ-bases with their sparsity graphs, and the best σ. Then the change of basis to the within-path eigenbases (G8b's nested basis: the
eigenvectors of the tridiagonal κ-blocks) — does the K₄ become a path there? (b) the same for the next three smallest sectors (by w) whose V₃-graph is not a path
(from G4's `graph_stats` over P_small; name them), and for one m < w sector with a cycle. (c) For these sectors, the simplest description you can find of the physical
copies (e.g. components as rationals or as roots of a quadratic; a pattern in the zeros).

**C4. (If time remains.) The exchange parity as a lattice symmetry.** For λ = λ′ sectors, G2/G3: the site exchange P satisfies PV₃P = −V₃ + c and the copies come in
P-pairs / P-eigenlines (signed permutations in the canonical basis, G3 §8). Is P, as a map on the slice, visible on the GT lattice in some σ-basis (a lattice reflection,
e.g. κ ↦ its complement in λ)? Test on [3,2,1,0]² and [4,2,1,0]²: compute P on M_μ from the canonical copies (G3 §8's signed permutation) and compare with the candidate
lattice involutions.

## 2. Deliverables and acceptance

* `relabelled_gl4.py` with a CLI `python -X utf8 relabelled_gl4.py check | bandwidth | cuts | small | parity | all` (and `run <stage,...>`), tables to
  `relabelled_gl4_tables.txt`, printed output of the final run to `relabelled_gl4_output.txt`; a public `sigma_mats(lam, lamp, mu, sigma)` returning the σ-patterns, the
  σ-coordinates and the exact matrices of the triple and of the strings on W^σ.
* `test_relabelled_gl4.py`: ≥ 6 exact regression tests (< 60 s): σ = id equals G3; the convention check (joint spectrum of M^σ = Spec(M_μ)) for two σ on one sector;
  the w₀ cut inequality; one coordinate-cut classification; the K₄ example's exact matrices; one bandwidth statistic.
* `RESULTS-G8c-relabelled-bases.md`: §0 headline (the per-σ path counts in one table; the cut rule; the verdict on "some σ makes the triple a path" with fractions
  by class; the K₄ anatomy in five lines; the parity result if done), §1 conventions (state the σ convention you fixed and how it was checked), §2–§6 per question with
  tables, §7 surprises and corrections (to this brief or to G4), §8 not done / limits, §9 usage. Every claim exact or "50 digits".
