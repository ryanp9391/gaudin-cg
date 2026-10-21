# Plan: the slice engine and the Gaudin missing-label operators for gl(4) (written 2026-10-21 night, for the next session; NOT started)

Paul: "the vertex operator approach may indeed extend to gl(4). I don't want to do it now, but set up a plan for the guy tomorrow to
implement." Context: for gl(3) the highest-weight (HW) Clebsch–Gordan problem of V_λ⊗V_λ′ collapses to a w-dimensional slice, the
multiplicity is resolved by ONE effective operator (the two-site Gaudin vertex operator V reduced to the slice), which is a Jacobi chain with
closed-form entries; see `Slice-engine-notes.{tex,md,pdf}`, `Paul/Python/gaudin_sl3/{slice_engine.py,jacobi_cg.py}` and their RESULTS, and the
literature check `Literature-check-missing-label.md`. For SU(3) a specialised recursion (Draayer–Akiyama; ndsu3lib, Herko et al. EPJA 62
(2026) 95) already exists and is of comparable speed; it is SU(3)-specific (one missing label, Biedenharn–Louck–Hecht zero pattern, closed-form
SU(3) seeds). For gl(N), N ≥ 4, the only general tool is the Munich null-space code (`Paul/Cpp/AlexCGC`, arXiv:1009.0437), which costs
O(n_W³) and dies at moderate size. The slice theorem is general, and the two-site Gaudin algebra supplies the right number of commuting
missing-label operators for every N — so gl(4) is where the method should give orders of magnitude AND a canonical multiplicity basis that
nobody has. This revives the project's original gl(4)/su(4) interest (CLAUDE.md and the memory still say gl(4) was dropped on 2026-10-16;
update both when this starts).

## 0. What carries over unchanged, what changes

| | gl(3) (done) | gl(4) (to do) |
|---|---|---|
| slice theorem | M_μ = {v ∈ V_λ[μ−λ′] : E₁₂^{λ′₁−λ′₂+1}v = E₂₃^{λ′₂−λ′₃+1}v = 0} (Kostant 1959 / PRV 1967) | same with THREE string conditions E_{i,i+1}^{λ′_i−λ′_{i+1}+1} v = 0, i = 1,2,3; dim M_μ = multiplicity; general theorem, nothing to prove |
| level recursion for the full vector | left inverses of the site-2 raising blocks, once per λ′ | same, with the gl(4) GT raising formulas |
| weight space V_λ[μ−λ′] | a chain in one GT entry (m₁₂), w ≤ a+1 | a 3-parameter region of the GT lattice (6 entries, 3 weight constraints); w grows faster (tens to hundreds for dim ~10⁴) but still ≪ n_W |
| missing labels needed | 1 | (N−1)(N−2)/2 = 3 |
| the operators | the (1,1,2) cubic invariant V = ½Σ sgn sgn e¹e¹e² | the mixed invariants of the two-site gl(4) Gaudin model: one cubic ((1,1,2) ≅ (1,2,2) on fixed μ) and two quartic ((1,1,1,2) ≅ (1,2,2,2), and (1,1,2,2)); the exact Gaudin combinations (quantum corrections!) are the ones that commute |
| reduction to the slice | one site-2 generator per term → one-line substitution (Cartan → λ′_c, lowering → 0, raising → −E¹) | terms with up to two site-2 generators: normal-order the site-2 factors (lowering left, Cartan, raising right) with the gl(4) commutation relations, then the same substitution, rightmost raising first (Zhelobenko extremal-projector calculus) |
| structure | tridiagonal, closed forms, three-term recurrence | banded in the GT lattice; joint eigenvectors of three commuting operators; closed forms and any recurrence structure are open |
| integrability cross-check | V = a₁ − η/2 (accessory parameter of the oper); chain char poly = Q-system polynomial | three accessory parameters of a 4th-order Fuchsian oper; the repo already has the two-site gl(4) Gaudin/Baxter solver (see §4) |
| benchmark | Munich 206 s → 0.006 s at [14,7,0]²; SU(3)-specific codes comparable | Munich's own SU(4) sweep exists: `Paul/Cpp/AlexCGC/results/hw_timings.txt`, [S,S,0,0]², S = 1..8 (S = 8: dim 7425, 80 s total) — multiplicity-free, so the slice alone (no operators) already settles it; multiplicities need another family |

## 1. Assets already in the repo (read before writing anything)

- gl(3) code to generalise: `Paul/Python/gl3_reps.py` (GT patterns, `_raise/_lower` formulas, `GL3Rep`), `Paul/Python/gaudin_sl3/slice_engine.py`
  (prv_space, level recursion with precomputed left inverses, V_red, Racah–Speiser), `jacobi_cg.py` (the chain), `gauss_dictionary.py` (PRV
  utilities), their SPEC/RESULTS files (conventions, pitfalls: GT-order traps, float precision of the recursion at large a, point choice).
- gl(4) material from the old programme (OLD conventions — port deliberately, check every sign): `Paul/Mathematica/old-conventions/old/Gaudin/Clean/Hamiltonians_SU4.wb`
  (gl(4)/su(4) GT basis and generators, Gaudin Hamiltonians), `FastQ.wb` (residue formulas), `Experiments/ResidueBuilder.wb`; the gl(4) L = 2
  Baxter/Gaudin work of July 2026 (HANDOFF.md blocks "State as of 2026-07-22": `baxter_SU4_L2_XXX.wb`, standalone `GaudinSolver.wb`, the
  untwisted τ₁, τ₂, τ₃ closed forms anchored at ΔCas = Cas(μ) − Cas(μ_top)); `Paul/Python/SU4_Paul_Normalisation_EXPLICIT.wl`; the gl(4)
  untwisting notes in memory `project-gl4-rect-untwist-why-clean`. Conventions of the gl(3) Gaudin oper: `Gaudin-sl3-CG-notes.md` §1–2.
- The Munich code and its SU(4) numbers: `Paul/Cpp/AlexCGC/README.md`, `results/hw_timings.txt`, `results/summary.csv`, `hwcgc.sh` (extracts the
  HW CG rows for any N — the reference data generator for small SU(4) cases; the binary is built on the Ubuntu laptop, not on Windows).

## 2. Work packages (Sonnet-sized, in dependency order; Python unless stated; new directory `Paul/Python/gaudin_gl4/` with SPEC-*.md / RESULTS-*.md beside each script; `*.pkl`/`*.npz` gitignored)

### WP-G1 — `glN_reps.py` and the gl(4) slice engine (pure representation theory; the SU(4) headline)
- GT basis and generators for gl(N) (the standard Gelfand–Tsetlin formulas for E_{k,k+1}, E_{k+1,k}, diagonal weights; others by commutators),
  sparse. Tests: N = 3 reproduces `gl3_reps.py` exactly (same pattern order or a documented map); N = 4: all [E_ij, E_kl] = δ_jk E_il − δ_li E_kj
  on random vectors for several irreps, Casimir eigenvalues, dimensions (Weyl formula), weight multiplicities (Kostant/Freudenthal or by
  counting patterns).
- Multiplicities of V_λ⊗V_λ′ for gl(4): Racah–Speiser with GT weight multiplicities (general N; cross-check Σ m dim = dim·dim′ and, for small
  cases, Littlewood–Richardson).
- `prv_space(lam, lamp, mu)` with the three string conditions; assert dim = multiplicity on every irrep of: [1,0,0,0]², [1,1,0,0]², [2,1,0,0]²,
  [2,1,1,0]² (adjoint-like), [2,2,0,0]², [3,2,1,0]², and λ ≠ λ′ cases ([2,1,0,0]⊗[1,1,0,0], [3,1,0,0]⊗[2,2,0,0]); then the rectangular family
  [S,S,0,0]², S = 1..10 (multiplicity-free: every M_μ is 1-dimensional — the slice IS the HW vector's germ).
- Level recursion `hw_from_slice` for gl(4) (three raising blocks; precomputed left inverses per λ′); HW residuals ≤ 1e-12, weights exact.
- Reference check against Munich: on the Ubuntu laptop (or wherever `clebsch`/`cgc_hw_bench` is built) run `hwcgc.sh 4 "S S 0 0" "S S 0 0" "mu"`
  for all μ of [S,S,0,0]², S ≤ 3, and compare the HW vectors up to normalisation (Munich uses the orthonormal GT basis; ours the plain GT
  basis — convert with the GT metric). 0 mismatches is the acceptance.
- Benchmark: all HW vectors of [S,S,0,0]² for S = 1..10 (and as far as feasible), t_a (slices) / t_b (full vectors), against Munich's table
  (S = 8: 80 s; S ≥ 9 not available there). Table of w_max, n_W,max, #irreps.
- Deliverable: `glN_reps.py`, `slice_engine_gl4.py`, RESULTS with the table; `ref/gl4-cg-slices-<lam>.txt` in the gl(3) file format.

### WP-G2 — the two-site gl(4) Gaudin missing-label operators (algebra; Python exact, or Mathematica K1 if the old notebooks are easier to reuse)
- Build the untwisted two-site gl(4) Lax matrix L(u) = E⁽¹⁾/(u−θ₁) + E⁽²⁾/(u−θ₂) and its quantum characteristic polynomial
  cdet(∂ − L(u)) = ∂⁴ − S₁∂³ + S₂∂² − S₃∂ + S₄ in the SAME conventions as the gl(3) notes §1 (column-ordered cdet) — the old gl(4) notebooks
  have this; port it and re-verify commutativity numerically on a small product (e.g. [1,0,0,0]⊗[1,1,0,0] and [2,1,0,0]²): all Laurent
  coefficients of S₂, S₃, S₄ at θ₁ commute with each other and with the diagonal gl(4).
- Identify the independent mixed invariants on a fixed irrep μ: remove site Casimirs and total Casimirs; expect exactly three independent
  operators: one cubic (the (1,1,2) one — check it equals the gl(3)-type ½Σ sgn sgn e¹e¹e² form generalised, up to Casimir terms) and two quartic.
  Write them as explicit polynomials in E⁽¹⁾_{ij}, E⁽²⁾_{kl} with the exact quantum ordering/corrections (store as term lists).
- Verify on test products with multiplicities ≥ 2 (e.g. [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [2,2,1,0]⊗[2,1,1,0]): the three operators commute
  on the HW space of each μ and have a simple JOINT spectrum (and record whether the cubic one alone already has simple spectrum — for gl(3) it
  does; if it does for gl(4) too, one operator suffices in practice, but keep the triple as the canonical definition). If the joint spectrum is
  degenerate anywhere, report the sector and the degeneracy; the fallback is the full shift-of-argument algebra (twist-direction dependent) —
  flag and stop rather than invent.
- Deliverable: `gaudin_gl4_ops.py` (term lists + builders + commutativity tests), RESULTS with the operator count, the simple-spectrum table.

### WP-G3 — reduction to the slice and the canonical copies
- Implement the general reduction: for a term with site-2 factors e²_{a₁b₁}…e²_{a_kb_k}, normal-order on site 2 (lowering operators to the left,
  then Cartan, then raising) using the gl(4) relations; acting on the top state of site 2 from the left, lowering → 0, Cartan → λ′_c, each
  raising e²_{ab} → −e¹_{ab} (transfer to site 1 with sign, rightmost first; the site-1 factors already present commute with site-2 operators).
  Result: V_red,k = polynomial in site-1 generators acting on V_λ[μ−λ′]. Verify against the direct route (reconstruct the m basis HW vectors with
  WP-G1's recursion, apply the full two-site operator, read the slice) on all sectors of the test products: identical to 1e-12; check the
  V_red,k preserve M_μ and commute on it.
- Joint diagonalisation on M_μ (m ≤ a few tens) → the canonical copies' slices; the full vectors by the recursion; G-orthogonality of copies
  (GT metric) as a check. Timing: n_W-free part vs full reconstruction, for the test products and for a family with growing multiplicity
  (candidates: [2a,a,0,0]² or [3a,2a,a,0]² — compute their multiplicity tables first and pick the one with m growing like a).
- Deliverable: `slice_ops_gl4.py`, RESULTS; `ref/gl4-cg-slices-*.txt` extended with the copies.

### WP-G4 — structure (the gl(4) analogue of the Jacobi chain; research, after G3)
- The GT lattice of V_λ[μ−λ′]: which entries move under each V_red,k (bandwidth per operator); matrix elements as rational functions of the GT
  labels (guess from the data as J did for gl(3): products of linear factors; verify on a held-out rep); whether M_μ is cut out by coordinate
  conditions (sub-lattice) or is a spectral subspace; whether the cubic operator alone is "tridiagonal along one direction" so that eigenvectors
  are products/sums of three-term recurrences (multivariate orthogonal-polynomial structure?). Deliverable: RESULTS with tables and whatever
  closed forms survive a held-out test.

### WP-G5 — the integrability cross-check (ties to the old gl(4) work)
- The joint eigenvalues of (V₃, V₄ᵃ, V₄ᵇ) on each HW sector vs the untwisted two-site gl(4) Baxter/Q-system solver of July 2026
  (`GaudinSolver.wb`; HANDOFF 2026-07-22: τ₁, τ₂, τ₃ closed forms on HW sectors): the dictionary eigenvalue ↔ accessory parameters of the
  4th-order oper (the gl(3) one was V = a₁ − η/2 with η = |Λ|² + cas₂(Λ) − cas₂(μ)/2). Report the map and the agreement (30 digits).

### WP-G6 — benchmark table and a section for the notes
- Munich (general N) vs slice engine vs operators for SU(4): [S,S,0,0]² (multiplicity-free: pure slice) and the multiplicity family; same
  machine where possible. One honest paragraph on what is new for N = 4 (canonical basis + closed-form-free but n_W-free algorithm) and what
  is classical (the slice theorem; Gaudin algebras as missing labels: Vinberg, Rybnikov, FFR, MTV, HKRW — see the literature check).

## 3. Acceptance for the whole programme
dim M_μ = multiplicity on every irrep of every test product (λ = λ′ and λ ≠ λ′); HW vectors agree with Munich on S ≤ 3 rectangular cases;
three commuting effective operators with simple joint spectrum on all multiplicity spaces of ≥ 3 test products; slice-reduced operators =
direct operators to 1e-12; a benchmark table with Munich's S = 8 (80 s) vs ours, and a multiplicity family to m ≥ 5.

## 4. Risks
- The two-site gl(4) Gaudin algebra might not separate every multiplicity space (MTV simple-spectrum results are for generic evaluation points;
  two points are affinely generic, but check); if it fails somewhere, the twisted (shift-of-argument) operators are the fallback — document.
- Quantum corrections in cdet: the naive trace-form invariants do NOT commute; use the exact Gaudin combinations from the old notebooks and
  verify commutativity numerically before anything else.
- GT generators for gl(4) in Python: port with tests (commutation relations, Casimirs); the old notebooks' normalisation of the GT basis may
  differ from the Munich code's orthonormal one — compare vectors only after converting with the GT metric.
- w grows: the exact PRV null space is still cheap (w ~ 10²), but the full reconstruction is output-bound; at dim ~10⁴ expect seconds to
  minutes for all vectors — report both numbers as for gl(3).
- Conventions: keep NEW conventions (CLAUDE.md), and keep this gl(4) work in its own directory; do not mix with the gl(3) files.
