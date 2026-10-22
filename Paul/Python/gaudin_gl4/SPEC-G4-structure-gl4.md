# Spec G4: structure of the reduced gl(4) missing-label operators on the GT lattice — stencils, closed-form matrix elements, the Weyl/BGG characterisation of the spectrum, chain structure (2026-10-22; final, written after RESULTS-G3)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram, no wolfbook tools, no git, no notebooks). Work package WP-G4 of
`Paul/Mathematica/Gaudin/gl4-slice-engine-plan.md` (research: the gl(4) analogue of the gl(3) Jacobi chain). Environment as in SPEC-G1/G2/G3 (`python -X utf8`, run from
`Paul/Python/gaudin_gl4/`, ≤ 2 heavy processes, `perf_counter`). **Never edit existing files; import G1/G2/G3/G5 read-only.** Files you own: `structure_gl4.py`,
`test_structure_gl4.py`, `RESULTS-G4-structure-gl4.md`, `structure_gl4_output.txt`, `structure_gl4_tables.txt`. Another agent (G6) is concurrently writing the notes
`Paul/Mathematica/Gaudin/gl4-slice-engine-notes.*` and `RESULTS-G6-*` — do not touch those. Time-box ≈ 2.5 hours; deliver what is finished and say what is not.

Read first: this file; the plan §0 (row "structure") and WP-G4; `RESULTS-G3-slice-ops-gl4.md` §3 (the universal reduced polynomials and their sparsity), §6 (scalars c, c′), §7
(three exact structural facts: BGG intertwining, commutation on whole weight spaces for every λ′, M_μ is a spectral subspace), §8 (copies, exchange structure), §13 (usage);
`RESULTS-G5-oper-gl4.md` §2 (the dictionary, T_k, τ_k) and §6; `RESULTS-G1-slice-engine-gl4.md` §0, §12; the gl(3) result you generalise: `Paul/Python/gaudin_sl3/jacobi_cg.py`
(docstring; `derive_symbolic`, `check_index_structure`, `chain`, `mult_closed`, `slice_from_V`) and `RESULTS-jacobi-cg.md` §7 (T1–T6: the weight space is a chain in m₁₂,
V_red tridiagonal with closed forms, the E₁₂ condition is a coordinate cut, E₂₃ the same cut in the Weyl-reflected chain, M_μ = sub-chain or spectral subspace of the m
LARGEST eigenvalues, slice = p_k(V)).

APIs: G3 — `so.reduce_terms(so.op_terms('V3'), lamp)` → {site-1 word: Fraction} (V₃: 39 universal words, V₄ᵃ: ~61, V₄ᵇ: 85), `so.reduce_terms_poly(terms, 4)` (the universal
polynomial in λ′, symbolic), `so.V_red(lam, lamp, mu, name)` (exact m×m in G1's PRV basis; names 'V3','V4a','V4b','S31','S41','S_k_j','S_k_j_sw'), `so.triple_mats`,
`so.copies(lam, lamp, mu, dps)` (canonical copies: `ev`, `slice` {global GT index: mpf}, `coords`), `so.normal_order`; G3 also has the whole-weight-space operators
(§7: the operators on W = V_λ[μ−λ′] in the GT basis — look for the function used in §7 (b)/(c), e.g. a `vred_matrix_on_weight_space`-type helper; if there is none, build the
w×w matrix from `reduce_terms` + `GLNRep.apply` yourself). G1 — `gr.GLNRep.of(lam)` (`.pats`, `.wsp`, `.E`, `.apply`), `se.prv_space(lam, lamp, mu)` (`.W .basis .free .m .w`),
`gr.racah_speiser`, `gr.candidates`. G5 — `oper_gl4.Oper`, `solve_sector` (not needed unless you want the oper side).

## 0. What is known (G3 §7, exact on the tested products) and the questions

Setting: site 1 = V_λ in the GT basis (pattern m_{ij}, rows (m₁₁), (m₁₂,m₂₂), (m₁₃,m₂₃,m₃₃), λ); a weight space W = V_λ[ν], ν = μ − λ′, fixes the row sums and leaves
THREE free lattice coordinates (e.g. m₁₂, m₁₃, m₂₃): a 3-parameter region R_ν of the GT lattice, w = |R_ν|. The reduced operators V_red,k(λ′), k ∈ {3, 4a, 4b}, are
weight-zero polynomials in the site-1 generators, universal and polynomial in λ′; M_μ = ∩_i ker(E_{i,i+1}^{k_i}) ⊂ W, k_i = λ′_i − λ′_{i+1} + 1 (PRV).
Known: (a) **BGG intertwining** E_i^{k_i} V_red(λ′) = V_red(s_i·λ′) E_i^{k_i} as maps V_λ[μ−λ′] → V_λ[μ−s_i·λ′], s_i·λ′ = λ′ − k_iα_i (dot action); (b) the three V_red,k(λ′)
**commute on the whole W for every λ′** (also non-dominant rational λ′); (c) **the joint spectrum on W is simple and M_μ is spanned by m of its w joint eigenlines** (the other
w − m eigenvectors are "spurious": they violate a string condition). Interpretation (G3, not proved): V_red(λ′) are the Gaudin operators on the singular vectors of
V_λ ⊗ M(λ′), M(λ′) the Verma module, and (a) is compatibility with the embeddings M(s_i·λ′) ↪ M(λ′).

Questions (answer each with data; closed forms only if they survive a held-out test; proofs are welcome but not required — say clearly what is verified and what is argued):

1. **Stencils.** For each operator the set of lattice displacements Δ(m₁₂, m₁₃, m₂₃) (equivalently the displacement of all six entries) occurring with nonzero matrix element,
   from the universal reduced words (which generators move which entries: E_{i,i+1} raises one entry of row i, E_{i+1,i} lowers one, longer generators move entries in rows
   i..j−1) and confirmed on data; the structural argument as in `check_index_structure`. Table of stencil sizes per operator; the gl(3) analogue was {0, ±(1,−1)} for V.
2. **Closed-form matrix elements.** Derive the matrix elements of V_red,3 (then V_red,4a, V_red,4b if feasible) as rational functions of the six GT labels, λ′ (symbolic from
   `reduce_terms_poly`) and λ, with sympy through the symbolic GT formulas (as `derive_symbolic` did; the off-diagonal elements should factor into products of linear factors),
   and **verify on held-out representations** (≥ 3 reps not used in the derivation, all weight spaces, exact against G3's matrices). Report the forms that survive and the ones
   that do not (size of the expressions, factorisation).
3. **The Weyl/BGG characterisation of the physical spectrum (the main question).** By (a), a joint eigenvector v of V_red(λ′) on V_λ[μ−λ′] with E_i^{k_i}v ≠ 0 is mapped to a
   joint eigenvector of V_red(s_i·λ′) on V_λ[μ−s_i·λ′] with the SAME joint eigenvalue. Conjecture (the Gaudin-spectrum version of the BGG resolution / Weyl character formula):
   **for every joint triple τ, Σ_{w∈S₄} sgn(w)·mult_τ(V_red(w·λ′) on V_λ[μ − w·λ′]) = [τ ∈ Spec(M_μ)]**, where w·λ′ = w(λ′+ρ) − ρ (dot action; V_λ[μ−w·λ′] = 0 when μ − w·λ′ is not
   a weight of V_λ), mult_τ = multiplicity of τ in the joint spectrum of the commuting triple at the (generally non-dominant) weight w·λ′, and [·] ∈ {0,1}. Test it exactly/at 50
   digits on every sector of [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]² (and [6,4,2,0]² if time permits): compute the joint spectra of V_red(w·λ′) on all 24 Weyl
   translates, match triples across spaces (tolerance 1e-30), and report (i) whether the signed count is exactly the indicator of Spec(M_μ), (ii) the simpler one-step version
   (spurious triples of W = ∪_i triples inherited from s_i·λ′, i = 1,2,3 — is the union disjoint, or do deeper translates matter?), (iii) the gl(3) sanity check of the same
   statement with `slice_engine`/`jacobi_cg` data on [4,2,0]² (there the one-step version is T2/T6 of RESULTS-jacobi-cg: the spectral subspace of the m largest eigenvalues).
   If the conjecture fails, report where and how (multiplicities, missing/extra triples) — that is a result too.
4. **Which eigenlines are physical, intrinsically?** In gl(3) the physical copies were the m LARGEST eigenvalues of V on W (E₂₃ sectors). For gl(4): is there an ordering or a
   sign rule on W that picks Spec(M_μ) (e.g. by V₃, by a fixed combination, by the sign of some entry of the eigenvector such as the component at the lattice point where the
   string cut acts)? Tabulate, for every sector of the products above, the rank of each physical triple in the orderings by V₃, V₄ᵃ, V₄ᵇ and by G2's generic combination, and look
   for a rule; test any candidate rule on held-out products ([3,2,1,0]⊗[4,2,1,0], [4,2,2,0]², [4,3,1,0]²).
5. **Chain structure and polynomiality.** (a) The sparsity graph of V_red,3 on R_ν (vertices = lattice points, edges = nonzero off-diagonal elements): connected? a path? a grid?
   same for the pair {V₃, V₄ᵇ}; when M_μ = W (m = w), is V_red,3 tridiagonal in some ordering of the lattice points (bandwidth after a Cuthill–McKee reordering)? (b) Is the k-th slice
   component of a copy a polynomial in its triple (v₃, v₄ᵃ, v₄ᵇ)? Fit G3's copies (50 digits, "first nonzero component = 1" and alternative normalisations) by polynomials of
   total degree ≤ m−1 with rational coefficients (least squares + rational reconstruction); report degrees or failure. (c) A cyclic vector that is a coordinate vector for the
   triple on M_μ or on W (then the copies are Lanczos polynomials in a single generic operator): test the coordinate vectors of R_ν.
6. **Closed-form multiplicity** (optional): m = w − #(spurious) — using question 3, express m as a signed count of lattice points Σ_w sgn(w)·|R_{μ−w·λ′}| (this is the Racah–Speiser/
   Kostant formula with the GT weight multiplicities; verify on all sectors as a consistency check of your implementation of the dot action).

## 1. Deliverables

`structure_gl4.py` (stages `stencil | symbolic | weyl | rule | chain | poly | all`), `structure_gl4_tables.txt` (every table), `test_structure_gl4.py` (a few-second regression
test: one stencil, one symbolic element on a held-out rep, the Weyl identity on one sector), `RESULTS-G4-structure-gl4.md`: the answers to 1–6 with numbers, the closed forms
that survived held-out tests (and the ones that did not), the Weyl/BGG result (the main table: per product, #sectors where the signed count is the indicator), the physical-rank
tables and any rule found, the chain/polynomiality findings, surprises, what is not done. Report failures with numbers; never paper over a failing test. Write Python files with
the file tools (Bash heredoc escape trap, hit by G1/G2/G3).
