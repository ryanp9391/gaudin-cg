# Spec W1: the sl(3) Gaudin chain, L = 2 — Hamiltonians, vacuum, localising operators, two-vacua untwist and CG (2026-10-15)

Delegated to a wolfbook-builder subagent (sonnet); planned by the main session. **New notebook:**
`Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain.wb` (self-contained; rep-generic with `Λ` as the only switch; first code
cell `Quiet[ClearAll["Global`*"]]`; [2,1,0] first, then [3,1,0] by changing the switch). Read first: `../Gaudin-sl3-CG-plan.md`
§0–§2 and §7 (the dictionary in §1 is binding). Reference cell sources (UTF-8) in `Paul/Mathematica/Gaudin/ref/`:
`ref-L2-two-vacua-untwist-all.txt` is the template for everything chain-side (cells 3, 5, 7, 9, 11, 13 …: generators and metric,
exact transfer matrices by interpolation, group data `expN`/`x0Of`/`ρdiag`/`rowExp`/`coverG`/`rGauss`/`nMat3`/`nInv3`/`ΩD`/`phiUof`/
`Pphys`/`ΩC`, HW kets and G-projectors `hwKets`/`PHWc`, Q machinery, `branchesAt`, K/S and CG cells). Copy the generic parts
verbatim; replace the XXX-specific parts (transfer matrices, twists, regime) by the Gaudin ones below. Everything exact where
possible (rational ζ, ω, θ, u); 200 digits where eigenstates enter; time one case before any loop; nothing over ~3 min without a
timed smaller run.

## 0. Conventions (from the plan §1; do not deviate silently — if something only works with a variant, say which)

- Λ = {2,1,0} (dim 8, 64 states) then {3,1,0} (dim 15, 225). θ = {2/7, −3/5}. Finite-twist point P: ζ = {3/2, 5/4, 2/3},
  ω = {5/2, 7/3, 3/5}, ω′ = {7/4, −2/5, 9/7}; control point P′: ζ″ = {2, 5/3, −1/2}, ω″ = {3, −1/3, 7/2}.
- Lax: `Lg[g][u]` = 3×3 table of d2×d2 sparse matrices, L_ij = g_ij·id2 + E1[i,j]/(u−θ₁) + E2[i,j]/(u−θ₂); derivatives
  L′_ij = −Σ_s E^{(s)}_ij/(u−θ_s)², L″_ij = 2Σ_s E^{(s)}_ij/(u−θ_s)³ (exact).
- cdet(∂ − L_g) = Σ_σ sgn σ · M_{σ(1)1}·M_{σ(2)2}·M_{σ(3)3}, M = ∂·1 − L_g, ∂ moved to the right with ∂f = f∂ + f′. Do the normal
  ordering ONCE symbolically on noncommuting symbols (e.g. `NonCommutativeMultiply` over symbols l[i,j], lp[i,j], lpp[i,j] with
  ∂ as a formal symbol `dd` and the rule dd ** f → f ** dd + f′; or expand by hand: each term is a product of up to three
  factors each equal to dd or −l[i,j]; move every dd to the right generating derivative terms), obtaining
  cdet = dd³ − S1·dd² + S2·dd − S3 with S1, S2, S3 noncommutative polynomials in {l, lp, lpp}. Print S2 and S3 in that symbolic
  form in a markdown cell (this is the "what are the Gaudin T's" record). Then `S[k][g][u]` = substitute the matrices.
- Checks (Task 0): S1 = (Σζ_k + |Λ|Σ_s 1/(u−θ_s))·id2 exactly; [S2(u₀), S2(u₁)] = [S2(u₀), S3(u₁)] = [S3(u₀), S3(u₁)] = 0 exactly at
  u₀ = 3/7, u₁ = −2/9 (and one more pair); S2, S3 conserve the weight for diagonal g; S2(u) − Σ_{i<j}(L_iiL_jj − L_jiL_ij)(u) is a
  combination of L′ terms and scalars (print it). **If the commutators are NOT zero**: try (i) the transposed placement
  L_ij → g_ji + Σ E^{(s)}_ji/(u−θ_s), (ii) the row-ordered product M_{1σ(1)}M_{2σ(2)}M_{3σ(3)}, (iii) ∂ + L; report which works and
  use it; the plan's expectation is the stated one.
- Laurent coefficients (the "T's"): `Hres[s]` = Res_{θ_s} S2; `S3c[s,k]` = coefficient of (u−θ_s)^{−k} in S3, k = 1,2,3; `S3inf` =
  lim_{u→∞} S3. Compute them exactly by evaluating S_k at several rational u and solving the partial-fraction ansatz (S2 has
  poles of order ≤ 2, S3 of order ≤ 3 at each θ_s, plus constants; degree in u of the numerator is bounded, so 7 resp. 9 sample
  points suffice — verify the ansatz by a further point). Check: S3c[s,3] and S3inf are scalars; Σ_s Hres[s] = Σ_kζ_k(E_tot)_kk
  + scalar; H₁ := Hres[1] equals Σ_{ij} ζ_i δ_ij E1[j,i] + Σ_{ik} E1[i,k].E2[k,i]/θ₁₂ + (print the difference).
- Bra vacuum: `gc[ζ]` = {{e1,−e2,e3},{1,0,0},{0,1,0}}; `PauxG[ζ]` rows (1, ζ_i − e1, ζ_i² − e1 ζ_i + e2); check PauxG.gc.PauxG⁻¹ =
  diag(ζ); g′ := PauxG[ζ].gc[ω].PauxG[ζ]⁻¹; Ω_{ζω} := `ΩD[ζ, ω]` with the reference `rGauss` formulas (z → ζ, w → ω, unchanged).
  Verify (Task 1) that Ω⊗Ω is the exact left eigencovector of S2^{(g′)}(u₀), S3^{(g′)}(u₀) (and at u₁). If not, test the variants
  PauxG⁻¹.gc.PauxG and the transposes of both; if still not, compute the exact common left eigencovector of S2^{(g′)}(u₀),
  S3^{(g′)}(u₀) in the top weight sector with the first choice of g′ and compare with ΩD (parallel?) — report. Also verify the
  one-site statement: with M the 3×3 matrix such that M.g′.M⁻¹ = diag(ω) (rows of M = left eigenvectors of g′ ordered ω₁, ω₂,
  ω₃, first entry 1), e₁ᵀρ(upper Gauss factor of M) ∝ ΩD[ζ,ω] (ρ of the Gauss factors as `rhoU`).
- Companion-frame vacuum (record only): P_phys^G(ζ) := the rep image of PauxG[ζ]⁻¹ (or PauxG, whichever makes the twist
  companion — mirror the reference `Pphys`/`x0Of` construction: Gauss-factorise the 3×3 matrix as P₋·diag(s)·P₊ and map each
  factor with `expN`/`ρdiag`); x0G := Flatten[KroneckerProduct[x, x]] with x = e₁ᵀP_phys^G⁻¹ (factorised inverse, never
  `Inverse` at small ε). Check x0G·Ψ_B ≠ 0 on all states at P.
- Regime: ζ = εφ, ω = εω̂, ω′ = εω̂′ with the plan's φ, grids, ε ∈ {1,2,4}·10⁻⁶ (`reg2`-type association). Note ΩD[εφ, εω̂] is
  exactly ε-independent (degree-0 homogeneity) — assert it.
- Eigenstates (Task 2): per weight sector (`secOf`), `Eigensystem` of N[S2(u₀) + π S3(u₁), pr] restricted to the sector, pr = 200;
  at the regime point the HW-limit branches = the dim(HW_μ) states of sector μ with the smallest raising residual; normalise
  each branch by G-norm 1 (`gNorm`) and sign (first nonzero component in GT order positive); copies ordered by Re of the
  eigenvalue. Store μ, copy, Ψ, eigenvalue, residual, ‖P_HWΨ‖/‖Ψ‖, and the reference-vacuum overlap ⟨ΩD[ζ,ω₀]⊗ΩD[ζ,ω₀]|Ψ⟩
  with ω₀ = ε·(2/5, 1/7) for the Q side. Also at the finite-twist point P: all 64 states (for Task 3 and for P1's use —
  print the 64 eigenvalue pairs (S2(u₀), S3(u₁)) to 30 digits in a markdown cell).

## 1. Tasks (each its own code cell + markdown result cell with a `Grid`; `checks` an Association; save after each)

**Task 0 — setup and the commuting family** (conventions above; symbolic S2, S3 printed; all exactness checks; timing of one
S3(u) evaluation).

**Task 1 — vacuum and frames** (eigencovector verification with the variant bookkeeping; the one-site Gauss statement; x0G;
the ε-independence assertion; `tiePts`-type checks at P, P′ and at a regime point).

**Task 2 — eigenstates, HW kets, projectors, HW-limit branches** (`hwKets`, `PHWc`, `branchesAt[ζ]` adapted; the 64-state table
at P; the branches at the regime point for ε = 10⁻⁶; multiplicity-2 sector [3,2,1]: the 2×2 matrix of S3c[1,2] (or of the full
S3(u₁)) on the exact HW space of the UNTWISTED model (ζ = 0): its eigenvectors are the Gaudin untwisted copies — print it
and them).

**Task 3 — the localising scan (which T's).** Exact, at P and P′. D₁ := independent subset of {id, W := Σ_kζ_k(E1+E2)[k,k],
H₁, S3c[1,2], S3c[1,1], S3c[2,2], S3c[2,1]} (report the rank; expect 7 incl. id); D₂ := D₁ ∪ {a_i.a_j, i ≤ j} (report rank;
≤ 28; also the rank of the cubic span, expected 64 = vacuous). Φ_X := reshape(⟨Ω⊗Ω|X) as d1×d1 (row = site 1, column = site 2).
Site-1 localisation ⇔ Φ_X·k = 0 for every ket k in the (d1−1)-dim annihilator of Ω₁ (NullSpace[{Ω₁}]): a homogeneous linear
system in the D₂ coefficients — `NullSpace` exactly; L₁ := the span of the resulting v = Φ_X·k₀/(Ω₁·k₀). Report ℓ₁ = dim L₁, a
basis, which site-1 weights of V are covered (rank of L₁ projected on each weight), and for each basis element the
coefficient vector in D₂ (sparse representatives: solve for the minimal number of nonzero coefficients among {D₁ only},
{D₁ + products with H₁}, {all}). Same for site 2 (Φ_X = Ω₁⊗u → L₂, ℓ₂). Then IDENTIFY: test whether each basis covector of L₁
lies in Ω₁·span{id, C_ζ := Σ_kζ_kE[k,k] (one-site), the one-site S2/S3 Laurent coefficients at θ₁ for the one-site Gaudin model
with twist diag(ζ) (build the one-site cdet the same way with L¹_ij = ζ_iδ_ij + E_ij/(u−θ₁)), Σ_{ij}(g − g′)_ij E[j,i]}; print the
decomposition. Expectation (from XXX, not guaranteed): ℓ₁ = 3 with L₁ = Ω₁·span{1, one-site S3 coefficients}, the analogue of
{Ω, ΩT₁₁, ΩT₂₁}; and the spectator factor of each localising X (the scalar f with ⟨Ω⊗Ω|X = f·v⊗Ω₁ when v is normalised
to first component 1) as a rational function of θ₁₂ (fit at P and P′ and one more θ; expect a polynomial in θ₁₂ of low degree —
the t → ∞ limit of f₍₁,₁₎ = (θ₁₂−h)(θ₁₂−3h) is θ₁₂²). Keep the scan ≤ ~5 min total (exact 64×64 with ≤ 28 operators is cheap;
for [3,1,0] use SparseArray and time one product first).

**Task 4 — the two-vacua family: exact undo, K/S, CG at leading order** (the Stage D pipeline, cells 15–33 of the reference,
with the Gaudin branches). Φ = ΩD[ζ,ω]⊗ΩD[ζ,ω′]. (a) Exact undo ⟨Φ|h⟩ = ⟨e₁ᵀρ(n(ζ,ω)n(ζ,ω′)⁻¹)⊗e₁ᵀ|h⟩ on all HW kets at P, P′ and
a regime point (exact; control: the naive pairing differs). (b) The undone unipotent n′ at the regime point is exactly
ε-free: print a′, c′, g′ = b′−a′c′ in closed form in (φ, ω̂, ω̂′) (`Factor`), the rank of (ω̂₁,ω̂₂,ω̂₁′,ω̂₂′) ↦ (a′,c′,g′) (expect 3), the
per-weight ranks of φ_u = e₁ᵀρ(n′) over the 24-point grid (expect rank = multiplicity on every weight of V). (c) Chain side at
the grid × ε: R_C = ⟨Φ|Ψ_B⟩, S = ⟨Φ|P_HWΨ_B⟩ (= ⟨φ_u⊗e₁ᵀ|P_HWΨ_B⟩, check), K = R_C − S; ε-exponents of |S|, |K| (log₂ of ratios
between consecutive ε), max |K/S| per branch over the grid at each ε. **Expectation to confirm or refute: K/S = O(ε) for EVERY
HW branch** (adjoint: 6, incl. the singlet and both octets). If a branch fails: the `showContam` breakdown by site-1 weight, no
fitting. (d) CG at leading order from the overlap alone: for 1-dim slices C_B := R_C/[φ_u]_slice — grid spread, error vs the
chain slice ⟨m_slice⊗e₁|P_HWΨ_B⟩ at two ε's (both ∝ ε), Richardson limits vs the slice of the exact normalised HW ket (the
truth: compute the exact HW kets' slices directly — these ARE the CG coefficients in the GT basis with G-norm 1); 2-dim slices
(octets): the 2×2 system from two grid points → both components, and the mixing R_c = (u·C_c)/(u′·C_c) with u = e₁ᵀE₁₃,
u′ = e₁ᵀE₂₃E₁₂ against the exact untwisted copies of Task 2. (e) Controls: ω′ fixed at (7/4, −2/5) (not untwisted) — expect
some slices to drop order and K/S ↛ 0 for the deep branches; common ω̂₁′ = ω̂₁ — expect a′ ≡ g′ ≡ 0 (reads one irrep). (f) Record
for the Q side: at one grid point and ε = 10⁻⁶ the four numbers ⟨Ω_ω⊗Ω_ω|Ψ⟩, ⟨Ω_{ω′}⊗Ω_{ω′}|Ψ⟩, ⟨Ω_ω⊗Ω_{ω′}|Ψ⟩, ⟨Ω_{ω′}⊗Ω_ω|Ψ⟩
per HW branch (G-normalised branches), with ε-exponents.

**Task 5 — summary cell** `### sl(3) Gaudin L = 2 chain result (2026-10-15)`: the symbolic S2, S3; the variant that works (if
any deviation); the eigencovector result; ℓ₁, ℓ₂ and the identified localisers with spectator factors; the K/S table; the CG
table (overlap-extracted vs exact); surprises. Then set Λ = {3,1,0}, rerun Tasks 0–4 (time the S3 construction first; reduce the
grid to 2×2×… if the Task 4 loop exceeds 10 min), append the [3,1,0] results to the summary. Save.

## 2. Rules and gotchas

- Create the notebook as JSON (cells: kind 1 markdown / 2 code, `value`, `languageId` "wolfram", `outputs` [], `metadata` {}),
  open it with `code <path>`, evaluate each code cell's exact saved source in K1 via `wolfbook_evaluateExpression`
  (`multiLine` false, generous `timeoutSeconds`), store the printed results as a markdown cell under the code cell, save by
  rewriting the JSON (or `wolfbook_getNotebookContext` action "save" if available). Say in the report that there are no native
  outputs. Do not evaluate cells of any other notebook; do not edit other notebooks, `paper/*.tex`, HANDOFF.md, or git.
- Shared kernel K1 is the only kernel; the first cell wipes it. Semicolons between statements. Never the two characters
  `*)` inside a comment. Do not use `D`, `E`, `N`, `I`, `C`, `K` as symbols (`KK`, `kk`, `dd`). `Max[{}]` = −∞ (floor at 0).
  Memoised definitions cache bad results — re-run the definition and everything downstream after a fix; `Do[DownValues[f]=…]`
  purges nothing. Exact zero tests with `===` on rationals; on numerics report `Precision`; `SetPrecision` before `Norm` of
  near-zero vectors; never `PossibleZeroQ` on numerics.
- The vacuum degenerates when ω_i = ζ_j; the plan's points avoid this; check any new point.
- Noncommutative normal ordering: keep it to the 3×3 symbolic level (27 symbols); never try to normal-order with 64×64 matrices
  inside symbolic expressions.
- Use `SparseArray` for E1, E2, L_ij and the products; convert to dense `N[…, pr]` only per weight sector.
- Report numbers, not adjectives; every residual with its `Precision`.
