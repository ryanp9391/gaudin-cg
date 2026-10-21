# Spec D3: the Gaudin B operator, the jet (SoV) basis, the Gaudin wavefunction rule, the Gaudin K tensor, and G from first principles (2026-10-16)

Delegated to a general-purpose subagent (sonnet), Python only, exact arithmetic (Fraction / sympy) wherever possible, high-precision mpmath
only as a fallback. Read first: `Paul/Mathematica/Gaudin/Gaudin-sl3-round3-plan.md` §0 Q2 (the question), §3, §6; the round-2 plan
`Gaudin-sl3-round2-plan.md` §1 with its "Added by D2" and "Added by D1" blocks (the jet formula you are deriving structurally) and
`Paul/Python/gaudin_sl3/RESULTS-twovac-limit.md` (D2: the formula, the 25 entries of G, the Plücker/ν bases, the role of minor_s(−1,0) = 0),
`RESULTS-twovac-fit.md` (D1: κ_s = a_s/3, jet ranks, the N₅ basis), `RESULTS-xxx-limit.md` (P2: the XXX port, x₀ ~ h per site → x₀^G,
notebook-order vs gl3_reps-order of the GT basis). The XXX SoV construction you are limiting: `Paul/Mathematica/Gaudin/ref/
ref-L2-rotation-sov-all.txt` (UTF-8 cell dump; line ~48 `qminor`/`tAt`, line ~231 `Bop2[v] := TMf[Id3][1,1][v].qminor[Id3][{1,2},{1,2}][v+hN]
+ TMf[Id3][2,1][v].qminor[Id3][{1,2},{1,3}][v+hN]` with TMf the UNTWISTED monodromy T = ℒ⁽¹⁾(u−θ₁)ℒ⁽²⁾(u−θ₂), the label rule `labRule`, the
SoV covectors `vecB`, the wavefunction test `wfTest`, the one-site monomials `siteMono`), and `Paul/Mathematica/new-conventions/
Two-site-rotation-brackets-notes.md` §1 ("SoV labels, wavefunctions, ladders"; "Two-site covectors, K tensor, blocks") and §5.

**Code to reuse (read-only, import or copy; never edit):** `xxx_chain.py` (P2: `Chain` with `lax`, `mono_untwisted`, `qmin`, `t_a`; the sparse
`SpOp` class; the formal-series class `Ser` in h), `xxx_limit.py` (`mp_x0_onesite`, `Paux_exact`, `Gc_exact`), `gaudin_chain.py` and `qside.py`
(P1: Gaudin states and dual/direct Q's at 60 digits; HW index 7 in gl3_reps order), `twovac_fit.py` (D1: `jets`, `plucker`, `x0_onesite`,
`onesite_inf`, `cov_inf`, `kappaG`, `cscal`, `d2_G` — the 25 entries of G as functions). **Files you own:** `gaudin_B.py` (runnable
top-to-bottom, `python -X utf8 gaudin_B.py`, stage flags), caches `gaudin_B_*.pkl`, `RESULTS-gaudin-B.md`. Environment: Python 3.13, sympy
1.14, mpmath 1.3; no `uv`. Λ = (2,1,0), L = 2, θ = (2/7, −3/5), ζ = (3/2, 5/4, 2/3) for the eigenstates (companion twist frame: Ψ_c =
(P_phys^{G −1}⊗P_phys^{G −1})ψ_d with P_phys^G = ρ((P_aux(ζ)⁻¹)ᵀ); x₀^c = e₁⊗e₁ so ⟨x₀^c|Ψ_c⟩ = ⟨x₀^G|ψ_d⟩). Mind the basis order throughout
(xxx_chain / notebook order has the HW ket first; gl3_reps order has it at index 7).

## 0. What is known (do not re-derive; verify where cheap)

XXX (new conventions, untwisted T): B(u) = T₁₁(u)·qmin_{12,12}(u+h) + T₂₁(u)·qmin_{12,13}(u+h), qmin_{I,J}(v) = Σ_σ sgn σ ∏_k T_{σ(i_k)j_k}(v−(k−1)h)
(rows antisymmetrised, column k at v−(k−1)h); [B(u),B(v)] = 0; B(u) is a polynomial of degree 6 in u; its spectrum is simple: eigenvalue
∏_{s=1,2}∏_{k=1}^{3}(u − θ_s − m_k^{(s)}h) with per-site labels (a|bc), roots m ∈ {−1−a, −b, −c}, λ₂ ≤ a ≤ λ₁, 0 ≤ c ≤ λ₂, c ≤ b ≤ a (8 labels
per site for the adjoint, 64 two-site labels; the labels are u-independent, the eigencovectors ⟨x| are common to all u). Wavefunction (companion
frame, companion-twisted eigenstates Ψ_c): ⟨x|Ψ_c⟩/⟨x₀|Ψ_c⟩ = N_x ∏_s M^{(s)}_x/M^{(s)}_vac, M_{(a|bc)} = minor_s(−1−a,−b)·Q³(θ_s−ch),
minor_s(m,m′) = Q²(θ_s+m′h)Q³(θ_s+mh) − Q³(θ_s+m′h)Q²(θ_s+mh), Q^j the dual Q's of the state, x₀ = the label with a = λ₁, b = c = λ₂ on both sites
(check: M_vac = minor(−3,−2)Q³(θ−h) for the adjoint), N_x normalisation constants absorbed into the K tensor. Two-site covectors Φ = Σ_x φ_x⟨x|/N_x
with φ = Φ·K (K the exact rational, z-independent 64×64 "K tensor"); ⟨Φ|Ψ⟩/⟨x₀|Ψ⟩·M¹_vacM²_vac = Σ φ_{x₁x₂}M_{x₁}M_{x₂}.
Gaudin (round 2): the overlap of the (∞,∞) family is 4c·Qnum_G·(P³p₀₂)(θ₁)(P³p₀₂)(θ₂) = Σ G_ij N_i(θ₁)N_j(θ₂), N = (P³p₀₂, P³′p₀₂, P³p₀₃, P³′p₀₃, P³p₂₃),
p_kl = P²^{(k)}P³^{(l)} − P³^{(k)}P²^{(l)} at θ_s, with p₀₁ = 0, p₁₂ = κp₀₂, κ_s = a_s/3; it is the h¹⁰ term of the XXX identity, and the one relation
Q^j(θ_s−h) = ρ_sQ^j(θ_s) ⇔ minor_s(−1,0) = 0 explains all cancellations. In the Gaudin limit all ladder points collapse to θ_s, so B(u) =
ν(u)³[1 + hb₁(u) + h²b₂(u) + …] with b₁ Cartan-like and degenerate (eigenvalue Σ_s(1+a+b+c)/(u−θ_s)): the Gaudin SoV "basis" must be a jet basis.

## 1. Tasks

**B0 — port and verify B(u).** Build B(u) from P2's `Chain` (untwisted, exact rationals at h = 1/2 first): [B(u₀),B(u₁)] = 0 exactly; degree 6 in u
(interpolate from 7 points); the 64 eigenvalues at u₀ = 3/7 equal the label-rule products (list the 8 one-site labels and the 64 pairs); the
left eigencovectors ⟨x| (exact nullspaces of B(u₀)ᵀ − λ_x) are eigencovectors of B(u₁) too (u-independence). Verify the wavefunction rule on
the companion-frame XXX eigenstates at h = 1/2 for a few states (needs XXX dual Q's: solve the dual Baxter equation ν₃(x)ν₃(x−h)ν₃(x−2h)Q(x+h) −
t₁(x)ν₃(x−h)ν₃(x−2h)Q(x) + t₂(x)ν₃(x−2h)Q(x−h) − t₃(x)Q(x−2h) = 0 for Q^j = z_j^{x/h}q^j, deg q^j = μ_j, at 60 digits with the companion twist G_c(z)
— or skip this check if it costs more than ~30 min and say so; the rule is established in the notebooks).

**B1 — the h-expansion of B.** With `Ser` (formal series in h, sparse rational matrix coefficients) at rational u: B(u) = ν(u)³[1 + hb₁(u) + h²b₂(u)
+ h³b₃(u) + h⁴b₄(u) + …] (as many orders as the structure needs; at least 4). Express b₁(u) in the Gaudin Lax entries L_ij(u) = Σ_s E^{(s)}_ij/(u−θ_s)
(expect a Cartan combination such as 2L₁₁ + L₂₂ — fit exactly in the basis {L_ij(u), L′_ij(u), 1}), and b₂(u) in {products L_ijL_kl, L′_ij, L″_ij, 1}
(exact linear algebra); report the commutators [b_i(u), b_j(v)] (which vanish, which combine to zero order by order as [B(u),B(v)] = 0 demands:
[b₁(u),b₁(v)] = 0; [b₁(u),b₂(v)] + [b₂(u),b₁(v)] = 0; …); the spectrum of b₁(u₀) (eigenvalues and multiplicities; compare with Σ_s(1+a+b+c)/(u₀−θ_s))
and whether b₂ restricted to each b₁-eigenspace splits it further (joint generalised eigenspaces of (b₁, b₂, b₃): dimensions; are the b_k
simultaneously triangularisable? is the joint structure a flag indexed by the labels?).

**B2 — the SoV covectors as functions of h and the jet basis.** For each of the 64 labels x, the exact left eigencovector ⟨x|(h) of B(u₀) at
rational h ∈ {1/2, 1/3, 1/4, …, 1/20} (nullspace of B(u₀)ᵀ − λ_x(h)·1 with the exact eigenvalue from the label rule; one-dimensional for h ≠ 0),
normalised by ⟨x|(h)·r = 1 for one fixed generic rational vector r (e.g. r_i = 1/(i+2)); then (a) reconstruct each component as a rational
function of h (Thiele continued fraction or linear solve for numerator/denominator coefficients with increasing degree bounds; validate on two
held-out h values; report the degrees), and expand at h = 0: ⟨x|(h) = ⟨x|₀ + h⟨x|₁ + h²⟨x|₂ + …; (b) the leading covectors ⟨x|₀: how many distinct
ones among the 64 (expect heavy coincidences — report the classes, e.g. by total excitation Σ(a+b+c) or by weight), the ranks of span{⟨x|₀},
span{⟨x|₀, ⟨x|₁}, … (the filtration), and the jet basis: a canonical basis of (V⊗V)* obtained by divided differences along the label lattice
(for a one-site ladder of labels x⁽⁰⁾, x⁽¹⁾, … whose leading covectors coincide, take ⟨x⁽¹⁾|−⟨x⁽⁰⁾| over h, etc.) or, equivalently, the joint
generalised eigenvectors of (b₁, b₂, b₃) ordered by the filtration — state the construction you used and verify that it is a basis (rank 64).
If the rational reconstruction is too hard (degrees > 30), fall back to the series directly: ⟨x|_k from degenerate perturbation theory
(the eigenvalue series is known exactly), or to 100-digit numerics with Richardson in h — say which.

**B3 — the Gaudin wavefunction rule.** In the jet basis {⟨y|} (y indexing the 64 jet-basis covectors), the Gaudin eigenstates (P1's states at
ζ = P, companion frame Ψ_c) have components ⟨y|Ψ_c⟩/⟨x₀^c|Ψ_c⟩. Hypothesis (from the XXX rule and the ladder → jet dictionary): each is a
product over the two sites of a jet monomial of the dual Q's at θ_s — a Plücker coordinate p_kl(θ_s) (or a jet of W = P²P³′ − P³P²′) times a
derivative P³^{(c)}(θ_s) — divided by M_vac^G = P³(θ_s)p₀₂(θ_s), with state-independent coefficients. Find the rule by exact fitting per
jet-basis element (64 states × the candidate monomials {P³^{(c)}p_kl, c ≤ 1, k < l ≤ 3} per site — 12 per site, products across sites 144; but a
jet-basis covector should be a short combination): report for each y the monomial(s) and coefficients, and the verification on all 64 states
(60 digits). This is the Gaudin SoV wavefunction: ⟨y|Ψ_c⟩ ∝ ∏_s (jet monomial)_s. If the basis of B2 is not adapted (each y needs many monomials),
re-choose the jet basis so that it is (the wavefunction rule defines the canonical jet basis) and report both.

**B4 — the Gaudin K tensor and the derivation of G.** Expand the companion-frame (∞,∞) covector Φ_c∞ = Ωc∞[ω₁]⊗Ωc∞[ω₁′], Ωc∞[ω₁] =
e₁ᵀρ(e^{E₁₂/ω₁}e^{E₂₃/ω₁}) (W2: gaussCov{1/ω₁,0,1/ω₁}), in the jet basis: Φ_c∞ = Σ_y φ_y(ω₁,ω₁′)⟨y| (exact rational in ω₁, ω₁′ — the Gaudin K tensor
is this change of basis, ζ-independent). Then ⟨Φ_c∞|Ψ_c⟩/⟨x₀^c|Ψ_c⟩ = Σ_y φ_y·(wavefunction rule)_y — rewrite in D2's basis N and compare with
D2's G (25 entries, `twovac_fit.d2_G`), including the prefactor 4c·(P³p₀₂)₁(P³p₀₂)₂ (c is the frame-tie scalar κ^Gκ^G′ with κ^G = Δ₁Δ₂ — check
how it arises here). Report: agreement entry by entry (exact if your φ_y are exact), the number of nonzero φ_y (the Gaudin analogue of the XXX
"ten products"), their structure in (ω₁, ω₁′) (homogeneity, the weights (2,3,3,4,5) of D2), and which jet-basis elements couple the two sites.
Also the equal-twist case Φ_c = Ωc^G[ω]⊗Ωc^G[ω] (finite ω₂): φ_y and the resulting Q-number vs P1's det N (should be proportional; report the
state-independent ratio).

**B5 (optional) — a single Gaudin B operator.** Is there a polynomial B_G(u) in the Gaudin L_ij(u), L′_ij(u) (degree ≤ 3) whose joint generalised
eigenspaces for all u reproduce the jet-basis filtration (i.e. B_G(u) is upper-triangular in the jet basis with a known diagonal)? Test the
natural candidates from B1: b₂(u) − (its Cartan/scalar parts), the "connected" part b₂ − b₁²/2, and the leading coefficient of B(u) − ν³(1 +
hb₁) after removing the commuting Cartan part; report what you find or that none of the candidates works. Also state the operator-valued
coefficients of u in b₁, b₂ (B(u) is a polynomial in u; the SoV basis diagonalises every u-coefficient simultaneously).

## 2. Reporting

`RESULTS-gaudin-B.md`: B0 (verification numbers), B1 (the b_k in Lax terms, commutators, spectra, joint eigenspace dimensions), B2 (degrees in h,
the leading-covector classes, the filtration ranks, the jet basis construction), B3 (the wavefunction rule: a table y → monomials, verification),
B4 (K tensor: nonzero φ_y count and structure, G reproduced entry by entry or where not, the equal-twist comparison), B5 (what was tried),
surprises, not done. Exact where possible; numbers with precision. No notebooks, no git, no Wolfram, no edits outside your files. If an
expectation is refuted and two attempts do not fix it, stop there, keep what works, report what you saw.
