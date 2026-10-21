# Clebsch–Gordan coefficients of gl(3) from separation of variables in the two-site Gaudin model: working notes on rounds 1–4

*Paul Ryan (project notes, compiled from the work of the 2026-10-15/16 agent rounds) — 16 October 2026*

*(Markdown twin of `Gaudin-sl3-CG-notes.tex`; same structure and numbering. Math is Unicode; displayed equations are numbered as in the PDF.)*

## Abstract

We want the Clebsch–Gordan coefficients of gl(3) for V_Λ⊗V_Λ, and we want them as explicit expressions in Q-functions, obtained from separation of variables (SoV) in the two-site rational Gaudin model. These notes collect everything that four rounds of computation (chain side in Wolfram, Q side and kernels in Python, about fourteen agents) have established, and say honestly what has not been established.

The main finding is that the overlap of a CG-complete two-vacua covector family with an eigenstate of the Gaudin Hamiltonians is a Q-number: a bilinear form, with a state-independent kernel, in the Plücker coordinates of the order-(λ₁+1) jets at θ₁ and θ₂ of the dual Q-functions. For the adjoint the kernel is a 5×5 matrix of polynomials in (ω₁,ω₁′,1/θ₁₂); it was found as the formal h→0 limit of the XXX ten-product formula, re-found by a direct fit, and derived from first principles from the Gaudin jet basis. Writing the Plücker coordinates through the Frobenius solution at θₛ shows that the overlap depends only on the Laurent data (Hₛ,aₛ) of the Hamiltonians at the two sites. Taking the untwisting limit ζ→0 of this expression gives, for all six highest-weight branches of the adjoint, the CG generating function of adj⊗adj as an exact polynomial identity, once the ε-series of the Laurent data are known; those series are exact, polynomial in the twist direction, and algebraic with a square root on the (4,1,1) and (3,3,0) branches. For general Λ we have closed forms for the one-site kernel, a block theorem reducing every two-row representation to the one-row master kernel, and the master kernel to first order in 1/θ₁₂ (all orders for [l,l,0]). No per-sector cleaning of the XXX kind is needed anywhere.

What is open: the master kernel beyond first order, the Laurent-data formulas and the limit for λ₁≥3, a derivation of the ε-series, and the norm of the untwisted eigenvectors from the Q side. Appendices hold all the data (kernels, series, CG tables, untwisted Laurent data) and a reproducibility map.

## Contents

- 0. Overview: what the program is, what we found, what is open
  - 0.1. The program in a nutshell
  - 0.2. What we found, in eight points
  - 0.3. What is open
  - 0.4. How to read these notes
- 1. Setting and conventions
  - 1.1. The model and its Hamiltonians
  - 1.2. Laurent data
  - 1.3. The dictionary to XXX
  - 1.4. Frames, vacua and the reference covector
  - 1.5. The XXX scaffold in one paragraph
  - 1.6. Representations, bases and numerical points
- 2. Q-functions, Wronskians, detN
  - 2.1. The oper and its Q-functions
  - 2.2. Local exponents and Frobenius constants
  - 2.3. The no-log conditions: Gaudin's Bethe equations
  - 2.4. Wronskians and quantum Weyl
  - 2.5. Orthogonality and the detN overlap
  - 2.6. The rank-4 rule
  - 2.7. Vacuum Q-functions and vacuum Laurent data
- 3. Which operators localise: the answer to “which T's?”
  - 3.1. The question
  - 3.2. The scan
  - 3.3. Spectators and closed forms
  - 3.4. Why Gaudin differs from XXX
- 4. The two-vacua family and the untwist: the chain side
  - 4.1. Why a two-vacua family
  - 4.2. The family, the exact undo and CG-completeness
  - 4.3. K/S and the CG data
  - 4.4. Multiplicity-2 sectors and universal mixing
  - 4.5. Duality [3,2,0]≅[3,1,0]*⊗det³
  - 4.6. The controls and the “no cleaning” answer
- 5. The jet formula: the overlap as a bilinear form in jets of the Q-functions
  - 5.1. The statement
  - 5.2. The adjoint in detail
  - 5.3. Route 1: the h→0 limit of the XXX ten-product formula
  - 5.4. Route 2: the direct fit
  - 5.5. Route 3: derivation from the jet basis
  - 5.6. General Λ: [3,1,0] and [3,2,0]
  - 5.7. Kernel laws
  - 5.8. What went wrong first
- 6. The Gaudin B operator and the jet basis
  - 6.1. The XXX operator and its eigenbasis
  - 6.2. The hierarchy b₁,b₂,…
  - 6.3. Where the degeneracy lives
  - 6.4. The canonical jet basis and its wavefunction
  - 6.5. The K tensor and the derivation of G
- 7. The overlap as a function of the Laurent data
  - 7.1. The Frobenius reduction (adjoint)
  - 7.2. Why this must hold for every Λ
- 8. The general kernel
  - 8.1. Conventions for this section
  - 8.2. The one-site kernel
  - 8.3. The block theorem
  - 8.4. Structure: polynomiality, weights, degrees
  - 8.5. The first-order term
  - 8.6. The single-index family [l,l,0] to all orders
  - 8.7. The top-state sum rule
  - 8.8. One-row representations are the master kernels
  - 8.9. The auxiliary sl₂ and what is open
- 9. The untwisting limit and the CG generating function
  - 9.1. The regime and the difficulty
  - 9.2. ε-order bookkeeping
  - 9.3. The leading polynomial identity
  - 9.4. The reduced form: only Laurent data
  - 9.5. The CG dictionary and φ-independence
- 10. The ε-series of the Laurent data
  - 10.1. The series
  - 10.2. How they were obtained
  - 10.3. The perturbative content: Ψ₁ and X_B
  - 10.4. Worked numbers for adj⊗adj
- 11. The untwisted model
  - 11.1. Laurent data from Casimirs
  - 11.2. Fuchsian at infinity: three polynomial solutions
  - 11.3. The ₃F₂ criterion
  - 11.4. Why the twisted monic Q's have poles in ε
- 12. Surprises and corrections along the way
- 13. Outlook: the round-5 plan
  - 13.1. WP-A: the master kernel in closed form
  - 13.2. WP-B: Laurent-data form and the limit for [3,1,0], [3,2,0]
  - 13.3. WP-C: the mechanism of the ε-series
  - 13.4. WP-D: the norm from the Q side
  - 13.5. WP-E, WP-F: the write-up and housekeeping
- A. Conventions and normalisations
  - A.1. The kernel normalisations
  - A.2. Frames and bases
  - A.3. Numerical points and grids
  - A.4. Vacuum eigenvalue functions
- B. The kernels
  - B.1. The adjoint kernel
  - B.2. [3,1,0] and [3,2,0]: base-point kernels
- C. The one-site kernel
- D. The ε-series and the CG tables
  - D.1. CG tables
- E. The untwisted data
- F. Reproducibility

## 0. Overview: what the program is, what we found, what is open

**Punchline.** For the two-site gl(3) Gaudin model, the overlap of a suitable pair of rotated vacua with an eigenstate is a Q-number that depends on the eigenstate only through the Laurent data of the Gaudin Hamiltonians at the two sites; sending the twist to zero turns this Q-number into the CG generating function. We have this fully for the adjoint, with the main structural pieces (kernel, Frobenius reduction, representation-independent rules) in place for [3,1,0], [3,2,0] and, at the level of kernels, many more.

### 0.1. The program in a nutshell

The goal is simple to state. Take two copies of an irreducible gl(3) representation V_Λ, Λ=(λ₁,λ₂,0), and ask for the Clebsch–Gordan coefficients ⟨e_ν⊗e₁|Ψ₀⟩ of the highest-weight (HW) vectors Ψ₀ in V_Λ⊗V_Λ. The idea is to realise V_Λ⊗V_Λ as the Hilbert space of a two-site Gaudin chain with inhomogeneities θ₁,θ₂ and a diagonal twist g=diag(ζ₁,ζ₂,ζ₃). The Gaudin Hamiltonians are the Laurent coefficients of the quantum characteristic polynomial of the Lax matrix L_g(u)=g+ΣₛE⁽ˢ⁾/(u-θₛ). Their eigenvectors for generic ζ are labelled by Q-functions, and as ζ→0 (in a controlled direction ζ=εφ) the eigenvectors that are regular become the HW vectors of the untwisted V_Λ⊗V_Λ. SoV tells us how to write the overlap of any covector with such an eigenvector as a determinant or bilinear form of Q-functions. If we pick covectors whose overlaps with the untwisted HW kets are exactly the CG coefficients (suitably rotated), the CG coefficients become Q-expressions.

This is the “XXX program” of the earlier notes (Lax matrix u+hE, companion twist, SoV bases, localising transfer matrices), transplanted to the rational Gaudin model, which is the h→0 limit of the XXX chain with twists z=e^(hζ). We use the XXX chain as a scaffold (Section 1.5) but the deliverable is Gaudin.[^1]

### 0.2. What we found, in eight points

1. **The chain side works without cleaning.** The family Φ_∞(ω₁,ω₁′)=Ω^∞[ζ,ω₁]⊗Ω^∞[ζ,ω₁′] of products of two rotated vacua is CG-complete, its exact undo is ⟨φᵤ⊗e₁| with ε-free Gauss parameters, and the contamination ratio K/S is O(ε) on every HW branch of [2,1,0], [3,1,0], [3,2,0], including the multiplicity-2 sectors that needed per-copy cleaning in XXX (Section 4).
2. **The overlap is a Q-number (the jet formula).** Up to a frame-tie scalar c and a vacuum monomial M, the overlap is a bilinear form, with a state-independent kernel K, in Plücker coordinates of the jets of the dual Q-functions at θ₁ and θ₂, equation (1) below (Section 5).
3. **Where the kernel comes from.** The adjoint kernel was obtained three independent ways: as the analytic h→0 limit of the XXX ten-product formula, by a direct fit, and as the Gaudin K tensor in the jet basis that replaces the SoV basis when h→0 (Sections 5 and 6).
4. **The overlap depends only on the Laurent data** (Hₛ,aₛ) at θ₁,θ₂, because every Plücker ratio is fixed by the Frobenius recursion at θₛ (Section 7).
5. **The kernel has closed forms.** One site: a decomposable bivector, explicitly a dressed Wronskian. Two sites: a block theorem relating every [λ₁,λ₂,0] to the one-row master kernel, the first-order term in 1/θ₁₂ for all pairs, and all orders for [l,l,0] (Section 8).
6. **The untwisting limit is established for the adjoint**, as a polynomial identity in the rotation parameters, for all six HW branches; the CG coefficients come out φ-independent to 77 digits (Section 9).
7. **The ε-series of the Laurent data are exact** and partly algebraic (Section 10), and the untwisted model has a rich structure of its own: Fuchsian at infinity, three polynomial solutions, Casimir closed forms (Section 11).
8. **Several of our own early expectations were wrong**, and we collect them in Section 12 because they are as instructive as the successes.

The formula of point 2, written out once so that the reader has it in front of them, is

```text
c (⟨Φ_∞|Ψ⟩)/(⟨x₀^G|Ψ⟩)  M(θ₁)M(θ₂)=Σ_((kl),(k′l′))K_((kl),(k′l′))(ω₁,ω₁′,θ₁₂) (P³pₖₗ)(θ₁) (P³p_(k′l′))(θ₂),        (1)
```

with M=P³p_(0,λ₂+1). Here Φ_∞ is the two-vacua covector, Ψ any eigenstate, x₀^G a reference covector, P²,P³ two of the three dual Q-functions of Ψ, pₖₗ=P²⁽ᵏ⁾P³⁽ˡ⁾-P³⁽ᵏ⁾P²⁽ˡ⁾ their Plücker coordinates at θₛ, and (kl) runs over the C(λ₁+2,2)-C(λ₂+1,2) pairs 0≤k<l≤λ₁+1 with l≥λ₂+1. Every symbol is defined properly in Sections 1 to 5.

### 0.3. What is open

The master kernel R^(λ₁) for all four jet indices beyond first order; the Laurent-data (Frobenius) form of the overlap and the untwisting limit for λ₁≥3; a derivation, not a fit, of the ε-series of the Laurent data (the degenerate solver for the top state and the octets); the norm ⟨Ψ₀|Ψ₀⟩ from the Q side, which would remove the last piece of chain input; the multiplicity-3 sector [4,2,0]; and the write-up as a paper. Section 13 turns this into the round-5 plan.

### 0.4. How to read these notes

Each section opens with its punchline in two or three sentences, then tells the story, then gives the data. The main text works one example, the adjoint [2,1,0], all the way through; bulk data live in the appendices (A conventions, B kernels, C one-site kernel, D Laurent series and CG tables, E untwisted data, F reproducibility). A reader who skips all appendices still gets the whole story. Every number carries a provenance footnote: a file and a section or cell. Where sources disagreed, the later round, or the explicit reconciliation in `HANDOFF.md`, wins, and we say so in a footnote. Table 1 is the map of what has been done for which representation.

**Table 1.** Status map. “Chain” = chain-side checks (K/S, CG at leading order, Section 4); “Q” = Q-side rules (oper, Wronskians, detN, Section 2); “jet” = jet formula verified; “kernel” = closed forms (Section 8); “limit” = analytic untwisting limit (Section 9); “exact” = exact untwisted reference data (Appendices D, E). All entries are ours to the precision stated in the text; “fit” = determined by a unique exact fit, not derived.

| Λ | dimV | states | HW branches | chain | Q | jet | kernel | limit | exact |
|---|---|---|---|---|---|---|---|---|---|
| [1,1,0] | 3 | 9 | 2 | – | – | fit | exact | – | yes |
| [2,1,0] (adjoint) | 8 | 64 | 6 | yes | yes | derived | exact | yes | yes |
| [2,2,0] | 6 | 36 | 3 | – | – | fit | exact | – | yes |
| [3,1,0] | 15 | 225 | 10 | yes | yes | fit | laws | – | yes |
| [3,2,0] | 15 | 225 | 10 | yes | yes | fit | laws | – | yes |
| [4,1,0] | 24 | 576 | 14 | – | – | base point | block thm | – | yes |
| [l,0,0], l≤4; [3,3,0], [4,4,0] | 3–15 | 9–225 | – | – | – | fit | exact | – | – |

## 1. Setting and conventions

**Punchline.** The commuting family is the set of Laurent coefficients of a column-ordered quantum determinant; on any state it is described by six numbers (Hₛ,aₛ,bₛ)_(s=1,2), with a scalar double pole c₂=λ₂(λ₁+1) and no triple pole. The bra vacua are e₁ᵀρ(n) with explicit Gauss parameters, and they live in a frame tied to the diagonal frame of the eigenstates by an explicit scalar κ^G.

### 1.1. The model and its Hamiltonians

We work with gl₃ and a two-site chain; both sites carry the same irreducible representation V_Λ, Λ=(λ₁,λ₂,λ₃) with λ₃=0 throughout.[^2] On site s the generators are E⁽¹⁾=E⊗1, E⁽²⁾=1⊗E, they obey [Eᵢⱼ,Eₖₗ]=δⱼₖEᵢₗ-δₗᵢEₖⱼ, and the inhomogeneities are θ=(θ₁,θ₂), θ₁₂=θ₁-θ₂. There is no h anywhere in the Gaudin model. The twist is a Lie-algebra element g=diag(ζ₁,ζ₂,ζ₃) and the Lax matrix is

```text
L_g(u)ᵢⱼ=ζᵢ δᵢⱼ+Σₛ₌₁²(E⁽ˢ⁾ᵢⱼ)/(u-θₛ).        (2)
```

The commuting family is generated by the column-ordered quantum determinant of ∂ᵤ-L_g(u), with the derivative moved to the right (∂f=f∂+f′):

```text
cdet(∂ᵤ-L_g)=Σ_(σ∈S₃)sgn σ  M_(σ(1)1)M_(σ(2)2)M_(σ(3)3)=∂³-S₁(u) ∂²+S₂(u) ∂-S₃(u),   M=∂ᵤ1-L_g.        (3)
```

Writing lᵢⱼ=L_g(u)ᵢⱼ and l′ᵢⱼ=∂ᵤlᵢⱼ, l″ᵢⱼ=∂ᵤ²lᵢⱼ, and keeping operator products in the order written, the coefficients are

```text
S₁ =l₁₁+l₂₂+l₃₃,
S₂ =Σ_(i<j)(lᵢᵢlⱼⱼ-lⱼᵢlᵢⱼ)-l′₂₂-2 l′₃₃,
S₃ =l″₃₃-l₁₁l′₃₃-l₂₂l′₃₃-l′₂₂l₃₃+l₃₁l′₁₃+l₃₂l′₂₃+l′₃₂l₂₃
+l₁₁l₂₂l₃₃-l₁₁l₃₂l₂₃-l₂₁l₁₂l₃₃+l₂₁l₃₂l₁₃+l₃₁l₁₂l₂₃-l₃₁l₂₂l₁₃.        (4)
```

We checked that [Sₐ(u),S_b(v)]=0 exactly (exact rationals, 3×3 grids of rational u,v, for [2,1,0], [3,1,0], [3,2,0]), that S₁ is a scalar and that the weights are conserved.[^3]

### 1.2. Laurent data

Every non-scalar information about a state is contained in the Laurent coefficients at the two sites:

```text
S₁ =σ+|Λ|Σₛ1/(u-θₛ),   σ=ζ₁+ζ₂+ζ₃,
S₂ =e₂(ζ)+Σₛ[Hₛ/(u-θₛ)+c₂/((u-θₛ)²) ],
S₃ =e₃(ζ)+Σₛ[bₛ/(u-θₛ)+aₛ/((u-θₛ)²) ],        (5)
```

where eₖ(ζ) are elementary symmetric polynomials, |Λ|=λ₁+λ₂+λ₃, and

```text
c₂=λ₂(λ₁+1)        (6)
```

is a scalar. S₃ has no triple pole. On an eigenstate the six operators H₁,H₂,a₁,a₂,b₁,b₂ become six numbers; this is all the state-dependent data that the Gaudin Hamiltonians carry. They are not independent: three no-log conditions per site (Section 2.3) cut the solution set down to a finite one.

The value c₂=λ₂(λ₁+1) deserves a remark, because we got it wrong twice. It equals |Λ| for [2,1,0] and [3,1,0], so those two representations cannot tell the two apart; the first guess, half the quadratic Casimir, gives the right 3 for [2,1,0] but 6 for [3,1,0], whose true value is 4. The scan over representations is in Table 2. For one-row representations (λ₂=0) the double pole is the scalar 0.[^4]

**Table 2.** The scalar double pole c₂ of S₂ and the induced exponents. Dual exponents at θₛ: {0,λ₂+1,λ₁+2}; direct exponents {0,λ₁-λ₂+1,λ₁+2}. Sources: P3 (`RESULTS-Qside-310.md`), P4 (`RESULTS-Qside-320.md`), W2/W4 (`Gaudin-L2-chain-2.wb`, `Gaudin-untwisted-exact.wb`).

| Λ | c₂=λ₂(λ₁+1) | ∣Λ∣ | dual exponents | direct exponents |
|---|---|---|---|---|
| [1,1,0] | 2 | 2 | {0,2,3} | {0,1,3} |
| [2,1,0] | 3 | 3 | {0,2,4} | {0,2,4} |
| [3,1,0] | 4 | 4 | {0,2,5} | {0,3,5} |
| [2,2,0] | 6 | 4 | {0,3,4} | {0,1,4} |
| [3,2,0] | 8 | 5 | {0,3,5} | {0,2,5} |
| [l,0,0] | 0 | l | {0,1,l+2} | {0,l+1,l+2} |

The first Laurent coefficient has an explicit form, and its sign cost us a round-1 correction. With C_ζ⁽ˢ⁾=ΣᵢζᵢE⁽ˢ⁾ᵢᵢ and W=Σₖζₖ(E⁽¹⁾+E⁽²⁾)ₖₖ,

```text
H₁=-[ΣᵢζᵢE⁽¹⁾ᵢᵢ+1/θ₁₂ Σ_(i,k)E⁽¹⁾ᵢₖE⁽²⁾ₖᵢ]+|Λ|σ+(|Λ|²)/θ₁₂,   H₁+H₂=-W+2|Λ|σ,        (7)
```

and H₂ is the same with 1↔2 and θ₁₂→-θ₁₂. The sign is - in front of the bracket, not + as the first plan had it: it follows from the residue at θ₁ of Σ_(i<j)(lᵢᵢlⱼⱼ-lⱼᵢlᵢⱼ).[^5] The same formula at ζ=0 gives H₁=(|Λ|²-½[C₂(μ)-2C₂(Λ)])/θ₁₂ on a state of total weight μ (Section 11.1).

### 1.3. The dictionary to XXX

The Gaudin family is the limit of the XXX transfer matrices at twist z=e^(hζ). Expanding t₁(u)/ν(u)=3+hA₁+h²A₂+h³A₃+… and t₂(u)/(ν(u)ν(u-h))=3+hB₁+h²B₂+… with ν(u)=Π(u)=(u-θ₁)(u-θ₂) one finds, exactly,

```text
A₁=S₁,   B₁=2S₁,   B₂=A₂+½ S₁²+½ S₁′,
S₂=-A₂+½ S₁²+½ S₁′,   S₃=A₃-B₃+⅙ S₁³+½ S₁S₁′+½ S₁″,        (8)
```

with Sₖ(u→∞)=eₖ(ζ). The relations A₁=S₁ and B₁=2S₁ are forced by the Baxter equations at orders h¹ and h²; nothing in the dictionary needs more than S₁ and its derivatives.[^6] This is why the cdet convention could be pinned independently of the chain-side commutation test: the plan's convention is the only one that fits.

### 1.4. Frames, vacua and the reference covector

Eigenstates Ψ_B of the diagonal-twist Hamiltonians live in what we call the *diagonal frame*. The SoV construction naturally lives in the *companion frame*, where the twist is replaced by its Lie-algebra companion matrix, built from eₖ=eₖ(ζ),

```text
g_c(ζ)=[[e₁, -e₂, e₃], [1, 0, 0], [0, 1, 0]],   P_aux(ζ) with rows (1, ζᵢ-e₁, ζᵢ²-e₁ζᵢ+e₂),        (9)
```

so that P_auxg_cP_aux⁻¹=diag(ζ). The physical-space image of the frame change is P^G_phys(ζ)=ρ((P_aux(ζ)⁻¹)ᵀ), and the reference covector is x₀^G=x(ζ)⊗x(ζ) with x(ζ)=e₁ᵀ(P^G_phys)⁻¹ (one site). For the adjoint, in notebook order,

```text
x(ζ)=(ζ₂-ζ₁, ζ₃-ζ₁, ζ₂-ζ₁, 2ζ₃-ζ₁-ζ₂, 3(ζ₂-ζ₁), 3(ζ₃-ζ₁), 2(ζ₃-ζ₂), 3(ζ₃-ζ₂)),        (10)
```

a covector whose components are homogeneous of degree λ₂ in ζ for every representation we have looked at.[^7] We expected x₀^G not to be parallel to the h→0 limit of the XXX x₀; it is, to 10⁻³¹ (norms ∼h per site).

The twist of the bra vacuum is conjugated into the diagonal frame, g′=P_aux(ζ) g_c(ω) P_aux(ζ)⁻¹, where ω=(ω₁,ω₂,ω₃) are the three eigenvalues of the vacuum's own twist. The bra vacuum is

```text
Ω_(ζω)=e₁ᵀρ(n),   n=exp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃),        (11)
```

```text
r₁=(ω₁-ζ₁)/(ω₁-ζ₂),   r₂=((ω₁-ζ₁)(ζ₁-ω₂)(ζ₂-ζ₃))/((ζ₁-ζ₂)(ω₁-ζ₃)(ω₂-ζ₃)),   r₃=((ω₁-ζ₂)(ζ₂-ω₂)(ζ₁-ζ₃))/((ζ₁-ζ₂)(ω₁-ζ₃)(ζ₃-ω₂)).        (12)
```

We verified that Ω_(ζω)⊗Ω_(ζω) is an exact left eigencovector of all coefficient operators of S₂,S₃ built with the twist g′, for [2,1,0], [3,1,0], [3,2,0] at three twist points; the four other orientations of the conjugation, the transposes, and the naive control diag(ω) all fail.[^8] The covector Ω_(ζω) is the vacuum of a *different* twist than the ket's; its overlaps with the ζ-eigenstates are the objects the SoV machinery turns into Q-functions.

At ω₂→∞ the vacuum has a finite limit Ω^∞[ζ,ω₁]=e₁ᵀρ(n^∞) with

```text
r₁=(ω₁-ζ₁)/(ω₁-ζ₂),   r₂=-((ω₁-ζ₁)(ζ₂-ζ₃))/((ζ₁-ζ₂)(ω₁-ζ₃)),   r₃=((ω₁-ζ₂)(ζ₁-ζ₃))/((ζ₁-ζ₂)(ω₁-ζ₃)),        (13)
```

and it depends on only two parameters in the two-site family Φ_∞=Ω^∞[ζ,ω₁]⊗Ω^∞[ζ,ω₁′]. In the companion frame the same vacuum has Gauss parameters (1/ω₁, -1/(ω₁ω₂), (ω₁+ω₂)/(ω₁ω₂)), i.e. exactly the XXX ones with w→ω, and at ω₂=∞ it is Ω_c^∞[ω₁]=e₁ᵀρ(e^(aE₁₂)e^(aE₂₃)) with a=1/ω₁. The two frames are tied by

```text
Ω_D[ζ,ω] P^G_phys(ζ)=κ^G(ζ,ω) Ω_c[ω],   κ^G=Δ₁^(λ₁-λ₂)Δ₂^(λ₂-λ₃)Δ₃^(λ₃),        (14)
```

with Δ₁=ω₁²/((ζ₂-ω₁)(ζ₃-ω₁)) and, at ω₂=∞, Δ₂=ω₁/((ζ₁-ζ₂)(ζ₃-ω₁)). This gives the scalar that appears throughout the kernel normalisation,

```text
c=((ω₁ω₁′)^(2λ₁-λ₂))/(κ^G(ζ,ω₁) κ^G(ζ,ω₁′))=(ζ₂-ω₁)^(λ₁-λ₂)(ζ₃-ω₁)^(λ₁)(ζ₂-ω₁′)^(λ₁-λ₂)(ζ₃-ω₁′)^(λ₁)(ζ₁-ζ₂)^(2λ₂).        (15)
```

It was derived for the adjoint (c=(ζ₂-ω₁)(ζ₃-ω₁)²(…)(ζ₁-ζ₂)²), [3,1,0] and [3,2,0], with the other exponent choices excluded by the fits.[^9]

### 1.5. The XXX scaffold in one paragraph

All of the above was first understood in the XXX chain with L=2 sites in the new conventions (Lax uδᵢⱼ+hEᵢⱼ, companion twist on the right, B and C built from the untwisted T). There, the SoV basis consists of the eigencovectors of B(u), with simple spectrum and labels (a|bc) per site, and the overlap of the vacuum Ω_(zw)=e₁ᵀρ(n(z,w)) with an eigenstate is a determinant or a sum of products of ladder brackets of the Q-functions. For the family of two rotated vacua Ω_w⊗Ω_(w′) the untwisting z,w,w′→1 was clean: K/S=O(ε) on every HW branch, all CG data at leading order, and for the adjoint (∞,∞) family the Q-number is a sum of ten products of 2×2 bracket determinants.[^10] The XXX program also needed, for the single-vacuum approach, a “cleaning” step for some sectors (`Clean-CG-combinations.md`). The Gaudin model is the h→0 limit, equivalently θ₁₂/h→∞ at fixed twist, because of the scaling symmetry (u-θ)→(u-θ)/λ, g→λg. The XXX companion matrix is not the exponential of a Lie-algebra companion matrix (at z=1 it is a Jordan block), which is why the Gaudin companion frame had to be defined separately; the Jordan block turned out to be harmless, since (G′-1)/h→g′ at O(h).

### 1.6. Representations, bases and numerical points

We label the GT basis of V_Λ in *notebook order*: index 1 is the highest weight, patterns sorted in descending order. The Python modules use the *`gl3_reps` order*, which starts at the lowest-weight vector, so that the HW ket has the last index (7 for the adjoint, 14 for [3,1,0] and [3,2,0]). The map is a permutation, not a reversal: for the adjoint, notebook index a=0,…,7 corresponds to `gl3_reps` index (7,5,6,4,2,1,3,0).[^11] The two-site product basis has index (i,j)↦(i-1)dimV+j, site 1 slow, site 2 fast.

The numerical points reused across all agents are these. Inhomogeneities θ=(2/7,-3/5), θ₁₂=31/35. Twist points P: ζ=(3/2,5/4,2/3) and P′: ζ=(2,5/3,-1/2), with extra points (7/5,-3/4,1/6), (-4/3,2/9,5/7), (1/3,-2/7,5/11). Vacuum twists ω=(5/2,7/3,3/5), ω′=(7/4,-2/5,9/7) and the two-parameter pair (ω₁,ω₁′)=(5/2,7/4). The untwisting regime is ζ=εφ, ω=εω̂, ω′=εω̂′ with φ=(1/3,-2/7,5/11) (control direction φ′=(-2/5,3/8,7/11)), and ε∈{1,2,4}·10⁻⁶ for the chain-side checks and εₖ=10⁻⁶2⁻ᵏ, k=0,…,11, for the series fits; precision is 200 digits on the chain side and 60 to 120 digits on the Q side. The full list of grids is in Appendix A.

## 2. Q-functions, Wronskians, detN

**Punchline.** Every eigenstate is a triple of polynomial-times-exponential solutions of one third-order ODE (the oper built from S₁,S₂,S₃); the ODE has local exponents {0,λ₂+1,λ₁+2} at each θₛ, which forces linear relations among the jets of the solutions and fixes the Frobenius constant κₛ=aₛ/c₂. The overlap of two states is detected by a 6×6 determinant of residue brackets whose columns come from a Lagrange identity.

### 2.1. The oper and its Q-functions

The dual oper is

```text
(∂³-S₁∂²+S₂∂-S₃)P=0,   Pʲ=e^(ζⱼu) pⱼ(u),  degpⱼ=μⱼ, pⱼ monic,        (16)
```

where μ is the weight of the state. For each of the three exponentials there is exactly one polynomial solution (kernel dimension 1), and Σdeg=|μ|. The XXX dual Baxter equation becomes exactly this oper at order h³, with Qʲ=zⱼ^(u/h)qʲ→e^(ζⱼu)pⱼ.[^12]

The direct Q-functions are a little subtler, and this was our first surprise. One expects them to solve the formal adjoint of (16). They do not: the plain adjoint has *no* solutions of the form e^(-ζu)×polynomial. What is true is that ψᵢ=Π⁻ᵖQᵢ solves the adjoint, with Π(u)=(u-θ₁)(u-θ₂),

```text
Qᵢ=e^(-ζᵢu)qᵢ,   degqᵢ=2λ₁-μᵢ,   direct oper=Πᵖ∘(formal adjoint)∘Π⁻ᵖ,  p=λ₁-λ₃.        (17)
```

The XXX Baxter equation at order h³ gives the adjoint acting on Qᵢ/ν², which is the same statement for the adjoint. For [2,1,0] the gauge power 2 is forced by the degrees; for [3,1,0] and [3,2,0] a scan over p=0,…,5 gives no solution at p=0,1,2, exactly the right degrees 6-μᵢ at p=3 in all 36 sectors, and degrees raised by 2(p-3) above. All our data have λ₃=0, so p=λ₁-λ₃ and p=λ₁ are not separated by data; the exponent argument below selects the former.[^13] Table 3 lists the degrees for the three representations at the top state.

**Table 3.** Degrees of the Q-functions. μ is the state weight; dual degrees are μⱼ, direct degrees 2λ₁-μᵢ. The vacuum has μ=(2λ₁,2λ₂,0). Sources: P1, P3, P4 (all 64/225/225 states: degrees as stated in every sector, kernel dimension 1).

| Λ | sectors | dual degpⱼ (vacuum) | direct degqᵢ (vacuum) | gauge power p |
|---|---|---|---|---|
| [2,1,0] | 19 | (4,2,0) | (0,2,4) | 2 |
| [3,1,0] | 36 | (6,2,0) | (0,4,6) | 3 |
| [3,2,0] | 36 | (6,4,0) | (0,2,6) | 3 |

### 2.2. Local exponents and Frobenius constants

At θₛ the operator has a regular singular point with indicial polynomial e(e-1)(e-2)-|Λ|e(e-1)+c₂e, whose roots are

```text
{0, λ₂+1, λ₁+2}  (dual),  {0, λ₁-λ₂+1, λ₁+2}  (direct).        (18)
```

This was confirmed in 450/450 (state, site) cases for each of [2,1,0], [3,1,0], [3,2,0] from the vanishing orders of the span of the three solutions; for [3,2,0] the direct exponents are {0,2,5}, not {0,3,5}, and duality swaps the two patterns.[^14] A consequence of the first exponent being λ₂+1 is that the jets of the dual Q-functions at θₛ obey linear relations. For every solution P of the oper,

```text
P′(θₛ)=κₛP(θₛ),   κₛ=aₛ/c₂,        (19)
```

with the same κₛ for all three solutions, because the coefficient of (u-θₛ)⁻² in (∂³-S₁∂²+S₂∂-S₃)P reads (I(1)c₁-aₛc₀) with I(1)=c₂, and no regular part enters. For λ₂≥2 the second derivative is also tied, P″(θₛ)=κ_(2,s)P(θₛ) with κ_(2,s)=(aₛ²-Hₛaₛ+8bₛ)/24 for [3,2,0]. We checked the relations on all states and both sites, to 2.7·10⁻⁷⁸ (κ, adjoint, 80 digits) and 3.1·10⁻⁵⁷, 3.6·10⁻⁵⁷ (κ,κ₂, [3,2,0]).[^15]

The jet matrix [Pʲ⁽ᵏ⁾(θₛ)/k!]_(j,k) then has rank 1 for k≤λ₂ and rank 2 up to k=λ₁+1, jumping to 3 at order λ₁+2 (measured for [3,2,0]: 3×3 rank 1, 3×5 rank 2, 3×6 rank 3, in 450/450 cases). Hence the Plücker coordinates

```text
pₖₗ=P²⁽ᵏ⁾P³⁽ˡ⁾-P³⁽ᵏ⁾P²⁽ˡ⁾  (at θₛ)        (20)
```

vanish for 0≤k<l≤λ₂ (“forced” pairs), and each non-forced pₖₗ is the 2×2 determinant Dₛ=c₂d₃-c₃d₂ of the two non-trivial solutions' leading components times a universal function of the Frobenius coefficients. For the adjoint p₀₁=W(θₛ)=0 to 1.9·10⁻⁷⁸.

### 2.3. The no-log conditions: Gaudin's Bethe equations

The three polynomial solutions must be free of logarithms at θₛ. Resonances between the exponents 0,λ₂+1,λ₁+2 give three independent conditions per site, six in all, on the six numbers (Hₛ,aₛ,bₛ). For λ₂=1 the first of them is

```text
aₛ²-Hₛaₛ+c₂ bₛ=0,  i.e.   bₛ=(aₛ(Hₛ-aₛ))/c₂,        (21)
```

(c₂=3 for the adjoint, 4 for [3,1,0]); the other two involve the Taylor coefficients of the regular parts of S₁,S₂,S₃ and therefore the other site's data. For [3,2,0] the first condition arises at order x³ and is cubic. All hold on every state and both sites to 10⁻⁵⁶ to 10⁻⁶⁰ ([3,1,0]: 1.3·10⁻⁵⁸; [3,2,0]: 3·10⁻⁵⁹).[^16] We call these the Gaudin Bethe equations in Laurent-data form. In the untwisted model they are what constrains the eigenvalue data in the Casimir language of Section 11.1, and perturbatively in ε they determine the series of Section 10 (except in two degenerate cases).

### 2.4. Wronskians and quantum Weyl

The three solutions are linked by Wronskian (QQ) relations with a representation-dependent power of Π:

```text
PᵃPᵇ′-PᵇPᵃ′ =(ζ_b-ζₐ) e^(σu) Π^(λ₂) Q_c,
QₐQ_b′-Q_bQₐ′ =(ζₐ-ζ_b) e^(-σu) Π^(λ₁-λ₂) Pᶜ,        (22)
```

for (a,b,c) a permutation of (1,2,3) and monic pⱼ,qᵢ, so that the constants are exactly ζ_b-ζₐ and ζₐ-ζ_b. The triple Wronskian is W=det[P,P′,P″]=Π_(i<j)(ζⱼ-ζᵢ) e^(σu) Π^(|Λ|), state-independent. Relative remainders of the polynomial divisions are ≤10⁻⁴⁸ to 10⁻⁵⁴ on all states.[^17] The rule is F=Π^(λ₁-λ₂) for the direct side, F′=Π^(λ₂)=Π^(|Λ|-p) for the dual side. In particular W(P²,P³)(θₛ)=0, and W′(θₛ)=(ζ₃-ζ₂)e^(σθₛ)Π′(θₛ)Q₁(θₛ) for the adjoint.

The oper itself is recovered from the three solutions by the quantum-Weyl formulas

```text
S₁=(det[P,P′,P″′])/W,   S₂=-(det[P,P″′,P″])/W,   S₃=(det[P″′,P′,P″])/W,        (23)
```

(rows labelled by j), verified on all states to 4·10⁻⁴⁷ ([3,1,0]) and 1.8·10⁻⁴⁵ ([3,2,0]).

### 2.5. Orthogonality and the detN overlap

The SoV scalar product of an eigenstate with a vacuum is detected by a 6×6 determinant, and this is how we first saw the Q-structure on the Gaudin side. The orthogonality relation is not of fixed-coefficient type: it follows from the Lagrange identity ψ D_B[P]-P D_B^†[ψ]=C(P,ψ)′ with P=Pʲ of ket B and ψ=Π⁻ᵖQᵢ of bra A. Because D_B^†ψ=ΔD^†ψ=-(ΔS₂ψ)′-ΔS₃ψ (with ΔSₖ=Sₖ^B-Sₖ^A; S₁ and all constants cancel), taking the residue at each θₛ (the residue of a derivative vanishes) gives

```text
ΣₘΔₘ N[(s,i,j),m]=0,   Δ=(ΔH₁,ΔH₂,Δa₁,Δa₂,Δb₁,Δb₂),        (24)
```

```text
N[·,Hₜ]=-Res_(θₛ)ψP′/(u-θₜ),   N[·,aₜ]=Res_(θₛ)ψP/((u-θₜ)²),   N[·,bₜ]=Res_(θₛ)ψP/(u-θₜ),        (25)
```

t=1,2. The six non-scalar Laurent data differences Δₘ are state-dependent, so no fixed linear combination of brackets works; that is why the spec's search for fixed coefficients could never have succeeded. The identity holds to 6·10⁻⁵⁷ ([2,1,0], 300 pairs), 5.5·10⁻⁵⁶ ([3,1,0], 600 pairs), 8.4·10⁻⁵⁷ ([3,2,0], 600 pairs).

The stacked 18×6 matrix (rows (s,i,j) for the nine pairs (i,j), two sites) has rank 5 for generic A≠B, never 6; it has rank 6 for A=B. The 6×6 determinant detN, built on the rows (s,(i,j)) with (i,j)∈{(1,3),(1,2),(2,3)}, vanishes for A≠B at the same twist (≤1.6·10⁻⁶⁵ over 100 pairs, 2.2·10⁻⁶² for [3,2,0]) and is non-zero for A=B (minimum 2·10⁻³⁷). For a vacuum of twist ω, with R_B(ω)=⟨Ω_(ζω)⊗Ω_(ζω)|Ψ_B⟩ for *any* normalisation of Ψ_B and D_B(ω) the detN built from the ket's dual Pʲ and the vacuum's direct Qᵢ^vac,ω,

```text
(R_B(ω))/(R_B(ω′))=(D_B(ω))/(D_B(ω′)) ·(D_(B₀)(ω′))/(D_(B₀)(ω)),   B₀=e₁⊗e₁,        (26)
```

on all 64 (6.4·10⁻⁵⁵), 225 (6.8·10⁻⁵²) and 225 (1.6·10⁻⁵¹) states of [2,1,0], [3,1,0], [3,2,0]; the identity is x₀-free. In the untwisting regime it holds on the HW branches to 5.4·10⁻⁶⁷, 6.6·10⁻⁵², 2.2·10⁻⁵² at ε=10⁻⁶ with |detN|∼ε¹² for (4,2,0),(4,1,1),(3,3,0), ε¹¹ for the octets and ε¹⁰ for the singlet in [2,1,0]; for [3,1,0], ε¹⁶ on (6,2,0),(6,1,1),(5,3,0),(4,4,0), ε¹⁵ on both (5,2,1) and both (4,3,1), ε¹⁴ on (4,2,2),(3,3,2), and the same pattern transported by duality for [3,2,0]. Of the 84 triples of pairs, 72 work, 6 are degenerate (all i or all j equal) and the 6 permutation triples fail; the XXX triple works. Controls (a₂ column replaced by a triple pole, swapped orientation, wrong gauge Π⁻²,Π⁻⁴) fail at deviation ∼1 to 2.[^18] Why the same six equal-twist columns also close the mixed-twist overlaps, with no ΔS₁,Δe₂,Δe₃ columns, was observed and not derived.

### 2.6. The rank-4 rule

About 20% of different-sector pairs have rank 4 instead of 5. The rule is clean: the rank is a function of Δμ=μ_B-μ_A only, and **rank 4 if and only if Δμ has exactly one positive and two negative components** (all 4032 ordered pairs of [2,1,0]: 3318 of rank 5, 714 of rank 4; 133 rank-4 pairs among 600 for [3,1,0], 130 for [3,2,0]). The mechanism: for every k with Δμₖ≤-1 the product ψₖ^AP_Bᵏ is a rational function, its exponentials cancelling, and it decays at infinity, so its residues at θ₁ and θ₂ sum to zero: row(1,k,k)+row(2,k,k)=0 in all six columns (violation ≤5.4·10⁻⁵⁸). Each site's nine rows span a 3-dimensional space Sₛ, the rank is 6-dim(S₁∩S₂), each negative component of Δμ adds one vector to S₁∩S₂, and two negative components give rank 4. The second null vector has no closed form (open, minor).[^19]

### 2.7. Vacuum Q-functions and vacuum Laurent data

The vacuum covector is itself an eigencovector, with Q-functions of twist ω: direct Qᵢ^vac=e^(-ωᵢu)qᵢ of degrees 2λ₁-μ^vacᵢ as in Table 3, with q₁=1 (the vacuum has μ^vac=(2λ₁,2λ₂,0)). Its Laurent data are known in closed form and independent of ζ:

```text
Hₛ^vac=Σᵢ(|Λ|-λᵢ) ωᵢ ± (2e₂(Λ))/θ₁₂,   aₛ^vac=c₂ ω₃,   bₛ^vac=ω₃(Hₛ^vac-|Λ|ω₃),        (27)
```

with the upper sign for s=1, e₂(Λ)=λ₁λ₂+λ₁λ₃+λ₂λ₃. The middle statement is where |Λ| and c₂ were conflated before [3,2,0]: a^vacₛ=8ω₃ for [3,2,0], not 5ω₃. The consistency check is κ^vac=a^vac/c₂=ω₃, the exponent of P³. The vacuum eigenvalue functions for the three representations are tabulated in Appendix A.[^20]

## 3. Which operators localise: the answer to “which T's?”

**Punchline.** In Gaudin, every Laurent coefficient of S₂,S₃ at θ₁ localises the two-vacua covector on site 1 by itself: {1,H₁,S3c[1,2]} each separately for every representation, plus one extra quadratic combination X₄ once λ₁=3. The number of localisers is ℓ₁=λ₁+1, as in XXX, but the fused shapes of XXX are not needed.

### 3.1. The question

Paul's question at the start of the project was “I don't know which T's to use”. In XXX the answer was specific: the shapes T_((1,1))(θ₁-2h) and T_((2,1))(θ₁-2h) (for the adjoint), and in general T_((k,k)) (k≤d) and T_((k,d)) (d<k≤λ₁) at θ₁-λ₁h with d=λ₁-λ₂, which acting on ⟨Ω⊗Ω| produce f_λ⟨ΩT⁽¹⁾_λ|⊗⟨Ω|: a covector on site 1 times the unchanged vacuum on site 2. The set of such site-1 covectors is the “localising span” L₁, and the CG data are read off from L₁ as soon as it reaches every weight with rank equal to the multiplicity.

In Gaudin there are no shifts θₛ-mh; the commuting family is the Laurent coefficients themselves. We test every operator X in the algebra they generate: X *localises* if ⟨Ω⊗Ω|X=v⊗Ω₁ for a covector v on site 1 (the rank-1 scan: ⟨Ω⊗Ω|X viewed as an dimV×dimV matrix has rank 1 and second factor Ω).

### 3.2. The scan

Let S3c[s,k] denote the coefficient of (u-θₛ)⁻ᵏ in S₃ (so S3c[s,2]=aₛ and S3c[s,1]=bₛ), and

```text
D₁=span{1, W, H₁, S3c[1,2], S3c[1,1], S3c[2,2], S3c[2,1]},   D₂=D₁+quadratic monomials,        (28)
```

with dimD₁=7 (independent basis), dimD₂=24 for [2,1,0] and 26 for [3,1,0], [3,2,0]; the cubic span has rank 48,68,68 (not 64, hence not vacuous). The result, exact at P and P′, is as follows.

**Adjoint.** ℓ₁=ℓ₂=3 already inside D₁, with no products needed. Site 1: {1,H₁,S3c[1,2]} localise, each separately,

```text
⟨Ω⊗Ω|H₁=[Ω₁(c_H-C_ζ)]⊗Ω₁,   ⟨Ω⊗Ω|S3c[1,2]=[Ω₁S3^1-site[-2]]⊗Ω₁,        (29)
```

where C_ζ=ΣₖζₖEₖₖ, S3^1-site[-2] is the (u-θ₁)⁻² coefficient of the one-site cdet with twist diag ζ, and c_H is a constant (below). Hence the localising span is L₁=Ω₁·span{1,C_ζ,S3^1-site[-2]}, of full rank on every weight of V. Site 2: {1, W+H₁(=-H₂+scalar), S3c[2,2]} with +C_ζ.[^21]

**[3,1,0] and [3,2,0].** ℓ₁=4=λ₁+1: three over D₁ (the same three), four over D₂ (and over D₁+H₁D₁ at site 1). No single quadratic monomial supplies the fourth; the minimal supports are {H₁²,H₁S3c[1,2]} and {H₁²,S3c[1,2]²}, and the localiser is

```text
X₄=H₁²-½ H₁ S3c[1,2]-γ/θ₁₂  W-S3c[1,1]  (mod 1,H₁,S3c[1,2]),        (30)
```

where γ=2 for [3,1,0] and γ=1 for [3,2,0], and with ζ- and ω-independent coefficients (checked at 9 values of θ₂ at P and 2 at P′ for [3,2,0]). Its covector v₄=Ω₁(c₀+c₁C_ζ+c₂S3[-2]-S3[-1]+C_ζ²+½ C_ζ S3[-2]) is θ₁₂-independent. The coefficient γ changed from 2 to 1; it could be λ₁-λ₂ or 2/λ₂, and two data points cannot decide. The three-dimensional L₁ already reaches every weight (multiplicities ≤2), so the fourth is needed to match the count ℓ=λ₁+1 and not for reach.[^22]

The quadratic subspace has exact relations that make D₂ smaller than 28 by two. For [3,1,0] they are S3c[s,1]=¼ S3c[s,2](Hₛ-S3c[s,2]) for s=1,2. For [3,2,0] there is no scalar c_Λ (none of 1/4,1/3,1/5,1/6,1/8,1/(λ₁+1),1/|Λ|) and the relation is an inhomogeneous quadratic: with a=S3c[1,2], s=Σζ,

```text
S3c[1,1]=3H₁²-7/4 H₁a+¼ a²-3/θ₁₂ W-βH₁-γ′a-δ,        (31)
```

β=78/θ₁₂+15s, γ′=-45/(2θ₁₂)-9/2 s, δ=-504/θ₁₂²-210s/θ₁₂-18s²-3e₂(ζ), an exact 225×225 operator identity at six combined (ζ,θ) points; site 2 is the mirror θ₁₂→-θ₁₂, H₁→H₂=-H₁-W+10s. Its quadratic part factorises, 12H²-7Ha+a²=(3H-a)(4H-a).[^23]

### 3.3. Spectators and closed forms

Two things accompany each localiser: the constant c_H in (29) and the *spectator* factor f(X), defined by ⟨Ω⊗Ω|X|e₁⊗e₁⟩ with v normalised to first component 1 (the analogue of f_λ in XXX). All of them have closed forms:

```text
c_H=|Λ|σ+(2e₂(Λ))/θ₁₂,   f(H₁)=|Λ|σ-ζ·Λ+(2e₂(Λ))/θ₁₂,   f(S3c[s,2])=c₂ ζ₃,        (32)
```

with f(H₂) having the opposite sign of the 1/θ₁₂ term, and e₂(Λ)=λ₁λ₂+λ₁λ₃+λ₂λ₃. The constant c_H has a one-line origin: contracting the Ω's one has Σᵢₖ(ΩEᵢₖ)⊗(ΩEₖᵢ)=(Σᵢλᵢ²) Ω⊗Ω (the Casimir-type identity on the vacuum), and with |Λ|²-Σλᵢ²=2e₂(Λ) the 1/θ₁₂ coefficient follows. We initially guessed 2λ₁ for it (from the two data points 4,6 of [2,1,0],[3,1,0]) and were wrong: it is 2e₂(Λ), equal to 2λ₁λ₂=4,6,12 for [2,1,0],[3,1,0],[3,2,0]. Table 4 lists spectators at the two standard points.[^24]

The last formula in (32) is the place where c₂=λ₂(λ₁+1) and not |Λ| appears: f(S3c[s,2])=8ζ₃ for [3,2,0], not 5ζ₃; the same correction as for the vacuum's a^vacₛ=c₂ω₃ (Section 2.7), and for the same reason, since ⟨Ω⊗Ω|S3c[1,2]=Ω₁S3^1-site[-2]⊗Ω₁ and the one-site double pole is c₂ζ₃.

**Table 4.** Spectator factors at the two standard twist points (P=(3/2,5/4,2/3), P′=(2,5/3,-1/2)). Sources: W1, W2, W3 summary cells.

| Λ | f(H₁) at P | f(H₁) at P′ | f(S3c[1,2]) at P | at P′ |
|---|---|---|---|---|
| [2,1,0] | 6+4/θ₁₂ | 23/6+4/θ₁₂ | 2 | -3/2 |
| [3,1,0] | 95/12+6/θ₁₂ | 5+6/θ₁₂ | 8/3 | -2 |
| [3,2,0] | 121/12+12/θ₁₂ | 13/2+12/θ₁₂ | 16/3 | -4 |

### 3.4. Why Gaudin differs from XXX

In XXX only the specific fused shapes T_((1,1)),T_((2,1)),… localise; in Gaudin every Laurent coefficient in D₁ at θ₁ localises on its own. This was the round-1 surprise, and we have no derivation of it. Our reading, which has not been tested, is that in the h→0 limit the shifted transfer matrices at θₛ-mh collapse onto jets at θₛ, so that the products needed for a fused shape reduce to combinations of Laurent coefficients, and the spectators become Laurent polynomials in θ₁₂ with point-dependent constants instead of ladder products.

## 4. The two-vacua family and the untwist: the chain side

**Punchline.** The covector Ω_(ζω)⊗Ω_(ζω′) with ω≠ω′, or its two-parameter (∞,∞) sub-family, is CG-complete: after an exact undo it becomes ⟨φᵤ⊗e₁| with ε-free Gauss parameters, it has K/S=O(ε) on every HW branch of [2,1,0], [3,1,0], [3,2,0], and every CG coefficient is read off at leading order. No cleaning step is needed anywhere.

### 4.1. Why a two-vacua family

With a *single* vacuum and localised transfer matrices the XXX approach had to be cleaned: the localising covector is pure top weight at ε=0, its O(ε) lower-weight part meets the O(ε⁻ᵏ) raised corrections of the HW-limit eigenstate, and the contamination ratio comes out K/S=κ-1=O(1) (for example κ=(t+1)/2 for [3,3,0]). “Sector-wide cleaning” removed K for a whole sector by one combination of insertions; “per-copy cleaning” was needed where the clean space lacked HW rank ([4,3,1] in [3,1,0]^(⊗2), [5,3,2] in [3,2,0]^(⊗2), [6,4,2] in [4,2,0]^(⊗2)).[^25] The two-vacua family removed the need already in XXX (Stage D). Paul asked whether Gaudin avoids cleaning too; the answer is yes, and it is the content of this section.

### 4.2. The family, the exact undo and CG-completeness

The family is Φ=Ω_D[ζ,ω]⊗Ω_D[ζ,ω′] in general and Φ_∞(ω₁,ω₁′) for ω₂,ω₂′→∞. On every HW ket h of the untwisted algebra the product of two rotated vacua reduces to a rotation on site 1 only, because ρ(n)⊗ρ(n) with n unipotent acts trivially on HW kets:

```text
⟨Φ|h⟩=⟨φᵤ⊗e₁ᵀ|h⟩,   φᵤ=e₁ᵀρ(n′),   n′=e^(a′E₁₂)e^(g′E₁₃)e^(c′E₂₃).        (33)
```

We checked this identity exactly on every HW ket of [2,1,0] (at 11 points), [3,1,0] and [3,2,0] (9 points). In the regime ζ=εφ, ω₁=εω̂₁, ω₁′=εω̂₁′ and for the (∞,∞) family the Gauss parameters are ε-independent,

```text
a′=((φ₁-φ₂)(ω̂₁-ω̂₁′))/((φ₂-ω̂₁)(φ₂-ω̂₁′)),  
c′=((φ₁-φ₃)(φ₂-φ₃)(ω̂₁-ω̂₁′))/((φ₁-φ₂)(φ₃-ω̂₁)(φ₃-ω̂₁′)),  
g′=((φ₃-φ₁)(φ₃-φ₂)(ω̂₁-ω̂₁′))/((φ₃-ω̂₁)(φ₂-ω̂₁′)(φ₃-ω̂₁′)),        (34)
```

equal to the XXX leading coefficients exactly (the XXX values are the same rational functions). The Jacobian of the map from the direction parameters to (a′,c′,g′) has rank 3 for the full family and rank 2 for (∞,∞), and the per-weight ranks of the family of covectors φᵤ over a grid of ω̂ equal the weight multiplicities on every weight: the family is CG-complete. The ε-independence of Ω_D[εφ,εω̂] and of φᵤ is exact (as exact rationals) and is a gift of Gaudin: in XXX the vacua were O(1) with O(ε) corrections.[^26]

### 4.3. K/S and the CG data

For each HW-limit branch Ψ_B=Ψ₀+εΨ₁+… we decompose ⟨Φ|Ψ_B⟩=S+K, with S=⟨φᵤ⊗e₁ᵀ|P_HWΨ_B⟩ the “signal” and K the contamination from the non-HW part of Ψ_B. The result is that K/S=O(ε) on every branch, with S=O(1): the maximum over the 24-point grid at ε=10⁻⁶ doubles at 2·10⁻⁶ (ratios 2.000) and K≡0 on the top branch. Table 5 lists the values; the multiplicity-2 sectors (the ones that needed cleaning in XXX) behave like all others.

**Table 5.** max|K/S| over the grid at ε=10⁻⁶ (exponent of ε exactly 1.00 everywhere, doubling ratio 2.000). Full family: 24 grid points; (∞,∞): its two-parameter sub-family. [2,1,0] and [3,1,0] from W1/W2; [3,2,0] from W3. Entries “–” were not printed separately.

| Λ | branch μ | full family | (∞,∞) |
|---|---|---|---|
| [2,1,0] | (4,1,1) | 4.9·10⁻⁷ | ≤2.5·10⁻⁶ |
|  | (3,3,0) | 2.4·10⁻⁶ | – |
|  | (3,2,1) copies 1, 2 | 2.4·10⁻⁶, 1.6·10⁻⁶ | – |
|  | (2,2,2) | 5.4·10⁻⁶ | – |
| [3,1,0] | all ten branches | ≤7.7·10⁻⁶ [(3,3,2)] | ≤4.8·10⁻⁶ |
| [3,2,0] | (6,3,1) | 2.44·10⁻⁷ | 1.27·10⁻⁶ |
|  | (6,2,2) | 9.76·10⁻⁷ | 5.09·10⁻⁶ |
|  | (5,5,0) | 2.42·10⁻⁶ | 2.42·10⁻⁶ |
|  | (5,4,1) copies 1, 2 | 2.44·10⁻⁶, 1.81·10⁻⁶ | 5.75·10⁻⁷, 2.53·10⁻⁶ |
|  | (5,3,2) copies 1, 2 | 1.79·10⁻⁶, 1.42·10⁻⁶ | 2.39·10⁻⁶, 1.30·10⁻⁶ |
|  | (4,4,2) | 5.53·10⁻⁶ | 1.37·10⁻⁶ |
|  | (4,3,3) | 4.16·10⁻⁶ | 3.65·10⁻⁶ |

The CG coefficients are extracted from the overlap alone: for each branch the overlaps at several grid points give linear equations on the unknown slice ⟨e_ν⊗e₁|Ψ₀⟩ (one equation per grid point, using the known φᵤ), solved exactly for one-dimensional slices and by 2×2 systems for two-dimensional ones. The raw error is exactly O(ε) (ratio 4.00 between 4·10⁻⁶ and 10⁻⁶), and after Richardson extrapolation in ε the extracted slices equal the exact G-norm-1 HW-ket slices to 10⁻¹⁷–10⁻²¹ (Table 6); the exact values are listed in Appendix D.1. For the adjoint the one-dimensional slices are -1/√2, -1/√2, 0.1178511302 for (4,1,1),(3,3,0),(2,2,2) and the octet slices are (0.04931160001,0.2088872897) and (0.3379867346,-0.07978784486).

**Table 6.** Precision of the chain-side CG extraction after Richardson. Sources: W1 summary (task 4d), W2 (task A), W3 (A4 cleaning verdict).

| Λ | full family | (∞,∞) |
|---|---|---|
| [2,1,0] | 10⁻¹⁷ to 10⁻¹⁹ | 10⁻¹⁸ to 10⁻²¹ |
| [3,1,0] | 10⁻¹⁷ to 10⁻¹⁹ | 10⁻¹⁸ to 10⁻²¹ |
| [3,2,0] | ≤1.5·10⁻¹⁸ (chain Richardson ≤1.6·10⁻²⁰) | ≤1.9·10⁻¹⁸; (4,4,2): 1.3·10⁻¹⁷ |

### 4.4. Multiplicity-2 sectors and universal mixing

In a multiplicity-2 sector the HW space is two-dimensional and the untwisted copies of the HW ket are fixed by the eigenvectors of S3c[1,2]=a₁ restricted to it. We expected this to depend on θ. It does not: θ₁₂ S3c[1,2]|_HW is a θ-independent matrix with no scalar part (exact at θ=(2/7,-3/5),(1/5,-4/9),(3/7,-1/4)), so the copies, and the mixing, are universal. For each sector we record the matrix N₀=θ₁₂ S3c[1,2]|_HW and the coordinate-ratio quadratic of its eigenvectors:

**Table 7.** Universality of the multiplicity-2 mixing. N₀=θ₁₂S3c[1,2]|_HW on the HW space (coordinates in the canonical rational HW basis of the data files), its trace and determinant, the eigenvalues α_±=(tr∓√disc)/2 with disc=tr²-4det, and the quadratic satisfied by the ratio r of the coordinates of the two eigenvectors. Sources: W2 task C, W3 cells 33, 35, 37; the adjoint matrix at θ₁₂=31/35 is (105/62)[[3,-1],[-1,7]].

| rep, sector | N₀ | tr, det | eigenvalues | r-quadratic |
|---|---|---|---|---|
| [2,1,0] (3,2,1) | 3/2 [[3, -1], [-1, 7]] | 15, 45 | (15∓3√5)/2 | r²+4r-1=0 |
| [3,1,0] (5,2,1) | ⅓ [[19, -40], [-2, 44]] | 21, 84 | (21∓√105)/2 | 40r²+25r-2=0 |
| [3,1,0] (4,3,1) | 4/3 [[4, -1], [-2, 14]] | 24, 96 | 12∓4√3 | r²+10r-2=0 |
| [3,2,0] (5,4,1) | [[13, -6], [-1, 22]] | 35, 280 | (35∓√105)/2 | 6r²+9r-1=0 |
| [3,2,0] (5,3,2) | 4/3 [[16, -1], [-2, 26]] | 56, 736 | 28∓4√3 | r²+10r-2=0 |

The octet mixing parameter R=(e₁E₁₃·C)/(e₁E₂₃E₁₂·C) on the weight-(1,1,1) slice satisfies R²-5R-5=0 for the adjoint, with roots 5.854101966250 and -0.854101966250; the Gaudin-extracted values agree to 1.5·10⁻¹⁸ and 5.7·10⁻²⁰. For [3,1,0] (5,2,1) the roots are 4.311737691490 and -0.811737691490 (agreement 5·10⁻¹⁹), and exactly the same numbers appear for [3,2,0] (5,4,1), as duality requires. Notice that (5,3,2), the sector that needed per-copy cleaning in XXX, has literally the same quadratic as [3,1,0]'s (4,3,1), namely r²+10r-2=0. The XXX values at finite t=θ₁₂/h were different (for example -1.694 and 1.096 for the adjoint octet); the comparison with the t→∞ limit of those numbers was proposed but not carried out. For [3,1,0] agent W1 did not compute the mixing in the (4,3,1) sector (the functionals e₁E₁₃ and e₁E₂₃E₁₂ do not suffice there); the [3,2,0] extractions of W3 solve the two-dimensional slices from pairs of grid points.

### 4.5. Duality [3,2,0]≅[3,1,0]*⊗det³

The HW sectors of [3,2,0]^(⊗2) are exactly the dual images μ*=(6-μ₃,6-μ₂,6-μ₁) of those of [3,1,0]^(⊗2), with the same multiplicities: (6,2,0)→(6,4,0), (6,1,1)→(5,5,0), (5,3,0)→(6,3,1), (5,2,1)×2→(5,4,1)×2, (4,4,0)→(6,2,2), (4,3,1)×2→(5,3,2)×2, (4,2,2)→(4,4,2), (3,3,2)→(4,3,3). The eigenvalues of θ₁₂S3c[1,2] on a dual sector differ from those of the original by a *sector-dependent scalar* while the traceless parts agree: (4,3,1)→(5,3,2): +16 (96-24x+x²→736-56x+x²); (5,2,1)→(5,4,1): +7 (84-21x+x²→280-35x+x²); multiplicity one: (3,3,2)→(4,3,3): 20→36, (4,2,2)→(4,4,2): 18→30, (4,4,0)→(6,2,2): 0→24, (5,3,0)→(6,3,1): 0→16, (6,1,1)→(5,5,0): 8→0, (6,2,0)→(6,4,0): 0→0. The scalar shift is new. The G-orthonormal slice ratios agree for (4,3,1)↔(5,3,2) (-3.6049 and 0.06935 in both) but not for (5,2,1)↔(5,4,1), because the GT bases of V and V* differ by a rotation in a two-dimensional weight space; the basis-independent mixing R above agrees. The exponents of detN, |R(ω)| and the ε-orders transport by duality too (ε^(16,15,14)).[^27]

### 4.6. The controls and the “no cleaning” answer

Three controls test what the family needs.

- ω₁′=ω₁ (common rotation). Then a′=g′=0, the family reads only the top branches: (4,2,0),(4,1,1) for [2,1,0] (full family), (6,4,0),(6,3,1),(6,2,2) for [3,2,0], and only the top branch in the (∞,∞) family. The others are “not read”. This is as in XXX.
- ω₁′ held fixed at a finite value, not untwisted. In XXX this broke the untwisting. In Gaudin with *finite* ω₂ it does not: Ω_D[ζ→0,ω′] has a finite direction-dependent limit, K/S is still O(ε) on every branch, and the per-weight ranks over six directions are full. We expected the opposite and recorded it as a surprise.
- ω₁′ held fixed in the (∞,∞) family (at ω₂=∞). This *fails*: φᵤ has rank ≤1 on the multiplicity-2 weights, singlet-type branches have K/S=O(1) (for [3,2,0]: 0.465 at (4,4,2) and 0.388 at (4,3,3), identical to [3,1,0]'s duals (4,2,2) and (3,3,2)), and all branches with a 2-dimensional slice have O(1) CG errors (0.39 to 1.5; σ₂=1.5·10⁻⁷,1.7·10⁻⁸,8.7·10⁻¹⁵). Branches with one-dimensional slices still pass.

So the answer to Paul's question is: no cleaning is needed anywhere, in either family, for the adjoint, [3,1,0] or [3,2,0]. Specifically (5,3,2) has max|K/S|=1.8·10⁻⁶,1.4·10⁻⁶ (full family) and 2.4·10⁻⁶,1.3·10⁻⁶ ((∞,∞)) at ε=10⁻⁶, with the same doubling ratio as every other branch. What is untested is [4,2,0] (dimension 27, multiplicity-3 sector (6,4,2)), the other XXX trouble spot. The only weakly conditioned 2×2 extraction is (4,4,2) in the (∞,∞) family (σ₂∼3·10⁻⁴), still accurate to 1.4·10⁻¹⁷ at 200 digits; conditioning is a property of the chosen pair of grid points (best pairs have σ₂≈0.4 to 2.8), not of the family.[^28]

## 5. The jet formula: the overlap as a bilinear form in jets of the Q-functions

**Punchline.** For any eigenstate, the normalisation-free overlap c ⟨Φ_∞|Ψ⟩/⟨x₀^G|Ψ⟩ times M(θ₁)M(θ₂) equals a bilinear form in the monomials P³pₖₗ at the two sites, with a state-independent kernel that is unique, rational, and polynomial in ω₁,ω₁′,1/θ₁₂. The number of monomials per site is C(λ₁+2,2)-C(λ₂+1,2), namely 5,9,7 for [2,1,0],[3,1,0],[3,2,0].

### 5.1. The statement

We write Q_num(B)=⟨Φ_∞|Ψ_B⟩/⟨x₀^G|Ψ_B⟩; it is independent of the normalisation of Ψ_B. Equation (1) of the overview says, with mₖₗ=P³pₖₗ evaluated at θₛ on the dual Q-functions of the state,

```text
c Q_num M(θ₁)M(θ₂)=Σ_((kl),(k′l′))K_((kl),(k′l′)) mₖₗ(θ₁) m_(k′l′)(θ₂),   M=P³ p_(0,λ₂+1).        (35)
```

Several features are worth spelling out.

- *Which pairs.* (kl) runs over 0≤k<l≤λ₁+1 with not both k,l≤λ₂. The excluded pairs are those forced to vanish by P⁽ᵏ⁾=κ_(k,s)P for k≤λ₂ (Section 2.2). Their number is C(λ₁+2,2)-C(λ₂+1,2).
- *The lone factor is always P³.* The monomials P³⁽ᶜ⁾pₖₗ with c>0 are redundant: P³′p₀ₗ=P³p₁ₗ and P^(3″)p₀ₗ=P³p₂ₗ exactly (Frobenius).
- *M is the first non-forced coordinate.* M=P³p_(0,λ₂+1)∝Π′(θₛ)Q₁(θₛ)(ζ₃-ζ₂)e^(σθₛ)P³ by the Wronskian relation (adjoint), and it is the vacuum's own monomial.
- *Normalisation c.* The frame-tie scalar (15). ζ enters the formula only through c and the Q-functions; the kernel is ζ-independent.
- *The kernel.* Unique (the design matrix of products has full rank), rational, and homogeneous. For the adjoint K=G/4 with G the 5×5 matrix of Table 17; in terms of the symbols used in the notes there are two equivalent normalisations (G=4K, because the vacuum monomial is M_vac=2h²W′P³ per site, and Θ=K/(ω₁ω₁′)^(2λ₁-λ₂) in the companion frame), reconciled in Appendix A.

The formula was found, fitted, and derived. We describe each route.

### 5.2. The adjoint in detail

For [2,1,0] the pairs are (02),(03),(12),(13),(23), and since p₁₂=κp₀₂, p₁₃=κp₀₃ on every state, a basis of five one-site monomials is

```text
N=(P³p₀₂, P³′p₀₂, P³p₀₃, P³′p₀₃, P³p₂₃),  weights (2,3,3,4,5),        (36)
```

(the order of D2; P³′p₀ₗ=P³p₁ₗ). The kernel Gᵢⱼ(ω₁,ω₁′,θ₁₂) is a 5×5 matrix of polynomials listed in full in Table 17; a few entries show the flavour, with t=1/θ₁₂:

```text
G₁₁=-4ω₁²ω₁′² (3t²+tω₁-tω₁′-ω₁ω₁′),   G₃₃=ω₁²ω₁′²,   G₅₅=1,   G₄₅=-2(2t+ω₁).        (37)
```

Properties, all exact: G has rank 5; it is symmetric under site exchange in the form Gⱼᵢ(ω₁,ω₁′,θ₁₂)=Gᵢⱼ(ω₁′,ω₁,-θ₁₂); each entry is homogeneous of degree 10-wᵢ-wⱼ in (ω₁,ω₁′,1/θ₁₂) jointly, with w=(2,3,3,4,5) (checked by the scaling u→2u, ζ,ω→·/2: the refitted kernel equals G at (ω/2,2θ₁₂) to 5·10⁻⁷²); every entry is a polynomial in (ω₁,ω₁′) of total degree ≤6 (no poles, none in ω₁-ω₁′ or ω₁-ζⱼ) and in 1/θ₁₂ of degree ≤4; it depends on θ₁,θ₂ only through θ₁₂. Table 18 gives the same kernel in jets of the twisted functions e^(-ωₛu)Pʲ: in those variables the two vacua interact only through δ=ω₁′-ω₁ and t, the Gaudin image of the XXX statement that they couple only through s=a′/a.

Numerically: on all 64 states the formula holds with ratio exactly 1.0 and maximal relative deviation 2.8·10⁻⁵⁶ at (ω₁,ω₁′)=(5/2,7/4), and 2·10⁻⁵⁶ to 3·10⁻⁵⁵ at four further parameter points (other (ω₁,ω₁′), other θ, other ζ).[^29]

### 5.3. Route 1: the h→0 limit of the XXX ten-product formula

Agent D2 took the XXX identity for the adjoint (∞,∞) family,

```text
(⟨Φ_c|Ψ_c⟩)/(⟨x₀ᶜ|Ψ_c⟩) M⁽¹⁾_vacM⁽²⁾_vac=Σ_(c,c′∈{0,1})Q³(θ₁-ch)Q³(θ₂-c′h)Σₖ₌₁¹⁰Bk2⁽¹⁾[Fₖ,Gₖ] Bk2⁽²⁾[F′ₖ,G′ₖ],        (38)
```

with ladder points θₛ+mh, m=-3,…,0, 2×2 bracket determinants Bk2 of the dual Q-functions, M⁽ˢ⁾_vac=minorₛ(-3,-2)Q³(θₛ-h), and kernels rational in s=a′/a and t=θ₁₂/h (parsed from cell 62 of `L2-210-two-vacua-brackets.wb`; the 40 kernel laws and ten products reproduce the printed values with 0 mismatches), and expanded it formally at h→0 in the diagonal frame. The accounting is as follows.

The XXX normalisation κκ′=h⁻⁸/c and ⟨x₀ᵈ|ψ⟩=h²⟨x₀^G|ψ⟩ exactly (since x₀ is exactly linear in z) give Q_num^XXX∼h⁶ c Q_num^G. Each minor is O(h²) (the O(h) term is -W(θₛ)·(m′-m) and vanishes), so the naive order of the right-hand side is h⁴; the actual leading order is h¹⁰. The cancellations at orders h⁴,…,h⁹ are not accidental: they follow from *one* exact relation. The dual Baxter equation at x=θₛ has ν₃(θₛ)=0 and t₃(θₛ)=χ₃ν₁ν₂(x-h)ν₃(x-2h)|_(x=θₛ)=0, so

```text
Qʲ(θₛ-h)=ρₛ Qʲ(θₛ),  hence  minorₛ(-1,0)=0,        (39)
```

with ρₛ independent of j, exactly and to all orders in h, for every XXX eigenstate: the XXX generalisation of W(θₛ)=0. With only this relation imposed and all other jets random, the right-hand side has valuation exactly 10 (three random trials). Consequently the four ladder values of (Q²,Q³) at θₛ lie in a 2-plane, the minors are D times (-α,ρβ,ρ,β,1,0) with leading values (α,β,ρ)=(-2,3,1), and the leading coefficient depends only on the n=0 jets J_(0,0),J_(0,2),J_(0,3) (and J_(0,1)=κJ_(0,0)), through the Plücker coordinates; the O(h) corrections q₁ of the XXX Q-functions enter first at order h¹¹, so no extra equation for q₁ was needed. Two corrections to the plan came out: the vacuum monomial is M_vac=2h²W′(θₛ)P³(θₛ), not (5/2)h²W′P³ (the weight-2 minor contains an n=1 combination which the exact relation fixes to p/2), and only five of the twelve candidate monomials {P³⁽ᶜ⁾pₖₗ} occur. The result is (35) with the table of G.[^30]

A natural shortcut failed, and it is worth recording. Taking t→∞ at fixed a=1/ω₁ in the XXX kernels is wrong: of the 40 kernel laws, 25 survive and all equal s u_c[mm′]u_(c′)[nn′], a rank-one U⊗U, and the naive right-hand side divided by the vacuum monomials equals (1-a)³(1-a′)³, state independent. The Gaudin limit is a=e^(-hω₁)→1, not fixed a, and all the state dependence sits in the cancellation hierarchy, invisible to a rank-one naive kernel.[^31]

### 5.4. Route 2: the direct fit

Agent D1 fitted the kernel directly from Gaudin data, without reference to XXX. With target 4 c Q_num (P³p₀₂)₁(P³p₀₂)₂ and the basis N (25 unknowns) on 320 rows (5 twist points), the residual is below the fixed-point floor 10⁻⁹⁰, the rank is 25 (all products independent), the kernel is unique, all entries are rational, and it equals D2's G to 2.8·10⁻⁷⁴; per-ζ fits agree with the stacked one to 10⁻⁷². A hypothesis ladder led to this: Table 8 lists the levels. The direct-Q₁ bases do *not* fit with a ζ-independent kernel because p₀₃=Aₛ[(Λₛ-κ)Q₁+2Q₁′] with Λₛ=2σ+2/Π′(θₛ) contains σ=Σζ; the Plücker form is the ζ-free one. The controls fail as they should: lone index P² or P¹ instead of P³ gives residuals 0.99; jets built from site-1 or site-2 alone 0.37, 0.58; jets at wrong points 0.48; no M: 0.97.[^32] The kernel laws were then fixed on a 3×3 grid of (ω₁,ω₁′) and five θ variants: each entry is a polynomial of total degree ≤6 in (ω₁,ω₁′) with 28 coefficients, identical (exact rational comparison) to D2's analytic Gᵢⱼ in 25/25 entries.

**Table 8.** Hypothesis ladder for the adjoint fit (D1 F2): relative residual ‖Ak-y‖/‖y‖, 0.0 meaning below the 10⁻⁹⁰ floor; 320 rows (5 twist points); M=P³p₀₂ unless stated.

| level | basis per site | unknowns | rank | residual |
|---|---|---|---|---|
| a | N₅={P³p₀₂,P³′p₀₂,P³p₀₃,P³′p₀₃,P³p₂₃} | 25 | 25 | **0.0** |
| a | direct AₛQ₁⁽ᵏ⁾P³⁽ᶜ⁾, k≤2,c≤1 (6) | 36 | 36 | 0.294 |
| a | direct, k≤1,c≤3 (8) | 64 | 64 | 0.562 |
| a | D2's nine direct monomials with κ | 81 | 81 | 0.363 |
| b | Plücker {pₖₗP³⁽ᶜ⁾}, k<l≤3,c≤1 (12) | 144 | 64 | 0.0 |
| c | Plücker k<l≤3,c≤3 (24) | 576 | 256 | 0.0 |
| c | direct k≤3,c≤3 (16) | 256 | 256 | 0.0643 |
| ctl | N₅ without M₁M₂ | 25 | 25 | 0.966 |

### 5.5. Route 3: derivation from the jet basis

The third derivation (Section 6) constructs the Gaudin K tensor of Φ_c^∞ in the jet basis that replaces the SoV basis when h→0, and finds exactly 25 non-zero entries on N₅×N₅, equal to G/(4ω₁³ω₁′³) as an exact polynomial identity in (1/ω₁,1/ω₁′). The “4” is the jet normalisation M_vac=2h²P³p₀₂ per site, and c arises from the frame tie (14) exactly.

### 5.6. General Λ: [3,1,0] and [3,2,0]

The same fit pipeline, now with the sparse hypothesis tested first, gave the generalisation. For [3,2,0] (dim 15, 225 states) with c=Π_ω(ζ₂-ω)(ζ₃-ω)³(ζ₁-ζ₂)⁴ and M=P³p₀₃, the target is exactly (residual 0.0 at 80 digits) a bilinear form in P³ times the seven non-forced coordinates {p₀₃,p₀₄,p₁₃,p₁₄,p₂₃,p₂₄,p₃₄} per site; the kernel is unique (49 unknowns, rank 49 on 1125 rows), rational, and equal for every ζ point to ∼10⁻⁷². Dropping any one of the seven leaves a residual (0.41 without p₀₃ or p₀₄, 0.16,0.25,0.078,0.14,0.061 for p₁₃,p₁₄,p₂₃,p₂₄,p₃₄). For [3,1,0] (c=(ζ₂-ω)²(ζ₃-ω)³(ζ₁-ζ₂)² per vacuum, M=P³p₀₂), greedy backward elimination from 18 features ends with the pure Plücker basis of the nine non-forced pairs, 81 unknowns, rank 81, residual 0.0, unique rational kernel, K₉₉=1/9; the adjoint-like five-element basis fails (0.316). The pattern is the count

```text
#monomials=C(λ₁+2,2)-C(λ₂+1,2)=6-1, 10-1, 10-3=5, 9, 7.        (40)
```

In the ν-variables of [3,2,0] the seven monomials are (1,κ,κ₂)×(1,π₄) and π₃₄ with π₄=p₀₄/p₀₃, π₃₄=p₃₄/p₀₃; monomials of ν-weight ≤3 (13 of them) fail (0.0046), weight ≤4 fit.[^33] The kernel is neither symmetric at fixed θ₁₂ nor does K+Kᵀ express site exchange; the exchange symmetry is Kᵢⱼ(ω,ω,1/θ₁₂)=Kⱼᵢ(ω,ω,-1/θ₁₂) (deviation 5.8·10⁻⁶³).

### 5.7. Kernel laws

Table 9 collects the laws that each kernel obeys. They were established on grids of (ω₁,ω₁′) (9 values, 72 ordered pairs, 5 twist points each) and θ₁₂ (12 distinct values) with exact rational polynomial reconstruction, and the structure was verified at held-out points; Appendix B gives the data. The two sites exchange under θ₁₂→-θ₁₂, and the first non-forced coordinate carries weight 0, with the monomial weight eₖₗ=k+l-λ₂-1.

**Table 9.** Laws of the jet kernel for the three representations: e are the monomial weights relative to M; the ω-degree of entry (i,j) is ≤ total -eᵢ-eⱼ; t=1/θ₁₂. Sources: D1/D2 (adjoint), P4 ([3,2,0]), G1 ([3,1,0], `RESULTS-general-kernel.md` K0).

|  | [2,1,0] | [3,1,0] | [3,2,0] |
|---|---|---|---|
| monomials per site | 5 | 9 | 7 |
| pairs (kl) | (02), (03), (12), (13), (23) | (02), (03), (04), (12), (13), (14), (23), (24), (34) | (03), (04), (13), (14), (23), (24), (34) |
| e=k+l-λ₂-1 | 0,1,1,2,3 | 0,1,2,1,2,3,3,4,5 | 0,1,1,2,2,3,4 |
| max ω-degree | 3 each | 5 each | 4 each |
| t-degree | 4 | 8 | 7 |
| c | (ζ₂-ω)(ζ₃-ω)² per vacuum, (ζ₁-ζ₂)² | (ζ₂-ω)²(ζ₃-ω)³, (ζ₁-ζ₂)² | (ζ₂-ω)(ζ₃-ω)³, (ζ₁-ζ₂)⁴ |
| last entry K_last | 1/4 | 1/9 | 1/9 |
| monomials in Θ | 25 entries | 998 | 458 |

### 5.8. What went wrong first

The first fits all failed, and the failures are worth a paragraph. With the plan's vacuum monomial M=Q₁P³ and the assumption that the kernel is ζ-independent, no level of the ladder fit, not even the full jet space {J²ₐJ³_bJ³_c} with 1295 unknowns (training residual 4·10⁻², hold-out error 8·10²). The cause was normalisation: p₀₂=AₛQ₁ with Aₛ=(ζ₃-ζ₂)e^(σθₛ)Π′(θₛ) is ζ-dependent, and so is the frame-tie scalar. Once one divides by c and uses M=P³p₀₂, the kernel is exactly ζ-free. Likewise ⟨x₀^G|Ψ_B⟩/n_B turns out to be, to 2·10⁻⁷⁷, a constant times Q₁(θ₁)P³(θ₁)Q₁(θ₂)P³(θ₂) (exponent 1 on each, 0 on the other ten numbers), so the “x₀-free” target T(B)=⟨Φ_∞|Ψ⟩D_B(ω₀)/⟨Ω₀⊗Ω₀|Ψ⟩ is not constant times Q_num but c_T(ζ) Q_numM₁M₂: the two targets are the same fit with one free scalar per ζ point. And the one-site-fixed-a limit is wrong, as said above. The vacuum monomial of the plan was off by a factor 5/4 (M_vac=2h²W′P³, not 5/2 h²W′P³), and the twelve candidate monomials reduce to five.

## 6. The Gaudin B operator and the jet basis

**Punchline.** There is no single Gaudin B operator with simple spectrum. The XXX operator B(u) becomes, as h→0, a hierarchy ν³[1+hb₁+h²b₂+…] whose first member is Cartan-like and degenerate; the SoV basis becomes a ζ-independent *jet basis*, defined by a wavefunction rule in terms of eight one-site monomials, and the kernel G of the jet formula is exactly the K tensor of the vacuum family in that basis. Where the degeneracy lives is the diagonal frame, not B.

### 6.1. The XXX operator and its eigenbasis

In the new XXX conventions (gl(3), untwisted) B(u)=T₁₁(u) qmin_(12,12)(u+h)+T₂₁(u) qmin_(12,13)(u+h). It has simple spectrum ΠₛΠₖ(u-θₛ-mₖh) with labels (a|bc) per site and roots θₛ+{-1-a,-b,-c}h, λ₂≤a≤λ₁, 0≤c≤λ₂, c≤b≤a; for the adjoint there are eight labels per site, (2,2,1) being the vacuum with roots {-3,-2,-1}. Its left eigenbasis is the SoV basis, in which the wavefunction is ⟨x|Ψ_c⟩/⟨x₀|Ψ_c⟩=NₓΠₛM⁽ˢ⁾ₓ/M⁽ˢ⁾_vac with M_((a|bc))=minorₛ(-1-a,-b) Q³(θₛ-ch). We verified, exactly at h=1/2 on the adjoint, that B has degree 6 in u and in h, [B(u₀),B(u₁)]=0, all 64 label-rule eigenvalues are simple and the left eigencovectors are u-independent.[^34]

As h→0 all roots collapse to θₛ. So the “Gaudin B” cannot be one operator with simple spectrum; it has to be a hierarchy whose joint generalised eigenbasis is a jet basis, the h-expansion coefficients (divided differences) of the SoV covectors.

### 6.2. The hierarchy b₁,b₂,…

With t(x)=T(x)/ν(x)=1+hL(x)+h²C(x), ν=Π, Lᵢⱼ=ΣₛEˢᵢⱼ/(x-θₛ) the Gaudin Lax matrix and Cᵢⱼ=ΣₖE¹ᵢₖE²ₖⱼ/((x-θ₁)(x-θ₂)) the cross-site term, the exact all-order structure is

```text
B/ν³=(1+hs₁+h²s₂){ t₁₁(u)[t₁₁(u+h)t₂₂(u)-t₂₁(u+h)t₁₂(u)]
+t₂₁(u)[t₁₁(u+h)t₂₃(u)-t₂₁(u+h)t₁₃(u)]},        (41)
```

s₁=ν′/ν, s₂=1/((u-θ₁)(u-θ₂)), checked against the chain for b₀,…,b₆ at three values of u (orders 7, 8 vanish). Expanding B/ν³=1+hb₁+h²b₂+…,

```text
b₁=2L₁₁+L₂₂+ν′/ν,   b₂=2C₁₁+C₂₂+L₁₁′+L₁₁²+2L₁₁L₂₂-L₂₁L₁₂+L₂₁L₂₃+ν′/ν(2L₁₁+L₂₂)+s₂.        (42)
```

The first, b₁, is diagonal (Cartan) with eigenvalue Σₛ(1+a+b+c)/(u-θₛ) and has 25 eigenspaces of dimensions 1 (4 times), 2 (12 times) and 4 (9 times): it does not resolve the 64 labels. b₂ needs the cross-site Cᵢⱼ: it lies outside the span of {LᵢⱼLₖₗ,L′,L″,L,1} (171 residual entries) and inside it once the Cᵢⱼ are added. Further facts, all exact:

- the bₖ(u) do not commute individually, [bᵢ(u),bⱼ(v)]≠0 for all (i,j) except (1,1) and (6,6), but Σᵢ₊ⱼ₌ₙ[bᵢ(u),bⱼ(v)]=0 for n=2,…,8;
- compressed to the b₁-eigenspaces, b₂ has as eigenvalues the h²-coefficient of the label-rule series (25/25 eigenspaces) and splits every eigenspace completely: the 64 pairs (b₁,b₂^compressed) are distinct; b₂ alone has 62 simple roots and one double root with a Jordan block (the bₖ are not individually diagonalisable);
- B commutes with the grading H_tot=Σₛ(2Eˢ₁₁+Eˢ₂₂), with level block sizes 1,4,8,12,14,12,8,4,1, and every SoV covector lives in the level D(x)=Σₛ(6-(1+a+b+c));
- the Lie algebra generated by all coefficient matrices B⁽ᵐ⁾ₖ is *solvable* on every level block (derived series: 80→67→45→10→0 on the 14-block), so the bₖ(u) are simultaneously triangularisable block by block, with the level flag.
- **No single polynomial B_G(u) of degree ≤3 in L,L′ triangularises the jet basis** (tried: b₁,b₂,b₃,b₂-b₁²/2 in the jet basis; b₁ has only 25 distinct values and is not triangular inside level blocks). The operator that does resolve all 64 labels is the b₁-block compression of b₂, which needs the projectors onto H₁,H₂ eigenspaces: not a polynomial of small degree in L.[^35]

### 6.3. Where the degeneracy lives

We expected degeneracy in B: at h→0 all 64 SoV covectors become proportional. That is refuted. In the XXX companion frame the 64 covectors ⟨x|(h), normalised by ⟨x|r=1, are rational in h of degree ≤5 over ≤5 (common denominator; degree tally 0:11, 1:23, 2:19, 3:7, 4:3, 5:1), validated at six unseen rational h (64/64), and their *leading terms at h=0 are pairwise distinct and linearly independent (rank 64)*. The degeneracy lives in the diagonal frame, where the Gaudin states live: there the covectors are fᵈₓ(h)=fₓ(h)(P⊗P) with P=P_phys⁻¹(e^(hζ)), every entry O(h) with singular values h(1,h,h,h²,h²,h³,h³,h⁴). All 64 have valuation h² and *all 64 leading covectors are proportional to x^G⊗x^G* (rank 1). The ranks of the span of Taylor coefficients of fᵈₓ/h² to order K=0,…,8 are 1,5,13,25,39,51,59,63,64. Saturating over ℚ[[h]] by echelon reduction (dividing by h when the leading term vanishes) yields the *jet basis* with defects per level 1,4,8,12,14,12,8,4,1, the coefficients of (1+2t+2t²+2t³+t⁴)²: per site eight jets of weights 0,1,1,2,2,3,3,4, the one-site depth of P_phys.[^36]

### 6.4. The canonical jet basis and its wavefunction

The jet basis is made canonical by the wavefunction rule. Its one-site monomials are, in terms of the dual Q-functions at θₛ,

```text
m₀ =P³p₀₂,  m₁ =P³′p₀₂,  m₂ =P³p₀₃,  m₃ =P³′p₀₃,
m₄ =P³′p₁₂,  m₅ =P³p₂₃,  m₆ =P³′p₁₃,  m₇ =P³′p₂₃,        (43)
```

with weights w=c+k+l=2,3,3,4,4,5,5,6: eight monomials, one per site label (=dimV). These are what the twelve P³⁽ᶜ⁾pₖₗ, c≤1, l≤3 reduce to (p₀₁=0, p₁₂=κp₀₂, p₁₃=κp₀₃). The 64 products mᵢ(θ₁)mⱼ(θ₂) are independent functions on the 64 states (σₘᵢₙ/σₘₐₓ=6·10⁻¹²). The canonical jet covectors ⟨i,j|ᶜ are defined by

```text
(⟨i,j|ᶜΨ_c⟩)/(⟨x₀|Ψ_c⟩)=(mᵢ(θ₁) mⱼ(θ₂))/M₁M₂,        (44)
```

in the Gaudin companion frame (Ψ_c=(P_phys⁻¹⊗P_phys⁻¹)|_(z=ζ)Ψ_d). The 64×64 matrix Y expressing them in the physical basis is *exactly ζ-independent and rational*: identical exact rationals from the states at P, P′, X₁ (833 non-zero entries, inverse 975). Only five of the 64 covectors are product covectors (the 8×8 coefficient matrix has rank 1): (0,0),(0,1),(0,2),(1,0),(2,0); the other 59 couple the two sites through 1/θ₁₂ (ranks 2 to 8). Examples: ⟨00|=e₀₀, ⟨01|=-e₀₁, ⟨02|=-2e₀₂+e₀₁-2t e₀₀ with t=1/θ₁₂. The rule (44) was verified to 6·10⁻⁷⁶ (P), 1.4·10⁻⁷⁴ (P′), 1.5·10⁻⁷³ (X₁) on all 64×64 with Y built at P only; the saturation basis is block-triangular with respect to it, the filtration being the jet-weight filtration (level d=wᵢ+wⱼ-4).[^37]

### 6.5. The K tensor and the derivation of G

Expanding the companion-frame family Φ_c^∞=Ω_c^∞[a]⊗Ω_c^∞[a′] in the canonical jet basis gives its K tensor Θ_μ=(Ω_c^∞[a]⊗Ω_c^∞[a′])Y⁻¹, with exactly *25 non-zero entries out of 64*, all on N₅×N₅ (in the saturation basis: 36 non-zero entries). Writing Θᵢⱼ=Θ(a,a′,1/θ₁₂) one has

```text
Θᵢⱼ(a,a′,t)=(Gᵢⱼ(1/a,1/a′,t))/(4 ω₁³ω₁′³)   (adjoint), i.e.  Θ=K/((ω₁ω₁′)^(2λ₁-λ₂)),        (45)
```

as an exact polynomial identity in (a,a′), 25/25 entries, at three values of θ₁₂ (31/35,29/45,5/6). This is the third derivation of the kernel. The kernel has rank 1 at 1/θ₁₂=0 (the sites decouple, G becomes a product) and rank 5 otherwise, so the whole site coupling is through 1/θ₁₂. For equal twist (finite ω₂) all 64 jet coefficients are non-zero (levels 1,4,8,12,14,12,8,4,1) and the jet expansion reproduces ⟨Φ|ψ⟩/⟨x₀|ψ⟩ to 6·10⁻⁷⁵; moreover Q_num^eqM₁M₂/detN_B(ω) is state-independent (spread 1.2·10⁻⁷⁵, value 12.5381319888976 at ζ=P, θ=(2/7,-3/5)), i.e. the 6×6 determinant of Section 2.5 is a Q-number of the same kind.[^38]

## 7. The overlap as a function of the Laurent data

**Punchline.** Every Plücker monomial is M times a rational function of the Laurent data (Hₛ,aₛ,bₛ) of the state at the two sites, so c Q_num=ΣᵢⱼGˢᵢⱼνᵢ(θ₁)νⱼ(θ₂) and neither the Q polynomials nor their normalisation enter: the Q-number is a function of eigenvalue data.

### 7.1. The Frobenius reduction (adjoint)

Near θₛ, with x=u-θₛ and exponents {0,2,4} for the adjoint, every solution is Pʲ=Aⱼf₀+Bⱼf₂+Cⱼf₄ with Frobenius series

```text
f₀=1+κx+γ₂x²+γ₃x³+…,   f₂=x²+d₃x³+….        (46)
```

To derivative order 3 only f₀,f₂ matter, and the Plücker coordinates are pₖₗ=[A,B]₂₃ wₖₗ with [A,B]₂₃=A₂B₃-A₃B₂ and w₀₂=2, w₀₃=6d₃, w₁₂=2κ, w₂₃=12(γ₂d₃-γ₃), w₀₁=0. Hence Nᵢ=Mνᵢ with

```text
ν=(1, κ, π, κπ, εₛ),   π=p₀₃/p₀₂=3d₃,   εₛ=p₂₃/p₀₂=6(γ₂d₃-γ₃),        (47)
```

and M cancels from the left-hand side of (35):

```text
c Q_num=ΣᵢⱼGˢᵢⱼ(ω₁,ω₁′,θ₁₂) νᵢ(θ₁) νⱼ(θ₂),   Gˢ=G/4.        (48)
```

The Frobenius recursion, with δ=θₛ-θₜ (t≠s), σ₁=Σζ, e₂=e₂(ζ), e₃=e₃(ζ), gives exactly (sympy)

```text
κ =aₛ/3,   π=2Hₛ-aₛ-2σ₁-6/δ,
εₛ =-⅔ [Hₜaₛ/δ-aₛbₛ+aₛe₂+3aₛ/δ²-3aₜ/δ²-3bₜ/δ-3e₃],   bₛ=(aₛ(Hₛ-aₛ))/3,        (49)
```

the last being the first no-log condition (21) with c₂=3. We verified (49) against Nᵢ/M built from the Q polynomials to 10⁻¹⁷⁴–10⁻¹⁹³ on all six regime branches, both sites, at ε=10⁻⁶/8, and the no-log relation to 10⁻¹⁹⁹. Thus Q_num is a function of (H₁,a₁,H₂,a₂), ζ, θ₁₂, ω₁,ω₁′ only.[^39] The reduction also improves the limit identity by 14 or more digits (the residuals of Section 9 drop from 10⁻⁴³–10⁻⁷⁵ to 10⁻⁵⁵–10⁻¹⁵⁷) and requires ε-orders ≤5 instead of 8.

### 7.2. Why this must hold for every Λ

The argument is local. For exponents {0,λ₂+1,λ₁+2} the jet of every solution up to order λ₁+1 lies in the 2-plane spanned by the exponent-0 and the exponent-(λ₂+1) Frobenius vectors; the coefficients of these vectors are rational functions of (Hₛ,aₛ,bₛ), the other site's data through the regular parts, and eₖ(ζ). Hence every non-forced pₖₗ=Dₛ νₖₗ with Dₛ the 2×2 determinant of the leading components of the two non-trivial solutions (the analogue of [A,B]₂₃) and νₖₗ rational, and M=P³Dₛν_(0,λ₂+1), so that P³Dₛ cancels between the sides. We have *not* derived the νₖₗ for λ₁≥3, and the untwisting limit for [3,1,0], [3,2,0] is not done; the general statement is an expectation based on the adjoint and on the verified exponent and jet-rank structure (Section 2.2). This is work package B of round 5.

## 8. The general kernel

**Punchline.** In the companion normalisation Θ=K/(ω₁ω₁′)^(2λ₁-λ₂) the kernel is a polynomial in (a,a′,t)=(1/ω₁,1/ω₁′,1/θ₁₂). At t=0 it factorises into one-site kernels θᵢ(a)θⱼ(a′), and the one-site kernel is a dressed Wronskian with a closed form valid for all Λ. Every two-row kernel is a block of the one-row master kernel R^(λ₁) times θ⊗θ′ (block theorem). R^(λ₁) is known to first order in t for all pairs and to all orders for [l,l,0]; the general four-index form is open.

### 8.1. Conventions for this section

Put n=λ₁+1 and let pairs (kl) run over 0≤k<l≤n with l≥λ₂+1. Let a=1/ω₁, a′=1/ω₁′, t=1/θ₁₂, x=at, y=a′t, and let the monomial weight be eₖₗ=k+l-λ₂-1 (relative to M). The jet formula in canonical normalisation is

```text
c Q_num M₁M₂=ΣKᵢⱼ mᵢ(θ₁)mⱼ(θ₂),   Θ:=K/((ω₁ω₁′)^(2λ₁-λ₂)),        (50)
```

with M=P³p_(0,λ₂+1) exactly, no factor 4. (D2's G=4K for the adjoint, from its vacuum monomial 2h²P³p₀₂; D3's identity (45) reads Θ=G/(4ω₁³ω₁′³).) The common exponential factor e^((ζ₂+2ζ₃)θₛ) of all monomials of a site is dropped, since the formula is homogeneous of degree one in the monomials, so every quantity is an exact rational. We use the closed forms below only after checking them at points not used in their discovery.

### 8.2. The one-site kernel

At t=0 the sites decouple and Θ→θ(a)⊗θ(a′) exactly (symbolic identity entry by entry, for all computed two-site representations). The one-site kernel θₖₗ(a) is exactly computable for any Λ: for one site the eigenstates are GT weight vectors, the Q's are exact rationals, and the same fit gives, with rank equal to the number of pairs and residual exactly zero (eight twist points, all diagonal states, a-dependence from four values of ω),

```text
θₖₗ(a)=cₖₗ a^(k+l-λ₂-1),   cₖₗ=(-1)^(k+l-λ₂-1)C(λ₁-λ₂,l-λ₂-1)C(λ₁+1,k)(l-k)/((λ₂+1)C(l,λ₂+1)),        (51)
```

independent of ω and ζ, with c_(0,λ₂+1)=1. The formula was guessed from Λ≤[4,2,0] and then correctly *predicted* the other eight representations [3,0,0],[3,3,0],[4,0,0],[4,3,0],[4,4,0],[5,1,0],[5,2,0],[5,3,0] (17 in all). Examples (pairs in lexicographic order): [2,1,0]: (02) 1,(03) -½,(12) -3/2,(13) 1,(23) -½; [1,1,0]: (02) 1,(12) -1; [3,1,0]: 1,-1,⅓,-2,8/3,-1,-2,1,-⅓. Appendix C has the table up to λ₁=4.

The structure is a decomposable bivector. For λ₂=0 the generating symbol is

```text
Σ_(k<l)θ^([λ₁,0,0])ₖₗ(x₁ᵏx₂ˡ-x₁ˡx₂ᵏ)=(x₂-x₁)[(1-ax₁)(1-ax₂)]^(λ₁),        (52)
```

exact for λ₁=1,…,5, and in jet language, for general λ₂,

```text
Σ_(k<l)θₖₗ(a) mₖₗ=ε_Λ a^(-λ₂) P³ W[(1-a∂)^(λ₁)P², (1-a∂)^(λ₁)P³](θ),   ε_Λ=(-1)^(λ₂)(λ₂!(λ₁-λ₂)!)/(λ₁!),        (53)
```

a Wronskian of the two dual Q's after dressing by (1-a∂)^(λ₁). In Taylor-normalised jets p̃ₖₗ=pₖₗ/(k! l!) the kernel is ĝ∧ĥ with ĝₖ=(-a)ᵏ/(n-1-k)!, ĥₖ=(-a)ᵏ/(n-k)!: it is a Plücker vector of a 2-plane, the image of the lowest-weight vector (ξ₁-ξ₂) of the spin-λ₁ irrep of an auxiliary sl₂ under exp(-aS₊). A one-site top-state sum rule fixes the normalisation (Section 8.7).[^40]

### 8.3. The block theorem

The second structural fact: the λ₂-dependence of the two-site kernel is trivial. For all non-forced pairs i,j,

```text
Θ^([λ₁,λ₂,0])ᵢⱼ(a,a′,t)=[(λ₂!(λ₁-λ₂)!)/(λ₁!) ]²(aa′)^(-λ₂) Θ^([λ₁,0,0])ᵢⱼ(a,a′,t),        (54)
```

so that, equivalently, Θᵢⱼ=θᵢ(a)θⱼ(a′)R^(λ₁)ᵢⱼ(at,a′t) with a universal, λ₂-independent Laurent polynomial R^(λ₁)_((kl),(k′l′))(x,y), R(0,0)=1, homogeneous of degree 0 in (a,a′,θ₁₂) (a function of x,y only; the argument is (x,y)=(at,a′t)). All the λ₂ content is in the one-site kernel θ. In the K normalisation this has a particularly clean form (a rewriting of (54) that we also checked on the stored exact kernels, see below):

```text
K^([λ₁,λ₂,0])ᵢⱼ=[(λ₂!(λ₁-λ₂)!)/(λ₁!) ]²K^([λ₁,0,0])ᵢⱼ  (non-forced pairs).        (55)
```

For λ₁=3 the prefactor is 1/9 for both λ₂=1,2, so the [3,1,0] and [3,2,0] kernels in the K normalisation coincide on the seven shared pairs. We confirmed this at two parameter points, (ω₁,ω₁′,θ₁₂)=(5/2,7/4,31/35) and (-2/5,3,29/45), entry by entry on the exact stored polynomial laws (0 mismatches in 49 entries each), and confirmed (54) for [3,1,0] against [3,0,0] at (a,a′,t)=(2/5,4/7,35/31) (0 mismatches in 81 entries). G1's own verification: symbolic for λ₁≤3 (all λ₂), symbolic for [4,4,0] against [4,0,0], direct independent fits for [4,1,0],[4,2,0],[4,3,0] at two θ₁₂ and two (ω,ω′) (deviation ≤3·10⁻⁵⁸). The corner block (pairs (k,λ₁+1)) equals N(aa′)^(λ₁-λ₂) times the single-index kernel Θ^([λ₁,λ₁,0]), N=[λ₂!(λ₁-λ₂)!/λ₁!]².

The block theorem is useful as a predictor: [4,2,0] has 12 pairs and Θ=1/36(aa′)⁻²Θ^([4,0,0])| (verified directly, deviation ≤8·10⁻⁶²), [5,1,0] has 20 pairs, Θ=1/25(aa′)⁻¹Θ^([5,0,0])|, a-degree 9 and t-degree 16. The [5,1,0] prediction was not tested directly: a corner-block check was started but the 1225-state build was too slow and was stopped.

### 8.4. Structure: polynomiality, weights, degrees

Every kernel computed is a polynomial in (a,a′,t) with rational coefficients; each monomial aᵖa′^qtⁿ of entry (i,j) obeys p+q-n=eᵢ+eⱼ with 0≤p,q≤2λ₁-λ₂ (maximum attained) and the exchange rule

```text
Θᵢⱼ(a,a′,t)=Θⱼᵢ(a′,a,-t)  exactly.        (56)
```

Equivalently, writing Θᵢⱼ=ΣₙtⁿΘ_(ij,n)(a,a′), each Θ_(ij,n) is homogeneous of degree eᵢ+eⱼ+n in (a,a′) (the “torus” law: ω₁,ω₁′ scale like ζ and θ₁₂ like 1/ζ). The t-degree is Dₜ=4(λ₁-1)-max(0,2λ₂-3) for λ₁≥2, Dₜ=1 for λ₁=1, and 2l-1 for [l,l,0]; the one-row sequence 1,4,8,12 is 4(l-1) from l=2 on. The earlier guess C(l+2,2)-2=13 for [4,0,0] is refuted (12). Predicted: [5,0,0], [5,1,0]: 16; [5,2,0]: 15; [5,3,0]: 13. Table 10 summarises.

**Table 10.** Structure of the two-site kernels (G1 `RESULTS-general-kernel.md` K0, K2, K3): number of pairs, a-degree (=2λ₁-λ₂), t-degree, generic rank of Θ, and ranks of the t-expansion coefficients Θₙ. All kernels satisfy Θ(t=0)=θ⊗θ′, the exchange rule, and the top-state sum rule exactly. Rank-deficient at a generic point: [1,0,0],[3,0,0],[3,1,0],[4,0,0],[4,1,0] (cause not understood).

| Λ | pairs | a-deg | t-deg | generic rank | ranks of Θₙ, n=0,1,… |
|---|---|---|---|---|---|
| [1,0,0] | 3 | 2 | 1 | 2 | 1,2 |
| [1,1,0] | 2 | 1 | 1 | 2 | 1,2 |
| [2,0,0] | 6 | 4 | 4 | 6 | 1,2,4,6,6 |
| [2,1,0] | 5 | 3 | 4 | 5 | 1,2,4,5,4 |
| [2,2,0] | 3 | 2 | 3 | 3 | 1,2,3,2 |
| [3,0,0] | 10 | 6 | 8 | 8 | 1,2,4,8,10,10,10,8,6 |
| [3,1,0] | 9 | 5 | 8 | 8 | 1,2,4,8,9,9,8,6,4 |
| [3,2,0] | 7 | 4 | 7 | 7 | 1,2,4,7,7,6,4,2 |
| [3,3,0] | 4 | 3 | 5 | 4 | 1,2,3,4,3,2 |
| [4,0,0] | 15 | 8 | 12 | 12 | 1,2,4,8,13,15,15,15,15,12,10,8,6 |
| [4,4,0] | 5 | 4 | 7 | 5 | 1,2,3,4,5,4,3,2 |
| [4,1,0] | 14 | 7 | 12 | 12 | 1,2,4,8,13,14,14,14,13,10,8,6,4 (block thm) |
| [4,2,0] | 12 | 6 | 11 | 12 | 1,2,4,8,12,12,12,11,9,6,4,2 (block thm) |
| [4,3,0] | 9 | 5 | 9 | 9 | 1,2,4,7,9,9,8,6,4,2 (block thm) |

The supports shrink with n: the nonzero entries of Θₙ per order are 25,23,19,12,4 of 25 for the adjoint and 49,47,43,37,29,19,9,2 of 49 for [3,2,0]; the top order has only two to four non-zero entries. The “two channels per order” doubling of the rank ends at n=3 (for [4,0,0]: 13 at n=4, not 16).

### 8.5. The first-order term

The coefficient of t¹ is known in closed form for all pairs and all representations. With ρ(p)=(n-k)(n-l), σ(p)=k+l-1, f₂(p)=2(n-1)-σ(p), N₀=n(n-1), and p=(kl), q=(k′l′),

```text
Θ_(1,pq)=t θₚ(a)θ_q(a′)[α_(pq)a+β_(pq)a′+γ_(pq)a²/a′+ε_(pq)a′²/a ],        (57)
```

```text
α_(pq)=-3/n² ρ(p)(ρ(q)-N₀)-2/n ρ(p)σ(q)+1/n f₂(p)(ρ(q)-N₀)+f₂(p)σ(q),  β_(pq)=-α_(qp),        (58)
```

```text
γ_(pq)=(k′l′ρ(p))/n²,  ε_(pq)=-(kl ρ(q))/n².        (59)
```

This formula was fitted on λ₁=1,2,3 and then reproduced the λ₁=4 single-index kernel and all λ₂-blocks exactly: 0 mismatches over 1062 random-point entry tests (λ₁≤4) and 675 entries of the held-out [4,4,0] and [4,0,0]. As a matrix on the pair space α,β have rank 2 and γ,ε rank 1; the factors (n-k)(n-l) and k′l′ are the raising and lowering amplitudes of the spin-n/2 representation of the auxiliary sl₂ acting on the jet index (double hops: γ raises both indices of site 1 and lowers both of site 2). On a single index (l=l′=n): α=(λ₁-k)k′, β=-(λ₁-k′)k, γ=ε=0. Θ₁ has rank exactly 2, with θ(a) in its column space and θ(a′) in its row space. Beyond first order we know the extreme terms. At order 2 the exponent pairs xᵘyᵛ (u+v=2) are (u,v)=(-2,4),(-1,3),(0,2),(1,1),(2,0),(3,-1),(4,-2) with coefficient ranks 1,2,4,5,4,2,1; the extreme ones are ½(n-k)₂(n-l)₂k′₂l′₂/(n(n-1))² and its mirror. General m (exact on all 225 entries of [4,0,0] for m=1,…,4): the coefficient of x²ᵐy⁻ᵐ in R_(pq) is

```text
1/(m!) ((n-k)ₘ(n-l)ₘ(k′)ₘ(l′)ₘ)/((n)ₘ²)=(m!)³(C(n-k,m)C(n-l,m)C(k′,m)C(l′,m))/(C(n,m)²),        (60)
```

with falling factorials (z)ₘ, and a sign (-1)ᵐ and exchanged roles for the mirror y²ᵐx⁻ᵐ.

For λ₁=1 the full kernel is an R-matrix: Θ^([1,λ₂,0])=(1+z) θ(a)⊗θ(a′)-z θ(a′)⊗θ(a) with z=t/(ω-ω′)=aa′t/(a′-a), i.e. a permutation exchange of the two coherent states with the rational R-matrix R(u)=u+tP.

### 8.6. The single-index family [l,l,0] to all orders

For [l,l,0] (one Plücker index k∈{0,…,l}) the kernel is

```text
Θ_(kk′)=(-1)ᵏ⁺ᵏ′C(l,k)C(l,k′)aᵏa′ᵏ′Σₘ₌₀ˡ⁻ᵏΣ_(m′=0)ˡ⁻ᵏ′C(l-k,m)C(l-k′,m′)F_(mm′)(k,k′) (at)ᵐ(-a′t)ᵐ′,        (61)
```

```text
F₀₀=1,   Fₘ₀=(k′)ₘ,   F_(0m′)=(k)_(m′),   F_(mm′)=(k+m+1)_(m′-1)(k′+m′+1)ₘ₋₁ (kk′+m′k+mk′),        (62)
```

here with rising Pochhammer symbols; uniformly F=Γ(k+m+m′)Γ(k′+m+m′) [(k+m)(k′+m′)-mm′]/(Γ(k+m+1)Γ(k′+m′+1)). It is exact for l=1,…,4 (all 4+9+16+25 entries) and *predicted correctly for l=5*: direct fits of [5,5,0] (441 states, 36 unknowns, rank 36, residual 0.0) agree to 1.5·10⁻⁵⁹–9.2·10⁻⁵⁸ at three parameter points. Its symbol is nice. With A=1-aξ, B=1-a′η and Θ̂(ξ,η)=ΣΘ_(kk′)ξᵏηᵏ′,

```text
Θ̂^([l,l,0])=AˡBˡ+t l²aa′(ξ-η)Aˡ⁻¹Bˡ⁻¹+…,        (63)
```

every tⁿ term with n≥1 is divisible by (ξ-η) (so on ξ=η, plane waves like the top state's P³=e^(ζ₃u), there is *no coupling*), and the top coefficient is t²ˡ⁻¹ times (-1)ˡ⁺¹((2l-2)!(2l)!)/(2(l-1)! l!)(aa′)ˡ(ξ-η), with values 1,12,720,100800,25401600,10059033600 for l=1,…,6 (l≤4 from fitted kernels, l=5,6 from (61)). In the variables Ξ=ω-ξ, H=ω′-η, ν=ω-ω′, δ=ξ-η=ν-Ξ+H, the rescaled K̂=(ωω′)ˡΘ̂ is ω-free: K₁=ΞH+tδ, K₂=Ξ²H²+4tδΞH+2t²δ(3δ-2ν)-12t³δ, K₃=Ξ³H³+9tδΞ²H²+18t²δΞH(2δ-ν)+…. In a Laplace–Borel form, Θ̂=Σ_(K,K′)C(l,K)C(l,K′){E[(u+vs)^K(u′+v′s′)^(K′)]-KK′vv′ E[(u+vs₂)^(K-1)(u′+v′s₂′)^(K′-1)]} with u=-aξ, v=at, u′=-a′η, v′=-a′t and independent Gamma variables s∼Γ(K′), s′∼Γ(K), s₂∼Γ(K′+1), s₂′∼Γ(K+1) (the Pochhammer symbols are the moments).

Its determinants have a pattern. With ν=ω-ω′=(a′-a)/(aa′), T=l(l+1)/2 and dₗ=⌈l/2⌉, detΘ^([l,l,0])=cₗ(aa′)^(sₗ)tᵀν^(T-dₗ)pₗ(ν,t) with pₗ homogeneous of degree dₗ:

```text
p₁=ν+t,   p₂=ν+3t,   p₃=ν²+6tν+3t²,   p₄=ν²+10tν+15t²,   p₅=ν³+15tν²+45t²ν+15t³.        (64)
```

The kernel loses rank at ν=0 (ω=ω′, where CG-completeness fails) and on the roots of pₗ. We note, as a curiosity, that Θ is rank-deficient at a generic point for [1,0,0],[3,0,0],[3,1,0],[4,0,0] (ranks 2/3, 8/10, 8/9, 12/15): lost directions of unknown origin; all other computed representations have full rank.

### 8.7. The top-state sum rule

A strong exact constraint comes from the top state e₁⊗e₁, whose Q-functions are explicit: P³=e^(ζ₃u), P²=e^(ζ₂u)p₂ with p₂=Σ_(k≥0)(ζ₃-ζ₂)⁻ᵏ(Π^(λ₂))⁽ᵏ⁾ (from the Wronskian relation with Q₁=e^(-ζ₁u); monic, p₂′-(ζ₃-ζ₂)p₂=-(ζ₃-ζ₂)Π^(λ₂)). The sum rule is

```text
ΣᵢⱼΘᵢⱼ mᵢ^top(θ₁) mⱼ^top(θ₂)=M^top(θ₁)M^top(θ₂) R(a)R(a′),   R(a)=(1-ζ₂a)^(λ₁-λ₂)(1-ζ₃a)^(λ₁),        (65)
```

exact for ten representations (two-site; fully symbolic in ζ₂,ζ₃,θ₁,θ₂,a,a′ for [1,0,0],[1,1,0],[2,2,0],[2,1,0], exact rational at 3–4 random points for the rest) and ten one-site ones. The right-hand side is not M^topM^top, as the plan had it: the companion-frame top state is not e₁⊗e₁ because P_phys⁻¹e_hw is not parallel to e_hw, and c a^(2λ₁-λ₂)/x_hw=R(a) with x_hw=(ζ₂-ζ₁)^(λ₂) exactly. At a=0 the rule gives R(0)=1 and each further power of a one linear relation between the top-state jets and Θ: every candidate kernel must satisfy it.

### 8.8. One-row representations are the master kernels

The λ₂=0 representations work, with the same monomials, and are the master objects of the block theorem. The oper has exponents {0,1,λ₁+2}, c₂=0 (a scalar, see Section 1.2), M=P³p₀₁ with nothing forced, and C(λ₁+2,2) monomials (3,6,10 for [1,0,0],[2,0,0],[3,0,0]). Nothing breaks, and in particular the plan's worry that S₂'s double pole would be non-scalar did not materialise.

### 8.9. The auxiliary sl₂ and what is open

The jet index space (jets of order ≤λ₁+1, dimension λ₁+2) carries an sl₂ of spin (λ₁+1)/2 acting on symbols by the Möbius map (Uₛf)(ξ₁,ξ₂)=Πᵢ(1-sξᵢ)^(λ₁+1)f(ξᵢ/(1-sξᵢ)); the one-site kernel obeys θ⁽ᵃ⁺ˢ⁾=Uₛθ⁽ᵃ⁾ exactly (checked symbolically, λ₁≤4), a coherent state of the spin-λ₁ irrep of Λ²V. The two-site kernel is *not* covariant under the diagonal shift a→a+s, a′→a′+s (tested, λ₁=1,2). Further tests that failed: Θ is not a compound matrix Λ²M (Θ₀, Θ₁ are Plücker, Θ₂,Θ₃,… and Θ are not); the “box exchange” hypothesis (the kernel as a sum over permuted parameters u=(ω+ct,ω′+c′t) with shifts up to ±3 in t/2 units) is refuted for λ₁=2 with the fused-R picture not literally a sum over parameter permutations; and the XXX R-matrix poles (ω-ω′+mt) appear only in determinants of Θ, not in its entries. The number of t-orders is only linear in λ₁, not quadratic as a fusion count of λ₁² R-matrix factors would suggest.

**Open:** the four-index master kernel R^(λ₁)_((kl),(k′l′))(x,y) beyond first order. We know the single-index family (all orders), the first-order term (all pairs), the extreme multi-hops, universality in λ₂, the rank pattern, and the sum rule; the full-λ₁=5 grids and the full (a′,t) polynomial grids for [4,1,0],[4,2,0],[4,3,0] (replaced by the block theorem) were not done.

## 9. The untwisting limit and the CG generating function

**Punchline.** In the regime ζ=εφ, ω=εω̂ the jet formula is a statement about the ε-jet of the Q-functions, with cancellations of depth up to six orders. Reduced to Laurent data it becomes a polynomial identity in (ω̂₁,ω̂₁′) whose solution is the CG generating function. For all six HW branches of the adjoint the identity holds to 10⁻⁵⁵–10⁻¹⁵⁷ and returns every CG coefficient, φ-independent to 2·10⁻⁷⁷.

### 9.1. The regime and the difficulty

Let ζ=εφ, ω₁=εω̂₁, ω₁′=εω̂₁′ and ε→0. Chain side (Section 4): ⟨Φ_∞|Ψ_B⟩=⟨φᵤ⊗e₁ᵀ|Ψ₀⟩(1+O(ε)) for every HW-limit branch Ψ_B=Ψ₀+εΨ₁+⋯ with Ψ₀ an exact HW vector of weight μ, G-norm 1; so

```text
⟨φᵤ(ω̂₁,ω̂₁′)⊗e₁|Ψ₀⟩=Σ_ν[e₁ᵀρ(n′)]_ν C_ν,   C_ν=⟨e_ν⊗e₁|Ψ₀⟩,        (66)
```

is the CG generating function: a polynomial in the rotation parameters whose coefficients are the CG coefficients of the GT basis. We want to write the left-hand side as a Q-expression using (35).

The difficulty is that the limit is delicate on the Q side. At ζ=0 the three dual Q-functions e^(ζⱼu)pⱼ lose their exponentials and their Wronskians vanish (they carry factors ζ_b-ζₐ), so the untwisted limits of the three branch Q's are linearly dependent: the limit of the solution space is a ζ-jet (p⁽⁰⁾ and its ε-derivatives), not just a triple of polynomials. Moreover the jet formula's total is ε^(N_B) smaller than its individual terms, with N_B=8,7,7,6,6,4 for (4,2,0),(4,1,1),(3,3,0),(3,2,1)×2,(2,2,2), while single terms are O(ε²…ε³): a cancellation of up to six orders.

### 9.2. ε-order bookkeeping

The ingredients scale as follows (for the adjoint, at generic φ; φ′ gives identical exponents). The frame scalar is c=ε⁸ĉ with ĉ=(φ₂-ω̂₁)(φ₃-ω̂₁)²(φ₂-ω̂₁′)(φ₃-ω̂₁′)²(φ₁-φ₂)²; the kernel splits by degree in (ω₁,ω₁′), Gˢ(εω̂)=Σ_dεᵈG^(s[d])(ω̂) (at fixed θ₁₂); M₁M₂∼ε² on all branches; ⟨Φ_∞|Ψ_B⟩∼ε⁰ at all 16 grid points; and

```text
⟨x₀^G|Ψ_B⟩∼ε^(2+k_B),   k_B=0,1,1,2,2,4  for  (4,2,0), (4,1,1), (3,3,0), (3,2,1)×2, (2,2,2).        (67)
```

Table 11 gives the exponents. The interpretation of k_B: on the exact HW kets of every non-top μ, ⟨x(φ)⊗x(φ)|v⟩=0 identically (≤4·10⁻²⁰² at 200 digits, both φ, the whole HW space and both octet copies); for the top it is (φ₂-φ₁)². So ⟨x₀^G|Ψ_B⟩ gets its first contribution from the k_B-th order vector Ψ_(k_B), and we define

```text
X_B:=⟨x(φ)⊗x(φ)|Ψ_(k_B)⟩,        (68)
```

which turns out to be a product of differences φᵢ-φⱼ (below). The consistency check is N_B=8+0-(2+k_B)+2=8-k_B on all six branches.

**Table 11.** ε-orders for the six adjoint branches (U2 L1; the same at φ′). Nᵢ are the exponents of the five jet monomials (P³p₀₂,P³′p₀₂,P³p₀₃,P³′p₀₃,P³p₂₃) (the same on both sites); N_B is the exponent of ΣGNN; “first n” is the first order at which ⟨x(φ)⊗x(φ)|Ψₙ⟩≠0.

| branch | k_B | ⟨x₀^G∣Ψ_B⟩ | M₁M₂ | Nᵢ exponents | N_B | first n |
|---|---|---|---|---|---|---|
| (4,2,0) | 0 | ε² | ε² | 1,2,1,2,3 | 8 | 0 |
| (4,1,1) | 1 | ε³ | ε² | 1,1,2,2,3 | 7 | 1 (odd n only) |
| (3,3,0) | 1 | ε³ | ε² | 1,2,1,2,3 | 7 | 1 (odd n only) |
| (3,2,1) a, b | 2 | ε⁴ | ε² | 1,1,1,1,2 | 6 | 2 |
| (2,2,2) | 4 | ε⁶ | ε² | 1,1,1,1,1 | 4 | 4 (n=0..3 vanish to 10⁻⁷⁴) |

### 9.3. The leading polynomial identity

All orders t<N_B of ΣᵢⱼGˢᵢⱼ(εω̂)NᵢNⱼ cancel, with residuals 4.5·10⁻⁵¹ ((4,2,0)), 1.2·10⁻⁵⁵, 1.1·10⁻⁵⁵, 6.6·10⁻⁶⁶, 2.9·10⁻⁸⁴ at t=N_B-1, against 0.2–0.4 at t=N_B. The leading polynomial is L_B(ω̂)=Σ_(d+n+m=N_B)G^(s[d])ᵢⱼN⁽ⁿ⁾ᵢN⁽ᵐ⁾ⱼ, where Nᵢ⁽ⁿ⁾ is the εⁿ coefficient of Nᵢ, and the identity is

```text
L_B(ω̂₁,ω̂₁′)=ĉ ⟨φᵤ⊗e₁|Ψ₀⟩ m_B/X_B,   m_B=limε⁻²M₁M₂,  κ_B=1,        (69)
```

as a polynomial in (ω̂₁,ω̂₁′), coefficient by coefficient (relative coefficient residual 10⁻³⁹ to 10⁻⁷⁵ over six branches and two directions; the constant κ_B between the two sides is exactly 1). Each octet copy pairs with its own HW ket (different generating functions, same constant R_B=m_B/X_B). For the top state m_B=-(φ₃-φ₂)²θ₁₂² (Q₁=1; checked to 13 digits), so R_B(top)=-(φ₃-φ₂)²θ₁₂²/(φ₂-φ₁)²; in general Mₛ=P³AₛQ₁ with Aₛ=(ζ₃-ζ₂)e^(σθₛ)(θₛ-θₜ) and m_B=-(φ₃-φ₂)²θ₁₂²Q₁⁰(θ₁)Q₁⁰(θ₂) (not separately checked beyond the top). Which (d,n,m) contribute at order N_B (terms above 10⁻³⁰): (4,2,0): d=2..6, n,m=1..4 (largest (4,2,2),(3,2,3),(5,1,2),(4,1,3)); (4,1,1) and (3,3,0): d=1..5, n,m≤5 (largest (3,2,2),(4,1,2),(2,2,3)); octets d=0..4, n,m≤4; singlet only six triples (1,1,2),(2,1,1),(1,2,1),(0,2,2),(0,1,3),(0,3,1).

**Table 12.** The constant R_B=m_B/X_B of (69) at the two directions (U2 table L2(c)): best-fit value, and relative residual of the polynomial identity.

| branch | R_B(φ) | residual (φ) | R_B(φ′) | residual (φ′) |
|---|---|---|---|---|
| (4,2,0) | -1.1217793595 | 2.4·10⁻⁴³ | -0.0892224658 | 7.7·10⁻³⁹ |
| (4,1,1) | -4.8470775935 | 2.7·10⁻⁴⁷ | 0.1598811435 | 2.8·10⁻⁴⁶ |
| (3,3,0) | -1.1349116102 | 1.8·10⁻⁴⁷ | 0.0721026407 | 1.6·10⁻⁴⁵ |
| (3,2,1) a, b | -7.5969660751 | 2.9·10⁻⁵⁷, 7.3·10⁻⁵⁸ | -0.2001611630 | 2.1·10⁻⁵⁴, 3.4·10⁻⁵⁵ |
| (2,2,2) | -90.949002471 | 7.0·10⁻⁷⁵ | -0.7937987684 | 3.8·10⁻⁷² |

### 9.4. The reduced form: only Laurent data

Using Section 7 the leading identity can be written without any Q-function or normalisation:

```text
ΣᵢⱼGˢᵢⱼ(εω̂₁,εω̂₁′) νᵢ(θ₁)νⱼ(θ₂)=ε^(6-k_B) ĉ (⟨φᵤ⊗e₁|Ψ₀⟩)/X_B+O(ε^(7-k_B)),        (70)
```

verified coefficientwise to 10⁻⁵⁵…10⁻¹⁵⁷ with lower orders vanishing to 10⁻⁶²…10⁻¹⁶⁵; no M, no κ_B and no Q-normalisation enters. The ν-series are needed only to order 3 (top), 5 ((4,1,1), (3,3,0)), 4 (octets) and 2 (singlet), i.e. the Laurent data to order ε³,ε⁵,ε⁴,ε². Table 13 lists the (i,n) pairs that enter.

**Table 13.** Depth needed in the ν-form (U2 L3): i=0,…,4 labels (1,κ,π,κπ,εₛ), n the ε-order of νᵢ.

| branch | ν exponents | n needed per site |
|---|---|---|
| (4,2,0) | 0,1,0,1,2 | only the leading coefficients: (0,0),(1,1),(2,1),(3,2),(4,3) |
| (4,1,1) | 0,0,1,1,2 | i=1: 0,1,2,4; i=2: 1,2,4; i=3: 1,2,3,5; i=4: 2,3,4 |
| (3,3,0) | 0,1,0,1,2 | i=1: 1; i=2: 0,1,2,4; i=3: 1,2,3,5; i=4: 2,3,4 |
| (3,2,1) a, b | 0,0,0,0,1 | i=1..3: 0..4; i=4: 1..4 |
| (2,2,2) | 0,0,0,0,0 | i=1..4: 0..2 |

### 9.5. The CG dictionary and φ-independence

Solving the polynomial identity (69) or (70) for the HW ket (least squares over the polynomial coefficients, unknown in the HW space H_μ, then G-normalise and fix the sign by the first non-zero GT component) returns the exact CG slices: C₇=1 for the top (gl3 index; notebook ν=1), -0.70710678119 for (4,1,1) and (3,3,0), the octet pairs (C₂,C₄)=(0.208887289734,0.049311600013) and (-0.0797878448606,0.337986734608), and C₀=0.117851130198 for the singlet; the full 64-component vector equals the chain's Ψ₀ to 3·10⁻⁷⁷ (octets), 6·10⁻¹⁰⁰ ((4,1,1)), 3·10⁻¹⁰² (singlet). The result is *φ-independent*: repeating the whole procedure with φ′=(-2/5,3/8,7/11) changes N⁽ⁿ⁾, ĉ, the Gauss parameters, R_B, X_B and m_B (Table 12) but the recovered C_ν differ by ≤2·10⁻⁷⁷ (0 for the one-dimensional sectors, 10⁻²⁰² for the singlet).[^41]

The recipe for the CG generating function in “CG from Q” form is thus: write GF(ω̂)=Σ_(d+n+m=6-k_B)G^(s[d])ᵢⱼνᵢ⁽ⁿ⁾(θ₁)νⱼ⁽ᵐ⁾(θ₂)/ĉ and expand in the basis ĉ [e₁ᵀρ(n′(ω̂))]_ν (polynomials; the weight-(0,2,1) component -2a′g′ carries (φ₂-ω̂₁′)⁻² that ĉ does not clear, harmless because C₃=⟨e₃⊗e₁|Ψ₀⟩=0 on every HW ket of the sectors involved); the overall constant is fixed by the G-norm and the sign convention. The final formula needs exactly: (i) the kernel G split by degree; (ii) the one-site numbers ν(θₛ;ε)=(1,aₛ/3,π,aₛπ/3,εₛ), that is (Hₛ,aₛ) per site and branch as ε-series to the depths above; (iii) the Gauss parameters (34) and ĉ. The Q-polynomials are not needed beyond these numbers.

## 10. The ε-series of the Laurent data

**Punchline.** The Laurent data of the HW-limit branches have exact ε-series, polynomial in the direction φ (degree n at order εⁿ): linear for the top state, algebraic with a square root for (4,1,1) and (3,3,0), given to n=6 for the singlet and to n=5 over ℚ(√5) for the octets. They obey H₁+H₂=εΣ_(i≠j)φᵢμⱼ and bₛ=aₛ(Hₛ-aₛ)/3 order by order.

### 10.1. The series

In units θ₁₂=1 (at general θ₁₂, H⁽ⁿ⁾,a⁽ⁿ⁾∼θ₁₂ⁿ⁻¹ and b⁽ⁿ⁾∼θ₁₂ⁿ⁻², verified at two θ), with ε written e where terse:

- *Top, (4,2,0), exact to all orders:* H₁=4+ε(φ₁+2φ₂+3φ₃), a₁=a₂=3εφ₃ (κ=ζ₃, P³=e^(ζ₃u)), and Ψ_B=e₁⊗e₁ exactly for all ε. With θ₁₂: H₁⁽⁰⁾=4/θ₁₂.
- *(4,1,1)*, with x=ε(φ₂-φ₃): H₁=6+ε/2(2φ₁+5φ₂+5φ₃)+(√(1+x²/4)-1), a₁=6+3/2 ε(φ₂+φ₃)+3(√(1+x²/4)-1) (series to n=6: x²/8-x⁴/128+x⁶/1024, resp. 3x²/8,…; odd orders n≥3 vanish).
- *(3,3,0)*, with y=ε(φ₁-φ₂): H₁=6+3/2 ε(φ₁+φ₂+2φ₃)+(√(1+y²/4)-1), a₁=a₂=3εφ₃ exactly.
- *Singlet (2,2,2)*, with s₁=φ₁+φ₂+φ₃, Q=Σφᵢ²-Σ_(i<j)φᵢφⱼ, D=(φ₁-2φ₂+φ₃)(φ₁+φ₂-2φ₃)(2φ₁-φ₂-φ₃):

```text
H₁ =12+2s₁ε+Q/6 ε²+0·ε³-Q²/288 ε⁴+0·ε⁵+(deg6) ε⁶,
a₁ =12+s₁ε+Q/4 ε²+D/72 ε³-Q²/192 ε⁴-11/10368 DQ ε⁵+(deg6) ε⁶.        (71)
```

- *Octets (3,2,1) copy 0* (the other copy has √5→-√5; even orders are rational, odd orders carry √5), in terms of pᵢ=φᵢ: H₁⁽⁰⁾=9,

```text
H₁⁽¹⁾ =((15+√5)p₁+(20-2√5)p₂+(25+√5)p₃)/10,
H₁⁽²⁾ =(2(2p₁²-3p₁p₂-p₁p₃+3p₂²-3p₂p₃+2p₃²))/75,        (72)
```

H₁⁽³⁾=-√5(p₁-p₃)²(p₁-2p₂+p₃)/375, a₁⁽⁰⁾=(15-3√5)/2, a₁⁽¹⁾=[(3-√5)p₁+(9-3√5)p₂+(18+4√5)p₃]/10, b₁⁽⁰⁾=3√5, b₁⁽¹⁾=(12p₁+11p₂+22p₃)/5; orders up to n=5 in the data (Appendix D).

The structure is exact, not guessed: (i) the *global rules* H₁+H₂=εΣ_(i≠j)φᵢμⱼ (it holds at order 1 only) and b₁+b₂=ε²(φ₁φ₂μ₃+φ₁φ₃μ₂+φ₂φ₃μ₁); (ii) bₛ=aₛ(Hₛ-aₛ)/3 order by order; (iii) the six no-log conditions hold to 10⁻¹⁴⁸. Algebraic degree: on (4,1,1) and (3,3,0) the √(1+x²/4) is read off the series to n=6 (not proved) and suggests a 2×2 mixing, the HW-limit branch mixing with one neighbour; on the singlet and octets we did not identify closed algebraic forms.

### 10.2. How they were obtained

Two independent routes agree exactly. Agent W4 did exact Rayleigh–Schrödinger degenerate perturbation theory on the weight-μ space, with isolating operator H₁+a₁/7 (all six Laurent operators commute, and every series eigen-equation holds exactly): series to order 7 at φ and φ′ (data file record `E4_SERIES`), symbolic φ to order 3, and Ψₙ to order 4 symbolically. Agent U1 solved the six no-log conditions plus the two global rules order by order: that works (exact symbolic φ to n=6) when the Jacobian of the eight equations in six unknowns has rank 6, which is the case for (4,1,1), (3,3,0), (2,2,2). For the top state and the octets the Jacobian has rank 5 (kernel dimension 1), the implemented degenerate-solvability step failed (a cubic with no rational root), and U1 used instead the closed form from the Wronskian structure (top) and a polynomial fit in φ over 29 generic φ at 200–300 digits followed by exact ℚ(√5) reconstruction by PSLQ (octets, n≤5; n=6 would need ≥30 directions). All series were validated against 200 and 400-digit numerics at φ=(1/3,-2/7,5/11) to 10⁻⁶⁰ (n≤5).[^42]

Tables 20 and 21 (Appendix D) list the series at φ. A few exact entries as a sanity check: top H₁=140/31+(260/231)ε, a₁=(15/11)ε, b_(1,2)=±(700/341)ε-(25/231)ε²; (4,1,1) at φ: H₁⁽¹⁾=349/462, a₁⁽¹⁾=39/154, H₁⁽²⁾=100719/1660120; singlet H₁⁽¹⁾=232/231 and a₁⁽¹⁾=116/231 (2s₁ and s₁ at φ, in units θ₁₂=31/35, since 2(1/3-2/7+5/11)=232/231).

### 10.3. The perturbative content: Ψ₁ and X_B

The exact perturbation theory also gives the irreducible content of the first-order vector. For (4,1,1): Ψ₁=31/140(φ₃-φ₂)(1,1) lies entirely in the (4,2,0) irrep (its G-norm² is (θ₁₂²/8)(φ₂-φ₃)² while Ψ₀=(1,-1) has G-norm² 2); (3,3,0) is the same with φ₂-φ₁, entirely in (4,2,0); the singlet Ψ₁ lies entirely in the (3,2,1) octet irrep (G-norm² (15376/1225)Q, no (3,3,0),(4,1,1),(4,2,0) component); the octets have components in all four irreps, proportional to (φ₂-φ₃) and (φ₁-φ₂) combinations. The leading coefficient of the reference overlap, closed-form X_B=⟨x(φ)⊗x(φ)|Ψ_(k_B)⟩, with Ψ₀ of G-norm 1:

```text
X_((4,2,0)) =(φ₁-φ₂)²,   X_((4,1,1))=-θ₁₂/2(φ₁-φ₂)(φ₁-φ₃)(φ₂-φ₃),
X_((3,3,0)) =-θ₁₂/2(φ₁-φ₂)³,   X_((2,2,2))=(θ₁₂/2 )⁴[(φ₁-φ₂)(φ₁-φ₃)(φ₂-φ₃)]²,
X_octet =-(5766(√5-3))/6125(φ₁-φ₂)²(φ₁-φ₃)(φ₂-φ₃)  and its conjugate,        (73)
```

(for (4,1,1), (3,3,0) the formulas are for Ψ₀ of G-norm² 2, divided by √g₀ in ⟨x₀^G|Ψ_B⟩=ε^(2+k_B)X_B/√g₀; for the octet g₀=1080(7-3√5), for the singlet g₀=288; the data file has E5\_KB with the exact strings). These are the closed forms of U2's numerical X_B up to normalisation: e.g. at φ the numbers are 0.383220,-0.017394,-0.074288,0.0013930,6.9931·10⁻⁶.

The generating functions in the data file are ĉ⟨φᵤ⊗e₁|Ψ_B⟩=Σ_ν W_ν s[ν] with W_ν=ĉ[e₁ρ(n′)]_ν polynomial in (ω̂₁,ω̂₁′) (ν=1..6,8; W₇ has a simple pole 1/(ω̂₁′-φ₂), which never carries a slice) and s[ν]=Ψ(ν,1), exactly the CG coefficients at order 0 and exact linear or quadratic forms in φ at orders 1,2 (record `E5_GF` of `untwisted-exact-210.txt`).

### 10.4. Worked numbers for adj⊗adj

The complete exact CG table of adj⊗adj is Table 22 (rows with Λ=(2,1,0)): C₁=1 for (4,2,0); C₂=-1/√2 for (4,1,1); C₃=-1/√2 for (3,3,0); C₈=2/√288=√2/12 for the singlet (2,2,2); and for the octet copy 1,

```text
C₄=(21-9√5)/(√(7560-3240√5)),   C₅=(-3+3√5)/(√(7560-3240√5)),        (74)
```

copy 2 with √5→-√5 (indices in notebook order: ν=4,5 are the two weight-(1,1,1) states). Their Laurent data at ζ=0 and the first corrections are in Tables 23 and 20. What the Q-side formula needs for each branch is: Gˢ, the series of (Hₛ,aₛ) to the depths of Table 13, the Gauss parameters and ĉ; the output is the polynomial in (ω̂₁,ω̂₁′) of which C_ν are the coefficients up to the overall constant m_B/X_B (Table 12), fixed by the G-norm. A single end-to-end script reproducing the whole table from the pieces (G, the series, the Gauss parameters) was not written; it is item F3 of the round-5 plan.

## 11. The untwisted model

**Punchline.** At ζ=0 the Laurent data of a HW state are Casimir eigenvalues: H₁ from C₂, a₁-a₂ from the cubic Casimir, and bₛ=0 except in multiplicity-2 sectors. The untwisted oper is Fuchsian at infinity with three polynomial solutions of degrees (μ₃,μ₂+1,μ₁+2), the Weyl-shifted weight; it is a ₃F₂ equation only on branches with a consecutive exponent pair. The monic twisted Q-functions have poles in ε; the finite objects are the Plücker monomials and the Laurent data.

### 11.1. Laurent data from Casimirs

On the HW kets of weight μ in V_Λ⊗V_Λ at ζ=0 all six Laurent operators act within the HW space, and we can compute their eigenvalues in closed form. With T=|Λ| and the Casimirs of the GT normalisation EᵢₖEₖᵢ and EᵢⱼEⱼₖEₖᵢ,

```text
C₂(l)=Σᵢlᵢ(lᵢ+4-2i),   C₃(l)=l₁³+l₂³+l₃³+4l₁²+l₂²-2l₃²-l₁l₂-l₁l₃-l₂l₃+4l₁-2l₂-2l₃,        (75)
```

one has H₁=-H₂=η/θ₁₂, a₁=-a₂=α/θ₁₂, b₁=-b₂=β/θ₁₂² with η,α,β independent of θ (verified at three θ points), and

```text
η=T²-½ [C₂(μ)-2C₂(Λ)],  
α=(4T(1+T)²-3(1+T)C₂(μ)+C₃(μ))/6=((μ₁+2)(μ₂+1)μ₃)/2,        (76)
```

(the first from ΣᵢₖE¹ᵢₖE²ₖᵢ=½(C₂^tot-C₂⁽¹⁾-C₂⁽²⁾), as in (7); the second found by an exact fit over 30 multiplicity-one sectors and then verified on every sector, and agreeing with U1's independently found product form over 23 sectors). The product is half the product of the entries of the Weyl-shifted weight (μ₁+2,μ₂+1,μ₃): α=0 iff μ₃=0, α=η for μ₃>0 except when λ₂≥2 where α-η=4,6,6,10 for [2,2,0] (4,2,2), [3,2,0] (6,2,2), (4,4,2), (4,3,3). For the adjoint η=4,6,6,9,12 and α=0,6,0,15/2,12 on (4,2,0),(4,1,1),(3,3,0),(3,2,1),(2,2,2). In terms of the pair (μ₃, cubic Casimir) this is a “determinant” of the shifted weight: a cubic Casimir quantity.

In *multiplicity-1* sectors β=0. In the multiplicity-2 sectors (the only ones that occur) the two copies are the eigenvectors of the 2×2 matrix of S3c[1,2]=a₁ (Table 7): α_±=(tr∓√disc)/2, β=±√disc with

```text
tr=η+2c₂,  disc=η²-4c₂²  (7 of the 8 multiplicity-2 sectors),        (77)
```

giving disc=45 for [2,1,0] (3,2,1); 105 and 192 for [3,1,0] (5,2,1), (4,3,1); 105 and 192 for [3,2,0] (5,4,1), (5,3,2); and β_copy=α_copy(η-α_copy)/c₂ (the no-log relation). **There is one exception: [3,2,0] (5,3,2)**, where tr=56=η+2c₂+18 and disc=192 (the prediction would have been 228); there the no-log relation in this form fails. No closed form for β(Λ,μ) is known: it is not a function of C₂,C₃ alone but needs the multiplicity-space operator. For the adjoint octet α_±=(15∓3√5)/2, β=±3√5, with eigenvalue equation λ²-15λ+45=0. At the standard point this gives a₁=4.68085… and 12.25463… for the two copies.

Two sum rules hold as operator identities on all of V⊗V at ζ=0 (exactly, five or six representations, 37 sectors):

```text
b₁+b₂=0,   a₁+a₂+θ₁b₁+θ₂b₂=0,   (H₁+H₂=0).        (78)
```

Hence S₂=O(u⁻²) and S₃=O(u⁻³) as u→∞ (no constants at ζ=0): the oper is regular singular at infinity.[^43] The exact-data file lists all of it (Appendix E).

### 11.2. Fuchsian at infinity: three polynomial solutions

The untwisted dual oper ∂³-S₁∂²+S₂∂-S₃ has, at θₛ, exponents {0,λ₂+1,λ₁+2} (Section 2.2), and at infinity solutions growing as uᵏ, where for the adjoint k solves k³-9k²+(14+η)k-2α=0. The indicial roots at infinity are, on all 45 branches,

```text
k=(μ₃, μ₂+1, μ₁+2),        (79)
```

the Weyl-shifted weight (e.g. top (0,3,6), (4,1,1): (1,2,6), singlet (2,3,4)). The Fuchs relation reads: twice the sum of the three exponents at one θₛ, minus Σk, equals 3 (for the adjoint 12-Σk=12-(|μ|+3)=3). No logarithms occur (all exponents differ by integers, the local monodromy at θₛ is trivial by the no-log conditions, and at infinity the exponents are integers), so *all three* solutions are polynomials, of degrees (μ₃,μ₂+1,μ₁+2), with Wronskian proportional to Π³ and apparent singularities at θ₁,θ₂. The exact dimension of polynomial solutions of degree ≤d is, for the top, 1,1,1,2,2,2,3 (d=0,…,6), for (4,1,1) 0,1,2,2,2,2,3, and always reaches 3. This refutes the expectation recorded in the round-4 plan that fewer than three polynomial solutions exist at ζ=0 (the Wronskians of the twisted Q's degenerate as ζ→0, which suggested a degenerate solution space; Section 11.4 explains where the degeneracy really sits).

The exact echelon bases in x=(u-θ₂)/θ₁₂ (site 2 at x=0, site 1 at x=1), all verified to solve the oper in sympy:

```text
(4,2,0):  {1, ½ x²(2x-3), ¼ x⁴(2x-3)²}={1,F,F²},  F=∫Π (a Sym²),
(4,1,1):  {½(2x-1), x², ½ x⁴(2x²-6x+5)},
(3,3,0):  {1, x²(x-1)², ¼ x²(4x³-10x+5)},
(2,2,2):  {¼(2x-1)², ⅛(2x-1)(4x²+2x-1), x⁴},
(3,2,1)_copy 0:  {1/10(10x-5+√5), ½ x²(2x-3-√5), ½ x⁴(2x-5+√5)},        (80)
```

(copy 1: √5→-√5; the copies are Galois conjugate). For the top state the basis is {1,F,F²} with F=∫Π (note the twisted top-state Q-functions are e^(ζ₃u) and e^(ζ₂u)(Π+Π′/d+Π″/d²) with d=ζ₃-ζ₂).[^44]

### 11.3. The ₃F₂ criterion

A natural hope was that the untwisted oper is a hypergeometric ₃F₂ equation, so that its polynomial solutions would be hypergeometric polynomials. This is true only sometimes. W4 compared the oper after Möbius maps and gauge factors with the ₃F₂ equation for all role assignments of the three singular points θ₂,θ₁,∞ and found: **the oper is a ₃F₂ equation if and only if some singular point carries a consecutive pair of exponents**, which happens on exactly 20 of the 45 branches. At θₛ the exponents {0,λ₂+1,λ₁+2} contain a consecutive pair iff λ₁=λ₂ (pair {λ+1,λ+2}: [1,1,0],[2,2,0], all five branches); at infinity iff μ₂=μ₃ or μ₁=μ₂. The parameters are b=(-λ₂,-λ₁-1) with witness roles s=0 at θ₂, s=1 at infinity, s=∞ at θ₁ and gauge (1-x)^… (for λ₁=λ₂: s=x, gauge (1-x)^(λ+1)), the parameters aⱼ being tabulated by witness. These are degenerate integer parameters, reducible, so “hypergeometric” is not the natural frame: the polynomial bases above are the closed forms. On the other 25 branches (all with μ₁>μ₂>μ₃ and λ₁≠λ₂, including every multiplicity-2 sector) there is no consecutive pair at any of the three points, hence no ₃F₂ form. U1's adjoint verdict (only (4,1,1), (3,3,0), (2,2,2), with X=(u-θ₂)/(u-θ₁) and gauge (u-θ₁)ᶜ, c=2,5,3 or 4) agrees with W4's.[^45]

The untwisted *direct* oper (our Πᵖ∘adjoint∘Π⁻ᵖ) has polynomial solutions of degrees {2λ₁+2-μ₃, 2λ₁+1-μ₂, 2λ₁-μ₁} on 43 of the 45 branches; the exceptions are [1,1,0] (2,2,0) (degrees {0,2,4}) and [2,2,0] (4,4,0) ({0,2,6}).

### 11.4. Why the twisted monic Q's have poles in ε

We expected the untwisting limit of the twisted Q's to be the three untwisted polynomials. The monic pⱼ(u;ε) have *no* coefficientwise limit: the coefficient cₘ of pⱼ behaves as $\epsilon^{-\min(\mu_j-\mu_3,\ \mu_j-m)}$ (top: P¹ coefficients ε^(-4..0), P² ε^(-2..0), P³ finite; (4,1,1): P¹ at ε⁻³; (3,3,0): P¹,P² at ε⁻³; (3,2,1): P¹ at ε⁻², P² at ε⁻¹; singlet finite). Expanding in the exact untwisted basis vₐ(ε) (solutions of the twisted oper with the same jets at u=0), Pʲ=ΣₐMⱼₐ(ε)vₐ(ε), the matrix M has pole order μⱼ-μ₃ on v₃ and O(ε^(-(Nⱼ-1))) on v₂,v₁, with detM=c_W Π(0)³∼ε³ exactly; the normalised rows (1,ρ₂,ρ₁)ⱼ have ρ=εr⁽¹⁾+⋯ with the *first-order direction common to all j*, differing only at ε² (octets, singlet) to ε⁴ (top: R₁-R₃=O(ε⁴), R₂-R₃=O(ε³)). So the degeneracy sits in the monic normalisation, not in the solution space: the limiting 3-plane is the span of the untwisted polynomials and the ζ-jet basis is the Smith (saturation) basis of M; its elementary-divisor valuations were not computed, and no closed form for r⁽¹⁾ was found. Limits of ε-differences and derivatives: the leading Laurent coefficients of Pʲ solve the untwisted oper, the lower ones solve D₀p⁽ⁿ⁺¹⁾=-D₁p⁽ⁿ⁾-⋯ (inhomogeneous), consistent with the numerics. The finite objects are the Plücker monomials and the Laurent data, which is why the formulation of Section 7 is the right one.

## 12. Surprises and corrections along the way

**Punchline.** About two dozen of our own expectations were wrong or had to be sharpened. We list them candidly, with where each was caught, because the pattern (conflating |Λ| with c₂, expecting degeneracy where there was none, expecting too many monomials) tells us where to be careful with λ₁≥3 and λ₃>0.

**Table 14.** Corrections to the plans, in the order they were found. “Round” is the plan round that carried the wrong expectation; the correcting agent is in the last column.

| # | expected | found | by |
|---|---|---|---|
| 1 | H₁=+[…] (plan, round 1) | sign is -: H₁=-[ΣζᵢEᵢᵢ+ΣE¹ᵢₖE²ₖᵢ/θ₁₂]+∣Λ∣σ+∣Λ∣²/θ₁₂, H₁+H₂=-W+2∣Λ∣σ | W1 |
| 2 | x₀ (XXX limit) not parallel to the Gaudin x₀^G | exactly parallel (10⁻³¹), x₀∼h per site, x linear in z at all orders | P2, D2 |
| 3 | twelve P³⁽ᶜ⁾pₖₗ monomials per site | five (adjoint); in general C(λ₁+2,2)-C(λ₂+1,2) with the lone factor always P³ | D1, D2, P4 |
| 4 | M_vac=5/2 h²W′P³ | 2h²W′P³ (the weight-2 minor contains an n=1 piece fixed by minor(-1,0)=0) | D2 |
| 5 | S₂ double pole =∣Λ∣ (or Casimir/2) | c₂=λ₂(λ₁+1) | P3, then P4 |
| 6 | aₛ^vac=∣Λ∣ω₃, f(S3c)=∣Λ∣ζ₃ | c₂ω₃, c₂ζ₃ (equal on [2,1,0], [3,1,0]) | P4, W3 |
| 7 | the direct oper is the formal adjoint | the adjoint has no e^(-ζu)×polynomial solutions; direct =Πᵖ∘adjoint∘Π⁻ᵖ, p=λ₁-λ₃ | P1, P2, P3 |
| 8 | dual Wronskian F′=Π for [3,2,0] | F′=Π^(λ₂)=Π²; direct exponents {0,2,5}, not {0,3,5} | P4 |
| 9 | c_H has 1/θ₁₂ coefficient 2λ₁ | 2e₂(Λ) (2λ₁λ₂ for λ₃=0) | W2 |
| 10 | the [3,1,0] fourth localiser is a single product, S3c[s,1]=c_Λ S3c(H-S3c) | two quadratic monomials plus a W term; for [3,2,0] no c_Λ at all, γ changed 2→1 | W2, W3 |
| 11 | ℓ₁=1 risk: no localiser inside D₂ | ℓ₁=3 already inside D₁; every Laurent coefficient localises | W1 |
| 12 | ω₁′ fixed control fails in Gaudin as in XXX | passes at finite ω₂, fails at ω₂=∞ | W1, W2 |
| 13 | fixed-coefficient orthogonality among brackets | the relation has state-dependent Δ (Lagrange identity); six residue columns | P1 |
| 14 | ζ-independent kernel with vacuum monomial M=Q₁P³ | M=P³p_(0,λ₂+1) and the scalar c are needed; direct-Q₁ bases carry a hidden σ | D1 |
| 15 | naive t→∞ at fixed a gives the Gaudin kernel | gives a rank-one state-independent (1-a)³(1-a′)³; wrong limit | D2 |
| 16 | the SoV covectors degenerate at h→0 in the companion frame | they are independent there; the degeneracy lives in the diagonal frame (P⁻¹_phys=O(h)) | D3 |
| 17 | one polynomial Gaudin B_G(u) | none; solvable hierarchy, no polynomial of low degree | D3 |
| 18 | cleaning needed for (5,3,2) etc. | no cleaning anywhere (K/S=O(ε) on all branches) | W1–W3 |
| 19 | top-state sum rule =M^topM^top | extra factor R(a)R(a′), R=(1-ζ₂a)^(λ₁-λ₂)(1-ζ₃a)^(λ₁) | G1 |
| 20 | λ₂=0: S₂ double pole non-scalar | scalar 0; one-row reps are the master kernels | G1 |
| 21 | ≤2 polynomial solutions of the untwisted oper | all three are polynomials, degrees (μ₃,μ₂+1,μ₁+2) | U1, W4 |
| 22 | the untwisted oper is ₃F₂ | only with a consecutive exponent pair (20 of 45 branches) | U1, W4 |
| 23 | the monic twisted Q's have untwisted limits | poles in ε; finite objects are Plücker monomials and Laurent data | U1 |
| 24 | ĉ clears all denominators of e₁ρ(n′) | not the weight-(0,2,1) component -2a′g′ with (φ₂-ω̂₁′)⁻² (harmless) | U2 |
| 25 | detN rank 5 for all A≠B | rank 4 iff Δμ has one positive and two negative components | P3 |
| 26 | gauge power fixed by λ₃=0 data | p=λ₁-λ₃ vs λ₁ not separated (all data λ₃=0) | P3, P4 |
| 27 | octet copies distinguished only at higher order | they split already at ε⁰ by aₛ⁽⁰⁾ (while Hₛ⁽⁰⁾ coincide) | U2 |

A few of these deserve a sentence more than the table gives. The first plan said that with the XXX-like ansatz the Gaudin formula would involve twelve monomials and M_vac=5/2 h²W′P³ (items 3, 4); both numbers came from the heuristic that the XXX ladder points θₛ+mh become the order-3 jet, and in fact the true count follows from the Frobenius structure of the oper. The *cleaning question* (item 18) was Paul's own: whether the per-sector, per-copy cleaning of the single-vacuum approach can be avoided; the plan's careful list of expected failures ((5,3,2), multiplicity 3) is a good example of a plan that was too pessimistic. The conflation of |Λ| and c₂ (items 5, 6) is the one to remember: on [2,1,0] and [3,1,0] one has λ₂=1 and hence c₂=λ₁+1=|Λ|, so the two quantities coincide and only [3,2,0] could tell them apart. Finally, item 21 is a good reminder that exponents integer-spaced with trivial monodromy mean polynomial solutions: the “degeneracy” of the Wronskians at ζ=0 does not reduce the solution space.

## 13. Outlook: the round-5 plan

**Punchline.** The adjoint works end to end. Round 5 aims to turn “the adjoint works” into “the formula”: a closed form for the master kernel R^(λ₁), the Laurent-data form and untwisting limit for [3,1,0] and [3,2,0], a derivation of the ε-series, the norm of the untwisted eigenvectors from the Q side, and the write-up. Nothing below is implemented; Paul asked for the plan, not for work.[^46]

### 13.1. WP-A: the master kernel in closed form

The goal is a four-index closed form of R^(λ₁)_((kl),(k′l′))(x,y), x=at, y=a′t, valid for all λ₁; then every Λ follows by the block theorem and θ⊗θ′. Data in hand: exact R for λ₁≤4 ([4,0,0]: 225 entries), the first-order formula, the single-index Γ-form, the extreme-hop coefficients, the exchange rule, homogeneity, the rank pattern 1,2,4,8 and the top-state sum rule. Three attempts are planned: (A1) a fit-guided Pochhammer ansatz generalising the Γ-form, validated on a new representation, [5,0,0] (1225 states), done at 40 digits through the ν-form, which needs only Laurent data and no Q-solving, so that Section 7 pays off computationally; (A2) a structural derivation of the jet basis for one-row representations V_((λ₁,0,0))=Sym^(λ₁)(ℂ³) in the polynomial realisation Eᵢⱼ=xᵢ∂ⱼ, where the SoV wavefunction rule becomes explicit differential or integral operators (start with λ₁=1,2); (A3) the auxiliary sl₂: is R covariant under a twisted shift (a,a′,t)→(a+s,a′+s,t′)? If so R is fixed by its value at a=0, that is by the one-site data plus finitely many coefficients; also whether R^(λ₁) is the spin-(λ₁+1)/2 fusion of the λ₁=1 R-matrix form. Acceptance: R^(λ₁) reproduces all exact kernels for λ₁≤4 and predicts [5,0,0] (and [5,1,0]) on ≥200 states to ≥30 digits.

### 13.2. WP-B: Laurent-data form and the limit for [3,1,0], [3,2,0]

(B1) Derive νₖₗ for general (λ₁,λ₂) from the Frobenius recursion at θₛ (exponents {0,λ₂+1,λ₁+2}; the order-(λ₁+1) jet lies in a 2-plane with coefficients rational in (Hₛ,aₛ,bₛ)), prove that P³Dₛ cancels between the sides, and verify on the 225 states of [3,1,0] and [3,2,0] to 10⁻⁵⁰. (B2) Run U2's pipeline for both representations: k_B per branch (10 each) and the depth; test whether k_B follows a weight rule (adjoint 0,1,1,2,2,4) and the exact CG tables of the data files. (B3) The exact ε-series of the Laurent data for the 10+10 branches (duality check: [3,2,0]'s series are the dual images of [3,1,0]'s). Acceptance: all 20 CG tables from the Q side to ≥30 digits, φ-independent.

### 13.3. WP-C: the mechanism of the ε-series

(C1) Implement the degenerate solver: Fredholm alternative at each order for the rank-5 cases (top state, octets), symbolic in φ, to n=8, from the Bethe side alone. (C2) Identify the 2×2 mixings behind the algebraic (4,1,1), (3,3,0) forms and the algebraic degrees on the others; conjecture for Hₛ(ζ) on a HW-limit branch as eigenvalues of a finite matrix built from the twisted gl₃ Casimirs. (C3) Closed form via Casimirs on the chain side: the Laurent operators on the μ-isotypic component are expressible through C₂, C₃, ΣζₖEₖₖ and products; write them in that algebra and diagonalise symbolically on each isotypic component (small matrices since the twist mixes only within weight spaces). Acceptance: an algorithm or formula producing the series of every branch of [2,1,0], [3,1,0], [3,2,0] exactly.

### 13.4. WP-D: the norm from the Q side

The CG formula fixes the generating function up to one constant per branch, the G-norm. Classical results say the norm of a Gaudin Bethe vector is the Hessian of the Yang–Yang function at the Bethe roots (Reshetikhin–Varchenko, Mukhin–Varchenko); here the Bethe equations are the no-log conditions. Hypothesis: ⟨Ψ|Ψ⟩_G/⟨x₀|Ψ⟩² is proportional to the determinant of the Jacobian of the six no-log conditions with respect to the Laurent data, or to the equal-twist detN suitably normalised (D3's observation that Q_num^eqM₁M₂/detN_B is state-independent is norm-like). Test numerically on the 64 states, then take the untwisting limit; if none of the candidates is state-independent the true Yang–Yang Hessian in terms of Q-zeros is needed, a different formalism, and the plan says to flag it and stop.

### 13.5. WP-E, WP-F: the write-up and housekeeping

The notes you are reading are E1. E2 is the paper: `paper/gaudinCG.tex` is a gl(4)-rectangular stub, to be replaced by the sl(3) Gaudin CG paper (recommended) or a new file; proposed structure: introduction; the two-site model, cdet, Laurent data, opers; vacua, frames and the two-vacua family (no cleaning); the jet formula (basis, K tensor, block theorem, one-site kernel, master kernel); the overlap as a function of Laurent data; the untwisting limit; worked example adj⊗adj; outlook; appendices on conventions, kernels, series. Housekeeping: rerun the four Gaudin notebooks in VS Code to store native outputs (the agents stored none), gitignore `*.pkl` (≈230 MB in `Paul/Python/gaudin_sl3/`) before committing, promote the notebooks to `Clean/` and run the docs agent, update `CLAUDE.md`'s outdated goal statement, a one-file demonstration `cg_from_Q.py`, and the backlog: [4,2,0] (dim 27, multiplicity 3) as the last cleaning check, the (5,3,2) exception of the multiplicity-2 law and β(Λ,μ), the second null vector of the rank-4 pairs, the λ₃>0 dependence of the gauge rule.

The decisions that are Paul's: the paper file; the priority between WP-A (the more striking result) and WP-C/D (the more complete one, no chain input); whether to spend a round on [4,2,0] before writing; and the update of `CLAUDE.md`.

## Appendix A. Conventions and normalisations

This appendix reconciles the notations of the different agents and fixes every factor. The main text uses one convention throughout (the canonical one below); each source file is mapped to it.

### A.1. The kernel normalisations

**Table 15.** Kernel normalisations in the different sources and their relations. c is the frame scalar (15); M=P³p_(0,λ₂+1).

| symbol | where | definition and relation |
|---|---|---|
| K | G1, P4, main text | c Q_numM₁M₂=ΣK m m, mₖₗ=P³pₖₗ, M exactly (no extra factor) |
| G | D1, D2, D3 | 4c Q_numM₁M₂=ΣG N N (D2's “4c”); for the adjoint G=4K, because M_vac=2h²P³p₀₂ per site in the XXX limit and 2×2=4 |
| Gˢ | U2 | Gˢ:=G/4=K (c Q_numM₁M₂=ΣGˢNN; U2 states “G of the spec = D2's printed G divided by 4”) |
| Θ | G1, D3 | Θ=K/(ω₁ω₁′)^(2λ₁-λ₂), a polynomial in (a,a′,t)=(1/ω₁,1/ω₁′,1/θ₁₂); adjoint: Θ=G/(4ω₁³ω₁′³) (D3, exact in 25/25 entries) |
| Q_num | all | ⟨Φ_∞∣Ψ⟩/⟨x₀^G∣Ψ⟩ in the diagonal frame; companion-frame version Q_num^comp=Q_num/(κ^Gκ^(G′))=c Q_num/(ω₁ω₁′)^(2λ₁-λ₂) (D1, D3) |
| T(B) | D1 | ⟨Φ_∞∣Ψ⟩D_B(ω₀)/⟨Ω₀⊗Ω₀∣Ψ⟩=c_T(ζ) Q_numM₁M₂ |
| ν | U2 | Nᵢ/M=(1,κ,π,κπ,εₛ); the main text writes ν for this vector, and ν(u)=Π(u) appears only in the XXX B-operator discussion (Section 6.1) |
| D₂ | W2, W3 | the quadratic localising algebra; not to be confused with D2 (agent) or D_B (detN) |
| γ₂,γ₃,d₃ | this text | Frobenius coefficients of f₀,f₂ (U2 calls γ₂,γ₃ “c₂,c₃”; renamed to avoid the clash with the double pole c₂) |

The relation Θ^adj=G/(4ω³ω′³) was found by D3 to be exact; G1 verified symbolically in 25/25 entries that Θ defined from K equals it, so no further factor hides in the conventions. The “4 of D2/D3” is a convention, not physics: with M=P³p_(0,λ₂+1) exactly one uses K=G/4 and no 4 appears anywhere.

### A.2. Frames and bases

- *Frames.* Diagonal frame: states of the diagonal-twist Hamiltonians; companion frame: Ψ_c=(P_phys⁻¹⊗P_phys⁻¹)|_(z=ζ)Ψ_d, x₀ᶜ=e₁⊗e₁ (D3), Φ^∞_c=Ω_c^∞[ω₁]⊗Ω_c^∞[ω₁′]. x₀^G in the diagonal frame is x(ζ)⊗x(ζ) and ⟨x₀ᶜ|Ψ_c⟩=⟨x₀ᵈ|ψ_d⟩.
- *Basis orders.* Notebook (W-agents, data files): GT patterns descending, index 1 or 0 is the HW ket; `gl3_reps` (P, D, G, U agents): lowest weight first, HW ket last. Adjoint map NB→gl3: (7,5,6,4,2,1,3,0) (0-based). Weights in notebook order for the adjoint: (2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,1,1),(1,0,2),(0,2,1),(0,1,2); the metric diagonal G=(1,1,1,2,6,9,4,9). U2's slices (`gl3` index ν): top C₇, (4,1,1) C₅, (3,3,0) C₆, octets (C₂,C₄), singlet C₀; W4's (notebook, 1-based): C₁,C₂,C₃,(C₅,C₄) resp. (C₄,C₅), C₈; Table 22 uses W4's.
- *Signs.* HW kets are G-normalised with the first non-zero component, in descending GT order, positive (D1, U2); W4's canonical HW kets are primitive integer vectors (first non-zero component positive). For (4,4,2) of [3,2,0] the exact slice has an exactly vanishing component and the chain branch's sign convention differs; the CG cells align signs and use norm-based errors (W3).
- *Copies.* In a multiplicity-2 sector copy 1 is the eigenvector of S3c[1,2] with the smaller eigenvalue (α₁⁽⁰⁾); U1's “copy 0” is W4's copy 1 and √5→-√5 exchanges the copies.
- *θ₁₂ units.* U1's series are in units θ₁₂=1, with H⁽ⁿ⁾,a⁽ⁿ⁾∼θ₁₂ⁿ⁻¹, b⁽ⁿ⁾∼θ₁₂ⁿ⁻²; W4's data are θ-scaled, H=η/θ₁₂, a=α/θ₁₂, b=β/θ₁₂², with θ=(2/7,-3/5), θ₁₂=31/35 for the numerical tables and series (so “140/31” is 4/θ₁₂).
- *Kernel pair orders.* G1: lexicographic pairs (kl), 0≤k<l≤λ₁+1, l≥λ₂+1. D2/D1: N=(P³p₀₂,P³′p₀₂,P³p₀₃,P³′p₀₃,P³p₂₃); the map to G1's order (02),(03),(12),(13),(23) is 1→(02), 2→(12), 3→(03), 4→(13), 5→(23) because P³′p₀ₗ=P³p₁ₗ.
- *c for general Λ:* (15); wrong exponents tested and excluded (Section 1.4).

### A.3. Numerical points and grids

- θ=(2/7,-3/5), θ₁₂=31/35; θ variants for kernel laws: (1/5,-4/9), (3/7,-1/4), (1/2,-1/3), (-1/7,-6/5), (4/5,-2/3) (θ₁₂=29/45,19/28,5/6,37/35,22/15); further θ₁₂ values 7/9,13/8,-9/7,-17/10,11/6,3/4,-5/3,19/13,-7/11,9/5,-13/9,25/14,-19/10,15/11,-11/7 (G1, 12–20 values); pairs with equal θ₁₂: (2/7,-3/5),(-1/3,-128/105),(5/4,51/140).
- Twist points: P=(3/2,5/4,2/3), P′=(2,5/3,-1/2), X₁=(7/5,-3/4,1/6), X₂=(-4/3,2/9,5/7), X₃=(1/3,-2/7,5/11); further X₄..X₆ and 26 real-spectrum random rationals (D1); some triples give a complex spectrum of S₂+πS₃ and are skipped.
- Vacuum twists: ω=(5/2,7/3,3/5), ω′=(7/4,-2/5,9/7); ω₀=(5/2,7/3,3/5) for T(B); (ω₁,ω₁′): base (5/2,7/4), second (-2/5,3), a 3×3 grid ω₁∈{5/2,-2/5,9/7}×ω₁′∈{7/4,3,-5/3}, a 9-value grid {5/2,-2/5,9/7,3,-5/3,7/4,11/3,-7/2,13/5} (72 ordered pairs), and G1's w′ list up to 23/9,-3/2,10/3.
- Regime: φ=(1/3,-2/7,5/11), control φ′=(-2/5,3/8,7/11); ω̂₁∈{3/4,-5/3,1/2,-7/3}, ω̂₁′∈{-12/5,7/2,1/8,5/4} (the 24-point Stage-D grid is {3/4,-5/3}×{-1/5,9/7,5/3} for ω̂ and {-12/5,7/2}×{1/8,-3/4} for ω̂′); ε∈{1,2,4}·10⁻⁶ (chain-side, three values) and εₖ=10⁻⁶2⁻ᵏ, k=0..11 (series fits, K=12; the extended set K=24, 400 digits, reliable to n=18).
- Generic points for operator identities: u₀=3/7, u₁=-2/9 (note: Laurent fits at θ₁=3/7=u₀ hit a pole).
- Precision: chain 200 digits; Q side 60 (finite twist), 80 (fits, kernels), 100–120 (regime), 200 (series), 400 (extended series); fixed-point 320-bit integer QR in D1; exact Fractions in D3, G1 (one-site).

### A.4. Vacuum eigenvalue functions

The exact eigenvalue functions of the vacuum Ω_(ζω)⊗Ω_(ζω) (exact rationals; xₛ=u-θₛ), with S₁=e₁(ω)+|Λ|(1/x₁+1/x₂), e₁(ω)=163/30 for ω and 369/140 for ω′, and the constants e₂(ω)=131/15, e₃(ω)=7/2 (e₂(ω′)=29/28, e₃(ω′)=-9/10). The table lists the Laurent data of S₂,S₃ (the double poles are c₂ and aₛ=c₂ω₃).

**Table 16.** Vacuum Laurent data (exact). Hₛ is the residue of S₂ at θₛ, bₛ that of S₃, a the common double pole of S₃. Sources: P1, P3, P4 section 0.

| Λ | twist | H₁ | H₂ | b₁ | b₂ | a |
|---|---|---|---|---|---|---|
| [2,1,0] | ω | 12539/930 | 4139/930 | 2173/310 | 493/310 | 9/5 |
| [2,1,0] | ω′ | 40463/4340 | 1263/4340 | 30501/4340 | -19899/4340 | 27/7 |
| [3,1,0] | ω | 5789/310 | 1589/310 | 3027/310 | 507/310 | 12/5 |
| [3,1,0] | ω′ | 54107/4340 | -4693/4340 | 40869/4340 | -34731/4340 | 36/7 |
| [3,2,0] | ω | 885/31 | 45/31 | 2376/155 | -144/155 | 24/5 |
| [3,2,0] | ω′ | 48341/2170 | -10459/2170 | 44217/2170 | -31383/2170 | 72/7 |

Check against (27) for [3,2,0], ω: H₁=2ω₁+3ω₂+5ω₃+12/θ₁₂=5+7+3+420/31=15+420/31=885/31.

## Appendix B. The kernels

### B.1. The adjoint kernel

**Table 17.** The adjoint kernel Gᵢⱼ(ω₁,ω₁′,t), t=1/θ₁₂, in the basis N=(P³p₀₂, P³′p₀₂, P³p₀₃, P³′p₀₃, P³p₂₃) (weights 2,3,3,4,5); G=4K. Source: `RESULTS-twovac-limit.md` section L2(d), parsed from the agent's printout (D2); equal to the independent fit of D1 in 25/25 entries.

| (i,j) | Gᵢⱼ |
|---|---|
| 1,1 | -4 ω₁² ω₁′² (3 t² + t ω₁ - t ω₁′ - ω₁ ω₁′) |
| 1,2 | 2 (6 t³ ω₁² + 24 t³ ω₁ ω₁′ + 6 t³ ω₁′² + 2 t² ω₁³ + 12 t² ω₁² ω₁′ - 2 t² ω₁′³ + 2 t ω₁³ ω₁′ - 6 t ω₁² ω₁′² - 2 t ω₁ ω₁′³ - 3 ω₁³ ω₁′²) |
| 1,3 | -2 ω₁² ω₁′² (2 t + ω₁) |
| 1,4 | 4 (-3 t⁴ - 4 t³ ω₁ - 2 t³ ω₁′ - t² ω₁² + 3 t² ω₁ ω₁′ + t² ω₁′² + 3 t ω₁² ω₁′ + t ω₁ ω₁′² + ω₁³ ω₁′) |
| 1,5 | -2 (6 t³ + 8 t² ω₁ + 4 t² ω₁′ + 4 t ω₁² + 2 t ω₁ ω₁′ + ω₁³) |
| 2,1 | -2 (6 t³ ω₁² + 24 t³ ω₁ ω₁′ + 6 t³ ω₁′² + 2 t² ω₁³ - 12 t² ω₁ ω₁′² - 2 t² ω₁′³ - 2 t ω₁³ ω₁′ - 6 t ω₁² ω₁′² + 2 t ω₁ ω₁′³ + 3 ω₁² ω₁′³) |
| 2,2 | -12 t² ω₁² - 48 t² ω₁ ω₁′ - 12 t² ω₁′² - 4 t ω₁³ - 12 t ω₁² ω₁′ + 12 t ω₁ ω₁′² + 4 t ω₁′³ + 9 ω₁² ω₁′² |
| 2,3 | 12 t⁴ + 16 t³ ω₁ + 8 t³ ω₁′ + 4 t² ω₁² - 12 t² ω₁ ω₁′ - 4 t² ω₁′² - 4 t ω₁² ω₁′ + 4 t ω₁ ω₁′² + 3 ω₁² ω₁′² |
| 2,4 | 2 (6 t³ + 8 t² ω₁ + 4 t² ω₁′ + 2 t ω₁² - 6 t ω₁ ω₁′ - 2 t ω₁′² - 3 ω₁² ω₁′) |
| 2,5 | 6 t² + 8 t ω₁ + 4 t ω₁′ + 3 ω₁² |
| 3,1 | 2 ω₁² ω₁′² (2 t - ω₁′) |
| 3,2 | 12 t⁴ - 8 t³ ω₁ - 16 t³ ω₁′ - 4 t² ω₁² - 12 t² ω₁ ω₁′ + 4 t² ω₁′² - 4 t ω₁² ω₁′ + 4 t ω₁ ω₁′² + 3 ω₁² ω₁′² |
| 3,3 | ω₁² ω₁′² |
| 3,4 | -2 (-6 t³ - 2 t² ω₁ + 2 t² ω₁′ + 2 t ω₁ ω₁′ + ω₁² ω₁′) |
| 3,5 | 6 t² + 4 t ω₁ + ω₁² |
| 4,1 | 4 (-3 t⁴ + 2 t³ ω₁ + 4 t³ ω₁′ + t² ω₁² + 3 t² ω₁ ω₁′ - t² ω₁′² - t ω₁² ω₁′ - 3 t ω₁ ω₁′² + ω₁ ω₁′³) |
| 4,2 | 2 (-6 t³ + 4 t² ω₁ + 8 t² ω₁′ + 2 t ω₁² + 6 t ω₁ ω₁′ - 2 t ω₁′² - 3 ω₁ ω₁′²) |
| 4,3 | -2 (6 t³ + 2 t² ω₁ - 2 t² ω₁′ - 2 t ω₁ ω₁′ + ω₁ ω₁′²) |
| 4,4 | -4 (3 t² + t ω₁ - t ω₁′ - ω₁ ω₁′) |
| 4,5 | -2 (2 t + ω₁) |
| 5,1 | -2 (-6 t³ + 4 t² ω₁ + 8 t² ω₁′ - 2 t ω₁ ω₁′ - 4 t ω₁′² + ω₁′³) |
| 5,2 | 6 t² - 4 t ω₁ - 8 t ω₁′ + 3 ω₁′² |
| 5,3 | 6 t² - 4 t ω₁′ + ω₁′² |
| 5,4 | 2 (2 t - ω₁′) |
| 5,5 | 1 |

**Table 18.** The adjoint kernel in twisted jets Ñ built from e^(-ωₛu)Pʲ (ωₛ=ω₁ on site 1, ω₁′ on site 2): it depends only on δ=ω₁′-ω₁ and t=1/θ₁₂. Same basis order as Table 17. Source: `RESULTS-twovac-limit.md`.

| (i,j) | G̃ᵢⱼ |
|---|---|
| 1,1 | -4 δ² t² (δ+ 3 t)(2 δ+ 3 t) |
| 1,2 | -4 δt² (2 δ+ 3 t)² |
| 1,3 | -4 δt² (δ² + 5 δt + 3 t²) |
| 1,4 | -4 t² (3 δ² + 8 δt + 3 t²) |
| 1,5 | -4 t² (2 δ+ 3 t) |
| 2,1 | 4 δt² (2 δ+ 3 t)² |
| 2,2 | 2 δt (δ+ 6 t)(2 δ+ 3 t) |
| 2,3 | 2 t² (5 δ² + 10 δt + 6 t²) |
| 2,4 | 4 t (δ² + 5 δt + 3 t²) |
| 2,5 | 2 t (2 δ+ 3 t) |
| 3,1 | 4 δt² (δ² + 5 δt + 3 t²) |
| 3,2 | 2 t² (5 δ² + 10 δt + 6 t²) |
| 3,3 | 2 δt² (δ+ 6 t) |
| 3,4 | 4 t² (2 δ+ 3 t) |
| 3,5 | 6 t² |
| 4,1 | -4 t² (3 δ² + 8 δt + 3 t²) |
| 4,2 | -4 t (δ² + 5 δt + 3 t²) |
| 4,3 | -4 t² (2 δ+ 3 t) |
| 4,4 | -4 t (δ+ 3 t) |
| 4,5 | -4 t |
| 5,1 | 4 t² (2 δ+ 3 t) |
| 5,2 | 2 t (2 δ+ 3 t) |
| 5,3 | 6 t² |
| 5,4 | 4 t |
| 5,5 | 1 |

For the adjoint the canonical Θ is G/(4ω₁³ω₁′³) with ω=1/a. A few entries in this normalisation (G1's K0 table, pairs (02),(03),(12),(13),(23), a′=a primed, t=1/θ₁₂): Θ₍₀₂₎₍₀₂₎=1+(a-a′)t-3aa′t², Θ₍₀₂₎₍₀₃₎=-a′(2at+1)/2, Θ₍₀₃₎₍₀₂₎=a(2a′t-1)/2, Θ₍₀₃₎₍₀₃₎=aa′/4, Θ₍₂₃₎₍₂₃₎=a³a′³/4, Θ₍₁₃₎₍₁₃₎=-a²a′²(3aa′t²-at+a′t-1), Θ₍₀₃₎₍₂₃₎=aa′³(6a²t²+4at+1)/4. The full table follows from G above or from the exchange rule Θᵢⱼ(a,a′,t)=Θⱼᵢ(a′,a,-t); all 25 entries are stored exactly in `general_kernel_rep_210.pkl`.

### B.2. [3,1,0] and [3,2,0]: base-point kernels

The exact kernel K at the base point (ω₁,ω₁′,θ₁₂)=(5/2,7/4,31/35) was evaluated from G1's stored exact polynomial laws (`general_kernel_rep_310.pkl`) in the basis of pairs (02),(03),(04),(12),(13),(14),(23),(24),(34), rows i and columns j; K₉₉=1/9. By the block theorem (55) the kernel of [3,2,0] is the sub-block with rows and columns (03),(04),(13),(14),(23),(24),(34) (we verified exact equality of all 49 entries, at this point and at (-2/5,3,29/45)); it reproduces the base-point values printed by P4 (for example K₁₁=-2477531875/3936256, K₁₂=-15649375/95232, K₁₇=1241317375/2770563, K₇₁=-3124391459/709264128, K₇₇=1/9 in the [3,2,0] order).

- Row (02): 315367650959375/30261936128, 718176615625/244047872, 6342284375/5904384, -6585110076593125/117265002496, 4028315502843156875/655020316754688, -92562651185431875/14086458424832, 37273401868703875/1760807303104, -84056493864875/113600471168, -6592621876625/2748398496.
- Row (03): -1163517096875/488095744, -2477531875/3936256, -15649375/95232, 56471003621623654875/3493441689358336, -102534604765838875/42259375274496, 201835815824125/113600471168, -1066229095827875/227200942336, 841817515875/1832265664, 1241317375/2770563.
- Row (04): 4287885875/23617536, 8703625/190464, 42875/4608, -59815795579096375/42259375274496, 213468283144375/1022404240512, -1720011432125/10993593984, 29908228000/85887453, -4528239625/88658016, -66781625/2144952.
- Row (12): 2958159365546875/117265002496, -26499034046453294375/1746720844679168, 38085651450599375/10564843818624, 142414335777875/3782742016, -409392543660729125/42259375274496, 917570019491125/227200942336, -1924121671270925/227200942336, 399981133475/458066416, 30381041525/44329008.
- Row (13): -45956030935735311125/5240162534037504, 92531938123379125/21129687637248, -1336404322859375/1022404240512, -1017838029413828875/84518750548992, 2973626853294125/1022404240512, -7269823449875/5496796992, 12484966188275/5496796992, -1902394025/5541126, -179699825/1072476.
- Row (14): 94733549284435075/112691667398656, -319234542990475/908803769344, 2843820794675/21987187968, 496844069229025/454401884672, -651463216075/2748398496, 7355513725/59105344, -5263492675/29552672, 16750475/476656, 147925/11532.
- Row (23): 4226124769018075/1760807303104, -347834314964675/454401884672, 1716815095975/5496796992, 1917549022254845/908803769344, -5119851097655/10993593984, 13720239855/59105344, -17517770585/59105344, 14281435/238328, 227285/11532.
- Row (24): -584809663796365/1817607538688, 1123221350545/14658125312, -14416082485/354632064, -462993234235/1832265664, 4168591945/88658016, -27137425/953312, 6967205/238328, -57785/7688, -365/186.
- Row (34): 2805654394403/87948751872, -3124391459/709264128, 60638627/17159616, 6925335571/354632064, -2906729/1072476, 100009/46128, -37121/23064, 203/372, ⅑.

Laws (all verified, Table 9 and Section 8.4): each entry is a polynomial in (ω₁,ω₁′). For [3,2,0] the table of total degrees is rows 8,7,7,6,6,5,4 / 7,6,6,5,5,4,3 (twice) / 6,5,5,4,4,3,2 (twice) / 5,4,4,3,3,2,1 / 4,3,3,2,2,1,0, i.e.\ dᵢⱼ=8-(eᵢ+eⱼ) with e=(0,1,1,2,2,3,4), at most degree 4 in each of ω₁,ω₁′; for [3,1,0], e=(0,1,2,1,2,3,3,4,5), bidegree (5,10) in (1/ω₁′,t) at fixed a, 998 monomials; the dependence on θ is through θ₁₂ only (equal-θ₁₂ pairs agree to 4.8·10⁻⁷², 8·10⁻⁷⁰, 2.4·10⁻⁶⁸ for [2,1,0],[2,0,0],[3,2,0]).

## Appendix C. The one-site kernel

Table 19 lists the coefficients of θₖₗ(a)=cₖₗa^(k+l-λ₂-1) for all two-row and one-row representations up to λ₁=4. They are computed from the closed form (51); we cross-checked the entries printed by G1 for [1,0,0],[1,1,0],[2,0,0],[2,1,0],[2,2,0],[3,1,0],[3,2,0],[4,2,0] (they agree).

**Table 19.** The one-site kernel θₖₗ(a)=cₖₗ a^(k+l-λ₂-1): coefficients cₖₗ from the closed form (eq. 51) for all two-row (and one-row) representations up to λ₁=4. Pairs (kl) in lexicographic order, 0≤k<l≤λ₁+1, l≥λ₂+1. The entries for Λ≤[4,2,0] coincide with the exact one-site fits of G1; the rest were predicted by the formula and then confirmed (17 representations in all).

| Λ | # pairs | cₖₗ |
|---|---|---|
| [1,0,0] | 3 | (01): 1, (02): -1, (12): 1 |
| [1,1,0] | 2 | (02): 1, (12): -1 |
| [2,0,0] | 6 | (01): 1, (02): -2, (03): 1, (12): 3, (13): -2, (23): 1 |
| [2,1,0] | 5 | (02): 1, (03): -½, (12): -3/2, (13): 1, (23): -½ |
| [2,2,0] | 3 | (03): 1, (13): -2, (23): 1 |
| [3,0,0] | 10 | (01): 1, (02): -3, (03): 3, (04): -1, (12): 6, (13): -8, (14): 3, (23): 6, (24): -3, (34): 1 |
| [3,1,0] | 9 | (02): 1, (03): -1, (04): ⅓, (12): -2, (13): 8/3, (14): -1, (23): -2, (24): 1, (34): -⅓ |
| [3,2,0] | 7 | (03): 1, (04): -⅓, (13): -8/3, (14): 1, (23): 2, (24): -1, (34): ⅓ |
| [3,3,0] | 4 | (04): 1, (14): -3, (24): 3, (34): -1 |
| [4,0,0] | 15 | (01): 1, (02): -4, (03): 6, (04): -4, (05): 1, (12): 10, (13): -20, (14): 15, (15): -4, (23): 20, (24): -20, (25): 6, (34): 10, (35): -4, (45): 1 |
| [4,1,0] | 14 | (02): 1, (03): -3/2, (04): 1, (05): -¼, (12): -5/2, (13): 5, (14): -15/4, (15): 1, (23): -5, (24): 5, (25): -3/2, (34): -5/2, (35): 1, (45): -¼ |
| [4,2,0] | 12 | (03): 1, (04): -⅔, (05): ⅙, (13): -10/3, (14): 5/2, (15): -⅔, (23): 10/3, (24): -10/3, (25): 1, (34): 5/3, (35): -⅔, (45): ⅙ |
| [4,3,0] | 9 | (04): 1, (05): -¼, (14): -15/4, (15): 1, (24): 5, (25): -3/2, (34): -5/2, (35): 1, (45): -¼ |
| [4,4,0] | 5 | (05): 1, (15): -4, (25): 6, (35): -4, (45): 1 |

## Appendix D. The ε-series and the CG tables

**Table 20.** Exact ε-series coefficients of the Laurent data of the six adjoint branches at φ=(1/3,-2/7,5/11), θ₁₂=31/35 (so H₁⁽⁰⁾=4/θ₁₂ etc.), orders n=0,1,2 (Hₛ=ΣₙεⁿHₛ⁽ⁿ⁾). Source: `ref/untwisted-exact-210.txt`, record `E4_SERIES` (agent W4; exact Rayleigh–Schrödinger).

| branch | series | n=0 | n=1 | n=2 |
|---|---|---|---|---|
| (4,2,0) | H₁ | 140/31 | 260/231 | 0 |
| (4,2,0) | a₁ | 0 | 15/11 | 0 |
| (4,2,0) | b₁ | 0 | 700/341 | -25/231 |
| (4,1,1) | H₁ | 210/31 | 349/462 | 100719/1660120 |
| (4,1,1) | a₁ | 210/31 | 39/154 | 302157/1660120 |
| (4,1,1) | b₁ | 0 | 1160/1023 | -107/462 |
| (3,3,0) | H₁ | 210/31 | 221/154 | 5239/123480 |
| (3,3,0) | a₁ | 0 | 15/11 | 0 |
| (3,3,0) | b₁ | 0 | 1050/341 | 5/154 |
| (3,2,1) copy 1 | H₁ | 315/31 | 82/77+157/1155 √5 | 4645474/140072625 |
| (3,2,1) copy 1 | a₁ | 525/62-105/62 √5 | 509/770+541/2310 √5 | 7649591/140072625-13875662/700363125 √5 |
| (3,2,1) copy 1 | b₁ | 3675/961 √5 | 76/31 | -1/11+1952029/20010375 √5 |
| (3,2,1) copy 2 | H₁ | 315/31 | 82/77-157/1155 √5 | 4645474/140072625 |
| (3,2,1) copy 2 | a₁ | 525/62+105/62 √5 | 509/770-541/2310 √5 | 7649591/140072625+13875662/700363125 √5 |
| (3,2,1) copy 2 | b₁ | -3675/961 √5 | 76/31 | -1/11-1952029/20010375 √5 |
| (2,2,2) | H₁ | 420/31 | 232/231 | 782347/11205810 |
| (2,2,2) | a₁ | 420/31 | 116/231 | 782347/7470540 |
| (2,2,2) | b₁ | 0 | 2320/1023 | -17/231 |

**Table 21.** The same series, numerically (6 significant digits), orders n=0,…,7; exact values to order 7 are in the data file. Odd orders vanish identically on the algebraic branches (4,1,1) and (3,3,0) beyond n=1 and on the singlet for H₁ (see text).

| branch | series | n=0 | n=1 | n=2 | n=3 | n=4 | n=5 | n=6 | n=7 |
|---|---|---|---|---|---|---|---|---|---|
| (4,2,0) | H₁ | 4.51613 | 1.12554 | 0 | 0 | 0 | 0 | 0 | 0 |
| (4,2,0) | a₁ | 0 | 1.36364 | 0 | 0 | 0 | 0 | 0 | 0 |
| (4,2,0) | b₁ | 0 | 2.05279 | -0.108225 | 0 | 0 | 0 | 0 | 0 |
| (4,1,1) | H₁ | 6.77419 | 0.755411 | 0.0606697 | 0 | -0.00163007 | 0 | 8.75937e-05 | 0 |
| (4,1,1) | a₁ | 6.77419 | 0.253247 | 0.182009 | 0 | -0.00489022 | 0 | 0.000262781 | 0 |
| (4,1,1) | b₁ | 0 | 1.13392 | -0.231602 | 0.0202232 | 0 | -0.000543358 | 0 | 2.91979e-05 |
| (3,3,0) | H₁ | 6.77419 | 1.43506 | 0.0424279 | 0 | -0.0007972 | 0 | 2.9958e-05 | 0 |
| (3,3,0) | a₁ | 0 | 1.36364 | 0 | 0 | 0 | 0 | 0 | 0 |
| (3,3,0) | b₁ | 0 | 3.07918 | 0.0324675 | 0.0192854 | 0 | -0.000362364 | 0 | 1.36173e-05 |
| (3,2,1) copy 1 | H₁ | 10.1613 | 1.36889 | 0.0331648 | -9.34224e-05 | -0.00028439 | 7.10857e-06 | 4.61811e-06 | -2.62555e-07 |
| (3,2,1) copy 1 | a₁ | 4.68085 | 1.18472 | 0.0103104 | 0.000820064 | -0.000148242 | -8.918e-06 | 3.12722e-06 | 7.11436e-08 |
| (3,2,1) copy 1 | b₁ | 8.55104 | 2.45161 | 0.127221 | 0.00973109 | -0.000715096 | -5.10438e-05 | 1.19737e-05 | 3.85103e-07 |
| (3,2,1) copy 2 | H₁ | 10.1613 | 0.760985 | 0.0331648 | 9.34224e-05 | -0.00028439 | -7.10857e-06 | 4.61811e-06 | 2.62555e-07 |
| (3,2,1) copy 2 | a₁ | 12.2546 | 0.137354 | 0.0989128 | 0.00171494 | -0.000786658 | -4.38194e-05 | 1.17683e-05 | 1.1896e-06 |
| (3,2,1) copy 2 | b₁ | -8.55104 | 2.45161 | -0.309039 | 0.00973109 | 0.000715096 | -5.10438e-05 | -1.19737e-05 | 3.85103e-07 |
| (2,2,2) | H₁ | 13.5484 | 1.00433 | 0.0698162 | 0 | -0.000539655 | 0 | 2.38687e-06 | 0 |
| (2,2,2) | a₁ | 13.5484 | 0.502165 | 0.104724 | -0.00635184 | -0.000809482 | 0.000180024 | 3.5803e-06 | -4.80708e-06 |
| (2,2,2) | b₁ | 0 | 2.26784 | -0.0735931 | 0.0403721 | 0 | -0.000607703 | 0 | 1.14447e-05 |

Symbolic data to higher order are in `untwisted_Q_pert_sym.pkl` (singlet to n=6, algebraic branches to n=6), `untwisted_Q_laurent_exact.pkl` (all branches n≤5), `untwisted_Q_octet_fit.pkl` and the data-file records `E4_SYMBOLIC` (orders 0..3), `E4_PSI` (to order 4), `E4_GNORM`, `E5_KB`, `E5_SLICE`, `E5_GF`. The numerical ε-series of the monomials Nᵢ,M, the Plücker coordinates and the polynomial coefficients are in `untwisted_Q_series_210.pkl` (12 points, 200 digits, n≤10), `_ext.pkl` (24 points, 400 digits, n≤20, reliable to 18) and `_phi2.pkl` (φ′).

### D.1. CG tables

**Table 22.** Exact CG slices C_ν=⟨e_ν⊗e₁|Ψ₀⟩/‖Ψ₀‖_G of the highest-weight kets of V_Λ⊗V_Λ, six representations (45 HW branches); ν is the site-1 Gelfand–Tsetlin index in notebook order (index 1 = highest weight) with its weight; the two copies in a multiplicity-2 sector are the eigenvectors of θ₁₂ a₁, copy 1 having the smaller eigenvalue; C_ν=(c_A+c_B√D₀)/√(g_A+g_B√D₀). Source: `ref/untwisted-exact-<Lambda>.txt`, records `BRANCH`.

| Λ | μ | copy | ν, weight | c_A+c_B√D₀ | g_A+g_B√D₀ | C_ν |
|---|---|---|---|---|---|---|
| (1,1,0) | (2,2,0) |  | 1, (1,1,0) | 1 | 1 | 1.0000000000 |
| (1,1,0) | (2,1,1) |  | 2, (1,0,1) | -1 | 2 | -0.7071067812 |
| (2,1,0) | (4,2,0) |  | 1, (2,1,0) | 1 | 1 | 1.0000000000 |
| (2,1,0) | (4,1,1) |  | 2, (2,0,1) | -1 | 2 | -0.7071067812 |
| (2,1,0) | (3,3,0) |  | 3, (1,2,0) | -1 | 2 | -0.7071067812 |
| (2,1,0) | (3,2,1) | 1 | 4, (1,1,1) | 21-9√5 | 7560-3240√5 | 0.0493116000 |
|  |  |  | 5, (1,1,1) | -3+3√5 |  | 0.2088872897 |
| (2,1,0) | (3,2,1) | 2 | 4, (1,1,1) | 21+9√5 | 7560+3240√5 | 0.3379867346 |
|  |  |  | 5, (1,1,1) | -3-3√5 |  | -0.0797878449 |
| (2,1,0) | (2,2,2) |  | 8, (0,1,2) | 2 | 288 | 0.1178511302 |
| (2,2,0) | (4,4,0) |  | 1, (2,2,0) | 1 | 1 | 1.0000000000 |
| (2,2,0) | (4,3,1) |  | 2, (2,1,1) | -1 | 4 | -0.5000000000 |
| (2,2,0) | (4,2,2) |  | 3, (2,0,2) | 1 | 12 | 0.2886751346 |
| (3,1,0) | (6,2,0) |  | 1, (3,1,0) | 1 | 1 | 1.0000000000 |
| (3,1,0) | (6,1,1) |  | 2, (3,0,1) | -1 | 2 | -0.7071067812 |
| (3,1,0) | (5,3,0) |  | 3, (2,2,0) | -1 | 4 | -0.5000000000 |
| (3,1,0) | (5,2,1) | 1 | 4, (2,1,1) | 264-24√105 | 1496880-136080√105 | 0.0564580374 |
|  |  |  | 5, (2,1,1) | -30+6√105 |  | 0.0983443476 |
| (3,1,0) | (5,2,1) | 2 | 4, (2,1,1) | 264+24√105 | 1496880+136080√105 | 0.2998902852 |
|  |  |  | 5, (2,1,1) | -30-6√105 |  | -0.0538008073 |
| (3,1,0) | (4,4,0) |  | 7, (1,3,0) | 1 | 12 | 0.2886751346 |
| (3,1,0) | (4,3,1) | 1 | 8, (1,2,1) | -168+96√3 | 5225472-2985984√3 | -0.0074430331 |
|  |  |  | 9, (1,2,1) | 24-24√3 |  | -0.0758903002 |
| (3,1,0) | (4,3,1) | 2 | 8, (1,2,1) | -168-96√3 | 5225472+2985984√3 | -0.1036680780 |
|  |  |  | 9, (1,2,1) | 24+24√3 |  | 0.0203347447 |
| (3,1,0) | (4,2,2) |  | 10, (1,1,2) | 18 | 233280 | 0.0372677996 |
|  |  |  | 11, (1,1,2) | -6 |  | -0.0124225999 |
| (3,1,0) | (3,3,2) |  | 14, (0,2,2) | -3 | 17280 | -0.0228217732 |
| (3,2,0) | (6,4,0) |  | 1, (3,2,0) | 1 | 1 | 1.0000000000 |
| (3,2,0) | (6,3,1) |  | 2, (3,1,1) | -1 | 4 | -0.5000000000 |
| (3,2,0) | (6,2,2) |  | 3, (3,0,2) | 1 | 12 | 0.2886751346 |
| (3,2,0) | (5,5,0) |  | 4, (2,3,0) | -1 | 2 | -0.7071067812 |
| (3,2,0) | (5,4,1) | 1 | 5, (2,2,1) | 22-2√105 | 18480-1680√105 | 0.0423435280 |
|  |  |  | 7, (2,2,1) | -3+√105 |  | 0.2037459499 |
| (3,2,0) | (5,4,1) | 2 | 5, (2,2,1) | 22+2√105 | 18480+1680√105 | 0.2249177139 |
|  |  |  | 7, (2,2,1) | -3-√105 |  | -0.0701153289 |
| (3,2,0) | (5,3,2) | 1 | 6, (2,1,2) | -168+96√3 | 5225472-2985984√3 | -0.0074430331 |
|  |  |  | 8, (2,1,2) | 24-24√3 |  | -0.0758903002 |
| (3,2,0) | (5,3,2) | 2 | 6, (2,1,2) | -168-96√3 | 5225472+2985984√3 | -0.1036680780 |
|  |  |  | 8, (2,1,2) | 24+24√3 |  | 0.0203347447 |
| (3,2,0) | (4,4,2) |  | 12, (1,2,2) | 4 | 2880 | 0.0745355992 |
| (3,2,0) | (4,3,3) |  | 13, (1,1,3) | -3 | 17280 | -0.0228217732 |
| (4,1,0) | (8,2,0) |  | 1, (4,1,0) | 1 | 1 | 1.0000000000 |
| (4,1,0) | (8,1,1) |  | 2, (4,0,1) | -1 | 2 | -0.7071067812 |
| (4,1,0) | (7,3,0) |  | 3, (3,2,0) | -1 | 6 | -0.4082482905 |
| (4,1,0) | (7,2,1) | 1 | 4, (3,1,1) | 450-90√21 | 5443200-1088640√21 | 0.0557300243 |
|  |  |  | 5, (3,1,1) | -42+18√21 |  | 0.0600589519 |
| (4,1,0) | (7,2,1) | 2 | 4, (3,1,1) | 450+90√21 | 5443200+1088640√21 | 0.2670185879 |
|  |  |  | 5, (3,1,1) | -42-18√21 |  | -0.0385423778 |
| (4,1,0) | (6,4,0) |  | 7, (2,3,0) | 3 | 360 | 0.1581138830 |
| (4,1,0) | (6,3,1) | 1 | 8, (2,2,1) | -150+60√6 | 6912000-2764800√6 | -0.0081097820 |
|  |  |  | 9, (2,2,1) | 18-12√6 |  | -0.0304894699 |
| (4,1,0) | (6,3,1) | 2 | 8, (2,2,1) | -150-60√6 | 6912000+2764800√6 | -0.0802785656 |
|  |  |  | 9, (2,2,1) | 18+12√6 |  | 0.0128118003 |
| (4,1,0) | (6,2,2) |  | 10, (2,1,2) | 40 | 4838400 | 0.0181848242 |
|  |  |  | 11, (2,1,2) | -10 |  | -0.0045462060 |
| (4,1,0) | (5,5,0) |  | 13, (1,4,0) | -1 | 144 | -0.0833333333 |
| (4,1,0) | (5,4,1) | 1 | 14, (1,3,1) | 690-150√21 | 1669248000-362880000√21 | 0.0010394114 |
|  |  |  | 15, (1,3,1) | -90+30√21 |  | 0.0188810654 |
| (4,1,0) | (5,4,1) | 2 | 14, (1,3,1) | 690+150√21 | 1669248000+362880000√21 | 0.0238611846 |
|  |  |  | 15, (1,3,1) | -90-30√21 |  | -0.0039407078 |
| (4,1,0) | (5,3,2) |  | 16, (1,2,2) | -120 | 261273600 | -0.0074239234 |
|  |  |  | 17, (1,2,2) | 30 |  | 0.0018559808 |
| (4,1,0) | (4,4,2) |  | 22, (0,3,2) | 4 | 1036800 | 0.0039283710 |

## Appendix E. The untwisted data

**Table 23.** Untwisted Laurent data and oper data on all 45 highest-weight branches: H₁=-H₂=η/θ₁₂, a₁=-a₂=α/θ₁₂, b₁=-b₂=β/θ₁₂²; exponents of the dual oper at infinity (μ₁+2,μ₂+1,μ₃) (the three polynomial solutions have exactly these degrees); degrees of the three polynomial solutions of the untwisted direct oper (see text for the two exceptions); last column: whether the oper is a ₃F₂ equation (consecutive exponent pair at a singular point). Source: `ref/untwisted-exact-<Lambda>.txt`.

| Λ | μ | copy | η | α | β | at ∞ | direct deg. | ₃F₂ |
|---|---|---|---|---|---|---|---|---|
| (1,1,0) | (2,2,0) |  | 2 | 0 | 0 | (4,3,0) | (0,2,4) | yes |
| (1,1,0) | (2,1,1) |  | 4 | 4 | 0 | (4,2,1) | (3,2,0) | yes |
| (2,1,0) | (4,2,0) |  | 4 | 0 | 0 | (6,3,0) | (6,3,0) | no |
| (2,1,0) | (4,1,1) |  | 6 | 6 | 0 | (6,2,1) | (5,4,0) | yes |
| (2,1,0) | (3,3,0) |  | 6 | 0 | 0 | (5,4,0) | (6,2,1) | yes |
| (2,1,0) | (3,2,1) | 1 | 9 | 15/2-3/2 √5 | 3√5 | (5,3,1) | (5,3,1) | no |
| (2,1,0) | (3,2,1) | 2 | 9 | 15/2+3/2 √5 | -3√5 | (5,3,1) | (5,3,1) | no |
| (2,1,0) | (2,2,2) |  | 12 | 12 | 0 | (4,3,2) | (4,3,2) | yes |
| (2,2,0) | (4,4,0) |  | 8 | 0 | 0 | (6,5,0) | (0,2,6) | yes |
| (2,2,0) | (4,3,1) |  | 12 | 12 | 0 | (6,4,1) | (5,2,0) | yes |
| (2,2,0) | (4,2,2) |  | 14 | 18 | 0 | (6,3,2) | (4,3,0) | yes |
| (3,1,0) | (6,2,0) |  | 6 | 0 | 0 | (8,3,0) | (8,5,0) | no |
| (3,1,0) | (6,1,1) |  | 8 | 8 | 0 | (8,2,1) | (7,6,0) | yes |
| (3,1,0) | (5,3,0) |  | 10 | 0 | 0 | (7,4,0) | (8,4,1) | no |
| (3,1,0) | (5,2,1) | 1 | 13 | 21/2-½ √105 | √105 | (7,3,1) | (7,5,1) | no |
| (3,1,0) | (5,2,1) | 2 | 13 | 21/2+½ √105 | -√105 | (7,3,1) | (7,5,1) | no |
| (3,1,0) | (4,4,0) |  | 12 | 0 | 0 | (6,5,0) | (8,3,2) | yes |
| (3,1,0) | (4,3,1) | 1 | 16 | 12-4√3 | 8√3 | (6,4,1) | (7,4,2) | no |
| (3,1,0) | (4,3,1) | 2 | 16 | 12+4√3 | -8√3 | (6,4,1) | (7,4,2) | no |
| (3,1,0) | (4,2,2) |  | 18 | 18 | 0 | (6,3,2) | (6,5,2) | yes |
| (3,1,0) | (3,3,2) |  | 20 | 20 | 0 | (5,4,2) | (6,4,3) | yes |
| (3,2,0) | (6,4,0) |  | 12 | 0 | 0 | (8,5,0) | (8,3,0) | no |
| (3,2,0) | (6,3,1) |  | 16 | 16 | 0 | (8,4,1) | (7,4,0) | no |
| (3,2,0) | (6,2,2) |  | 18 | 24 | 0 | (8,3,2) | (6,5,0) | yes |
| (3,2,0) | (5,5,0) |  | 14 | 0 | 0 | (7,6,0) | (8,2,1) | yes |
| (3,2,0) | (5,4,1) | 1 | 19 | 35/2-½ √105 | √105 | (7,5,1) | (7,3,1) | no |
| (3,2,0) | (5,4,1) | 2 | 19 | 35/2+½ √105 | -√105 | (7,5,1) | (7,3,1) | no |
| (3,2,0) | (5,3,2) | 1 | 22 | 28-4√3 | 8√3 | (7,4,2) | (6,4,1) | no |
| (3,2,0) | (5,3,2) | 2 | 22 | 28+4√3 | -8√3 | (7,4,2) | (6,4,1) | no |
| (3,2,0) | (4,4,2) |  | 24 | 30 | 0 | (6,5,2) | (6,3,2) | yes |
| (3,2,0) | (4,3,3) |  | 26 | 36 | 0 | (6,4,3) | (5,4,2) | yes |
| (4,1,0) | (8,2,0) |  | 8 | 0 | 0 | (10,3,0) | (10,7,0) | no |
| (4,1,0) | (8,1,1) |  | 10 | 10 | 0 | (10,2,1) | (9,8,0) | yes |
| (4,1,0) | (7,3,0) |  | 14 | 0 | 0 | (9,4,0) | (10,6,1) | no |
| (4,1,0) | (7,2,1) | 1 | 17 | 27/2-3/2 √21 | 3√21 | (9,3,1) | (9,7,1) | no |
| (4,1,0) | (7,2,1) | 2 | 17 | 27/2+3/2 √21 | -3√21 | (9,3,1) | (9,7,1) | no |
| (4,1,0) | (6,4,0) |  | 18 | 0 | 0 | (8,5,0) | (10,5,2) | no |
| (4,1,0) | (6,3,1) | 1 | 22 | 16-4√6 | 8√6 | (8,4,1) | (9,6,2) | no |
| (4,1,0) | (6,3,1) | 2 | 22 | 16+4√6 | -8√6 | (8,4,1) | (9,6,2) | no |
| (4,1,0) | (6,2,2) |  | 24 | 24 | 0 | (8,3,2) | (8,7,2) | yes |
| (4,1,0) | (5,5,0) |  | 20 | 0 | 0 | (7,6,0) | (10,4,3) | yes |
| (4,1,0) | (5,4,1) | 1 | 25 | 35/2-5/2 √21 | 5√21 | (7,5,1) | (9,5,3) | no |
| (4,1,0) | (5,4,1) | 2 | 25 | 35/2+5/2 √21 | -5√21 | (7,5,1) | (9,5,3) | no |
| (4,1,0) | (5,3,2) |  | 28 | 28 | 0 | (7,4,2) | (8,6,3) | no |
| (4,1,0) | (4,4,2) |  | 30 | 30 | 0 | (6,5,2) | (8,5,4) | yes |

The sectors of each V_Λ⊗V_Λ and their multiplicities: [1,1,0]: (2,2,0),(2,1,1). [2,1,0]: (4,2,0),(4,1,1),(3,3,0),(3,2,1)×2,(2,2,2). [2,2,0]: (4,4,0),(4,3,1),(4,2,2). [3,1,0]: (6,2,0),(6,1,1),(5,3,0),(5,2,1)×2,(4,4,0),(4,3,1)×2,(4,2,2),(3,3,2). [3,2,0]: (6,4,0),(6,3,1),(6,2,2),(5,5,0),(5,4,1)×2,(5,3,2)×2,(4,4,2),(4,3,3). [4,1,0]: (8,2,0),(8,1,1),(7,3,0),(7,2,1)×2,(6,4,0),(6,3,1)×2,(6,2,2),(5,5,0),(5,4,1)×2,(5,3,2),(4,4,2). Multiplicity ≤2 everywhere (45 branches in all). One-site Casimirs (C₂,C₃)(Λ) from the data files: [1,1,0]: (4,8); [2,1,0]: (9,30); [2,2,0]: (12,36); [3,1,0]: (16,72); [3,2,0]: (19,77); [4,1,0]: (25,140). The Laurent operators in the canonical HW basis, the Gram matrices and the HW kets themselves are in the files (records `LAURMAT`, `GRAM`, `HWBASIS`).

## Appendix F. Reproducibility

All Python work lives in `Paul/Python/gaudin_sl3/` (run with `python -X utf8`, Python 3.13, sympy 1.14, mpmath 1.3, numpy; no `uv`), all notebooks in `Paul/Mathematica/Gaudin/Experiments/`, reference data in `Paul/Mathematica/Gaudin/ref/`. Nothing was committed by the agents; about 230 MB of `*.pkl` caches in `Paul/Python/gaudin_sl3/` should be gitignored before committing. The four notebooks store no native outputs (results are in markdown cells below each code cell); Paul reruns them in VS Code. The common pitfall: run heavy stages separately (an all-stage run stalled under CPU contention from parallel agents).

**Table 24.** Which file produces what. Times are those reported by the agents (warm = with caches, cold = from scratch). Specs are `SPEC-*.md`/`*-spec.md` next to each file.

| agent | files | produces | time, precision |
|---|---|---|---|
| W1 | `Experiments/Gaudin-L2-chain.wb` (48 cells, 16 code) | Hamiltonians, vacuum, localisers, two-vacua untwist and CG for [2,1,0], [3,1,0] (switch Λ in the first code cell) | one clean run, 200 digits |
| W2 | `Experiments/Gaudin-L2-chain-2.wb` (80 cells, 21 code) | (∞,∞) family, [3,1,0] fourth localiser, closed forms, A5 reference tables | every code cell run twice, 200 digits |
| W3 | `Experiments/Gaudin-L2-chain-320.wb` (75 cells, 32 code) | [3,2,0] chain side, cleaning verdict, A7 tables | one clean run, 8 min; scan cell 59 about 6 min cold |
| W4 | `Experiments/Gaudin-untwisted-exact.wb` (39 cells, 14 code); `ref/untwisted-exact-{110,210,220,310,320,410}.txt` | exact untwisted reference data, perturbation theory, CG tables | one clean run, exact |
| P1 | `gaudin_chain.py`, `qside.py`, `RESULTS-Qside.md`, `qside_output.txt`, `qside_cache.pkl` | Q side for [2,1,0] (opers, Wronskians, detN, regime) | about 75 s; 60/100 digits |
| P2 | `xxx_chain.py`, `xxx_limit.py`, `RESULTS-xxx-limit.md` | XXX→Gaudin limit, dictionary, x₀, localiser prediction | about 2 min (X4 sympy 95 s); exact series |
| P3 | `qside_310.py`, `RESULTS-Qside-310.md`, `qside_310_*.pkl` | [3,1,0] Q side; rank-4 rule | 4 min warm, 15 min cold; 60/100 digits |
| P4 | `qside_320.py`, `twovac_fit_320.py` (stages T0..T7), `RESULTS-Qside-320.md`, `twovac_fit_320_*.pkl` | [3,2,0] Q side and jet formula, [3,1,0] comparison | 4 min warm / 8 min cold; 60/80/100 digits |
| D1 | `twovac_fit.py`, `RESULTS-twovac-fit.md`, `twovac_fit_*.pkl` | direct fit of the adjoint kernel, kernel laws, regime Q-side CG | 2.5 min warm; 80 digits (F4: 120) |
| D2 | `twovac_limit.py`, `RESULTS-twovac-limit.md`; input `ref/xxx-210-infinf-formula-cells.txt` | h→0 limit of the ten-product formula, G | about 1 min; exact, 60 digits for L3 |
| D3 | `gaudin_B.py`, `RESULTS-gaudin-B.md`, `gaudin_B_*.pkl` | B hierarchy, jet basis, K tensor, G from first principles | 6 min warm / 15 min cold; exact Fractions, 80 digits for states; `GB_THETA=1/5,-4/9` repeats at another θ |
| G1 | `general_kernel.py` (stages K1, REP, ZC, BLK, SGL, CORNER, ANALYZE), `general_kernel_bld.py`, `RESULTS-general-kernel.md`, `general_kernel_rep_<Lambda>.pkl` | general-Λ kernels, one-site closed form, block theorem | 80 digits; per point ∼35 s ([3,1,0], [3,2,0]), [4,1,0] about 2 min, [4,2,0] about 4 min |
| U1 | `untwisted_Q.py` (stages T0, T1, EXPORT), `RESULTS-untwisted-Q.md`, `untwisted_Q_*.pkl` | untwisted Laurent data, oper, ε-series | 200/400 digits |
| U2 | `untwist_limit.py` (stages R, N, NU), `untwist_limit_frob.py`, `RESULTS-untwist-limit.md`, `untwist_limit_reg_<k>.pkl` | analytic untwisting limit for the adjoint | about 1 min with caches; 200 digits, 12 ε points ∼1 s each |
| data | `ref/ref-*.txt` | cell dumps of the XXX notebooks used by the Python agents | – |

Other pitfalls recorded by the agents: a numerically vanishing feature (p₀₁∼10⁻⁷⁸) must be zeroed explicitly before QR or it is counted as an independent feature; the 100-digit regime solver needs a kernel threshold of 10⁻⁵⁰ (the default 10⁻³⁵ gives spurious lower degrees); a threshold-based degree reading of q fails when the leading coefficient is ε-suppressed (read the smallest degree admitting a kernel); the memoised `laurent` cache holds bad values after a fit at a pole; the mp eigensolver returns 10⁻⁵⁹ imaginary noise in some sectors; the all-stage Python run stalled under contention. The notebooks were built from JSON and run cell by cell in the shared kernel (wolfbook new-notebook workaround), hence no native outputs.

The sources of this document are the plans `Gaudin-sl3-CG-plan.md`, `-round2-`, `-round3-`, `-round4-`, `-round5-plan.md`, the ten `RESULTS-*.md` files, the four notebooks and the six data files named above, and the repository's `HANDOFF.md` (entries from “sl(3) GAUDIN” onwards). The generator of this document (hand-written prose plus tables built from the data files) lives in `Paul/Mathematica/Gaudin/notes-src/`.

---

[^1]: Paul's steer of 2026-10-16: “the goal of the project is Gaudin anyway, not XXX; I assumed XXX would be simpler”; and later “project target is gl(3) CG coefficients from SoV in the Gaudin model”. `HANDOFF.md`, entries “ROUND 2” and “ROUND 4”. The repository's `CLAUDE.md` still states a su(4) goal; that sentence is outdated and was deliberately not edited by the notes writer.

[^2]: All representations are written in the Gelfand–Tsetlin (GT) basis with plain, non-Hermitian rational generators Eᵢⱼ and the invariant diagonal metric G (EᵢⱼᵀG=GEⱼᵢ); this is the convention of `Paul/Python/gl3_reps.py` and of the notebooks' `JJ`. Raising operators are Eᵢⱼ with i<j; e₁ is the highest-weight ket and covectors are rows acted on from the right, so e₁ᵀEᵢⱼ≠0 only for i≤j. Sources: `Gaudin-sl3-CG-plan.md` section 1; `RESULTS-Qside.md` section 0.

[^3]: Agents W1 (`Gaudin-L2-chain.wb`, summary cell), P1 (`RESULTS-Qside.md` section 0), P3, P4. The convention was not tuned: it worked at the first attempt. The row-ordered product with Eᵢⱼ and the column-ordered product with Eⱼᵢ do *not* commute (274 non-zero entries in [S₂,S₃]); the doubly flipped variant is just the transposed chain.

[^4]: Agent P3 (`RESULTS-Qside-310.md` section 0) found c₂=λ₂(λ₁+1) on (1,1,0),(2,1,0),(3,1,0),(2,2,0),(3,2,0) and reported the double pole as “not scalar” for (n,0,0); agent G1 (`RESULTS-general-kernel.md` K2) later checked [1,0,0],[2,0,0],[3,0,0] on 9/36/100 states and found the scalar 0 every time. We take G1's statement, as the later round. (2,1,1) has c₂=8 but an S₃ triple pole and was not used.

[^5]: Agent W1, `Gaudin-L2-chain.wb` summary cell (“Sign difference to the plan”), exact on [2,1,0]; agents P3 and P4 report the same explicit form for [3,1,0] and [3,2,0] (exact fits, 0 mismatches). The plan `Gaudin-sl3-CG-plan.md` section 1 carries the corrected row.

[^6]: Agent P2, `RESULTS-xxx-limit.md` section X2: exact rational series in h, checked at four values of u and three twist sets, and at finite h to O(h³) residuals. P2 also found the vacuum, frame and Baxter-to-oper statements used below.

[^7]: Exactly linear: agent D2, `RESULTS-twovac-limit.md` L1(b), with the XXX x₀(z) exactly linear in z at every order in h; x₀(e^(hζ))=h x(ζ)+O(h²). Degree λ₂ in general: agent W4, `Gaudin-untwisted-exact.wb` E5 and the `XVEC` records of the data files ([1,1,0]: (z₂-z₁,z₃-z₁,z₃-z₂); [2,2,0]: quadratic).

[^8]: Agents P1, P2 (independent), W1; `RESULTS-Qside.md` section 0, `RESULTS-xxx-limit.md` X1, `Gaudin-L2-chain.wb`. The rₖ are the same degree-zero homogeneous functions as in XXX: rₖ(z,w)→rₖ(ζ,ω) with O(h) corrections, and that is why the Gaudin vacuum is *exactly* ε-independent in the regime ζ=εφ, ω=εω̂. At ζ=(3/2,5/4,2/3), ω=(5/2,7/3,3/5) the Gauss parameters are r=(4/5,-7/11,65/44), rep-independent.

[^9]: Agents W2/W3 (frame tie, `Gaudin-L2-chain-2.wb` task A, `Gaudin-L2-chain-320.wb` A2) and P4 (`RESULTS-Qside-320.md` T0): with the correct c the per-ζ fits reproduce the stacked kernel to 10⁻⁷¹ to 10⁻⁷⁴; the wrong exponents (1,2,2),(1,3,1),(0,3,2),(2,3,2),(1,3,3),(1,4,2) give residuals 0.18 to 0.69.

[^10]: `Paul/Mathematica/new-conventions/Two-site-rotation-brackets-notes.md`, sections 0 and 5 (Stages D and C), notebooks `L2-two-vacua-untwist.wb` and `L2-210-two-vacua-brackets.wb`; the XXX Q-number reproduces the chain overlap to 7.9·10⁻¹¹¹ and the six adjoint CG data to 10⁻¹⁸.

[^11]: Agent D3, `RESULTS-gaudin-B.md` header (`NB2G`), also D1 `RESULTS-twovac-fit.md` section 0 (the x₀^G direction of P2 mapped with `NB2G` is reproduced to 10⁻⁸¹). Agent W4's data files use 1-based notebook order; P1/P3/P4, D1, U2 use `gl3_reps` order. In this text, slices and indices are in notebook order unless said otherwise.

[^12]: Agents P1 (`RESULTS-Qside.md` Q1; all 64 states of [2,1,0], ODE residual 3.9·10⁻⁵⁵) and P2 (`RESULTS-xxx-limit.md` X4, sympy, symbolic h).

[^13]: Agents P1 (gauge 2, `RESULTS-Qside.md` section 0), P3 (`RESULTS-Qside-310.md` Q1, gauge scan table), P4 (`RESULTS-Qside-320.md` Q1), P2 (X4: the direct Baxter equation is ν²∘adjoint∘ν⁻² exactly). Wrong-sign and wrong-gauge controls have no solutions to degree 8–10.

[^14]: Agents P3, P4. The duality [3,2,0]≅[3,1,0]*⊗det³ exchanges λ₂+1↔λ₁-λ₂+1 between dual and direct opers.

[^15]: Agents D1 (`RESULTS-twovac-fit.md` F1(b): the constant is aₛ/3 and involves neither Hₛ nor the regular parts; the spec's expectation that it would was refuted) and P4 (`RESULTS-Qside-320.md` Q1).

[^16]: Agents P3 (`RESULTS-Qside-310.md` Q1, `frobenius_polys()` in `qside_310.py`) and P4 (`RESULTS-Qside-320.md` Q1); for the adjoint, U1 (`RESULTS-untwisted-Q.md` T2, to 10⁻¹⁴⁸ at 200 digits) and U2. The number of solutions of the system versus the number of states was not counted.

[^17]: Agents P1 (adjoint, F=Π both sides), P3 ([3,1,0]: direct F=Π², dual F′=Π), P4 ([3,2,0]: direct F=Π, dual F′=Π²); plan tables of rounds 2–3 had the dual power wrong for [3,2,0] and were corrected. The constants and exponential factors are the h→0 images of the XXX (1/z_b-1/zₐ)χ₃^(-u/h) at z=e^(hζ).

[^18]: Agents P1 (`RESULTS-Qside.md` Q3), P3 (`RESULTS-Qside-310.md` Q3), P4 (`RESULTS-Qside-320.md` Q3).

[^19]: Agent P3, `RESULTS-Qside-310.md` section R.

[^20]: Agents P1/P3/P4 (exact eigenvalue functions from the exact rational vacua, `RESULTS-Qside*.md` section 0), W2 task C (closed forms, 10 representations) and W3/P4 (correction for [3,2,0]: aₛ^vac=8ω₃ verified at ω₃=3/5,9/7: 24/5, 72/7).

[^21]: Agent W1, `Gaudin-L2-chain.wb`, task G3 (cells “Which T's localise” and the summary cell). The prediction by agent P2 (`RESULTS-xxx-limit.md` X3) from the XXX limit was T_((1,1))(θ₁-2h)=3h²θ₁₂²[1+h X₁/…] with X₁=-4631/124 1+3H₁-S3c[1,2] exactly (the constant at one numerical point only) and ⟨Ω⊗Ω|X₁=u⊗Ω₁, u=Ω₁(431/62-3C_ζ-S3^1-site[-2]); W1 found this to be a combination of the three separate localisers.

[^22]: Agents W2 (`Gaudin-L2-chain-2.wb` task B; X₄ for [3,1,0]), W3 (`Gaudin-L2-chain-320.wb` summary: γ=1, “not c_Λ”). The XXX analogue is presumably the third fused shape T_((3,2)); the identification is a guess from degrees.

[^23]: Agent W3, `Gaudin-L2-chain-320.wb` cells 63–66 and the summary table (verdict “a different relation” from the plan's c_Λ).

[^24]: c_H: agent W2 task C, ten representations (the conjecture 2λ₁ refuted); W3 cell 69. Spectators: W1, W2, W3 (Laurent polynomials in θ₁₂ fitted on six values with three checks).

[^25]: `Paul/Mathematica/new-conventions/Clean-CG-combinations.md`; recapitulated in `Gaudin-sl3-round3-plan.md` section 0 (question Q1).

[^26]: Agents W1 (`Gaudin-L2-chain.wb` task 4 and summary), W2 (task A: (∞,∞) family, closed forms, `Gaudin-L2-chain-2.wb` summary), W3 (`Gaudin-L2-chain-320.wb` A4); also U2, `RESULTS-untwist-limit.md` “Chain-side facts”. The (∞,∞) vacua span only a 4-dimensional subspace of the 8-dimensional one-site space (D1 `RESULTS-twovac-fit.md` Surprise 4), yet the family is CG-complete on HW slices.

[^27]: Agents W3 (`Gaudin-L2-chain-320.wb`, duality row of the summary table) and P4 (`RESULTS-Qside-320.md` Q4).

[^28]: Agent W3, `Gaudin-L2-chain-320.wb`, summary “Cleaning verdict (A4)”, cells 41–57. The exact slice of (4,4,2) has an exactly vanishing component and the sign convention of the HW kets differs between the exact ket and the chain branch there, so the CG cells align signs and use norm-based errors.

[^29]: Agent D2, `RESULTS-twovac-limit.md` L3 (60 digits); the constant between left and right sides is exactly 1, the one predicted from the normalisations, so no leftover scalar exists.

[^30]: Agent D2, `RESULTS-twovac-limit.md` L0–L2; reproduced exactly for the top state at finite h=1/2 (e₁⊗e₁, exact rationals), and independently (h⁶ c Q_num for the top state converges to 102245/41472, to 10⁻¹⁰ by Richardson, against c Q_num^G(top)=2.4653983410).

[^31]: Agent D2, `RESULTS-twovac-limit.md` L4 (a plan item requested by W2).

[^32]: Agent D1, `RESULTS-twovac-fit.md` F2; fixed-point 320-bit column-pivoted QR, a 512×400 problem in 30 s. Larger bases (144, 576 unknowns) fit with null spaces whose sparse representative is the N₅ kernel.

[^33]: Agent P4, `RESULTS-Qside-320.md` section J (hypothesis ladder table, [3,1,0] comparison T7), points P,P′,X₁,X₂,X₃ built at 80 digits; the wrong exponents of c ((1,3,2) for [3,1,0]) give 0.32.

[^34]: Agent D3, `RESULTS-gaudin-B.md` B0 (Fractions, u₀=3/7, u₁=-2/9; nnz of B₀,…,B₆ are 64,64,232,251,247,178,75). The XXX wavefunction check with dual Q's at h=1/2 was skipped as allowed; the rule is not used by the derivation, and the formula built on it is reproduced independently.

[^35]: Agent D3, `RESULTS-gaudin-B.md` B1, B5.

[^36]: Agent D3, `RESULTS-gaudin-B.md` B2, at ζ=P,P′,X₁.

[^37]: Agent D3, `RESULTS-gaudin-B.md` B3. The unipotent gauge freedom inside each filtration level is irrelevant for G (a bilinear form on jets, not on the basis). A determined solve at P is the weakest evidence; the independent evidence is rationality, exact triangularity, ζ-independence and the match with G.

[^38]: Agent D3, `RESULTS-gaudin-B.md` B4. The reproduction of the plan's expectation that the kernel is the “K tensor times wavefunction monomials” was the main goal of round 3; the origin of c is the frame tie Ω_c^∞[1/ω]P_phys⁻¹(z=ζ)=Ω_D^∞[ζ,ω]/κ^G exactly.

[^39]: Agent U2, `RESULTS-untwist-limit.md` L3 and `untwist_limit_frob.py`; independently re-verified by U1 (`RESULTS-untwisted-Q.md` T2: the ν-formulas hold to 10⁻¹⁸⁵ for all six branches, both sites). The result needs no use of the Q polynomials; it was reached by “the Frobenius mechanism”, which is representation-generic.

[^40]: Agent G1, `RESULTS-general-kernel.md` headline 1, K1 (exact Fractions; GT states with multiplicity one only: 6/8, 9/15, 12/24, \dots of the states; the others are not GT eigenvectors and were not used) and section K3(d). The symbol identity holds for λ₁=1,…,5 and c^((λ₂))ₖₗ/c⁽⁰⁾ₖₗ=ε_Λ for all pairs.

[^41]: Agent U2, `RESULTS-untwist-limit.md` headlines 1–4, L2(d). Normalisation-free check: ρ_B(ω̂;ω̂⁰)=[ĉ(ω̂⁰)/ĉ(ω̂)]L(ω̂)/L(ω̂⁰) against the exact ratio from ⟨φᵤ⊗e₁|Ψ₀⟩ at all 16 grid points: maximal deviation 2.5·10⁻¹⁵³, 3.1·10⁻⁵², 4.4·10⁻⁵¹, 3.0·10⁻⁶⁴, 1.2·10⁻⁶⁴, 3.4·10⁻⁸² (φ).

[^42]: Agent U1, `RESULTS-untwisted-Q.md` T2(b)/T3 (`untwisted_Q_laurent_exact.pkl`, n≤5; `untwisted_Q_pert_sym.pkl`, n≤6; `untwisted_Q_octet_fit.pkl`); agent W4, `Gaudin-untwisted-exact.wb` E4, data file `ref/untwisted-exact-210.txt`. U1 and W4 agree exactly to the depth where both have data; the degenerate solver is the missing piece (work package C1 of round 5).

[^43]: Agents U1 (`RESULTS-untwisted-Q.md` T0, Fractions on exact HW kets, five representations, two θ points) and W4 (`Gaudin-untwisted-exact.wb` E2, all 37 sectors of six representations; the 45 branches in Table 23). W4's exact tables: α_copy copy 1 is the smaller, “(3,2,1)”: α_±=(15∓3√5)/2.

[^44]: Agent U1, `RESULTS-untwisted-Q.md` T1 and `v0_bases()` in `untwisted_Q.py`; W4 E3 (exponents and degrees on all 45 branches, `Gaudin-untwisted-exact.wb`). The explicit statement for the top Q's was checked numerically at 200 digits by U1.

[^45]: Agent W4, `Gaudin-untwisted-exact.wb` E3; Table 23 last column. We re-counted the criterion from the exponents in the 45 branches of the data files and obtained 20.

[^46]: `Paul/Mathematica/Gaudin/Gaudin-sl3-round5-plan.md` (written 2026-10-16 about 18:00, not dispatched). Suggested allocation is sonnet agents only, ≤4 Python agents at once, the wolfbook builder alone on the kernel K1.
