# Untwisting the two-rotation det₆×det₂ formula and extracting CG coefficients

gl(3) XXX, L=2, new conventions — research notes of the 2026-10-10 thread

Research notes on the 2026-10-10 thread “untwisting the two-rotation det₆×det₂ formula and extracting CG coefficients” (adjoint [2,1,0]^(⊗2) in full; [3,1,0]^(⊗2), point A, in full; the untwisting of the Q-determinants themselves, adjoint, partially). Audience: Paul Ryan and future Claude sessions. This file is the Markdown companion of `Untwist-two-rotation-notes.tex`/`.pdf`; both have the same content (they are generated from one source).

**Primary sources.** `HANDOFF.md` (entry “NEWEST (2026-10-10)”); the write-ups `L2-210-untwist-CG.md`, `L2-310-untwist-CG.md`, `L2-210-untwist-Qdets.md`; the specs `L2-210-untwist-CG-spec.md`, `L2-310-untwist-CG-spec.md`, `Adjoint-untwist-perturbative-spec.md`, `L2-210-untwist-Qdets-spec.md`; and the four notebooks `L2-210-untwist-CG.wb` (45 cells), `Adjoint-untwist-perturbative.wb` (22), `L2-310-untwist-CG.wb` (90), `L2-210-untwist-Qdets.wb` (45). Background only (not re-documented): `Two-rotations-det-overlap.md` §5.6b, §5.8 (the family), `Undo-site2-rotation-notes.md` §1–§4, §5.2–5.4 (the undo, point A), `Clean-CG-combinations.md` §2–3 (signal and contamination). All paths are under `Paul/Mathematica/new-conventions/` unless stated. Cell numbers are 1-based positions in the `.wb` file. **None of the four notebooks carries a native output**: every number quoted below was read from the markdown result cells that the agents stored under (or after) the code cells, or from the cell sources; no number was recomputed in Wolfram.

**Evidence tags.** Every claim carries one of:

- **[exact]**: a symbolic identity, an exact-rational computation (ranks, `Series`, exact `Solve`), or agreement at 200-digit working precision with the stated residual;
- **[numeric]**: true at the stated numeric point(s) with the stated residual, not proved;
- **[fitted]**: a closed form recognized from a scan and validated on held-out points (the held-out residual is quoted);
- **[conjecture]**: pattern or prediction not (fully) verified; the text says what is verified.

Combined tags such as **[exact/numeric]** mean that part of the statement is exact and the rest numerical; the text says which.

**Scope of this write-up.** No notebook, no existing `.md`/`.tex`, `HANDOFF.md`, `CLAUDE.md` or `paper/*.tex` was edited, nothing was evaluated in Wolfram, nothing was committed. Independent checks made by hand or with plain arithmetic while writing (not notebook cells) are marked “re-derived here” and listed in Appendix B.


---

## 0. Summary

**The question.** The two-rotation family of `Two-rotations-det-overlap.md` §5.6b/§5.8 gives product covectors ⟨ v₁⊗ v₂| whose overlap with a twisted eigenstate Ψ_B of the L=2 transfer matrix is a product of two Q-determinants, det₆×det₂. The previous thread (`Undo-site2-rotation-notes.md`) showed that, after undoing the site-2 rotation, the covectors span every site-1 weight space at finite twist. But a CG coefficient is a property of the highest-weight (HW) ket h_ν, and the eigenstate Ψ_B at finite twist is not h_ν. This thread asks: in which limit of the twist can the HW CG data be read off from the Q-formula alone, with what contamination, and can the contamination be controlled?

**The answer in one paragraph.** Untwist everything together, z=1+εφ, w₁=1+εω₁, w₂=1+εω₂, with the ladder-1 insertion slaved to w₂ (§2). In this regime the undone covector has O(1) components on four site-1 weights and one direction of a fifth, the E₁₃ direction degenerates by exactly one order, and the ratio of contamination to signal K/S is O(ε) for the branches whose CG slice sits at order ε⁰ (“clean”) and a constant independent of ε, ω, φ for those at higher order (“contaminated with a structural constant”). All six HW CG data of adj⊗adj and all ten of [3,1,0]⊗[3,1,0] (point A) are then recovered from the Q-formula, the top-killing determinant, one Q-number τ_c and closed forms in t=θ₁₂/h (§3, §5), with O(ε) errors at finite ε and Richardson limits agreeing with the chain side to 10⁻¹⁸ or better. For the adjoint the two structural constants are derived by first-order perturbation theory (§4); for [3,1,0] they are fitted (§5). The Q-determinants themselves are untwisted only for det₂ (§6).

**What is extracted, by representation.** “Leading” = read at the lowest order of the covector expansion from Qₙᵤₘ alone; “O(ε) layer” = needs the structural constant listed; R = octet-type mixing of a two-dimensional multiplicity space.


### 0.1 Adjoint [2,1,0]^(⊗2) (6 HW states, 5 irreps)

| branch | site-1 weight | order of covector | extraction | tag |
|---|---|---|---|---|
| [4,2,0] (27) | (2,1,0) | ε⁰ | trivial: Qₙᵤₘ/δ^Λ=[Ψₜₒₚ]₁₁ | **[numeric]** |
| [4,1,1] | (2,0,1) | ε⁰ | leading, C=Qₙᵤₘ/(δ^Λ c₀') | **[numeric]** |
| [3,3,0] | (1,2,0) | ε⁰ | leading, C=Qₙᵤₘ/(δ^Λ a₀') | **[numeric]** |
| [3,2,1] octet, 2 copies | (1,1,1) (dim 2) | one direction ε⁰, second ε¹ | functional (u+u')· C_c at leading order; mixing R_c^Q=(Δ_c-slope_c^Q)/(t-2) with Δ_c=(1-t²)+tτ_c/(3(t+3)) | **[fitted]** and **[exact]** (derived, §4) |
| [2,2,2] singlet | (0,1,2) | ε¹ | C=Qₙᵤₘ/(κ φ_((0,1,2))) with κ=2(t-3)/(3(t-2)) | **[fitted]** and **[exact]** (derived) |


### 0.2 [3,1,0]^(⊗2), point A (10 HW states, 8 irreps)

| branch | site-1 weight | order of covector | extraction | tag |
|---|---|---|---|---|
| [6,2,0] | (3,1,0) | ε⁰ | trivial (top) | **[numeric]** |
| [6,1,1] | (3,0,1) | ε⁰ | leading, clean | **[numeric]** |
| [5,3,0] | (2,2,0) | ε⁰ | leading, clean | **[numeric]** |
| [5,2,1]×2 | (2,1,1) (dim 2) | ε⁰ and ε¹ | u'· C at leading order; R_c=3/2-2t+τ_c/(2(t+4)) from the ω₂-slope with signal slopes (t-5)/4, (7-3t)/2 and remainder Δ_c | **[fitted]** (10⁻¹⁸ for R_c; 10⁻⁹ for Δ_c) |
| [4,4,0] | (1,3,0) | ε¹ | κ=(t+1)/2 | **[fitted]** |
| [4,3,1]×2 | (1,2,1) (dim 2) | ε¹ and ε² | κ_c=A+Bτ_c; R_c=A_R+B_Rτ_c; second direction needs α(t) | **[fitted]**, second direction only to 10⁻⁸ |
| [4,2,2] | (1,1,2) (dim 2) | ε¹ and ε² | κ=(11t-31)/(5(3t-7)), R=-2(t-3)/(3(3t-7)), slope constant rational fit | **[fitted]** |
| [3,3,2] | (0,2,2) | ε² | κ=(t+1)(19t-49)/(40(t-2)) | **[fitted]** |

**Closed forms obtained in this thread** (t=θ₁₂/h; τ_c = untwisted eigenvalue of t₂(θ₁+h)/h⁴ of the copy c, a Q-number):

- b'/(a'c')=1+ε(t-2)(φ₂-ω₂)+O(ε²) **[exact]**;
- adjoint singlet: K/S=-t/(3(t-2)), κ=2(t-3)/(3(t-2)) **[exact]**, derived in §4; per copy K_(1,c)/S₁=A+Bτ_c **[exact]**;
- adjoint octets: Δ_c=(1-t²)+tτ_c/(3(t+3)) **[exact]**; τ₁+τ₂=3(2t-1)(t+2), τ₁τ₂=9t(t-1)(t+1)(t+3) **[exact]**; t₂(θ₁+h) on HW_((3,2,1)) is an explicit 2×2 matrix **[exact]**;
- adjoint clean branches, top-killed: K₁ᵗᵏ/S₀=(t-1)(φ₂-φ₃)/4 for [4,1,1] and (t-1)(φ₁-φ₂)/4 for [3,3,0] **[exact]**;
- [3,1,0]: all closed forms of the table above, plus t₂(θ₁+h) on HW_([4,3,1]): trace 8t²+16t-8, determinant 16t(t-1)(t+1)(t+4), and on HW_([5,2,1]): trace 8t²+19t-8, determinant 4(t+4)(t-1)(t+1)(4t+3) **[exact, rational data of the untwisted chain]**;
- det₂⁽¹⁾/minor₁=343000/1161 F(y)·(1,τ₁₀,τ₁₁,τ₂₀,τ₂₁) with F explicit **[numeric]** (residual 10⁻¹⁴⁰); τ⁰_λ = untwisted quantum-Wronskian formula **[numeric]** (residual 10⁻¹⁹⁷);
- det₆ ratios: D₆^([3,3,0])/D₆ᵗᵒᵖ=(w₁-z₁)/(w₁-z₂) and D₆^([4,1,1])/D₆ᵗᵒᵖ=(w₁-z₂)(w₂-z₂)/((w₁-z₃)(w₂-z₃)) exactly in ε **[numeric]** (deviation 10⁻¹²²).

**What is not done.** (i) Analytic untwisting of det₆ for the non-top branches; τ¹ in Q-language (§6). (ii) Closed form of α(t) for the [4,3,1] second direction, and proofs of the [3,1,0] constants (§5). (iii) [4,1,0] at k=2 and [3,2,0]-type representations. (iv) Stored native outputs: all four notebooks must be rerun in VS Code (§7).


---

## 1. Setting, conventions, the family, the numeric points

**Conventions (new, 2026-10-01).** Lax Lᵢⱼ(u)=(u-θₛ)δᵢⱼ+hEᵢⱼ; T(u)=L⁽¹⁾L⁽²⁾ with twist on the right, T^M=T· M; tₐ the sum of principal quantum minors with column k at u-(k-1)h; CBR fusion T_λ(x)=det_(i,j)t_(λ'ⱼ+i-j)(x+(i-1)h), t₀=1, tₐ=0 outside 0..3. Generators: the plain GT generators `JJ` (exact rationals, **not Hermitian**), so the invariant inner product is a diagonal metric G (below) and the HW projector P_(HW) is G-orthogonal, which is the orthogonal projector of the Hermitian EH frame. The GT basis has the HW vector e₁ first, covectors are rows acted on from the right, and Eᵢⱼ with i<j raise (checked in each notebook: `L2-210-untwist-CG.wb` cell 3, `L2-310-untwist-CG.wb` cells 3–4).

**Notation for t.** In this file t:=θ₁₂/h everywhere (so t=62/35 at P1, θ₁₂-h=27/70, θ₁₂-2h=-8/70 in absolute units). This differs from `Undo-site2-rotation-notes.md`, which writes θ₁₂ (31/35 at P1) in its [3,1,0] formulas; see Appendix C. The scans use h=1, θ=(t,0), where the two agree.

**Numeric points.** P1: h=1/2, θ=(2/7,-3/5), t=62/35, y₀=(t-1)/(t-2)=-27/8. P2: h=1/3, θ=(1/5,-4/9), t=29/15, y₀=-14. P3: h=1/2, θ=(3/7,-1/4), t=19/14, y₀=-5/9 (in no fit set). P2 lies in the fit set of the adjoint Δ_c (see Appendix C). Scans in t: h=1, θ=(t,0). Directions φ=(1/3,-2/7,5/11) (primary), φ'=(-2/5,3/8,7/11) (control), and six further directions for the octet scan; ω₁∈{3/4,-5/3,7/2}, ω₂∈{-1/5,9/7,-12/5,5/3,2/9,-7/4} (the 3×6=18 grid, “ω A” =(3/4,-1/5), “ω B” =(7/2,-7/4)); ε∈{1,2,4}·10⁻⁶ (four-point ladder 1,2,4,8·10⁻⁶ where stated); working precision 200 digits. Never let ωᵢ equal a φₖ (poles of Δ and the Gauss parameters).

**Diagonal frame, vacuum, family.** Ket twist diag(z). The one-site vacuum is Ω_w=e₁ᵀe^(r₁E₁₂)e^(r₂E₁₃)e^(r₃E₂₃)=e₁ᵀρ(n) with the representation-independent Gauss parameters (w=(w₁,y,w₃), w₃ drops out)

    r₁=(w₁-z₁)/(w₁-z₂),   r₂=((w₁-z₁)(z₁-y)(z₂-z₃))/((z₁-z₂)(w₁-z₃)(y-z₃)),   r₃=((w₁-z₂)(z₂-y)(z₁-z₃))/((z₁-z₂)(w₁-z₃)(z₃-y))

(`oneSiteVac`, `u12rs`: `L2-210-untwist-CG.wb` cells 5, 11; `L2-310-untwist-CG.wb` cell 7), and n=[[1, r₁, r₂+r₁r₃], [0, 1, r₃], [0, 0, 1]]∈ N₊. The family of `Two-rotations-det-overlap.md` §5.6b is

    v₂=Ω_((w₁,y))=e₁ᵀρ(n),     v₁=v₂ ρ(Δ),   Δₖ=1/(w₂-zₖ),     w₂=(y (t-2))/(t-1) ⇔ y=w₂ y₀,

with overlap, state-independent at finite twist,

    (⟨ v₁⊗ v₂|Ψ_B⟩)/(⟨ x₀|Ψ_B⟩)=c·(det₆[Q_B,Q_(vac,(w₁,w₂))])/(Ψ_Q^B(x₀))·(det₂⁽¹⁾[Q_B;y])/(minor₁(B)) .

For [3,1,0] point A the site-1 covector carries a second-order dressing, v₁∝ e₁ᵀ(1+β E₁₂+γ(E₁₃+1/2E₁₂E₂₃))ρ(n_A)ρ(Δ), y=y_A=w₂y₀, and the overlap is det₆×det₂^((0,1,2))[y_A]/(Ψ_QW(3,4)) (closed forms of β,γ and the identification: `Undo-site2-rotation-notes.md` §5.3–5.4; regression of both here: `L2-310-untwist-CG.wb` cells 16–18, §5). The Q-side objects are those of the background notes: dual Qʲ(x)=zⱼ^(x/h)qʲ(x) with deg qʲ=μⱼ; vacuum direct Qᵢᵛᵃᶜ(x)=wᵢ^(-x/h)qᵢᵛᵃᶜ(x) with degrees (0,2,4) (adjoint) resp. (0,4,6) ([3,1,0]); W(a,b)=Q²(θ₁-ah)Q³(θ₁-bh)-Q³(θ₁-ah)Q²(θ₁-bh); minor₁=W(2,3) (adjoint), W(3,4) ([3,1,0]); Ψ_Q with the vacuum SoV triple. det₆ is the 6×6 bracket determinant over all poles (rows (k,(i,j)) with k=1,2 the insertions 1 and e^(2π ix/h) and (i,j)∈{(1,3),(1,2),(2,3)}; for [3,1,0] an 8-pole version, poles n=0..3 on both ladders); the fixed-degree monic Baxter solvers (overdetermined least squares, no `Chop`) are those of `L2-diagonal-frame-overlap-EH.wb` cell 95, and every Baxter residual is evaluated in exact arithmetic on the stored numbers (`L2-210-untwist-CG.wb` cell 8).

**The undo.** For n∈ N₊, ρ(n)⊗ρ(n) fixes every HW ket, so ⟨ v₁⊗ v₂|h_ν⟩=s⟨φᵤ⊗ e₁ᵀ|h_ν⟩ with the undone covector φᵤ=v₁ρ(n)⁻¹=δ^Λ e₁ᵀg'ρ(n'), n'=Δ⁻¹nΔ n⁻¹=[[1, a', b'], [0, 1, c'], [0, 0, 1]], δₖ=1/(w₂-zₖ), s=(v₂)₁=1, g' the site-1 dressing (g'=1 for the adjoint). Here δ^Λ=δ₁²δ₂ (adjoint) and δ₁³δ₂ ([3,1,0]), and

    a'=(r₁(δ₂-δ₁))/(δ₁),   c'=(r₃(δ₃-δ₂))/(δ₂),   b'=(r₁r₃(δ₁-δ₂))/(δ₁)-((r₂+r₁r₃)(δ₁-δ₃))/(δ₁),   (b')/(a'c')=(y-w₂)/(y-z₂).

All of this is **[exact]** (finite twist; `Undo-site2-rotation-notes.md` §3, re-checked in this thread: adjoint exact undo to 10⁻¹⁵² at 324 points, `L2-210-untwist-CG.wb` cell 25; [3,1,0] to 4.9·10⁻¹³⁸ at 540 points, `L2-310-untwist-CG.wb` cell 32).

**The CG target, signal and contamination.** The twisted eigenstate of the HW-limit branch of ν behaves as Ψ_B→ε⁻ᵏ(h_ν+εψ₁+…) with ⟨ x₀|Ψ_B⟩=1. Define R_C:=⟨ v₁⊗ v₂|Ψ_B⟩, S:=⟨ v₁⊗ v₂|P_(HW)Ψ_B⟩=s⟨φᵤ⊗ e₁ᵀ|P_(HW)Ψ_B⟩ (signal, the CG slice ⟨ m⊗top|h_ν⟩ weighted by the covector), K:=R_C-S (contamination), κ:=R_C/S=1+K/S. The invariant metric is G=diag(1,1,1,2,6,9,4,9) for the adjoint and diag(1,1,2,3,24,32,4,12,24,64,288,576,36,128,576) for [3,1,0] (GT order of `ps`; `L2-210-untwist-CG.wb` cell 25, `L2-310-untwist-CG.wb` cell 4). **Top-killing** (exact): for ν≠ top, ⟨ v₂⊗ v₂|h_ν⟩=0, so T:=δ^Λ s⟨ v₂⊗ v₂|Ψ_B⟩/⟨ x₀|Ψ_B⟩ is pure contamination, and it is itself a Q-number, the det₆ of the (w₁,y)-vacuum over Ψ_Q (verified to 1.6·10⁻¹¹⁶ for the adjoint, `L2-210-untwist-CG.wb` cell 25; 2.5·10⁻⁸² for [3,1,0], `L2-310-untwist-CG.wb` cell 50). The pure Q-side number is

    Qₙᵤₘ(B)=(R_Q(B))/(R_Q(top))·(⟨ v₁⊗ v₂|e₁e₁⟩)/(⟨ x₀|e₁e₁⟩),     R_Q(B)=(det₆[B,F_(vac)])/Ψ_Q·(det₂[B;y])/(minor₁),

and Qₙᵤₘ=R_C is the finite-twist identity at small ε (ε=10⁻⁶ etc.): see §3.2, §5.4.

**Branches.** Adjoint: five HW sectors μ∈{(4,2,0),(4,1,1),(3,3,0),(3,2,1),(2,2,2)} of sizes 1,2,2,6,10 (21 states), multiplicities 1,1,1,2,1: six HW-limit branches, selected as the dimHW_μ states with the smallest raising residual (below 4·10⁻⁷ versus at least 1.3 for the others at ε=10⁻⁶), copies ordered by the real part of the eigenvalue of t₁+π t₂ at u=3/7. [3,1,0]: eight HW sectors of sizes 1,2,2,6,3,10,14,18 (56 states), multiplicities 1,1,1,2,1,2,1,1, ten branches; the raising residual is ∼10⁻⁷ versus at least 1.41 for the others (`L2-310-untwist-CG.wb` cell 27).


---

## 2. The untwisting regime and why


### 2.1 The problem, and why fixed y does not work

At finite twist the covector ⟨ v₁⊗ v₂| sees the eigenstate Ψ_B of the twisted chain, not the HW ket h_ν of the untwisted one. As z→1 the HW-limit branch behaves as Ψ_B=ε⁻ᵏ(h_ν+εψ₁+…) and the overlap splits exactly as R_C=S+K (signal plus contamination, §1; mechanism and measured κ values in `Clean-CG-combinations.md` §2). The question is therefore not whether the covectors span (the previous thread showed they do, at finite twist) but at which order in ε each CG slice ⟨ m⊗top|h_ν⟩ appears in the undone covector, and how it compares with the contamination that comes from the O(ε) corrections ψ₁,ψ₂,… of the eigenstate.

At **fixed y** (the ladder-1 insertion and w₁ fixed, z→1) the resolvent torus Δ=diag(1/(w₂-zₖ)) becomes a scalar at z=1, the two rotations v₁,v₂ coincide in the limit (spec §4: in the joint regime they “stay distinct” because the torus is non-scalar), and the undone covector has components O(εᵈ) where d is the depth of the weight. For the adjoint the leading orders at fixed y are ε⁰,ε¹,ε¹,ε² (rank 2), ε³,ε³,ε⁴ over the seven weights (`Undo-site2-rotation-notes.md` §3, computed in `L2-rotation-sov.wb` cell 119 with a',c'=O(ε), b',a'c'=O(ε²)) **[exact]**. The contamination, which starts at relative order ε, is then at or above the signal for every non-top state: “at fixed y nothing is clean” (HANDOFF; spec `L2-210-untwist-CG-spec.md` §4). That table is not used in this thread; I did not re-derive it here.


### 2.2 The joint regime

The regime was chosen by the main session before the notebooks were written (spec §0, §4: sympy-checked at P1, then verified in `L2-210-untwist-CG.wb` cells 10–12): untwist **everything together**,

    z=1+εφ,     w₁=1+εω₁,     w₂=1+εω₂,     y=w₂ y₀,   y₀=(t-1)/(t-2),

with fixed rational directions φ,ω₁,ω₂ and ε→0 (rational paths; the direction at first order is what matters). The product condition w₂(y)=w₂ then holds exactly. The resolvent torus is now Δ=ε⁻¹diag(1/(ω₂-φₖ)), a genuine non-scalar torus element at ε=0 after removing the overall factor ε⁻¹, so the two rotations stay distinct. w₃ drops out of the determinant identity; it is set to 1+5ε/13 (the spec's 1+2ε/9 coincides with w₂ at ω₂=2/9, see §7).


### 2.3 The mechanism: a', c' and the degeneration of the E₁₃ direction

With δ₂/δ₁=(w₂-z₁)/(w₂-z₂) one has δ₂/δ₁-1=(z₂-z₁)/(w₂-z₂), hence, in the regime,

    a'=((w₁-z₁)(z₂-z₁))/((w₁-z₂)(w₂-z₂))=((ω₁-φ₁)(φ₂-φ₁))/((ω₁-φ₂)(ω₂-φ₂)),

**exactly independent of ε** **[exact]** (re-derived here; confirmed symbolically in `L2-210-untwist-CG.wb` cell 11, with a'₀=-2275/783 at φ, ω A). Likewise c'=r₃(ε)(φ₃-φ₂)/(ω₂-φ₃)=O(1), with r₃(0)=((ω₁-φ₂)(φ₁-φ₃))/((φ₁-φ₂)(ω₁-φ₃)) and c'₀=2755/3549 at φ, ω A (`L2-210-untwist-CG.wb` cells 11, 45). The decisive quantity is the ratio

    (b')/(a'c')=(y-w₂)/(y-z₂)=(w₂)/(w₂(t-1)-(t-2)z₂)=(1+εω₂)/(1+ε((t-1)ω₂-(t-2)φ₂))=1+ε (t-2)(φ₂-ω₂)+O(ε²),

**[exact]** (the middle expression is re-derived here from the finite-twist identity b'/(a'c')=(y-w₂)/(y-z₂) of the background notes and y=w₂(t-1)/(t-2); the first-order coefficient is confirmed in `L2-210-untwist-CG.wb` cell 11 as 24/1225=(t-2)(φ₂-ω₂) at t=62/35, φ₂-ω₂=-3/35, and symbolically in `Adjoint-untwist-perturbative.wb` cell 22). So b'₀=a'₀c'₀ and at ε=0 the undone covector is

    φᵤ(0)/δ^Λ=e₁ᵀe^(a'₀E₁₂)e^(c'₀E₂₃)    (adjoint),

**[exact]** (`L2-210-untwist-CG.wb` cell 11: exact equality; the reversed order of the exponentials fails). The E₁₃ direction (b'-a'c') is thus **one order down**: it is the ε=0 shadow of the old blind spot b'=0⇔ω=y. The main session's reading (spec §4) is that the vacuum's r₂/(r₁r₃) and the torus Δ are both built from the same φ, so that n₀Δ n₀⁻¹ has b'=a'c'; I take that as an interpretation, the all-orders formula above is the proven statement. (Ordering of exponentials: as a matrix n'=e^(c'E₂₃)e^(b'E₁₃)e^(a'E₁₂), the Undo notes' form; at b'=a'c' this equals e^(a'E₁₂)e^(c'E₂₃) because [E₁₂,E₂₃]=E₁₃. The two statements are consistent.)


### 2.4 What degenerates: per-weight orders

Adjoint, leading order and ranks of the coefficient vectors over the 3×6 (ω₁,ω₂) grid (exact, `L2-210-untwist-CG.wb` cells 11–12; identical for φ and φ') **[exact]**:

| site-1 weight | dim | leading order | ranks (ε⁰, ε¹, both) |
|---|---|---|---|
| (2,1,0) | 1 | ε⁰ | (1,0,1) |
| (1,2,0) | 1 | ε⁰ | (1,0,1) |
| (2,0,1) | 1 | ε⁰ | (1,1,1) |
| (1,1,1) | 2 | ε⁰ (one direction) | (1,2,2) |
| (1,0,2) | 1 | ε⁰ | (1,1,1) |
| (0,2,1) | 1 | ε¹ | (0,1,1) |
| (0,1,2) | 1 | ε¹ | (0,1,1) |

The (1,1,1) component is b'u+a'c'u' with u=e₁ᵀE₁₃=(-1,3) and u'=e₁ᵀE₂₃E₁₂=(2,0) on the two GT indices; at ε=0 only u+u'=(1,3) is seen. The (0,1,2) component is 3b'(a'c'-b')δ^Λ=-3(t-2)(φ₂-ω₂)(a'₀c'₀)²δ^Λ ε+O(ε²) **[exact]** (the first-order coefficient re-derived here from the ratio above; notebook value -72/1225 (a'₀c'₀)² at ω A, cell 11). So in the adjoint the singlet direction and the second octet direction sit at O(ε), at the same order as the first contamination.

[3,1,0] point A (`L2-310-untwist-CG.wb` cells 20–24, exact, symbolic ε; again identical for φ,φ') **[exact]**:

| site-1 weight | dim | leading order | ranks (ε⁰, ε¹, ε²; cumulative 0..1, 0..2) |
|---|---|---|---|
| (3,1,0) | 1 | ε⁰ | 1,0,0; 1,1 |
| (3,0,1), (2,2,0), (2,0,2) | 1 each | ε⁰ | 1,1,1; 1,1 |
| (2,1,1) | 2 | one direction at ε⁰ (u'=e₁ᵀE₁₂E₂₃=(2,8)), second at ε¹ (u=e₁ᵀE₁₃=(-1,8)) | 1,2,2; 2,2 |
| (1,3,0) | 1 | ε¹ | 0,1,1; 1,1 |
| (1,2,1), (1,1,2) | 2 each | one direction at ε¹, second at ε² | 0,1,2; 1,2 |
| (1,0,3), (0,3,1), (0,2,2) | 1 each | ε² | 0,0,1; 0,1 |
| (0,1,3) | 1 | beyond ε² (ε³) | 0,0,0; 0,0 |

The slice weight carrying the CG vector of an irrep ν is ν-Λ; for [3,1,0] these are (3,1,0), (3,0,1), (2,2,0), (2,1,1), (1,3,0), (1,2,1), (1,1,2), (0,2,2) for [6,2,0], [6,1,1], [5,3,0], [5,2,1]×2, [4,4,0], [4,3,1]×2, [4,2,2], [3,3,2] respectively. The weight (0,1,3) is not a slice weight.


### 2.5 Why three orders for [3,1,0] and two for the adjoint

The dressing of point A, e₁ᵀ(1+β E₁₂+γ(E₁₃+1/2E₁₂E₂₃)), behaves in the regime as follows (`L2-310-untwist-CG.wb` cells 8, 21, 24):

- γ=O(ε) (numerator carries three ε's, denominator two; γ(10⁻⁶)/γ(2·10⁻⁶)=0.5000003; ε¹ coefficient -2565/56056 at φ, ω A) **[exact]**;
- β(0)=((φ₁-φ₂)(ω₁-φ₁))/(2(ω₁-φ₂)(ω₂-φ₁)) **[exact]**; with the torus rescaling of the E₁₂ coefficient, β₀ δ₂/δ₁=β₀(ω₂-φ₁)/(ω₂-φ₂)=-a'₀/2 (re-derived here), so that the undone covector at ε=0 is

    φ_A(0)/δ^Λ=e₁ᵀ(1-(a'₀)/2E₁₂)e^(a'₀E₁₂)e^(c'₀E₂₃)=e₁ᵀ(1+(a'₀)/2E₁₂)e^(c'₀E₂₃)

**[exact]** (cell 21: exact equality; the no-dressing variant and the reversed order fail). The second equality uses e₁ᵀE₁₂³=0: (1-a/2X)e^(aX)=1+a/2X-(a³)/12X³+… (re-derived here). **The E₁₂² term cancels exactly**, which is why the weight (1,3,0) drops from ε⁰ to ε¹.

Consequently the CG slices of [3,1,0] enter the covector at orders ε⁰ (four 1-dimensional weights and one direction of (2,1,1)), ε¹ ((1,3,0), one direction each of (1,2,1),(1,1,2), the second direction of (2,1,1)) and ε² ((0,2,2) and the second directions of (1,2,1),(1,1,2)). The Q-number must therefore be expanded to three orders, against two for the adjoint (orders ε⁰ and ε¹; the highest slice weight is (0,1,2) at ε¹). The sources do not derive why the contamination sits at the same order as the signal for the contaminated branches; the data are as follows (§3.1, §5.3). A branch whose slice weight is at ε⁰ is clean (K/S=O(ε)). A branch whose slice weight (first direction, for the dim-2 weights) is at ε¹ or ε² has K/S equal to a constant independent of ε, ω, φ. For the dim-2 weights the mixing R of the two directions is read from the next layer of the expansion (§3.3, §5.5).

For [4,1,0] at k=2 (dressing of E₁₃-degree 2) the [3,1,0] write-up expects four orders **[conjecture]**; it was not tried.


### 2.6 Bookkeeping of ε-exponents

For every branch the exponent of S equals (exponent of δ^Λ) + (exponent of ‖Ψ‖, as a power of ε⁻¹) + (leading order of the covector at the slice weight), with δ^Λ∼ε⁻³ (adjoint) or ε⁻⁴ ([3,1,0]). Adjoint: [4,2,0]: -3-2+0=-5; [4,1,1], [3,3,0]: -3-3+0=-6; octets: -3-4+0=-7; singlet: -3-6+1=-8, matching the measured exponents of §3.1. [3,1,0]: [4,4,0]: -4-4+1=-7; [3,3,2]: -4-7+2=-9 (stated in `L2-310-untwist-CG.wb` cell 32), and likewise [6,1,1]: -4-3+0=-7, [4,2,2]: -4-6+1=-9 **[exact]** (consistency of the tables, re-derived here).


---

## 3. Adjoint [2,1,0]^(⊗2)

Notebook: `L2-210-untwist-CG.wb` (45 cells, no stored outputs; the agent ran each cell's source in the shared kernel K1, results in the markdown cells 25, 33, 40, 45). Implementation by a sonnet wolfbook-builder from `L2-210-untwist-CG-spec.md`; planned and reviewed by the main session.


### 3.1 Chain side: signal and contamination

Task 1 (cells 13–16), three ε and 18 ω at φ, P1. The exact undo S=s⟨φᵤ⊗ e₁ᵀ|P_(HW)Ψ⟩ holds to better than 10⁻¹⁵² at 324 points **[exact/numeric]**. ε-exponents (identical for all 18 ω) of S, K, ‖Ψ‖ and the contamination ratio:

| branch | S | K | \|Ψ\| | max\| K/S\| at ε=1,2,4·10⁻⁶ | κ-1 at ω A |
|---|---|---|---|---|---|
| [4,2,0] | -5 | K≡0 | -2 | 0 | 0 |
| [4,1,1] | -6 | -5 | -3 | 9.6·10⁻⁷, 1.9·10⁻⁶, 3.8·10⁻⁶ | 1.1·10⁻⁷, 2.2·10⁻⁷, 4.4·10⁻⁷ |
| [3,3,0] | -6 | -5 | -3 | 9.4·10⁻⁷, 1.9·10⁻⁶, 3.7·10⁻⁶ | 8.6·10⁻⁸, 1.7·10⁻⁷, 3.5·10⁻⁷ |
| octet 1 | -7 | -6 | -4 | 2.8·10⁻⁶, 5.5·10⁻⁶, 1.1·10⁻⁵ | 1.5·10⁻⁷, 3.0·10⁻⁷, 6.0·10⁻⁷ |
| octet 2 | -7 | -6 | -4 | 1.1·10⁻⁶, 2.3·10⁻⁶, 4.5·10⁻⁶ | -1.0·10⁻⁷, -2.0·10⁻⁷, -4.0·10⁻⁷ |
| singlet | -8 | -8 | -6 | 2.58 at every ε and ω (=31/12) | κ=3.58 constant (=43/12) |

**[numeric]** The spec's expectation is confirmed: K/S=O(ε), linear in ε, for [4,1,1], [3,3,0] and both octets, and κ↛1 for the singlet. The contamination carriers (undone frame, ε=10⁻⁶, ω A; Σ_ν c_ν=K to 10⁻¹⁶ or better): singlet, site-1 weight (1,0,2) at ε⁻⁸ (fraction 1.0), (0,1,2),(0,2,1),(1,1,1) at ε⁻⁷ (10⁻⁷ each), (1,2,0),(2,0,1) at ε⁻⁶, (2,1,0) at ε⁻⁵; octet 1, weights (1,1,1),(1,2,0),(2,0,1) all at ε⁻⁶ (fractions 0.94,2.31,0.37, signed) and (2,1,0) at ε⁻⁵ (7·10⁻⁸) **[numeric]**. In words: the singlet's contamination is carried by the weight (1,0,2), which is O(1) in the covector, meeting the O(ε) correction ψ₁ of the singlet eigenstate, while its own slice weight (0,1,2) is O(ε).


### 3.2 Q side at small ε and the leading-order CG numbers

Task 2 (cells 17–21) **[numeric]**. Dual Q's of the branches (deg qʲ=μⱼ) and vacuum direct Q's of the twist (w₁,w₂,w₃) (degrees 0,2,4) by fixed-degree least squares; exact Baxter residuals 8.4·10⁻¹⁴⁹ (dual) and 6.9·10⁻¹⁹⁰ (vacuum). The finite-twist identity at small ε: R_C/R_Q is state-independent over the six branches to 1.56·10⁻¹¹⁶ (ε=10⁻⁶) and 1.34·10⁻¹²⁵ (4·10⁻⁶) (max over 18 ω), and over all 21 states of the five HW sectors to 2.4·10⁻¹¹⁶, 1.0·10⁻¹²⁶; R_Q keeps 126 digits (the limit is the working precision, not the identity). Qₙᵤₘ=R_C to the same numbers.

Leading-order extraction. With the ε=0 covector of §2.3, C^Q(ω,ε):=Qₙᵤₘ/(δ^Λ fac(ω)) with fac=c'₀ ([4,1,1], weight (2,0,1)), a'₀ ([3,3,0], weight (1,2,0)) and a'₀c'₀ (octets, which gives the functional (u+u')· C_c on the (1,1,1) slice, u+u'=(1,3) independent of ω); for the 27 the check is Qₙᵤₘ/δ^Λ=[Ψₜₒₚ]₁₁ (deviation below 10⁻¹⁴⁴). The chain side is C_B=f·[P_(HW)Ψ_B]_((ν-Λ),e₁). Spread of C^Q over the 18 ω and error against the chain side, both O(ε):

| branch | ε=10⁻⁶: spread / error | ε=4·10⁻⁶: spread / error |
|---|---|---|
| [4,1,1] | 8.5·10⁻⁷ / 7.9·10⁻⁷ | 3.4·10⁻⁶ / 3.2·10⁻⁶ |
| [3,3,0] | 8.5·10⁻⁷ / 9.4·10⁻⁷ | 3.4·10⁻⁶ / 3.7·10⁻⁶ |
| octet 1 | 3.5·10⁻⁶ / 3.4·10⁻⁶ | 1.4·10⁻⁵ / 1.4·10⁻⁵ |
| octet 2 | 1.6·10⁻⁶ / 1.8·10⁻⁶ | 6.3·10⁻⁶ / 7.3·10⁻⁶ |

(scaling 4.00±0.00001 in ε). Second-order Richardson on f(ε)=C(ε)ε^(-p_B) with p_B= the exponent of ‖Ψ_B‖ (2,3,3,4,4), ε=1,2,4·10⁻⁶ at ω A (cell 21): chain and Q limits 2.609467456 (27), 46.67453891 ([4,1,1]; agreement 1.2·10⁻¹⁹), 10.92853892 ([3,3,0]; 1.6·10⁻¹²¹), 184.3305612 (octet 1, functional; 3.2·10⁻²⁰), -238.5587523 (octet 2; 9.5·10⁻²⁰) **[numeric]**. The two octet functionals are (u+u')· C_c only: the mixing R_c requires the O(ε) layer.


### 3.3 The O(ε) layer

**Top-killing** (cells 22–24, 30). T=δ^Λ s⟨ v₂⊗ v₂|Ψ_B⟩/⟨ x₀|Ψ_B⟩ is subtracted: Qₙᵤₘᵗᵏ:=Qₙᵤₘ-T=R_C-T. The identity ⟨ v₂⊗ v₂|Ψ⟩/⟨ x₀|Ψ⟩=det₆[F_(vac)(w₁,y)]/Ψ_Q (top-normalised) holds to 1.6·10⁻¹¹⁶ **[numeric]**.

**Measurements before the fits** (Task 3, cell 25, with ρ_B:=Kᵗᵏ/S): for the octets ρ_B=O(ε), so the ε-independent quantity is g=ρ_B/ε (ε-independence 3.0·10⁻⁶ and 1.7·10⁻⁶ for copies 1, 2, O(ε)); for the singlet ρ_B is ε-independent (6.1·10⁻⁶). The singlet K/S=2.5833=31/12 is independent of ε, ω and φ (g(φ')/g(φ)=1-3·10⁻⁸). The octet g₀=lim(4g(10⁻⁶)-g(4·10⁻⁶))/3 takes six distinct values over the 18 ω: it depends on ω₂ only (the ω₁-only ratios equal 1), and it is **not** φ-independent. These are the facts the later fits are built on.

**Singlet** **[fitted]**, then **[exact]**. K/S is constant in ε, ω, φ. A scan in t at h=1, θ=(t,0) over eight values t∈{29/15,5/2,3/2,7/3,11/4,13/5,9/4,17/6} (fresh kernel at t=5/2 and at P2; φ' at t=7/3), residual 3.1·10⁻¹¹ (cell 33), gives

    K/S=-t/3(t-2),     κ_(singlet)=1+K/S=2(t-3)/3(t-2)  (=43/12 at P1).

The sources give the structural reason (main session): the first-order correction ψ₁ of the untwisted eigenvector comes from the twist derivative, whose off-diagonal part within a weight sector is h²ΣⱼE⁽¹⁾ₖⱼE⁽²⁾ⱼₖ, independent of θ, while the singlet signal carries the explicit factor (t-2) of §2.3. This is made rigorous in §4. The Q-side CG number: C^Q=Qₙᵤₘ/(κ φ_((0,1,2))(ε,ω)) with φ_((0,1,2)) the exact (0,1,2) component; spread over 18 ω / error against the chain side 1.27·10⁻⁶ / 1.47·10⁻⁶ at ε=10⁻⁶ and 5.06·10⁻⁶ / 5.89·10⁻⁶ at 4·10⁻⁶; Richardson limit of C^Qε⁶: chain and Q -5118.53831440474756, agreement 4.5·10⁻¹⁹ **[numeric]**.

**Octets** (cells 29–32, 37–39). Define F:=Qₙᵤₘᵗᵏ/(δ^Λ a'₀c'₀)=F₀+ε F₁+… (computed as f=Fε⁴ from three ε by an exact 3×3 inverse; two points would leave an O(ε) error in F₁). Measured at P1 over eight directions φ **[numeric]**:

- F₀ is ω-independent (4.8·10⁻¹⁷, max over 8 φ); F₀ itself (=(u+u')· C_c) depends on φ (e.g. 184.33, 38.96, 120.66, -37.60 for copy 1 at the first four φ; F₀⁽²⁾/F₀⁽¹⁾=-1.29418);
- g=F₁/F₀ is exactly a straight line in ω₂ (residual 2.1·10⁻¹¹ copy 1, 9.6·10⁻¹² copy 2; the residual is O(ε²)), ω₁-independent (slope and intercept differences below 10⁻¹⁶⁸), with a **φ-independent slope**: copy 1 -1.5687193188 (spread over 8 φ 8.4·10⁻¹²), copy 2 0.7163668808 (2·10⁻¹²), and a φ-dependent intercept;
- the chain-side mixing R_c=(u· C_c)/((u+u')· C_c): ≈-1.694367 (copy 1; φ-spread 4·10⁻⁷), 1.0955643 (copy 2).

The signal contributes exactly -(t-2)R_c to the slope (proved in §4: the signal part of the slope is exactly -(t-2)R_c, remainder 0); the rest is a structural constant

    Δ_c:=slope_c+(t-2)R_c,

independent of ω and φ: Δ₁∈[-1.18143555,-1.18143546], Δ₂∈[0.46595216,0.46595223] at P1 (spread 9·10⁻⁸, 6·10⁻⁸, O(ε)) **[numeric]**. **A scan over 22 values of t** (15 in the fit set, 7 held out; h=1, θ=(t,0); copies sorted by τ; τ_c the eigenvalue of t₂(θ₁+h) at u=θ₁+h read from the branch's own t₂ polynomial, Richardson-extrapolated, and divided by h⁴, the h-independent scale because t₂ has degree 4 in (u,θ,h)) gives (cells 38–40):

| series | fitted form | held-out relative residual |
|---|---|---|
| τ₁+τ₂ | 6t²+9t-6=3(2t-1)(t+2) | 2.4·10⁻²¹ |
| τ₁τ₂ | 9t⁴+27t³-9t²-27t=9t(t-1)(t+1)(t+3) | 2.1·10⁻²² |
| Δ₁+Δ₂ | 3(2-t²)/(t+3) | 4.8·10⁻¹² |
| Δ₁Δ₂ | (t-3)(t-1)(t+1)/(t+3) | 1.6·10⁻¹¹ |
| A in Δ_c=A+Bτ_c | 1-t² | 2.2·10⁻¹¹ |
| B in Δ_c=A+Bτ_c | t/(3(t+3)) | 3.7·10⁻¹¹ |

**[fitted]** so that

    Δ_c=(1-t²)+(t τ_c)/3(t+3),     τ²-3(2t-1)(t+2) τ+9t(t-1)(t+1)(t+3)=0 .

The fit-set data are in cell 40 (e.g. t=29/15: Δ=-1.46172486, 0.40496811, τ=9.76840506, 24.0582616; t=5/2: -2.50844796, 0.19026615, 18.0942434, 35.9057566; t=17/6: -3.16347254, 0.06347254, 23.8677676, 43.7988990). Consistency re-derived here: 2A+B(τ₁+τ₂)=3(2-t²)/(t+3) and Δ₁Δ₂=(t-3)(t-1)(t+1)/(t+3) follow algebraically from the quoted τ-sums and A,B. The hypothesis Δ₁=Δ₂ is false (-2.508 versus 0.190 at t=5/2).

**The Q dictionary for τ_c.** τ_c is a Q-number: t₂(θ₁+h)=χ₃ ν₁(θ₁+h) ν₂(θ₁) Qʲ(θ₁-h)/Qʲ(θ₁), for any j (j-spread 0), with χ₃=z₁z₂z₃ and νₙ(x)=Πₛ(x-θₛ+hΛₙ) (cell 42; `Clean-CG-combinations.md` §4); the τ_c from this dictionary equal the eigenvalue-based ones to relative difference 0 at P1, P2, P3 **[numeric]**. Hence the mixing is read from the Q side as

    R_c^Q=(Δ_c-slope_c^Q)/(t-2),

with slope_c^Q the ω₂-slope of F₁/F₀ computed from Qₙᵤₘ-T alone.


### 3.4 Final extraction at P1, P2, P3

Cells 41–45 (`finalRun`). Per branch: error =max_ω|C^Q/C-1| at ε=10⁻⁶ (the expected O(ε) of the leading-order formula), then the Richardson limit of Cεᵖ (p=2,3,3,4,4,6) from ε=1,2,4·10⁻⁶, chain-side limit, and |Q/chain-1| **[numeric]**:

| branch | P1 | P2 | P3 (out of sample) |
|---|---|---|---|
| 27 | 0 / 2.609467456 / 0 | same | same |
| [4,1,1] | 7.9·10⁻⁷ / 46.67453891 / 1.2·10⁻¹⁹ | 1.1·10⁻⁶ / 38.57793522 / 6.0·10⁻²⁰ | 7.6·10⁻⁷ / 100.8170040 / 1.8·10⁻¹⁹ |
| [3,3,0] | 9.4·10⁻⁷ / 10.92853892 / 1.6·10⁻¹²¹ | 1.1·10⁻⁶ / 9.032771962 / 9.7·10⁻¹²¹ | 4.3·10⁻⁷ / 23.60564406 / 3.5·10⁻¹²¹ |
| octet 1, (u+u')· C | 3.4·10⁻⁶ / 184.3305612 / 3.2·10⁻²⁰ | 3.6·10⁻⁶ / 122.3424387 / 2.8·10⁻²⁰ | 3.1·10⁻⁶ / 937.1257171 / 3.8·10⁻¹⁹ |
| octet 2, (u+u')· C | 1.8·10⁻⁶ / -238.5587523 / 9.5·10⁻²⁰ | 1.2·10⁻⁶ / -245.4884068 / 1.6·10⁻²⁰ | 3.6·10⁻⁶ / -254.7138085 / 2.9·10⁻¹⁹ |
| singlet | 1.5·10⁻⁶ / -5118.538314 / 4.5·10⁻¹⁹ | 2.0·10⁻⁶ / -3805.202777 / 6.4·10⁻¹⁹ | 1.1·10⁻⁶ / -20997.70633 / 2.2·10⁻¹⁹ |

(Per-cell entries: error at ε=10⁻⁶ / Q-side limit (equal to the chain limit to the digits shown) / |Q/chain-1|. The octet limits are of the functional (u+u')· C_c.) Octet mixing, Q side against chain side:

| point | τ_Q copy 1, 2 | Δ_c | slope^Q | R^Q=Rᶜʰᵃⁱⁿ | \| R^Q-Rᶜʰᵃⁱⁿ\| |
|---|---|---|---|---|---|
| P1 | 7.72932830, 21.0412839 | -1.18143552, 0.465952203 | -1.56871932, 0.716366881 | -1.69436661, 1.09556421 | 4.9·10⁻¹¹, 3.7·10⁻¹¹ |
| P2 | 9.76840506, 24.0582616 | -1.46172486, 0.404968107 | -1.58305272, 0.476746417 | -1.81991789, 1.07667465 | 8.4·10⁻¹¹, 5.6·10⁻¹¹ |
| P3 | 3.18098696, 14.0843192 | -0.511570329, 0.620469626 | -1.37387704, 1.36685128 | -1.34136599, 1.16103812 | 1.2·10⁻¹¹, 5.3·10⁻¹² |

The differences are at the Richardson O(ε³) level. The chain ε-exponents of the slices are -2,-3,-3,-4,-4,-6 at all three points. **All six HW CG data of adj⊗adj are thus obtained from the det₆×det₂ formula, the top-killing det₆, one Q-number τ_c and the two closed forms κ(t), Δ_c(t,τ_c).** Caveat on “from the Q-formula alone”: the closed forms were fitted from chain-side scans (and later derived, §4); at the final stage the chain side enters only the comparison columns. P2 lies in the Δ_c fit set (cell 43 states this), so only P3 and the seven held-out t are independent tests of Δ_c.


---

## 4. The perturbative derivation of the adjoint constants

Notebook: `Adjoint-untwist-perturbative.wb` (22 cells; spec `Adjoint-untwist-perturbative-spec.md`; no native outputs, the cells were run in K1, results in the markdown cell 22). The model is `Adjoint-clean-analytic.wb` (first-order perturbation theory of the eigenvector, Ψ=h+εψ₁, for the T-insertion family). All results of this section are exact symbolic identities in (t,φ₁,φ₂,φ₃,ω₁,ω₂) (octets: in (m,φ,ω), m the conic parameter below) **[exact]**, except the finite-ε cross-check (80 digits, **[numeric]**). Setting: h=1, θ=(t,0) with **t symbolic**, plain GT generators, exact arithmetic.


### 4.1 Setup

**Covectors.** V₁:=v₁/δ^Λ=v₂ ρ(diag(Dₖ))/(D₁²D₂) with Dₖ=1/(ω₂-φₖ); the only ε-dependence is in v₂, through the Gauss parameters rₖ(ε) (cell 5). The setup checks (listed in the results cell 22): gl(3) relations, positivity of G=diag(1,1,1,2,6,9,4,9), HW multiplicities (1,1,1,2,1), sector sizes (1,2,2,6,10); [tₐ(u),t_b(v)]=0 and the first-order identity [tₐ₀,t_b']+[tₐ',t_(b0)]=0 (symbolic in t,φ); exact undo at a rational point and ρ(n)⊗ρ(n) HW=HW; φ̃(0)=e₁ᵀe^(a₀'E₁₂)e^(c₀'E₂₃), a' ε-free, b'/(a'c')=1+ε(t-2)(φ₂-ω₂)+O(ε²), all symbolic.

**The twisted transfer matrix to first order.** Since tₐ(u;z)=Σ_(|I|=a)z_I m_I(u) with z_I=Π_(i∈ I)zᵢ and the untwisted principal quantum minors m_I (the principal minors preserve the total weight, so every term stays inside a weight sector), the ε-derivative is exact:

    tₐ'=Σ_(|I|=a)(Σ_(i∈ I)φᵢ)m_I,     tₐ(u;ε)=tₐ₀+ε tₐ'+O(ε²)

(cell 4). **The family member.** A generic member M=t₁(3/7)+5/3t₂(1/5) of the commuting family restricted to a weight sector (a second member t₁(2/9)-7/4t₂(-3/5) is used as a check). **Perturbation theory.** M₀ is gl(3)-invariant but not symmetric in the plain basis, so use the left eigenvector ℓ of the eigenvalue m₀ of the HW ket h:

    (M₀-m₀)ψ₁=-(M₁-m₁)h,     m₁=(ℓ M₁h)/(ℓ h),     ℓ·ψ₁=0,

solved as the nonsingular system (M₀-m₀+hℓ/(ℓ h))ψ₁=-(M₁-m₁)h (cell 7, `pert1`; ψ₁=Σᵢφᵢψ₁⁽ⁱ⁾ over the unit directions; the component of ψ₁ along h is irrelevant). Irrep projectors Π_ν on a weight sector are built from C₂+C₃/7 (`PiIrrep`). The overlap to first order is

    ⟨ v₁⊗ v₂|Ψ⟩=⟨ V₁(ε)⊗ v₂(ε)|h⟩+ε ⟨ V₁⁰⊗ v₂⁰|ψ₁⟩+O(ε²),

with the scalar normalisation series of Ψ cancelling in the ratios below. Define S₁ as the order-ε coefficient of the signal ⟨ V₁(ε)⊗ v₂|h⟩ (zero at ε⁰ for the singlet) and K₁:=⟨ V₁⁰⊗ v₂⁰|ψ₁⟩.


### 4.2 Singlet (Claim 1)

Sector (2,2,2) (dim 10), h the gl(3)-invariant singlet ket (checked: annihilated by all Eᵢⱼ, i≠ j, and Eᵢᵢ↦2; cell 9).

    S₁=(3(t-2)(ω₁-φ₁)²(φ₁-φ₃)²(φ₂-φ₃)²)/((ω₂-φ₂)(ω₁-φ₃)²(ω₂-φ₃)²)=-3(t-2)(φ₂-ω₂)(a₀'c₀')² C_(singlet),     (K₁)/(S₁)=-t/3(t-2),

**[exact]**, for all φ,ω₁,ω₂ (C_(singlet)=1 in the normalisation used; the equality of the two expressions for S₁ is re-derived here from a₀',c₀' of §2.3), independent of the family member (second member: K₁-K₁'=0 and ψ₁-ψ₁'∥ h). Hence κ_(singlet)=1+K/S=2(t-3)/(3(t-2)), =43/12 at t=62/35, reproducing the fitted law of §3.3.

**Which states contribute.** The twist derivative t'=Σ_I(Σ_(i∈ I)φᵢ)m_I is an adjoint tensor operator (the sources' argument; φ↦ t' is linear in the Cartan element diag(φ)), so it takes the singlet only to octets: the projections of ψ₁ onto the 27, [4,1,1] and [3,3,0] vanish identically **[exact]**. **100% of K₁ comes from the two octets** (each contributing a 2-dimensional weight-(2,2,2) space), with copy shares

    (K_(1,c))/(S₁)=A(t)+B(t) τ_c,     A=(3t⁴+2t³-9t²-11t+6)/(3(t-2)(5t²+4)),   B=-(t²-2)/(3(t-2)(5t²+4)),

also independent of φ,ω **[exact]** (cell 16: the spectral projectors of M₀ inside the octet subspace of the (2,2,2) sector). The octet pair enters only in the symmetric, rational-in-t combination. Re-derived here: 2A+B(τ₁+τ₂)=-t/(3(t-2)) with τ₁+τ₂=3(2t-1)(t+2) (the numerator 2(3t⁴+2t³-9t²-11t+6)-3(t²-2)(2t²+3t-2)=-t(5t²+4) cancels against the denominator). At P1 the two shares are 0.7307 and 1.8527 (using τ=7.7293,21.0413), summing to 31/12.


### 4.3 Octets (Claim 2)

Sector (3,2,1) (dim 6, two HW kets). **The τ_c.** The two octet HW kets are the eigenvectors of t₂(θ₁+h) on HW_((3,2,1)); restricting the operator to the HW kets, t₂(θ₁+h)H=HA(t) with (cell 12, symbolic in t)

    A(t)=[[3/2(2t²+5t-1)3/2(t-1)], [3/2(t-3)3/2(t-1)(2t+3)]],     tr A=3(2t-1)(t+2),   det A=9t(t-1)(t+1)(t+3),

**[exact]**. The discriminant is (trA)²-4det A=9(5t²+4) (re-derived here), a conic in (t,τ), rationally parametrized by t=4m/(5-m²), for which √(5t²+4)=(10+2m²)/(5-m²)=2+mt (re-derived here). So both τ_± are rational functions of the single symbol m, and **every octet computation is exact and symbolic**, with no `RootReduce` and no sampling: an identity in m holds on the whole curve τ²-(trA)τ+det A=0, i.e. for both copies at all t. The deck involution is m↦-5/m (it fixes t and flips the sign of the square root; re-derived here).

**The top-killed overlap.** With Sig₀=⟨ V₁⁰⊗ v₂⁰|h_c⟩=a₀'c₀' (u+u')· C_c (u=(-1,3), u'=(2,0) on the (1,1,1) indices),

    g_c=(F₁)/(F₀)=(⟨∂_ε(V₁⊗ v₂-v₂⊗ v₂)|h_c⟩+⟨ V₁⁰⊗ v₂⁰-v₂⁰⊗ v₂⁰|ψ₁⟩)/(Sig₀)

(the factors δ^Λ, a'₀c'₀ and the scalar normalisation series cancel up to an ω₂-independent constant; cells 11, 13). Results **[exact]**:

- g_c is affine in ω₂ (second derivative 0); the slope is independent of ω₁ and of φ;
- the slope is -(t-2)R_c+Δ_c with R_c=(u· C_c)/((u+u')· C_c); the **signal part** contributes exactly -(t-2)R_c (remainder 0), and the **contamination part** ∂_(ω₂)⟨ V₁⁰⊗ v₂⁰-v₂⁰⊗ v₂⁰|ψ₁⟩/Sig₀ equals

    Δ_c=(1-t²)+(t τ_c)/3(t+3)

exactly, for both copies, symbolic in m, and for the second family member as well.

**Which irreps carry the slope** (cell 15): ψ₁=Σ_νΠ_νψ₁ with ν∈{27, 10, 10̅, other octet}. The two decuplets contribute Δ_c⁽¹⁰⁾=Δ_c^((10̅))=1/2Δ_c=-1/2(t²-1)+tτ_c/(6(t+3)) each; the 27 and the other octet contribute ω₂-independent, φ-dependent (intercept-only) terms with zero slope; the sum over ν equals Δ_c symbolically. In the notation of §3.3 the intercept of g_c is where the φ-dependence of ρ_B sits.


### 4.4 The clean branches (bonus)

For [4,1,1] and [3,3,0] (sectors of dimension 2: the decuplet HW ket and a 27 descendant), cell 18, exact in (t,φ,ω) **[exact]**:

    
    [4,1,1]:   S₀=((ω₁-φ₂)(φ₁-φ₃)(φ₂-φ₃))/((φ₁-φ₂)(ω₁-φ₃)(φ₃-ω₂)),   (K₁)/(S₀)=((t-1)(φ₂+φ₃-2ω₂))/4,   (K₁ᵗᵏ)/(S₀)=((t-1)(φ₂-φ₃))/4,
    
    [3,3,0]:   S₀=((ω₁-φ₁)(φ₁-φ₂))/((ω₁-φ₂)(φ₂-ω₂)),   (K₁)/(S₀)=((t-1)(φ₁+φ₂-2ω₂))/4,   (K₁ᵗᵏ)/(S₀)=((t-1)(φ₁-φ₂))/4.
    

The top-killed contamination is ω-independent. The 27 ([4,2,0]) has a 1-dimensional sector, so ψ₁=0. These are not needed for the CG extraction. **Cross-check against the numerics (re-derived here):** K/S=ε K₁/S₀, and at φ, ω A the formulas give K₁/S₀=0.10970 ([4,1,1]) and 0.08633 ([3,3,0]), i.e. κ-1=1.1·10⁻⁷ and 8.6·10⁻⁸ at ε=10⁻⁶, exactly the values measured in `L2-210-untwist-CG.wb` cell 25 (§3.1).


### 4.5 Cross-checks of the perturbation theory itself

Cells 19–21, 25. **Independent finite-ε evaluation** **[numeric]** at the rational point m=2 (t=8, where τ_±=252,198 are rational), at ε=10⁻⁸ and 80 digits, with no Taylor pieces and no symbolic m: the exactly twisted sector operator M(ε)=Σ_Iz_I[m_I(3/7)+5/3m_I(1/5)], its true eigenvector near h_c, and the exact finite-ε covectors. The true eigenvector reproduces ψ₁ to 2.6·10⁻⁹ and 4.6·10⁻⁹ relative (O(ε)); the numerical slopes of g in ω₂ give Δ_c=-1.90909088 against the claim -1.90909091 (τ=252), and -15.0000024 against -15 (τ=198). The two agree to 3·10⁻⁸ (τ=252) and 2.4·10⁻⁶ (τ=198) in absolute terms; the write-up and HANDOFF quote “3·10⁻⁸” for the pair, which is true only of the first (see Appendix C). At P1 (t=62/35) the closed form gives Δ=-1.18143552, 0.46595220 against the earlier fitted -1.1814355, 0.4659522 (re-derived here from the closed form and the τ's of cell 45: -1.18143552, 0.46595220). Independence of the family member: cell 21.

What the derivation does **not** give: it does not prove the [3,1,0] constants (second-order perturbation theory is needed for the ε² layer), and it proves Δ_c and κ_(singlet) in the regime where the eigenvector is the first-order perturbation of a generic family member (non-degenerate within the sector after the irrep decomposition; the two octets are split by the inhomogeneity through the τ_c).


---

## 5. [3,1,0]^(⊗2), point A, through the same pipeline

Notebook: `L2-310-untwist-CG.wb` (90 cells; no native outputs: VS Code never opened the file, the sonnet agent ran each cell's source in K1 via `Get` and stored the printed result as a markdown cell directly under each code cell; summary cells 18, 24, 32, 38, 50, 51 for round 1 and 64, 84, 89, 90 for round 2). Spec `L2-310-untwist-CG-spec.md` (round 1); round 2 (closed forms, ε² layer, final extraction) was ordered by the main session during the day. Representation Λ=(3,1,0) on both sites (dim 15, 225 states, 36 distinct total weights), 10 HW states in 8 irreps. The pipeline is that of §3; the point-A data are those of `Undo-site2-rotation-notes.md` §5.2–5.4, §6.


### 5.1 Regression at finite twist

Cells 15–18. Before any small-ε work, point A was rebuilt from the closed-form dressing at z=(3/2,5/4,2/3), w=(5/2,7/3,3/5) (y_A=-63/8): (i) ⟨ v₁⊗ v₂|HW⟩=s⟨φ_A⊗ e₁ᵀ|HW⟩ exactly on the 10 HW kets (and, at three regime points, also (ρ(n)⊗ρ(n))HW=HW, cell 10) **[exact]**; (ii) ⟨ v₁⊗ v₂|Ψ⟩/⟨ x₀|Ψ⟩=const×det₆/Ψ_Q×det₂^((0,1,2))[y_A]/W(3,4) over **all 225 states** (36 weight sectors, 200 digits): max|ratio/ratio₁-1|=1.3·10⁻¹⁷⁹ (the constant has precision 187); the controls fail (wrong y=2: 10.3, y=3: 44.1, and the pole set {0..3} at y_A: 44.9) **[numeric]**. Dual-Q Baxter residual 3.8·10⁻¹⁸⁴ (225 states), vacuum 5.3·10⁻¹⁹⁵ (degrees 0,4,6). Timings: 225 states 13.4 s, det₆×det₂ on 225 states 3.7 s, one R_Q 0.014 s. The closed forms of the spec are right as given: K₁₂=β(w₂-z₁)/w₂=2/25 and K₁₃=γ(w₂-z₁)(y_A-z₃)/w₂²=-189/704; `Undo-site2-rotation-notes.md` §5.3 (md line 247, tex line 265) quotes +189/704, a **sign typo** in the notes (the formula passes the test; cell 17 prints -189/704) (Appendix C).


### 5.2 The ε-structure of the undone covector

Task 0 (cells 19–24; exact, symbolic ε, rational φ and ω at φ, ω A). The per-weight orders and ranks are in the table of §2.4 and the mechanism (the E₁₂² cancellation, γ=O(ε), b'/(a'c')) in §2.5; the control that the factorized form φ_A/δ^Λ=[e₁ᵀg diag(Πδₖ^(nₖ))]ρ(n') equals φ_A/δ^Λ holds exactly at ε=1/1000 (15 components) **[exact]**. Here δ^Λ=δ₁³δ₂ and φ_A(0)/δ^Λ=e₁ᵀ(1+(a'₀)/2E₁₂)e^(c'₀E₂₃) with b'₀=a'₀c'₀ (b₀'=a₀'c₀' and β̃₀=-a'₀/2 verified exactly, cell 21). The leading coefficients of the twelve weight components of φ_A/δ^Λ (cell 21, exact) are, at φ, ω A **[exact]**:

| weight | order | leading coefficient(s) |
|---|---|---|
| (3,1,0) | ε⁰ | 1 |
| (3,0,1) | ε⁰ | 2755/3549 (=c'₀) |
| (2,2,0) | ε⁰ | -2275/783 (=a'₀) |
| (2,0,2) | ε⁰ | -26172500/3737097 |
| (2,1,1) | ε⁰ | (-2375/1053, -9500/1053)∝(1,4)∝ u'=(2,8): one direction |
| (1,3,0) | ε¹ | -295750/204363 |
| (1,2,1) | ε¹ | (-204250/147987, -396625/147987) |
| (1,1,2) | ε¹ | (-38807500/18110547, -451250/223587) |
| (1,0,3) | ε² | -3978220000/969845877 |
| (0,3,1) | ε² | -2470000/12874869 |
| (0,2,2) | ε² | -72200000/363604059 |
| (0,1,3) | ε³ | -1097440000/4700022327 |

(γ=ε·(-2565/56056)+… at the same point; the identification of (3,0,1) and (2,2,0) with c'₀, a'₀ is re-checked here against the values of §2.3: a'₀=-2275/783, c'₀=2755/3549.) The ranks over the 3×6 grid are identical for φ and φ', and the only weight with a component beyond ε² is (0,1,3) **[exact]**.


### 5.3 Chain side: clean and contaminated branches

Task 1 (cells 25–32): exact undo to 4.9·10⁻¹³⁸ at 540 points (3 ε × 18 ω × 10 branches), minimum Precision of K 137.6; the sum of carrier contributions equals K to 3.3·10⁻¹⁴⁰. Raising residual of the branches ∼10⁻⁷ versus at least 1.41 for the other states. ε-exponents, κ and classification **[numeric]**:

| branch | slice weight | S | K | ‖Ψ‖ | κ at ω A | class |
|---|---|---|---|---|---|---|
| [6,2,0] | (3,1,0) | -6 | K≡0 | -2 | 1 | exact |
| [6,1,1] | (3,0,1) | -7 | -6 | -3 | 1+O(10⁻⁶); max\| K/S\|=9.6·10⁻⁷ at 10⁻⁶ | clean |
| [5,3,0] | (2,2,0) | -7 | -6 | -3 | 1+O(10⁻⁶); 2.6·10⁻⁷ | clean |
| [5,2,1] #1, #2 | (2,1,1) | -8 | -7 | -4 | 1+O(10⁻⁶); 8.3·10⁻⁷, 1.8·10⁻⁶ | clean |
| [4,4,0] | (1,3,0) | -7 | -7 | -4 | 1.3857144 (=97/70), constant | contaminated |
| [4,3,1] #1, #2 | (1,2,1) | -8 | -8 | -5 | 1.4105702; 0.98671310 | contaminated |
| [4,2,2] | (1,1,2) | -9 | -9 | -6 | 1.3661018 (=1+108/295) | contaminated |
| [3,3,2] | (0,2,2) | -9 | -9 | -7 | 4.6508004 | contaminated |

The clean set is exactly the set of branches whose slice weight is at ε⁰ in the table of §2.4; the contaminated set is the branches whose slice weight is at ε¹ or ε², where K (which comes from ε⁰ weights next door) has the same order as S. The κ of the contaminated branches are constant in ε (e.g. 1.3857144,1.3857144,1.3857146 at ε=1,2,4·10⁻⁶) with spreads ≤10⁻⁵ over ε, ω, φ (cell 50). Contamination carriers (ε=10⁻⁶,2·10⁻⁶, all c_ν scaling like S): [4,4,0] from (2,2,0) (100%); [4,3,1] #1 from (2,1,1) (100%); #2 from (2,1,1) (≈100%, plus (1,2,1): 1.6·10⁻⁵, (1,3,0): 4.1·10⁻⁵, (2,2,0): 4.3·10⁻⁶ of K); [4,2,2] from (2,0,2) (≈100%); [3,3,2] from (1,1,2) (85%) and (2,0,2) (15%) **[numeric]**.


### 5.4 Q side and the leading-order extraction

Task 2 (cells 33–38) **[numeric]**. Vacuum direct Q's (degrees 0,4,6) of the 54 twists (w₁,w₂,w₃): Baxter residual ≤2.6·10⁻¹⁸⁹; dual Q's of the 56 HW-sector states: 2.4·10⁻¹⁴⁰ (exact arithmetic). R_C/R_Q is state-independent over the 10 branches to 2.5·10⁻⁸², 3.7·10⁻⁸⁶, 8.9·10⁻⁹³ (ε=1,2,4·10⁻⁶, max over 18 ω; limited by Precision(R_Q)≥103.7 digits, not by the identity), and over all 56 HW-sector states at ω A to 6.4·10⁻⁷⁷, 2.0·10⁻⁸⁵, 2.8·10⁻⁸⁷. Leading-order extraction for the five clean branches, C^Q=Qₙᵤₘ/(δ^Λ fac) with fac the component of φ_A(0)/δ^Λ at the slice weight: 1, c'₀, a'₀ (the GT matrix element of e₁ᵀE₁₂ is 2, so the coefficient a'₀/2 of E₁₂ gives a'₀), and (a'₀c'₀/2)u' for (2,1,1) with u'=(2,8):

| branch | spread / error at ε=10⁻⁶ | at 4·10⁻⁶ | Richardson limit (Q = chain) | \| Q/chain-1\| |
|---|---|---|---|---|
| [6,2,0] | 0 / 0 | 0 | 2.609467455621302 | 0 |
| [6,1,1] | 8.5·10⁻⁷ / 7.9·10⁻⁷ | 3.4·10⁻⁶ / 3.2·10⁻⁶ | 46.67453891138102 | 1.2·10⁻¹⁹ |
| [5,3,0] | 2.0·10⁻⁶ / 1.9·10⁻⁶ | 7.8·10⁻⁶ / 7.4·10⁻⁶ | -36.88381884387802 | 5.4·10⁻²⁰ |
| [5,2,1] #1 (u'· C) | 4.2·10⁻⁶ / 4.0·10⁻⁶ | 1.7·10⁻⁵ / 1.6·10⁻⁵ | -1325.402366533307 | 1.1·10⁻¹⁹ |
| [5,2,1] #2 (u'· C) | 2.1·10⁻⁶ / 2.5·10⁻⁶ | 8.6·10⁻⁶ / 9.8·10⁻⁶ | -334.3316171624337 | 1.2·10⁻¹⁹ |

(Richardson of Cεᵖ, p=2,3,3,4,4, three points; scaling ×4.00.) The top value 2.6094674556 and the [6,1,1] value 46.6745389114 coincide with the adjoint's 27 and [4,1,1] values of §3.4. For the [5,2,1] pair only the functional u'· C is read at leading order; C itself needs the ε¹ layer.


### 5.5 The O(ε) layer

**Top-killing** (cells 40–41). T=δ^Λ s⟨ v₂⊗ v₂|Ψ_B⟩/⟨ x₀|Ψ_B⟩ equals the (w₁,y_A)-vacuum det₆/Ψ_Q (no det₂), top-normalised, to 2.5·10⁻⁸² over the 10 branches (vacuum Baxter residual 9.5·10⁻¹⁹⁴) **[numeric]**. T∼ε⁻⁶ for every branch (pure contamination: ⟨ v₂⊗ v₂|HW⟩=0 for non-top HW kets); T/K=2.30137 for [6,1,1] and -1.24138 for [5,3,0] are O(1), and ≤4·10⁻⁸ for all other branches, so the top-kill matters at leading order only for the two clean decuplet-type branches [6,1,1], [5,3,0] and at relative O(ε) for the [5,2,1] pair.

**The two-direction weight (2,1,1) and the [5,2,1] pair** (cells 42–47, 54–64). Convention: u=e₁ᵀE₁₃=(-1,8) and u'=e₁ᵀE₁₂E₂₃=(2,8) on the two GT indices (the order E₁₂E₂₃ is the one that appears in e₁ᵀ(1+(a'₀)/2E₁₂)e^(c'₀E₂₃); the adjoint's u' is e₁ᵀE₂₃E₁₂), φ_((2,1,1))/δ^Λ=Aᵤu+A_(u')u', R_c:=(u· C_c)/(u'· C_c), C_c=[P_(HW)Ψ_c]_(((2,1,1),e₁)) Richardson-extrapolated. Exact (rational, at φ, ω₁∈{3/4,-5/3,7/2}): Aᵤ(0)=0, A_(u')(0)=fac₀=a'₀c'₀/2, and the first-order coefficients are exact straight lines in ω₂, ω₁-independent: A⁽¹⁾_(u')/fac₀=-113/140ω₂+const and A⁽¹⁾ᵤ/fac₀=59/70ω₂+const **[exact]**. With F=ε⁴Qₙᵤₘᵗᵏ/(δ^Λfac₀)=F₀+ε F₁+… and g=F₁/F₀ **[numeric]**:

| quantity | copy 1 | copy 2 |
|---|---|---|
| F₀ spread over 12 ω | 2.2·10⁻¹⁷ | 4.7·10⁻¹⁷ |
| g(ω₂) a straight line: residual | ≤2.5·10⁻¹¹ | ≤2.5·10⁻¹¹ |
| ω₁-dependence of the slope | 0 (10⁻¹⁶⁶) | 0 (10⁻¹⁶⁷) |
| slope (8 φ; spread) | -1.90107065557 (1.0·10⁻¹¹) | 0.97351289979 (6·10⁻¹²) |
| R_c=(u· C)/(u'· C), φ-independent | -0.8921886416 | 1.2040698297 |
| signal slope =-113/140+59/70R_c | -1.5591304265 | 0.2077159993 |
| Δ_c:=slope-slope_(sig) (spread over 8 φ) | -0.34194022909 (6·10⁻¹²) | 0.76579690047 (6·10⁻¹²) |

F₀ itself (=u'· C) and the intercept depend on φ (-1325.40, -334.33 at φ; -280.13, -70.66 at φ'); slope, R_c and Δ_c are structural. Δ₁+Δ₂=0.42385667137. Note that the exact signal slopes are p'₁=(t-5)/4 for the u'-part and q'₁=(7-3t)/2 for the u-part (times R_c), found as exact straight lines at all 20 t (-113/140 and 59/70 at P1), and that the signal's coefficient of R_c is **not** -(t-2) as for the adjoint; so the combination of the spec, slope+(t-2)R, is not the signal-free constant here (see Appendix C).

**The branches with F₀≡0** (leading covector order ε¹ or ε²; cells 48–50) **[numeric]**. K/S is a constant over ε, ω, φ (3 ε × 18 ω; 8 φ at ω A; the variations are O(ε)):

| branch | slice weight, order | K/S (ε→0) | spread over ε,ω | spread over 8 φ | R_ν=C[idx2]/C[idx1] (φ / φ') |
|---|---|---|---|---|---|
| [4,4,0] | (1,3,0), ε¹ | 0.38571428571=27/70 (κ=97/70) | 8.5·10⁻⁶ | 3.1·10⁻⁷ | - |
| [4,3,1] #1 | (1,2,1), ε¹ | 0.41056992897 | 2.2·10⁻⁵ | 2.7·10⁻⁷ | -0.0789657625 / same |
| [4,3,1] #2 | (1,2,1), ε¹ | -0.01328651383 | 7.3·10⁻⁴ | 2.2·10⁻⁵ | -8.7671880836 / same |
| [4,2,2] | (1,1,2), ε¹ | 0.36610169492=108/295 | 1.8·10⁻⁵ | 2.0·10⁻⁶ | -1/3 exactly / -1/3 |
| [3,3,2] | (0,2,2), ε² | 3.65080357143 | 1.8·10⁻⁵ | 6.6·10⁻⁷ | - |

(R_ν in the GT-index convention u= unit covector on the first index, u' on the second; φ-independent to 10⁻²⁰.) [4,4,0] at P1 is rational, 27/70=θ₁₂-h in absolute units. The two [4,3,1] copies have copy-dependent constants, as the adjoint octets (where τ_c entered). The clean branches [6,1,1], [5,3,0] also have an O(ε) layer, K/S=λε with λ ω- and φ-dependent (at ω A, φ: λ=0.10970, -0.00789); it is not needed for the extraction and was not pursued. Observation (re-derived here): λ=0.1097031447 for [6,1,1] equals, to 8·10⁻⁸ relative, the adjoint [4,1,1] value K₁/S₀=(t-1)(φ₂+φ₃-2ω₂)/4=5913/53900=0.1097031540 of §4.4; no explanation is given in the sources **[numeric]**.

**Round 2, Task A: closed forms by a t-scan** (cells 53–65; h=1, θ=(t,0); 20 values of t: 29/15, 5/2, 3/2, 7/3, 11/4, 13/5, 9/4, 17/6, 19/7, 8/3, 12/5, 7/2, 23/8, 16/5, 11/3, 14/5, 25/9, 13/4, 9/5, 21/8; fit on 14, held out the 3rd, 6th, ..., 18th; rational fits by the null vector of P(t)-fQ(t)=0 at 40 digits; φ, ω₁=3/4, six ω₂, three-ε Richardson; 30 s per t). At t=62/35 (h=1) the round-1 constants are reproduced. Exact untwisted data of the chain (rational, all 20 t, held-out exact): t₂(θ₁+h) on HW_([4,3,1]): tr=8t²+16t-8, det=16t(t-1)(t+1)(t+4); on HW_([5,2,1]): tr=8t²+19t-8, det=4(t+4)(t-1)(t+1)(4t+3) **[exact]**; the copy eigenvalues of the branches' own t₂ polynomials agree with these to 9·10⁻¹⁹; t₁(θ₁+h) is degenerate on the [4,3,1] pair (difference 4·10⁻¹⁸), so the affine-in-t₁ alternative cannot separate the copies. Closed forms (τ_c=eig(t₂(θ₁+h))/h⁴, D=t³+3t²-17t+5):

| constant | closed form | residual (fit / held-out) |
|---|---|---|
| κ_([4,4,0]) | (t+1)/2 | 1.4·10⁻¹⁸ / 2.6·10⁻¹⁸ |
| κ_([4,2,2]) | (11t-31)/(5(3t-7)); pole at t=7/3 (chain value -2.5·10⁶ there, ε-regulated; point skipped) | 3.4·10⁻¹⁵ / 1.9·10⁻¹⁷ |
| κ_([3,3,2]) | (t+1)(19t-49)/(40(t-2)) | 7.2·10⁻¹⁷ / 4.9·10⁻¹⁷ |
| κ_c^([4,3,1])=A+Bτ_c | A=(-7t⁴+26t³-166t+75)/(9D), B=(7t²-22t+23)/(36D) | 3.0·10⁻¹⁰ (data amplified near the root t≈2.65 of D) / 9.1·10⁻¹⁶ |
| R_c^([5,2,1]) | 3/2-2t+τ_c/(2(t+4)); Rₐ+R_b=(16-7t)/(2(t+4)), RₐR_b=-7t/(2(t+4)) | 9.7·10⁻¹⁹ / 2.3·10⁻¹⁸ |
| signal slopes [5,2,1] | p'₁=(t-5)/4, q'₁=(7-3t)/2 | exact (all 20 t) |
| Δ_c^([5,2,1])=slope-p'₁-q'₁R_c | (111+44t-92t²)/105+((23t-13) τ_c)/105(t+4) | 7.1·10⁻¹⁰ / 1.6·10⁻⁹ (noise-limited) |

**[fitted]** (κ's and R_c) and **[fitted]** at the 10⁻⁹ level (Δ_c^([5,2,1]): the coefficients with the common 1/105 are rationalized from a 10⁻⁸-accurate least-squares fit and are not exact-proven; Δ₁+Δ₂=(-3t²+143t/105+992/105)/(t+4) and the product are rational fits of degrees (2,1), (3,1) with held-out 5.6·10⁻¹¹, 1.5·10⁻⁹; re-derived here: the sum follows from the quoted Δ_c and tr=8t²+19t-8). The worst line residual (1.3·10⁻⁹) is at t=29/15, close to the singularity of y₀ at t=2. The sources note that the first automatic fit of q₁'(t) reported “NO FIT” only because q₁' vanishes at t=7/3 (relative error 0/0); the linear form is exact. The cross-check against round 1 at t=62/35: the closed form gives Δ=0.765797 (τ=37.479) and -0.341940 (τ=13.282) (re-derived here); round 1 measured 0.76579690 and -0.34194023. (Cell 64 prints “0.7655” and cell 61 prints a “sanity” line with other numbers; see Appendix C.)

Re-derived here from the quoted trace and determinant: the discriminants of the two 2×2 problems are (tr)²-4det=64(3t²+1) on HW_([4,3,1]) and 105t²+256 on HW_([5,2,1]) (the adjoint octets had 9(5t²+4)). Each is a conic in (t,τ), so the rational-parameter trick of §4.3 would apply to the first-order perturbation theory of [3,1,0] as well **[conjecture]** (not attempted). At P1 these give τ₂-τ₁=25.82 and 24.20, matching the printed τ's (35.6315-9.8150 and 37.4788-13.2820).


### 5.6 The O(ε²) layer

Round 2, Task B (cells 66–84): the weights (1,2,1) ([4,3,1]×2) and (1,1,2) ([4,2,2]), whose first direction is at ε¹ and whose second direction enters at ε². **Exact covector structure** **[exact]** (all rational): φ_ν/δ^Λ=ε c₁+ε²c₂+… with c₁=fac₁(ω,φ) u', u'=(1,x_U) (normalised so its first GT component is 1), u=(0,1) and c₂=a₂u'+b₂u. Convention R:=(u· C)/(u'· C)=C₂/(C₁+x_UC₂). a₂/fac₁ and b₂/fac₁ are exactly linear in ω₂ (exact interpolation of degree 1 through 10 points, both weights) and ω₁-independent:

| weight | x_U(t) | ω₂-slope of a₂/fac₁ | ω₂-slope of b₂/fac₁ | at P1 |
|---|---|---|---|---|
| (1,2,1) | -(t+3)/(2(t-3)) | t(2-t)/(t-3) | not fitted | x_U=167/86, a-slope -496/1505, b-slope 525/3698 |
| (1,1,2) | -3(t-1)/(2(t-3)) | -19/35 at P1 (zero at an isolated t; t-form not fitted) | 3(t-1)/(t-3) | x_U=81/86, a-slope -19/35, b-slope -81/43 |

(The t-forms of x_U at P1 are re-derived here: 167/86 and 81/86.) **The Q-number structure at P1** (four-ε Richardson with ε=1,2,4,8·10⁻⁶; ten ω₂, ω₁=3/4 unless stated; F=εᵖ Qᵗᵏₙᵤₘ/(δ^Λfac₁), p=4,5) **[numeric]**: F₀ is ω-independent (spread 3·10⁻²¹, 1.3·10⁻²¹, 4·10⁻²¹); g=F₁/F₀ is a straight line in ω₂ (degree-1 residual ≤4.3·10⁻¹⁵; degrees 2–4 do not improve beyond the 10⁻¹⁵ truncation level); its slope is ω₁-independent (ω₁=7/2: identical to 10⁻¹⁵) and φ-independent (φ', φ₃, φ₄: identical to ≤5·10⁻¹⁶); the intercept and F₀ depend on φ. Write g=s₁+rem with the signal part s₁=a₂/fac₁+(b₂/fac₁)R (the R-carrying coefficient is the b-slope):

| branch | R | s₁ slope | g slope | remainder slope | κ₀ | kslope=κ₀×rem |
|---|---|---|---|---|---|---|
| [4,3,1] #1 (τ=9.815) | -0.0932674373 | -0.3428091568 | -1.0162824575 | -0.6734733008 | 1.4105699290 | -0.9499811860 |
| [4,3,1] #2 (τ=35.631) | 0.5471061670 | -0.2518961924 | 0.7940852803 | 1.0459814726 | 0.9867134862 | 1.0320840253 |
| [4,2,2] | -0.4858757062 | 0.3723970944 | -0.1527827012 | -0.5251797956 | 1.3661016949 | -0.71743 |

The natural structural constant is kslope, the ω₂-slope of the O(ε) correction κ₁ of R_C/S=κ₀+εκ₁+… (F₁/F₀=s₁+κ₁/κ₀). **t-scan** (20 t, 14 fit, 6 held out, after the memo fix of §7; [4,2,2] excludes t=7/3) **[fitted]**:

- R_c^([4,3,1])=A_R+B_Rτ_c, A_R=(-10+28t/3+4t²-4t³+2t⁴/3)/D, B_R=(-3/2+t-(t²)/6)/D (held-out 1.9·10⁻¹⁶, 5.7·10⁻¹⁸; Rₐ+R_b and RₐR_b rational (3,3) with the same D: 8.8·10⁻¹⁸, 1.8·10⁻¹⁷);
- [4,2,2]: R=-2(t-3)/(3(3t-7)) (exact in the fit), kslope=(-28/45t³+2912/675t²-448/45t+4948/675)/(t-7/3)² (held-out 1.5·10⁻¹¹; the coefficients are 10⁻⁹-noise-limited rationals; re-derived here: at P1 this evaluates to -0.717449 against the measured -0.71743);
- [4,3,1]: kslope_c=α(t)+β(t)τ_c with β=(3331-7425t+7876t²-2262t³-481t⁴+247t⁵-6t⁶)/(216D²) (held-out 1.7·10⁻¹¹) and α=P₈(t)/D² with P₈ of degree 8 (leading coefficient ≈1/9), only a **numerical** fit (held-out 10⁻⁹ relative; the coefficients are not rationalizable at the 10⁻⁶ noise level and live in the notebook variable `caC`). At P1 the closed form reproduces the measured kslope: -0.94998118599 (copy 1) against -0.949981186023.

What did **not** fit: (i) the remainder slope rem_c=kslope_c/κ_c as an affine-in-τ_c form with polynomial/Dᵐ coefficients (degrees up to 6–8, fit-set residual ≥10⁻¹; one has to divide by κ_c=A+Bτ_c, which is not affine); (ii) symmetric functions of rem_c and of kslope_c in low degree; (iii) α(t) is not an identified rational; (iv) the t-forms of the (1,2,1) b-slope and the (1,1,2) a-slope.


### 5.7 Final pure-Q extraction of all ten branches

Round 2, Task C (cells 85–89). Inputs from the Q side only: Qₙᵤₘᵗᵏ=Qₙᵤₘ-T_q (T_q the top-kill det₆ with the (w₁,y_A)-vacuum) at 6 ω₂ × 4 ε; τ_c=t₂(θ₁+h)/h⁴ from the dual-Q dictionary t₂(θ₁+h)=χ₃ν₁(θ₁+h)ν₂(θ₁)Qʲ(θ₁-h)/Qʲ(θ₁) (single-ratio form, Richardson in ε; agrees with the exact untwisted eigenvalues of Task A to 10⁻¹⁸); the closed forms of Tasks A–B; the exact covector functions (facₖ, x_U, the a- and b-slopes are analytic in the twists). Extraction: first direction u'· C=Qᵗᵏₙᵤₘ/(δ^Λεᵏfacₖ κ) (k=0,1,2 the leading covector order); second direction (dim-2 weights): R^Q from the ω₂-slope of F₁/F₀ minus the known signal and the closed-form remainder; the full slice from {u· C=Rc₁, u'· C=c₁}. The chain side enters only the comparison columns. Columns: error =max_(6 ω) at ε=10⁻⁶ of the first-direction extraction against the chain; Richardson =|Q/chain-1| of the ε→0 limits at ω A; dR=|R^Q-Rᶜʰᵃⁱⁿ|; slice = max componentwise relative error of the ε→0 slice **[numeric]**.

| branch (P1) | extracted | error | Richardson | τ_Q | R^Q | dR | slice |
|---|---|---|---|---|---|---|---|
| [6,2,0] | full | 0 | 0 | - | - | - | - |
| [6,1,1] | leading | 3.1·10⁻⁷ | 9.9·10⁻²⁰ |  |  |  |  |
| [5,3,0] | leading | 2.1·10⁻⁶ | 5.4·10⁻²⁰ |  |  |  |  |
| [5,2,1] #1 | u'C+R | 4.0·10⁻⁶ | 1.1·10⁻¹⁹ | 13.28200213 | -0.892188641588 | 6.2·10⁻¹⁷ | 1.6·10⁻¹⁶ |
| [5,2,1] #2 | u'C+R | 2.5·10⁻⁶ | 1.5·10⁻¹⁹ | 37.47881420 | 1.20406982971 | 2.6·10⁻¹⁶ | 1.3·10⁻¹⁵ |
| [4,4,0] | leading, κ=(t+1)/2 | 9.5·10⁻⁷ | 6.3·10⁻¹⁹ |  |  |  |  |
| [4,3,1] #1 | u'C+R | 2.5·10⁻⁶ | 9.6·10⁻²⁰ | 9.815049286 | -0.0932674374557 | 1.7·10⁻¹⁰ | 1.8·10⁻⁹ |
| [4,3,1] #2 | u'C+R | 1.7·10⁻⁶ | 1.5·10⁻¹⁹ | 35.63148133 | 0.547106166786 | 2.4·10⁻¹⁰ | 7.4·10⁻⁹ |
| [4,2,2] | u'C+R | 5.6·10⁻⁷ | 1.2·10⁻¹⁸ |  | -0.485875706215 | 1.9·10⁻¹⁶ | 3.8·10⁻¹⁶ |
| [3,3,2] | leading, ε² covector | 2.7·10⁻⁷ | 1.9·10⁻¹⁹ |  |  |  |  |

| branch (P3, out of sample) | extracted | error | Richardson | τ_Q | R^Q | dR | slice |
|---|---|---|---|---|---|---|---|
| [6,2,0] | full | 0 | 0 | - | - | - | - |
| [6,1,1] | leading | 5.4·10⁻⁷ | 2.8·10⁻¹⁹ |  |  |  |  |
| [5,3,0] | leading | 2.1·10⁻⁶ | 4.2·10⁻²⁰ |  |  |  |  |
| [5,2,1] #1 | u'C+R | 3.8·10⁻⁶ | 3.7·10⁻¹⁹ | 5.660760026 | -0.685948111847 | 4.7·10⁻¹⁷ | 2.5·10⁻¹⁶ |
| [5,2,1] #2 | u'C+R | 4.7·10⁻⁶ | 3.1·10⁻¹⁹ | 26.85964814 | 1.29261477851 | 3.1·10⁻¹⁷ | 1.0·10⁻¹⁶ |
| [4,4,0] | leading | 6.4·10⁻⁷ | 4.4·10⁻¹⁹ |  |  |  |  |
| [4,3,1] #1 | u'C+R | 1.8·10⁻⁶ | 3.1·10⁻¹⁹ | 4.006458515 | -0.0492541950103 | 1.6·10⁻⁸ | 3.3·10⁻⁷ |
| [4,3,1] #2 | u'C+R | 4.7·10⁻⁶ | 3.6·10⁻¹⁹ | 24.44252108 | 0.865785999106 | 2.3·10⁻⁸ | 2.1·10⁻⁷ |
| [4,2,2] | u'C+R | 2.4·10⁻⁶ | 1.6·10⁻¹⁸ |  | -0.373983739837 | 1.8·10⁻¹⁷ | 5.0·10⁻¹⁷ |
| [3,3,2] | leading | 2.9·10⁻⁶ | 9.5·10⁻¹⁹ |  |  |  |  |

(P3 is h=1/2, θ=(3/7,-1/4), t=19/14.) “First direction only” is never the limit: every branch is extracted completely (one-component weights from the leading formula with κ(t); the dim-2 weights through R). The [4,3,1] second direction is limited to 10⁻⁸ (10⁻¹⁰ at P1) by the noise of the fitted α(t); the [4,2,2] and [5,2,1] ones by nothing but the arithmetic (10⁻¹⁶). The ε=10⁻⁶ errors are the expected O(ε) of the leading-order formulas, and the ε→0 Richardson limits agree to ≤1.6·10⁻¹⁸. The same caveat as for the adjoint holds: the closed forms are fits (only κ's and R_c are at the 10⁻¹⁶ level; Δ_c^([5,2,1]) and the [4,2,2] kslope at 10⁻⁹–10⁻¹¹; [4,3,1] α(t) at 10⁻⁹).


### 5.8 What remains for [3,1,0]

Identify α(t) (needs either higher accuracy, four-ε Richardson at smaller ε, or second-order perturbation theory as in §4); prove the fitted κ's, Δ_c and kslope by first- and second-order perturbation theory; the [5,2,1] remainder and the [4,2,2] kslope are rationalized from 10⁻⁹–10⁻¹¹ fits; [4,1,0] at k=2 (dressing of E₁₃-degree 2: four orders expected **[conjecture]**); poles of the structural constants at t=7/3 ([4,2,2]) and t=2 ([3,3,2]) must be excluded from fits.


---

## 6. Untwisting the Q-determinants themselves

Notebook: `L2-210-untwist-Qdets.wb` (45 cells; no native outputs, results in the markdown cells 15, 22, 40, 44, 45); spec `L2-210-untwist-Qdets-spec.md` (“Spec C1”); write-up `L2-210-untwist-Qdets.md`. Setting as in §3 (adjoint, P1, joint regime), cells 2–8 copied from `L2-210-untwist-CG.wb` cells 3–9. ε∈{1,2,4,8}·10⁻⁶; **48 direction points**: φ and φ' over the 18 (ω₁,ω₂) grid (36 points), plus φ₃=(1/2,1/5,-3/7) and φ₄=(2/3,-1/4,1/9) over six ω₂ at ω₁=3/4 (12 points); 200 digits; branch order 27, [4,1,1], [3,3,0], octet 1, octet 2, singlet. Goal (spec §0): the two pieces of the CG pipeline that so far enter as black boxes, D₂:=det₂⁽¹⁾[Q_B;y]/minor₁ and D₆:=det₆[Q_B,F_(vac)(w₁,w₂,w₃)]/Ψ_Q (and the top-kill partner D₆ₜ:=det₆[Q_B,F_(vac)(w₁,y,w₃)]/Ψ_Q), to leading and next-to-leading order in ε, in terms of **untwisted** objects: the untwisted Q-polynomials of the HW ket with staggered degrees, the untwisted transfer eigenvalues and the directions. The spec cites the only existing untwisting technology, `old-conventions/old/XXX/Experiments/SU3_L2_Untwisting.{md,tex}` §3.2–3.5 (old conventions, and for the norm, not for this overlap; I did not read it): bare rows at infinity, a Laplace split into a bare block that becomes a Vandermonde-type monomial in the twist directions and a periodic block that becomes an untwisted minor, and staggered untwisted degrees deg qʲ=μⱼ+(3-j), deg qᵢ=(4-μᵢ)+(i-1).


### 6.1 Exponents and direction monomials

Task 1 (cells 9–15) **[numeric]**. Exponents are log₂ ratios over the doubling ladder (all within 10⁻⁶ of an integer, except D₂ where the non-top exponent is 0 to 10⁻⁶); coefficients c₀=f/εᵖ by a cubic extrapolation in ε (error ∼10⁻²²); monomial test by least squares in the logs of the nine differences φ₁₂,φ₁₃,φ₂₃,ω₁-φₖ,ω₂-φₖ followed by a constant-spread test at all 48 points.

| factor | ε-exponent per branch |
|---|---|
| D₂=det₂⁽¹⁾/minor₁ | 1,0,0,0,0,0 |
| D₆=det₆/Ψ_Q | 10,10,10,9,9,8 |
| D₆ₜ (top-kill partner) | 6 for all |
| ⟨ x₀\|Ψ̂⟩ (=1/‖ P_(HW)Ψ‖_G; so ‖Ψ‖_G∼ε⁻ᵖ) | 2,3,3,4,4,6 |
| Qₙᵤₘ | -5,-6,-6,-7,-7,-8 |

Consistency: p(Qₙᵤₘ^B)=p(D₆^B)+p(D₂^B)-11-5 (top-state values D₂∼ε, D₆∼ε¹⁰, D₆ₜ∼ε⁶, topFac=⟨ v₁₂|e₁e₁⟩/⟨ x₀|e₁e₁⟩∼ε⁻⁵; re-checked here against the table).

**Leading coefficients** (constants exact where rational; spread ≤8·10⁻²² over 48 points unless stated):

- D₂ is independent of φ and ω except for the top state: top 893025/1024(ω₂-φ₂)ε (893025=945²); [4,1,1] and [3,3,0]: -1157625/512 (=-105³/2⁹); singlet -1157625/128; octets -7116.62856099…, -3487.22204765… (algebraic, the two copies) **[numeric]**.
- D₆ₜ is a pure monomial in ω₁ alone: D₆ₜ∼const·ε⁶Πₖ(ω₁-φₖ)^(4-μₖ), i.e. exponents (0,2,4) (top), (0,3,3) ([4,1,1]), (1,1,4) ([3,3,0]), (1,2,3) (octets), (2,2,2) (singlet); no φᵢⱼ and no ω₂ dependence. The constant is the same complex number 1.47784527·10⁶+1.41296475·10⁶ i for the three clean branches (it contains the e^(2π ix/h) insertion).
- **D₆ itself is not a monomial** (fit residual O(1), spread ∼10): the top branch carries a branch-independent, twist-dependent complex normalisation constant of the det₆ formula (the same value for 27, [4,1,1], [3,3,0]). The branch **ratios** D₆^B/D₆ᵗᵒᵖ are exact monomials, listed in the table below.
- ⟨ x₀|Ψ̂⟩ is a function of φᵢⱼ only (absolute value): top |φ₁₂|², [4,1,1] |φ₁₂φ₁₃φ₂₃|, [3,3,0] |φ₁₂|³, octets φ₁₂²|φ₁₃φ₂₃|, singlet |φ₁₂φ₁₃φ₂₃|²; constants 1, 0.272741187029, 0.272741187029, 0.0398430144518, 0.0773784594992, 0.00746241357134.
- **Qₙᵤₘ is an exact monomial for every branch** (constant spread ≤10⁻²²; integer deviation of fitted exponents ≤10⁻²²), see the table below.

| branch | D₆^B/D₆ᵗᵒᵖ (monomial; constant) | D₆ₜ^B/D₆ₜᵗᵒᵖ (constant) | Qₙᵤₘ^B (constant) |
|---|---|---|---|
| top | 1 | 1 | φ₁₂⁻²(ω₂-φ₁)⁻²(ω₂-φ₂)⁻¹ (1, =δ^ΛΨₜₒₚ(1,1)) |
| [4,1,1] | ((ω₁-φ₂)(ω₂-φ₂))/((ω₁-φ₃)(ω₂-φ₃)) (1; spread 10⁻¹²⁷) | (ω₁-φ₂)/(ω₁-φ₃) (1) | φ₁₂⁻²(ω₁-φ₂)/(ω₁-φ₃)(ω₂-φ₁)⁻²(ω₂-φ₂)⁻¹(ω₂-φ₃)⁻¹ (-70/27) |
| [3,3,0] | (ω₁-φ₁)/(ω₁-φ₂) (1) | (ω₁-φ₁)/(ω₁-φ₂) (1) | φ₁₂⁻²(ω₁-φ₁)/(ω₁-φ₂)(ω₂-φ₁)⁻²(ω₂-φ₂)⁻² (-70/27) |
| octets | (ω₁-φ₁)/((ω₁-φ₃)(ω₂-φ₃)) (-0.7767209220418339, 2.0514347617920940; algebraic; spread 2·10⁻²⁰) | (ω₁-φ₁)/(ω₁-φ₃) (1.17753621, 0.53110063) | φ₁₂⁻²(ω₁-φ₁)/(ω₁-φ₃)(ω₂-φ₁)⁻²(ω₂-φ₂)⁻²(ω₂-φ₃)⁻¹ (6.33834161515, -8.20301776029) |
| singlet | ((ω₁-φ₁)²)/((ω₁-φ₃)²(ω₂-φ₃)²) (-9800/2619=-3.7418862161130; spread 10⁻²¹) | ((ω₁-φ₁)²)/((ω₁-φ₃)²) (2107/2619=0.80450554) | φ₁₂⁻²((ω₁-φ₁)²)/((ω₁-φ₃)²)Πₖ(ω₂-φₖ)⁻² (38.8047459449) |

Reading: Qₙᵤₘ^B/Qₙᵤₘᵗᵒᵖ=[D₆^B/D₆ᵗᵒᵖ][D₂^B/D₂ᵗᵒᵖ] (the state-independent chain-side constants cancel); the top branch's D₂ carries the factor (ω₂-φ₂), which turns the ω₂-exponent of the D₆ ratio into the one listed. So the ω₁-dependence of Qₙᵤₘ is exactly that of D₆^B/D₆ᵗᵒᵖ. These exact monomials are the reason the leading-order extraction of §3.2 works with the simple factors c'₀, a'₀, a'₀c'₀.


### 6.2 det₂⁽¹⁾: the τ-polynomial and the untwisted Wronskian

Task 2 (cells 16–22). **(a) The τ-polynomial identity at small ε** **[numeric]**:

    D₂=c F(h,θ₁₂,y)·(1,τ₁₀,τ₁₁,τ₂₀,τ₂₁),
    
    F=(2h²(θ₁₂-h)(θ₁₂-2h)(θ₁₂-3h)y⁴, 2h(θ₁₂-h)(θ₁₂-3h)y³, (θ₁₂-h)y², (θ₁₂-3h)y², 2/3hy),

(θ₁₂ absolute here, as in `Two-rotations-det-overlap.md` §5.4; the dot product pairs the vector F with (1,τ_λ) componentwise-summed), τ_λ the eigenvalue of the fused T_λ at v₀=θ₁-2h read from the branch's own t₁,t₂,t₃ polynomials through the CBR (T₁₀=t₁, T₁₁=t₂, T₂₀=t₁t₁^([+])-t₂^([+]), T₂₁=t₂t₁^([+])-t₃^([+])). The constant c is the same number for all six branches, all 18 ω, all four ε, both φ and φ': max|c_B/cₜₒₚ-1|=0 at Precision 140 (the working precision of D₂), and

    c=343000/1161=(70³)/27·43=295.4349698535745…

(independent of ε, ω, φ, branch; the write-up adds “so the constant depends only on (h,θ)”, but only P1 was tested, so that part is untested). The convention check (no mismatch) passed.

**(b) Untwisted τ⁰ and first-order shifts.** HW kets of the untwisted chain (octet sector diagonalised in HW_((3,2,1)), copies ordered by the real part of the eigenvalue of t₁+π t₂ at u=3/7; eigen-residuals 10⁻¹⁹⁶). τ⁰ at P1 **[numeric]**:

| branch | T₁₀ | T₁₁ | T₂₀ | T₂₁ |
|---|---|---|---|---|
| 27 | -0.078571429 | -0.17770408 | 0.015153061 | 0.034271501 |
| [4,1,1] | -0.57857143 | -0.17770408 | 0.40086735 | 0.12312354 |
| [3,3,0] | -0.57857143 | -0.63841837 | 0.11158163 | 0.12312354 |
| octet 1 | -1.3285714 | -0.51334555 | 1.3150802 | 0.50813271 |
| octet 2 | -1.3285714 | 0.31865168 | 0.48308302 | -0.11586521 |
| singlet | -2.0785714 | 0.40086735 | 2.4794388 | -0.47817748 |

(27: T₁₀=-11/140.) The branch τ_λ(ε) obey τ=τ⁰+ετ¹+O(ε²): max|(τ(ε)-τ⁰)/ε-τ¹|/|τ¹|=2.67·10⁻⁴,5.34·10⁻⁴,1.07·10⁻³,2.14·10⁻³ at ε=1,2,4,8·10⁻⁶ (successive ratios exactly 2.000: O(ε)), with τ¹=⟨ℓ|T'_λ|r⟩/⟨ℓ|r⟩, T'_λ=∂_ε T_λ|_(ε=0) built from the principal minors (∂_εΠ_(i∈ I)zᵢ=Σ_(i∈ I)φᵢ), ℓ the left eigenvector in the weight sector. **τ¹ is linear in φ**, τ¹=Σₖφₖaₖ with aₖ=∂τ/∂ zₖ|_(z=1) (linearity checked at φ, φ' to 10⁻¹⁹⁶), and the scaling sum rule Σₖaₖ=|λ| τ⁰ holds to 10⁻¹⁹⁶ (z→ cz multiplies T_λ by c^(|λ|)) **[exact/numeric]**. The aₖ are rational for the non-degenerate branches, e.g. 27: T₁₀: (0,-27/140,4/35), T₁₁: (0,-3483/19600,-3483/19600) (the sum rule is re-checked here: -27/140+16/140=-11/140=τ⁰₁₀ and 2·(-3483/19600)=2τ⁰₁₁); singlet a₁=a₂=a₃=|λ|τ⁰/3 (symmetric weight); [4,1,1]: a₁=0 for T₁₀,T₁₁; octet aₖ are quadratic irrationals. No uniform closed formula was found beyond linearity and the sum rule.

**(c) τ⁰ through the untwisted dual Q's.** At z=1 the dual Baxter equation has a **3-dimensional polynomial kernel** in degrees ≤μ₁+2 for every branch; the staggered basis (elimination from the top degree) has degrees (μ₃,μ₂+1,μ₁+2)=(0,3,6),(1,2,6),(0,4,5),(1,3,5),(1,3,5),(2,3,4) for 27, [4,1,1], [3,3,0], octet 1, octet 2, singlet (checked against the expected degrees); Baxter residual of each vector ≤2·10⁻¹⁹⁹ (exact arithmetic on the stored numbers). **The untwisted quantum-Wronskian formula holds** **[numeric]**:

    T_λ(x)=P_λ(x) (det[v(x+λ₁h), v(x+(λ₂-1)h), v(x-2h)])/(det[v(x), v(x-h), v(x-2h)]),     P_λ(x)=Π_(boxes (r,c))ν₃(x+(c-r)h),

with v the three staggered untwisted dual solutions (any basis of the kernel: the ratio is basis-independent), ν₃(x)=Πₛ(x-θₛ+hΛ₃)=(x-θ₁)(x-θ₂) and P₁₀=ν₃(x), P₁₁=ν₃(x)ν₃(x-h), P₂₀=ν₃(x)ν₃(x+h), P₂₁=ν₃(x)ν₃(x+h)ν₃(x-h) (derived by Cramer from the 3×4 system): max|τ/(PW)-1|=0 at Precision 197 for all 6 branches × 4 shapes × 4 points x=θ₁-2h, 3/7, -2/9, 11/13; the control P=1 fails with O(1) deviation. Hence **at leading order D₂=cF(y)·(1,τ⁰) with τ⁰=PW[untwisted staggered dual Q's] is an explicit expression in untwisted Q's**, and at next order it is the same polynomial with τ→τ⁰+ετ¹ (F also depends on ε through y=(1+εω₂)y₀). The only non-Q input left at next order is τ¹ (a first-order eigenvalue shift of the untwisted chain, linear in φ).


### 6.3 det₆: what is and is not untwisted

Task 3 (cells 23–40) **[numeric]** unless stated.

**(a) Entry exponents.** In the raw 6×6 matrix (rows (k,(i,j)): k=1 bare rows (1,3),(1,2),(2,3), then k=2 periodic rows; columns c₁..c₆=-⟨ QF⟩,-⟨ xQF⟩,⟨ QF(x+h)⟩,⟨ xQF(x+h)⟩,⟨ Q(x-h)F⟩,⟨ xQ(x-h)F⟩) the exponent of every entry is p_(ab)=rₐ+s_b with s=(0,-1,0,-1,0,-1) (the x-weighted brackets are one power of ε larger) and row exponents (bare | periodic): 27: (5,3,3 | 0,-2,-2); [4,1,1]: (4,4,2 | 0,0,-2); [3,3,0]: (5,2,3 | 0,-3,-2); octets: (4,3,2 | 0,-1,-2); singlet: (3,3,1 | 0,0,-2). All 216 exponents are integers to 2.7·10⁻⁵. The raw k=2 rows contain divergent pieces (the twisted monic q² tends to a multiple of ε⁻¹q³ plus the finite untwisted solution of degree μ₂+1); the determinant-preserving row operations (row (1,2) `+=` y·row (1,3), row (2,3) `+=` x·row (1,3), both k; y=-Q²(a)/Q³(a), x=-F₂(a')/F₁(a') at a=1/3, a'=-1/5; y,x∼ε⁻²; the determinant is unchanged to 1.6·10⁻¹¹⁶) bring the periodic rows to exponents (0,1,1) for all branches.

**(b) The column basis: a rational one exists for the H directions.** Filtration of the column space by the weighted order in the three bare rows (levels Fⱼ={v:ordₐ(Mv)≥ rₐ-1+j}): dimensions (j=0..4) (6,5,3,1,0) for **all six branches**. The level F₂ (dim 3) is the same subspace for all six branches and for two directions (stacked rank 3, differences 10⁻¹⁵), with the rational reduced echelon basis

    H₁=c₁+c₅,     H₂=c₃-c₅,     H₃=c₄-5/2c₅-c₆  (5/2=5h at h=1/2; the 5h is untested),

and the complement O=(c₂,c₅,c₆) (any complement works). Change of basis CC (columns [H₁H₂H₃O₁O₂O₃]): det CC=1, det(M CC)/det CC=det M to 3·10⁻¹¹⁶. Entry-exponent profile of M CC in the bare rows: H columns (n+2,n+2,n+2), O (n,n+1,n) (e.g. 27, row 1: (6,6,6 | 4,5,4)). The combinations that vanish at leading order in the periodic rows (null space of the 3×6 leading periodic block, dim 3) are state- and direction-dependent and complex: **there is no rational O of that type; only the H space is rational**. The old thread's (c₁+c₃, c₅-c₆+5c₃, c₂-c₃) is the same space in its labelling.

**(c) The Laplace split holds for the top state only.** The block exponents add exactly: the ε-exponent of det₆ is 12,12,12,11,11,10 (for det₆ itself; D₆=det₆/Ψ_Q has 10,10,10,9,9,8; the difference of 2 is presumably the order of Ψ_Q, not checked in the sources) = BO + PH with BO (bare rows × O columns) =10,9,9,8,8,6 and PH (periodic rows × H columns) =2,3,3,3,3,4. The split ratio r_B=limdet₆/(-BO·PH): for 27 it is 1 (1+2.2·10⁻⁷ at ε=10⁻⁶, O(ε); limit 1+3·10⁻²³); for the other five it does **not** tend to 1: at (φ,ω A) it is -66/53=-1.2453 ([4,1,1]), 0.32794 ([3,3,0]), 0.72267 (both octets), 0.38888 (singlet), and it is direction dependent (spread over 48 points 2.8,1.8,0.75,0.75,6.7). Diagnosis: for the five non-top branches all 20 Laplace terms (subsets S of three columns for the bare rows) are of the same ε-order; the terms with S containing one H column and two O columns are of order 1 relative to the determinant (e.g. [4,1,1]: seven terms with |term/det| up to 15), because both BO and PH are anomalously suppressed (O₁,O₃ proportional at leading order; PH one order below the naive count). The bare-times-periodic product is therefore **not** the leading term in this basis (the determinant and the blocks are invariant under the row regularisation, so regularisation does not help).

**(d) The periodic block against the untwisted minor: fails.** Untwisted objects: dual q³,q² = the staggered kernel solutions (degrees μ₃,μ₂+1); the vacuum direct kernel at z=1 (the w→1 limit; the vacuum eigenvalues at w=1 equal the top state's, difference 0 at 198 digits): polynomial kernel dimensions 1,1,1,2,2,2,3 for maximal degree 0..6, i.e. staggered degrees (0,3,6)=(4-λᵢ)+(i-1) for λ=(4,2,0) (the finite-ε monic qᵢᵛᵃᶜ have degrees 0,2,4). The untwisted periodic minor PHᵤ (rows (1,3),(1,2),(2,3), insertion e^(2π ix/h), columns H₁,H₂,H₃, with F₁ (deg 0), F₂ (deg 3), q³,q²): **top 13453 (nonzero); the five others are numerically zero** (1.7·10⁻²⁰⁷, 0, 6.9·10⁻²⁰⁸, 7.2·10⁻²⁰⁸, 2.4·10⁻²⁰⁸; orthogonality of different irreps at zero twist). So PH(ε)∼ε^((3,3,3,3,4)) is controlled by first-order Q corrections and cannot tend to the untwisted minor. Even for the top state c₀(PH)/PHᵤ is not a monomial in the directions (spread 5): a direction-dependent factor common to PH, BO and D₆ of all branches (the state-independent constant r₆ of the det₆ formula; LHSₜₒₚ=⟨ΩΩ|e₁e₁⟩/⟨ x₀|e₁e₁⟩=ε⁻²/φ₁₂² exactly, and r₆=LHSₜₒₚ/D_(6, top) is not a monomial; D₆ is independent of the third vacuum twist w₃, checked). Best alternative tried: the top state's own untwisted direct solutions in place of the vacuum Q's change nothing (they coincide); no alternative reproduces PH.

**(e) The bare block.** BO itself is not a monomial (the common factor above), but the branch ratios are exact monomials with rational constants (spread ≤2.7·10⁻²⁰, exponent deviations ≤8·10⁻²⁰):

    BO_B/BOₜₒₚ=-5/3(ω₁-φ₂)/((ω₁-φ₃)(ω₂-φ₃)) ([4,1,1]),     -1/(ω₁-φ₂) ([3,3,0]),
    
    (5/2)/((ω₁-φ₃)(ω₂-φ₃)) (octets),     (5/2)/((ω₁-φ₃)²(ω₂-φ₃)²) (singlet).

So E_ν/E₂₇=-5/3,-1,5/2,5/2,5/2 for [4,1,1],[3,3,0], octet 1, octet 2, singlet; the old thread's E_λ/E₂₇ (for the norm; E_λ=-3/160,-1/36,-1/36,-2/3,-2/3,-1) would be 40/27,40/27,320/9,320/9,160/3 (re-derived here from the quoted E_λ): no resemblance (different object; our absolute E₂₇ is hidden in the common factor). PH_B/PHₜₒₚ is **not** a monomial (spread 4 to 170).

**Exact results for the full D₆ ratios** (not from the Laplace split) **[numeric]**:

    (D₆^([3,3,0]))/(D₆ᵗᵒᵖ)=(w₁-z₁)/(w₁-z₂),     (D₆^([4,1,1]))/(D₆ᵗᵒᵖ)=((w₁-z₂)(w₂-z₂))/((w₁-z₃)(w₂-z₃)),

hold **exactly in ε** (deviation 6·10⁻¹²² and 6·10⁻¹³⁰ at ε=10⁻⁶, Precision ∼140): ratios of Gauss parameters of the vacuum covector, with no Q input at all (this suggests an algebraic rather than SoV-determinant origin; the write-up lists it as a lead). For the octets and the singlet the leading monomial of the table above holds with a numerical O(ε) correction ρ_B=c₁/c₀ that is not a simple 9-term rational function of the directions (fit residuals 0.05 to 1; cell 43).


### 6.4 Reconstruction of Qₙᵤₘ to O(ε)

Task 4 (cells 41–44): Qₙᵤₘ^B=[D₆^B/D₆ᵗᵒᵖ][D₂^B/D₂ᵗᵒᵖ] topFac, with (1) D₂^B/D₂ᵗᵒᵖ=F(y)·(1,τ_B)/F(y)·(1,τₜₒₚ), τ_B=τ⁰_B+εΣₖφₖaₖ (untwisted τ⁰=PW of the untwisted staggered dual Q's, τ¹ from the untwisted chain), τₜₒₚ **exact** (the top state is the reference state; eigenvalues `auditVac`(z)), y=(1+εω₂)y₀ exact, c=343000/1161 cancels; (2) D₆^B/D₆ᵗᵒᵖ=ε^(p_B-10)c₀(1+ερ_B) with c₀ the Task-1 monomial (exact rational constant) and ρ_B=c₁/c₀ **numerical**, fitted from the four-ε ladder of the same direction; (3) topFac exact. Test at ε=3·10⁻⁶,5·10⁻⁶ (outside the ladder), direction A (φ, ω=(3/4,-1/5)) and direction B (φ', ω=(-5/3,9/7)); residual Qₙᵤₘʳᵉᶜ/Qₙᵤₘ-1 **[numeric]**:

| direction | branch | ε | full (c₁,τ¹) | leading only | c₁ only | τ¹ only |
|---|---|---|---|---|---|---|
| A | [3,3,0] | 3·10⁻⁶ | 6.2·10⁻¹³ | 7.2·10⁻⁷ | 7.2·10⁻⁷ | 6.2·10⁻¹³ |
| A | [3,3,0] | 5·10⁻⁶ | 1.7·10⁻¹² | 1.2·10⁻⁶ | 1.2·10⁻⁶ | 1.7·10⁻¹² |
| A | singlet | 3·10⁻⁶ | 8.9·10⁻¹³ | 1.7·10⁻⁶ | 7.2·10⁻⁷ | 1.0·10⁻⁶ |
| A | singlet | 5·10⁻⁶ | 2.5·10⁻¹² | 2.9·10⁻⁶ | 1.2·10⁻⁶ | 1.7·10⁻⁶ |
| B | [3,3,0] | 3·10⁻⁶ | 8.7·10⁻¹³ | 8.2·10⁻⁷ | 8.2·10⁻⁷ | 8.7·10⁻¹³ |
| B | [3,3,0] | 5·10⁻⁶ | 2.4·10⁻¹² | 1.4·10⁻⁶ | 1.4·10⁻⁶ | 2.4·10⁻¹² |
| B | singlet | 3·10⁻⁶ | 6.9·10⁻¹³ | 2.1·10⁻⁶ | 8.7·10⁻⁷ | 1.2·10⁻⁶ |
| B | singlet | 5·10⁻⁶ | 1.9·10⁻¹² | 3.5·10⁻⁶ | 1.5·10⁻⁶ | 2.0·10⁻⁶ |

The full reconstruction scales as ε² (ratio 2.7–2.8 for 5/3 in ε, i.e. (5/3)²), the leading-only one as ε (ratio 1.7). For [3,3,0] the D₆ ratio is exact, so only the τ¹ term matters. The fitted c₀ agrees with the monomial to 2·10⁻¹²¹ ([3,3,0]) and 10⁻²⁴ (singlet, limited by the cubic ladder extrapolation). **Caveat on “held out”:** the direction-B residuals are held out in ε and in the direction of the test, but ρ_B is fitted from the ladder of the same direction (`d6Coef[lab, io, ib]`, cell 42), so the O(ε²) agreement checks the structure, not a prediction of ρ_B.

**What is an untwisted-Q formula and what is still numerical** (cell 44):

- untwisted-Q formulas, verified: D₂ at leading order =cF(y)·(1,τ⁰) with τ⁰=P_λ W evaluated on the untwisted staggered dual Q's (residual 10⁻¹⁹⁷); topFac, τₜₒₚ (exact vacuum eigenvalues), F(y), y exact; D₆^([4,1,1])/D₆ᵗᵒᵖ and D₆^([3,3,0])/D₆ᵗᵒᵖ as exact rational functions of (w₁,w₂,z) (Gauss parameters), no Q input;
- untwisted-chain (not Q) formula: τ¹=Σₖφₖaₖ with rational aₖ;
- empirical (fitted monomial with exact rational constant, not derived from untwisted Q's): the leading coefficients of D₆ᵒᶜᵗᵉᵗ/D₆ᵗᵒᵖ and D₆ˢⁱⁿᵍˡᵉᵗ/D₆ᵗᵒᵖ (the periodic-block derivation failed);
- numerical: ρ_B (the first correction of D₆ for the octets and singlet; it vanishes for [4,1,1], [3,3,0]) and, for the octets, the individual copies' τ (quadratic irrationals).


---

## 7. Bugs and tooling

**The `Do`/`DownValues` purge bug (adjoint, cells 26–28, found in cells 34–36).** The first version of the t-scan for the singlet reported K/S=-(62/105)/(t-2). Cause: `DownValues` has attribute `HoldAll`, so the line `Do[DownValues[f] = ..., {f, {Tm, tBasisC, tZof, wgt, branchesAt}}]` of `purge[]` assigned to the `Do` iterator `f`, not to `Tm`, `tBasisC`, ...; the memoized transfer matrices and states of the first point (P1) were never purged and the scan kept P1's transfer matrices while y₀ and θ₁₂ changed. The result was K/S=-t₁₂(P1)/(3(t-2)) with 62/105=t_(P1)/3. **Symptom:** a clean-looking law whose numerator was a P1 number: at t=5/2 the scan gave -1.180952381 whereas the correct law gives -5/3. Everything at P1 itself (Tasks 1–3, 5 and the P1 part of Task 4) was unaffected. Fix (cell 36): purge with literal symbols, `DownValues[Tm] = Select[DownValues[Tm], !FreeQ[#[[1]], Pattern] &]` and so on, wrapped in `setPoint[h, theta]`; verified in fresh kernels (restart, cell 3 with the line `hN = 1/2; thN = {2/7, -3/5};` replaced by the test point, then cells 4–9, 14, 35): (h,θ)=(1,(5/2,0)): -1.6666666666672 against -5/3 (difference 5·10⁻¹³); (1/3,(1/5,-4/9)): 9.6666666666674 against 9.6666666666667 (7·10⁻¹³; the old scan gave 8.857142857); (1,(7/3,0)): φ: -2.3333333333339, φ': -2.3333333333343 against -7/3. The corrected eight-point scan agrees to 3.1·10⁻¹¹. The bug cost one wrong law.

**The (φ,ε,ω)-keyed memo bug ([3,1,0] round 2, Task B, cell 84).** `qtkC`, a memoized Q-side quantity depending on (h,θ) but keyed only on (φ,ε,ω), was not purged by `setPoint`; the first t-scan reused the Q-numbers of the first t, which made the “g slope” t-independent and produced a spurious (φ₃,ω₁=7/2) discrepancy. Fix: `qtkC` is now purged with literal symbols in `purge[]`; the P1 tables were computed before any scan (fresh cache) and the t-scan was redone. Task A used the non-memoized `qtk` and is unaffected. **Rule:** purge every (h,θ)-dependent memo with literal symbols.

**The ω₂=2/9 collision.** The spec's w₃=1+ε·2/9 coincides with w₂ at ω₂=2/9 (two equal vacuum twists, vacuum Baxter residual 4·10⁻¹⁰ instead of 10⁻¹⁹⁰). Use w₃=1+ε·5/13; w₃ drops out of the identity.

**Other gotchas found or re-met.** (1) `` Quiet[ClearAll["Global`*"]] ``: the front-end wrapper symbol `ClearGlobals` is Protected and aborts a bare `ClearAll`. (2) The naive Baxter residual prints an exact 0 with Precision ∞ (the trap of `Clean-CG-combinations.md` §7); evaluate it in exact arithmetic on the stored numbers (`L2-210-untwist-CG.wb` cell 8). (3) Precision limits: R_Q keeps 126 digits at working precision 200 (adjoint) and about 104–112 digits ([3,1,0]), so identities are quoted as 10⁻¹¹⁶ and 10⁻⁸², limited by the working precision, not by the identity. (4) Poles of structural constants at t=7/3 ([4,2,2]) and t=2 ([3,3,2]; y₀ also singular at t=2): exclude such t from fits. (5) t₁(θ₁+h) does not split the [4,3,1] pair; use t₂(θ₁+h). (6) τ scales as h⁴ (t₂ has degree 4 in (u,θ,h)): quote τ/h⁴. (7) Ω_((w₁,y)) degenerates at y=zᵢ; the Gauss parameters have poles at w₁=zₖ; Δ at w₂=zₖ. (8) Inherited: the string `y*)` closes a Wolfram comment; the notebooks write `y-star`. (9) Order of exponentials in n' (§2.3).

**No native outputs (all four notebooks).** VS Code did not execute the cells of new files (controller not associated, as for `L2-410-two-rotations.wb`): for each notebook an agent ran each cell's source in the shared kernel K1 (`ToExpression`/`Get` of the saved source, cells in order) and stored the printed results as markdown cells (`L2-210-untwist-CG.wb` cells 25, 33, 40, 45; `L2-210-untwist-Qdets.wb` cells 15, 22, 40, 44, 45; `Adjoint-untwist-perturbative.wb` cell 22; `L2-310-untwist-CG.wb`: a markdown cell with the printed output directly under every code cell, plus the summaries). `Adjoint-untwist-perturbative.wb` was created from the cell sources by a script (no create-notebook tool). **A rerun in VS Code is needed to store native outputs**; run times from the notebooks: adjoint, 21 states at one ε 1.5 s, one parameter point 0.02 s, one t point ≈3 s, a full run a few minutes; [3,1,0], 56 states 3.1 s per ε, 225 states 13.4 s, one t point ≈30 s (a 20-point scan ≈9 min); `Adjoint-untwist-perturbative.wb` a few seconds. Dependencies for a partial rerun of the adjoint notebook (its write-up §6): cells 3–9 always; Task 5 needs 11 and 14; the final cells need 3–9, 11, 14, 28, 36.

**Shared kernel.** Every one of these notebooks starts with a `ClearAll` that wipes the live state of the other new-conventions notebooks (files intact). The main session sequenced the agents; `L2-310-untwist-CG.wb` cell 52 records that “the kernel was cleared by another agent” between rounds and that setup cells 1–6 and the helper cells (10, 20, 21, 23, 24, 30–32) had to be re-run before round 2. Never evaluate two of these notebooks concurrently. Nothing was committed; `img/L2-210-untwist-CG/` is untracked tooling output (HANDOFF also records an unrelated 12-line change in `.vscode/settings.json`).


---

## 8. Open problems and suggested next steps

1. **Analytic untwisting of det₆ for the non-top branches** (§6). Needed: a column basis, possibly state- or direction-dependent or with O(ε) admixtures, that makes a bare×periodic term dominant; or a different route. The exact Gauss-parameter ratios for [4,1,1], [3,3,0] suggest an algebraic origin; the common non-monomial factor r₆ (state-independent, w₃-independent) is unidentified; τ¹ in Q-language (currently from the untwisted chain). Until then the leading coefficients of the octet and singlet D₆ ratios (constants -0.7767, 2.0514, -9800/2619) and ρ_B are empirical.
2. **Proofs and the missing constants for [3,1,0].** First- and second-order perturbation theory as in §4 for the fitted κ's, Δ_c^([5,2,1]), kslope; identify α(t) for [4,3,1]. Suggestion (mine, re-derived from the quoted traces and determinants): the conic parametrization of §4.3 should apply, since the discriminants are 64(3t²+1) ([4,3,1]) and 105t²+256 ([5,2,1]) **[conjecture]**.
3. **[4,1,0] at k=2** through the same pipeline (dressing of E₁₃-degree 2; four orders expected **[conjecture]**); then the [3,2,0]-type representations, where point A does not exist (`Undo-site2-rotation-notes.md` §7) and the undone covector is one direction short.
4. **Why is K/S of order one exactly for the branches whose slice weight is at ε^(≥1)?** The data are in §3.1 and §5.3; the sources give the exponent bookkeeping but no derivation. A derivation would also predict the constants (the perturbative proof of §4 is the adjoint instance).
5. **Relations to the earlier clean combinations.** κ_([4,4,0])=(t+1)/2 here coincides with the number (t+1)/2=97/70 quoted in `Clean-CG-combinations.md` §2.1 for the adjoint [3,3,0] (a different setting: T-insertion, not this overlap); and the [3,1,0] clean branch [6,1,1] has the same O(ε) constant as the adjoint [4,1,1] (§5.5). Whether these are structural or accidental is not known.
6. **Housekeeping.** Rerun the four notebooks in VS Code to store native outputs; decide with Paul whether the adjoint pipeline is promoted to a clean notebook; the HANDOFF NEXT list (1) analytic untwisting of det₆×det₂ (this section: det₂ done, det₆ open), (2) [3,1,0] k=2 then [4,1,0] (first done, second open), (3) derive κ_(singlet)(t), Δ_c(t,τ_c) (done for the adjoint), (4) rerun in VS Code (open).


---

## Appendix A. Notebook and cell maps

All notebooks are in `Paul/Mathematica/new-conventions/`; cells are 1-based. None has native outputs; “result cell” means a markdown cell holding the printed numbers. The first code cell of each notebook is `` Quiet[ClearAll["Global`*"]] ``.


### A.1 `L2-210-untwist-CG.wb` (45 cells)

| cells | content |
|---|---|
| 1–2 | title, setup note (everything a function of the twists; G makes P_(HW) the EH-frame orthogonal projector) |
| 3 | kernel reset; GT basis, plain generators `JJ`, relations, nilpotency, invariant metric G, Λ, `Ee` |
| 4 | monodromy and transfer polynomials: minors at the seven interpolation points `uS`, basis polynomials `tBasisC`, `tZof`; quantum determinant check |
| 5 | Gauss factors, x₀ with the factorized inverse (`x0Of`), `oneSiteVac`, `rhodiag` (the Greek letter rho in the notebook), `Dy`, the `u12` machinery (`u12rs`, `u12phi`, ...) |
| 6 | exact HW kets of V⊗ V, weights, G-orthogonal projectors per HW sector |
| 7 | Q machinery: poles, brackets `det6`, `su2det1`, `minor1`, `psiQ`; fixed-degree Q solvers (`auditVac`, `auditMonicQ`); state = association with dual Q's and twist |
| 8 | exact-arithmetic Baxter residual (`baxRes`) |
| 9 | all 21 eigenstates of the five HW sectors at a ket twist, ⟨ x₀\|Ψ⟩=1, HW-limit branches (`branchesAt`), t₁,t₂ polynomials |
| 10–12 | Task 0: ε-structure of the undone covector (symbolic ε), per-weight ranks over the 3×6 grid |
| 13–16 | Task 1: R_C, S, K, exact undo, ε-exponents, κ tables, contamination carriers (`showContam`) |
| 17–21 | Task 2: R_Q, state independence, Qₙᵤₘ, leading-order extraction, Richardson |
| 22–24 | Task 3: top-killing, ρ_B diagnostics |
| 25 | result cell, Tasks 0–3 |
| 26–28 | Task 4: singlet κ(t) scan (first version, with the purge bug), law check |
| 29–32 | Task 5: octet O(ε) layer, slope and intercept over 8 φ |
| 33 | result cell, Tasks 4–5 |
| 34–36 | Task 6: leak test; `purge[]`/`setPoint` fixed (cell 36), corrected scan |
| 37–39 | Task 7: Δ_c, τ_c scan over 22 t; rational fits |
| 40 | result cell, Tasks 6–7 |
| 41–45 | FINAL: `finalRun` (42), P2 (43), P3 (44), result cell (45) |


### A.2 `Adjoint-untwist-perturbative.wb` (22 cells)

| cells | content |
|---|---|
| 1–2 | title, regime, method (generic family member, conic parametrization of the octets, evidence tags), setup note |
| 3 | kernel reset; GT generators, metric, HW kets per weight sector |
| 4 | two-site transfer matrices at h=1, θ=(t,0) with t symbolic; principal quantum minors m_I; tₐ' |
| 5 | regime and covectors (exact in ε,t,φ,ω) |
| 6–7 | first-order perturbation theory (`pert1`), irrep projectors (`PiIrrep`) |
| 8–10 | singlet: ψ₁, K₁, S₁ (cell 9); Claim 1, second family member, per-irrep decomposition (cell 10) |
| 11–16 | octets: τ_c and the 2×2 matrix (12), ψ₁ and g_c (13), Claim 2 (14), per-irrep slope contributions (15), octet shares of the singlet K₁ (16) |
| 17–18 | bonus: clean branches [4,1,1], [3,3,0] |
| 19–21 | cross-checks: independent finite-ε evaluation at t=8 (20), second member, 27 and other octet, P1 (21) |
| 22 | result cell |


### A.3 `L2-310-untwist-CG.wb` (90 cells)

Every code cell is followed by a markdown cell “Output of the cell above ...” holding the printed output (in round 1; in round 2 some of these hold only `Null`).

| cells | content |
|---|---|
| 1–2 | title, setup note |
| 3–14 | setup (code, output): reset and generators with G (3–4); sparse transfer polynomials (5–6); Gauss factors, `x0Of`, point-A machinery (`dressBeta`, `dressGamma`, `vA1`, `vA2`, `phiA`), regime `reg` (7–8); HW kets and projectors (9–10); Q machinery (11–12); eigenstates `branchesAt` (13–14) |
| 15–18 | finite-twist regression (code 16, output 17, result 18) |
| 19–24 | Task 0: ε-structure (20–21), ranks (22–23), result (24) |
| 25–32 | Task 1: machinery (26–27), loop and classification (28–29), contamination carriers (30–31), result (32) |
| 33–38 | Task 2: Q side (34–35), leading extraction and Richardson (36–37), result (38) |
| 39–50 | Task 3: top-kill and ρ (40–41), exact (2,1,1) covector (42–43), [5,2,1] pair (44–45), 8 directions (46–47), F₀≡0 branches (48–49), result (50) |
| 51 | round-1 summary |
| 52–65 | round 2, Task A: `setPoint`, `scanT` (54), fits (56), κ closed forms (58–59), [5,2,1] pair over 20 t (60–61), rationalized closed forms (62–63), result (64), driver (65) |
| 66–84 | round 2, Task B: machinery (67), P1 structure (69–72), driver (73), t-scan analysis (74–75), affine fits (76–77), `kslope` fits (78–79), symmetric functions (80–81), α,β (82–83), result (84) |
| 85–90 | round 2, Task C: `taskC` (86), drivers (88; the P3 call is commented out in the source), result (89), summary (90) |


### A.4 `L2-210-untwist-Qdets.wb` (45 cells)

| cells | content |
|---|---|
| 1 | title, tooling note |
| 2–8 | setup copied from `L2-210-untwist-CG.wb` cells 3–9 |
| 9–15 | Task 1: machinery (10), data over 48 points (11), coefficients (12), monomial table (13), D₆ ratios (14), result (15) |
| 16–22 | Task 2: τ-polynomial identity (17), untwisted τ⁰ and τ¹ (18, 21), staggered kernel and quantum Wronskian (19–20), result (22) |
| 23–40 | Task 3: matrix and exponents (24–25), filtration and H space (26–28), regularisation (29–31), O space (32), block determinants and split (33–36), untwisted vacuum Q's and periodic block (37–38), normalisation r₆ (39), result (40) |
| 41–44 | Task 4: reconstruction (42–43), result (44) |
| 45 | summary |


---

## Appendix B. Checks summary

**What the notebooks check** (counts are literal `checks["..."] =` assignments in the sources; no notebook states a total pass count): `L2-210-untwist-CG.wb` 57, `Adjoint-untwist-perturbative.wb` 40 (the HANDOFF quotes 42), `L2-310-untwist-CG.wb` 56, `L2-210-untwist-Qdets.wb` 20. Key residuals:

- Adjoint (`L2-210-untwist-CG.wb`): gl(3) relations, nilpotency, positivity of G, quantum determinant 10⁻¹⁹⁹, [t₁,t₂] 10⁻¹⁹⁷, factorized x₀= Inverse-based x₀ exactly, top branch =e₁⊗ e₁ to 10⁻¹⁹¹; Baxter residuals 8.4·10⁻¹⁴⁹ (dual), 6.9·10⁻¹⁹⁰ (vacuum); exact undo 10⁻¹⁵² at 324 points; R_C/R_Q state-independence 1.6·10⁻¹¹⁶; Richardson Q against chain ≤6.4·10⁻¹⁹; octet mixing ≤8.4·10⁻¹¹; singlet law 3.1·10⁻¹¹ (8 t); Δ_c fit held-out ≤4·10⁻¹¹ (22 t).
- Perturbative (`Adjoint-untwist-perturbative.wb`): all symbolic identities of §4; finite-ε cross-check 2.6·10⁻⁹, 4.6·10⁻⁹ for ψ₁, 3·10⁻⁸/2.4·10⁻⁶ for Δ_c.
- [3,1,0]: regression 1.3·10⁻¹⁷⁹ (225 states); exact undo 4.9·10⁻¹³⁸ (540 points); R_C/R_Q 2.5·10⁻⁸²; Richardson ≤1.6·10⁻¹⁸; closed forms: see the table of §5.5.
- Qdets: τ-polynomial identity 0 at Precision 140; quantum Wronskian 0 at Precision 197; det under row operations and column change 10⁻¹¹⁶; D₆ ratios 6·10⁻¹²²; reconstruction 10⁻¹² at ε∼3·10⁻⁶.

**Independent checks made while writing these notes** (by hand, or with plain floating-point arithmetic in a shell; not Wolfram, not notebook cells; nothing in the repository was touched):

1. a' is exactly ε-free; the all-orders formula for b'/(a'c') and its first-order coefficient (§2.3); the first-order coefficient of the (0,1,2) component, -3(t-2)(φ₂-ω₂)(a'₀c'₀)² (§2.4).
2. The E₁₂² cancellation at ε=0 for [3,1,0] and β₀δ₂/δ₁=-a'₀/2 (§2.5).
3. S₁ of the singlet equals -3(t-2)(φ₂-ω₂)(a'₀c'₀)² (§4.2); 2A+B(τ₁+τ₂)=-t/(3(t-2)) (algebra), and the two copy shares at P1, 0.7307 and 1.8527, sum to 31/12 (arithmetic).
4. Δ₁+Δ₂ and Δ₁Δ₂ from Δ_c=A+Bτ_c and the τ-sum/product; the discriminant 9(5t²+4) and the conic parametrization (√(5t²+4)=2+mt, the deck involution) (§4.3, §3.3).
5. Clean-branch K₁/S₀ at P1 (0.10970, 0.08633) against the measured κ-1 (§4.4).
6. [3,1,0] closed forms evaluated at P1 against the printed numbers: κ_([4,4,0])=97/70, κ_([4,2,2])=1.366101695, κ_([3,3,2])=4.650803571, κ_c^([4,3,1])=1.410570, 0.986713 (with the printed τ_c), R_c^([5,2,1])=-0.89218864, 1.20406983, Δ_c^([5,2,1])=-0.34194023, 0.76579690, R^([4,2,2])=-0.485875706, R_c^([4,3,1])=-0.0932674, 0.5471062, kslope^([4,2,2])=-0.717449, x_U=167/86, 81/86; the traces and determinants against the printed τ's (trace 50.76082 and determinant 497.79369 for [5,2,1]; 45.44653 and 349.72475 for [4,3,1]); the discriminants 64(3t²+1) and 105t²+256 (§5.5).
7. Qdets: the τ¹ sum rule for 27; E_λ ratios; the exponent consistency p(Qₙᵤₘ)=p(D₆)+p(D₂)-16 (§6); the exponent bookkeeping of §2.6.

Not checked: any printed number was recomputed in Wolfram (none was); any fitted form was re-fitted; any claim about the old `SU3_L2_Untwisting` technology (not read).


---

## Appendix C. Source inconsistencies found while compiling

Reported, not resolved silently. “Notebook wins” means the notebook's printed number was used in the text.

1. **K₁₃ sign (Undo notes).** `Undo-site2-rotation-notes.md` §5.3 (md line 247; tex line 265) writes K₁₃=γ(w₂-z₁)(y_A-z₃)/w₂²=+189/704; `L2-310-untwist-CG.wb` cell 17 prints -189/704, and the point-A regression (cells 16–18) passes with the spec's formula. The notes carry a sign typo; the notebook wins.
2. **t=θ₁₂ versus t=θ₁₂/h.** The Undo notes write t:=θ₁₂=31/35 at P1 in the [3,1,0]/[4,1,0] formulas; this thread, the specs, HANDOFF and `Two-rotations-det-overlap.md`/`Clean-CG-combinations.md` use t=θ₁₂/h=62/35 (e.g. b'/(a'c')=1+ε(t-2)(φ₂-ω₂) and κ_([4,4,0])=(t+1)/2=97/70 need the latter). The scans use h=1, where they agree. This file writes t=θ₁₂/h throughout.
3. **Two different Δ_c.** For the adjoint octets Δ_c:=slope_c+(t-2)R_c (the signal contributes exactly -(t-2)R_c); for the [3,1,0] [5,2,1] pair Δ_c:=slope-p'₁-q'₁R_c with q'₁=(7-3t)/2, because the signal's coefficient of R_c is not -(t-2). The [3,1,0] spec (§2, Task 3) expected the adjoint combination “Δ:=slope+(t-2)R can be read off”; the notebook (cell 64) states that this combination is also affine in τ_c (Δ'_c=(216-647t/2+13t²)/105+(263-13t)/420(t+4)τ_c, the “coordinator's combination”) but is not the one that removes the signal. HANDOFF uses the symbol Δ_c for both without saying so. This file keeps the two definitions apart.
4. **[5,2,1] sanity values.** In `L2-310-untwist-CG.wb` the result cell 61 prints a line “t = 62/35 sanity” with the headings “Delta, Delta-prime”, namely {{0.71618427,-0.57000985},{0.57719571,-1.8534932}}, which agrees neither with the round-1 values (-0.34194023, 0.76579690) nor with the closed form (re-derived here: -0.341940, 0.765797); cell 64 prints “Δ=0.7655 for τ=37.479” (closed form: 0.7658; round 1: 0.76580) and calls the difference “within the quoted 10⁻⁹ noise”, which 3·10⁻⁴ is not. Likely a stale print or a typo; the closed form and round 1 agree with each other, and cell 61's line is unresolved.
5. **The t=8 check of Δ_c (adjoint perturbative notebook).** The write-up (§4.1) and HANDOFF say “Δ_c to 3·10⁻⁸”; the results cell prints -1.90909088 against -1.90909091 for τ=252 (difference 3·10⁻⁸) but -15.0000024 against -15 for τ=198 (difference 2.4·10⁻⁶).
6. **“Validated at P2, P3”.** `L2-210-untwist-CG.md` §0 item 3 and HANDOFF say the constants are validated at P2 and P3; cell 43 states that t=29/15 (P2) is in the Δ_c fit set (it is the first of the 22 values). Only P3 and the held-out t's are independent tests of Δ_c (for κ_(singlet), P2 has the same t as the first scan point but a different h, so it tests the h-independence).
7. **Rounded bounds in HANDOFF.** “Richardson limits Q vs chain ≤6·10⁻¹⁹”: the largest printed value is 6.4·10⁻¹⁹ (P2 singlet); “octet mixing ≤8·10⁻¹¹”: 8.4·10⁻¹¹ (P2); “≤3.6·10⁻⁶”: 3.61·10⁻⁶ (P2 octet 1). Rounding only.
8. **“42 passing checks”** (HANDOFF, perturbative notebook): the results cell states no count; the source has 40 literal `checks[...] =` assignments (59 occurrences of the string). Unverifiable.
9. **Cell count of the [3,1,0] notebook.** HANDOFF's earlier bullet says “51 cells”; the file has 90 (cell 51 is the round-1 summary; round 2 was added the same day). Same type as item 6 of the corresponding appendix of the Undo notes.
10. **Round-2 task message.** `L2-310-untwist-CG.md` header cites “the round-2 task message recorded in HANDOFF”; the HANDOFF entry records results, not the task message.
11. **Vacuum Baxter residual at small ε ([3,1,0]).** `L2-310-untwist-CG.md` §1 gives a single vacuum residual, 5.3·10⁻¹⁹⁵, next to the dual residuals of both the regression and small ε; that value is the regression point (cell 12), while cell 38 gives ≤2.6·10⁻¹⁸⁹ for the 54 small-ε twists. Minor.
12. **Per-weight orders in HANDOFF.** The [3,1,0] bullet lists the weights at ε¹ without the second direction of (2,1,1), which cell 24/51 puts at ε¹. Not a contradiction, an omission.
13. **Specs against notebooks (expectations that did not hold; all recorded by the agents).** `L2-210-untwist-CG-spec.md` Task 3 expected ρ_B=Kᵗᵏ/S to be ε-independent and φ-independent; for the octets ρ_B=O(ε) (ρ_B/ε is the structural object) and its intercept is φ-dependent. `L2-310-untwist-CG-spec.md` Task 0 expected (1,3,0) and one direction each of (1,2,1),(1,1,2) at ε⁰, and (0,2,2) at ε¹ (§0: “one direction of each dim-2 weight and the weight (0,2,2) are O(ε)”); found (1,3,0) and one direction each of (1,2,1),(1,1,2) at ε¹ (the E₁₂² cancellation), and (0,2,2) at ε². The Qdets spec expected the bare×periodic Laplace split to hold; it holds for the top state only. The Qdets spec asked for a 2×3 (ω₁,ω₂) grid; the notebook used the 18-point grid for two directions plus two further directions on six ω₂ (48 points). The adjoint spec's w₃=1+ε·2/9 collides with ω₂=2/9 (§7).
14. **“Depends only on (h,θ)” (c=343000/1161).** `L2-210-untwist-Qdets.md`'s statement is based on P1 only.
15. **“Held out” reconstruction (Qdets, Task 4).** Direction B is held out in the test point and in ε, but ρ_B is fitted from the ladder of the same direction (§6.4).
16. **Order of the exponentials in n'.** The Undo notes write n'=e^(c'E₂₃)e^(b'E₁₃)e^(a'E₁₂), `L2-210-untwist-CG.wb` cell 11 verifies φᵤ(0)=e₁ᵀe^(a'₀E₁₂)e^(c'₀E₂₃) (reversed order fails); these are consistent at b'₀=a'₀c'₀ because [E₁₂,E₂₃]=E₁₃ (§2.3). Likewise u'=e₁ᵀE₂₃E₁₂ (adjoint) and e₁ᵀE₁₂E₂₃ ([3,1,0]) differ because the dressing changes the ordering.
17. **Exponents of det₆ versus D₆.** `L2-210-untwist-Qdets.md` quotes 10,10,10,9,9,8 (for D₆=det₆/Ψ_Q) in §0 and block exponents summing to 12,12,12,11,11,10 (for det₆) in the notebook; the sources do not say that the difference is the order of Ψ_Q.
18. **Not found in the sources (not guessed).** The values of the α(t) polynomial coefficients (variable `caC`); the P3 Δ_c numbers of the [5,2,1] pair (only R^Q is tabulated); the P3 code call in cell 88 (commented out; its results are in cell 89); any check of c=343000/1161 or of the Qdets structure at P2/P3; why the 5h in H₃ is 5h; the content of `SU3_L2_Untwisting` (not read).
