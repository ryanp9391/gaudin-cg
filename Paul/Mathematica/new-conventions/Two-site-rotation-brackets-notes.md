# Two rotations on the two sites, as ladder brackets: site-2 mirror, two-site K tensor, untwisting of <Omega_w(x)Omega_(w')|, and the adjoint (inf,inf) bracket formula - gl(3) XXX, L=2, new conventions -- research notes of the 2026-10-11 thread (Stages A, B, D, C)

Paul Ryan (with Claude), 11 October 2026

Research notes on the 2026-10-11 thread executing `Two-site-rotation-brackets-plan.md` (Stages A, B, D and C; adjoint [2,1,0]^((x)2) in full for all four stages, [3,1,0]^((x) 2) in Stages A, B, D and as a cost note in C, [3,2,0] only in Stage A). Audience: Paul Ryan and future Claude sessions. These notes are a compiled reading of the four notebooks, the specs and the write-ups listed below; nothing was recomputed in Wolfram.

**Sources.** All paths are under `Paul/Mathematica/new-conventions/`. The four notebooks (cell numbers are 1-based positions in the `.wb` JSON file; every claim below is cited as [X:n] with the letter fixed by this table):

| **tag** | **notebook** | **cells** | **stage** |
|---|---|---|---|
| A | `L2-site2-rotation-brackets.wb` | 36 | Stage A: site-2 rotation, vacuum-family identification |
| B | `L2-two-site-rotation-brackets.wb` | 58 | Stage B: two-site K, R-matrix, blocks, families; follow-up Task 4b |
| D | `L2-two-vacua-untwist.wb` | 40 | Stage D: untwisting, reach, K/S, CG at leading order |
| C | `L2-210-two-vacua-brackets.wb` | 70 | Stage C: adjoint (inf,inf) bracket formula, Q-side CG |

The stage specs `L2-site2-rotation-brackets-spec.md`, `L2-two-site-rotation-brackets-spec.md`, `L2-two-vacua-untwist-spec.md`, `L2-210-two-vacua-brackets-spec.md`, the plan `Two-site-rotation-brackets-plan.md`, and the background `Two-rotations-det-overlap.md` S1, S4, `L2-210-untwist-CG.md` S0-S1, `L2-210-localising-span-spec.md` S6, `L2-210-site1-rotation-brackets.md`, `Site1-rotation-brackets-notes.tex` were read for context; the main session's running summary (`HANDOFF.md`, entry "NEWEST (2026-10-11)") was used only to locate things -- numbers are cited to notebook cells, and the interpretations marked [reading] are the main session's.

**None of the four notebooks carries a native output**: VS Code did not execute cells of the newly written files, every code cell was evaluated in the shared kernel from its saved source, and the printed results sit in markdown cells (immediately below the code cell in B and C; in the following "result" cell in A and D). Two numbering conventions occur inside the notebooks' own text (Appendix C, item 7); here all cell numbers are absolute positions in the file.

**Evidence tags.** [exact]: symbolic identity or exact rational computation (ranks, Grobner/`Factor`, exact rationals); [numeric, 10^(-n)]: true at the stated points to the stated deviation (60-digit finite-twist data in A, B, C part A; 200-digit chain data in D and C part B); [fitted]: closed form recognised from exact data at several parameter points, fitted on some and validated on the rest; [conjecture]: not verified; [reading]: interpretation by the main session, not tested by a notebook. Combined tags mean that part of the statement is exact and the rest numerical.

**Scope.** No notebook, no other `.md`/`.tex`, `HANDOFF.md`, `CLAUDE.md` or `paper/*.tex` was edited; nothing was evaluated in Wolfram; nothing was committed. Arithmetic done by hand while writing is marked "re-derived here". Places where sources disagree are collected in Appendix C and not resolved silently.

## 0. Summary

**What was asked** (Paul). Put two rotations, one on each site, on the SoV vacuum (or use a different site-2 reference than e_1^T), keep the number of ladder-bracket determinants polynomial in the number of boxes of the representation, and untwist with the two bra twists w,w'->1 independently (the old det_6xdet_2 family slaved the second vacuum to the first and degenerated).

**What was found.**

1. **The [3,1,0] five-bracket locus is the one-site vacuum family** [exact] ([A:15], [B:19]). The companion-frame vacuum Omega_c[w]=e_1^Texp(E_(12)/w_1)exp(-E_(13)/(w_1w_2))exp((w_1+w_2)/w_1w_2E_(23)) equals e_1^Trho(g) with (alpha,beta,gamma)=(1/w_1,-1/(w_1w_2),(w_1+w_2)/(w_1w_2)), first row (1,1/w_1,1/w_1^2), so beta+alpha*gamma=alpha^2. The adjoint quadric 2(t-2)(beta+alpha*gamma)=(t-3)alpha^2 is a torus image of this family ([A:31]).

2. **The SoV construction is site-asymmetric, and the asymmetry is the Yangian R-matrix** [exact] ([A:18], [A:20], [B:36]). e_1^T(x) e_nu^T expands on 35/64 (resp. 103/225) SoV states with all dim V site-1 labels, per-label bracket count 17/40 against 5/7 for the site-1 rotation. The R-matrix R^(Lambda*Lambda)(theta_(12)) on V(x) V (one-dimensional solution space, both reps) satisfies ~t=RtR^(-1), ~BR=RB; the dressed covector (e_1^T(x) c_2)R has the site-1 structure (support on (lambda_1|lambda_1c_1), count 5/7, adjoint conic constant cc_1(-t), [3,1,0] best locus beta+alpha*gamma=alpha^2), while the undressed site-2 conic is cc_2(t)=cc_1(1-t) [fitted, 15 (h,theta) points].

3. **The two-vacuum family Omega_w(x)Omega_(w') untwists cleanly** ([D:19]-[D:40]). The undone site-1 covector e_1^Trho(n(w)n(w')^(-1)) has all three Gauss parameters O(1) at eps=0 (rank 3 in the four direction parameters) [exact]; every weight of V sits at order eps^0 with rank equal to the multiplicity [exact]; K/S=O(eps) for every highest-weight branch (6 for the adjoint, 10 for [3,1,0]), the first gl(3) family with this property [numeric, 200 digits]; all CG data (6/6 and 10/10) are read at leading order from the overlap alone, errors ~eps, Richardson overlap-versus-chain agreement 10^(-16)-10^(-20) [numeric]. Pinning w_2' or fixing w' entirely does *not* degrade this; only w_1'=w_1 does (it reads one extra irrep and nothing else).

4. **Bracket cost is generic, not cheap, for [3,1,0]** ([B:41], [B:46], [B:49]). Pivot-decomposition counts sum P (upper bounds; sum r the matrix-rank lower bound), adjoint / [3,1,0]: generic pair 23/54, det_6 vacuum Omega_c[w](x)Omega_c[w] 18/40, Omega_c[w](x)Omega_c[w'] 22/54, R-dressed two-rotation candidates 23/54. Cheapness correlates with loss of reach (common w_1, global rotations), with one exception: the (inf,inf) family w_2,w_2'->inf, w_1!= w_1', which is CG-complete with only two parameters and costs 12 (adjoint; 10 after a full pivot scan, [C:32]) but 45 for [3,1,0].

5. **Adjoint (inf,inf) family in closed form, six CG data from the Q side** ([C:32], [C:57], [C:70]). Ten bracket products =20 Bk2 (=sum r, attained), exact scaling law A^((cc'))=a^(m+m'+e)a'^( n+n'+e')kappa^((cc'))(s), s=a'/a [exact]; kernels rational in t with poles only at t in {0,+-1} [fitted]; no non-trivial u-form [exact]; the Q-number equals the overlap on all 21 HW-sector states to 7.9*10^(-111), and the six CG data of adj(x)adj are reproduced from the Q side alone (Richardson versus chain <=1.2*10^(-18)) [numeric].

6. **Open decision**: for [3,1,0] no family found is both CG-complete and near the planned target (~25 products); the best searched decomposition of the (inf,inf) family costs 45 (lower bound 24, minimality unproved; [C:69]). Whether to accept the block formula, change the cost metric (bracket determinants of size >2; det_6 packs 36 bracket entries), or look for a determinant form of <Omega_w(x)Omega_(w')| with non-vacuum bra data is Paul's call (S6).

## 1. Setting and conventions

**Chain and frames.** gl(3), L=2 sites, representation Lambda=(lambda_1,lambda_2,0): [2,1,0] (adjoint, dim V=8, 64 states on two sites) and [3,1,0] (dim V=15, 225 states); [3,2,0] (dim V=15) only in Stage A. New conventions: L_(ij)(u)=(u-theta_s)delta_(ij)+hE_(ij), T(v)=L^((1))(v-theta_1)L^((2))(v-theta_2), twist on the right, plain exact Gelfand-Tsetlin generators. Numeric point throughout: h=1/2, theta=(2/7,-3/5), hence t:=theta_(12)/h=62/35; ket twist z=(3/2,5/4,2/3) (controls z'=(3/2,5/4,3/7), z''=(2,5/4,2/3)); bra twists w=(5/2,7/3,3/5), w'=(7/4,-2/5,9/7) in B ([B:16]); 60-digit eigenvectors in A, B, C part A, 200 digits in D and C part B ([A:1], [B:1], [D:1], [C:1]). The *diagonal frame* has twist diag(z). The *companion frame* is reached by (P_( phys)^(-1))^((x)2) with P_( phys)=P_-P_dP_+; there the SoV vacuum is x_0=e_1^T(x) e_1^T ([A:8]: < x_0^C|ket> equals the diagonal-frame value to 10^(-54)/10^(-49)). The bracket formulas live in the companion frame; the chain-side CG extraction lives in the diagonal frame. The dictionary ([D:15], [C:57]) is the frame tie Omega_c[w]|| Omega_D[z,w] P_( phys)(z) up to the scalar kappa(z,w)=prod_kDelta_k(nP_3)^(lambda_k-lambda_(k+1)) (not simply equal; Omega_c is normalised to first component 1), and <Phi_c|Psi_c>/< x_0^c|Psi_c>=<Phi_d|psi_d>/(kappa*kappa') with kappa,kappa'~eps^(-4) each.

**SoV labels, wavefunctions, ladders.** The untwisted B(u) has a simple spectrum; its left eigencovectors <x| carry a label (a|bc) on each site with roots theta_s+{-1-a,-b,-c}h and the rule lambda_2<= a<=lambda_1, 0<= c<=lambda_2, c<= b<= a (reproduced for both sites, [A:8]; dim V labels per site). The wavefunction is < x|Psi>/< x_0|Psi>=N_xprod_s M^((s))_x/M^((s))_(vac) with the monomial M^((s))_((a|bc))=minor_s(-1-a,-b) v_3(theta_s-ch), minor(n,n')=Q^2(n')Q^3(n)-Q^3(n')Q^2(n), v=(Q^1,Q^2,Q^3) the dual Q's (verified on all states to 10^(-43)-10^(-46) for [2,1,0] and 10^(-39)-10^(-43) for [3,1,0] in three row/k conventions, [A:8], [B:14]). Ladder s has the points x_m=theta_s+mh, m=-1-lambda_1,...,0 (lambda_1+2 points); < Q^jF>_s=sum_mF_mQ^j(x_m) and

    Bk2^((s))[F,G]=< Q^3F>_s< Q^2G>_s-< Q^2F>_s< Q^3G>_s=sum_(n<n')(F^ G)_(nn') minor_(23)(n,n')

(verified, [A:20]). A 2-form on the ladder points is *decomposable* (one wedge, one Bk2) iff its Pfaffian vanishes; a 2-form of rank 2k costs k wedges.

**Two-site covectors, K tensor, blocks.** Every covector Phi on V(x) V is expanded as Phi=sum_(x_1,x_2)phi_(x_1x_2)<x_1x_2|/N_(x_1x_2) with phi=Phi* K, where K_(nu_1nu_2;x_1x_2) is the exact rational *two-site K tensor* built in B ([B:24]); <Phi|Psi>/< x_0|Psi>* M^((1))_vac M^((2))_vac=sumphi_(x_1x_2)M^((1))_(x_1)M^((2))_(x_2). Each label carries on each ladder a *minor* (m,m') and a *lone point* (site 1: Q^3(theta_1-ch), site 2: Q^3(theta_2-c'h)), c,c' in {0,...,lambda_2}. The *blocks* A^((c,c')) are the n_(min)x n_(min) matrices (n_(min)=C(lambda_1+2,2)=6 resp. 10) of phi with rows the ladder-1 minors and columns the ladder-2 minors; there are (lambda_2+1)^2 of them, so the count below is "polynomial in the boxes by construction". Two counts are reported throughout ([B:1], [B:41]): r=sum_(cc')rank A^((cc')) and P=sum_(cc')P_(cc'), where P_(cc') is the cost sum_k(rank u_k/2)(rank w_k/2) of a pivot decomposition A=sum_ku_k(x) w_k into ladder-1 and ladder-2 2-forms (the smaller of the two pivot orientations). Since every product Bk2^((1))Bk2^((2)) contributes a rank-one matrix, r<=(minimal number of products)<= P [exact]; P depends on the pivot choice and is not monotone in r ([B:58], item 10(i)), and the minimal count is a tensor-rank problem that was not solved except where r=P is attained. Stage A also uses a different, coarser *per-label* count (ladder-2 brackets with the ladder-1 monomial of each site-1 label as a scalar factor).

**Gauss chart and vacua.** g=e^(alpha E_(12))e^(beta E_(13))e^(gamma E_(23)), first row of the defining-representation matrix b=(1,alpha,b_3), b_3=beta+alpha*gamma; "genuine" means alpha*beta*gamma!=0. One-site vacua: companion Omega_c[w] (above), polynomial in (a,s)=(1/w_1,1/w_2) as gaussCov{a,-as,a+s} ([B:46]); diagonal frame Omega_D[z,w]=e_1^Trho(n(z,w)), n=exp(r_1E_(12))exp(r_2E_(13))exp(r_3E_(23)) with the rational r_k(z,w) ([D:1]); w_2=inf is the limit, gaussCov{a,0,a} ([B:49], [D:15]). Omega(x)Omega is an exact left eigencovector of the transfer matrices with the companion-matrix twist G(w_1,w_2,w_3) for any w_3 ([A:15]; control G^T fails), in the diagonal frame with G'(w)=P_( aux)G(w)P_( aux)^(-1), not diag(w) ([D:15]).

**Regime.** z=1+eps*phi, w_k=1+eps*om_k, w'_k=1+eps*om'_k independent, eps in {1,2,4}*10^(-6) (S4). K is exactly z-independent ([A:34], [B:57]), depends on (h,theta) only through t (rational) and is homogeneous of degree 0 in (h,theta) ([C:65]).

## 2. Stage A: the vacuum family and the site-2 mirror

Notebook A (36 cells; summary [A:36]); Stage A asks two things: what the site-1 five-bracket locus is, and what the mirror construction Phi=e_1^T(x) e_1^Trho(g_2) (rotation on site 2, vacuum on site 1) looks like.

### 2.1 The vacuum family is the [3,1,0] locus

[exact] ([A:15]; same results in [B:19], [B:16]). For both reps and four resp. three rational w: Omega_c[w] (dense `MatrixExp` of the rep generators) equals gaussCov{1/w_1,-1/(w_1w_2),(w_1+w_2)/(w_1w_2)}; at the 3x3 level the first row of the Gauss product is (1,1/w_1,1/w_1^2), hence beta+alpha*gamma=alpha^2 symbolically. Converse: six rational points (alpha,gamma) with beta=alpha^2-alpha*gamma (including the points (-3,3,-2),(-3,6,-1),(-3,12,1),(-3,15,2) of the site-1 table) equal Omega_c[{1/alpha,1/(gamma-alpha)}] exactly. So the [3,1,0] best locus is {Omega_c[(w_1,w_2)]:w_1=1/alpha, w_2=1/(gamma-alpha)}. On the adjoint the best site-1 locus has cc=(t-3)/(2(t-2))!=1 and the vacuum family costs the generic 5 brackets. Omega_c[w](x)Omega_c[w] is an exact left eigencovector of t_1,t_2,t_3 at u=3/7,19/23 with the companion twist G(w_1,w_2,w_3)=([[chi_1,-chi_2,chi_3],[1,0,0],[0,1,0]]) for w_3=3/2,-4,7/3 and a second w (any w_3); the control G^T fails for a=1,2.

**The two [3,1,0] lambda_2-group measures with alpha=1/w_1** [exact wedge identities; numeric 10^(-43.7)-10^(-46.6) on all 225 states] ([A:15], cell 14 of the notebook's code): with P_3',P_4' the node-polynomial derivatives (supports {-4,-3,-2}, {-4,...,-1} of ladder 1),

    column (0,1):    F(u)=(h^4 w_1^(-(u-theta_1)/h))/(P_3'(u)(u-theta_1)(u-theta_2+h)),  G(u)=(h^3 w_1^(-(u-theta_1)/h))/(P_3'(u)(u-theta_2+h)),  amplitude (48(t-3)(t+1))/(t w_1^8),
    column (0,0):    F(u)=(-h^4 w_1^(-(u-theta_1)/h) om(u))/(P_4'(u)(u-theta_2+h)),  G=F (u-theta_1)/h,  amplitude (-12(t-3)(w_2+3w_1))/(w_2 w_1^9),
with om=1 except om(theta_1-h)=3w_1/(w_2+3w_1) (symbolically 3(s-1)/(3s-2), s=gamma/alpha=1+w_1/w_2; no t-dependence). So both read as (bra twist w_1^(-(u-theta_1)/h))x(SoV measure of the support)x 1/(u-theta_2+h)x(low-degree polynomial) [exact]. The amplitudes at t=62/35 are -100104/(1085 w_1^8) and 516(3w_1+w_2)/(35 w_1^9w_2).

### 2.2 The site-2 rotation is not the mirror image: support and counts

[exact support; numeric 10^(-41)-10^(-46) identities] ([A:18], [A:20]). For Phi=e_nu^T(x) e_1^T (site 1 rotated) only the site-2 labels (lambda_1|lambda_1c_2) occur (12 resp. 24 SoV states). For Phi=e_1^T(x) e_nu^T the expansion involves *all* dim V site-1 labels and all site-2 labels on a "triangular" set of pairs (Table 1); the vacuum label on either site pairs with every label on the other, the set shrinks monotonically towards the corner label (1|00), and no additive weight rule exists: a naive rule wt(x_1)+wt(x_2)=Lambda+wt(nu) is violated by 11/51 nonzero K entries ([2,1,0]/[3,1,0]) and a linear solve for any additive assignment is inconsistent ([A:18]; the support rule could not be characterised more sharply in two attempts).

**Table 1.** Site-2 mirror: supports and counts ([A:18], [A:20], [A:35]); site-1 numbers from [A:12].

|  | [2,1,0] | [3,1,0] | [3,2,0] |
|---|---|---|---|
| SoV states in the support of e_1^T(x) e_nu^T (of dim V^2) | 35 of 64 | 103 of 225 | 122 of 225 |
|   of which with x_1 in {(lambda_1/lambda_1c_1)} | 15 | 29 | 41 |
| site-1 rotation: support | 12 | 24 | -- |
| per-label automatic count, site 2 (generic = random covector) | **17** | **40** | **54** |
| site-1 automatic count, generic / reduced | 5 / 4 | 7 / 5 | 11 / 9 |
| block form: matrix ranks r(c_1,c') | 3,3,3,2 | 6,6,6,4 | sum r=30 |
| block form: products (explicit; random covector) | 15 (16; 16) | 35 (38; 35) | 49 (52; 50) |
| z-independence of the maps K^((x_1)) (exact, z',z'') | True, True | True, True | True, True |

The per-label count 17/40 exceeds dim V=8/15 in the [3,1,0] case: the ladder-2 bracket formula for a site-2 rotation on the vacuum is not a few-bracket expression, because the ladder-1 factor is the full monomial minor_1(-1-a,-b)v_3(theta_1-ch) of every label rather than a single Q^3 ratio ([A:18]). It is still set by ladder lengths and column ranks (<=4), and the block-form 35-38 products at [3,1,0] are of the order of the det_6 Laplace count ([A:20]). [exact]

### 2.3 Loci, reach and the best site-2 locus

[exact] ([A:24]; all Grobner solves <1 s). *Adjoint.* Only two columns are non-automatic: (x_1=(2|20),c'=0) with Pfaffian ~ g_B^2(g_B+g_Ag_C)^2, genuine part b_3=0; and (x_1=vac,c'=0) with Pfaffian ~ g_B^2(66g_A^2-97b_3) at t=62/35, genuine part the conic b_3=cc alpha^2, cc=66/97. The joint locus of genuine rotations is empty. The constant is cc_2(t)=(t+2)/(2(t+1))=cc_1(1-t) with cc_1(t')=(t'-3)/(2(t'-2)) the site-1 adjoint constant (re-derived here: cc_1(1-t)=(-t-2)/(2(-t-1))) [fitted, exact on 15 (h,theta) points]; of the shifts t'=t+k, k-t (k=-6..6) only t'=1-t matches; t'=-t, t+-1, 2-t do not. The constant is not 1: the site-2 conic is not the vacuum family. *[3,1,0].* 11 non-automatic columns; conic-type factors b_3=cc alpha^2 with exact t-laws (2t^2+5t-4)/(2(t-1)(t+2)), ((t-1)(8t^2+13t-7))/(6(t-2)t(t+1)), (t^2+3t+6)/(3t(t+1)), (t+2)/2t, (2(t+3))/(3(t+1)) (values 2273/1188, -453483/96224, 8852/9021, 33/31, 334/291 at t=62/35) [fitted, 15 points]; cc=1 occurs nowhere. Single 5-point columns ((3|30),0) and (vac,0) have empty genuine part; of 57 minimal empty sets, 25 maximal sets are jointly decomposable: 24 pairs (each a dim-1 curve; rational points found for 4) and one 4-set {(6,0),(6,1),(8,0),(11,0)} with locus b_3=0.

**Table 2.** Site-2 loci: per-label counts and finite-twist reach ([A:24], [A:35]).

| locus | [2,1,0] | [3,1,0] | [3,2,0] | reach |
|---|---|---|---|---|
| generic | 17 | 40 | 54 | all weights, rank 2 on doubles |
| vacuum family b_3=alpha^2 | 17 | 40 | 54 | all |
| adjoint-type site-1 quadric | 17 | 40 | 54 | all |
| b_3=0 | **13** | **28** | **43** | **incomplete**: 5/7, 7/12, 9/12 weights |
| conic 66alpha^2=97b_3 (CG-complete) | 16 | -- | -- | all |
| dim-1 curves with rational points (CG-complete) | -- | 38 | 52-53 | all |
| non-reduced curve | -- | -- | 42 | 9/12 weights |

The undo identity < e_1^T(x) g_2|h>=< e_1^T(x) e_1^Trho(g_2L_0D_0)|h> holds exactly on all HW kets (6 resp. 10 HW states, 5 resp. 8 irreps; g_2=(1,2,3),(-3,6,-1),(2/3,-3/7,5/4)), because rho(U)(x)rho(U) fixes the HW kets ([A:24]) [exact]. **Key answer** ([A:24]): the site-2 best locus is *not* the vacuum family: the cheapest locus b_3=0 is not CG-complete, and the CG-complete reducing loci save only 1-2 brackets (16 against 17; 38 against 40).

### 2.4 Closed forms on b_3=0 and the exponent law

[exact closed forms; numeric 10^(-44)-10^(-46) on all states at three fresh points] ([A:29]). On b_3=0 (g_B=-g_Ag_C, parametrised g_2=e^(gamma E_(23))e^(alpha E_(12))) the columns decompose exactly into 13 wedges (12 columns, [2,1,0]) resp. 28 wedges (21 columns, [3,1,0]), residual 0 for every column; the wedges are sparse (|supp G|=1 for all but two, i.e. a point value of Q; |supp F|<=4); the explicit adjoint formula is printed in [A:29]. Exponent law (with g_C=s_Rg_A): A_(mm')=g_A^(m+m'+k)a_(mm')(s_R,t), a of degree <=1 in s_R, offset k one value per column (k=8,...,5 for the adjoint, 11,...,7 for [3,1,0]); each coefficient is a rational function of t, fitted on 15 (h,theta) points and reproduced exactly (poles at t in {0,+-1} for the adjoint, {0,+-1,+-2} for [3,1,0]; numerator factors t+k, k=-6..4, i.e. the values theta_1-x_m=h(t-m) of the ladder-1 pole polynomial at the ladder-2 nodes) [fitted]. u-form tests ([A:29]): no column is natural in the global SoV-measure, polynomial or cross-pole classes for degree <=1; the rank-2 three-point columns are delta_(-4)^(2-point vector) with a t-dependent annihilator ratio (e.g. r=-8342/13527 at t=62/35, (2|20),1) rather than the site-1 form n_m=-m(t+m+1). So the expected w_1'^(-(u-theta_2)/h)xmeasurexcross polexpolynomial form exists only as the exponent law [exact].

### 2.5 The adjoint quadric is a torus image of the vacuum family

[exact] ([A:31]). e_1^Trho(g)rho(D)=rho(D)_(11)e_1^Trho(D^(-1)gD) and D^(-1)gD has Gauss parameters (alpha d_2/d_1,beta d_3/d_1,gamma d_3/d_2), so b_3 scales by d_3/d_1 and alpha^2 by (d_2/d_1)^2: the torus with d_1d_3/d_2^2=cc maps {b_3=alpha^2} onto {b_3=cc alpha^2}. Verified on the rep (t=62/35, cc=43/16, d=(1,4,43); three w; converse with d=(1,1,cc)). Hence the adjoint quadric is {Omega_c[w]rho(D):d_1d_3/d_2^2=(t-3)/(2(t-2))} and the [3,1,0] best locus is the vacuum family modulo the unimodular torus d_1d_3=d_2^2. With kappa_( singlet)=2(t-3)/(3(t-2)): cc=(3/4)kappa_( singlet), i.e. d_1d_3/d_2^2=3/4kappa (noted, nothing more). The site-1 factor of the det_6xdet_2 product point of `Two-rotations-det-overlap.md` S5 has b_1b_3/b_2^2=1.5013, 1.0196, 0.9933 at three w: a non-constant symmetric function of (w_1,w_2), on the quadric only along a curve, and not on the vacuum family; the diagonal-frame resolvent ratio (y-z_2)^2/((y-z_1)(y-z_3)) equals cc(t) only at y=0.43438, 1.53476, not at the values used there: no identification [numeric].

### 2.6 [3,2,0] (optional, [A:35])

Same picture [exact/numeric]: support 122 of 225 (41 with x_1 in  F_1), per-label count 54 against site-1 11/9; the vacuum-x_1 conic constants are 2(t+2)/(3(t+1)) (c'=0) and (t+2)/(2(t+1)) (c'=1), the latter the same function as the adjoint's (exact over 15 points); b_3=0 gives 43 (35 columns, 43 wedges, exponent offsets 12 to 8, poles t in {0,+-1,+-2}) but reaches only 9 of 12 weights; CG-complete curves cost 52-53.

## 3. Stage B: the two-site K tensor, the R-matrix, block counts

Notebook B (58 cells; summary [B:58]; follow-up Task 4b in [B:47]-[B:49]).

### 3.1 The K tensor and supports

[exact] ([B:24]). K is a 64x64 resp. 225x225 exact rational array (all d^2 GT covectors e_(nu_1)^T(x) e_(nu_2)^T expanded in the SoV basis; expansion error 1.4*10^(-89) / 7.6*10^(-150), rationalisation error 3.4*10^(-74) / 3.6*10^(-73) against the 60-digit data, denominators <=8 resp. 12 digits). It is very sparse (278 of 4096 resp. 2051 of 50625 entries nonzero), of full exact rank, the state-label correspondence is bijective, and it is z-independent exactly (rebuilt at z',z'', [B:57]) and rational in t (rebuilt at two more (h,theta)). The slice nu_2=1 equals the one-site site-1 maps K^((c_2)) recomputed by an independent solve (support 12/24 states, ranks 4,8 / 9,15), and the slice nu_1=1 equals Stage A's site-2 maps (35/103 states, 15/29 with x_1 in  F_1, ranks per site-1 label 1,2,4,2,4,6,4,8 / 1,2,4,2,4,6,5,9,2,5,9,12,4,9,15). The overlap identity from phi holds to 10^(-45) / 10^(-42).

**Supports of two-site covectors** ([B:29]; exact zero tests). Every c_1(x) c_2 with generic factors, and also Omega_c[w](x)Omega_c[w'] (including w=w'), has the *full* support, 64 resp. 225 label pairs, every (c,c') block fully populated (slots 25,15,15,9 resp. 81,54,54,36). There is no dependence of the site-1 label set on the site-2 label: "the measure couples the two ladders" of `Two-rotations-det-overlap.md` S4 means that the support is everything and all information is in the coefficients (phi has rank 8 resp. 15 as a dx d matrix; the vacuum Omega_c[w](x)Omega_c[w] has rank 8 resp. 12). The generic-pair support contains the union of the two one-site supports but is strictly larger.

### 3.2 The R-matrix: the site asymmetry is the Yangian conjugation

The reading to be tested (main session, [reading]): the reversed-order chain L^((2))L^((1))=RTR^(-1) with the Yangian R^(Lambda*Lambda)(theta_(12)), so that ~t=RtR^(-1), ~Psi=RPsi and the SoV wavefunctions are identical, and mirror symmetry makes the cheap site-2 covectors the R-dressed ones.

**Construction and checks** [exact] ([B:36]). R solves R T_(ab)(u)=~T_(ab)(u)R for all nine (a,b), ~T=L^((2))(u-theta_2)L^((1))(u-theta_1), inside the commutant of the diagonal gl(3), i.e. the algebra generated by 1, C_(12)=sum E^1_(ab)E^2_(ba), the cubic couplings D_(12),D_(21) and two quartic couplings Q_1,Q_2; its dimension is 8 resp. 14=sum(HW multiplicity)^2 on V(x) V (checked), and the quartic invariants are needed (words in 1,C_(12),D_(12),D_(21) alone give 6 of 8 for the adjoint). The Gram matrix of the flattened residuals at u=3/11,2/9 has a one-dimensional null space for both reps. Then (all True for both reps): RT_(ab)=~T_(ab)R at a third u=19/23; ~t_a(u)R=Rt_a(u) with twist diag(z); ~B(u)R=RB(u); R commutes with rho(g)(x)rho(g); (e_1^T(x) e_1^T)R=e_1^T(x) e_1^T after normalisation; R R^(-1)=1. R has 322 resp. 2285 nonzero entries (denominators <=7 / 10 digits) and R->1 as theta_(12)->inf (max|R-1|=0.86/1.02 at t=62/35 and 0.020/0.044 at t=1404/7, i.e. O(1/t), [B:36], [B:57]).

**Dressed test** [exact support and counts] ([B:36]). The correct side is Phi=(e_1^T(x) c_2(g_2))R (row vector times R):

|  | [2,1,0] | [3,1,0] |
|---|---|---|
| (e_1(x) c_2)R: support; all with x_1 in  F_1? | 12 of 64; True | 24 of 225; True |
| per-label count, generic g_2 / g_2=Omega_c[w'] (dressed) | 5 / 5 | 7 / 5 |
| the same, undressed e_1(x) c_2 (Stage A) | 17 | 40 |
| (e_1(x) c_2)R^(-1) (wrong side): support; count | 36 pairs; 17 | 106 pairs; 41 |
| control (c_1(x) e_1)R (dressing the cheap site-1 rotation): support; counts | 35 pairs; 17, 17 | 103 pairs; 40, 41 |

So the R-dressed site-2 rotation has the site-1 structure (support only on x_1 in  F_1, factors v_3(theta_1-c_1h)/v_3(theta_1-lambda_2h), count 5/7), and dressing the cheap site-1 rotation makes it expensive: the asymmetry of Stage A is the R-matrix conjugation theta_(12)->-theta_(12) [exact]. *Loci of the dressed covector* (symbolic Gauss chart, [B:36]): adjoint, only the column (x_1=(2|20),c'=0) is non-automatic, Pfaffian =const* g_B^2(167g_A^2-264 b_3); the conic constant is cc=167/264=cc_1(-t)=(t+3)/(2(t+2)) at t=62/35 (re-derived here), also 37/59 at (h,theta)=(1/3,(1/5,-4/9)), 61/94 at (1/2,(3/7,-1/4)), 1425/2836 at theta_2=-100 ([B:36], [B:57]); count 4 on the conic (also on g_B=0). [3,1,0]: conic constants {1,404/501,167/264,3371/7348,101/66}; cc=1, the vacuum family, gives the best count 5, the mirror of the site-1 result (the constant 167/264 also appears, count 7).

### 3.3 Block counts: the main table and the det_6 regression

[exact ranks and pivot counts; numeric identities 10^(-39) or better] ([B:41], [B:58]). Blocks (c,c')=(0,0),(0,1),(1,0),(1,1). sum r / sum P per row:

**Table 3.** Block counts sum r / sum P ([B:41]; summary [B:58], item 4).

| Phi | [2,1,0] sum r | [2,1,0] sum P | [3,1,0] sum r | [3,1,0] sum P |
|---|---|---|---|---|
| (a) generic pair c_1(g_1)(x) c_2(g_2) (= random covectors) | 14 | **23** | 27 | **54** |
| (c) Omega_c[w](x)Omega_c[w] (=det_6) | 13 | **18** | 26 | **40** |
| (d) Omega_c[w](x)Omega_c[w'] (and reverse) | 13 | **22** | 27 | **54** |
| (e) Omega_c[w](x) e_1 | 4 | 5 | 4 | 5 |
| (e') e_1(x)Omega_c[w'] | 11 | 15 | 22 | 35 |
| (f) c_1(g_1)(x) e_1 generic | 4 | 5 | 4 | 7 |
| (f') e_1(x) c_2(g_2) generic (= Stage A block form) | 11 | 15 | 22 | 35 |
| (e3) c_1(adjoint quadric point)(x) e_1 | 4 | **4** | -- | -- |
| (g) (e_1(x)Omega_c[w'])R (R-dressed vacuum) | 4 | **5** | 4 | **5** |
| (g2) (e_1(x) c_2(dressed conic))R | 4 | **4** | -- | -- |
| (h) R-dressed two-rotation candidates (e_1(x)Omega_c[w'])R(rho(g_1)(x)1), (Omega_c[w](x) e_1)R(1(x)rho(g_2)), (e_1(x) c_2)R(rho(g_1)(x)1) | 14 | **23** | 27 | **54** |

Per-block ranks: generic pair (5,3,3,3) / (9,6,6,6); det_6 vacuum (5,3,3,2) / (8,6,6,6) with P per block (10,3,3,2) / (22,6,6,6); Omega_c[w](x)Omega_c[w'] (5,3,3,2) with P=(10,5,5,2) / (9,6,6,6) with P=(18,12,12,12).

**Regression for the one-site limits** [exact] ([B:41]): for site 2 =e_1 every block has rank <=1 and P(c,c') equals half the rank of column c of the one-site tensor of the group c_2=c' (checked for Omega_c[w](x) e_1 and generic c_1(x) e_1 against the one-site tensors rebuilt from K^((c_2))); totals 5/5 resp. 5/7 are the general notebook's generic and vacuum-family counts. All four blocks are populated (slots 3,5,1,3 / 6,9,3,6), not only c'=lambda_2.

**det_6 regression** (plan S7, [B:19], [B:41]). Omega_c[w](x)Omega_c[w] is proportional to det_6[F_(vac,w)]/Psi_Q on all 64 / 225 states at 5.7*10^(-46), 4.5*10^(-44) / 4.4*10^(-43), 5.2*10^(-40) (for w, w'; vacuum direct q degrees (0,2,4) / (0,4,6)), with Psi_Q=M^((1))_vac M^((2))_vac to 10^(-46)/10^(-43) [numeric]. The plan's expectation that sum r over the four blocks is bounded by the 20 Laplace terms of the single det_6 is *not* a valid bound: each Laplace term feeds all four (c,c') blocks, so the sum is bounded by 80 (four blocks of at most 20 each; per block r<=min(20,n_(min))=n_(min), non-binding), and indeed sum r=26>20 for [3,1,0]. What is visible, and passes, is that the vacuum has strictly lower ranks and counts than a generic pair (13 vs 14, 26 vs 27; 18 vs 23, 40 vs 54; rank phi=12<15 for [3,1,0] versus 15 for the generic pair and for Omega_c[w](x)Omega_c[w']). The termwise comparison with the 20 Laplace terms was not done (the per-ladder 3x3 minors are not polynomials in the monomials without the Baxter/Plucker relations of the eigenstates). Nothing had to be fixed in the bookkeeping [exact/numeric].

**The decisive comparison** (plan S6, [B:41]). One-site counts: site 1 5 generic / 4 (adjoint quadric) / 5 ([3,1,0] vacuum family) / 7 ([3,1,0] generic); site 2 dressed 5/4 resp. 7/5, undressed 15/35. Products of one-site counts: genericxgeneric dressed 25 / 49, bestxbest 16 / 25, undressed 75 / 175. The count of Omega_c[w](x)Omega_c[w'] is 22 / 54: just below the generic pair (23) for the adjoint, *exactly* the generic pair's counts for [3,1,0], far above bestxbest. Two different vacua do not decouple at the bracket level, and the R-dressed candidates (h) are exactly generic [exact].

### 3.4 The four-parameter family Omega_c[w](x)Omega_c[w']

[exact] ([B:46]). In (a_1,s_1,a_2,s_2)=(1/w_1,1/w_2,1/w_1',1/w_2') the blocks are matrices of polynomials (entries of degree <=3 resp. <=5 in each variable; <=12 resp. <=22 terms); symbolic phi equals the numeric one. The w_2->inf limit is gaussCov{a,0,a}=e_1^Texp(aE_(12))exp(aE_(23)) (computed, and the guess was right). Generic ranks (40 of 40 random points): (5,3,3,2) / (9,6,6,6), sum P=22 / 54.

**Table 4.** Sub-families of Omega_c[w](x)Omega_c[w'] ([B:46]): sum r, sum P.

| sub-family | [2,1,0] sum r | sum P | [3,1,0] sum r | sum P | reach ([D:23]) |
|---|---|---|---|---|---|
| generic (4 parameters) | 13 | 22 | 27 | 54 | complete |
| w'=w (vacuum) | 13 | 18 | 26 | 40 | (not a two-rotation family) |
| common w_1'=w_1 (3 parameters) | 13 | **18** | 26 | **40** | 2 weights: no |
| w_2'=w_2; w_1'=w_2; w_2'=w_1; swap; w_1=w_2; w_1'=w_2' | 13 | 22 | 27 | 54 | generic |
| w_2->inf | 13 | 20 | 27 | 54 | -- |
| w_2'->inf | 10 | 14 | 24 | 45 | complete |
| w_2,w_2'->inf | 10 | 12 | 24 | 45 | complete (w_1!= w_1') |
| w_1'=w_1, w_2->inf | 12 | 14 | 24 | 27 | 2 weights: no |
| w_1'=w_1, w_2'->inf | 9 | **10** | 18 | **21** | 2 weights: no |
| w_1'=w_1, w_2,w_2'->inf (1 parameter) | 8 | **8** | 14 | **14** | global rotation |
| w_1,w_1'->inf (both e_1) | 4 | 4 | 4 | 4 | trivial |

The count drops on w_1'=w_1 and on the w->inf limits, not on the other coincidences; for [2,1,0] the coincidence w_1'=w_1 is not a matrix-rank drop (the P drop is a decomposability locus), for [3,1,0] the (0, 0) rank drops 9->8 there. No sub-family count is a product of one-site counts. Rank-drop hypersurfaces (common irreducible factors of det(PAQ) for random small-integer P,Q; Grobner bases not used, codimension >=2 probed only through the table): coordinate hyperplanes (a_1,a_2,s_2), a degree-14 hypersurface (163 terms) in block (0, 0) for the adjoint, and for [3,1,0] block (1, 1) the t-dependent ratio 9a_1=44a_2 (w_1'/w_1=44/9), 2134a_1^2-1188a_1a_2+279a_2^2, 210a_1^2-665a_1a_2-93s_2(a_1-a_2) and a degree-4 factor; block (0, 0) of [3,1,0] only through a random line section and not reliable.

### 3.5 Follow-up Task 4b: the (inf,inf) family with w_1!= w_1'

[exact] ([B:49]). Phi=gaussCov{a,0,a}(x)gaussCov{a',0,a'}, a=1/w_1, a'=1/w_1' (the limit is confirmed symbolically; numerically Omega_c[{3/2,10^6}] differs from gaussCov{2/3,0,2/3} by 10^(-6) relative); three pairs (2/3,-5/7),(3/2,1/4),(-4/5,7/3) give identical numbers; overlap identity with 12 / 45 explicit bracket products 10^(-44) / 10^(-39) (see the caution on the first two points in Appendix C, item 6).

**Table 5.** Task 4b ([B:49]); blocks (0,0),(0,1),(1,0),(1,1).

|  | [2,1,0] r | P; sum r/sum P | [3,1,0] r | P; sum r/sum P |
|---|---|---|---|---|
| **(inf,inf), a!= a'** | 3,3,2,2 | 4,4,2,2; 10 / **12** | 6,6,6,6 | 12,12,9,12; 24 / **45** |
| gaussCov{a,0,a}(x)Omega_c[w'] | 5,3,3,2 | 8,5,5,2; 13 / 20 | 9,6,6,6 | 18,12,12,12; 27 / 54 |
| Omega_c[w](x)gaussCov{a',0,a'} | 3,3,2,2 | 5,5,2,2; 10 / 14 | 6,6,6,6 | 12,12,9,12; 24 / 45 |
| vacuum Omega_c[w](x)Omega_c[w] | 5,3,3,2 | 10,3,3,2; 13 / 18 | 8,6,6,6 | 22,6,6,6; 26 / 40 |
| generic pair | 5,3,3,3 | 10,5,5,3; 14 / 23 | 9,6,6,6 | 18,12,12,12; 27 / 54 |
| global rotation a=a' | 2,2,2,2 | 2,2,2,2; 8 / 8 | 2,3,3,6 | 2,3,3,6; 14 / 14 |

On the (a,a') plane: generic ranks (3,3,2,2) / (6,6,6,6); sum P on lines (two members each): a=a' 8 / 14, a=0 8 / 19, a'=0 4 / 4 (trivial). Rank-drop components found: [2,1,0]: (0, 1): a, a-a', a'; (1, 0), (1, 1): a' (the drop of block (0, 0) from 3 to 2 on a=a' is visible in the line table but not detected by the common-factor heuristic); [3,1,0] blocks (0,0),(0,1): a, a-a', a'; (1, 0): also 31a-66a', 66a-31a', 53a-22a' and a quartic; (1, 1): also 9a-44a', 6a-19a', 2134a^2-1188aa'+279a'^2 and a quartic. **Verdict** ([B:49]): for [3,1,0] the family costs 45 against 54 (generic) and 40 (vacuum), lower bound 24, no block with P=r; not cheap, no wedge decomposition listed. For the adjoint it costs 12 (generic 23, vacuum 18), blocks (1,0),(1,1) exact; the explicit decomposition at (a,a')=(2/3,-5/7) is printed in the output cell of [B:47] (cell 48): block (0, 0): 3 terms with ladder-1 wedges 2+1+1, ladder-2 one wedge each; (0, 1) the same; (1, 0) and (1, 1) 2 terms. The 8/14 line is the global rotation a=a' only; the earlier Task 4 listing of w_1'=w_1, w_2,w_2'->inf as "1 parameter, minimal" ([B:46]) refers to the same line.

### 3.6 One-site-fixed slices and the type-pair scan

[exact] ([B:54]). *(A) site 2 fixed at Omega_c[w'], g_1=(g_A,g_B,g_C) free.* Blocks are polynomial in three variables; generic ranks (5,3,3,2) / (9,6,6,6) (the (1, 1) block already has the vacuum value 2 for every g_1). The site-1 vacuum family, b_3=0, and the adjoint-type conics cc_1(t),cc_1(-t),cc_1(1-t) all give 22 / 54, the same as generic g_1: the reduction 23->22 relative to the generic pair is due to the fixed vacuum on site 2, not to g_1. Reducing loci are cross loci depending on w': alpha=alpha' (common w_1; 20 / 49), with the vacuum condition 18 / 40, and beta=0 (20 / 54). *(B) site 1 fixed at Omega_c[w], g_2 free.* Generic 23 / 54; the site-2 vacuum family is a genuine rank-drop locus of block (1, 1) for the adjoint (3->2; count 22) but not for [3,1,0] (54); b_3=0 gives 22 / 50; beta=0 gives **16** / 48 (ranks 3,3,3,3 resp. 6,6,6,6); alpha=alpha_w and the vacuum family 18 / 40; the three conics 23 / 54. (B) is not the mirror of (A). Best count at fixed g_2: 18 / 40; at fixed g_1: 16 / 40. *Type-pair scan* over simple loci per site {generic, vacuum, beta=0, gamma=0, alpha=0, b_3=0, three conics}: minimum 14 for the adjoint at (beta_1=0,beta_2=0) and (alpha_1=0,beta_2=0), and 40 for [3,1,0] at (alpha_1=0,beta_2=0); CG-completeness of these loci was not tested. The coincidence loci of the table above (10/21, 8/14) are cheaper than any pair of independent simple loci.

**The sites do not decouple at the bracket level** (plan S3 item 5): the one-site best loci are not count-reducing when the other site is a vacuum, and the cross terms appear exactly as loci depending on w'.

### 3.7 Driver and conclusion of Stage B

[exact] ([B:57]). Rebuilt at (h,theta) points t=62/35 (also with z',z''), 29/15, 19/14, 1404/7: K identical at z',z''; the R solution space is 1-dimensional; cc_1(-t) is among the dressed conic constants; dressed/undressed generic counts 5/17 and 7/40; sum P for (generic, vacuum, Omega_c[w](x)Omega_c[w']) =(23,18,22) and (54,40,54), independent of (h,theta,z) at all points; 17-35 s per point for [3,1,0].

**Conclusion** [exact counts; reading for the interpretation] ([B:58], items 6, 7b, 12). (i) Omega_c[w](x)Omega_c[w'] is as expensive as a generic pair for [3,1,0] (54 = 54), so the plan's success target (~25) is not met by Omega_w(x)Omega_(w') itself; reductions come from coincidences (w_1'=w_1) and w->inf limits. (ii) [reading] Cheapness correlates with the pair becoming a near-global rotation (common w_1 differs only in E_(23), so it reads one irrep; w_2=w_2'=inf with w_1=w_1' is a global rotation and reads the top irrep only), so no cheap *and* CG-complete two-vacuum family exists in this cost metric -- *except* the (inf,inf) family with w_1!= w_1', which Stage D shows to be CG-complete (S4) and which costs 12 (adjoint, cheap) / 45 ([3,1,0], not cheap). The agent's initial Stage C recommendation (common-w_1 family, 18/40) was overruled by the reach argument and then by Stage D. (iii) The pivot count P is not monotone in r (vacuum [3,1,0] block (0, 0): r=8, P=22; generic r=9, P=18), so r is the robust number. (iv) Skipped: greedy/tensor-rank minimisation (done in Stage C for the (inf,inf) family), termwise Laplace comparison, Grobner bases for the rank-drop ideals, [3,2,0], reach tests (Stage D).

## 4. Stage D: the untwisting of <Omega_w(x)Omega_(w')| (chain side)

Notebook D (40 cells; summary [D:40]). Diagonal frame, 200 digits on the chain, factorised x_0 (never `Inverse` at small eps). The family is Phi=Omega_D[z,w](x)Omega_D[z,w']. For each highest-weight (HW) ket of V(x) V there is an HW-limit branch Psi_B of the twisted chain at z (normalised < x_0|Psi_B>=1 with x_0=(e_1^TP_( phys)^(-1))^((x)2)); P_(HW) is the G-orthogonal projector onto the HW space of the sector, R_C=<Phi|Psi_B>, S=<Phi|P_(HW)Psi_B>, K=R_C-S. A HW state of weight mu has its CG slice at the one-site weight nu=mu-Lambda (site 2 =e_1). The regime is z=1+eps*phi with phi=(1/3,-2/7,5/11) (control phi'=(-2/5,3/8,7/11)), w=(1+eps*om_1,1+eps*om_2), w'=(1+eps*om_1',1+eps*om_2'), 24 grid points om_1 in {3/4,-5/3}, om_2 in {-1/5,9/7,5/3}, om_1' in {-12/5,7/2}, om_2' in {1/8,-3/4}, eps in {1,2,4}*10^(-6); w_3,w_3' drop out of the covector ([D:1]).

### 4.1 Setup, frame tie, exact undo

[exact unless stated] ([D:15]). Gelfand-Tsetlin relations, invariant positive metric, quantum determinant (200 digits, 10^(-200)); HW kets/sectors 6/5 (multiplicities 1,1,1,2,1) for the adjoint and 10/8 (1,1,1,2,1,2,1,1) for [3,1,0]; slice weights nu=mu-Lambda: (2,1,0)_1,(2,0,1)_1,(1,2,0)_1,(1,1,1)_2,(0,1,2)_1 resp. (3,1,0)_1,(3,0,1)_1,(2,2,0)_1,(2,1,1)_2,(1,3,0)_1,(1,2,1)_2,(1,1,2)_2,(0,2,2)_1. (a) The frame tie holds at 6 rational points (2 with w_2=inf, 1 regime point) up to the scalar kappa (the spec's wording "equal" is corrected: Omega_c is normalised to first component 1). (b) The exact undo

    <Omega_D[w](x)Omega_D[w']|h>=< e_1^Trho(n(w)n(w')^(-1))(x) e_1^T|h> on all HW kets

holds at 7 points (2 with w_2=inf, 2 in the regime), with no scalar since e_1^Trho(n)e_1=1; the naive e_1(x) e_1 pairing differs. (c) Omega_D[w](x)Omega_D[w] is a left eigencovector of t_a(2/9) with the bra twist G'(w)=P_( aux)G(w)P_( aux)^(-1) (control diag(w) fails for t_1,t_2). Continuity at w_2=10^(10),10^(30) against w_2=inf: 2.2*10^(-10),2.2*10^(-30) (adjoint). At eps=10^(-6) there are 21 / 56 states and 6 / 10 branches, raising residual <=4.1*10^(-7) / 10^(-6) against >=1.3 for the other states, HW fraction 1-1.7*10^(-14) / 1-8.2*10^(-14), exact dual Baxter residual 1.1*10^(-147) / 2.4*10^(-140) [numeric].

### 4.2 The undone unipotent

[exact, rep-independent] ([D:19]). With n(w)=([[1,a,b],[,1,c],[,,1]]), a=r_1, c=r_3, b=r_2+r_1r_3=(w_1-z_1)/(w_1-z_3) *independent of w_2*, put n'=n(z,w)n(z,w')^(-1)=e^(a'E_(12))e^(g'E_(13))e^(c'E_(23)) with g'=b'-a'c', hence g'=(b-b_(w'))-(a-a_(w'))c: it depends on w' only through w_1'. All three are O(1) at eps=0 (orders 0,0,0), with (phi=(f_1,f_2,f_3), om=(o_1,o_2), om'=(p_1,p_2))

    a'_0 =((f_1-f_2)(o_1-p_1))/((o_1-f_2)(p_1-f_2)),  
g'_0=((f_1-f_3)(f_2-f_3)(o_1-p_1)(o_2-p_1))/((f_3-o_1)(f_3-o_2)(f_3-p_1)(p_1-f_2)) (independent of p_2),
    c'_0 =r_3^0(om)-r_3^0(om'),   r_3^0(om)=((o_1-f_2)(f_2-o_2)(f_1-f_3))/((f_1-f_2)(o_1-f_3)(f_3-o_2)),
(c'_0 has a long closed form with numerator N=f_2f_3(o_1+o_2-p_1-p_2)-(f_2+f_3)(o_1o_2-p_1p_2)+o_1o_2(p_1+p_2)-p_1p_2(o_1+o_2) in [D:19]). The map (om_1,om_2,om_1',om_2')|->(a'_0,c'_0,g'_0) has rank 3 (Jacobian at three generic points, one with phi'; one-dimensional fibre). The w_2->inf limits used in Stage C: g'_0=(f_1-f_3)(f_2-f_3)(o_1-p_1)/((f_3-o_1)(f_2-p_1)(f_3-p_1)) (the limit om_2->inf of the closed form above; [C:57]).

**Controls** [exact] ([D:19]): (i) w'=w: n'=1. (ii-a) the old regime's pinning w_2'=y_0(1+eps*om_2'), y_0=(t-1)/(t-2)=-27/8: orders (0,0,0), g_0' the same function (only c_0' loses om_2'); (ii-b) w_2'=7/4 fixed: same; (ii-c) w'=(7/4,-2/5) fixed entirely: orders (0,0,0), a'_0=(f_1-f_2)/(f_2-o_1), g'_0=-(f_1-f_3)(f_2-f_3)/((f_3-o_1)(f_3-o_2)), and (om_1,om_2)|->(a_0',c_0',g_0') has rank 2; (iii) common w_1'=w_1: a'== g'==0 exactly, n'=exp(c'E_(23)) with c_0'=(f_1-f_3)(f_2-f_3)(f_2-o_1)(o_2-p_2)/((f_1-f_2)(f_3-o_1)(f_3-o_2)(f_3-p_2)) (the localising-span result); (iv) common w_2'=w_2, w_1' untwisted: orders (0,0,0), c_0'~(o_1-p_1)(f_2-o_2). **Surprise** (spec expectation not realised): the spec expected g' to lose an order when w_2' is pinned (the old det_6xdet_2 degeneration). It does not: the old family had the second vacuum slaved to a resolvent of the first (y pinned, g'=a'c'(z_2-w_2)/(y-z_2) with a pinned denominator), whereas here the two vacua are independent and g' has no pinned denominator. The only way to lose a direction is w_1'=w_1.

### 4.3 Finite-twist reach

[exact] ([D:23]; three rational (z,w,w') points per family, z also varied). Rank of the components of phi_u=e_1^Trho(n(z,w)n(z,w')^(-1)) on each one-site weight space; "CG-complete" means every slice weight with full rank.

|  | (i) generic | (ii) w_1'=w_1 | (iii) w_2=w_2'=inf, w_1'!= w_1 | (iii-c) same with w_1'=w_1 | (iv) w_2'=inf only | (v) w_2'=w_2 |
|---|---|---|---|---|---|---|
| adjoint: sum ranks (of 8) | 8 | 2 | **8** | 1 | 8 | 8 |
| adjoint: HW sectors read (of 5) | 5 | 2 | 5 | 1 | 5 | 5 |
| [3,1,0]: sum ranks (of 15) | 15 | 2 | **15** | 1 | 15 | 15 |
| [3,1,0]: HW sectors read (of 8) | 8 | 2 | 8 | 1 | 8 | 8 |
| CG-complete | yes | no | **yes** | no | yes | yes |

Rank 2 on every double weight in the complete families. (ii) reads exactly e_1^T and e_1^TE_(23): the top irrep and one further irrep ([4,1,1] resp. [6,1,1]). (iii-c) is a global rotation (n'=1): top weight only. **(iii) with w_1'!= w_1 is not a global rotation: it reads every weight of V with full rank, with only two free parameters (w_1,w_1')** (in the inf limit r_2,r_3 and hence n(z,w) depend on w_1 only); the spec's expectation "top irrep only" holds only for the control (iii-c).

### 4.4 Regime, per-weight orders; K/S

[exact] ([D:27]; exact in symbolic eps at the 24 grid points, phi and phi' identical). In the main family every weight of V (7 resp. 12 slice-relevant weights) has its leading order at eps^0, the rank of the eps^0 coefficient vectors over the grid equals the multiplicity (1 on single weights; 2 on (1,1,1) for the adjoint, on (2,1,1),(1,2,1),(1,1,2) for [3,1,0]) and nothing new appears at eps^1. The same holds for the control w'=(7/4,-2/5) fixed (6 points), w_2=w_2'=inf and w_2'=inf (12 points each); for common w_1'=w_1 only (2,1,0),(2,0,1) [(3,1,0),(3,0,1)] are nonzero, all others identically zero. For comparison the old det_6xdet_2 family had the depth table eps^0,eps^1,eps^1,eps^2,eps^3,eps^3,eps^4 (adjoint). The control with w' fixed does not reproduce a depth table (spec expectation not realised): n(z,w)n_1^(-1) with n_1 the finite limit of n(z,w') is a generic unipotent.

**Signal versus contamination** [numeric, 200 digits] ([D:31]). *Decisive expectation: K/S=O(eps) for every HW branch -- confirmed* for both reps, including the adjoint singlet and both octets and all ten [3,1,0] branches ([4,4,0], [4,3,1]x2, [4,2,2], [3,3,2] included, which had eps-constant kappa in the old family); max|K/S| doubles with eps (ratios 2.000, 2.000) in every branch; the top branch has K=0 exactly; the exact undo holds to 10^(-157) (adjoint) / 10^(-149) ([3,1,0]) at 200 digits. The control w' fixed, w_2=w_2'=inf and w_2'=inf are equally clean; common w_1'=w_1 reads only the top and the first irrep (|S|<10^(-100)|R_C| elsewhere).

**Table 6.** K/S, main family (24 points) ([D:31]): eps-exponents of S and K, and max|K/S| at eps=10^(-6),2*10^(-6),4*10^(-6).

| [2,1,0] | S exp | K exp | max/K/S/ (10^(-6)) | x2 | [3,1,0] | S exp | K exp | max/K/S/ (10^(-6)) | x2 |
|---|---|---|---|---|---|---|---|---|---|
| [4,2,0] | -2 | -- (K=0) | 0 | -- | [6,2,0] | -2 | -- | 0 | -- |
| [4,1,1] | -3 | -2 | 4.25*10^(-7) | 2.000 | [6,1,1] | -3 | -2 | 4.25*10^(-7) | 2.000 |
| [3,3,0] | -3 | -2 | 2.1*10^(-6) | 2.000 | [5,3,0] | -3 | -2 | 3.12*10^(-7) | 2.000 |
| [3,2,1] (copy 1) | -4 | -3 | 1.5*10^(-5) | 2.000 | [5,2,1] (1) | -4 | -3 | 5.65*10^(-7) | 2.000 |
| [3,2,1] (copy 2) | -4 | -3 | 8.05*10^(-6) | 2.000 | [5,2,1] (2) | -4 | -3 | 2.0*10^(-5) | 2.000 |
| [2,2,2] | -6 | -5 | 2.57*10^(-6) | 2.000 | [4,4,0] | -4 | -3 | 4.21*10^(-6) | 2.000 |
|  |  |  |  |  | [4,3,1] (1) | -5 | -4 | 3.69*10^(-6) | 2.000 |
|  |  |  |  |  | [4,3,1] (2) | -5 | -4 | 1.09*10^(-5) | 2.000 |
|  |  |  |  |  | [4,2,2] | -6 | -5 | 4.51*10^(-6) | 2.000 |
|  |  |  |  |  | [3,3,2] | -7 | -6 | 4.44*10^(-6) | 2.000 |

With phi' at eps=10^(-6): <=1.6*10^(-5) (adjoint), <=4.4*10^(-5) ([3,1,0]). The largest values in the other families at eps=10^(-6) are 2.89*10^(-5) (w_2'=inf, [4,1,1] and [6,1,1]) and 6.0*10^(-6) (control w' fixed). The exponents of K are those of S plus one (K exp of [3,3,2] is printed as -6.000,-5.997 in the main family, -6.000 elsewhere).

### 4.5 CG data at leading order from the overlap alone

[numeric] ([D:35]). For 1-dimensional slices C:=R_C/phi_u(0)[slice] with phi_u(0) the exact eps^0 coefficient vector (verified equal to e_1^Trho(e^(a'_0E_(12))e^(g'_0E_(13))e^(c'_0E_(23))) with the closed forms at all 24 points); for 2-dimensional slices the 2x2 system from two grid points with independent phi_u(0)[slice] (12 pairs). Errors and spreads are O(eps) (error ratio 4.00 between 10^(-6) and 4*10^(-6) everywhere). Second-order Richardson on f=Ceps^p, p=-(eps-exponent of the slice), at eps=10^(-6),2*10^(-6),4*10^(-6): overlap-extracted limit and chain limit agree to 10^(-16)-10^(-20) ("<10^(-N)" = identical to the working precision for the 1-dimensional slices of the first irreps). **All CG data are read at leading order from the overlap: 6 of 6 HW data of adj(x)adj and 10 of 10 of [3,1,0]**, for the main family, the control w' fixed, and the sub-families w_2=w_2'=inf and w_2'=inf (common w_1'=w_1 reads only the first two irreps).

**Table 7.** CG slices, main family ([D:35]): chain Richardson limit of the slice and overlap-vs-chain agreement.

| [2,1,0] | dim | chain limit | agree | [3,1,0] | dim | chain limit | agree |
|---|---|---|---|---|---|---|---|
| [4,2,0] | 1 | 2.60946745562 | <10^(-191) | [6,2,0] | 1 | 2.60946745562 | <10^(-191) |
| [4,1,1] | 1 | 46.6745389114 | <10^(-184) | [6,1,1] | 1 | 46.6745389114 | <10^(-184) |
| [3,3,0] | 1 | 10.9285389167 | <10^(-184) | [5,3,0] | 1 | -36.8838188439 | <10^(-183) |
| [3,2,1] (1) | 2 | (248.327054711, -21.3321644927) | 1.2*10^(-18) | [5,2,1] (1) | 2 | (-835.970434496, 43.3173128073) | 1.4*10^(-19) |
| [3,2,1] (2) | 2 | (11.3988399453, -83.3191974257) | 7.3*10^(-18) | [5,2,1] (2) | 2 | (22.7423320599, -47.4770351603) | 9.7*10^(-19) |
| [2,2,2] | 1 | -5118.53831440 | 1.5*10^(-19) | [4,4,0] | 1 | -115.853020728 | 1.7*10^(-19) |
|  |  |  |  | [4,3,1] (1) | 2 | (-3129.70920645, 247.139873999) | 3.7*10^(-19) |
|  |  |  |  | [4,3,1] (2) | 2 | (-65.0505836460, 570.310701774) | 4.5*10^(-17) |
|  |  |  |  | [4,2,2] | 2 | (4647.35409165, -1549.11803055) | 3.3*10^(-18) |
|  |  |  |  | [3,3,2] | 1 | 23369.8056448 | 7.0*10^(-19) |

The slice normalisation is < x_0|Psi>=1, so chain limits are family- and rep-independent for shared irreps, and they reproduce the earlier notebooks' numbers: adjoint [4,1,1] 46.67453891, [3,3,0] 10.92853892, singlet -5118.538314, octet functionals (u+u')* C=C_1+3C_2 =184.3305612 and -238.5587523, octet mixing R=u* C/((u+u')* C)=-1.6944, 1.0956 (earlier -1.694, 1.096); [3,1,0]: 46.6745, -36.8838, u'* C=2C_1+8C_2 =-1325.402367 and -334.3316172 (earlier -1325.40..., -334.33...). New [3,1,0] slices: [4,4,0] -115.8530207, [3,3,2] 23369.80564, [4,3,1] (-3129.709206, 247.1398740) and (-65.05058365, 570.3107018), [4,2,2] (4647.354092, -1549.118031), [5,2,1] (-835.9704345, 43.31731281) and (22.74233206, -47.47703516) (GT-index components of the slice).

### 4.6 The four overlaps and the exchange parity

[numeric] ([D:39]; one grid point (om,om')=((3/4,-1/5),(-12/5,1/8)), eps=10^(-6), exponents from 10^(-6)->2*10^(-6); w_3=1+5eps/13, w_3'=1-4eps/13). Vacuum direct Q's of w,w' have degrees (0,2,4) / (0,4,6) with exact Baxter residuals 4.5*10^(-191),2.5*10^(-190) / 7.4*10^(-190),9.2*10^(-190); dual Q's of all HW-sector states of degree mu_j (21 resp. 56 states). The diagonal overlaps are det_6: <Omega_w(x)Omega_w|Psi>=c_wdet_6[B,F_vac(w)]/Psi_Q with a state-independent constant, deviation 3.2*10^(-116),1.8*10^(-117) [2,1,0] and 6.4*10^(-77),1.5*10^(-85) [3,1,0] (limited by the ~100-digit Q-determinants) [numeric]. The mixed overlaps are not det_6 and are the Stage C targets.

| adjoint, main family | < ww/Psi> | < w'w'/Psi> | < ww'/Psi> | < w'w/Psi> | exponents (ww,w'w',ww',w'w) | parity |
|---|---|---|---|---|---|---|
| [4,2,0] | 2.609e12 | 2.609e12 | 2.609e12 | 2.609e12 | -2,-2,-2,-2 | + |
| [4,1,1] | -1.198e12 | -2.409e12 | -4.241e18 | 4.241e18 | -2,-2,-3,-3 | - |
| [3,3,0] | 1.05e12 | 3.374e12 | -9.732e18 | 9.732e18 | -2,-2,-3,-3 | - |
| [3,2,1] (1) | 4.367e18 | 5.889e18 | -1.515e26 | 2.066e26 | -3,-3,-4,-4 | (2-dim) |
| [3,2,1] (2) | -1.153e19 | -1.555e19 | -1.585e26 | 1.411e26 | -3,-3,-4,-4 | (2-dim) |
| [2,2,2] | -4.533e25 | -8.244e25 | 5.02e39 | 5.02e39 | -4,-4,-6,-6 | + |

The diagonal overlaps scale like eps^(-2) (top) to eps^(-4) ([2,2,2], [4,2,2], [3,3,2]), the mixed ones carry the CG-slice exponent (-2 top, -3 for [4,1,1],[3,3,0],[6,1,1],[5,3,0], -4 for the octets and [4,4,0], -5 for [4,3,1], -6 for [2,2,2],[4,2,2], -7 for [3,3,2]). For the one-dimensional slices <Omega_w(x)Omega_(w')|Psi>=+-<Omega_(w')(x)Omega_w|Psi> to the leading digits, the sign being the exchange parity of the irrep: - for [4,1,1],[3,3,0],[6,1,1],[5,3,0],[3,3,2]; + for the top irrep, [2,2,2],[4,4,0],[4,2,2] ([4,2,2] is a 2-dimensional slice with one copy, symmetric); for the two-copy 2-dimensional slices (octets, [5,2,1], [4,3,1]) the two mixed overlaps differ. For w_2=w_2'=inf (also tabulated in [D:39]) the diagonal exponents are all -2 (the inf-vacua are top-heavy) and the mixed exponents are again the slice exponents.

**Surprises recorded in [D:40]**: (1) no untwisting of w' is needed, pinned and fixed w' are as clean as the full regime, the family needs only w_1'!= w_1; (2) the sub-family w_2=w_2'=inf (with w_1!= w_1') is CG-complete at finite twist, at eps^0 in the regime, in K/S and in the CG extraction; [D:40] attributes to it the Stage B block counts 8/14 ("cheap") and concludes that the HANDOFF reading "no cheap and CG-complete two-vacuum family" is contradicted at the chain level -- but per [B:49] the 8/14 line is the global rotation a=a', and the family with w_1!= w_1' costs 12 / 45 (Appendix C, item 2); only two direction parameters (om_1,om_1') remain, and the O(eps) layer is open; (3) the frame tie holds up to kappa, and the eigen-covector twist is G'(w), not diag(w); (4) undo accuracy 10^(-157)-10^(-192) ([3,3,2] 10^(-149)) at 200 digits. Not done in D: the O(eps) layer, and the Q side (Stage C).

## 5. Stage C: the adjoint (inf,inf) family as a bracket formula, and the Q-side CG

Notebook C (70 cells; summary [C:70]); adjoint [2,1,0] only, plus a [3,1,0] cost note. Tasks in the order 0,1,3,2,4 because the Task 2 and Task 4 drivers wipe `Global`` ([C:1]). The family is Phi(a,a')=gaussCov{a,0,a}(x)gaussCov{a',0,a'}=Omega_c[(1/a,inf)](x)Omega_c[(1/a',inf)] (first rows (1,a,a^2), (1,a',a'^2)); diagonal frame Omega_((w_1,inf))(x)Omega_((w_1',inf)) with w_1=1/a, w_1'=1/a'. The identity to realise is

    (<Phi|Psi>)/(< x_0|Psi>) M^((1))_vac M^((2))_vac=sum_(c,c' in {0,1})Q^3(theta_1-ch) Q^3(theta_2-c'h)sum_kBk2^((1))[F_k,G_k] Bk2^((2))[F'_k,G'_k],

M^((s))_vac=minor_s(-3,-2)Q^3(theta_s-h), ladder points m=-3,...,0.

### 5.1 Setup and the count: ten, not twelve

Part A (60 digits) reproduces Stage B's setup (64 states, largest sector 10, eigen residual 10^(-52), dual Baxter 10^(-49), wavefunction rule 10^(-43), K: 278 nonzero rationals, rank 64; bracket identity verified symbolically) [exact/numeric] ([C:21]). The four blocks A^((cc'))(a,a') (6x6) are matrices of polynomials of total degree <=3 with 11,11,9,9 nonzero slots; generic ranks (3,3,2,2), sum r=10; at (a,a')=(2/3,-5/7),(3/2,1/4),(-4/5,7/3) the default pivots give P=(4,4,2,2), sum P=12 (rows/columns orientation 4,4,2,2 / 5,5,2,2), reproducing Stage B exactly; the identity with the 12 explicit products on all 64 states is <10^(-44) [exact/numeric] ([C:24]).

**Resolution of 12 versus ten** [exact] ([C:27], [C:32]). The 12 is the count of the *default* `RowReduce` pivot choice, in which one of the three terms in blocks (0, 0) and (0, 1) has a ladder-1 factor of rank 4 (two wedges); the "3+3+2+2=10" of Stage B's wedge list is the number of rank-one *terms* u_k(x) w_k, i.e. the lower bound sum r. The scan of *all* independent pivot sets in both orientations at the generic point (a,a')=(7/5,-4/3) gives: rank-3 blocks, rows orientation one admissible set (P=4), columns orientation 7 sets with P=3 (twice), 4 (four times), 5 (once); rank-2 blocks, all 3 sets P=2 in both orientations. The minimum is 3,3,2,2: **ten products =sum r, attained exactly**, because for the rank-3 blocks every one of the 3+3 2-forms is decomposable; the Pfaffians of all 20 ladder 2-forms vanish identically in (a,a') (the code returns `$Failed` otherwise and reconstructs the blocks exactly from the products). Winning choices: blocks (0,0),(0,1): columns orientation, pivot rows {2,3,4}; (1, 0): rows orientation, pivot columns {1,3}; (1, 1): rows, {2,3}. The identity of the ten products holds on all 64 states at (5/7,-4/9),(-3/5,7/4),(9/5,2/7) to <10^(-45) with all 64 overlaps nonzero ([C:32]). *Caution*: (a,a')=(2/3,-5/7) has w_1=1/a=z_1 and (3/2,1/4) has w_1=z_3 -- degenerate vacua for which 22 of the 64 states are orthogonal to the covector (42 nonzero overlaps); the identity holds there trivially on those states, and non-vacuity is claimed only at the generic points ([C:24], [C:70]).

### 5.2 Closed forms, the scaling law and where a-a' enters

[exact] ([C:32]; printed in [C:29] (s-frame, t=62/35), [C:31] (raw rational functions of (a,a'), last row), [C:62] (exact in (s,t))). Every block entry is homogeneous in (a,a') and

    A^((cc'))_((mm'),(nn'))(a,a')=a^( m+m'+e) a'^( n+n'+e') kappa^((cc'))_((mm'),(nn'))(s),   s=a'/a,

with kappa a polynomial in s and (e,e')=(7,5),(7,4),(6,5),(6,4) for (c,c')=(0,0),(0,1),(1,0),(1,1) (e+e'=12,11,11,10; verified in [C:29]). So on ladder 1 every wedge vector is F_m=a^mf_m(s), G_m=a^(m+e)g_m(s) (the bra twist a^((u-theta_1)/h)), on ladder 2 the same with a' and (F'_n,G'_n), and **the two ladders couple only through the single variable s=a'/a**; the whole (a,a') plane is a one-parameter family of kernels. (This is torus covariance: Phi(a)=e_1^Trho(e^(aE_(12))e^(aE_(23))) is Phi(1) conjugated by the diagonal torus.) Gauge: pivot (i,j) a nonzero slot of the antisymmetric matrix, F_i=0, F_j=-1, G_j=0. Example, product 1 of block (0, 0) (s-frame, t=62/35, [C:29]; m=-3,...,0):

    f =(0, -((s-1)(349071s-169136))/(9s(11579s-8184)), -1, 0),   g =(1,0,0,0),
f' =((35(s-1)s)/264,0,-1,0),   g' =(-(s(35s-132))/97, -132s/31,0,0),
which in the raw variables is G=(a^4,0,0,0), F=(0,((169136a-349071a')(a-a'))/(9a^2(8184a-11579a')a'),-1/a,0) ([C:32]; re-derived here that the two agree). *Structure* ([C:32]): of the 160 vector entries of the ten products, 92 vanish identically, 33 are constants, 31 are monomials c s^k (jointly monomials in a,a'), and only 4 are non-monomial in s (f'_(-2) of product 4, f'_(-3) of product 5, g_(-2) of product 6, f_(-2) of product 8); many vectors are delta functions; a-a'=a(1-s) appears through (s-1) factors, the pole 1/(s-1) occurring once (f'_(-2) of product 7, 132/(97(s-1)) at t=62/35); the other poles in s sit at roots of pivot minors. *Rank-drop loci* (gcd of all rx r minors of the kernels): (0, 0): (s-1)s^2; (0, 1): (s-1)s^3; (1,0),(1,1): s. On the lines: a=a' (s=1) ranks (2,2,2,2), P=8 (global rotation, Stage B's 8); a'=0 (s=0) ranks (1,1,1,1), P=4; a=0 (s=inf) ranks (2,2,2,2), P=8; generic (3,3,2,2), P=10 with optimal pivots.

### 5.3 t-dependence

[fitted; exact at all sampled points] ([C:60], [C:62], [C:65]). The kernels kappa^((cc'))(s) at three points with t=62/35 (base; (h,theta)->(2h,2theta)=(1,(4/7,-6/5)); other z=(5/3,-3/7,4/5)) are *identical*, as are (e,e'); K is z-independent and, in the canonical column order of the label pairs, invariant under (h,theta)->(2h,2theta): the entry ratio is a single value per column equal to 1 for all 64 columns, i.e. homogeneity degree 0 in the normalised-wavefunction units (monic Q's): K depends on (h,theta) only through t. The scaling law holds exactly at all 15 sampled points and the minimal-P closed forms exist at every point (blocks reconstructed exactly with the base pivots and gauge forced). The 15 points have 13 distinct values of t (62/35 three times, 29/15, 6/35, 16/15, 5/2, 7/3, 3/5, 10/7, 82/21, -2/21, -16/21, 57/14, 26/21). Every s-coefficient c_q(t) of every kernel entry kappa_(ij)(s,t)=sum_qc_q(t)s^q is fitted as a rational function of t of lowest total degree through 11 fit values (one null vector, 0 failures); the held-out t=57/14,26/21 and the two duplicates are reproduced exactly. The 40 kernel laws are printed in [C:62]; e.g.

    kappa^((1,1))_((-3,-2),(-3,-2)) =s,   kappa^((1,1))_((-3,-2),(-3,-1)) =-(2s(t+2))/t,
kappa^((1,1))_((-3,-1),(-3,-1)) =(4s(t-2)(t+2))/((t-1)(t+1)),   kappa^((0,0))_((-3,-2),(-3,-2)) =((s-1)(2s-1))/(t-1).
**Poles: the kernel coefficients have poles only at t=0,1,-1** (denominators t,t^2,t-1,t+1 and products; numerators contain t+-2, t+-3, 2t+1, 3t+1, 2t-3, 2t^2-2t-3, t^2-8t+6, 4t^2+t+1); there is no pole at the integers 2, 3 where the ladders meet -- those enter only the pivot gauge: the closed forms F,G,F',G' in (s,t) have denominators in t in {-2,-1,0,1,2}, the s-only factors s,s-1,s-2 and four mixed (s,t) polynomials (vanishing curves of pivot minors; gauge artefacts). Compare the site-2 locus b_3=0 of Stage A, whose kernels have poles at {0,+-1} (adjoint) and {0,+-1,+-2} ([3,1,0]) -- different objects.

### 5.4 The u-form: only trivial ones

[exact] ([C:64], [C:65]). Exact class tests on the ladder planes P=span{f,g} of the ten products (20 planes). Two planes have a two-point support (product 3 and product 6 on ladder 2: the plane is all functions on {-3,-2}; F,G are delta functions, Bk2 is a Q-minor at two points); 18 planes have a three-point support. For the 18 the dimension of the intersection with each class ((a) wgt_S(m)x polynomial; (b) the global six-pole SoV measure x polynomial; (c) polynomial; (d) cross-ladder poles x polynomial; (e) (d)xwgt_S; degrees 0,1,2) is always the generic dimension d (0,1,2) at three generic (s,t): the plane is never contained in a natural class (dimension 2 at degree 1 would be needed); the single exception is product 10, ladder 1, whose plane contains the constant vector (class (c), degree 0). Each plane contains exactly one natural vector per class, measure(m)*(m-m^*), with root m^*_( class)(s,t) an explicit rational function of (s,t) of high degree; irrational roots do not occur (the obstruction is dimensional, not arithmetic); 8 of the 18 planes contain a delta function. So the ladder kernels are a^((u-theta_1)/h) (resp. a'^((u-theta_2)/h)) times coefficients rational in (s,t) with delta- or two-point-supported wedge vectors, and there is no explicit F(u),G(u) beyond the trivial cases. Mobius/natural-measure tests beyond the dimension counts were not run.

### 5.5 The Q-number at small eps: dictionary and verification

Part B (200 digits) reuses Stage D's machinery verbatim (frame tie, exact undo, HW-limit branches with dual Q's, closed forms of n': [C:34]-[C:47], summary [C:48]); regime z=1+eps*phi, phi=(1/3,-2/7,5/11), w_1=1+eps*om_1, w_1'=1+eps*om_1', w_2=w_2'=inf, a=1/w_1, a'=1/w_1'; direction grid om_1 in {3/4,-5/3}, om_1' in {-12/5,7/2}, eps in {1,2,4}*10^(-6) ([C:33]). **Dictionary** [numeric, 4.9*10^(-123)] ([C:50], [C:57]): <Phi_c|Psi_c>/< x_0^c|Psi_c>=<Phi_d|psi_d>/(kappa*kappa') with kappa=kappa_( Tie)(z,w), kappa'=kappa_( Tie)(z,w'), Psi_c=(P^(-1)(x) P^(-1))psi_d and < x_0^c|Psi_c>=< x_0^d|psi_d>=1; at (om_1,om_1',eps)=(3/4,-12/5,10^(-6)): kappa=-1.79*10^(25), kappa'=9.38*10^(22) (~eps^(-4) each), all 21 HW-sector states. **Verification** [numeric] ([C:50], [C:57]): Q_( num)(B)=sum_(cc')Q^3Q^3sum_kBk2^((1))Bk2^((2))/(M^((1))_vac M^((2))_vac) with the ten closed-form products at the exact rationals a=1/(1+eps*om_1), a'=1/(1+eps*om_1') and the 200-digit dual Q's of the branch equals <Phi_c|Psi_c>/< x_0^c|Psi_c> on all 21 states with maximum relative deviation **7.9*10^(-111)** (non-vacuity: |Q_( num)| between 1.3*10^(-36) and 4.7*10^(-10); the deviation is limited by the ~160-digit Q's and the ~37-digit cancellation). At every one of four grid points, three eps and six branches (72 checks, [C:54]) the full R_Q=kappa*kappa'Q_( num) equals the chain overlap <Omega_D[w](x)Omega_D[w']|psi_B> to 1.6*10^(-113).

### 5.6 eps-exponents: where the cancellations sit

[exact entries; numeric exponents] ([C:52], [C:57]; point (om_1,om_1')=(3/4,-12/5); entries exact, the rest log_2 ratios between eps=10^(-6) and 2*10^(-6), integer to <=8*10^(-4), confirmed by 2*10^(-6)->4*10^(-6)).

| branch mu | brackets < Q^jF> | Bk2's | product terms | block sums | total Q_( num) | digits lost (in worst Bk2 / in the sum) |
|---|---|---|---|---|---|---|
| (4,2,0) | -3...1 | 0...2 | 0,1 | 0,0,0,0 | +6 | 19.9 / 36.7 |
| (4,1,1) | -1...1 | 0...2 | 0,1 | 0,0,0,0 | +5 | 6.7 / 29.5 |
| (3,3,0) | -4...1 | 0...2 | 0,1 | 0,0,0,0 | +5 | 26.2 / 30.9 |
| (3,2,1), copy 1 | -2...1 | 0...2 | 0,1 | 0,0,0,0 | +4 | 13.2 / 23.3 |
| (3,2,1), copy 2 | -2...1 | 0...2 | 0,1 | 0,0,0,0 | +4 | 14.2 / 23.4 |
| (2,2,2) | -1...1 | 0...2 | 0,1 | 0,0,0,0 | +2 | 6.8 / 10.7 |

Entries of F,G,F',G' have orders -1...1 (only f'_(-2) of product 7 has eps^(-1), from 1/(s-1); eps^1 entries come from (s-1) factors). Of the ten product terms, products 3 and 6 carry eps^1 and the other eight eps^0. The brackets are as large as eps^(-3)-eps^(-4) (the dual Q's of the HW-limit branches are large at the near-degenerate twist) but the 2x2 determinants are O(1)... O(eps^2): the first cancellation sits *inside the Bk2's* (up to 26 digits at eps=10^(-6)); the second in the *final sum over the four blocks*: every product and every block sum is O(eps^0), the total is O(eps^6...eps^2) (36.7, 29.5, 30.9, 23.3, 23.4, 10.7 digits). The total exponent equals the CG-slice exponent: R_Q=kappa*kappa'Q_( num)~eps^(-8+total)=eps^(-2),eps^(-3),eps^(-3),eps^(-4),eps^(-4),eps^(-6), exactly the mixed-overlap exponents of the Stage D table (S4.6). The brackets are not small; the whole eps-dependence of the CG datum is produced by cancellations between O(1) terms, and the 200-digit precision is what makes the formula usable at eps=10^(-6) (about 40 digits of headroom needed).

### 5.7 All six CG data from the Q side alone

[numeric] ([C:54], [C:56], [C:57]; 4 grid points, eps=1,2,4*10^(-6); the chain slice is used only as comparison). phi_u(0)=e_1^Trho(e^(a'_0E_(12))e^(g'_0E_(13))e^(c'_0E_(23))) with the w_2,w_2'->inf limits of the Stage D closed forms, equal to the direct leading Gauss parameters of n(w)n(w')^(-1) at w_2=w_2'=inf (exact, symbolic) and to the eps^0 series coefficient at all four points; every one-site weight is at order eps^0 with rank equal to the multiplicity.

**Table 8.** CG data from R_Q ([C:57]).

| mu | slice dim | spread (10^(-6),4*10^(-6)) | error vs chain (1,2,4x10^(-6)) | Richardson limit from R_Q only | Q vs chain | reference |
|---|---|---|---|---|---|---|
| (4,2,0) | 1 | 0 | 0 | 2.609467455621301775 | 0 | (none earlier) |
| (4,1,1) | 1 | 1.68*10^(-6), 6.71*10^(-6) | 2.22,4.43,8.86*10^(-6) | 46.67453891138101664 | 0 (10^(-131)) | 46.67453891 |
| (3,3,0) | 1 | 1.61*10^(-6), 6.44*10^(-6) | 2.10,4.21,8.42*10^(-6) | 10.92853891670459718 | 0 (10^(-112)) | 10.92853892 |
| (3,2,1) copy 1 | 2 | 5.02*10^(-6), 2.01*10^(-5) | 5.31,10.6,21.2*10^(-6) | functional (u+u')* C=184.33056123341339 | 1.2*10^(-18) | 184.3305612 |
| (3,2,1) copy 2 | 2 | 1.44*10^(-5), 5.75*10^(-5) | 1.77,3.54,7.08*10^(-5) | functional -238.55875233168809 | 5.2*10^(-19) | -238.5587523 |
| (2,2,2) | 1 | 2.24*10^(-6), 8.95*10^(-6) | 2.27,4.53,9.07*10^(-6) | -5118.5383144047475541 | 7.3*10^(-19) | -5118.538314 |

All errors and spreads scale as eps (ratio 4.00 between 10^(-6) and 4*10^(-6)); the second-order Richardson limits agree with the chain limits to <=1.2*10^(-18) (to working precision for the 1-dimensional slices of the first irreps); the normalisation is Stage D's (< x_0|psi>=1), no rescaling is needed. **Octet mixing** ([C:55], [C:56]): R_c=u* C/((u+u')* C) with u=e_1E_(13), u'=e_1E_(23)E_(12) on the (1,1,1) slice, (u+u')=(1,3): -1.694366610178479773 and 1.095564214968898936 from the Q side, equal to the chain values to 5*10^(-19) / 4*10^(-20) and to the earlier -1.694 / 1.096.

### 5.8 Cost statement

[C:56], [C:70] (item 7). Ten products = twenty 2x2 bracket determinants Bk2 and four block sums (c,c'). The ladder vectors have 18 (ladder 1) + 15 (ladder 2) distinct projective classes among the 2x20 vectors, so 2x33=66 brackets < Q^2F>,< Q^3F> (each a four-term sum; 80 without sharing); the Q-values needed are Q^2,Q^3 on the four points of each ladder, 16 numbers (the Q^3(theta_s-ch) and M_vac use the same points). Comparison: det_6 is 36 bracket entries (each an 8-pole sum) and one 6x6 determinant; the direct sum over states is dim V^2=64. The count is fixed by the ladder length and the block ranks, not by dim V; it is not a determinant, and for the adjoint it is of the same order as det_6 (66 four-term sums plus 20 2x2 determinants against 36 eight-term sums plus one 6x6 determinant); its content is the closed-form (a,a')-dependence of a CG-complete family.

### 5.9 The [3,1,0] note

[exact, one point] ([C:68], [C:69]). Regression of Stage B reproduced on the symbolic blocks: ranks (6,6,6,6), default pivots give 45 at three points with the explicit-bracket identity on all 225 states to 10^(-38)-10^(-39). At (a,a')=(7/5,-4/3): pivot-set scan: blocks (0,0),(0,1): rows orientation one admissible set (P=12), columns orientation 84 sets with P between 17 and 24 (best 17), so 12; block (1, 0): 9 (rows, single set; columns 12); block (1, 1): 12 in both orientations; sum of the best =mathbf45 (the Stage B default). Greedy elimination (cheapest term first, 400 random restarts per block): best 12,12,9,15, i.e. **48** (the greedy never beats the scan: (0, 0) and (0, 1) are 12 in all 400 runs, (1, 0) between 9 and 11, (1, 1) between 15 and 21). *The lower bound 24 is not attained by any searched decomposition*: every decomposition found needs on average one rank-4 (two-wedge) factor per term; best total 45 (det_6 40, generic 54). The searched family is pivot decompositions and greedy Schur-complement eliminations; a general decomposition with decomposable factors (the column space meets G(2,5) in a 2-dimensional family, but six independent decomposable columns with decomposable duals would be needed) was not excluded by proof [exact for the searched family; minimality unproved].

## 6. What it means and what is open

This section is short and factual. Statements marked [reading] are the main session's interpretations (from `HANDOFF.md`); the others restate notebook results.

1. **The cost metric (Paul's decision).** In the metric "number of products of two 2x2 bracket determinants" the best CG-complete two-vacuum family costs 10 (adjoint; =sum r, [C:32]) and, for [3,1,0], 45 (searched; lower bound 24; [C:69]), against the generic pair 23 / 54, the det_6 vacuum 18 / 40 ([B:41]) and the plan's target ~25. The det_6 itself is 36 bracket entries (each an 8-pole sum) and one 6x6 determinant ([C:56]), so the two metrics (products of 2x2 determinants versus determinants of brackets) are not commensurable: "det_6=40 products" counts only what the product-of-2x2 formula would need to reproduce it. Whether to (a) accept the generic block formula (polynomial in the boxes by construction: (lambda_2+1)^2 blocks of size C(lambda_1+2,2), 54 products for [3,1,0]), (b) change the metric to bracket determinants of size >2, or (c) search further (the minimal count is a tensor-rank problem; lower bound 24 for [3,1,0] unattained by all searched decompositions) is not decided by the notebooks.

2. **A determinant form of <Omega_w(x)Omega_(w')|Psi> with non-vacuum bra data.** The mixed overlaps are not det_6 ([D:39]); the diagonal ones are ([B:19], [D:39]). [reading] The exact two-site K tensor ([B:24]), z-independent and sparse, makes a linear-algebra search for a determinant form (bra data F on the 6 poles, insertions, rows/columns of the det_6 pattern) possible; not attempted.

3. **R-dressed two-site covectors.** The R-dressed candidates tested, (e_1(x)Omega_c[w'])R(rho(g_1)(x)1) and the other order, are exactly generic (23 / 54, [B:41]); the R-dressed (inf,inf) pair and dressings on both sites were not tested.

4. **cc_1(1-t) versus cc_1(-t).** The undressed site-2 adjoint conic constant is cc_1(1-t)=(t+2)/(2(t+1)) ([A:24], [fitted], 15 points; same function for [3,2,0], [A:35]); the R-dressed one is cc_1(-t)=(t+3)/(2(t+2)) ([B:36], checked at four (h,theta) points). [reading] The main session reads the difference as the R shift; no notebook derives it, and the shift of the argument by one unit of t is unexplained by the cells.

5. **The [3,1,0] negative.** No family found is both CG-complete and near ~25: CG-complete ones cost 45-54, the cheap ones (common w_1: 18/40; beta_2=0: 16/48; w_1'=w_1 with w_2'->inf: 10/21; the type-pair (alpha_1=0,beta_2=0): 40) lose reach or were not tested for reach ([B:54], [D:23]). In particular the type-pair loci and the line beta_2=0 were never tested for CG-completeness.

6. **The s-dependence.** The (inf,inf) formula is a one-variable (s=a'/a) deformation of the global rotation ([C:32]); whether kappa(s,t) is a known object (for instance a ratio of the two one-site vacuum Q's) is open ([reading]: raised by the main session; nothing in the cells).

7. **Not done:** the O(eps) layer of the Q-number; a minimality proof for [3,1,0]; Mobius/natural-measure tests ([C:70]); re-running the four notebooks in VS Code to store native outputs; [3,2,0] in Stages B-D.

## Appendix A Notebook and cell maps

Absolute cell positions; "p-q, r" means code cells p-q with their result cell r (in B and C the printed output of each code cell is the markdown cell right below it, listed only in the ranges).

### A.1 `L2-site2-rotation-brackets.wb` (A, 36 cells)

| **cells** | **content** |
|---|---|
| 1 | title, setup, one switch (`Lambda` line of cell 2), notation |
| 2-7, 8 | setup (generators, eigenstates, dual Q's, P_( phys)/x_0, B, wavefunction rule); result 8 |
| 9-11, 12 | site-1 regression (supports, tensors, counts 5/7, reduced 4/5); result 12 |
| 13-14, 15 | vacuum-family identification, lambda_2-group measures; result 15 |
| 16-17, 18 | Task 1: the mirror, support, weight rule; result 18 |
| 19, 20 | Task 2: ladder-2 tensors, counts (per-label 17/40, block form); result 20 |
| 21-23, 24 | Task 3: loci, joint loci, reach; result 24 |
| 25-28, 29 | Task 4: closed forms on b_3=0, verification, u-structure; result 29 |
| 30, 31 | Task 5: torus statement, section-5 comparison; result 31 |
| 32-33, 34 | driver (z-independence, (h,theta)-dependence; wipes `Global``) and analysis; result 34 |
| 35 | [3,2,0] (optional) |
| 36 | summary |

### A.2 `L2-two-site-rotation-brackets.wb` (B, 58 cells)

| **cells** | **content** |
|---|---|
| 1 | title, notation, P and r |
| 2-13, 14 | setup (generators, eigenstates, Q's, kets, B, wavefunction rule; printed outputs in 3,5,...,13); result 14 |
| 15-18, 19 | Task 0: Omega_c[w], det_6 regression; result 19 |
| 20-23, 24 | Task 1: two-site K, slice regressions; result 24 |
| 25 | tools (2-forms, wedges, pivot decompositions) |
| 27, 29 | Task 2: supports; result 29 |
| 30-35, 36 | Task 2b: R-matrix, dressed test, dressed loci; result 36 |
| 37-40, 41 | Task 3: block ranks and counts, det_6 regression, decisive comparison; result 41 |
| 42-45, 46 | Task 4: four-parameter family, sub-families, rank-drop loci; result 46 |
| 47-48, 49 | Task 4b: (inf,inf) family; result 49 |
| 50-53, 54 | Task 5: one-site-fixed slices, type-pair scan; result 54 |
| 55-56, 57 | driver (wipes `Global``); result 57 |
| 58 | summary |

### A.3 `L2-two-vacua-untwist.wb` (D, 40 cells)

| **cells** | **content** |
|---|---|
| 1-2 | title, notation, regime, cell map |
| 3-14, 15 | setup, group data, vacua, HW kets, Task 0 checks, Q machinery and branches; result 15 |
| 16-18, 19 | Task 1: undone unipotent, closed forms, controls; result 19 |
| 20-22, 23 | Task 1b: finite-twist reach; result 23 |
| 24-26, 27 | Task 2: per-weight eps-orders; result 27 |
| 28-30, 31 | Task 3: signal versus contamination, K/S; result 31 |
| 32-34, 35 | Task 4: CG at leading order; result 35 |
| 36-38, 39 | Task 5: det_6-type data and the four overlaps; result 39 |
| 40 | summary |

### A.4 `L2-210-two-vacua-brackets.wb` (C, 70 cells)

| **cells** | **content** |
|---|---|
| 1-2 | title, task map, notebook state |
| 3-20, 21 | Part A setup (Stage B cells + tools); result 21 |
| 22-23, 24 | Task 0: blocks in (a,a'), ranks, count 12, identity; result 24 |
| 25-31, 32 | Task 1: pivot scan (26-27), closed forms and scaling law (28-29), numeric identity and structure (30-31); result 32 |
| 33-47, 48 | Part B setup (Stage D machinery, closed forms of n'); result 48 |
| 49-56, 57 | Task 3: dictionary and Q_( num) (49-50), eps-exponents (51-52), six CG data (53-54), octet mixing and cost (55-56); result 57 |
| 58-64, 65 | Task 2: driver (59-60, wipes `Global``), fit (61-62), u-form (63-64); result 65 |
| 66-68, 69 | Task 4: [3,1,0] driver (67-68); result 69 |
| 70 | summary |

## Appendix B Numeric points and grids

| **item** | **values (source)** |
|---|---|
| base point | h=1/2, theta=(2/7,-3/5), t=62/35, z=(3/2,5/4,2/3); controls z'=(3/2,5/4,3/7), z''=(2,5/4,2/3) ([A:1], [B:1]) |
| bra twists (B) | w=(5/2,7/3,3/5), w'=(7/4,-2/5,9/7) ([B:16]) |
| (h,theta) points (A driver) | (1/3,(1/5,-4/9)) [t=29/15]; (1/2,(3/7,-1/4)) [19/14]; h=1/2, theta_2=1/5,-26/105,-27/28,-37/42,-1/70,-3/7,-5/3,1/3 [t=6/35,16/15,5/2,7/3,3/5,10/7,82/21,-2/21] and, for [3,1,0] only, theta_2=2/5,-4/5,-2/3,-6/5 [-8/35,76/35,40/21,104/35] ([A:34]) |
| (h,theta) points (B driver) | base, z', z'', (1/3,(1/5,-4/9)), (1/2,(3/7,-1/4)), theta_2=-100 [t=1404/7] ([B:57]) |
| (h,theta) points (C driver) | 15 points, 13 distinct t: 62/35 (base; (2h,2theta); z=(5/3,-3/7,4/5)), 29/15, 6/35, 16/15, 5/2, 7/3, 3/5, 10/7, 82/21, -2/21, -16/21, 57/14, 26/21; held out: 57/14, 26/21 ([C:60], [C:65]) |
| Gauss points | undo/reach (A): g_2=(1,2,3),(-3,6,-1),(2/3,-3/7,5/4); closed-form verification on b_3=0: (g_A,g_C)=(5/2,-3/7),(-4/3,2/5),(7/3,11/4); site-1 adjoint quadric point (4,23,5); [3,1,0] site-1 locus points (-3,3,-2),(-3,6,-1),(-3,12,1),(-3,15,2) ([A:15], [A:24], [A:29]) |
| reach data at t=62/35 (A) | L_0=([[1,0,0],[-23/12,1,0],[5/6,-2/3,1]]), D_0=diag(1,-1/4,35/72), U_0=([[1,1,1],[0,1,10/3],[0,0,1]]) ([A:24]) |
| w points (A, B) | Omega_c[w] identity: (2/3,-5/7),(-3/2,4/5),(5/3,7/2) (A), four points in B; measure verification at (2/3,-5/7),(-3/2,4/5) ([A:15], [B:19]); transfer-matrix tests at u=3/7,19/23; R solve at u=3/11,2/9, check at 19/23 ([B:36]) |
| torus data (A) | d=(1,4,43) (d_1d_3/d_2^2=43/16), (2,3,9/2), (3/2,-5/4,7/3), converse (1,1,cc); product-point samples w=(7/3,-5/7),(-3/2,4/5),(5/3,7/2) ([A:31]) |
| (a,a') points (B 4b, C) | B: (2/3,-5/7),(3/2,1/4),(-4/5,7/3) (the first two degenerate vacua); C pivot scan (7/5,-4/3); C identity: (2/3,-5/7),(5/7,-4/9),(-3/5,7/4) ([C:24]) and (5/7,-4/9),(-3/5,7/4),(9/5,2/7) ([C:32]) |
| Stage D grid | phi=(1/3,-2/7,5/11), phi'=(-2/5,3/8,7/11); om_1 in {3/4,-5/3}, om_2 in {-1/5,9/7,5/3}, om_1' in {-12/5,7/2}, om_2' in {1/8,-3/4} (24 points); w_3=1+5eps/13, w_3'=1-4eps/13; eps=1,2,4*10^(-6); sub-families inf2, inf1, common1 (12 points each), fixed-w' control (6 points) ([D:1]) |
| Stage D det_6/four-overlap point | (om,om')=((3/4,-1/5),(-12/5,1/8)), eps=10^(-6) ([D:39]) |
| Stage C regime | same phi; om_1 in {3/4,-5/3}, om_1' in {-12/5,7/2}; exponent table at (3/4,-12/5) ([C:33], [C:57]) |
| precision | 60 digits (finite twist: A, B, C part A); 200 digits (D, C part B); exact rationals elsewhere |

## Appendix C Source notes: inconsistencies and ambiguities

Reported, not resolved silently. Items 1-6 were named in the task; 7-13 were found while compiling.

1. **The spec's 229/342.** [B:36] and [B:58] (item 3) state that the dressed adjoint conic constant at t=62/35 is 167/264=cc_1(-t)=(t+3)/(2(t+2)) and that the spec's value 229/342 is not equal to this (re-derived here: 167/264~=0.6326, 229/342~=0.6696) and is not matched by any shift tested; the formula, not the number, is what was confirmed (also at t=29/15,19/14,1404/7).

2. **Stage D attributes the 8/14 counts to the wrong family.** [D:23] (last bullet) and [D:40] (item 6.2) call the sub-family w_2=w_2'=inf "cheap (Stage B block counts 8/14)" and conclude that the "no cheap and CG-complete family" reading is contradicted. Per [B:46], [B:49] and [C:32] the exact 8/14 is the global-rotation line w_1'=w_1 (a=a'); the CG-complete family w_2=w_2'=inf, w_1!= w_1' costs 12 (10 after the pivot scan) for the adjoint and 45 for [3,1,0]. The reach statement of Stage D is correct; the cost attribution is not.

3. **Wedge count 3+3+2+2 versus sum P=12.** [B:49] lists "3 terms" in blocks (0,0),(0,1) and 2 in (1,0),(1,1) (ten terms) next to P=4,4,2,2 (twelve products). [C:32] resolves it: ten is the number of rank-one terms =sum r; twelve is the default-pivot count with one rank-4 (two-wedge) ladder-1 factor in each of the two rank-3 blocks; the full pivot scan attains ten.

4. **The det_6 Laplace bound.** The plan (S1, S7) expected sum r<=20 (20 products of per-ladder 3x3 determinants). [B:41] derives the per-block bound r_(cc')<=20 and the sum bound 80, [B:58] (item 5) states each block has r<=min(20,n_(min))=n_(min) (non-binding); both say sum r<=20 is not valid (sum r=26 for the det_6 vacuum at [3,1,0]). The two statements in B ("bounded by 80" in [B:41] versus "non-binding per block" in [B:58]) are different bounds for the same sum; neither is <=20. Also the units differ: the plan's 20 are products of 3x3 bracket minors, the notebooks' 18 / 40 are products of 2x2 minors, and [C:56] counts det_6 as 36 bracket entries plus a 6x6 determinant.

5. **The [3,1,0] greedy result.** [C:68], [C:69] and [C:70] (item 8) give: pivot scan best 45 (12,12,9,12); greedy elimination best 48 (12,12,9,15); default pivots 45. The running summary (`HANDOFF.md`) states that the scan *and* the 400-restart greedy give 45 (12,12,9,12). The notebook is followed here.

6. **Degenerate regression points.** [C:24] and [C:70] (item 10(iii)) note that (a,a')=(2/3,-5/7) has w_1=z_1 and (3/2,1/4) has w_1=z_3: degenerate vacua with 22 of 64 states orthogonal to the covector. Stage B's Task 4b used these two points among its three ([B:49]) without remarking on it, and states the overlap identity at (2/3,-5/7); that identity is partly trivial there. The block ranks and counts are K-based and are not affected (identical numbers at the third, generic pair). Moreover [C:24] and [C:32] check the identity at different point sets (the first set contains the degenerate (2/3,-5/7), the second does not).

7. **Two cell-numbering conventions.** A and C cite absolute positions in the file. The result cells of B and D cite the ordinal of the *code* cell: [B:14] "cells 1-6 above" (absolute 2,4,...,12), [D:15] "cells 1-6", [D:19] "cell 7" (absolute 17), [D:23] "cell 8" (21), [D:27] "cell 9" (25), [D:31] "cell 10" (29), [D:35] "cell 11" (33), [D:39] "cell 12" (37). The stage specs and this write-up use absolute positions (the Stage D pointers "reach cell 21/23, Gauss parameters 17, K/S 29, CG 33, overlaps 37-39" are absolute: code cell then result cell).

8. **Adjoint conic fit points.** [A:24] and [A:36] state that cc_2(t)=(t+2)/(2(t+1)) was fitted over 15 (h,theta) points; [A:34] lists 11 points for [2,1,0] (base, two (h,theta) points, eight theta_2 values) and 15 for [3,1,0] (four more added after a first run in which one [3,1,0] coefficient needed numerator degree 4 over denominator degree 5), while its table says "all 15 sampled values" for both reps. Whether the adjoint was re-run with the 15 points is not stated.

9. **Two meanings of "generic" in B.** "Generic pair" ([B:41]) is a pair of independent rotations (23 / 54); "generic (4 parameters)" in [B:46] is the generic member of Omega_c[w](x)Omega_c[w'] (22 / 54); the running summary's "Omega_c[w](x)Omega_c[w'] 13/22, 27/54" is the latter.

10. **The 8/14 line inside B.** [B:46] introduces w_1'=w_1, w_2,w_2'->inf as a "1-parameter" minimal family (8 / 14); [B:49] shows it is the global rotation a=a' only. Same fact, restated in the follow-up; the earlier text was not amended.

11. **Pivot-dependent counts that differ between sources.** The Stage A block-form product count is 15 (explicit decomposition 16, random covector 16) for the adjoint and 35 (38; 35) for [3,1,0] ([A:20]); [B:41] quotes 15 / 35. The per-label counts (17 / 40) and the block counts (15 / 35) are different upper bounds for the same covector.

12. **"Any w_3".** [A:15] and [B:19] state that Omega_c[w](x)Omega_c[w] is an eigen-covector of t_a with G(w_1,w_2,w_3) for any w_3; the tests are at w_3=3/2,-4,7/3 and one second w (w_3=5/3) ([A:15]) and at two w in B: exact at those points.

13. **Minor.** [C:57] says the top-irrep slice 2.6094674556 has "no earlier reference number"; [D:35] lists 2.60946745562 for [4,2,0] and [6,2,0]. [D:31] prints the K exponent of [3,3,2] as -6.000,-5.997 in the main family and -6.000 elsewhere (within the stated tolerance). [A:29] and [A:24] refer to the same t-fit as "cell 33" (code cell 33, analysis of the driver data).

**Not reconciled / not found.** (a) The support rule of e_1(x) e_nu (triangular set of pairs, no additive weight rule) was not characterised ([A:18]). (b) No cell derives the relation between cc_1(1-t) and cc_1(-t) (open item 4 of S6). (c) The termwise comparison of the 20 det_6 Laplace terms with the blocks was not done ([B:41]). (d) Whether the type-pair loci of [B:54] (e.g. (alpha_1=0,beta_2=0), 40 for [3,1,0]) and the line beta_2=0 are CG-complete was not tested. (e) The explicit 12-wedge decomposition of the adjoint (inf,inf) family at (2/3,-5/7) printed in the output cell of [B:47] (cell 48) was not read cell by cell for this write-up (only its block structure, quoted in [B:49]); the closed forms used here are those of [C:29], [C:31], [C:62].
