# Brief: literature check on fast Clebsch–Gordan coefficients for SU(4) / SU(N ≥ 4), the outer-multiplicity labels, and the claims of the gl(4) notes (2026-10-22)

For a general-purpose subagent (SONNET) with web search/fetch. **No repo code, no git, no notebooks**; a throwaway script in the session scratchpad for reading arXiv
metadata is fine (as yesterday). Output: `Paul/Mathematica/Gaudin/Literature-check-SU4-CG.md` in the style of yesterday's `Literature-check-missing-label.md` (read it first, whole:
the access legend `[FT]` full text / `[AB]` abstract / `[CIT]` via a citing paper / `[MEM]` memory-unverified, the per-claim verdict table in §0, one section per claim, "premises
that did not check out", "recommended framing and must-cite", "what I could not access", an appendix of verified bibliographic data with DOIs/arXiv ids). Everything not
tagged `[FT]` is "I know what it is about, not what it proves". Do not repeat yesterday's SU(3)-specific findings; reuse its general-N references (Kostant/PRV, Munich,
Gaudin-algebra literature) and extend them. Time-box ≈ 3 hours; deliver what is finished and list what is not.

Read first, in this order: yesterday's check (whole); `Paul/Mathematica/Gaudin/gl4-slice-engine-notes.md` §1 (punchline), §4.3, §5.4, §6, §7, §9 (what we claim and with what
status); `Paul/Python/gaudin_gl4/benchmark_gl4.md` (the numbers we would quote); `Paul/Cpp/AlexCGC/README.md` (the Munich code as we use it).

## 0. The question

For SU(3) the state of the art turned out to be the Draayer–Akiyama lineage (SU3lib, ndsu3lib), SU(3)-specific, as fast as our engine for small reps, and the cubic
missing-label operator was known (Pluhař–Smirnov–Tolstoy, Crampé–Poulain d'Andecy–Vinet). **Is there an SU(4) (or general SU(N), N ≥ 4) analogue of any of this?** Concretely:
(i) which codes and algorithms compute SU(4) / SU(N) Clebsch–Gordan coefficients numerically, how fast, up to what size, with which multiplicity convention; (ii) what is known
about the outer-multiplicity (missing-label) problem for SU(4)⊗SU(4) ⊃ SU(4) and for U(n) in general, including Gaudin/Bethe-algebra labels; (iii) which of the structural
statements of our gl(4) notes are known, partially known, or new. The deliverable is a per-claim verdict with citations, a table of codes, and a must-cite list for the gl(4)
notes/paper. **Search broadly before judging** (nuclear-physics SU(4) "Wigner supermultiplet" literature, particle-physics SU(4) tables of the 1970s, the tensor-network/quantum-
simulation SU(N) literature, the U(n) representation-theory literature, the Gaudin/Bethe-algebra literature), and record for each item how you accessed it.

## 1. Claims to check (give each a verdict: KNOWN / PARTIALLY KNOWN / NOT FOUND / CONTRADICTED, with the evidence)

**K1 — Numerical SU(4) and SU(N) CGC codes (the benchmark risk).** Candidates to locate and characterise (algorithm, basis, multiplicity resolution, general N or not, largest
reported case, reported timings, availability): Alex–Kalus–Huckleberry–von Delft 2011 (the Munich code; does the paper or anything citing it report SU(4) timings or limits?);
Weichselbaum's QSpace tensor-network library (Ann. Phys. 2012; Phys. Rev. Research 2020 "X-symbols"; later papers), which computes SU(N) CGC tensors for arbitrary N — what algorithm,
what sizes, what multiplicity ("outer multiplicity") treatment, any timings; Gliske–Klink–Ton-That (Acta Appl. Math. 2007?, "algorithms for computing U(N) Clebsch–Gordan
coefficients") and Klink–Ton-That (symbolic / generalized Wigner operators); Rowe–Repka vector-coherent-state U(n) Wigner coefficients (J. Math. Phys. 1997 and later), Rowe's "U(n)
Wigner coefficients via … " ; Vlasii–(?) or other recent arXiv "SU(N) Clebsch–Gordan via Young tableaux / Gelfand–Tsetlin" (2015–2026; e.g. search "Clebsch-Gordan coefficients
SU(N) algorithm", "arbitrary SU(N) Clebsch-Gordan", "Gelfand-Tsetlin Clebsch-Gordan numerical"); Mathematica/other packages (LieART, GroupMath, Susyno, SymPy, LiE, Sage — which of
them actually compute CGCs for SU(4) and how); nuclear SU(4) (Wigner supermultiplet SU(4) ⊃ SU(2)×SU(2)): Hecht–Pang 1969 J. Math. Phys. 10, 1571, Partensky–Maguin 1978, Draayer
and collaborators' SU(4) ⊃ SU(2)×SU(2) codes (1980s–2000s), Kota's SU(4) Wigner–Racah work, Alisauskas (SU(4) isofactors, U(n) isoscalar factors); particle-physics SU(4) tables
(Haacke–Moffat–Savaria 1976 J. Math. Phys. 17, 2041 and related 1970s SU(4) flavour tables; "SU(4) Clebsch–Gordan coefficients" 1973–1978); lattice-gauge/quantum-simulation uses
of SU(N) CGCs (Hamiltonian lattice gauge theory, "SU(3)/SU(4) Clebsch–Gordan for qudit simulations", 2020–2026); any code that computes the HIGHEST-WEIGHT coefficients as a
separate step (Munich does, Eq. 35) and any that uses a slice / string-condition / extremal-projector reduction for general N. **Verdict wanted:** is anything faster than the
Munich code for SU(4) highest-weight CGCs or for full tables, and does anything pre-empt the slice engine for N ≥ 4? Produce a table (code | algorithm | basis | multiplicity
convention | N | largest case | timing | source | access).

**K2 — The outer-multiplicity (missing-label) problem for SU(4)⊗SU(4) ⊃ SU(4), and U(n) in general.** How many labels (½(N−1)(N−2) = 3 for N = 4) and which operators are
proposed: Biedenharn–Louck canonical unit tensor operators / "operator patterns" for U(n) (the canonical resolution; Biedenharn, Louck, Chacón, Ciftan 1972–1985; Louck–Biedenharn
on U(n) canonical multiplicity labels); Alisauskas's labelling operators for U(n) tensor products; Sharp and others on "fourth-order external labeling operators"; Moshinsky's
complementary-group approach (Pan–Draayer 1998 has SU(4) closed multiplicity formulas); extremal-projector constructions for U(n) (Tolstoy's review arXiv:1010.4054 and the
Asherova–Smirnov–Tolstoy papers; Pluhař–Smirnov–Tolstoy generalisations beyond SU(3)?); and the Gaudin/Bethe-algebra labels for gl_N: Vinberg's problem, Rybnikov (2006, 2020),
Feigin–Frenkel–Rybnikov (2010), Mukhin–Tarasov–Varchenko (2009 and the "Bethe algebra of gl_N Gaudin model" papers), Halacheva–Kamnitzer–Rybnikov–Weekes (2020), Chervov–Falqui–
Rybnikov 2010 "Limits of Gaudin algebras, quantization of bending flows, Jucys–Murphy elements and Gelfand–Tsetlin bases" (Lett. Math. Phys.; the two-site Gaudin algebra and GT
bases — read it, it may contain our V_red commutation-on-weight-spaces / GT-lattice structure), Chervov–Talalaev / Talalaev "quantum spectral curve" (the cdet we use), and any
explicit two-site sl_N or sl_4 treatment (Campoamor-Stursberg–Musso did sl(3); is there sl_N?), and CPV-type follow-ups on the diagonal centralizer of U(sl_N) (Crampé, Poulain
d'Andecy, Vinet and coauthors 2021–2026: search "missing label sl_N", "centralizer diagonal embedding U(sl_n)", "Lehrer–Racah"). **Verdict wanted:** is "one cubic + two quartic
Gaudin Laurent coefficients as the three SU(4) missing labels, with simple joint spectrum" stated anywhere; is the exchange-parity degeneracy of the cubic label known.

**K3 — The slice / PRV reduction and the level recursion as an algorithm for U(n).** Beyond yesterday's findings (Kostant 1959, PRV 1967, extremal projectors, Braunschweig–
Williams for SU(3), Rowe–Bahri VCS): any U(n)/SU(N) numerical implementation of "highest-weight vector from its component at the second factor's highest weight" (Zhelobenko's
Gauss-decomposition approach to U(n) Wigner coefficients; Klimyk; Gould's U(n) Wigner coefficient formulas via characteristic identities — e.g. M. D. Gould 1980s J. Math. Phys.;
Alisauskas's "U(n) isoscalar factors via recursive methods"); any mention of the "dangling"/string conditions in a numerical recursion. **Verdict wanted:** prior art for G1.

**K4 — Structural statements of the notes (§5.4, §6, §7).** For each, search and give a verdict: (a) the Gaudin operators on V_λ[μ−λ′] for every (also non-dominant) λ′ commute
and obey the BGG intertwining E_i^{k_i}V(λ′) = V(s_i·λ′)E_i^{k_i} — Gaudin/Bethe algebras on tensor products with Verma modules, singular vectors and the dot action (MTV, Varchenko
school; "Bethe algebra Verma module", "Gaudin model Verma modules BGG resolution"); (b) the Weyl/BGG identity for the physical spectrum Σ_w sgn(w)·mult_τ(V(w·λ′) on V_λ[μ−w·λ′]) =
[τ physical] (search "Euler characteristic" + "Bethe eigenvalues", "Weyl group" + "Gaudin spectrum" + "Verma"; Kostant's multiplicity formula is the dimension shadow); (c) the forced
dictionary: Laurent coefficients of the quantum characteristic polynomial on a joint eigenvector are fixed by the central characters of λ, λ′, μ plus the free accessory
parameters (this is standard oper theory — find the cleanest statement: FFR 2010 "opers with regular singularities", Frenkel's reviews, Chervov–Talalaev); polynomial solutions
of degrees μ_j + N − j and the Wronskian ∝ Π(u−θ_s)^{|λ_s|} (MTV "spaces of polynomials with prescribed exponents", Schubert calculus; the "reproduction"/direct–dual relation =
formal adjoint + gauge Π^{λ₁−λ_N}); the exact Q-system/Bethe-ansatz solve being complete and needing all N degree conditions (MTV completeness, "Bethe ansatz is complete for
gl_N Gaudin"); (d) the GT-lattice stencil and product closed forms for matrix elements of mixed invariants (Gould's characteristic identities give matrix elements of invariant
polynomials in the GT basis; Perelomov–Popov Casimirs; Alisauskas); (e) CG coefficients / eigenvector components polynomial in the labelling eigenvalues (Hahn-type structure
for SU(2), CPV for SU(3); anything for N ≥ 4; "Krylov"/cyclic-vector arguments); (f) the correspondence of our oper with the h → 0 limit of an XXX Baxter equation (standard;
cite the cleanest reference for gl_N: e.g. "Gaudin model as the rational limit of XXX", Mukhin–Tarasov–Varchenko "XXX Bethe algebra / Gaudin limit").

**K5 — Speed claims in print.** Any published timings or complexity statements for SU(4)/SU(N) CGC computation (Munich paper figures/tables; Weichselbaum; Gliske et al.;
nuclear SU(4) codes); the cost model of Munich (null space O(n_W³)) as stated by the authors; our claims (×860 vs Munich as distributed, two machines; S = 20 in 4.6 s) — what
can honestly be said against the literature.

**K6 — Descendants.** Standard references for building the full CG table from the highest-weight vector by lowering in the GT basis (Munich's Eq. 35 ff.; Alex et al.;
Louck–Biedenharn pattern calculus); any statement that the highest-weight step dominates / does not dominate the cost.

**K7 — Anything else** that bears on the gl(4) notes: SU(4) basis conventions in the physics tables (SU(2)×SU(2) or SU(3)×U(1) chains vs the canonical GT chain; conversion
literature), symmetry properties of SU(4) CGCs, the "canonical" multiplicity basis debates (Biedenharn–Louck vs Gaudin/Bethe vs Gram–Schmidt), and recent (2023–2026) arXiv
papers on SU(N) CGCs or on Gaudin models and multiplicity labels we have not seen.

## 2. Deliverable

`Literature-check-SU4-CG.md`: §0 bottom line + the verdict table (K1–K7, one line each); §1 the codes table (K1); §§2–8 one section per claim with the evidence and access
tags; §9 premises that did not check out; §10 recommended framing for the gl(4) notes/paper and the must-cite list (ordered by how much each reference constrains our claims),
with the specific sentences in `gl4-slice-engine-notes.md` §9.2 that should change; §11 what could not be accessed; Appendix: verified bibliographic data (DOI / arXiv id /
journal, as yesterday). Be concrete and honest: "not found" after a stated search is a result; "UNVERIFIED" where it is memory. Reply with a ≤ 40-line summary: the verdict
table, the three most important references, and whether anything pre-empts (i) the engine, (ii) the triple, (iii) the Weyl identity.
