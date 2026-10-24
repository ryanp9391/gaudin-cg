# RESULTS G8a: block-Jacobi structure of the gl(4) triple on lattice gradings, and the forward recursion `slice_from_triple` (2026-10-24)

Brief: `SPEC-G8a-lattice-recursion.md`; plan: `Paul/Mathematica/Gaudin/gl4-jacobi-plan.md` §1 H-A. Code `lattice_recursion_gl4.py` (CLI `python -X utf8 lattice_recursion_gl4.py gradings|recursion|charpoly|n0|separability|all`), tests `test_lattice_recursion_gl4.py` (8 tests, 3 s, all pass), raw output `lattice_recursion_gl4_output.txt`, full tables `lattice_recursion_gl4_tables.txt` (per grading, per weight space, per class A, per sector). Whole `all` run: 21 s. Everything is exact (Fractions, exact rref / exact left inverses, sympy QQ polynomials) except where "50 digits" is said.

## 0. Headline

**H-A is true in the following precise sense, and false in the strong sense hoped for.** Data set: the 222 weight spaces with w >= 2 of the 12 products P_small + P_big (sizes w = 2..15; 171 of them have a void E12-cut, W1 = W; 1 has m = 0).

1. **Forward recursion exists for every weight space.** Taking the best grading *per weight space*, the triple is block-banded and fully forward determined (every stacked inter-level block has full column rank, exact) on **222 of 222** weight spaces of W1 (and of W), with bottom-level size n0 = 1 on **188 / 222 (84.7 %)** of W1 (180 / 222 of W), n0 = 2 on 33 (40), n0 = 3 on 1 (2), never more. **No single grading works everywhere:** no grading of the 124 is fully determined on all 222 spaces. Best single gradings (W1): **+m12** determined on 218/222 (98.2 %; the 4 failures have rank deficit 1 at the first step: two m = 1 sectors of [4,2,0,0]^2, level sizes [1,2,1] and [1,3], and two w = 12 sectors of [6,4,2,0]^2, level sizes [2,6,3] with rank 5 of 6) but n0 = 1 only on 41; **-m23** determined on 214 with n0 = 1 on 82 (37 %); **+m33 (=(0,-1,-1))** 216 with n0 = 1 on 37; **-m12** 215 / 73. For the large spaces (w >= 9, 23 spaces, mostly [6,4,2,0]^2) n0 = 1 is possible in only **6 of 23** (all w = 10): w = 9, 12, 15 need n0 = 2 (w = 12: 11 of 12; one needs 3). So n0 = 1 is **not** a property of the whole product; the counting obstruction n_{d+1} <= 3 n_d (necessary for full column rank of the stacked block) and genuine rank drops produce n0 >= 2.
2. **The recursion is implemented and validated:** `slice_from_triple(lam, lamp, mu, trip, grading, space)` reproduces **358 / 358** copies (137 sectors of P_small ∪ [6,4,2,0]^2 whose W1 has an n0 = 1 grading; the first-choice grading works for all) to **<= 2.9e-43 at 50 digits**; the 161 remaining copies of [6,4,2,0]^2 / P_small lie in n0 > 1 sectors and are covered by item 4: **203 / 203** copies of all 34 sectors with minimal n0 = 2 or 3 (incl. w = 12, 15) verified through the polynomial-matrix formulation. The unused (characteristic) equations vanish at the true triples (<= 4e-45 absolute) and not at a perturbed triple (>= 0.023). Recursion 0.3 ms per copy (8x faster than the 50-digit joint eigen-decomposition of the sector; exact plan 0.1 ms per sector).
3. **Characteristic equations (n0 = 1).** ψ_x = P_x(v3,v4a,v4b) with **total degree of P_x = level of x exactly** (14/14 sectors; max degree 2 for m <= 10 in [6,4,2,0]^2), and **equal to G4's minimal Krylov degree on the PRV coordinates in 14/14**. The unused equations C_j(v) (9 to 21 per sector, degrees 2–3) vanish at every joint triple of W1 (both physical and spurious; the bottom component was nonzero in all 14 sectors' spectra) and generate a zero-dimensional ideal J with **dim Q[v]/J = dim W1 = number of distinct W1 triples in 14 / 14**, i.e. **J is the radical (vanishing) ideal of the joint spectrum of W1**: the gl(4) continuant is the Gröbner basis of the C_j. The E23/E34 conditions that cut W1 down to M_mu are *not* among them (spectral): e.g. [4,2,2,0]^2 (6,4,4,2): dim W1 = 6, dim M~ = m = 3.
4. **n0 > 1.** ψ = T(v)ψ0 with polynomial T; characteristic equations C(v)ψ0 = 0, C a polynomial matrix (entry degrees 1–2 for n0 = 2). Verified at 50 digits at **203 / 203** copies of **34 / 34** sectors (n0 = 2: 33, n0 = 3: 1): rank C(v) = n0 − 1 (σ_{n0}/σ_1 < 1e-40, σ_{n0−1}/σ_1 > 1e-15), ker C(v) = the copy's bottom-level component (1 − cos <= 1.6e-61), T(v)·ker = the copy (error <= 4.2e-44); at a perturbed triple C has full rank. For n0 = 2 the 2×2 minors (41–66 distinct, degrees 2–4) generate, in the 5 sectors tested symbolically, a zero-dimensional ideal of dimension exactly dim W1 (4,4,4,5,5), i.e. the vanishing ideal again, not just a multiple.
5. **Vandermonde fact corrected** (§2). **Separability (A5): false** — Ψ[κ][m12] has rank 1 in only 4 of the 146 (of 304) physical copies with a non-degenerate support; the other two cuts never.

## 1. Conventions

W = V_λ[μ − λ'] with the plain GT basis of `glN_reps.GLNRep` (global indices ascending); lattice coordinates x = (m12, m13, m23) = `structure_gl4.coords`; the three operators = `slice_ops_gl4.weight_space_mats(...)`, M[row][col] = coefficient of e_row in V e_col. Spaces: **W** (whole weight space) and **W1** = {m12 <= μ1 − λ'_2} (`Sector.in1`). Verified exactly on all 222 spaces: W1 = {c : E12^{k1} e_c = 0} (`power_images`) and W1 is invariant under the triple (0 violations). A grading g = a m12 + b m13 + c m23 is *directed*: levels ascending = recursion direction; g and −g are both analysed (124 directed gradings = 62 pairs; **only 98 directed (49 pairs) are primitive**, the rest are multiples with identical level structure and results). Levels are classified by the level index; **band** = max level-index distance of a nonzero entry (the value-band equals the structural bound |a| + max(|b|,|c|,|b−c|) for all 124 gradings: the bound is attained). Forward determinacy (greedy, = "every F_d has full column rank" for band 1): process levels ascending; at level d the 3 n_d rows determine all not yet known unknowns of levels <= d + B at once; full determinacy iff the stacked block on them has full column rank at every step (exact rank over Q); failure = (level, rank deficit). The bottom level L_0 supplies the initial data ψ0 (n0 numbers).

The recursion (docstring of `slice_from_triple`): F_d ψ_{d+1} = [ v_k ψ_d − Σ_{c known} M_k[L_d, c] ψ_c ]_{k=3,4a,4b}, F_d = [M_k[L_d, L_{d+1}]]_k (3 n_d × n_{d+1}); solved with the exact left inverse L_d (computed once over Q); the unused equations are T_d·rhs with T_d F_d = 0 (3 n_d − n_{d+1} per level; all 3 n_top at the top; in total 3n − (n − n0) = 2n + n0 unused scalar equations).

## 2. A1: gradings (data: `lattice_recursion_gl4_output.txt`, `..._tables.txt`)

Weight spaces by size: w = 2: 66, 3: 17, 4: 54, 5: 26, 6: 13, 7: 23, 9: 4, 10: 6, 12: 12, 15: 1. Products (count): [2,1,0,0]^2 3, [2,1,1,0]^2 1, [3,2,1,0]^2 8, [4,2,0,0]^2 19, [4,2,1,0]^2 18, [3,2,1,0]x[4,2,1,0] 10, [4,2,1,0]x[3,2,1,0] 14, [4,2,2,0]^2 8, [4,3,1,0]^2 22, [6,4,2,0]^2 62, [5,3,1,0]^2 39, [5,3,2,1]^2 18.

Named gradings, space W1 (det = fully forward-determined, of 222; n0=1 = determined with n0 = 1; same for w >= 5, 85 spaces):

| grading | det | n0=1 | det (w>=5) | n0=1 (w>=5) | worst deficit |
|---|---|---|---|---|---|
| +m12 (bottom far from the cut) | 218 | 41 | 83/85 | 3 | 1 |
| −m12 (bottom at the cut) | 215 | 73 | 80/85 | 11 | 1 |
| +m13 | 211 | 90 | 74/85 | 24 | 1 |
| −m13 | 212 | 32 | 81/85 | 8 | 2 |
| +m23 | 202 | 87 | 72/85 | 20 | 2 |
| −m23 | 214 | 82 | 78/85 | 20 | 2 |
| +m33 | 216 | 37 | 82/85 | 7 | 2 |
| −m33 | 217 | 63 | 80/85 | 11 | 1 |
| +(m12+m13) | 172 | 117 | 38/85 | 14 | 4 |
| −(m12+m13) | 162 | 98 | 36/85 | 10 | 4 |
| +(m12−m13) | 163 | 100 | 34/85 | 8 | 8 |
| −(m12−m13) | 199 | 81 | 65/85 | 10 | 4 |
| +(m12+m13+m23) | 181 | 108 | 46/85 | 10 | 4 |
| −(m12+m13+m23) | 170 | 109 | 38/85 | 9 | 4 |

(Space W: same pattern, +m12 218 / 41; −m12 only 164 (worst deficit 6); −m23 210 / 74; +m13 207 / 77.) Band-1 gradings (±m12, ±m13, ±m23, ±m33) are by far the best for determinacy; the band-2 diagonal gradings give more n0 = 1 (thinner bottom levels) but are determined on only 34–65 of the 85 spaces with w >= 5. Full 124-row tables per space in the tables file.

Failures of +m12 on W1 (all at level 0, deficit 1): [4,2,0,0]^2 (6,4,1,1) levels [1,2,1] and (5,3,2,2) [1,3]; [6,4,2,0]^2 (9,7,6,2) [2,6,3] and (9,7,4,4) [2,6,3] (rank 5 of 6 unknowns). Failures of −m12: (7,3,1,1), (6,4,1,1), [4,3,1,0]^2 (7,4,3,2) [1,4] (a counting failure 3 < 4), and others in the tables.

**Best grading(s) per weight space** (smallest n0, then band): in the tables file, one line per weight space. Smallest n0 over all determined gradings: W1: n0 = 1: 188, 2: 33, 3: 1; W: 180 / 40 / 2. By size (W1): w <= 3 always 1; w = 4: 51 × 1, 3 × 2; w = 5: 21 × 1, 5 × 2; w = 6: 12 × 1, 1 × 2; w = 7: 15 × 1, 8 × 2; w = 9: 4 × 2; w = 10: 6 × 1; w = 12: 11 × 2, 1 × 3; w = 15: 1 × 2.

**Rule "A -> grading"** (W1, A = active strings; counts = spaces of the class where the grading is determined; classes: A=- n=63, 1: 22, 2: 45, 3: 42, 12: 8, 13: 18, 23: 21, 123: 3):

| class | n | primitive directed gradings determined on EVERY space of the class (count; first ones) |
|---|---|---|
| A = ∅ | 63 | 8: −m12, (−1,+1,0), +m33, −m13, ±m23, … (+m12 and −m12 both 63/63) |
| 1 | 22 | 10: (−2,−1,0), (−1,−1,0), −m12, +m33, −m13, −m23, … (+m12 22/22) |
| 2 | 45 | **3: +m12, −m33 = (0,1,1), (−1,1,1)** (−m12 fails once, 44/45) |
| 3 | 42 | **2: +m33 = (0,−1,−1), −m13** (+m12 40/42, −m12 37/42) |
| 12 | 8 | 5: −m12, −m23, −m33, (−1,1,1), (1,0,−1) (+m12 7/8) |
| 13 | 18 | 5: −m12, (−1,1,0), +m33, −m13, (1,0,1) (+m12 17/18) |
| 23 | 21 | 3: **+m12**, +m33, (−1,1,1) (−m12 20/21) |
| 123 | 3 | 8 |

So there is a partial rule: **+m12 is determined on every space of the classes {−,1,2,23,123}, −m12 on {−,1,12,13,123}; the class A = {3} (E34 string active alone) needs +m33 or −m13** (and no single grading is n0 = 1 on all spaces of any class; the best is 19 of 22 for A = {1} with +m33). The rule is empirical (no proof); it is *not* "the cut end is the bottom": both ±m12 work in most classes.

**Vandermonde side check (corrected).** For the three Δm12 = +1 root moves of the row 3 with distinct spectators k = 1, 2, 3 and *coherent orientation* — {(Δm13,Δm23,Δm33)} = {(0,−1,1), (−1,0,1), (−1,1,0)} resp. its negative — 38 sampled (source, sector) pairs (all 38 that exist in the products with w >= 3, seed 7): det[(V3,V4a,V4b) entries] = **+1 · Π_α c_α (y2−y1)(y2−y3)** for the first set (19 of 19) and **−1 ·** the same for the second set (19 of 19) (y_k = m_{k3} − k of the source, c_α the V3 entries). For mixed orientations (e.g. {(0,−1,1),(1,0,−1),(1,−1,0)}) the ratio is *not* ±1 (values 2, 7/4, −3, −2, −5/2, 3 …; 24 of 52 triples over [4,2,1,0]^2, [6,4,2,0]^2). So the plan's formula holds exactly for the two coherent sets, up to the sign −1 for the opposite one.

## 3. A2: recursion validation

`slice_from_triple` on W1 with the best n0 = 1 grading (mostly −m12 ×55, −m23 ×29, +m23 ×23, +m13 ×20, −m13 ×9): per product, sectors with m >= 1, w >= 2 (copies = all copies of the sector):

| product | sectors validated | no n0=1 grading (copies untested) | copies reproduced (first grading) |
|---|---|---|---|
| [2,1,0,0]^2 | 3 | 0 | 4/4 |
| [2,1,1,0]^2 | 1 | 0 | 2/2 |
| [3,2,1,0]^2 | 8 | 0 | 22/22 |
| [4,2,0,0]^2 | 15 (1 with m=0) | 3 (5) | 21/21 |
| [4,2,1,0]^2 | 18 | 0 | 47/47 |
| [3,2,1,0]x[4,2,1,0] | 10 | 0 | 28/28 |
| [4,2,1,0]x[3,2,1,0] | 13 | 1 (4) | 28/28 |
| [4,2,2,0]^2 | 8 | 0 | 15/15 |
| [4,3,1,0]^2 | 19 | 3 (13) | 47/47 |
| [6,4,2,0]^2 | 42 | 20 (139) | 144/144 |
| total | 137 | 27 (161) | 358/358 |

Max error 2.9e-43 (G3 normalisation: first nonzero GT component = 1), max unused residual at true triples 4.0e-45, smallest residual at the relatively perturbed triple (1e-3) 0.023. The 27 n0 > 1 sectors (161 copies) of this table are among the 34 verified in §5. Timing: 0.112 s for 358 copies (0.3 ms/copy) vs 0.85 s for the joint eigen-decompositions of the same 137 sectors; exact plans 0.013 s in total. A band-2 grading ((1,1,0)) also reproduces copies exactly (test).

## 4. A3: polynomials and characteristic equations

14 sectors (11 of P_small with m = 3..5, 3 of [6,4,2,0]^2 with m = 10, 7, 6; all with an n0 = 1 band-1 grading on W1; `charpoly` stage). Per sector (tables file): ψ_x = P_x(v3,v4a,v4b) exact; **deg P_x = level index of x exactly** (levels 0..1 or 0..2, so total degree <= 2 even for m = 10), equal to G4's minimal graded-Krylov degree max_y deg(y) on the PRV coordinates in **14 / 14** (and to d_gen(m) where G4 has it); the number of unused scalar equations is 2 w1 + 1 (9 for w1 = 4, 13 for w1 = 6, 21 for w1 = 10), 7–21 of them non-zero (12 or 13 of the 13), of degree mostly 2–3 (degree 1 also appears in [4,2,0,0]^2 (6,4,2,0) and [4,2,2,0]^2 (6,4,4,2)); the nonzero C_j, e.g. for [4,2,1,0]^2 μ = (6,4,3,1) (g = −m12): `3 v3^2 − 2 v3 v4a − 378 v3 + 155 v4a − 15 v4b + 12240`, `18 v3^2 − 11 v3 v4a − 2268 v3 − 13 v4a^2 − 6 v4a v4b + 2334 v4a + 50112`, …. Evaluated at all joint triples of W1 (physical and spurious): all C_j vanish (<= 1e-42 absolute, 1e-45 for the small ones, coefficient sizes up to 1e7) in 14/14 sectors, and at the physical triples in 14/14; the bottom component of every W1 eigenvector was nonzero (0 exceptions). Gröbner (grevlex, sympy): dim Q[v]/J = dim W1 = #distinct W1 triples in 14/14. In the language of RESULTS-jacobi-cg §4: the **gl(4) continuant is the set of the C_j; they generate the vanishing ideal of the joint spectrum of W1**. M~ and M_mu: in all 14 sectors dim M~ = m (so M~ = M_mu here) and the physical triples are exactly the M~ triples; W1 has additional spurious triples in 4 of the 14 sectors (1, 1, 3, 1 extra; e.g. [4,2,2,0]^2 (6,4,4,2): dim W1 = 6, dim M~ = m = 3). The C_j vanish on those spurious triples as well, so none of them is "needed" for M~: the E23 string condition that removes them is spectral (not a polynomial consequence of the recursion on W1). For the 10 sectors with W1 = M_mu (no spurious triples) the C_j alone give the physical spectrum.

## 5. A4: n0 > 1

Minimal n0 over gradings by (w, m, A): table in the tables file (A4 section). Summary (W1): by w: n0 = 1 always for w <= 3, nearly always for w = 4 and 6, 10; w = 5: 21/26; w = 7: 15/23; **w >= 9: 6/23**; by A: A=− 58/63 (n0=1), 1: 19/22, 2: 43/45, 3: 35/42, 12: 4/8, 13: 13/18, 23: 15/21, 123: 1/3. **Is there a grading with n0 = 1 on every sector of a product?** For [2,1,0,0]^2, [2,1,1,0]^2, [3,2,1,0]^2, [3,2,1,0]x[4,2,1,0], [4,2,2,0]^2: yes (single grading, e.g. (−2,−2,−1) ~ (−2,−2,−1)·: n0 = 1 in 3/3, 8/8, 10/10, 8/8); for [4,2,0,0]^2 16/19 (some grading per sector gives 16), [4,2,1,0]^2 15–16/18 by one grading (18/18 with per-sector choice), [4,3,1,0]^2 19/22, [4,2,1,0]x[3,2,1,0] 12/14, [5,3,1,0]^2 26/39 (per-sector 32), [5,3,2,1]^2 15–16/18 (18), [6,4,2,0]^2 29/62 by one grading (per-sector 42/62). No for the larger products.

Verification (50 digits) on all **34** sectors (m >= 1) with minimal n0 >= 2 on W1 (33 with 2, 1 with 3; includes w = 12 and 15 of [6,4,2,0]^2): 203 copies; rank C(v) = n0 − 1 at all 203; bottom component of the copy nonzero in all; the kernel of C(v) equals it (max 1 − cos = 1.6e-61); T(v)·ker reproduces the copy (<= 4.2e-44); at a perturbed triple σ_{n0}/σ_1 is O(1e-2..1) (full rank). Symbolic C(v) (n0 = 2, 5 sectors): 10–12 nonzero rows, entries of degree 1–2, 41–66 distinct nonzero 2×2 minors (degrees 2, 3, 4), ideal of minors zero-dimensional with dimension 4,4,4,5,5 = dim W1.

## 6. A5: separability

304 physical copies of P_small (all sectors with m >= 1, weight spaces of any size w >= 1), ranks with the 1e-30 relative threshold at 50 digits. Pairs = (rank over rows with the complete m12-range, rank over all rows; missing patterns = 0).

* Ψ[κ][m12]: (1,1) 158, (1,2) 111, (1,3) 7, (2,2) 12, (0,1) 4, (0,2) 8, (0,3) 4. Almost all of the 158 (1,1) are trivial supports (one nonzero row or column). Restricting to copies with >= 2 nonzero rows and >= 2 nonzero columns (146 copies): **rank 1 in only 4**, rank 2 in 131, rank 3 in 11.
* Ψ[m13][(m23,m12)]: (1,1) 119 (trivial), (1,2) 85, (1,3) 3, (0,1) 3, (0,2) 74, (0,3) 20; the 182 non-degenerate copies have rank 2 (159) or 3 (23), **never 1**.
* Ψ[m23][(m13,m12)]: (1,1) 138 (trivial), (1,2) 87, (1,3) 12, (0,1) 2, (0,2) 58, (0,3) 7; the 164 non-degenerate copies have rank 2 (145) or 3 (19), **never 1**.

So ψ(m12, κ) = f(m12) g(κ) fails for essentially every copy with a non-trivial support, and likewise for the two other cuts: slices are not separable (A5 negative).

## 7. Surprises and corrections

* **Brief / plan:** (i) the plan says "a path in 29 of 44" etc.; here the right invariant is the counting condition n_{d+1} <= 3 n_d plus rank; (ii) of the 62 gradings only 49 are distinct; (iii) the Vandermonde formula holds only for coherent orientation sets (+1 / −1), not for arbitrary choice of one target per spectator; (iv) the E12 cut is void (W1 = W) on 171 of 222 spaces, so the "E12-cut W1" is a nontrivial restriction only for 51.
* n0 = 1 is **not** the generic situation for large spaces (6 of 23 for w >= 9), against the hope; the generic case needs the polynomial-matrix formulation, which works exactly (§5).
* The degree of P_x equals the level exactly and equals G4's minimal Krylov degree, so the recursion reproduces the optimal polynomials, not just some polynomials.
* J (the C_j) is radical in all tested sectors; dim Q[v]/J = dim W1 exactly.
* No clean bottom = cut-end rule: for the same weight space both +m12 and −m12 are usually fine; failures of ±m12 are counting failures at the first step.

## 8. Not done / limits

* M~ was only used as a spectral comparison in A3 (counts of triples), not graded in A1 (as agreed in the brief); no recursion on M~.
* Symbolic C(v) minors only for 5 n0 = 2 sectors with w <= 8; no symbolic n0 = 3 minors; Gröbner dimension for w1 <= 10 only (w1 = 12, 15 numerically only).
* The "A → grading" rule is empirical on 222 spaces; no proof, no closed-form coefficients for the recursion (the entries of F_d are G4's closed forms restricted; this is for G8d).
* A3 uses sectors with an n0 = 1 grading (biased towards small w); the large-w sectors are covered by A4 numerically.
* Products beyond P_small ∪ P_big not run. The residual magnitudes quoted for A3 are absolute with unnormalised primitive-integer polynomials (coefficients up to 1e7).

## 9. Usage

```
cd Paul/Python/gaudin_gl4
python -X utf8 lattice_recursion_gl4.py all            # 21 s: writes lattice_recursion_gl4_output.txt and lattice_recursion_gl4_tables.txt
python -X utf8 lattice_recursion_gl4.py run recursion,charpoly
python -X utf8 test_lattice_recursion_gl4.py           # 3 s
```
```python
import lattice_recursion_gl4 as L, slice_ops_gl4 as so
cp = so.copies(lam, lamp, mu)[0]
r = L.slice_from_triple(lam, lamp, mu, cp['ev'], grading='best', space='W1')   # r['slice'] == cp['slice'] to 1e-43; r['resid'] = unused equations
sec = L.sector(lam, lamp, mu); sy = L.sym_recursion(sec, 'W1', (-1, 0, 0))    # sy['psi'][pos][0] = P_x ; sy['unused'] = C_j
```
