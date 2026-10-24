# Questions on SPEC-G9-fast-cg-ubuntu.md (2026-10-24, Ubuntu laptop)

Prerequisites check passed on Ubuntu: HEAD = 05dc470 (round G8), environment matches the brief (i5-1345U, 12 threads,
15 GB, g++ 11.4.0; Python 3.10.12, numpy 2.2.6, scipy 1.15.3, sympy 1.14.0, mpmath 1.3.0), every referenced function
exists where stated, the G7 Munich cache is present (190 `hw_*.txt`, orig/eig/compact timings, 7 `bench_eig` files).
Three points to clarify before dispatching:

1. **Weight-space counts.** §0 says G8 was established on "222 weight spaces" but that `closed_form_mats` is "exact on
   all 328 weight spaces"; T5 says "W₁ = W, 171 of 222 spaces". Which test set is the 222, which the 328? (One sentence
   would stop the agent flagging it as an inconsistency.)

2. **Basis of the W₁ matrices.** Route B applies the reduced words with the *orthonormal* CSR generators (`float_ops()`),
   while the slice output is in the plain basis (u_c/√g_c) and T1 compares against G3's `weight_space_mats` (plain basis?).
   Please state explicitly: are V₃, V₄ᵃ, V₄ᵇ on W₁ and the string test s(u) taken in the orthonormal basis, and is the
   T1 twin comparison to be done after conjugating by √g? (Otherwise a metric-factor mismatch could look like a failure.)

3. **Untracked binaries.** `Paul/Cpp/AlexCGC/cgc_hw_prof` and `cgc_hw_prof_ob` are untracked build products on Ubuntu.
   Should they be gitignored (presumably), and should the G9 agent leave them alone?

---

# Answers (2026-10-24, Windows session; the brief has been patched in the same places — §0, §1 S1/S2, T1, T3 — so the agent can read the brief alone)

1. **Counts.** Both numbers refer to the same 12 test products (P_small ∪ P_big, with [3,2,1,0]⊗[4,2,1,0] in both orders). Verified just now with
   `glN_reps.candidates`: the products have **328 weight spaces W = V_λ[μ−λ′] with w ≥ 1, of which 222 have w ≥ 2**. G8a/G8c/G8b's structural statements
   (gradings, recursion, path tests, "W₁ = W in 171 of 222") are about the 222 with w ≥ 2, where they are non-trivial; G8d's `closed_form_mats` check ran on all 328
   including the 106 one-point spaces. (G8c reports 221: it dropped one space with m = 0.) No inconsistency.

2. **Basis.** Fixed as follows, and now written into the brief:
   * G3's `weight_space_mats` (and everything exact in G1–G8) are in the **plain** GT basis: M[row][col] = coefficient of e_row in V e_col.
   * `GLNRep.float_ops()` are the generators in the **orthonormal** basis ẽ_c = e_c/√g_c (Eh[r,q] = E[r,q]·√(g_r/g_q)). So Route B produces the operators on W₁ in
     the orthonormal basis: M̃ = D M D⁻¹ with D = diag(√g_c), c ∈ W₁ (entrywise M̃[r,c] = M[r,c]·√(g_r/g_c)). **Work in the orthonormal basis throughout S1–S2**
     (matrices, the eigenproblem, the string test); Route A (closed forms, plain by construction) must be conjugated by D before use or comparison.
   * **T1 twin comparison: conjugate G3's exact plain matrix by √g** (M̃_ref[r,c] = float(M[r,c])·√(g_r/g_c), g from `R.g`) and compare with Route B in the
     orthonormal basis, relative to max |M̃|. Eigenvalues are basis-independent, so the triples compare directly.
   * **Slices:** an eigenvector ũ of the orthonormal-basis matrix has plain components u_c = ũ_c/√g_c (v = Σ u_c e_c = Σ ũ_c ẽ_c). G3's canonical normalisation
     (first nonzero component = 1, ascending GT index) is applied to the PLAIN components, and `se.hw_from_slice` takes the PLAIN slice (it multiplies by `sg1[W]`
     internally, i.e. it converts back to ũ) — so pass u, not ũ.
   * **String test:** compute s(u) = max_i ‖E_{i,i+1}^{k_i} ũ‖₂ / ‖ũ‖₂ in the orthonormal basis with the orthonormal CSR powers (apply E_{i,i+1} k_i times to ũ; the
     intermediate weight spaces are different from W₁, so do not restrict to W₁ during the power). Note that G4/G8b measured their "spurious violate by ≥ 16, physical
     by 0" in the PLAIN basis with the max norm at 50 digits (`structure_gl4.string_violations`), so **re-measure the gap in your basis**: report min over spurious and
     max over physical of s per product; the acceptance is a gap ≥ 10⁶ with τ = 1e-8 placed inside it (adapt τ and say so if the gap is smaller). The physical
     count must equal Racah–Speiser in every sector regardless of τ.

3. **Binaries.** Yes: `Paul/Cpp/AlexCGC/cgc_hw_prof` and `cgc_hw_prof_ob` are now in `.gitignore` (this commit; `clebsch`, `cgc_hw_bench` were already ignored).
   Leave them alone — do not rebuild, delete or commit them; the G9 benchmark re-runs Munich only for a missing cached file (G7 rule) and uses these binaries as they are.
