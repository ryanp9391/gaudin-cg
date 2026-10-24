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
