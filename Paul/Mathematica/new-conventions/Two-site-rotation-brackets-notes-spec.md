# Spec: write-up of the two-site rotation-bracket thread (Stages A, B, D, C of `Two-site-rotation-brackets-plan.md`) — `Two-site-rotation-brackets-notes.{md,tex}` (2026-10-11)

Delegated to a sonnet general-purpose agent (write-only within the two files named below; read-only everywhere else;
no notebook edits, no `paper/*.tex`, no HANDOFF, no git). Model: the earlier thread write-ups
`Site1-rotation-brackets-notes.tex` and `Untwist-two-rotation-notes.tex` (same structure: setting and conventions, results
with evidence tags [exact]/[numeric, 10⁻ⁿ], notebook + cell for every claim, an appendix of source inconsistencies, a
one-page summary first). Compile the .tex with `pdflatex` if available (twice), else leave the .tex and say so.

## Sources (read all; cell numbers must be quoted from the notebooks, not from memory)

- The plan `Two-site-rotation-brackets-plan.md` and the four stage specs: `L2-site2-rotation-brackets-spec.md`,
  `L2-two-site-rotation-brackets-spec.md`, `L2-two-vacua-untwist-spec.md`, `L2-210-two-vacua-brackets-spec.md`.
- Notebooks (read via the JSON files or `wolfbook_getNotebookContext`; markdown result cells hold the numbers):
  `L2-site2-rotation-brackets.wb` (36 cells, Stage A), `L2-two-site-rotation-brackets.wb` (58 cells, Stage B + follow-up
  Task 4b in cells 47–49, summary 58), `L2-two-vacua-untwist.wb` (40 cells, Stage D, summary 40),
  `L2-210-two-vacua-brackets.wb` (Stage C, adjoint; summary = last cell). UTF-8 dumps of the first three are in the
  scratchpad (`ref-site2-all.txt`, `ref-two-site-all.txt` (55-cell version), `ref-two-vacua-untwist-all.txt`).
- HANDOFF.md top section "NEWEST (2026-10-11)" (the main session's running summary and readings — treat its
  interpretations as the main session's, the numbers as the notebooks').
- Background: `L2-210-site1-rotation-brackets.md`, HANDOFF's "2026-10-10, late night" section (site-1 general results),
  `Two-rotations-det-overlap.md` §1, §4, `L2-210-untwist-CG.md` §1–§3, `L2-210-localising-span-spec.md` §6.

## Structure (sections; keep the maths in real math mode, Greek as \alpha etc., never Unicode in the .tex)

0. **Summary (one page).** What was asked (Paul: two rotations on both sites / a different site-2 reference, bracket count
   polynomial in the boxes, untwist with w, w′ → 1 independently), what was found: (i) the [3,1,0] 5-bracket locus is the
   one-site vacuum family (Ωc = e₁ᵀρ(g), α = 1/w₁); (ii) the SoV construction is site-asymmetric and the asymmetry is the
   Yangian R-matrix conjugation (undressed site-2 rotation 17/40, R-dressed 5/7 with cc₁(−t)); (iii) the two-vacuum family
   Ω_w⊗Ω_{w′} untwists cleanly — every CG datum of adj⊗adj and [3,1,0]⊗[3,1,0] at leading order from the overlap alone,
   K/S = O(ε) for all branches (first gl(3) family with this property); (iv) its bracket cost is generic (22/54), the det₆
   costs 18/40 in product units, the CG-complete (∞,∞) sub-family costs 12 (adjoint, cheap) / 45 ([3,1,0], not cheap);
   (v) the adjoint (∞,∞) family in closed form and its Q-side CG (Stage C results); (vi) the open decision: cost metric /
   determinant form for [3,1,0].
1. **Setting and conventions** (companion vs diagonal frame, labels, monomials, blocks A^{(c,c′)}, the two counts r and P
   with the lower/upper-bound statement, regime).
2. **Stage A** — vacuum-family identification, closed forms in w, torus image of the adjoint quadric, the site-2 mirror
   (supports, counts, loci, cc₂(t) = cc₁(1−t), reach, closed forms on b₃ = 0).
3. **Stage B** — two-site K tensor; R-matrix (construction in the commutant, checks, dressed support/count/loci); block
   tables (a)–(h) and Task 4b; det₆ regression (and why Σr ≤ 20 was not a valid bound); 4-parameter family and
   sub-families; one-site-fixed loci; the conclusion that cheapness correlates with loss of reach, EXCEPT the (∞,∞) family
   with w₁ ≠ w₁′ (complete, 12/45).
4. **Stage D** — exact undo, closed-form ε = 0 Gauss parameters, controls (pinned w₂′ does not degrade; common w₁ reads
   one irrep), reach table of the sub-families, per-weight table, K/S table, CG at leading order with the reference numbers,
   the four-overlap table and the exchange-parity sign.
5. **Stage C (adjoint (∞,∞))** — closed-form F, G, count, t-laws, u-forms, Q-number verification, ε-exponents, six CG data
   from the Q side; [3,1,0] greedy result. FILL FROM THE STAGE C NOTEBOOK.
6. **What it means and what is open** — for Paul: the cost-metric question (products of 2×2 bracket determinants vs
   determinants of brackets: det₆ = 36 entries vs 40 products), the determinant-form question for ⟨Ω_w⊗Ω_{w′}| with
   non-vacuum bra data (the two-site K tensor makes a linear-algebra search possible), the R-dressed two-site covectors,
   the interpretation of cc₁(1−t) vs cc₁(−t), and the [3,1,0] negative. Keep this section short and factual; mark the main
   session's readings as such.
- **Appendix A: cell maps** per notebook. **Appendix B: numeric points and grids.** **Appendix C: source inconsistencies**
  (e.g. the spec's wrong 229/342; Stage D's misattribution of the 8/14 line to w₁ ≠ w₁′; the Stage B wedge-count 3+3+2+2 vs
  ΣP = 12, as resolved by Stage C; any others you find — list them, don't fix them).

## Rules

Every number in the notes comes from a notebook cell (cite notebook + cell); when two sources disagree, report both in
App. C. Evidence tags on every claim. No speculation beyond §6, and there flagged. Length: aim for 15–25 pages of .tex.
Write the .md as a faithful plain-text companion (same sections, ASCII math as in the earlier .md write-ups).
