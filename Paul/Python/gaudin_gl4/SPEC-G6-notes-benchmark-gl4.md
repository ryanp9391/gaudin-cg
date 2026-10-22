# Spec G6: the benchmark table and the write-up — hep-th-style notes on the gl(4) slice engine and the two-site Gaudin missing labels (2026-10-22)

Delegated to a general-purpose subagent (SONNET). **No Wolfram kernel, no wolfbook tools, no notebooks, no git.** Work package WP-G6 of
`Paul/Mathematica/Gaudin/gl4-slice-engine-plan.md`. Run anything Python from `Paul/Python/gaudin_gl4/` with `python -X utf8` (≤ 2 heavy processes; you should not need more
than a few re-runs of existing CLI stages — the numbers exist in the RESULTS files). LaTeX: MiKTeX is on PATH (`pdflatex`); compile with `pdflatex` twice (+ `bibtex` if you use
a .bib; a `thebibliography` environment is fine); if MiKTeX wants to install a package on the fly and hangs, use only standard packages (amsmath, amssymb, booktabs, hyperref,
geometry, graphicx). **Never edit existing files** (in particular do not touch `Paul/Mathematica/Gaudin/Slice-engine-notes.*`, `Gaudin-sl3-CG-notes.*`, `paper/*.tex`, the plan, or any
RESULTS/SPEC). Files you own: `Paul/Mathematica/Gaudin/gl4-slice-engine-notes.tex`, `.md`, `.pdf` (and the LaTeX by-products), `Paul/Python/gaudin_gl4/RESULTS-G6-notes-benchmark-gl4.md`,
`Paul/Python/gaudin_gl4/benchmark_gl4.md` (the standalone benchmark table). Another agent (G4, structure) is running concurrently and may finish before you: check for
`RESULTS-G4-structure-gl4.md` before you write the structure section; if it is absent, write that section from G3 §7 and mark the G4 items as pending.

Read first (in this order): this file; the plan (all of it); the gl(3) notes `Paul/Mathematica/Gaudin/Slice-engine-notes.md` (whole; the style, structure and tone to match — plain
hep-th prose, short sections, every claim tied to a check, honest caveats, a reproduction appendix) and its `.tex` (for the LaTeX conventions used in this repo: document class,
macros, how tables and references are done); `Paul/Mathematica/Gaudin/Literature-check-missing-label.md` §0, §2, §3, §10 (what is classical, the must-cite list, the recommended
framing); then the four reports: `RESULTS-G1-slice-engine-gl4.md` (§0, §1, §4.3, §5, §6, §7, §9, §10, §11), `RESULTS-G2-gaudin-gl4-ops.md` (§0, §2, §4, §5, §6, §9),
`RESULTS-G3-slice-ops-gl4.md` (§0, §1, §3, §6, §7, §8, §9, §11, §12), `RESULTS-G5-oper-gl4.md` (§0, §2.1, §3, §4, §5, §6, §7); and `RESULTS-G4-structure-gl4.md` if present. Also
`Paul/Python/gaudin_sl3/RESULTS-slice-engine.md` §5 (the gl(3) benchmark numbers), `Paul/Cpp/AlexCGC/README.md` and `results/hw_timings.txt` (Munich), `HANDOFF.md` top block
(session summary and the ndsu3lib block for the gl(3) context: "orders of magnitude" holds only against general-N codes).

## 0. What the document is for

Paul will read it to decide what goes into the paper; a reader who knows representation theory and integrable systems but has not followed the project should be able to
follow it without the code. It is the gl(4) companion of `Slice-engine-notes` (gl(3)). It must be honest about three things: (i) what is a theorem (Kostant/PRV slice theorem;
the dictionary of G5 given G2's identities; the reduction rule), what is an exact verification on tested cases (everything else: commutativity, simple joint spectrum,
BGG intertwining, Weyl/BGG spectrum identity if G4 confirms it, closed forms), and what is a benchmark on one machine; (ii) that the speed comparison is against the only
general-N code (Munich, arXiv:1009.0437) and on a different machine than Munich's own table, and that the gl(3)-specific lineage (Draayer–Akiyama, ndsu3lib) does not apply
to N = 4; (iii) what is classical (slice theorem; Gaudin algebras as missing labels: Vinberg, Rybnikov, Feigin–Frenkel–Rybnikov, Mukhin–Tarasov–Varchenko, Halacheva–Kamnitzer–
Rybnikov–Weekes; the oper/Bethe dictionary; for N = 3 the operator is Pluhař–Smirnov–Tolstoy / Crampé–Poulain d'Andecy–Vinet's X) and what is new here (the explicit
engine with N−1 string conditions and the level recursion; the exact cdet triple with its t-independence and structure; the failure of the cubic operator alone and the
parity argument; the slice reduction rule with normal ordering; the canonical copies; the forced dictionary with T_k/τ_k; the exact Q-system solve showing all four degree
conditions are needed; the BGG intertwining and spectral-subspace facts; the benchmark; the adjoint+gauge identification of the July-2026 ODE).

## 1. Deliverables

**G6.1 `benchmark_gl4.md`** — one standalone page: (A) the SU(4) table [S,S,0,0]², S = 1..12 (+14, 16, 20): Munich total hw-solve (its machine, i5-1345U; n/a beyond S = 8),
our t_a / t_b (G1, i7-12700KF), ratio, entries written, worst residual; with Munich's per-irrep singlet row at S = 8 (18.7 s vs 4.5 ms warm / 45 ms cold). (B) the multiplicity
families (G3 §9): [2a,a,0,0]², a = 1..6 (+7, 8 slices only) and [3a,2a,a,0]², a = 1..4: #irreps, Σ copies, max m, max w, max n_W, t_a (exact) and t_a with float64 eigen-decomposition,
t_b; the breakdown (reduction once per λ′ ≈ 4 ms; exact operator application; eigen; recursion). (C) the gl(3) context in three lines (from `RESULTS-slice-engine.md`: [14,7,0]² in
0.10/0.17 s vs Munich 206 s / 10.9 s repaired; ndsu3lib caveat). (D) a paragraph "how to read it": what each column measures, which steps are n_W-free, what is exact and what is
float, the machine difference, the Python-overhead regime (S ≤ 2), scaling exponents (Munich ×5–7 per unit S, ours ×1.4–2). Every number must be traceable to a RESULTS file
(cite section numbers in a final "sources" line).

**G6.2 `gl4-slice-engine-notes.tex` + `.md` + compiled `.pdf`** — the notes, ≈ 12–20 pages, structure (adapt, do not pad):
1. Punchline (half a page): the HW CG problem of gl(4) is a w-dimensional slice problem; the three missing labels are the Laurent coefficients of the two-site Gaudin cdet; the
   engine does Munich's 80 s benchmark in 0.09 s and reaches sizes Munich cannot; the cubic label alone fails for N = 4, the triple works; the oper dictionary is forced.
2. The problem and the standard algorithm (Munich's null space, cost O(n_W³)), with the SU(4) numbers.
3. The slice theorem for gl(N) (Kostant 1959 / PRV 1967, N−1 string conditions, site roles), the level recursion (A_a, left inverses, exact twin, the dangling conditions needed for
   N = 4 — G1 §4.3), the brute-force cross-check, dim M_μ = multiplicity on 26 products.
4. The missing-label operators: the column-ordered cdet, Laurent data, homogeneity in t and the triple; commutativity (what was checked); the structure identities (V₃ = ΣD(2,1) +
   ((N−2)/2)η, V₄ᵃ, V₄ᵇ, the trace forms, naive traces fail); self-adjointness and real spectrum; **V₃ alone is degenerate for N = 4** with the exchange-parity argument
   (PV₃P = −V₃ + c, ≥ |p₊−p₋| zero modes; forced vs accidental per G3 §8.1), {V₃,V₄ᵇ} separates everything tested, the triple as the canonical definition; the N = 3 specialisation
   V₃ = V + η/2 (so the notes are consistent with the gl(3) notes).
5. Reduction to the slice: the rule (normal ordering; lowering kills, Cartan → λ′, raising → −e¹ reversed), the N = 3 test, the direct-route verification, universality in λ′;
   the structural facts of G3 §7 (commute on the whole weight space for every λ′, BGG intertwining, M_μ a spectral subspace) with the Verma-module interpretation stated as an
   interpretation; the canonical copies (joint diagonalisation, ordering, normalisation, G-orthogonality, exchange pairs).
6. The oper side (G5): the dictionary table with T_k, τ_k (state it as a theorem given G2's identities), exponents and polynomial solutions of degrees μ_j + 4 − j, the Wronskian,
   the exact Q-system solve (Σ = joint spectrum, 180 sectors, m ≤ 12; all four degree conditions needed — the gl(4) pencil-bug lesson), μ₄ = 0 ⇒ V₄ᵃ = V₄ᵇ = 0, the July-2026
   `gaudinODE` = 8Π²Π^{(S,P)}D†Π^{−(S,P)} with the Casimir dictionary.
7. Structure on the GT lattice (G4 if available, else G3 §7 + "pending"): stencils, closed forms, the Weyl/BGG characterisation of the physical spectrum, chain/polynomiality.
8. Benchmark (the table of G6.1 in condensed form) and caveats (machines, Python, exact vs float, what was not run: the literal Munich vector comparison — the GT bases coincide,
   G1 §9).
9. What is new, what is classical, what is next (proofs: commutativity for all λ, λ′; the Weyl identity; closed forms; a compiled operator application; gl(N); the paper).
Appendices: A. conventions (GT order, metric, sites, slice, the term-list format); B. the reduced universal words of V₃ (39 words) in a table if it fits, else the counts;
C. reproduction (every CLI stage of G1–G5 with its runtime, the reference files and what each contains).
References: use the must-cite list of the literature check (Kostant 1959; PRV 1967; Alex–Kalus–Huckleberry–von Delft 2011; Vinberg 1991; Rybnikov 2006/2020; Feigin–Frenkel–Rybnikov
2010; Mukhin–Tarasov–Varchenko 2009; Halacheva–Kamnitzer–Rybnikov–Weekes 2020; Talalaev 2006 / Chervov–Talalaev for the cdet; Crampé–Poulain d'Andecy–Vinet 2023 and Pluhař–
Smirnov–Tolstoy 1986 for N = 3; Asherova–Smirnov–Tolstoy extremal projectors; Zhelobenko's book) — cite only what you can attribute from the literature check or the gl(3) notes;
do not invent references.
The `.md` is a faithful Markdown rendering of the `.tex` (same text, Unicode math where readable), as for `Slice-engine-notes.md`.

**G6.3 `RESULTS-G6-notes-benchmark-gl4.md`** — short: files, how the numbers were sourced (section pointers), what you re-ran, the compile log summary (0 errors; warnings listed),
a list of statements in the notes that are verifications rather than theorems (so Paul can see the status of each claim at a glance), a 10-line paragraph proposing how the
gl(4) material could enter `Slice-engine-notes` or the paper (Paul's decision; do not edit them), and anything inconsistent you found between the RESULTS files (report, do not fix).

## 2. Rules of writing

Plain sentences, short sections, every number with its source and its machine; no superlatives; say "verified on N sectors" rather than "proved" unless it is a proof; one idea
per paragraph; tables for numbers; the notation of the RESULTS files (V₃, V₄ᵃ, V₄ᵇ, S_{k,j}, t = θ₁₂, M_μ, w, m, n_W, T_k, τ_k). No code in the prose; the CLI goes to the appendix.
If RESULTS files disagree on a number, quote both and flag it. Write the `.tex` with the file tools; compile; read the log; fix; recompile.
