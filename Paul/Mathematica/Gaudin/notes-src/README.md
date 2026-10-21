# notes-src: generator of `Gaudin-sl3-CG-notes.{tex,md}`

The notes `../Gaudin-sl3-CG-notes.tex` / `.md` / `.pdf` were written on 2026-10-16 from the results of rounds 1-4 of the sl(3) Gaudin project.
The `.tex` is ordinary, hand-editable LaTeX and may be edited directly; this directory is only the (optional) pipeline that produced it.

* `src/*.src` - the hand-written prose in a light markup (Markdown-like headings and lists, inline LaTeX math, `^[...]` footnotes,
  `$$ ... $$` equations, `@table label | caption` ... `@endtable` tables, `@include gen/...`). Edit these if you want to regenerate both twins.
* `gen/*.src` - the data tables (CG slices, untwisted Laurent data, epsilon-series, one-site kernel, adjoint kernel `G`, base-point kernel `K` of [3,1,0]).
  They are **produced, not typed**, by `gen_data.py` from the exact data files `Paul/Mathematica/Gaudin/ref/untwisted-exact-*.txt`,
  from `Paul/Python/gaudin_sl3/RESULTS-twovac-limit.md` (the adjoint kernel, parsed) and from the exact polynomial laws stored in
  `Paul/Python/gaudin_sl3/general_kernel_rep_310.pkl` (read-only; the `.pkl` caches are not in git).
* `build.py` - renders `src/*.src` to `Gaudin-sl3-CG-notes.tex` and the Markdown twin (LaTeX math converted to Unicode; formulas that do not convert cleanly stay as inline LaTeX).

## Regenerate

    cd Paul/Mathematica/Gaudin/notes-src
    python -X utf8 gen_data.py            # refresh gen/*.src from the data files (needs the .pkl for the [3,1,0] kernel table)
    python -X utf8 build.py ..            # writes ../Gaudin-sl3-CG-notes.tex and .md
    cd .. && pdflatex Gaudin-sl3-CG-notes.tex && pdflatex Gaudin-sl3-CG-notes.tex   # twice for references

`build.py <outdir>` overwrites the `.tex`/`.md` in `<outdir>`, so do not run it after editing the `.tex` by hand unless you want to discard those edits.
Python 3.13, standard library only for `build.py`; `gen_data.py` uses only the standard library plus `pickle` (the `.pkl` contains mpmath objects, so mpmath must be installed to load it).
Note: write `.src` files with an editor, not with shell heredocs (backslashes).
