# The Magenta Gap — September revision

The longer July JMP manuscript is available separately in [`jmp/`](jmp/), with its verified recovered TeX and original 51-page PDF. The files in this directory remain the shorter September revision.

This is the September 28, 2026 author-review draft, imported on November 11, 2026. It is work in progress, not a claim of submission readiness or a frozen publication.

- [TeX source](Magenta_Gap_JOSAA_2026-09-28.tex): imported byte for byte from the author's September draft.
- [PDF](Magenta_Gap_JOSAA_2026-09-28.pdf): the compiled copy supplied in the September source archive.
- [Editorial notes](EDITORIAL_NOTES.md): the dated account of the September corrections and outstanding review decisions, preserved unchanged.

The TeX contains the complete text and bibliography. The three figure PDFs and the original layout support files are included. The recent Hilbert-geometry comparison, matrix-effect interpretation, and Berthier–Provenzi citations have not been inserted.

## Build

Run from this directory with pdfLaTeX and latexmk installed:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error Magenta_Gap_JOSAA_2026-09-28.tex
```

The import was checked for source identity and figure dependencies. A fresh build was attempted but the local TeX installation lacked `pdflatex.fmt`; no fresh compilation is claimed. The included PDF is preserved from the source archive.

The package retains the prior Optica length-check layout, including local compatibility shims `algorithm.sty`, `algpseudocode.sty`, and `multibib.sty`. These cover unused interfaces loaded by the class, not full replacements for the standard packages. The original source-archive documentation describes their provenance; they may be replaced by the standard packages on a complete TeX installation.

## Calculations and scope

The corrected calculation archive is already available on [`claude/corrected-bound-calculation`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/a0282ae23d11c102550c88a9f44a8e2a99f99211/recompute), including its code, input arrays, results, and numerical figure. It has not been merged into `main` by this import. The top-level scripts on `main` are the older pipeline.

The manuscript distinguishes:

- a conditional null-boundary compatibility theorem;
- a numerical illustration using a conventionally scaled Fisher metric;
- the additional calibration and control assumptions needed for a physiological interpretation.

The September editorial notes retain their historical wording, including statements about what had not yet been published at that time.

## Provenance

Imported from `Magenta_Gap_JOSAA_Source_and_Calculations_2026-09-28.zip`. Its TeX is identical to the separately supplied `Magenta_Gap_JOSAA_2026-09-28.tex`.

TeX SHA-256:

```text
1198503e699c84206342c52c2561ee9d9ba6e081e3481ff5c647e53364f3dabe
```

The dated filename identifies the imported draft. Git history records subsequent edits.
