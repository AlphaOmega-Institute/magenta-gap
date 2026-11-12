# The Magenta Gap — September revision

The longer July JMP manuscript is available separately in [`jmp/`](jmp/), with its verified recovered TeX and original 51-page PDF. The files in this directory remain the shorter September revision.

This is the September 28, 2026 author-review draft, imported on November 11 and consolidated on November 12, 2026. It remains work in progress, not a claim of submission readiness or a frozen publication. The dated filename is retained so the shared TeX has one version history.

- [TeX source](Magenta_Gap_JOSAA_2026-09-28.tex): the revised self-contained text and bibliography.
- [PDF](Magenta_Gap_JOSAA_2026-09-28.pdf): a fresh build of the revised source.
- [Editorial notes](EDITORIAL_NOTES.md): the original September notes followed by the November consolidation record.

The TeX contains the complete text, mathematical arguments and bibliography. The three existing figure PDFs and layout support files are included. The effect-bicone algebra, three height conventions, conditional altered-colour constraint, paired-carrier readout limitation and certified object-colour obstruction are now incorporated. The exploratory Hilbert-gamut comparison and the separate cubic object-colour theorem are not inserted. No new perceptual encoding or geometry is selected.

## Build

Run from this directory with pdfLaTeX and latexmk installed:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error Magenta_Gap_JOSAA_2026-09-28.tex
```

The November draft is freshly compiled and its pages are visually checked. The original import's missing-format build limitation was resolved locally using the installed TeX sources and standard font maps; no local runtime paths are required by the manuscript source. A complete TeX installation can use the command above.

The package retains the prior Optica length-check layout, including local compatibility shims `algorithm.sty`, `algpseudocode.sty`, and `multibib.sty`. These cover unused interfaces loaded by the class, not full replacements for the standard packages. The original source-archive documentation describes their provenance; they may be replaced by the standard packages on a complete TeX installation.

## Calculations and scope

The corrected calculation archive is already available on [`claude/corrected-bound-calculation`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/a0282ae23d11c102550c88a9f44a8e2a99f99211/recompute), including its code, input arrays, results, and numerical figure. It has not been merged into `main` by this import. The top-level scripts on `main` are the older pipeline.

The manuscript distinguishes:

- a conditional null-boundary compatibility theorem;
- a numerical illustration using a conventionally scaled Fisher metric;
- the additional calibration and control assumptions needed for a physiological interpretation.

The November consolidation also records a different, continuum bounded-reflectance model. Its exact certificate is in `research/object-colour/`, and the effect-algebra checks are in `research/colour-interface/`. The original September numerical table and figure are unchanged. The separate independent reproduction archive remains on `claude/independent-colour-checks` at `581ff68`; all 24 executions were reviewed and their scientific output reproduced. Most of those scripts print numerical checks rather than assert formal pass/fail conditions.

The September editorial notes retain their historical wording, including statements about what had not yet been published at that time; the appended November section records the present changes.

## Provenance

Imported from `Magenta_Gap_JOSAA_Source_and_Calculations_2026-09-28.zip`. At import, its TeX was identical to the separately supplied `Magenta_Gap_JOSAA_2026-09-28.tex`; the current source includes the November consolidation and review corrections recorded in the editorial notes.

Original imported TeX SHA-256 (historical, not the revised file):

```text
1198503e699c84206342c52c2561ee9d9ba6e081e3481ff5c647e53364f3dabe
```

The original imported source and PDF remain recoverable at commit `cb997e8`. `PROVENANCE.json` records that import and is not a checksum manifest for the revised files. Git history records the November edits. The July manuscript, research certificate and input data are unchanged by this revision.

