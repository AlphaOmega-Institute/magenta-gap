# The Magenta Gap — canonical working manuscript

[`Magenta_Gap.tex`](Magenta_Gap.tex) is the sole active manuscript, with a matching [`Magenta_Gap.pdf`](Magenta_Gap.pdf). It uses the standard LaTeX `article` class in a single-column format and remains a public work in progress. A venue-specific layout can be applied later without maintaining a second scientific text.

The reviewed mathematical revision was merged at `51b37d5`. The format conversion retains its complete arguments, three appendices, 16 bibliography entries, numerical table, three figures and assistance disclosure. It restores explicit equation-name prefixes previously supplied by the Optica class and clarifies one forward reference. No perceptual encoding or new empirical claim is added.

The November 13 conceptual-restoration draft expands the Introduction and Discussion around the transformation problem: relating stimulus distributions, observer state and context to organized color appearance. It restores the motivation in corrected form, with sourced blackness/afterimage examples and a qualitative luminosity-threshold comparison. It keeps stimulus gamut, appearance target and regime-specific image distinct. The initial draft at `cee2c7b` had 22 pages and 22 references, with its mathematical content and figure files unchanged from `5c4d9e6`.

The follow-up restores the separate color-body/achromatic-quotient figure, the bee-purple comparison, and the ten historical and philosophical references omitted from July, adding Chittka's color-hexagon paper for the bee example. The current draft has 24 pages, four figures and 33 references. All 50 displayed equations, theorem/proof environments and the numerical table are unchanged. It incorporates Claude's standalone figure repairs from `e07eb71`, removing cropped-away July text, and retains the drawing sources. It introduces no new psychophysical data or completed encoding model. This revision is prepared for review before merging.

The [July manuscript in `jmp/`](jmp/) is a frozen historical baseline, not a second active paper. Its source and original PDF are unchanged. The reviewed Optica-format September/November text, PDF and former class/support files remain recoverable in [Git history at `51b37d5`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/51b37d506dee66de9228a1ee2c0e81aa7a193488/manuscript).

## Build

Run from this directory with a complete TeX installation providing pdfLaTeX, latexmk and the standard packages named in the preamble:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error Magenta_Gap.tex
```

For Overleaf, upload `Magenta_Gap.tex` and the four figure PDFs below, set `Magenta_Gap.tex` as the main document, and select pdfLaTeX. The arguments and bibliography are in the TeX source; no external companion note, BibTeX database, custom class or local compatibility shim is required.

- `fig1_magenta_gap.pdf`
- `fig2_color_body.pdf`
- `fig3_spindle.pdf`
- `fisher_bound_convex_hull.pdf`

The three conceptual figures also have standalone TikZ sources beside their PDFs. Those sources are optional for the manuscript build; rebuilding them requires the standard `standalone` class and the packages listed in each figure source. Keep `\noindent` at the start of the abstract in subsequent revisions.

The former Optica class and its unused support/shim files are removed from the current working tree. They remain in Git history. The canonical PDF is freshly built and visually inspected; the verification record is in [EDITORIAL_NOTES.md](EDITORIAL_NOTES.md).

## Calculations and scope

The paper distinguishes conditional null-boundary compatibility, an exact effect-algebra connection, and a certified obstruction for one separately specified bounded-reflectance model. Its Fisher illustration uses an uncalibrated convention. The stimulus encoding, report procedures and physical controls needed for a perceptual realization remain open.

The corrected calculation archive is now in [`recompute/`](../recompute/), merged from [`claude/corrected-bound-calculation`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/a0282ae23d11c102550c88a9f44a8e2a99f99211/recompute), including code, input arrays, results and the numerical figure. The source commit remains the pinned reproduction target cited in the manuscript. The repository's top-level scripts are still the older pipeline.

The exact continuum certificate is in [`research/object-colour/`](../research/object-colour/), and the 11 effect-algebra checks are in [`research/colour-interface/`](../research/colour-interface/). The original September numerical table and figure remain unchanged. The independent reproduction archive is preserved on [`claude/independent-colour-checks`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/581ff6880effe636dfd32b781375b7a82d2d741c/research/independent-checks); all 24 executions were reviewed and their scientific output reproduced. Most of those scripts print numerical checks rather than assert formal pass/fail conditions.

## Provenance

The working manuscript descends from `Magenta_Gap_JOSAA_Source_and_Calculations_2026-09-28.zip`, imported November 11, followed by the November 12 consolidation and review. The original imported TeX and PDF remain recoverable at `cb997e8`.

Original imported TeX SHA-256 (historical, not the current file):

```text
1198503e699c84206342c52c2561ee9d9ba6e081e3481ff5c647e53364f3dabe
```

`PROVENANCE.json` records that import and is not a checksum manifest for the current files. The editorial notes retain their historical wording and append subsequent decisions and verification. Git history preserves the rename from the dated source to the canonical filename.
