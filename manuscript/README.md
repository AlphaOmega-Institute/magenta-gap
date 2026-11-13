# Magenta Gap — working manuscript

[`Magenta_Gap.tex`](Magenta_Gap.tex) contains the complete text, arguments and bibliography; [`Magenta_Gap.pdf`](Magenta_Gap.pdf) is its current 24-page build. The paper uses the standard LaTeX `article` class, with four figures and 33 references.

## Build

With a complete TeX installation providing pdfLaTeX, latexmk and the standard packages named in the preamble, run from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error Magenta_Gap.tex
```

For Overleaf, upload `Magenta_Gap.tex` and these four PDFs, set `Magenta_Gap.tex` as the main document, and select pdfLaTeX:

- `fig1_magenta_gap.pdf`
- `fig2_color_body.pdf`
- `fig3_spindle.pdf`
- `fisher_bound_convex_hull.pdf`

No companion note, BibTeX database, custom class or local compatibility shim is required. The three conceptual figures also have optional standalone TikZ sources beside their PDFs; rebuilding those requires the standard `standalone` class and the packages listed in each source. Keep `\noindent` at the start of the abstract.

## Reproduction

- [`../recompute/`](../recompute/): numerical calculation, pinned dependencies, input arrays and stored outputs. The manuscript cites the code and outputs at commit `a0282ae`.
- [`../research/object-colour/`](../research/object-colour/): exact-arithmetic continuum certificate and its declared input.
- [`../research/colour-interface/`](../research/colour-interface/): effect-algebra analysis and 11 executable checks.

The Fisher illustration is not psychophysically calibrated. A perceptual realization requires a specified stimulus encoding, observer conditions and readout.

A concise [version summary](../CHANGELOG.md) records substantive changes. Earlier drafts, import records and detailed editorial notes remain accessible in [snapshot `5176130`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/51761306772c3b234ba65099afe46afc051cd394).
