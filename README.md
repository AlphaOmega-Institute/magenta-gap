# Magenta Gap

A geometric reachability model of nonspectral hue and color-gamut boundaries. This repository is a public work in progress, with one active manuscript.

## Working manuscript

- [Magenta_Gap.tex](manuscript/Magenta_Gap.tex) — complete text, arguments and bibliography.
- [Magenta_Gap.pdf](manuscript/Magenta_Gap.pdf) — the reviewed 24-page build.
- [Build instructions](manuscript/README.md) — standard LaTeX article, four figure PDFs, Overleaf instructions.
- [Editorial history](manuscript/EDITORIAL_NOTES.md) — corrections, reviews and provenance.

The conceptual restoration and figure repairs were reviewed and merged at `44ee356`. The working paper has four figures and 33 references. It retains its AI-assistance disclosure. No submission venue or completed perceptual encoding is selected.

## Calculations and research

| Location | Contents and status |
| --- | --- |
| [`recompute/`](recompute/) | Corrected calculation, pinned dependencies, input arrays and stored numerical outputs used by the paper. |
| [`research/object-colour/`](research/object-colour/) | Exact-arithmetic continuum certificate, declared input table, and historical cubic-contact proof and audit. |
| [`research/colour-interface/`](research/colour-interface/) | Effect-algebra analysis and 11 executable checks. |
| [Independent checks at `a3a947c`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/a3a947c0c27d76dcb4b5e57d7b1329791fc771e7/research/independent-checks) | Claude's reproductions and exploratory luminosity comparison; preserved on a separate branch, not merged into `main`. |

Run the current algebra checks and certificate from the repository root:

```sh
python research/colour-interface/check_colour_interface.py
python research/object-colour/certify_continuum_rim.py
```

For the corrected numerical calculation, follow the dependency and reproduction instructions in [`recompute/README.md`](recompute/README.md).

The paper distinguishes a conditional null-boundary compatibility theorem, an exact effect-algebra connection, and a certified obstruction for one specified bounded-reflectance model. Its Fisher illustration remains uncalibrated. Neither the numerical illustration nor the historical MacAdam comparison establishes psychophysical validation.

## Archive

The [archive index](archive/README.md) separates historical material from the working paper:

- [July 2026 manuscript](archive/july-2026/) — original 51-page PDF and recovered TeX, frozen with superseded claims intact.
- [Legacy pipeline](archive/legacy-pipeline/) — the nine superseded root scripts and their historical documentation.
- [Earlier versions and retired branch tips](archive/README.md#versions-and-branch-history) — immutable commit links and recovery information.

Historical research notes retain their dated status statements. The current manuscript and the index above identify what has since been incorporated. Git history is preserved.

## License

[MIT](LICENSE).
