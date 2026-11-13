# Magenta Gap

A geometric reachability model of nonspectral hue and color-gamut boundaries. This repository is a public work in progress, with one active manuscript.

## Working manuscript

- [Magenta_Gap.tex](manuscript/Magenta_Gap.tex) — complete text, arguments and bibliography.
- [Magenta_Gap.pdf](manuscript/Magenta_Gap.pdf) — the reviewed 24-page build.
- [Build instructions](manuscript/README.md) — standard LaTeX article, four figure PDFs, Overleaf instructions.
- [Version summary](CHANGELOG.md) — substantive changes between versions.

The working paper has four figures and 33 references. Its assumptions, empirical scope and AI-assistance disclosure are stated in the manuscript.

## Calculations and research

| Location | Contents and status |
| --- | --- |
| [`recompute/`](recompute/) | Numerical calculation, pinned dependencies, input arrays and a saved result summary used by the paper. |
| [`research/object-colour/`](research/object-colour/) | Exact-arithmetic continuum certificate, declared input table, and cubic-contact proofs and endpoint analysis. |
| [`research/colour-interface/`](research/colour-interface/) | Effect-algebra analysis and 11 executable checks. |
| [Independent checks at `a3a947c`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/a3a947c0c27d76dcb4b5e57d7b1329791fc771e7/research/independent-checks) | Claude's reproductions and exploratory luminosity comparison; preserved on a separate branch, not merged into `main`. |

Run the current algebra checks and certificate from the repository root:

```sh
python research/colour-interface/check_colour_interface.py
python research/object-colour/certify_continuum_rim.py
```

For the numerical calculation, follow the dependency and reproduction instructions in [`recompute/README.md`](recompute/README.md).

The paper distinguishes a conditional null-boundary compatibility theorem, an exact effect-algebra connection, and a certified obstruction for one specified bounded-reflectance model. Its Fisher illustration remains uncalibrated. A perceptual realization requires an independently specified encoding and readout.

## Earlier versions

Earlier drafts, legacy scripts, import records and the development log remain available in [snapshot `5176130`](https://github.com/AlphaOmega-Institute/magenta-gap/tree/51761306772c3b234ba65099afe46afc051cd394).

## License

[MIT](LICENSE).
