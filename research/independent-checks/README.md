# Independent colour checks — November 2026 audit

These are the second reviewer's scripts from the November 2026 audit of the Magenta Gap manuscripts. They were written
independently of the committed certificate and of `research/colour-interface/check_colour_interface.py`, and they reproduce
the numbers quoted during the review. They are preserved as evidence. They change no manuscript and select no geometry.

Nothing here is certified. The one rigorous result is the
[continuum rim certificate](../object-colour/Continuum_Rim_Certificate.md); `certificate_review.py` reproduces it in
floating point.

## Running

```sh
pip install -r research/independent-checks/requirements.txt
python research/independent-checks/run_all.py        # about two minutes; writes outputs/*.txt
python research/independent-checks/continuum_rim.py 830   # any single script can be run on its own
```

Tested with Python 3.11, NumPy 2.4 and SciPy 1.17. The stored `outputs/` are from that run. Lines reporting run times
differ between runs; `cubic_law_and_september_numbers.py` also prints a solver warning that does not affect its ratios.

Most scripts print quantities and do not assert. `run_all.py` reports only whether each script ran.

## Input

Every data-dependent script reads `research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv` through `common.py`:
the 441 Stockman–Sharpe 2° LMS samples, 390–830 nm, bit-identical to the arrays used during the audit. No other data file is
used. Scripts marked "none" below need no data.

## Status labels

- **Exact** — an identity or statement with a short proof; the script checks it numerically.
- **Numerical** — a computed estimate that depends on sampling, bin width or solver tolerance.
- **Reproduction** — recomputes a number stated elsewhere in this repository.

## Index

### Manuscript conventions and the bound

| Script | Input | What it checks | Status |
| --- | --- | --- | --- |
| `hull_and_purple_edge.py` | table, 390–700 nm | Hull of the 311 chromaticities in both conventions of the September Sec. 5; the 410–700 nm edge; `eq:edge-derivative` equals tan(θ − θ₀) | Reproduction (188 vertices, four non-adjacent edges); identity exact |
| `full_block_fisher.py` | table, 390–700 nm | Scale c\* of the manuscript's angular-only test; whether a full-block Fisher identification could make the boundary null | Numerical; reproduces c\* = 0.3761, 0.4028 |
| `chart_comparison.py` | table, 390–700 nm | max (u/β)² in the appendix chart and in the July opponent chart | Numerical; both reviewers agree to nine digits |
| `linear_channel_symmetry.py` | none | Symmetry of the image of independent bounded activations under the July opponent formulas with any linear white/black channel | Numerical illustration on chosen models |

### The symmetric body and its geometry

| Script | Input | What it checks | Status |
| --- | --- | --- | --- |
| `effect_body_and_vector_identity.py` | none | Bicone ⇔ 0 ≤ E ≤ 1 for real symmetric 2×2; determinant sheets; the vector identity X² + Y² + N² = Z² | Exact; prior art Berthier (2020) Eq. (80) |
| `hue_projectors_holonomy_boost.py` | none | Hue projector decomposition; sign holonomy; boost action on pure hues | Exact; boost map is Berthier (2020) Eq. (89) |
| `white_point_and_chord_region.py` | none | Boosts move white; fraction of the disc beyond a chord | Exact |
| `symmetric_body_and_order_interval.py` | table, 390–700 nm | `eq:baseline` = −4 det; central versus mirror symmetry; order-interval body of the measured hull; Hilbert angular scale at the boundary | Identities exact; ridge heights and volumes numerical |
| `span_collinearity_hilbert.py` | table, 390–700 nm | Spectral signals have rank 3; collinear S = 0 stretch; Hilbert metric on the disc and on the measured hull | Numerical; Hilbert construction is Berthier (2020) Sec. 7.3 |

### The bounded-reflectance signal model

| Script | Input | What it checks | Status |
| --- | --- | --- | --- |
| `solid_lib.py` | table | Helper: finite-bin body by linear programming | Exact for the finite-bin model, to solver tolerance |
| `bounded_reflectance_body.py` | table, 311 bins | Central symmetry, containment, volumes and widest-point heights of the 1 nm body | Numerical; see corrections in its header |
| `widest_point_types.py` | table, 0.25 nm bins | Reflectance type at the widest point; scaling of the slope jump with step size | Numerical evidence |
| `crease_slopes_and_near_black.py` | table, bins to 0.05 nm | Crease slopes predicted from two tangent planes; near-black behaviour of bins; deficit exponents | Numerical; exponents at 550 nm are solver noise |
| `continuum_rim.py` | table, PCHIP, no bins | Widest points of the continuum body on 390–700 nm, or 390–830 nm with argument `830`; supporting-normal criterion; crease certificates; hue-weighted shares | Numerical on 1,440 sampled normals |
| `three_transition_points.py` | table, PCHIP, 390–830 nm | Three-transition rim points satisfy the criterion; weighting by normals versus hue angle | Numerical |
| `crease_fan_and_550nm.py` | table, PCHIP and spline | Width of the supporting-normal fan at edge colours; the 550 nm contact at several scales | Numerical |
| `cubic_law_and_september_numbers.py` | table | Deficit / (A z³) with exact band integrals at 510 and 550 nm; the long and short contacts and τ\* of the September note | Reproduction (703.417797 nm, τ\* = 0.0674706704); single contacts only |
| `certificate_review.py` | certificate script and CSV | Core routines of the certificate; floating-point reproduction of its roots and support point; white as a positive combination of the transition vectors | Reproduction of a certified result |
| `flat_face_counterexample.py` | table, PCHIP, 390–830 nm | The S = 0 face through black and the catch vector (0.06566, 0.01201, 0) | Numerical verification of an exact counterexample |
| `timelike_direction_step.py` | none | A conormal annihilating a timelike direction is spacelike, on 200,000 random metrics | Numerical check of a short proof |

### The connection to LG

| Script | Input | What it checks | Status |
| --- | --- | --- | --- |
| `ground_plane_interval.py` | none | Squared Krein form of a ground-plane vector equals Z² − X² − Y²; the neutral circle is the null circle | Exact; pure vectors only |
| `lg_effect_family.py` | none | The bicone from LG's R₀, Q₀; its three symmetries; the ground-plane characterisation of the family | Exact |
| `two_encoding_routes.py` | none | Sector-phase statistics do not depend on the input; paired effects do; defect extensions agree on paired states | Exact, on one random finite model |
| `defect_extensions_and_altered_points.py` | none | Range-mixing dynamics separates defect extensions; illustrative altered-colour points against a bicone | Exact; the points are illustrative, not measured |

## Exploratory work kept separately

[`exploratory/luminosity-thresholds/`](exploratory/luminosity-thresholds/README.md) compares published luminosity-threshold
measurements (Morimoto et al. 2021) with bounded-reflectance limits computed from the same cone table. It is exploratory,
needs third-party files that are not stored here, is not run by `run_all.py`, and supports no manuscript claim.

## Statements corrected during the review

Each affected script carries the correction in its header. In short:

- A finite-bin body is a zonotope, so it cannot show a smooth maximum; exact contact with the cone near black is a bin
  artifact.
- The exponents near one at 550 nm were linear-programming error. 550 nm is an ordinary quadratic contact.
- The claim that no boundary point other than black has a tangent plane containing its direction from black is false for
  this table (the S = 0 face).
- Central symmetry of an attainable body needs hypotheses the July manuscript does not supply; the sampled percentages
  describe chosen models.
- The region beyond a mixture chord is a convex restriction, not a linear defect.
- The Hilbert metric on a general gamut is prior art (Berthier 2020, Sec. 7.3).
- Shares of the hue circle are estimates, weighted by hue angle of the boundary point.
- The altered-colour incompatibility is conditional on a calibration and does not show that the Berthier–Provenzi approach
  excludes altered colours.

## Scripts not preserved

- A purple-chord script using the historical 390–700 nm chord with raw channels centred at 1/3, a convention the
  September manuscript rejects. Superseded by `hull_and_purple_edge.py`.
- A conic fit to the spectral locus, whose residual was withdrawn as coordinate-dependent. The coordinate-free statement
  is the collinear stretch in `span_collinearity_hilbert.py`.
- An earlier corner classifier whose "smooth maximum" label was wrong. Superseded by `widest_point_types.py` and
  `continuum_rim.py`.
- An earlier version of the July symmetry check that used a white/black formula recalled from memory. Superseded by
  `linear_channel_symmetry.py`, which assumes no formula.
- Two text-extraction utilities used to read papers.

## Sources

- M. Berthier (2020), *Geometry of color perception. Part 2: perceived colors from real quantum states and Hering's rebit*,
  J. Math. Neurosci. 10:14, [doi:10.1186/s13408-020-00092-x](https://doi.org/10.1186/s13408-020-00092-x). Sections 5 to 7
  were read for this record.
- M. Berthier and E. Provenzi (2022), *Quantum measurement and color perception: theory and applications*,
  [arXiv:2504.12304](https://arxiv.org/abs/2504.12304).
- The September object-colour notes and the certificate in `research/object-colour/`.
- LG, `manuscripts/lg/LG_Master_Unified.tex` in AlphaOmega-Theory, for the ground-plane labels cited in the headers.
