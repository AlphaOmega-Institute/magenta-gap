# Luminosity thresholds against bounded-reflectance limits — exploratory

**Status: exploratory and numerical. Nothing here is certified, and none of these numbers is used in a manuscript.**

Morimoto, Numata, Fukuda and Uchikawa (2021) measured the luminance at which a coloured patch, seen among other coloured
patches, stops looking like an illuminated surface and starts to look self-luminous. They compared those thresholds with
the largest luminance a reflecting surface can have at the same chromaticity. That upper limit is a boundary of the
bounded-reflectance body studied in `research/object-colour/`, reached by scaling outward along a chromaticity ray from
black. These scripts ask how closely the published thresholds follow that limit, in shape and in absolute level.

This is one small, testable piece of the question of how a physical preparation class (reflecting surfaces) relates to a
perceptual judgment. It does not bear on the manuscript's null-boundary hypothesis, on the certified obstruction to it, on
the cylinder and spindle picture, or on any calibration of a perceptual boundary.

## Inputs

None of the third-party files is stored here. Fetch them and pass their locations on the command line.

| Input | Source | Licence | Pinned by |
| --- | --- | --- | --- |
| Cone table | `research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv` in this repository | see that folder | SHA-256 `144fc58d…b9d02c40` |
| Raw thresholds | Zenodo, [doi:10.5281/zenodo.5590120](https://doi.org/10.5281/zenodo.5590120), `Data_v1.zip` | CC BY 4.0 | SHA-256 `7bcf39304749d98dd1861fa117b7f801f99ce0e25c628b4dba21b88215b54065` |
| Authors' analysis files | [github.com/takuma929/LuminosityThreshold](https://github.com/takuma929/LuminosityThreshold) | MIT | commit `f857dc7480e0a5d3c202d946f536537f8a8ec010` |

Files read from the authors' repository, with SHA-256:

```
5783ee8543edaa58a3f022718984b62c64de2b3e420125156216f24ce777f83d  LT_Data.mat
88ecaada8be1a6102b4d1e237e5d8b4929ece8e080b1a0ea99bffe82e456fd93  SpectrumLocus_MB.mat
bc662c026ad74f88c1af1e690ddc9228400f4120670ca603d6de2367dfe9b697  ill6500K.mat
c09334f76fe3a05a39ac9d50b5c2cf43b36a5470b2ab495fecb876015ea947cf  MB_bbl_500to25000with500step.mat
3729bbc5c8781644923fb6f97f58dcd32dad8a169d4f9d522e66166d66d49635  OP_MB/OP_6500_MB.mat
```

Only data files are read. The authors' MATLAB code is not run; its prediction rule is re-implemented from
`LT_saveCorrelationCoefficient.m`.

Reference: T. Morimoto, A. Numata, K. Fukuda and K. Uchikawa, "Luminosity thresholds of colored surfaces are determined
by their upper-limit luminances empirically internalized in the visual system," *Journal of Vision* 21(13):3 (2021),
[doi:10.1167/jov.21.13.3](https://doi.org/10.1167/jov.21.13.3).

## Running

```sh
python first_pass.py              TABLE_CSV ZENODO_DATA_DIR [S_AXIS_FACTOR]
python exp1_authors_comparison.py TABLE_CSV ZENODO_DATA_DIR AUTHORS_REPO_DIR
```

Needs NumPy and SciPy. The stored outputs are from Python 3.11, NumPy 2.3.5 and SciPy 1.17.0. Neither script is part of
`run_all.py`, because both need the external files above.

| Script | What it checks | Output |
| --- | --- | --- |
| `first_pass.py` | Experiments 1 and 2: correlation of each observer's settings with the limit computed from our table, and setting ÷ limit | `outputs/first_pass.txt`, and the two `…_s_axis_x….txt` files |
| `exp1_authors_comparison.py` | Experiment 1 at 6500 K: A conventions; B the authors' own prediction rule; C our optimisation against their stored optimal colours, point by point; D one factor per observer | `outputs/exp1_authors_comparison.txt` |

## What the outputs show

All of this is provisional: four observers per experiment, one laboratory, one judgment criterion.

**Conventions agree.** The authors' spectrum-locus file reproduces our cone table to within 0.4% at every wavelength from
400 to 700 nm once a global scale is removed, and their black-body locus agrees with ours to three parts in 10⁴ on the
S axis. Their stored white is reproduced to 10⁻⁸ by spline interpolation of their 10 nm illuminant file.

**The authors' rule and reported average are reproduced.** Their rule picks the nearest stored optimal colour in a
weighted chromaticity distance. Re-implemented, it gives a mean correlation of 0.9013 over four observers and five
surrounds; the paper reports 0.901.

**The two numerical procedures agree at the test chromaticities.** At the nine Experiment 1 chromaticities our constrained
optimisation and their nearest-sample lookup differ by −0.1% to +0.5% under their conventions, and by at most +0.6% under
ours (390–830 nm, Planck illuminant). At each selected sample's own chromaticity the optimisation returns the stored
luminance.

**Away from the test points the stored set is not everywhere the upper limit.** Of 374 stored optimal colours drawn at
random, 96.8% agree with the optimisation to 10⁻³; ten are beaten by more than 1%, all with S/(L+M) below 2 × 10⁻⁵. That
is the region where the S fundamental is tabulated as zero, the same flat face recorded in
`../../flat_face_counterexample.py`. It does not touch Experiment 1.

**Thresholds lie above the computed limit.** Median setting ÷ limit is 1.38 in Experiment 1 and 1.16–1.23 in Experiment 2,
with per-observer medians from about 0.95 to 1.7.

**One multiplicative factor per observer is not supported uniformly** (Experiment 1, authors' predictions, weighted by
the reported standard errors):

| Observer | Factor | χ²/dof, factor | χ²/dof, additive offset | χ²/dof, affine | Affine slope, intercept (cd/m²) | Setting ÷ limit across chromaticities |
| --- | --- | --- | --- | --- | --- | --- |
| KK | 1.62 | 11.5 | 20.8 | 11.5 | 1.73, −2.4 | 1.43–1.87, rising with L/(L+M) |
| MI* | 1.37 | 3.4 | 5.5 | 3.4 | 1.42, −1.2 | 1.30–1.47 |
| TM | 1.35 | 17.8 | 9.3 | 8.1 | 0.79, +13.4 | 1.14–1.46, lowest near white |
| YK* | 1.52 | 24.8 | 12.3 | 12.4 | 0.94, +13.9 | 1.28–1.61, lowest near white |

A factor fits better than an offset for two observers and worse for the other two. Systematic dependence on chromaticity
remains for three of the four. Every χ²/dof is far above one, so no one-parameter model describes the data within the
reported repetition error; the standard errors do not include differences between surrounds, which move the fitted factor
by 8–22% depending on the observer.

**The two public sources label two observers differently.** The Zenodo CSV and the repository's `LT_Data.mat` hold
identical numbers for Experiment 1, but the second and fourth observer columns carry the labels MI and YK in one and YK
and MI in the other. Observers above carry the Zenodo label, starred where the repository key differs.

## Corrections to earlier readings

- *"At 6500 K the ratio is roughly constant, so the thresholds look like the surface limit scaled outward from black."*
  Withdrawn. That came from averaging over observers and surrounds; the individual residuals have opposite trends.
- *"The offset reflects criterion and observer."* Too strong. Observers were asked for the midpoint between two perceptual
  endpoints, and neither endpoint has been identified with the physical limit. Criterion and observer differences may
  contribute; illuminant estimates, calibration and model mismatch remain possible contributors.
- *"I inferred the S-axis scaling."* Now confirmed against the authors' files. The ±10% runs deliberately miscalibrate the
  model against the data; they show that the 3000 K and 20,000 K correlations depend on that relative calibration, not
  that the physical limit depends on a coordinate choice.
- High correlation does not establish proportionality: adding a constant to every prediction leaves the correlation
  unchanged.

## What this does not show

It does not identify a measured perceptual boundary with the bounded-reflectance surface, establish an enlarged
appearance space, or support any claim about altered or chimerical colours. Seeing a source as luminous is ordinary
vision. Experiments 2 and 3 have not been compared against the authors' saved predictions.
