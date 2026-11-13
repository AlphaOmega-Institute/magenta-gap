# Numerical gamut-bound calculation

This folder computes the paper's Fisher-metric compatibility illustration for a specified chromaticity gamut. It contains the calculation script, pinned dependencies, input arrays and the reference summary [`bound_results.json`](bound_results.json). It fits no empirical discrimination data and does not calculate the separate cubic-contact asymptotic.

## Run

From the repository root, in a Python environment compatible with the pinned packages:

```sh
python -m pip install -r recompute/requirements.txt
python recompute/recompute_bound.py
```

The script reads `inputs/LMS.npy` and `inputs/wavelengths.npy`. It writes the result summary, ten CSVs, four diagnostic NPZs, and three figures in PDF, SVG and PNG formats into this folder. The generated tables and figure variants are ignored by Git; the reference JSON is tracked. A rerun may change its software-version metadata.

The paper uses [`manuscript/fisher_bound_convex_hull.pdf`](../manuscript/fisher_bound_convex_hull.pdf). The other two figures, `fisher_bound_comparison` and `fisher_bound_sensitivity`, are spectral-arc diagnostics, not the primary gamut result. Rendered image bytes can vary by platform even when the plotted data agree.

## Inputs and conventions

The inputs are the Stockman–Sharpe 2-degree fundamentals supplied by Colour Science 0.4.7, sampled at 1 nm from 390 to 830 nm. Their SHA-256 hashes and software versions are recorded in the JSON; the [input license](inputs/COLOUR_SCIENCE_LICENSE.txt) is retained.

The raw equal-energy channel integrals over that interval are

```
L = 115.97840823; M = 94.82117898; S = 58.41945060.
```

Raw equal-energy white therefore has chromaticity `(0.4307957163, 0.3522082976, 0.2169959861)`. The E-balanced convention first divides each channel by its own integral, giving white at `(1/3,1/3,1/3)`. The raw-channel comparison uses the actual integrated white as its polar origin.

For each convention, the primary boundary is the convex hull of the 311 normalized responses at 390–700 nm. This is the gamut of nonnegative mixtures of those tabulated stimuli; piecewise-linear interpolation of the LMS responses gives the same hull. It is not a measured surface-color solid or a model of the entire continuous human gamut.

Both hulls have 188 vertices. Nonadjacent wavelength edges join 410–700, 453–458, 599–616 and 616–700 nm. The last edge lies on the tabulated S=0 boundary; the preceding edge meets it at 616 nm. The 390–409 nm samples lie inside the hull, so the spectral arc closed by its endpoint chord is a different boundary.

## Method and reference results

On each open hull edge, the radial slope is evaluated from analytic ray/segment formulas in floating point. For `x=a+td`, centered on white,

\[
u=|\alpha'/\alpha|=\frac{|x\cdot d|}{|\det(x,d)|},\qquad
\beta^2=\frac{\sin^2\theta}{l}+\frac{\cos^2\theta}{m}
+\frac{(\sin\theta-\cos\theta)^2}{s}.
\]

There are 101 interior samples per edge, 18,988 in total. The fixed reporting subset is finite `beta(c=1)<10`. Corners have no unique angular derivative, and sampled extrema are not certified continuous extrema.

| Quantity | E-balanced | Raw, E-centered |
| --- | ---: | ---: |
| Retained points | 14,436 | 12,267 |
| Points satisfying `u < beta` at c=1 | 14,436 | 12,267 |
| Minimum `beta-u` | 1.008471 | 0.882087 |
| Sampled critical scale `max(u/beta)^2` | 0.3760819070 | 0.4028317754 |

For `g_disc -> c g_disc`, beta becomes `sqrt(c) beta`. On the fixed sample set, strict compatibility holds exactly when `c` exceeds the sampled critical scale. The choice `c=1` and the cutoff `beta<10` are uncalibrated model conventions. S is tabulated as zero from 616 nm onward; the script represents the Fisher divergence as infinity without an epsilon floor.

Assigning the simplex metric separately after channel balancing changes the metric model. The two columns are not coordinate transforms of one fixed physical tensor.

## Spectral-arc control

An additional control closes the entire spectral arc by its 390–700 nm endpoint chord. A coarse 1-degree angular grid misses a narrow raw-channel feature. Analytic derivatives of a PCHIP interpolation of the LMS responses resolve it: at 403.4 nm, `u=5.9852316075`, `beta=5.2311889036`, and `(u/beta)^2=1.3090647100`. These points are inside the mixture hull, so this is not a failure of the primary hull calculation.

## Scope

The calculation establishes compatibility on the declared samples under the stated metric convention. It does not calibrate human discrimination, certify an all-hue inequality, or identify a perceptual boundary. The exact-arithmetic continuum certificate uses a different bounded-reflectance model and is documented in [`research/object-colour/`](../research/object-colour/).
