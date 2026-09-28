# Corrected original gamut-bound illustration

This calculation repairs the numerical illustration in Appendix A of the
original Magenta Gap manuscript. It does **not** analyze the later cubic
object-color-boundary theorem. No empirical discrimination data are fitted.

## Primary model and scope

Let `LMS.npy` contain the supplied Stockman–Sharpe 2-degree fundamentals from
Colour Science 0.4.7, on the saved 390–830 nm, 1 nm grid. The raw E channel
integrals, using the trapezoidal rule over that full grid, are

```
L = 115.97840823; M = 94.82117898; S = 58.41945060.
```

Thus E has raw chromaticity `(0.4307957163, 0.3522082976, 0.2169959861)`,
not `(1/3,1/3,1/3)`. For the primary E-balanced convention, each sensitivity
channel is divided by its own integral before chromaticities are formed.
For the raw-channel control, channels remain unchanged and the polar origin
is their actual integrated E chromaticity.

For each convention the primary boundary is the convex hull of the 311
normalized tabulated responses at 390–700 nm. This is the chromaticity gamut
of arbitrary nonnegative mixtures of those tabulated stimuli. It is also the
gamut for piecewise-linear interpolation of the *LMS responses* between those
wavelengths: normalization maps each interpolated response to the corresponding
chromaticity segment. It is a declared discretized stimulus model, not a claim
about a measured surface-color solid or a full continuous human gamut.

The hull has 188 vertices in both conventions. Besides edges joining adjacent
wavelength samples, it contains edges joining 410–700, 453–458, 599–616, and
616–700 nm. The last two lie on the tabulated S=0 boundary. In particular,
390–409 nm spectral samples lie inside this hull. The original assertion that
the entire spectral locus itself was the gamut boundary was incorrect for
these data.

The radial function is computed by exact ray/segment intersection. At a point
`x = a + t d` on an open hull edge, with origin at E,

\[
u=|\alpha'/\alpha|=\frac{|x\cdot d|}{|\det(x,d)|},\qquad
\beta^2=\frac{\sin^2\theta}{l}+\frac{\cos^2\theta}{m}
+\frac{(\sin\theta-\cos\theta)^2}{s}.
\]

No numerical differentiation is required for this open-edge calculation.
The polygon has corners: no unique derivative is asserted there. Interior
edge samples approach each endpoint from its own side, but do not constitute
a theorem about corner smoothing or continuous extrema.

## Results to use in the manuscript

There are 101 points per open edge (`t=0.00001,…,0.99999`), 18,988 total.
The *fixed* reporting subset is `beta(c=1)<10`. Table entries are sampled
extrema, not certified continuous extrema.

| Quantity | E-balanced channels | Raw channels, E center |
|---|---:|---:|
| Retained points | 14,436 | 12,267 |
| Points satisfying `u < beta` at c=1 | 14,436 | 12,267 |
| Minimum `beta-u` | 1.008471 | 0.882087 |
| Angle of minimum margin | −152.6914° | −156.7225° |
| Maximum squared ratio `(u/beta)^2` | 0.3760819070 | 0.4028317754 |
| Angle of maximum squared ratio | −147.1657° | −153.4432° |
| Maximum retained u | 2.179950 | 1.925549 |
| Retained beta range | 1.414425–9.999271 | 1.418761–9.999092 |

The 360-point angular grid also passes at c=1: 282/282 retained E-balanced
samples and 255/255 raw-E samples. These counts describe the full hull, not
the spectral subset alone. Primary results above use the more resolved
open-edge checks, avoiding the original finite-difference issue.

For `g_disc -> c g_disc`, beta becomes `sqrt(c) beta`. On the fixed retained
sample set the strict inequality holds exactly when

\[
c>c_*:=\max(u/\beta)^2.
\]

Thus the sampled critical scales are 0.37608 and 0.40283. The choice c=1 is
an **uncalibrated simplex-model convention**. The result shows compatibility
under that convention; it does not establish a physiological inequality.
The subset is held fixed when discussing c: otherwise rescaling would also
change which points are reported. The cutoff beta<10 is a convenience for
excluding singular tails; it is not a measured discrimination threshold.
S is tabulated as exactly zero from 616 nm onward, and the script treats the
corresponding Fisher divergence as infinity, without an arbitrary epsilon floor.

E balancing changes the chromaticity model. Reassigning the canonical simplex
metric after balancing is **not** merely transforming the same physical metric
under a coordinate change. The two columns are distinct modeling conventions.
The fixed Euclidean benchmark beta=1 also has a coherent tensor transformation;
its numerical failure must not be called a failure of coordinate invariance.

## Original-locus countercontrol

The script also reproduces the original choice of the entire spectral arc
plus the 390–700 nm endpoint chord, with corrected white points and direct
metric evaluation. That radial curve is not the mixture-gamut hull.
Coarse 1° central differences pass all retained spectral samples (170/170
E-balanced and 145/145 raw-E), but hide a narrow violet feature.

A 0.01 nm wavelength grid with PCHIP interpolation of the original LMS
responses and analytic PCHIP derivatives finds a raw-E failure between
401.60 and 404.78 nm. At 403.4 nm,

```
u = 5.9852316075; beta = 5.2311889036; (u/beta)^2 = 1.3090647100.
```

The corresponding E-balanced maximum squared ratio is 0.7217116474 and
remains below one. The polygonal-locus control independently finds the raw
failure. These points are inside the convex hull, so this does not contradict
the primary hull result. It does show why the original table and a claim of
normalization-independent, all-hue verification should not be retained.

## Figures and reproducibility

Run `python3 recompute_bound.py` from any directory. It reads the bundled
unaltered inputs and writes `bound_results.json`, grid CSVs, diagnostic NPZs,
and figures. Exact package versions and input SHA-256 hashes are in the JSON.

Use `fisher_bound_convex_hull.pdf` or `.svg` as the main vector figure. It shows
the primary E-balanced angular coefficients and the metric-scale ratios for
both conventions. Pink shading marks hull edges joining nonadjacent sampled
wavelengths; it is not a hue-perception classification. Values beta>=10 are
omitted from the beta trace and the scale-ratio panel. Curves are displayed
at 1° spacing; table maxima use the denser open-edge samples. The dotted line
in the upper panel is the fixed Euclidean beta=1 benchmark.

`fisher_bound_sensitivity.pdf` documents the original-locus narrow violet
failure and is a diagnostic, not the primary gamut result.
`fisher_bound_comparison.pdf` is an additional original-locus diagnostic.

The corrected script and plots were prepared with ChatGPT/Codex code
assistance. Figures are deterministic Matplotlib plots of computations from
the supplied data; no generative raster imagery or synthesized measurements
were used. Model/version metadata beyond the tool name should not be guessed.
The author should inspect the code, assumptions, and output before submission.

Repository source inspected read-only on 2026-09-28:

- <https://github.com/AlphaOmega-Institute/magenta-gap/blob/main/step2_spectral_locus.py>
  (blob `7f583c09efb72cfbb80abd3e03d336357c77832a`): incorrect unbalanced E center.
- <https://github.com/AlphaOmega-Institute/magenta-gap/blob/main/step3_close_and_sample.py>
  (blob `8cc513fb1abaa81efd1c707d8354691e8379604f`): 60-point endpoint chord and polar interpolation.
- <https://github.com/AlphaOmega-Institute/magenta-gap/blob/main/step4_alpha_prime.py>
  (blob `b94b7a4b1e18e2741d346327dc531b93c5f865cb`): 1° central differences.
- <https://github.com/AlphaOmega-Institute/magenta-gap/blob/main/step5_fisher_beta.py>
  (blob `a2cfd58ee4f30fe873ffb45f70d5288871de47ce`): metric interpolation and epsilon floor.
- <https://github.com/AlphaOmega-Institute/magenta-gap/blob/main/fix1_purple_beta.py>
  (blob `f8b0c43bf026a33c4a601760f76770248a959676`): direct purple metric correction.

No changes were made to that repository. The old MacAdam count is not carried
forward as independent validation: no shared psychophysical metric scale has
been established.
