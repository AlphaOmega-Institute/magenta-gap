# magenta-gap

Numerical verification of the Lorentzian gamut bound, companion code for:

> Yiannopoulos, A. (2026). *The Magenta Gap: A Geometric Reachability Model of Nonspectral Hue and Color Gamut Boundaries.*

## What this code does

The paper derives a necessary constraint on any admissible psychophysical model of color perception. If the hue-dependent saturation ceiling α(θ) of the visual gamut is modeled as the boundary of a piecewise null reachability envelope in a (2,1) Lorentzian metric, then the following inequality must hold at every hue angle θ:

    |α'(θ) / α(θ)| < β(θ)

where β(θ) is determined by the local Fisher information metric of the cone fundamentals. This inequality is coordinate-invariant: under any smooth reparameterization θ → f(θ), both sides scale by 1/|f'|, so the bound is preserved. It is strictly stronger than the naive flat-space version |α'/α| < 1, which depended on an arbitrary 2π normalization of the hue coordinate.

This repository verifies that bound numerically against:
- The Stockman-Sharpe 2-degree cone fundamentals (via colour-science 0.4.7)
- The MacAdam (1942) empirical discrimination ellipses

## Results

| Quantity | Value |
|---|---|
| α(θ) range | [0.1082, 0.6577] |
| Flat-space violations \|α'/α\| > 1 | 90/360 (25.0%) |
| Non-degenerate angles (β < 10) | 297/360 |
| Bound satisfied (non-degenerate) | 297/297 (100%) |
| Bound satisfied (spectral arc) | 188/188 (100%) |
| Binding constraint — tightest margin | 0.542 at θ = -149° (\|α'/α\| = 1.933, β = 2.475) |
| max \|α'/α\| (degenerate wing) | 2.480 at θ = -23° (β → ∞; bound holds trivially) |
| MacAdam ellipse check | 25/25 (100%) |

The flat-space bound fails 25% of the time — as predicted by the theory, since it is not coordinate-invariant. The Fisher-corrected bound holds at every non-degenerate angle. The test is non-trivial: the 65 angles where the flat bound fails but β < 10 are rescued by genuine Fisher curvature (β ranging from 2.13 to 7.45), not by degenerate blow-up. The binding constraint is the *smallest margin* (0.542 at θ = -149°, the violet shoulder), not the largest derivative: |α'/α| peaks at 2.480 near the red terminus (θ = -23°), but β diverges there into the degenerate wing, so the bound holds there with astronomical room rather than binding.

## Environment setup

    conda create -n color python=3.12 numpy scipy matplotlib -y
    conda activate color
    pip install colour-science
    python3 -c "import colour, numpy, scipy; print('OK')"

All scripts must be run from this directory with the color env active.

## Pipeline

Run scripts in order. Each step saves .npy files consumed by the next.

| Script | Purpose |
|---|---|
| step1_cmfs.py | Load Stockman-Sharpe cone fundamentals |
| step2_spectral_locus.py | Compute spectral locus and raw α(θ) |
| step3_close_and_sample.py | Close gamut with purple line, resample uniformly |
| step4_alpha_prime.py | Differentiate α, compute |α'/α| |
| step5_fisher_beta.py | Compute Fisher β, run the bound test |
| step5b_diagnostic.py | Non-degenerate region analysis |
| step5c_plot_clipped.py | Publication-quality plot |
| fix1_purple_beta.py | Recompute β in purple-chord region from mixture chromaticities |
| fix2_macadam.py | MacAdam ellipse cross-check |

Run fix1_purple_beta.py after step5_fisher_beta.py to get the corrected 297/297 non-degenerate result. Run fix2_macadam.py after fix1_purple_beta.py for the MacAdam cross-check. step5c_plot_clipped.py automatically uses the purple-corrected β when fix1 has already been run; re-run it after fix1 to regenerate the publication plot over the full 297-angle non-degenerate region (run before fix1, it plots the 188-angle spectral arc).

## Key technical notes

**Purple region β correction (fix1):** The original pipeline interpolated β in the purple-chord region from adjacent spectral locus values, producing physically meaningless β values (6751-1,315,047) from the simplex singularity. fix1_purple_beta.py recomputes β at purple-chord grid angles from the actual chromaticity of the purple mixture at each angle. This expands the non-degenerate region from 188 to 297 angles.

**MacAdam scale caveat (fix2):** The MacAdam β values (284-3910) are much larger than the Fisher simplex β values (1.4-9.4). Both are consistent with the bound, but direct quantitative comparison requires a luminance-normalized photon-noise Fisher metric. This is left for future work.

**Coordinate invariance:** The bound |α'/α| < β is coordinate-invariant. Under θ → f(θ), α'/α scales by 1/|f'| and β = sqrt(g_θθ)/r scales by 1/|f'|, so the inequality is preserved.

## Dependencies

- Python 3.12
- colour-science 0.4.7
- numpy >= 2.0
- scipy
- matplotlib

## License

MIT

## Citation

    @article{Yiannopoulos2026,
      author  = {Yiannopoulos, Alexander},
      title   = {The Magenta Gap: A Geometric Reachability Model of Nonspectral Hue and Color Gamut Boundaries},
      journal = {Preprint},
      year    = {2026}
    }
