# A fixed linear map from spectra to a colour readout — exploratory

**Status: exploratory and numerical. Nothing here is certified, fitted to appearance data, or used in a manuscript.**

The question is whether one fixed, injective linear map from stimulus states to observer states can reproduce the
cone-response calculation, and what a two-level (opponent) colour factor then requires. Both scripts read only
`research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv`.

```sh
python rebit_embedding.py          research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv
python cone_register_and_trine.py  research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv
```

Needs NumPy and SciPy. Stored outputs are from Python 3.11, NumPy 2.3.5 and SciPy 1.17.0.

| Script | What it checks | Output |
| --- | --- | --- |
| `cone_register_and_trine.py` | One explicit map with three positive cone channels and a two-level opponent factor: injectivity, wavelength-diagonal pull-backs, the exact obstruction to a positive two-level S channel, and where lights sit in the state disc | `outputs/cone_register_and_trine.txt` |
| `rebit_embedding.py` | Minimum-area ellipses enclosing the spectral chromaticities in two chromaticity planes, with free centre and with centre at white | `outputs/rebit_embedding.txt` |

## What the outputs show

**A linear, injective map reproducing the cone catches exists.** Send each wavelength to an output that keeps its
wavelength label and carries its three cone amplitudes on a three-level register. The three register projections pull back
to the cone sensitivities, and the map is injective because every total response in the table is positive. Different
spectra stay different complete output states even when the three measurements cannot tell them apart. (This form of the
construction is the first reviewer's.)

**The cone channels cannot be positive operators on a two-level factor alone.** A nonzero positive operator on two levels
has a kernel of dimension at most one. The table has S = 0 from 616 to 830 nm with L/(L+M) running from 0.80654 at 616 nm
to 0.94167 at 700 nm, so those lights would have to share one colour state and one L:M ratio. This depends on the exact
tabulated zeros. (Also the first reviewer's observation.)

**Both can be had in one map.** Attach to each register level one of three real two-level states at 60° spacing. The
reduced two-level state of a light then has Bloch vector

    x = a_L − (a_M + a_S)/2,   y = (√3/2)(a_M − a_S),   divided by a_L + a_M + a_S,

the standard opponent combination of white-balanced catches. Equal-energy white is the maximally mixed state. Every
selected observable still pulls back to a wavelength-diagonal operator, so for these observables and incoherent light the
model is classical; the two-level factor adds states and symmetry, not new statistics.

**Three nested convex sets result.** In the state disc the three register levels are three rim points; the triangle they
span holds every cone-excitation pattern and covers 41.3% of the disc; mixtures of 390–700 nm lights cover 55.8% of that
triangle, 23.1% of the disc. No monochromatic light is a pure state: the largest purity radius is 0.941, at 416 nm. The
least pure light is at 494 nm (radius 0.176), with a local minimum of 0.500 at 569 nm.

**Enclosing ellipses do not settle the embedding.** Across two chromaticity planes and two centring rules the least pure
light falls at 490, 570, 585 or 600 nm and lights cover 33–54% of the disc. Reproducing a measured saturation minimum by
choosing among these is a constraint that could select an embedding, not a prediction.

## Scope and corrections

- The construction with three rim states is fixed by the white balance. A different channel scaling gives a different disc.
- Mixing is preserved by the map on unnormalised states. Normalising the output changes mixture weights.
- Coherent superposition and incoherent mixture agree for the selected measurements only because the pulled-back
  observables are wavelength-diagonal here. An arbitrary linear map does not have that property.
- A two-level colour factor must sit inside a larger output space, and mixed two-level states come from a vector map only
  through that larger space. Ellipse containment makes the cone channels signed observables on the two-level factor, not
  positive ones.
- Treating dimming as mixing with a dark state is a preparation rule that needs its own justification; it does not follow
  from a ground axiom. With the dark state at the disc centre the radius grows with the light's weight and hue is fixed;
  with it off-centre hue can shift but the radius need not grow monotonically.
- No calibrated relation between purity radius and perceived saturation is established. Hue as disc angle and saturation
  from the state's entropy are Berthier's proposals (2020); he gives no map from spectra.
- Nothing here supplies independently selected grounds, dynamics, context dependence, or a perceptual boundary.
