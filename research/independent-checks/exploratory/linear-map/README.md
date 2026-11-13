# A fixed linear map from spectra to a colour readout — exploratory

**Status: exploratory and numerical. Nothing here is certified, fitted to appearance data, or used in a manuscript.**

The question is whether one fixed, injective linear map from stimulus states to observer states can reproduce the
cone-response calculation, what a two-level (opponent) colour factor then requires, and how an adapting condition could
act on it. The first two scripts read only `research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv`. The third
reads the corresponding-colour data of Luo and Rhodes (1999), [doi:10.5281/zenodo.3270903](https://doi.org/10.5281/zenodo.3270903):
37 `.dat` files, 39,896 bytes. That record states no licence, so the files are not stored here.

```sh
python rebit_embedding.py          research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv
python cone_register_and_trine.py  research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv
python adaptation_rules.py             DIRECTORY_WITH_THE_DAT_FILES
python adaptation_strength_heldout.py  DIRECTORY_WITH_THE_DAT_FILES
python map_family_structure.py         research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv
```

Needs NumPy and SciPy. Stored outputs are from Python 3.11, NumPy 2.3.5 and SciPy 1.17.0.

| Script | What it checks | Output |
| --- | --- | --- |
| `cone_register_and_trine.py` | One explicit map with three positive cone channels and a two-level opponent factor: injectivity, wavelength-diagonal pull-backs, the exact obstruction to a positive two-level S channel, and where lights sit in the state disc | `outputs/cone_register_and_trine.txt` |
| `rebit_embedding.py` | Approximate minimum-area ellipses enclosing the spectral chromaticities in two chromaticity planes, with free centre and with centre at white | `outputs/rebit_embedding.txt` |
| `adaptation_rules.py` | Two ways an adapting condition could act on the colour state, their algebra, and their predictions for 654 corresponding-colour pairs | `outputs/adaptation_rules.txt` |
| `adaptation_strength_heldout.py` | The same rules with one strength parameter each, fitted on some experiments and scored on others | `outputs/adaptation_strength_heldout.txt` |
| `map_family_structure.py` | Linear algebra of the condition-dependent family of maps: injectivity, ranges, consistency of the matching rule, and the constraint on a ground correspondence | `outputs/map_family_structure.txt` |

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

## Reweighting the register versus filtering the two-level factor

An adapting condition can act on the colour state in two ways that give the same three cone outputs for every spectrum
(the distinction and its first example are the first reviewer's):

- **Rule R, reweight the cone register.** Multiply the three catches by gains; the three rim states stay fixed. With
  gains equal to the reciprocal catches of the adapting white this is the classical von Kries rule.
- **Rule F, filter the two-level factor.** Apply a positive operator to the two-level state and renormalise. The weights
  change in the same way, and the rim states attached to the cone channels move as well.

For given gains the positive filter is unique; other filters with the same gains differ from it by an orthogonal factor.
It exists exactly when the reweighted white has radius at most 1/2, and is invertible, as the matching predictions need,
only strictly below 1/2. It then sends white exactly twice as far as reweighting does. So knowing the cone gains does not
determine the colour state.

The filter that sends an adapting white with cone proportions p_i to the disc centre has cone gains proportional to
1 − p_i, whereas reweighting uses 1/p_i. The two agree only for balanced white. So a filter cannot both reproduce the
reciprocal cone gains and make the adapting white neutral; the two filtering rules below are different proposals.

A corresponding-colour pair is two stimuli that look alike under two adapting conditions. Each rule predicts the second
chromaticity from the first by equating colour states, which needs no calibration of hue or saturation. Mean error in
CIE 1976 u′v′ over the 25 data sets with an illuminant change (654 pairs), assuming complete adaptation:

| Rule | Mean | Median | Notes |
| --- | --- | --- | --- |
| No adaptation | 0.0633 | 0.0625 | |
| R, von Kries gains | 0.0158 | 0.0117 | closer than F-white on 67% of pairs and in 20 of 25 data sets |
| F, filter sending the adapting white to the centre | 0.0222 | 0.0179 | 101 predicted matches have negative CIE Z, which no nonnegative light spectrum can produce |
| F, filter reproducing the von Kries gains | 0.0806 | 0.0795 | 614 pairs; no such filter exists for one data set; no adaptation gives 0.0619 on the same pairs |
| CAT02 (benchmark) | 0.0128 | 0.0092 | |

**Conclusion as agreed in review.** Under the specified encoding and matching rule, cone-register reweighting describes
these corresponding-colour data better than the two tested filtering rules. The white-centring filter also predicts
physically unavailable matches for numerous samples. These results constrain the proposed realization without determining
the full adaptation mechanism. The poor result for the gain-matching filter rejects the extra identification that
reproducing the cone gains also implements the intended white adaptation, not filtering in general.

### Strength fitted on some experiments, scored on others

Each rule gets one strength parameter between no adaptation and the full transformation, and the encoding is held fixed.
The viewing conditions recorded with the data define the held-out designs.

| Held-out design | Pairs | Reweighting, linear gains | Reweighting, geometric gains | Centring filter |
| --- | --- | --- | --- | --- |
| 1000 lux, 20% background: leave one of five experiments out | 413 | 0.0103 | 0.0099 | 0.0186 |
| Breneman: leave one of nine phases out | 107 | 0.0107 | 0.0115 | 0.0217 |
| Braun and Fairchild: leave one of four phases out | 66 | 0.0096 | 0.0093 | 0.0108 |
| McCann: leave one of four phases out | 68 | 0.0225 | 0.0267 | 0.0305 |
| All held-out pairs | 654 | 0.0116 | 0.0119 | 0.0196 |

Linear-gain reweighting has the lower error in 20 of the 22 held-out units. At the fitted strengths the centring filter
still predicts a negative CIE Z for 78 pairs. In-sample, the best single strength is 0.89 (linear gains), 0.85 (geometric
gains) and 0.85 (filter), with errors 0.0129, 0.0121 and 0.0201.

The fitted strength for linear-gain reweighting is 0.92–0.94 for the 1000 lux experiments, about 0.8 for Breneman and for
Braun and Fairchild (129 lux), and 0.56–0.68 for McCann (14–40 lux, strongly coloured illuminants). That ordering follows
the recorded illuminance, but illuminance is confounded here with medium, method and illuminant.

Two qualifications. The viewing conditions decide which experiments are grouped for fitting; they do not determine the
strength through an established illuminance law. And the five-experiment comparison at 1000 lux is the strongest transfer
test; the other three hold out phases within one study.

The 25 phases are not 25 independent studies; several share observers, samples or apparatus. Corresponding-colour matches
constrain which stimuli are equivalent across conditions. They do not identify the internal representation or where
adaptation occurs physiologically.

This is a first pass. Every rule assumes complete adaptation, which several of these experiments did not have. Cone
signals are Hunt–Pointer–Estévez with equal-energy normalisation, and CAT02 shows that the choice of sensors matters about
as much as the difference between rule R and the centring filter. Rule F depends on the disc centre and on the choice of
rim states; rule R does not. Only chromaticity is compared.

## The working realization and its place in LG

**Geometric gains are the form that fits a fixed map per condition.** The linear-gain interpolation is kept as a
benchmark, but it does not compose: going from condition 1 to 2 and back does not return the starting colour (for strength
1/2 and gain ratios 2, 1, 1 the round trip multiplies the first channel by 1.125). So it cannot come from equating the
outputs of one fixed map assigned to each condition. The geometric rule composes exactly, and it is what the family

    A_c = ( I ⊗ diag(w_{c,i}^(−D/2)) ⊗ I ) A

gives when normalised colour outputs are equated across two conditions, with w_c the cone catches of the adapting white.
(The distinction and this form of the family are the first reviewer's.) Its held-out error is 0.0119 against 0.0196 for
the centring filter. It is a candidate realization; equating model states with appearance matches is an assumption.

**What the family supplies**, stated against LG (`manuscripts/lg/LG_Master_Unified.tex` in AlphaOmega-Theory):

- For each condition a linear, injective map, so `ax:physical-presentability` (Ker A = 0) holds. A_c*A_c is
  wavelength-diagonal and positive, and the polar isometry, hence the final support, depends on the condition.
- A nonzero Ker A_c*, so `ax:strict-phenomenal-excess` holds formally. Per wavelength the output space has dimension 6
  and the range dimension 1. From 390 to 615 nm the joint span of the ranges of all conditions has dimension 3, so the
  5 remaining dimensions split into 2 inside that span and 3 outside. From 616 to 830 nm, where the table has S = 0, the
  joint span has dimension 2 and the split is 1 and 4: no gain activates a zero cone component. "Inside the joint span"
  means a linear combination of preparations, not a direction that one allowed preparation reaches.
- Three positive cone observables and two opponent observables whose compressions to the joint span commute, so the
  static statistics of these measurements on every output are classical. The opponent observables themselves do not
  preserve that span. Whether measurement-update operations are available, and what they would do, is additional.

**What remains additional:**

- *Ground data.* The construction has no generators, algebras or ground vectors, so `ax:physical-vacuum`,
  `ax:phenomenal-ground` and `ax:ground-bridge` are not addressed. A vector whose image is to lie on one ray for every
  condition must be a common eigenvector of all the gain operators. A vector with equal amplitudes on the three cone
  levels is not, and is moved by the gains. But an extra register level is not needed merely to obtain a fixed grey ray:
  two wavelengths with S = 0, both on the S level and carrying orthogonal two-level vectors, give a fixed ray whose
  reduced colour state is maximally mixed and which is orthogonal to every preparation. That ray has S-channel
  probability one, so it is not a no-absorption state. (An earlier statement here that a single register level forces a
  rim colour state overlooked entanglement with wavelength; the counterexample is the first reviewer's.)
- *Dynamics.* None is specified, so a dynamically accessible defect is not defined. The 2 + 3 split above is a split by
  conditions, not by dynamics.
- *LG's indefinite structure.* The two-level factor here is an ordinary positive-definite one. No relation to LG's
  ground-plane operators or sector grading is established.
- *The bridge to appearance*, the strength D, and any context beyond the adapting white.
- *Complete output states.* The map keeps wavelength labels, so metamers differ as complete states. Whether those
  differences are phenomenal is an interpretation.

### Candidate extension: a no-absorption pair

The three cone projections sum to the identity on the output space, so every normalised output state registers in some
channel and none has all three cone probabilities zero. An explicit no-absorption outcome therefore needs an added
sector. The smallest one adds a dark input ray |d⟩ and an output ray |Ω_d⟩, left unchanged by adaptation:

    Â_c = a₀ ⊕ A_c,   F̂_i = 0 ⊕ F_i,   a₀ > 0.

Then Â_c|d⟩ = a₀|Ω_d⟩ and Â_c*|Ω_d⟩ = a₀|d⟩ for every condition, each Â_c is injective, all three cone readouts vanish
on the added ray, and the added pair creates no further defect because it is paired. (This formulation is the first
reviewer's; the script checks it.)

What this does not supply: generators and cyclic algebras that would make the pair LG's physical and luminous grounds;
an identification of "no incident light" with LG's physical vacuum, which LG does not itself make; or any statement
about what darkness looks like. Extending the opponent readouts by zero gives neutral selected readings on the dark ray;
that is not the same as a maximally mixed two-level state.

## Scope and corrections

- White balance alone does not select the three rim states. Requiring three pure conditional states whose equal mixture is
  the disc centre does, up to rotation and reflection. That is a modelling choice. A different channel scaling gives a
  different disc.
- The enclosing ellipses are approximate: the solver stops within about 10⁻⁵ of enclosing and the ellipse is then enlarged
  to enclose. "On the rim" means within 10⁻⁴. An earlier remark that the collinearity obstruction also holds for a merely
  locally straight locus was wrong and is removed; the obstruction uses the exact tabulated zeros.
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
