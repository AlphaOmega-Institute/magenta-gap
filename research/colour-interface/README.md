# Colour interface: algebra checks

`check_colour_interface.py` contains 11 standard-library checks supporting the effect-body and paired-carrier analysis in the self-contained [working manuscript](../../manuscript/Magenta_Gap.tex). They check finite algebraic identities and counterexamples; they are not psychophysical measurements or substitutes for the proofs.

## Run

From the repository root:

```sh
python research/colour-interface/check_colour_interface.py
```

No input data or third-party packages are required. The script uses exact rational arithmetic where specified and floating-point checks for complex phase identities.

## What is checked

| Check group | Content |
| --- | --- |
| Effect body | Positivity against principal minors; illustrative eigenvalues; uniform recalibration and reference points |
| State coordinates | Pure-vector interval identity; determinant term for mixed states |
| Symmetries | Conjugation and complementation of the bicone |
| Fixed bridge | Scalar active-graph compression; input-independent sector-phase statistics; a counterexample to extending the ground-plane characterization to larger carriers |
| Readout | Defect extensions indistinguishable on paired preparations but distinguishable after range mixing; internal paired effects with input-dependent statistics |

## Algebra and interpretation

On the ground plane, let \(J=\sigma_z\), \(R_0=\sigma_x\), and \(Q_0=\sigma_y\). The effects

\[
E(x,y,z)=\tfrac12((1+z)I+xR_0+yQ_0)
\]

satisfy \(0\le E\le I\) exactly when \(\sqrt{x^2+y^2}+|z|\le1\). This three-dimensional slice is unitarily equivalent to the real two-level effect body; the unrestricted complex-qubit effect space has four real coordinates. Effects assign probabilities to states. Their zero and unit operators are not the two sector ground vectors, and an effect-body model does not identify the optical gamut with its boundary.

For an isometry \(u_A\), define \(W_\theta\psi=2^{-1/2}(\psi\oplus e^{i\theta}u_A\psi)\). Extending the sector operators over the active carrier gives

\[
W_\theta^*E(x,y,z)W_\theta
=\tfrac12(1+z+x\cos\theta+y\sin\theta)I.
\]

Thus this family has input-independent statistics on a fixed graph. Internal effects can instead carry input-dependent content while preserving the bridge:

\[
W_0^*\operatorname{diag}(F,u_AF u_A^*)W_0=F.
\]

A stimulus encoding and justified choice of readout remain necessary for a colour realization. On larger carriers, the condition \(E+JEJ\) scalar does not characterize only the three-parameter family above.

With \(Q=I-u_Au_A^*\), every effect \(u_AF u_A^*+D\), \(0\le D\le Q\), has the same paired statistics. Dynamics confined to the range preserves that indistinguishability; dynamics mixing the range with its complement can distinguish the extensions. A proper convex preparation image need not produce any linear adjoint-kernel defect.

## Height calibration

Opponent contrast, effect magnitude, and a white-normalized cone sum are different height coordinates. Using the same symbol does not identify their endpoints or scales.

The illustrative points \((r,z)=(1/2,\pm1)\) and \((13/10,0)\) lie outside the unit effect bicone. This constrains a model only if it assigns those coordinates to appearances and identifies ordinary black and white with the bicone tips. The points are not measured locations of altered colours.

For the recalibrated family

\[
E_s=\tfrac12((1+z/s)I+(x/s)R_0+(y/s)Q_0),\qquad s>0,
\]

the effect condition \(0\le E_s\le I\) holds exactly when \(s\ge |z|+r\). Recalibration changes where reference points map; it does not relax probability positivity. An ordinary preparation image can also be a proper subset of an effect space. Neither altered-colour descriptions nor these checks select a calibrated perceptual target.

## Sources

- [LG manuscript](https://github.com/AlphaOmega-Institute/AlphaOmega-Theory/blob/main/manuscripts/lg/LG_Master_Unified.tex): ground Pauli operators, active factorization and paired effects.
- M. Berthier (2020), *Geometry of color perception. Part 2: perceived colors from real quantum states and Hering's rebit*, Sections 5.1 and 7.3–7.4, [doi:10.1186/s13408-020-00092-x](https://doi.org/10.1186/s13408-020-00092-x).
- M. Berthier and E. Provenzi (2022), *Quantum measurement and colour perception: theory and applications*, [author preprint](https://arxiv.org/abs/2504.12304).

The two-level colour geometry is prior art. The separate [continuum boundary certificate](../object-colour/Continuum_Rim_Certificate.md) concerns a specified bounded-reflectance signal model, not a perceptual boundary.
