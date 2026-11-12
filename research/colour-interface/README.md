# Colour interface: effect algebra and calibration

Research record, 12 November 2026. This records the agreed mathematical results of the colour-interface audit. It is not a completed physiological realization, a manuscript revision, or a claim of new two-level/Jordan algebra. The July and September manuscript sources and the [continuum boundary certificate](../object-colour/Continuum_Rim_Certificate.md) remain unchanged.

## 1. The exact algebraic connection

On LG's ground plane, write the Pauli triple as

\[
J=\sigma_z,\qquad R_0=\sigma_x,\qquad Q_0=\sigma_y.
\]

The three-dimensional family

\[
E(x,y,z)=\tfrac12\big((1+z)I+xR_0+yQ_0\big)
\]

has eigenvalues \((1+z\pm r)/2\), where \(r=\sqrt{x^2+y^2}\). Consequently,

\[
0\le E\le I\quad\Longleftrightarrow\quad r+|z|\le1.
\]

This is precisely the symmetric bicone. Its lower and upper sheets are, respectively, \(\det E=0\) and \(\det(I-E)=0\). A unitary change of basis sends the span of \(I,R_0,Q_0\) to the real symmetric two-by-two matrices. This is the real two-level effect body described by Berthier (2020), Section 5.1, especially Equations (79)–(80); see also Berthier–Provenzi (2022). A full complex-qubit effect space has four real coordinates. The three-dimensional slice is part of the construction, not a property of every unrestricted two-level effect space.

An effect is an operator assigning a probability to each state. Identifying these effects with perceived colours is a further model choice. The zero and unit effects are not LG's physical and phenomenal ground vectors. Nor does the existence of this effect body identify the independently calculated optical gamut with its boundary.

On the ground plane, the family is characterized among effects by \(E+JEJ\) being a scalar multiple of the identity. Conjugation by \(J\) sends \((x,y,z)\) to \((-x,-y,z)\); complementation \(E\mapsto I-E\) reverses all three coordinates. Their composition sends \(z\) to \(-z\) at fixed \(x,y\). These exact symmetries of the ideal effect body do not establish corresponding symmetries of a physiological preparation image.

## 2. Vectors, states and the fixed bridge

For a vector \(v\) in the ground plane, put

\[
Z=v^*v,\quad X=v^*R_0v,\quad Y=v^*Q_0v,\quad T=v^*Jv.
\]

Then \(Z^2-X^2-Y^2=T^2\). Normalizing \(Z=1\) gives a disk after omitting the sign of \(T\). The neutral pure rays form its boundary circle; the two normalized sector grounds both project to its centre. This quadratic construction is not the restriction of the Krein form to LG's fixed neutral graph: that restriction vanishes identically.

For a positive two-by-two matrix \(M\), with the corresponding trace coordinates,

\[
Z^2-X^2-Y^2=T^2+4\det M.
\]

Equal sector weights therefore do not by themselves imply saturation. For example, \(M=I/2\) has \(T=X=Y=0\), \(Z=1\), and nonzero determinant.

Let \(u_A\) be the polar isometry between active supports. LG factorizes the active carrier as a sector coordinate times the active physical carrier, with

\[
R_A=\sigma_x\otimes I,\quad Q_A=\sigma_y\otimes I,
\quad J=\sigma_z\otimes I.
\]

The same effect formula, using the active identity, has the same bicone. For normalized \(\psi\), define

\[
W_\theta\psi=2^{-1/2}(\psi\oplus e^{i\theta}u_A\psi).
\]

Its compression is

\[
W_\theta^*E(x,y,z)W_\theta
=\tfrac12(1+z+x\cos\theta+y\sin\theta)I.
\]

For one fixed graph, this entire family has statistics independent of the physical input \(\psi\). Operational access to the sector phase alone therefore does not supply stimulus-dependent colour. Interpreting hue as that phase additionally needs a justified preparation or control rule. On a larger active carrier, the condition \(E+JEJ\) scalar permits arbitrary operator-valued off-diagonal blocks; it no longer characterizes only this three-parameter family.

A different route encodes colour in internal physical state content. For any admissible physical effect \(F\), the paired effect satisfies

\[
W_0^*\operatorname{diag}(F,u_AF u_A^*)W_0=F.
\]

These statistics can depend on \(\psi\) while the bridge remains fixed. An internal two-Pauli colour-effect family can itself have a bicone. The choice of internal observables, stimulus encoding, adaptation variables and empirical readout remains additional. This is a representation route, not an established colour realization.

## 3. What paired data say about the defect

Set \(P=u_Au_A^*\), \(Q=I-P\). Static paired data fix only the compression of a phenomenal effect:

\[
u_A^*E_\Phi u_A=F.
\]

The canonical extension \(u_AF u_A^*\) vanishes on the defect, but the data do not force that choice: every

\[
E_\Phi=u_AF u_A^*+D,\qquad 0\le D\le Q,
\]

has the same paired statistics. More general extensions can also have cross terms, subject to positivity. In particular, adding \(tQ\) preserves the effect condition for \(0\le t\le1\), not for arbitrary real \(t\).

These extensions remain indistinguishable for preparations and dynamics confined to the paired range. Independently justified dynamics that reaches the defect can distinguish them; permanent invisibility does not follow. The accompanying checks include a finite range-mixing example.

An ordinary preparation image can be a proper convex subset while spanning the full carrier. That does not produce an adjoint-kernel defect. If actual ordinary and fatigue-induced appearances are both represented through one fixed \(u_A\), both have zero defect weight, even if their preparation families differ. The mere fact that an observer is physically embodied does not independently establish that representation assumption. LG's identification of an aligned ray as one cognitive event must not be read, without a scope stipulation, as its converse for every phenomenal determination.

## 4. Three different height coordinates

| Construction | Meaning of height | What is fixed |
| --- | --- | --- |
| July manuscript, opponent chart | \(z=A_{W/B}-50\), a centred white/black opponent activation; its discussion concerns focal-versus-surround contrast | Opponent-coordinate bounds; these do not automatically calibrate ordinary reference black and white in every context |
| Two-level effect body | \(z=2a_3-1\), where \(a_3=\operatorname{Tr}E/2\) is effect magnitude | Order bounds \(0\le E\le I\); identification with luminance, brightness or opponent activation is additional |
| Certified bounded-reflectance signal model | \(z=\tfrac23(a_L+a_M+a_S)-1\), for white-balanced cone catches | Perfectly absorbing and perfectly reflecting surfaces under the specified illuminant map to \(-1,+1\) |

These are not identified by using the same letter. The July manuscript already separates an ambient cylindrical response model from an ordinary spindle inside it. The cylinder permits chromatic variation at its height extrema. That separation is existing manuscript content, not a new outcome of this audit. The symmetric effect bicone is a different candidate geometry unless a map relating the two is specified.

## 5. Conditional altered-colour constraint

The following points are illustrative coordinates, not measured locations of altered percepts:

| Assumed point \((r,z)\) | Effect eigenvalues | Minimum uniform scale about mid-grey |
| --- | --- | --- |
| \((1/2,-1)\) | \(-1/4,1/4\) | \(3/2\) |
| \((1/2,1)\) | \(3/4,5/4\) | \(3/2\) |
| \((13/10,0)\) | \(-3/20,23/20\) | \(13/10\) |

These points are incompatible with the simultaneous assumptions that ordinary black and white are the zero/unit effect tips, that the same coordinate calibration assigns these points to altered appearances, and that all appearances inhabit this effect bicone.

Uniform enlargement recalibrates the map into effects:

\[
E_s=\tfrac12\big((1+z/s)I+(x/s)R_0+(y/s)Q_0\big),
\qquad s>0.
\]

It is an effect exactly when \(s\ge |z|+r\). At \(s=3/2\), ordinary achromatic reference points \(z=-1,+1\) map to \(I/6,5I/6\), rather than \(0,I\). This changes calibration, not probability positivity. A different target shape could instead permit chromatic variation at unchanged height extremes. Neither construction is selected by the illustrative points.

The statement “the full perceptual target is the effect space” does not imply that ordinary preparations fill that target. It permits an ordinary image \(C_{\rm ord}\subsetneq\mathcal E\) with altered responses in \(\mathcal E\setminus C_{\rm ord}\). Berthier (2020), Sections 7.3–7.4, explicitly distinguishes the ideal space from observer-specific convex subsets and leaves contextual effects for further work. Therefore this audit does **not** establish that the Berthier–Provenzi approach leaves no room for altered colours. Conversely, it supplies no calibrated account of them.

Descriptions such as “self-luminous” or “as dark as black” do not alone assign exact values to an effect's trace coordinate. Shared calibration and an operational response map are needed for a cross-regime comparison. No CFT, measured physiological metric, or cognitive final defect is inferred here.

## 6. Where this belongs in the manuscripts

- **September revision:** include the exact effect-bicone connection with prior-art attribution; distinguish the three heights and the ordinary image from the full target; state the conditional altered-colour incompatibility; incorporate the separate invariant null-boundary obstruction with its signal-model scope. Do not include the withdrawn assertion that Berthier–Provenzi excludes altered colours.
- **July source:** preserve its imported historical state pending a deliberate substantive revision. Its cylinder/spindle distinction is already present; its stronger physiological and ridge claims still need the earlier corrections.
- **LG/TAUM:** no core theorem is repaired by this note. Any colour-realization discussion should distinguish content from sector phase, preparation restrictions from linear defects, and paired events from general phenomenal possibilities.

The manuscripts are to remain self-contained: mathematical statements and explanations needed by a reader should be written into their own TeX, not deferred to this research record. Nothing in this commit applies those manuscript revisions.

## Reproducible checks and sources

Run from the repository root:

```sh
python research/colour-interface/check_colour_interface.py
```

The checks verify finite algebraic identities, exact rational examples and counterexamples. They are regression checks supporting the displayed arguments, not empirical tests or numerical substitutes for the proofs. Claude's reported 200,000-point tests and scratch scripts have not been imported; this is a separate reproducible check set.

Sources checked for this record:

- LG, repository blob `615da8db989cec99099f85bcd8431b634bf51e9e`, especially the ground Pauli triple, active factorization, paired effects, operational structural-phase question and phenomenal-excess question: [source](https://github.com/AlphaOmega-Institute/AlphaOmega-Theory/blob/main/manuscripts/lg/LG_Master_Unified.tex).
- July and September Magenta sources at magenta-gap commit `3734ac9e9abf5528b155c7dd36fde0a659755776`; the imported July source already contains the ambient cylinder and ordinary spindle.
- M. Berthier (2020), *Geometry of color perception. Part 2: perceived colors from real quantum states and Hering's rebit*, Sections 5.1, 7.3–7.4, [doi:10.1186/s13408-020-00092-x](https://doi.org/10.1186/s13408-020-00092-x).
- M. Berthier and E. Provenzi (2022), *Quantum measurement and colour perception: theory and applications*, [author preprint](https://arxiv.org/abs/2504.12304). The two-level colour geometry is prior art, not a novelty claim for this programme.
