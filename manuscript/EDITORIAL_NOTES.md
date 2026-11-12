# Editorial notes: JOSA-format Magenta Gap review draft

28 September 2026

## Status and provenance

This package contains an eight-page JOSA-format **author-review draft**, not a certification of submission readiness, scientific acceptance, or compliance with every publisher requirement. The source is the author's 51-page `Magenta_Gap_MSS(1).pdf`, read in full. The earlier JOSA architecture mockup informed the organization, but its provisional tables and editorial placeholders were not treated as established results. The original manuscript remains unchanged.

The condensation preserves the original argument: distinguish the full spectral stimulus domain from its monochromatic locus; represent nonspectral hue by a purple interval inserted into an auxiliary compactified spectral arc; collapse hue on the achromatic axis; introduce an indefinite metric as a proposed reachability gauge; and derive a conditional compatibility constraint for a prescribed gamut boundary. The paper remains *The Magenta Gap*, with its original geometric and explanatory aim.

The September cubic object-color theorem has **not** been inserted. It concerns bounded-reflectance attainability, uses different signal and mixture coordinates, and does not establish the original manuscript's phenomenal or indefinite-metric claims. Incorporating it would require a substantive new scientific section, not ordinary condensation.

## Ordinary editorial work

Repeated introductions, repeated definitions, extended philosophical restatements, and speculative physiological bridges were compressed. The response map, nonspectrality criterion, purple gluing, achromatic quotient, split-complex section, null-boundary algebra, and distinction between compatibility and confinement remain. The abstract, captions, section transitions, and conclusion were edited to describe the actual scope of the mathematical statements and numerical illustration. The manuscript now distinguishes colorimetric radius, modeled appearance coordinates, and measured discrimination.

The section mapping is:

| Original material | Revised location |
|---|---|
| Sections 1–3: motivation, spectral measures, compactification and resolution | Sections 1–2 |
| Section 4: quotient body, indefinite metric, split-complex sections | Section 3 |
| Section 5: spindle, null condition, bound, non-egress | Section 4 |
| Appendix A: cone-fundamental consistency calculation | Section 5, with recalculated methods and results |
| Sections 6–7: interpretation, limitations and conclusion | Sections 6–7 |
| Section 5.4: regular interior construction | Appendix A, with restricted scope |

## Substantive mathematical corrections

These changes go beyond wording and require the author's explicit scientific review.

**General versus diagonal bound.** The perfect-square identity yields the non-strict condition in the general radial–angular spatial block. The strict condition and formula for the radial coefficient apply only to the specified diagonal ansatz. Choosing that ansatz is not presented as a coordinate change that automatically preserves the boundary and calibration. A finite cross-term example shows that equality can be regular. Failure of the diagonal formula at equality is no longer called a universal geometric horizon or removable coordinate singularity.

**Ridge hypotheses.** The obstruction to a smooth null radial envelope now states its radial-graph, endpoint, and metric hypotheses, including the absence of spatial–luminosity cross-terms. The argument forces some failure of differentiability under those assumptions. The isoluminant equatorial location follows from the assumed symmetric bicone, not from nullity alone.

**Optical non-egress.** The control condition is imposed at the boundary on physically admissible one-sided velocities. It is not imposed throughout a collar on the entire linear image of a differential, which would include both a vector and its negative. The identity showing that the metric cancels remains central. Full invariance still needs a response/control model and suitable evolution assumptions.

**Regular extension.** The smooth extension claim is local on compact patches of open boundary sheets. A Cartesian interior core is distinguished from the boundary tips. For nonconstant angular ceiling, the prescribed luminosity coefficient has direction-dependent limits at those tips, obstructing a continuous extension there. No global smoothness claim conceals that obstruction.

## Substantive numerical revision

The new numerical analysis is not merely a reformatted version of the old appendix. It uses the convex hull of the selected 390–700 nm tabulated response vectors as the primary additive mixture boundary. That hull has 188 vertices; its long mixture edge joins the retained 410 and 700 nm vertices. Treating every monochromatic sample as an outer boundary point would include interior portions of the locus.

The calculation uses analytic derivatives on open polygon edges and finite sampling within them. Corner derivatives are one-sided. The summaries retain finite points satisfying the declared filter beta < 10 at metric scale one; this is a computational convention, not a physiological definition of nondegeneracy. All retained primary samples satisfy the strict condition, but the finite grid is not an all-point certificate.

The original locus-plus-chord construction remains a control. Finer interpolation resolves a narrow raw-normalization failure near 403 nm that the coarse angular grid missed. It lies inside the mixture hull, so it does not contradict the new outer-boundary calculation. Both results remain documented.

Equal-energy balancing is explicit. The Fisher coefficient is the fixed-radius angular component, not the pullback along the changing-radius boundary. Its scale remains a modeling choice. Reimposing the Fisher formula after channel normalization changes the metric identification as well as coordinates. The manuscript therefore distinguishes hue-coordinate covariance from discrimination-scale independence. The earlier MacAdam comparison is not retained as independent validation.

## References and author review

The bibliography retains sources supporting claims that remain in the draft. The short-cone study by Stockman, Sharpe and Fach (1999) accompanies the Stockman–Sharpe (2000) middle/long-cone citation; the latter alone does not document all three cone channels. Shevell's 2003 edited volume is distinguished from the different 2008 review that appeared under an inconsistent key in the prior mockup. Removed speculative bridges no longer borrow authority from unrelated empirical citations. Bibliographic correction does not establish novelty of the model.

Two DOI errors in the earlier material were identified: Resnikoff (1974) is `10.1007/BF00275798`, and Krantz (1975) is `10.1016/0022-2496(75)90026-7`. The revised bibliography uses the verified publication details and omits the incorrect DOI strings. Bujack et al. is cited for the limits of a global Riemannian account, not as evidence for an indefinite metric; Churchland's chimerical-color examples are distinguished from Crane and Piantanida's forbidden-hue reports.

Before submission, the author should accept or revise the substantive corrections, check the code and results, verify every retained claim and citation, and confirm the data-access statement. Funding and conflicts must be supplied from actual facts; none should be inferred from silence. Verify the reported AI tool names, versions, model identifiers, sources and uses, including the original manuscript's metadata. The assistance included audit and code development as well as editing, and should not be described as grammar assistance alone. The figure caption records code assistance; no new physiological measurements are claimed. Final disclosure and publisher-policy decisions remain the author's responsibility.

The current session identifies the assistance as OpenAI ChatGPT/Codex; an exact deployment/model version was not exposed and has not been guessed. The revised code and results accompany this package; they have not been published to the original project repository. The affiliation statement is retained as a fact, not treated as a completed declaration of all financial interests or funding.

---

## Consolidation of 12 November 2026

This section updates the historical notes above. The source remains the shorter September manuscript, now revised for author and Claude review. The July JMP source is unchanged. The revision retains the conditional angular compatibility theorem, equality counterexample, diagonal construction, original numerical table and three figures. It does not select a new physiological geometry or claim an empirical colour realization.

### Incorporated mathematics and scope

- **Section 3.C:** the standard real two-level effect bicone, with its eigenvalue proof and Berthier (2020)/Berthier–Provenzi (2022) attribution. Effect magnitude, opponent contrast and bounded-reflectance signal height are separated. The existing cylinder/spindle distinction is preserved. The altered-colour tip incompatibility is stated conditionally, with an illustrative recalibration; no claim that the effect approach excludes altered colours is retained.
- **Section 5.C:** chromatic-only rescaling at fixed height calibration is distinguished from whole-metric rescaling. The Poisson calculation shows why photon number scales both information components. A possible common calibration of measured angular and height sensitivities is stated as an additional bridge, with its constant-ratio condition and resulting inequality.
- **Section 6:** a coordinate-free conormal obstruction including metric cross-terms; a supporting-normal proposition for at least two independent genuine reflectance transitions; and the bounded-reflectance model with certified root and point enclosures. Differentiability is asserted at the certified point, not on a smooth neighbourhood. The 390–830 nm continuum reflectance model is explicitly different from the earlier 390–700 nm additive hull. Sampled hue percentages are not used.
- **Appendix B:** the exact polynomial obstruction to identifying the whole spatial block with a common multiple of the simplex Fisher tensor on a nontrivial open straight edge segment in the simplex interior. Singular edges where that tensor is undefined are outside the proof. The angular-only construction remains intact.
- **Appendix C:** the paired-carrier definitions and compression calculation are stated directly. Fixed sector effects have input-independent statistics on one aligned graph; internal effects can carry varying content. Paired statistics do not fix an effect on the orthogonal complement. A convex reachability restriction is distinguished from a linear defect.

The abstract, introduction, discussion and conclusion now reflect both the established algebraic connection and the certified negative result. The stimulus domain is described as an admissible set, since bounded reflectances are not a cone. Hyperbolic rotations are not asserted to preserve a radial half-plane. Four prior-art references were added, with the supporting claims checked in the primary texts. The three-dimensional effect family is not presented as new algebra.

The separate September cubic object-colour theorem and the exploratory Hilbert-gamut comparison are not inserted. Their distinct model assumptions would require a separate decision about scope. No LG/TAUM manuscript changes are made by this draft.

### Verification and reproducibility

- The exact continuum certificate rerun is byte-identical to the committed JSON at `cb997e8`.
- All 11 colour-interface regression checks pass.
- The preserved independent branch `claude/independent-colour-checks` at `581ff68` changes only its new research directory. All 24 executions complete, and their scientific output matches the preserved output at printed precision. Most scripts print results rather than assert pass/fail conditions; they are not certified numerical proofs. Runtime and absolute-warning-path differences were excluded from the comparison.
- The revised TeX compiles with pdfLaTeX/latexmk to 11 pages. It has no undefined references or citations, no LaTeX errors and no overfull boxes. Four underfull-line notices remain; page images were inspected for layout. The local incomplete TeX configuration was repaired using installed format sources and standard font maps. The manuscript requires no runtime-specific paths.
- All existing figure and layout-support files match their repository blobs. The original numerical table is unchanged. The text and bibliography remain in one TeX file; no companion note is needed to follow the arguments. Code and data locations are supplied for computational reproduction.

The original import source, PDF and provenance remain in Git history at `cb997e8`. The import checksum in `PROVENANCE.json` is historical, not the checksum of the newly revised source. The working filename retains the September date for continuity. This is an author-review draft, not a new submission or a frozen DOI version.

### Follow-up to Claude's review of `748df23` (12 November 2026)

The transition-normal proposition now requires **two of** the transition response vectors to be independent, allowing additional transitions in the same supporting plane. The complement map `rho -> 1-rho` is stated explicitly: it sends `a -> W-a` and, in the declared affine chart, `(x,y,z) -> (-x,-y,-z)`. A common ceiling on both polar tangent cones therefore requires `alpha(theta+pi)=alpha(theta)`. This is necessary for matching those cones, not a derivation of the separable bicone or same-hue mirror symmetry.

All four citation locators are moved into the prose, bibliography entries follow first-citation order, and redundant equation-name prefixes are removed to respect the class's definition of `eqref`. The cylinder/spindle sentence now cites both the ambient quotient and the ordinary image.

The full-Fisher appendix clarifies that `p`, `r` and `d` belong to the selected slice. Its formula already applies at any fixed open-sheet height because `R_z^2/C=1`. Moving to another slice generally changes the Fisher tensor through `p`; it is not merely multiplication of the entire expression by a constant. The nonconstancy proof applies independently on each regular straight slice.

The updated PDF remains 11 pages. A fresh build has no errors, unresolved references/citations or overfull boxes; the same four cosmetic underfull-line notices remain. Checks on extracted PDF text find no doubled equation names or locators inside citation brackets. The bibliography order and all labels were checked, and the rendered pages inspected. No numerical values, certificate code or inputs were changed by this follow-up.

