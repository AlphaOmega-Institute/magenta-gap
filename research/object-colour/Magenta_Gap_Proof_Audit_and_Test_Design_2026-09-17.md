# The Magenta Gap: proof audit and candidate measurement design

Research follow-up — 17 September 2026

## Findings

The conditional theorem survived this audit. I checked the proof's logical steps, supplied a separate dual derivation of its lower bound, and constructed a polynomial example that tests both the interior relation and the necessity of qualifying endpoint claims. This is a further audit by the same assistant, not an external referee's review or a formal proof-assistant verification.

The additional paper search found a public transcription of all three pages of West and Brill (1983). The 2015 color-cone paper and the book chapters remain unavailable in full text. Novelty remains unresolved.

The numerical screening identifies concrete test designs in the declared analytic spectral models. It does not establish that an available instrument can implement or resolve them. The PCHIP model needs much smaller gaps for a comparable test of both terms; the practical difference between models is substantial.

## 1. Proof audit

The claim being checked is

\[
d_\tau(z)=C(\tau)z^3+o(z^3),\qquad
C(\tau)=A_0(1-\tau)^3+A_1\tau^3,
\qquad A_i=\frac{q_i}{24b_i^3}.
\]

Its hypotheses specify a continuous response density on a compact spectral interval, a strictly positive white-normalized signal functional, and exactly two distinct interior quadratic supporting contacts. The local response is twice continuously differentiable. The limit holds for each fixed interior mixture; endpoints require the stated E1 or E2 condition. No uniform joint endpoint/signal limit is claimed.

| Step | Audit result |
|---|---|
| Existence of an optimum | The bounded reflectance class is weak-* compact; its three response integrals are continuous. The fixed-signal, fixed-ray feasible set is nonempty and compact. |
| Reduction to two moment constraints | White and the two contact chromaticities form a basis. Signal and one transverse linear moment specify the ray, leaving the deficit as the objective. |
| Change to local signal mass | With \(dt=b\,d\lambda\), the supporting cost is \(12A_i t^2+o(t^2)\). The factor \(1/24\) follows from integrating over a centered band. |
| Interior upper bound | The two-band mass Jacobian is invertible. Exact masses differ from \((1-\tau)z,\tau z\) by \(O(z^3)\); they remain positive at any fixed interior mixture. |
| Lower bound | The concentration and bounded-density filling argument is valid. The independent dual construction below gives the same coefficient without assuming a two-band optimizer. |
| Endpoint E1 | The averaged hue function can be written as \(\int_{-1/2}^{1/2}h_i(P_i(s+zu))du\), which is \(C^2\). Its nonzero derivative in \(s\) gives the required \(O(z^2)\) shift. |
| Endpoint E2 | The implicit correction mass is \(\Gamma_i z^3+o(z^3)\). Strict positivity of \(\Gamma_i\) makes it admissible; its own supporting cost is \(O(z^9)\). |
| Scope of invariance | Consistent linear coordinates and white-normalized intensity changes preserve the quantities. A change of spectral variable must be a smooth local diffeomorphism with the density transformed as well. |

No missing hypothesis or invalid inference was found in the theorem as currently qualified. Its application to fitted spectral data still depends on those models satisfying the hypotheses; the existing numerical contact checks are not interval-arithmetic certificates or measurements of physiology.

### A separate dual lower bound

Write \(d\mu=b\,d\lambda\), \(p=c/b\), \(g=n(p)\), and \(h=h_\tau(p)\). For any real multipliers \(\alpha,\beta\), every feasible reflectance obeys

\[
d_\tau(z)\geq \alpha z+
\int \min\{0,g-\alpha-\beta h\}\,d\mu.
\tag{A1}
\]

This follows directly from the two moment constraints and \(0\leq\rho\leq1\). Strong duality is unnecessary.

Set \(w_0=1-\tau\), \(w_1=\tau\), and choose multipliers so that

\[
\alpha+\beta h_\tau(p_i)=3A_iw_i^2z^2.
\]

The two contact values of \(h_\tau\) differ, so this system is solvable. Their weighted mean is zero, hence \(\alpha=3C(\tau)z^2\). Outside contact neighborhoods, \(g\) is bounded away from zero and the negative part in (A1) eventually vanishes. Rescaling each local mass coordinate by \(z\), its negative-part integral is

\[
\int \min\{0,12A_i t^2-3A_iw_i^2z^2\}\,dt
=-2A_iw_i^3z^3.
\]

Taylor remainders contribute \(o(z^3)\). Thus (A1) gives \(d_\tau(z)\geq C(\tau)z^3+o(z^3)\), matching the feasible upper construction. This separately checks the coefficient and the direction of the optimization inequality.

## 2. An explicit test of the endpoint qualification

Use a dimensionless spectral parameter \(\lambda\in[-2,2]\), and define

\[
u(\lambda)=\frac{(2-\lambda)(\lambda+1)^2}{4},\qquad
g(\lambda)=(\lambda^2-1)^2+(\lambda^2-1)^4,
\qquad c(\lambda)=\frac14(u(\lambda),g(\lambda),1).
\]

All three response components are nonnegative. Put

\[
G=\frac14\int_{-2}^{2}g(\lambda)d\lambda=\frac{2606}{315},
\quad W=(1/2,G,1),\quad \ell(X,Y,Z)=Z,\quad n(X,Y,Z)=Y/G.
\]

The two supporting contacts are \(-1,+1\), with chromaticities \((0,0,1)\) and \((1,0,1)\). Both are stationary in projective coordinates, and both supporting zeros are quadratic. The local coefficient is

\[
A_0=A_1=\frac{16}{3G}=\frac{840}{1303}
\approx0.6446661550.
\]

At the first endpoint, \(h_0=X-Y/(2G)\). E1 fails, and
\(\Gamma_0=-883/1303<0\), so E2 fails as well. The unqualified claim \(d_0(z)\sim A_0z^3\) is false. In fact,

\[
\boxed{d_0(z)=2z^3+o(z^3).}
\tag{A2}
\]

To see the lower bound, note that feasible endpoint responses satisfy \(n(x)=2X\). For any sequence with cost \(O(z^3)\), the original supporting cost confines the leading mass to the two contacts; the endpoint moment then forces that leading mass to \(\lambda=-1\). In the local mass coordinate \(t=(\lambda+1)/4\),

\[
2u(-1+4t)=24t^2-32t^3.
\]

The filling bound therefore gives leading coefficient at least 2.

For the upper bound, take a centered main band of mass \(z-e\) at \(-1\), plus a band of mass \(e\) ending at \(\lambda=2\). The latter has negative \(h_0\) and supplies the correction unavailable at the other supporting contact. The implicit function theorem gives

\[
e=\frac{G}{45}\left(1-\frac{A_0}{2}\right)z^3+o(z^3)>0.
\]

The main band's \(X\)-response is exactly \((z-e)^3\); the correcting band's is \(O(e^2)\). Since \(d_0=2X\), this attains (A2). The parameter interval can be affinely mapped to positive physical wavelengths; its dimensionless choice is immaterial.

This example satisfies the interior hypotheses and disproves the endpoint statement only when its additional conditions are removed. It is not a counterexample to the current theorem.

### Numerical checks independent of the spectral optimizer

The script `proof_stress_audit.py` uses polynomial antiderivatives, exact moment constraints, and the global weak-duality bound (A1). It does not call the spectral boundary optimizer.

For \(\tau=0.1,0.5,0.9\) at \(z=0.01,0.003,0.001\), all nine feasible-upper/dual-lower brackets are valid. At \(z=0.001\), every bracket lies within 0.0011% of the predicted leading value; the widest bracket is about 0.000243% of that value. The tested curve has stationary tangents at both contacts, so it also checks that the interior theorem does not silently rely on E1.

The endpoint constructions give \(d/(2z^3)=0.9999627141,0.9999966382,0.9999996263\). Their ratios to the incorrect coefficient tend instead to \(1303/420\approx3.10238\). The proof of (A2), rather than numerical extrapolation alone, establishes the limit.

## 3. Paper retrieval

The [public West–Brill transcription](https://www.academia.edu/54927465/Conditions_under_which_Schr%C3%B6dinger_object_colors_are_optimal) contains the text and captions from pp. 1223–1225. I read the barycentric argument, convexity/wavelength-ordering condition, and general support-line construction. No local cubic mixing relation was found. The original PDF was not obtained; OCR interleaves columns, and figures were not visually verified. Its assertion about convexity of the human locus must be read alongside the later Burns and Davis analyses.

Additional repository/index checks still returned only restricted publisher or institutional records for [Logvinenko (2015)](https://doi.org/10.1364/JOSAA.32.000314) and [the book](https://doi.org/10.1002/9781119885955). The unresolved full-text comparison is the 2015 paper and chapters 14–16, pp. 357–388. This search did not establish novelty of the two-contact coefficient.

## 4. Identify mixtures that test both terms

The earlier \(\tau=0.25,0.5,0.75\) tests are strongly dominated by \(A_1\tau^3\). They are useful numerical checks but weak experimental tests of the short-contact term. For example, at \(\tau=0.5\) in the equal-energy analytic model, the short-contact share of the leading coefficient is only about 0.0073%.

The contributions are equal when

\[
\tau_*=
\frac{(A_0/A_1)^{1/3}}{1+(A_0/A_1)^{1/3}}.
\tag{A3}
\]

| Model | \(\tau_*\) |
|---|---:|
| Stockman–Sharpe / PCHIP | 0.0674706704 |
| Govardovskii / equal energy | 0.0400748224 |
| Govardovskii / red tilt | 0.0627040413 |

These are signal-share coordinates, not radiant-power shares or perceptually uniform hue positions. A diagnostic design uses \(\tau_*/2,\tau_*,2\tau_*\), so the two contributions trade dominance. Coefficients or endpoint deficits must be calibrated independently of the boundary data used to test the formula; fitting both amplitudes to those same data would weaken the test.

## 5. Candidate finite-signal designs

I screened 72 model/mixture/signal combinations. Sixty-five satisfy the numerical feasibility and full-spectrum support checks. Seven solver failures are retained in the raw records and remain unresolved outside the selected designs. Failure of a numerical solve is not evidence that the theorem fails.

| Model | Tested mixture coordinates | Selected signal levels | Worst cubic-coefficient error across the selected design |
|---|---|---|---:|
| PCHIP | 0.0337353, 0.0674707, 0.1349413 | 0.00003, 0.00005 | 0.9552% |
| Analytic / equal energy | 0.0200374, 0.0400748, 0.0801496 | 0.001, 0.002, 0.003 | 4.5755% |
| Analytic / red tilt | 0.0313520, 0.0627040, 0.1254081 | 0.001, 0.002, 0.003 | 2.0027% |

All 24 selected points pass the numerical checks. These are discrete tested designs, not certified error bounds over every intervening signal or mixture. At the equal-contribution point alone, both analytic models remain within 1% across \(z=0.001,0.002,0.003\).

The three equal-energy analytic mixtures at \(z=0.003\) illustrate the practical tradeoff:

| \(\tau\) | Cubic error | Relative cone shortfall \(d/z\), ppm | Short reflecting band, nm | Long reflecting band, nm | Allowed \(\delta d/z\) for 10% gap error, ppm |
|---:|---:|---:|---:|---:|---:|
| 0.0200374 | 0.0360% | 9.7128 | 0.6019 | 5.5674 | 0.9713 |
| 0.0400748 | 0.6890% | 16.2318 | 0.5896 | 11.0537 | 1.6232 |
| 0.0801496 | 4.5755% | 69.2483 | 0.5650 | 21.5377 | 6.9248 |

The corresponding absolute normalized deficits are \(2.914\times10^{-8}\), \(4.870\times10^{-8}\), and \(2.077\times10^{-7}\). A 10% measurement budget for the gap alone is therefore \(2.914\times10^{-9}\), \(4.870\times10^{-9}\), and \(2.077\times10^{-8}\) in the same white-normalized units.

The PCHIP design is much harder. At \(z=0.00005\), its three relative gaps are approximately **4.09, 6.61 and 28.83 parts per billion**. The corresponding 10% gap-error budgets are 0.409, 0.661 and 2.883 parts per billion of signal. Its short reflecting bands are only about 0.032–0.036 nm wide. Numerical agreement in this model cannot reasonably be presented as an already practical measurement protocol.

### Avoid isolated apparent agreement

At the equal-contribution mixture, the PCHIP coefficient ratio is about 0.7079 at \(z=0.003\), then returns to about 0.9534 at \(z=0.01\). The latter isolated agreement within 5% is not evidence of a stable asymptotic regime. The exported screening figure shows this explicitly. Connecting its markers is only a visual guide; no unsampled interval is certified.

This model dependence is a substantive limitation. The easier analytic-model designs do not establish the same practical window for the measured-fundamental interpolation, or for actual observers. A continuous spectral calibration and interpolation-sensitivity study would be needed before fixing a human-color experiment's numbers.

### Numerical accuracy versus measurement uncertainty

The selected analytic designs were checked by comparing 16- and 64-point Gaussian quadrature on identical band endpoints at their largest signal. The largest relative deficit difference was below \(1.5\times10^{-33}\). This checks integration error in these represented models; it says nothing about sensor calibration error. Global support checks for the analytic model remain numerical, not interval-certified.

To measure the gap itself to relative error \(\eta\), a necessary budget is

\[
|\delta d|\lesssim\eta d.
\]

For coefficient inference, the signal and mixture errors matter as well. A first-order worst-case budget includes

\[
\frac{|\delta(d/z^3)|}{d/z^3}
\lesssim \frac{|\delta d|}{d}+3\frac{|\delta z|}{z},
\qquad
\frac{|\delta C|}{C}\approx\left|\frac{C'(\tau)}{C(\tau)}\right||\delta\tau|.
\tag{A4}
\]

The gap-error columns are standalone measurement budgets. Finite-signal bias, uncertainty in the supporting functional/contact locations, spectral band shape, stray light, and reference calibration still have to be added. They are not claims of 10% total coefficient accuracy.

The JSON records also give a conservative per-channel budget. In unit-white sensor coordinates \(y\), the calculation uses \(x=My\); thus \(\delta d=n^TM\delta y\). If every channel has absolute error at most \(\epsilon\), then \(|\delta d|\leq\|n^TM\|_1\epsilon\), before uncertainty in \(n\) or \(M\). This provides a concrete specification to compare with measured apparatus performance.

Small normalized \(z\) does not fix absolute photon flux: common illumination scaling cancels in the normalized quantities. Increasing flux could reduce counting noise, but it does not remove fractional calibration or spectral leakage errors. No instrument performance has been measured or assumed here.

## 6. Decision supported by this work

The mathematical claim can proceed to external review with its hypotheses intact. The endpoint qualification now has an explicit counterexample showing why it matters. The empirical direction has a concrete design and quantitative precision requirements, but experimental feasibility remains unestablished.

The next experimental decision requires actual spectral-response, filter-resolution and uncertainty data. The most informative comparison would use mixtures on both sides of (A3), several independently measured signal levels, and pre-specified finite-signal error tolerances. The remaining literature access gap and the novelty question remain separate from those measurement requirements.

## Reproducibility

The existing calculation archive now includes this note, `proof_stress_audit.py`, its results, `finite_window_audit.py`, the three model result files, `finite_window_summary.py`, the selected-design summary, and the screening figure in PNG and PDF. The earlier result files and conditional theorem are retained; the prior-art note updates the West–Brill access status.

Install `requirements.txt`, then run the polynomial benchmark. For each model (`pchip`, `govardovskii`, `govardovskii_tilted`), run the finite-window script once normally, once with `--refine`, then once with `--design`. Finally run the summary script:

```bash
python proof_stress_audit.py
python finite_window_audit.py --model govardovskii
python finite_window_audit.py --model govardovskii --refine
python finite_window_audit.py --model govardovskii --design
python finite_window_summary.py
```

The normal finite-window run replaces only its own audit JSON; follow-up modes append their missing cases. The earlier spectral result files are not overwritten. The summary requires the audit JSON for all three models. Error records are deliberately retained.