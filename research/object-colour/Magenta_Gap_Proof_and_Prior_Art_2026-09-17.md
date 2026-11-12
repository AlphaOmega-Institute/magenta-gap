# The Magenta Gap: proof and prior-art assessment

Research note — 17 September 2026; updated after the literature follow-up and proof audit

## Finding

The cubic mixing relation has a complete proof under explicit assumptions on a continuous spectral model. The proof covers every fixed interior mixture. It also covers the endpoints when either of two sufficient conditions below holds. The endpoint conditions are satisfied numerically by the models already tested, including the stationary long-wave endpoint of the PCHIP model.

The literature check does **not** establish novelty. The full-text follow-up strengthens the prior-art finding: an exact circle example in Koshevoy and Mosler (1997) already implies a cubic deficit in the same bounded-density geometric framework. Their general mixture construction is also established. Burns (2021) supplies directly relevant boundary optimization and long-wave reversal precedents. The particular two-contact coefficient in this note was not located in the works read, but important full texts remain unchecked. The defensible description is **a conditional asymptotic relation, derived within established geometry and checked in spectral models**. It is not yet a new empirical law.

This note concerns the physical attainability of passive object colors. It supplies no evidence for a Lorentzian perceptual metric or for a particular topology of phenomenal hue.

## 1. Definitions and theorem

Let \(J=[a,b]\) be a compact wavelength interval, and let \(c:J\to\mathbb R^3\) be a continuous response density. A reflectance is any measurable function \(0\leq\rho\leq1\). Define

\[
\mathcal O=\left\{\int_J c(\lambda)\rho(\lambda)\,d\lambda:
0\leq\rho\leq1\right\},\qquad
W=\int_J c(\lambda)\,d\lambda.
\]

Choose a linear functional \(\ell\) such that

\[
b(\lambda):=\ell(c(\lambda))>0\quad(\lambda\in J),
\qquad \ell(W)=1.
\]

The scalar \(z=\ell(x)\) is the signal level relative to white. It is not assumed to equal perceptual brightness.

Suppose a linear functional \(n\), normalized by \(n(W)=1\), satisfies:

1. \(f(\lambda):=n(c(\lambda))\geq0\) throughout \(J\).
2. Its only zeros are distinct interior points \(\lambda_0,\lambda_1\).
3. The density \(c\) is \(C^2\) in neighborhoods of these points, and
   \(q_i:=f''(\lambda_i)>0\). Since the zeros are interior minima, \(f'(\lambda_i)=0\).
4. The normalized contact vectors \(p_i=c(\lambda_i)/b_i\), where \(b_i=b(\lambda_i)\), are distinct.

These assumptions specify an exposed cone face with two isolated quadratic contacts. Smoothness at unrelated interpolation knots is unnecessary.

For fixed \(\tau\in[0,1]\), set

\[
p_\tau=(1-\tau)p_0+\tau p_1.
\]

The limiting cone mixture uses signal shares \(1-\tau\) and \(\tau\). This is neither a wavelength fraction nor a perceptually uniform hue coordinate.

At signal \(z\), define the largest attainable radial parameter

\[
r_\tau(z)=\max\{r:zW+r(p_\tau-W)\in\mathcal O\},
\qquad d_\tau(z)=z-r_\tau(z).
\]

The cone permits \(r=z\). The actual object-color solid generally falls inside that limit. Because \(n(p_\tau)=0\) and \(n(W)=1\),

\[
d_\tau(z)=n\!\left(zW+r_\tau(z)(p_\tau-W)\right).
\]

Thus the deficit is affine and requires no choice of Euclidean chromatic distance. If a radius \(R_\tau=\alpha_\tau r_\tau\) is introduced, then \(d_\tau=z-R_\tau/\alpha_\tau\).

**Interior theorem.** Under assumptions 1–4, for every fixed \(0<\tau<1\),

\[
\boxed{\displaystyle
d_\tau(z)=
\left[A_0(1-\tau)^3+A_1\tau^3\right]z^3+o(z^3),
\qquad A_i=\frac{q_i}{24b_i^3},
\qquad z\downarrow0.}
\tag{1}
\]

No nonzero projective tangent is required for this interior statement. Section 3 supplies the additional endpoint conditions. All limits in this note are pointwise in fixed \(\tau\); a uniform joint limit with \(\tau\to0\) or \(1\) is not asserted.

## 2. Proof of the interior theorem

### 2.1 Existence and the exact constraints

The set \(0\leq\rho\leq1\) is compact in the weak-* topology of \(L^\infty(J)\). Integration against each component of \(c\in L^1(J)\) is continuous in this topology. Consequently \(\mathcal O\) is compact and convex, and the constrained optimum exists.

The vectors \(W,p_0,p_1\) are linearly independent: \(p_0,p_1\) are independent because they are distinct and both have \(\ell\)-value one; \(W\) lies outside their span because its \(n\)-value is one.

Choose a nonzero linear functional \(h_\tau\) whose kernel is \(\operatorname{span}\{W,p_\tau\}\). The two equations

\[
\ell(x)=z,\qquad h_\tau(x)=0
\tag{2}
\]

are exactly the condition \(x=zW+r(p_\tau-W)\) for some \(r\). Minimizing \(n(x)\) under (2) therefore gives \(d_\tau(z)\). The constant reflectance \(\rho=z\), for \(0\leq z\leq1\), is feasible and gives \(r=0\). Nonnegativity of \(f\) gives \(d_\tau(z)\geq0\).

### 2.2 Use signal mass as the local coordinate

Put \(d\mu=b(\lambda)\,d\lambda\) and \(p(\lambda)=c(\lambda)/b(\lambda)\). Then

\[
\mu(J)=1,\qquad x=\int_J p\rho\,d\mu,\qquad
z=\int_J\rho\,d\mu.
\]

In a neighborhood of \(\lambda_i\), use the signed mass coordinate

\[
t_i(\lambda)=\int_{\lambda_i}^{\lambda}b(u)\,du.
\]

Since \(b_i>0\), this is a smooth local coordinate. Write \(P_i(t)=p(\lambda(t))\). The supporting cost in these coordinates is

\[
g_i(t):=n(P_i(t))
=\frac{q_i}{2b_i^3}t^2+o(t^2)
=12A_i t^2+o(t^2).
\tag{3}
\]

In particular, a band of mass \(m\) centered at \(t=0\) has response

\[
H_i(m):=\int_{-m/2}^{m/2}P_i(t)\,dt
=mp_i+\frac{m^3}{24}\ddot P_i(0)+o(m^3),
\tag{4}
\]

and supporting cost

\[
n(H_i(m))=A_i m^3+o(m^3).
\tag{5}
\]

Dots indicate derivatives in the mass coordinate, not wavelength.

### 2.3 Upper bound: satisfy both moments exactly

Consider two disjoint centered bands with masses \(m_0,m_1\). Require

\[
m_0+m_1=z,\qquad
h_\tau(H_0(m_0)+H_1(m_1))=0.
\tag{6}
\]

Extend \(H_i\) to signed \(m\) by the oriented integral in (4), solely to apply the implicit function theorem at the origin. The Jacobian with respect to \((m_0,m_1)\) there is

\[
\begin{pmatrix}
1&1\\h_\tau(p_0)&h_\tau(p_1)
\end{pmatrix}.
\]

It is invertible. Otherwise \(h_\tau(p_0)=h_\tau(p_1)\); together with \(h_\tau(p_\tau)=0\), this would make \(h_\tau\) vanish on the basis \(W,p_0,p_1\).

Equation (4) and the implicit function theorem give an exact solution of (6) with

\[
m_0=(1-\tau)z+O(z^3),\qquad
m_1=\tau z+O(z^3).
\tag{7}
\]

For fixed interior \(\tau\), both masses are positive for sufficiently small \(z\); the signed extension is no longer used. The corresponding indicator reflectance is physically admissible. Equations (5) and (7) yield

\[
\limsup_{z\downarrow0}\frac{d_\tau(z)}{z^3}
\leq A_0(1-\tau)^3+A_1\tau^3.
\tag{8}
\]

The construction need not be the exact optimizer at finite \(z\). A feasible construction suffices for this bound.

### 2.4 A scalar bound

If \(0\leq\eta\leq1\) on \(\mathbb R\) and \(\int\eta=m\), then

\[
\int_{\mathbb R}t^2\eta(t)\,dt\geq\frac{m^3}{12}.
\tag{9}
\]

To see this, set \(a=m/2\) and \(I=[-a,a]\). Pointwise,
\((t^2-a^2)(\eta-\mathbf1_I)\geq0\). Integrating and using equal masses cancels the \(a^2\) term. The remaining integral over \(I\) is \(m^3/12\).

This elementary filling argument is the bounded-density rearrangement step. It does not assume in advance that an optimizer has two bands.

### 2.5 Lower bound: concentration fixes the leading masses

Let \(\rho_z\) be an optimal reflectance. The upper bound gives \(d_\tau(z)=O(z^3)\). Define a probability measure

\[
d\nu_z=\frac{\rho_z\,d\mu}{z}.
\]

Its mean supporting cost is

\[
\int n(p)\,d\nu_z=d_\tau(z)/z=O(z^2).
\]

The continuous nonnegative function \(n(p)\) vanishes only at the two contacts. Every weak limit of \(\nu_z\) is therefore supported on \(\{\lambda_0,\lambda_1\}\). Moreover,

\[
\int p\,d\nu_z
=\frac{x_\tau(z)}z
=p_\tau+\frac{d_\tau(z)}z(W-p_\tau)\longrightarrow p_\tau.
\]

Distinct \(p_0,p_1\) force the unique limiting measure

\[
\nu_z\rightharpoonup
(1-\tau)\delta_{\lambda_0}+\tau\delta_{\lambda_1}.
\tag{10}
\]

For fixed disjoint small contact neighborhoods \(U_i\), let
\(m_i(z)=\int_{U_i}\rho_z\,d\mu\). Equation (10) gives

\[
m_0(z)/z\to1-\tau,\qquad m_1(z)/z\to\tau.
\tag{11}
\]

For every \(\varepsilon>0\), the neighborhoods can be chosen so that (3) implies
\(g_i(t)\geq(12A_i-\varepsilon)t^2\) inside them. Extend each restricted density by zero to the real line and apply (9). The cost outside the neighborhoods is nonnegative. Hence

\[
d_\tau(z)\geq
\sum_{i=0}^1(A_i-\varepsilon/12)m_i(z)^3.
\]

Divide by \(z^3\), take the lower limit using (11), and then let \(\varepsilon\downarrow0\). This matches (8), proving (1).

## 3. Endpoints: the missing qualification made explicit

Quadratic contact alone supplies the local cost, but an upper-bound construction must also satisfy the hue constraint with nonnegative band masses. The earlier calculation note did not state this endpoint requirement precisely enough.

Fix \(i\in\{0,1\}\), let \(j=1-i\), and choose \(h_i\ne0\) with
\(h_i(W)=h_i(p_i)=0\). Necessarily \(h_i(p_j)\ne0\).

Either of the following conditions is sufficient for

\[
\boxed{d_i(z)=A_i z^3+o(z^3).}
\tag{12}
\]

### E1. Regular projective tangent

Assume \(\dot P_i(0)\ne0\). Since
\(\ell(\dot P_i(0))=n(\dot P_i(0))=0\), this tangent is parallel to \(p_1-p_0\). Thus \(h_i(\dot P_i(0))\ne0\).

Use one band of mass \(z\), centered at mass coordinate \(s\). Define

\[
F(z,s)=\frac1z\,h_i\left(
\int_{s-z/2}^{s+z/2}P_i(t)\,dt\right),
\qquad F(0,s)=h_i(P_i(s)).
\]

The averaging expression extends smoothly enough through \(z=0\). We have
\(F(0,0)=0\), \(F_s(0,0)\ne0\), and \(F(z,0)=O(z^2)\).
The implicit function theorem gives \(s(z)=O(z^2)\) with \(F(z,s(z))=0\). By (3), its cost is \(A_i z^3+o(z^3)\). The concentration and scalar lower-bound argument above, now with limiting mass one at contact \(i\), proves (12).

### E2. A positive correction using the other contact

This alternative does not require a nonzero tangent. Define

\[
\Gamma_i:=
-\frac{h_i(\ddot P_i(0))}{24h_i(p_j)}.
\tag{13}
\]

Assume \(\Gamma_i>0\). Use a centered main band of mass \(z-e\) at contact \(i\), and a centered correction band of mass \(e\) at contact \(j\). The hue equation is

\[
G(z,e):=h_i(H_i(z-e)+H_j(e))=0.
\]

Its derivative \(G_e(0,0)=h_i(p_j)\ne0\). Expansion (4) therefore gives an exact solution

\[
e(z)=\Gamma_i z^3+o(z^3)>0.
\tag{14}
\]

For sufficiently small \(z\), also \(e(z)<z\). Both bands are admissible. Their total cost is

\[
A_i(z-e)^3+o(z^3)+O(e^3)=A_i z^3+o(z^3).
\]

The other band's own supporting cost is of order \(z^9\). The same lower-bound argument proves (12). Condition E2 also works at a regular endpoint; it is particularly useful when E1 fails.

If neither E1 nor E2 holds, this note makes no endpoint coefficient claim. Failure of these sufficient conditions is not itself a proof that (12) is false.

### Endpoint audit for the tested models

In the computation coordinates, \(W=(0,0,1)\). We used
\(h_i(v)=p_{i,x}v_y-p_{i,y}v_x\).

| Model | \(\Gamma_0\) | \(\Gamma_1\) | Endpoint conclusion |
|---|---:|---:|---|
| Stockman–Sharpe / PCHIP | 0.6761819442 | 1496.8526984 | E2 at both; E1 fails at the long contact |
| Govardovskii / equal energy | 0.1315458991 | 6049.0121365 | E1 and E2 at both |
| Govardovskii / red tilt | 0.1404559182 | 1252.6962552 | E1 and E2 at both |

These are high-precision numerical evaluations of the represented models. They do not constitute a proof that a fitted spectral model exactly describes physiology, or an interval-arithmetic certification of all global contact assumptions.

For PCHIP, the long contact is at approximately 703.417797 nm. The normalized spectral curve is stationary there. The positive value of \(\Gamma_1\) resolves the endpoint construction without incorrectly assuming a nonzero tangent.

The accompanying script constructs the E2 reflectance exactly to numerical tolerance. The following are feasible upper-bound witnesses, not asserted finite-level optima:

| Signal \(z\) | Constructed deficit / \((A_1z^3)\) | Correction mass / \(z^3\) |
|---|---:|---:|
| \(10^{-5}\) | 1.006943227 | 1507.245686 |
| \(3\times10^{-6}\) | 1.000624227 | 1497.787074 |
| \(10^{-6}\) | 1.000069352 | 1496.956508 |

Signal and hue residuals are below \(10^{-50}\) in this calculation. The limiting correction coefficient is 1496.852698. The proof, rather than these three numerical samples, supplies the asymptotic lower bound.

## 4. The resulting mixing relation and its scope

If each endpoint satisfies E1 or E2, define
\(C(\tau)=\lim_{z\downarrow0}d_\tau(z)/z^3\). Equations (1) and (12) give

\[
\boxed{C(\tau)=C(0)(1-\tau)^3+C(1)\tau^3,\quad 0\leq\tau\leq1.}
\tag{15}
\]

In particular,

\[
C(1/2)=\frac{C(0)+C(1)}8,
\qquad
\lim_{z\downarrow0}
\frac{d_\tau(z)}
{(1-\tau)^3d_0(z)+\tau^3d_1(z)}
=1.
\tag{16}
\]

The second expression removes the need to fit a separate amplitude to each interior mixture. It still requires a specified signal functional, white, contact pair, and mixture coordinate. The previous numerical study predicted coefficients from spectral derivatives; it did not independently measure the two endpoint deficits.

The cubic exponent comes from a one-dimensional bounded-density problem with quadratic supporting cost. It is not unique to magenta. A regular isolated spectral contact has the same local mechanism. The two-contact face supplies the two weighted cubic contributions.

The result concerns a *continuum* of wavelengths. Replacing that continuum by a fixed set of monochromatic atoms changes the asymptotic problem: small enough masses placed exactly at two contact atoms can remain on the cone face, giving zero gap until their capacities bind. A finite spectral grid is not, by itself, evidence for the continuum limit.

The theorem also excludes contacts at imposed spectral cutoffs, higher-order or flat contacts, extra supporting zeros, and changes to the reflectance bound. The finite-level optimizer need not always consist of two bands; the lower bound deliberately ranges over all admissible reflectances.

## 5. Exactly what is scale independent

An invertible linear coordinate change \(T\) sends \(c,W,p\) to \(Tc,TW,Tp\), and the covectors to \(\ell T^{-1},nT^{-1}\). It leaves \(z,d,\tau,A_i\) unchanged. This includes consistent changes in individual sensor gains.

A common illuminant-intensity change \(c\mapsto kc\) is canceled by white normalization: use \(\ell\mapsto\ell/k\) and \(n\mapsto n/k\). Again \(z,d,A_i\) are unchanged.

For a smooth change of wavelength variable \(\lambda=\phi(u)\), the response *density* must transform as
\(\widetilde c(u)=c(\phi(u))|\phi'(u)|\). At a quadratic contact,
\(\widetilde q_i=q_i|\phi'(u_i)|^3\) and
\(\widetilde b_i=b_i|\phi'(u_i)|\), so
\(\widetilde q_i/(24\widetilde b_i^3)=A_i\).

There is no discrimination metric in the construction, so its arbitrary overall scale does not enter.

Changing the illuminant's spectral shape, observer functions, or the signal functional generally changes coefficients and the mixture coordinate. An arbitrary nonlinear perceptual transformation need not preserve the relation. “Scale independent” should refer to the specified unit and coordinate changes, not to invariance under every change of viewing condition.

## 6. Prior-art assessment

### What was checked

The initial search followed color-science references and the general convex construction. The follow-up obtained and read the complete **published** Burns article and the complete Koshevoy–Mosler paper, closing two previously listed full-text gaps. It also checked both versions of Burns's author-hosted supplementary discussion and Davis's current mathematical vignette. These were compared with the definitions, fixed-mixture limit, and endpoint construction in this note. Failure to find the exact coefficient is not a novelty certificate.

| Primary source | Material actually inspected | Relevance and consequence |
|---|---|---|
| Logvinenko, *The color cone* (2015), [DOI](https://doi.org/10.1364/JOSAA.32.000314) | Publisher abstract and [university record](https://researchonline.gcu.ac.uk/en/publications/the-color-cone/); full text not obtained | The cone and effective visible spectrum are prior work. The university record supplies a publisher link, not an accessible manuscript. |
| Logvinenko and Levin, *Foundations of Colour Science*, chapters 14–16, [publisher contents](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119885955) | Contents, chapter summaries and metadata; chapter full texts not obtained | The book comparison remains open and must cover more than chapter 14. The cited p. 383 lies in chapter 16. Online publication is 2022; the 2025 article cites the book as 2023. |
| Logvinenko, Funt and Bastani, *The object-color solid* (2025), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11801397/) | Full article available; support representation, pole discussion, interpolation caveats, and Appendix A checked closely | Establishes the immediate framework: optimal reflectances, cone-like poles, two-contact geometry, and possible higher transition counts. The displayed cubic mixing coefficient was not located here. |
| Cox, *On the Mathematical Foundations of Colorimetry* (2021 preprint), [record](https://doi.org/10.31234/osf.io/afdr3) | Full PDF obtained; object-color and optimal-color sections, especially pp. 11–21 and equation (15), read | A measure-theoretic formulation and sign-based support rule predate this analysis. Its claim of at most two transitions is not adopted; later work shows why that restriction is unsafe. |
| West and Brill, *Conditions under which Schrödinger object colors are optimal* (1983), [DOI](https://doi.org/10.1364/JOSA.73.001223) | Public OCR text and captions of all three pages read; original PDF unavailable | Convexity, wavelength ordering, barycentric proof, and general support-line optimization are prior work. No local cubic mixing relation was located. |
| Burns, *The location of optimal object colors with more than two transitions* (2021), [published PDF](http://scottburns.us/wp-content/uploads/2025/03/The-location-of-optimal-object-colors-with-more-than-two-transitions-Burns-2021.pdf), [DOI](https://doi.org/10.1002/col.22693) | Complete published article, *Color Research & Application* 46, 1180–1193; preprint v7 and author supplements also checked | Boundary optimization, extra transitions and the modern long-wave reversal are prior work. Details below. |
| Davis, *Convexity and Transitions*, [current PDF](https://cran.r-project.org/web/packages/colorSpec/vignettes/convexity.pdf) | Complete eight-page vignette dated 8 June 2025, especially §§2–6 and Theorem 6.1 | A discrete formulation of the transition issue, with worked counterexamples. This is the currently available revision, not a verified copy of the 2020 version cited by Burns. It supplies no continuum two-contact asymptotic coefficient. |
| Huang and Slomka, *Approximations of convex bodies by measure-generated sets* (2017 preprint), [full text](https://arxiv.org/abs/1706.07112) | Definition 1.1 and Proposition 2.1 with its proof | The bounded-density mean construction and support selection are known as metronoids. They are not a new geometric object introduced by this project. |
| Huang, Slomka and Werner, *Ulam floating bodies* (2018 preprint), [full text](https://arxiv.org/html/1803.08224v2) | Definitions, main theorem statements, and relevant support construction | Related asymptotic convex geometry exists. Its principal volume-density assumptions differ from a measure supported on a spectral curve; its results cannot simply be substituted for (1). |
| Koshevoy and Mosler, *Zonoid trimming for multivariate distributions* (1997), [author-hosted published PDF](https://wisostat.uni-koeln.de/sites/dep_econometrics/pdf_publikationen/1997_KoshevoyM97.pdf), [DOI](https://doi.org/10.1214/aos/1069362382) | Complete paper, *Annals of Statistics* 25, 1998–2017; key equations on printed pp. 2008–2009 also checked against page images | Direct prior framework and an exact example already implying a cubic deficit. Details below. |

### The exact connection to established convex geometry

Our inference from those definitions can be written directly. Push the probability measure \(d\mu=b\,d\lambda\) forward by the normalized curve \(p=c/b\), and call the resulting measure \(\nu\). The attainable normalized colors at level \(z\) are

\[
\left\{\frac{x}{z}:x\in\mathcal O,\ \ell(x)=z\right\}
=
\left\{\int y\,u(y)\,d\nu(y):
0\leq u\leq1/z,\ \int u\,d\nu=1\right\}.
\tag{17}
\]

If \(p\) is not injective, conditional expectation of \(\rho\) given \(p\) gives the right-hand density; conversely every right-hand density pulls back to an admissible reflectance. Thus (17) is an equality, not just an analogy. It is a zonoid trimmed region, or equivalently the metronoid of \(\nu/z\).

This is precisely the structure of Koshevoy and Mosler's Definition 1.1 and Proposition 2.2. Equation (1) is therefore a local two-contact asymptotic calculation within a known class. Its proof uses Taylor expansion, a bounded-density filling inequality, concentration, and the implicit function theorem.

### The strongest new prior-art finding: the 1997 circle example

Koshevoy and Mosler give the following exact example in §6, example 2, printed p. 2009: for the uniform probability measure on a unit circle, the trimmed region is a disk of radius

\[
R(z)=\frac{\sin(\pi z)}{\pi z}.
\tag{18}
\]

**Our deduction from their formula:** in the lifted coordinates used here, the cone slice has radius \(z\), while the attainable slice has radius \(zR(z)\). Thus its radial deficit is

\[
d(z)=z-zR(z)=z-\frac{\sin(\pi z)}\pi
=\frac{\pi^2}{6}z^3+O(z^5).
\tag{19}
\]

The expansion in (19) is ours; the exact radius in (18) is theirs. The example has a single supporting contact, so it does not state our two-contact formula. It does establish a direct mathematical precedent for cubic scaling in the general construction.

Their Theorem 5.8, equation (5.6), printed p. 2008, gives the exact mixture identity. Written with two probability measures \(\mu_0,\mu_1\), it is

\[
D_z\bigl(\beta\mu_0+(1-\beta)\mu_1\bigr)
=\bigcup_{\substack{0\leq\delta_0,\delta_1\leq1\\
\beta\delta_0+(1-\beta)\delta_1=z}}
\left[
\frac{\beta\delta_0}{z}D_{\delta_0}(\mu_0)
+\frac{(1-\beta)\delta_1}{z}D_{\delta_1}(\mu_1)
\right].
\tag{20}
\]

Here \(\beta\) is a weight in the underlying measure, not automatically the selected mixture coordinate \(\tau\). This identity organizes the allocation of mass between components; it does not by itself prove our fixed-ray feasibility or endpoint correction. Neither the abstract nor a keyword search exposed this comparison as clearly as reading the paper did.

The paper's Theorem 3.1 assumes a density with respect to ambient Lebesgue measure. Our spectral measure lies on a curve, so that theorem cannot be imported without checking its hypothesis. The general definitions and the circle example do cover singular measures. The paper's empirical convergence results concern increasing sample size, a different limit from \(z\downarrow0\).

### What the complete Burns paper changes

The published version matters: it includes modern color-matching-function analysis absent from [preprint v7](https://arxiv.org/abs/2103.06997). Section 3, equation (3), optimizes bounded reflectances along a ray; §11 also discusses slices through a point on the gray axis. These methods require attribution.

Sections 5–6 compare unrestricted optima with two-transition candidates. That finite-level comparison is not our deficit from the cone face. Section 10, pp. 1191–1192, reports a reversal of the modern spectral locus near **703 nm**. Our stationary PCHIP contact near 703.417797 nm should therefore not be presented as discovery of that qualitative feature. The precise interpolation-dependent value and the E2 feasibility calculation are separate results.

No counterpart of (15), or a near-black cubic derivation, was located in the complete article or the supplementary discussions inspected. The [current journal supplement](http://scottburns.us/ocss_suppl_docs/) is dated 12 March 2025 and explicitly distinguishes its updateable author commentary from peer-reviewed material.

### West–Brill follow-up

The [public transcription](https://www.academia.edu/54927465/Conditions_under_which_Schr%C3%B6dinger_object_colors_are_optimal) closes the unread-text gap. Its OCR interleaves columns and garbles some notation; the figures were not visually verified. The original PDF remains desirable for those checks. Later Burns and Davis analyses qualify its claim about the human locus.

### Remaining full-text gaps

| Work still needed | Exact material | Access result |
|---|---|---|
| Logvinenko (2015), *The color cone* | *JOSA A* 32(2), pp. 314–322 | Publisher abstract and institutional metadata accessible; original full text unavailable in this session. |
| Logvinenko and Levin, chapter 14 | *Object-Colour Solid*, pp. 357–364; [chapter DOI](https://doi.org/10.1002/9781119885955.ch14) | Summary and contents checked; full chapter unavailable. |
| Same book, chapter 15 | *Trichromatic Regular Object-Colour Solid*, pp. 365–374; [chapter DOI](https://doi.org/10.1002/9781119885955.ch15) | Summary and contents checked; full chapter unavailable. |
| Same book, chapter 16 | *Object-Colour Stimulus Manifold*, pp. 375–388; [chapter DOI](https://doi.org/10.1002/9781119885955.ch16) | Summary and contents checked; full chapter unavailable. The 2025 article's cited p. 383 falls here. |

The earlier reading list named only chapter 14. The publisher's pagination shows that chapters 14–16 are needed to check the cited object-color boundary discussion. Checking only chapter 14 would leave that comparison incomplete. Public catalog, repository, publisher and preview checks did not yield these full texts; none is counted as read. Negative book-search snippets cannot rule out equations or differently named results. Copies supplied through the user's existing access would allow this remaining comparison.

### Novelty verdict

**Established background:** the physical solid and cone, supporting-plane optimization, complementary symmetry, bounded-density convex representation, allocation between component measures, and the need to permit extra transitions. Cubic behavior already follows from an exact example in the 1997 literature.

**Potential contribution, not verified as new:** the explicit fixed-mixture coefficient \(A_0(1-\tau)^3+A_1\tau^3\), its endpoint-normalized version, the endpoint feasibility treatment, and the continuum numerical validation in this colorimetric setting. These may amount to a specialized application or extension of known mathematics. Absence of the displayed formula in the texts read does not establish a new theorem's priority.

**Not supported:** a claim that the general geometry or cubic scaling mechanism is new; a claim that exponent three uniquely distinguishes magenta; a claim that the relation is already empirically established; or a claim that it validates the original Lorentzian construction.

The highest-priority remaining comparison is the book's chapters 14–16, followed by the complete 2015 cone paper. The 1997 zonoid and 2021 Burns full-text gaps are closed; West–Brill's text has now been read with the OCR limitations stated above. This remains a bounded review, not an exhaustive citation search or an expert determination of priority.

## 7. Implication for the paper

“The Magenta Gap” remains a suitable project title if the paper's central question is made precise: how a passive object-color boundary falls inside the nonspectral cone face, and how that deficit varies across mixtures.

The current result can support a mathematical section with explicit assumptions and a model-validation section. An empirical-law claim additionally needs independently calibrated spectral responses and boundary measurements, with uncertainty small enough to resolve the deficit. Predictions computed from the same spectral model being optimized are a consistency test, not independent experimental evidence.

The immediate research claim should be the conditional relation in (15), accompanied by its scope and the prior-art qualifications above. The historical title does not remove the need to revise unsupported claims about perception or metric geometry.

Following the full-text check, the introduction should credit the known convex construction and the earlier cubic example before presenting the two-contact calculation. Its scientific value must be assessed through the precision, applicability and validation of that calculation; novelty of the exponent is not a defensible premise.

## Reproducibility and revision record

The companion calculation package now includes this proof note, the clarified original calculation note, and endpoint_audit.py with endpoint_audit_results.json. Run the endpoint_audit.py script after installing the package requirements. It recomputes the endpoint conditions for all three model/illuminant combinations and the three PCHIP endpoint witnesses.

For reproducibility, wavelength derivatives are converted using

\[
\dot P=\frac{p'}b,\qquad
\ddot P=\frac{p''}{b^2}-\frac{p'b'}{b^3},
\]

\[
p'=\frac{c'}b-\frac{cb'}{b^2},\qquad
p''=\frac{c''}b-\frac{2c'b'}{b^2}
-\frac{cb''}{b^2}+\frac{2c(b')^2}{b^3}.
\]

The earlier interior numerical results are unchanged. The endpoint statement has been clarified: a concentration lower bound by itself does not establish an endpoint upper bound. Conditions E1 and E2 and their constructions provide the missing justification. These changes are confined to the separate research notes and calculation package.

The subsequent literature revision updates the finding and §§6–7 after full-text inspection of Burns (2021), Koshevoy and Mosler (1997), and Davis's current vignette. Equations (18)–(20) record the new prior-art comparison. It also corrects the remaining book scope from chapter 14 to chapters 14–16. The theorem, its proof, scripts and numerical results were not changed by this literature revision. Source links, editions and printed page references above identify what was inspected; downloaded publisher papers are not redistributed in the calculation package.

The latest follow-up updates the West–Brill access status and adds `Magenta_Gap_Proof_Audit_and_Test_Design_2026-09-17.md` to the package. That companion note records a further audit by the same assistant, a separate dual lower bound, an explicit counterexample to dropping the endpoint qualifications, and new finite-signal calculations with measurement budgets. No fatal gap was found in the qualified theorem; this is not an external referee report or empirical confirmation. The theorem and the earlier numerical result files remain unchanged.