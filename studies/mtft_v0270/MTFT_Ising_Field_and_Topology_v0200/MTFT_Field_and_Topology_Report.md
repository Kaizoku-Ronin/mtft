# MTFT Ising companion v0.2.0: fields, topology, and exact sampling

Built from the supplied **MTFT v0.26.0** source, 5 September 2026.

The first instrument counted the zero-field Ising states of the dual Manin graph. This extension resolves their magnetization, measures selected surface twists, reverses every antiferromagnetic bond in turn, couples two small copies, finds finite-polynomial zeros, and draws equilibrium configurations with exact integer conditional weights.

Three results stand out. One reversed bond can force at least **nine spins** to change. A small field changes the sign of the finite ensemble's **Fisher curvature**, with three independent numerical assemblies agreeing at high precision. And at three ferromagnetic coupling ratios, exact polynomial certificates place **all 56 Lee–Yang zeros on the unit circle**.

| Instrument | Completed scope | Evidence class |
|---|---|---|
| Joint energy and magnetization | N=6, 11, 15, 35, 143; all 2^56 states at N=143 | Exact integer tables; numerical observables |
| Surface twists | All 26 basis twists at N=143, plus lower-level controls; two retained homology coordinates | Exact counts, flatness/rank/gauge checks, independent cycle census |
| Single-bond defects | All 84 AF-to-ferro reversals at N=143 | Exact ground counts and closest witnesses |
| Two coupled fields | N=6, 11, 15, 35; largest case is two 16-spin layers | Exact joint tables and limiting identities |
| Complex zeros | Temperature polynomial at zero field; field polynomials for edge-weight ratios 2, 4, 8 | Exact factorization and field-circle certificates; numerical temperature roots |
| Equilibrium sampler | 20,000 independent draws at each of four rational weight ratios | Exact integer conditional weights; sampling agreement is diagnostic |

The delivered `verification.json` records 165 passing checks: 142 direct exact checks, ten recorded exact construction checks, seven numerical diagnostics, and six command-line integration checks. “Exact” here identifies finite integer or rational computations. It does not upgrade numerical curvature or complex-root residuals to interval bounds.

## The common finite model

At N=143, the dual graph has n=56 spin vertices, E=84 edges, one self-loop, and a surface embedding of genus 13. Spins sit on Manin triangles. The graph's drawing positions carry no physical metric.

Use J=kB=1 and natural parameters β and η=βh. For ordinary ferromagnetic bonds,

\[
S(\sigma)=\sum_{e=(u,v)}\sigma_u\sigma_v=E-2k,
\qquad M(\sigma)=\sum_v\sigma_v=2m-n,
\]

where k counts disagreeing edges and m counts up spins. A self-loop contributes +1 to S and never contributes a disagreement. The model is

\[
p_{\beta,\eta}(\sigma)=Z^{-1}\exp(\beta S+\eta M),\qquad
Z=\sum_{k,m}D_{k,m}\exp\{\beta(E-2k)+\eta(2m-n)\}.
\]

Thus every real-field observable below is a finite sum from the same counting table. The exactly zero mean magnetization at η=0 follows from global spin reversal; absolute magnetization and susceptibility expose the finite system's ordering behavior.

## 1. Joint counting and Fisher geometry

`joint_143.json` contains the complete 85×57 array D[k,m], including zero entries. Its entries sum to

\[
2^{56}=72,057,594,037,927,936.
\]

The extension preserves the previous contraction strategy. A min-fill elimination order has width at most 10 and uses at most 2,048 spin assignments in its largest intermediate table. This is an upper bound from one order, not a proof of optimal treewidth. Encoding two polynomial variables uses the integer marker B=2^57 and exponent k+85m. Nonnegative partial counts stay below B, preventing carries between coefficients. Larger packed integers explain why joint counting costs more than the earlier univariate calculation. The relation between contraction cost and graph width is established background; see [Markov and Shi, *Simulating Quantum Computation by Contracting Tensor Networks*](https://arxiv.org/abs/quant-ph/0511069).

The recorded N=143 joint contraction kernel took **5.223 seconds**. This excludes graph construction, unpacking, and verification. Timings describe this run, not a portable performance guarantee.

Checks included:

- Summing over m reproduces the earlier energy density of states coefficient by coefficient.
- Summing over k gives binomial(n,m), and D[k,m]=D[k,n−m].
- Cubic incidence gives the parity gate k≡m mod 2, including the self-loop's degree contribution.
- Full joint enumeration agrees on the smaller graphs through 16 spins.
- A different elimination order gives the identical N=143 integer table; direct enumeration checks its m=0,1,2,3 columns.
- Three evaluations with integer edge and field weights agree with separate scalar integer contractions.

For θ=(β,η), the Fisher metric is the full, unscaled covariance matrix

\[
g_{ab}=\partial_a\partial_b\log Z=\operatorname{Cov}(T_a,T_b),
\qquad T=(S,M).
\]

The covariance identity is a standard exponential-family result; see [Geyer's exponential-family notes](https://www.stat.umn.edu/geyer/5421/notes/expfam.html). We report Gaussian curvature K; scalar curvature is R=2K. Dividing the metric by n would multiply K by n, so the normalization matters.

For ψ=log Z, the Hessian-metric formula used is

\[
K=-\frac{1}{4(\det g)^2}
\det\begin{pmatrix}
\psi_{11}&\psi_{12}&\psi_{22}\\
\psi_{111}&\psi_{112}&\psi_{122}\\
\psi_{112}&\psi_{122}&\psi_{222}
\end{pmatrix}.
\]

Third derivatives are third cumulants of (S,M). `analyze_geometry.py` also assembles MTFT's Brioschi formula and a Christoffel–Riemann route using third and fourth cumulants. A three-outcome categorical model gives K=+1/4 and fixes the sign convention.

The three assemblies share the counting table and moment/cumulant inputs. Their agreement checks the curvature formulas and sign convention; it is not three independent derivations of the underlying probability distribution. The double-precision grid separately evaluates the moments with NumPy.

| β | η | Mean M/n | Susceptibility β Var(M)/n | Gaussian K |
|---:|---:|---:|---:|---:|
| 0.30 | 0 | 0 exactly | 0.89704 | +0.0669672 |
| 0.64 | 0 | 0 exactly | 20.07547 | +0.0871757 |
| 0.64 | 0.03 | 0.669864 | 6.89138 | −0.0553079 |
| 1.00 | 0.08 | 0.987432 | 0.07229 | −0.7162337 |

These are numerical finite-ensemble values. At all four points, the three curvature routes agree to absolute differences below 1.2×10^−59 with 60-digit working precision. This measures computational agreement, not a rigorous rounding enclosure. Double-precision curvature differs by less than 3.1×10^−13 at these points.

The grid covers β∈[0.05,1], η∈[−0.08,0.08], with 96×81 points. Every metric determinant is positive. Magnetization is odd and curvature even under η→−η within the recorded numerical tolerances. The field-induced sign change is a property of this family of probability distributions. It does not identify curvature of X0(143) or of spacetime. Finite-system ordering and peaks alone establish no continuum critical point.

Figure: `MTFT_Field_Geometry.png`.

## 2. Surface twists and two homology coordinates

The ordinary zero-field graph partition function omits the embedding. The earlier control becomes more useful here: N=6 and N=11 have the same labeled edge multiset and identical joint spin counts, but genera 0 and 1. The sphere admits no nontrivial flat sign-twist class; the torus admits two independent coordinates.

A bond preference b_e∈{0,1} assigns frustration

\[
k_b(\sigma)=\sum_{e=(u,v)}[(s_u\mathbin{\mathrm{xor}}s_v)\mathbin{\mathrm{xor}}b_e],
\]

where spin bits s=0,1 represent −1,+1. The Hamiltonian is H_b=2k_b−E. A spin gauge change shifts b by a dual coboundary and preserves the zero-field partition function.

The MTFT tree/cotree basis supplies 2g flat dual-edge sign cocycles via primal/dual identification. The calculation checks zero flux around dual cusp faces and independence modulo spin-gauge changes. At N=143 all 26 representatives were evaluated, each with a randomly chosen gauge control producing the identical integer density of states. Every nontrivial representative has positive minimum frustration.

| Minimum frustrated bonds | Number of the 26 chosen basis representatives |
|---:|---:|
| 1 | 16 |
| 2 | 9 |
| 3 | 1 |

These are statistics of the recorded lexicographic basis. A different basis can change this table. They are not intrinsic basis-free invariants or Atkin–Lehner eigen-sectors.

Numerical free-energy costs use

\[
\beta\Delta F_b=\log Z_0-\log Z_b.
\]

At β=0.64, basis representatives 0, 6, and 25 have ΔF/J≈0.87261, 1.03211, and 1.38377. At low temperature their energy penalties tend to twice their minimum frustration counts.

There is an independent exact route. The high-temperature expansion is a sum over even subgraphs, with a flat sign twist contributing a parity character. Retaining the first two cocycle pairings and performing a four-term Fourier transform produces four nonnegative histograms A_h[j], indexed by occupied-edge count j and two parity bits h. The relationship between surface Ising expansions and homology is discussed in [Chelkak, Cimasoni, and Kassel, *Revisiting the Combinatorics of the 2D Ising Model*](https://unige.ch/math/folks/cimasoni/Ising-Combinatorics.pdf).

The C++ program independently enumerated **all 2^29=536,870,912 even subgraphs**, using a spanning-tree cycle basis and Gray-code updates. Every one of the **4×85=340 coefficients** agrees with the transform of the spin counts. Each parity class contains 2^27 subgraphs when every subgraph has equal weight. The recorded enumeration took 4.002 seconds.

The plotted class probability uses weights tanh(β)^j. Only two of the 26 homology coordinates are retained; this is a four-class projection, not a computation of all 2^26 sectors. Flat bond-sign characters are linear parity data; they should not be confused with the quadratic spin-structure weights in a surface Pfaffian formula. For that distinct construction, see [Cimasoni and Reshetikhin, *Dimers on Surface Graphs and Spin Structures*](https://arxiv.org/abs/math-ph/0608070).

Figure: the first two panels of `MTFT_Twists_and_Layers.png`.

## 3. Every one-bond antiferromagnetic defect

The baseline antiferromagnet prefers disagreement on all 84 edges. Its ground states have nine frustrated edges and energy −66. There are exactly two ground configurations, related by global spin reversal.

For each edge, reverse that one preference from antiferromagnetic to ferromagnetic, holding every bond magnitude at one. Compute the entire new density of states, its ground degeneracy, and the closest ground configuration to a fixed baseline coloring. “Closest” minimizes Hamming distance over all new ground states and therefore includes the global spin-flipped choices.

| Minimum required spin changes | Number of bond reversals |
|---:|---:|
| 0 | 37 |
| 1 | 33 |
| 3 | 8 |
| 5 | 2 |
| 9 | 4 |

The new minimum frustration count is eight in nine cases, nine in 47 cases, and ten in 28 cases. Ground-state degeneracies range from two to 60.

Edges **4, 22, 34, and 48** each force at least nine spins to change. For edge 4, joining vertices 4 and 3 in the source's zero-based indexing, the ground energy stays −66 and the ground degeneracy stays two. The two new ground states lie at Hamming distances 9 and 47 from the reference.

This last statement has two exact computations behind it. Min-plus elimination lexicographically minimizes frustration and distance. Separately, gauging the reference coloring to all-zero spins turns up-spin count into Hamming distance. An exact joint counting contraction then gives one ground configuration at distance 9, one at 47, and zero at every other distance. Thus a closer state is ruled out by counting, not just by finding one optimizer. Min-plus contraction is related to the tropical tensor-network approach to spin ground states; see [Liu et al., *Tropical Tensor Network for Ground States of Spin Glasses*](https://arxiv.org/abs/2008.06888).

The self-loop is included in the 84-edge census. Reversing it changes a constant energy contribution and does not rearrange spins. The nine-spin examples concern ordinary bonds. A ground-state rearrangement gives an equilibrium comparison; no transition path or physical response time was computed.

Figure: `MTFT_One_Bond_Defect.png`. All witnesses and counts are in `defect_census.json`.

## 4. Two explicitly coupled copies

Put two spin fields σ and τ on the same graph and introduce the coupling κ through

\[
p(\sigma,\tau)\propto
\exp\left[\beta S_{\rm intra}+\kappa Q\right],
\qquad
Q=\sum_i\sigma_i\tau_i=n-2j.
\]

Here S_intra is the sum of the two within-layer interaction statistics, and j counts mismatched sites between layers. This is an explicit model extension: κ is a chosen coupling. It is not derived from an Atkin–Lehner action or identified with a product manifold.

`bilayer_counts` uses a four-state variable at each site and computes B[k,j], where k is the total disagreement count across both layers. Calculations completed at N=6,11,15,35. The largest case has **32 total spins**, with all 2^32 assignments counted; its recorded contraction kernel took 0.297 seconds.

Four coefficient-level identities hold:

\[
\sum_j B_{k,j}=(D*D)_k,
\qquad
\sum_k B_{k,j}=2^n\binom nj,
\qquad B_{k,j}=B_{k,n-j}.
\]

The j=0 column equals D[k/2] for even k and vanishes for odd k. Consequently,

\[
Z_{\rm double}(\beta,0)=Z_{\rm single}(\beta)^2,
\qquad
\lim_{\kappa\to+\infty}e^{-\kappa n}Z_{\rm double}(\beta,\kappa)
=Z_{\rm single}(2\beta).
\]

Full brute-force comparison also passes for the cases with at most 16 total spins.

At N=35 and κ=0.1, the numerical mean overlap Q/16 is approximately 0.11258, 0.61169, and 0.87053 at β=0.2,0.64,1.0. At κ=0 the overlap mean is exactly zero by flipping one entire layer. The full 112-spin N=143 bilayer was not attempted.

Figure: the third panel of `MTFT_Twists_and_Layers.png`.

## 5. Field zeros and temperature zeros

These are zeros of two different finite polynomials. The variables must be kept distinct.

### Lee–Yang field zeros: exact finite certificates

For integer satisfied/unsatisfied edge-weight ratio a/d with d=1, define

\[
P_a(y)=\sum_{m=0}^{56}\left[\sum_k D_{k,m}a^{84-k}\right]y^m,
\qquad y=e^{2\eta},\quad\beta=\tfrac12\log a.
\]

Spin reversal makes P reciprocal. Dividing by y^28 and substituting u=y+y^−1 reduces it to a degree-28 real polynomial. The exact substitution identity is checked coefficient by coefficient.

For a=2,4,8, rational Sturm isolation finds 28 simple real roots strictly inside (−2,2). Every such u yields two conjugate unit-circle roots of y²−uy+1=0. The shipped verifier independently checks 28 disjoint rational intervals with exact opposite signs at their endpoints. Their count exhausts the polynomial degree. This certifies all 56 field zeros on the unit circle in each of the three cases.

| a/d | β | Nearest positive-axis root angle, radians |
|---:|---:|---:|
| 2 | 0.346574 | 0.41934394 |
| 4 | 0.693147 | 0.06761055 |
| 8 | 1.039721 | 0.05684181 |

The angle displays are numerical; the unit-circle conclusion is algebraically certified. This supplies explicit certificates for these finite graphs, consistent with the established ferromagnetic Lee–Yang theorem, rather than a new general theorem. See [Fröhlich and Rodriguez, *Some Applications of the Lee–Yang Theorem*](https://arxiv.org/abs/1205.6643).

### Fisher temperature zeros: exact factors, numerical remaining roots

At zero field, write x=e^−2β and D(x)=Σ_k D_k x^k. Exact factorization gives

\[
D(x)=2(1+x)^{28}(1+x^2)Q_{45}(x).
\]

Therefore x=−1 has multiplicity 28 and x=±i each has multiplicity one. The factor (1+x)^28 has a graph explanation: an even subgraph of this cubic 56-vertex graph occupies at most 56 of its 84 edges. The high-temperature transform forces at least E−n=28 factors of (1+x). Its degree-56 even-subgraph coefficient is 1,784, so that multiplicity is exactly 28. This factor alone is not evidence for an additional arithmetic symmetry. Q45 has signed coefficients and is not automatically a smaller physical partition function.

The 45 Q roots were computed at 40 and 65 digits. Matching roots across precisions gives a maximum nearest-root displacement of 4.44×10^−41. At 80-digit working precision, the largest coefficient-scaled polynomial residual is 4.61×10^−67. These are numerical diagnostics, not certified complex isolating disks.

The pair nearest the positive real axis among roots with positive real part is approximately

\[
x=0.27188432269\pm0.10008127579\,i.
\]

The finite partition polynomial is positive for real x>0. The displayed complex zeros do not establish a thermodynamic-limit singularity or a result about Riemann zeta zeros.

Figure: `MTFT_Complex_Zeros.png`.

## 6. Exact-weight equilibrium sampling

Elimination can retain the conditional tables it would otherwise discard. Traversing them in reverse order assigns one spin at a time. If the two conditional integer weights are w0 and w1, a uniform integer draw in [0,w0+w1) chooses the spin with exactly the corresponding rational probability. There is no floating-point probability normalization.

With edge weights a,d and site weights v0,v1, the sampled configuration has weight

\[
W(\sigma)=a^{E-k}d^k v_0^{n-m}v_1^m.
\]

The sampler records fixed pseudorandom seeds. The algorithm yields independent equilibrium draws under ideal independent uniform integer draws; the implementation uses Python's pseudorandom generator. No Markov-chain burn-in or physical time interpretation is involved.

All configuration probabilities were checked as exact fractions for the four-spin N=6 graph with edge weights (3,2) and site weights (2,3). The product of conditional probabilities equals W/Z for every configuration and sums exactly to one.

At N=143, 20,000 configurations were drawn for each a/d=1,2,4,8 at zero field. Across all four runs, the largest standardized error of the sample mean interaction or magnetization is **1.182 standard errors**. This is a finite-sample diagnostic; it is separate from the exact conditional-weight calculation. The first 256 configurations and the full histograms are stored for each coupling.

`MTFT_Equilibrium_Samples.png` compares the sampled magnetization histograms with the finite distribution computed from the exact table. `MTFT_Equilibrium_Draws.gif` plays the first 32 draws at each coupling. Its playback speed has no physical meaning.

## Reproduction and scope

`README.md` contains the six playground commands, dependency setup, and the complete reconstruction sequence. `provenance.json` records the runtime. `SHA256SUMS` covers every delivered file other than itself.

The source archive SHA256 is:

`758f29a7a0c2ca08896a15964b3a607c7143ac1fd335bcc8f4d84a4880e53cdf`

The previous targeted release check is carried forward as `prior_release_checks.log`: 22 passed, one skipped because PARI/GP was unavailable. That was a targeted check of four test files, not the complete MTFT suite.

The main bridge to further work is now explicit: energy, magnetization, gauge twists, selected homology labels, and coupled-copy overlap can all be computed within one finite counting framework. The new statistical models do not change the separate AF09 bimodule result about vanishing one-forms. Extending the topology projection, certifying complex temperature-root disks, or increasing bilayer size would be additional computations beyond this delivered scope.
