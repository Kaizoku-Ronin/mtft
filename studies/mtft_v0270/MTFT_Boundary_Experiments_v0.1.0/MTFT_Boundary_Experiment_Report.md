# Boundary experiments on MTFT 0.26.0

Experiment bundle v0.1.0 · 5 September 2026

The first boundary experiment is complete: we can compute how much eight
boundary spins reveal about a fourteen-spin interior, choose which boundary
spins to observe, and measure what changes when one boundary bond is reversed.
Every configuration contributes through exact integer counts; probabilities,
entropies, and numerical optimizations are evaluated in floating point.

The main result is an intermediate-temperature information peak. With all
bonds ferromagnetic, the boundary reveals about **1.571 bits at β ≈ 0.612**.
At lower temperatures the interior becomes much easier to predict, while the
information added by the boundary tends toward one shared orientation bit.
Sensor placement and the assumed bond model both matter.

## The experiment we fixed

We use the 56-spin, 84-edge dual graph generated from `cell_complex(143)` in
the supplied MTFT 0.26.0 source. Its single self-loop contributes a constant
energy and remains in the partition-function normalization. There is no
external magnetic field. Coupling magnitudes and Boltzmann constant are 1,
so increasing β means cooling the system.

| Region | Vertex labels | Count |
|---|---|---:|
| Interior A | 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 35, 51, 52 | 14 |
| Boundary B | 6, 18, 22, 24, 36, 40, 50, 53 | 8 |
| Exterior C | All other vertices | 34 |

There are 17 A–A, 8 A–B, 16 B–C, and 43 C–C edges; no A–C or B–B edges.
The partition was chosen before the temperature and defect calculations from
a graph-Laplacian spectral sweep. It is one connected example, without a
certificate of optimal separator size. The compact modular curve itself has
no preferred boundary: B is a chosen separator in its finite graph.

The eight cut edges match B, in the order above, to contact spins
S = [7, 14, 17, 51, 35, 13, 10, 52] inside A. The other six interior spins are
[8, 9, 11, 12, 15, 16]. The defect is also fixed in advance: reverse the first
A–B edge in source order, edge 7, joining interior 7 to boundary 6. We compare
that model with the all-ferromagnetic model.

## What a boundary guarantees

For this local Gibbs model, absence of A–C edges gives the factorization

\[
P(A,C\mid B)=P(A\mid B)P(C\mid B).
\]

Once the boundary is known, observing the exterior supplies no additional
information about the interior. This is the classical separator property of
a local Gibbs distribution. Mutual-information area bounds for local
statistical models have an established literature; see
[Wolf, Verstraete, Hastings and Cirac](https://arxiv.org/abs/0704.3906).

For this eight-bit separator, data processing gives

\[
I(A:C)\leq I(A:B)\leq H(B)\leq 8\ \text{bits}.
\]

The experiment computes I(A:B), interior reconstruction, and every observed
subset of B. It does not directly enumerate I(A:C). The factorization is
structural, with an independent small-graph control described below.

At any finite β, all interior configurations retain positive conditional
probability. Knowing the boundary therefore determines an interior probability
distribution, not a unique microstate. One chosen cut also cannot establish
how information scales with boundary size across a family of regions.
These are finite classical experiments inspired by the holography image;
they do not establish gravitational holography or the black-hole entropy
coefficient. The relevant gravitational framework is reviewed by
[Bousso](https://arxiv.org/abs/hep-th/0203101).

## Exact compression makes the experiment practical

Let D_s[k] count interior configurations with contact pattern s and k
disagreements among the 17 internal edges. Each contact pattern has exactly
64 compatible interior states. Let G_b[k] count exterior configurations at
fixed boundary b, with k disagreements among the 59 exterior and B–C edges.
Each exterior coefficient row sums exactly to 2^34.

For q = exp(−2β), the entire calculation reduces to

\[
P(s,b)\propto D_s(q)G_b(q)q^{\operatorname{popcount}(s\oplus b\oplus d)}.
\]

Here d is zero for clean bonds and one for the selected reversed bond, in the
declared boundary bit order. The physical partition function includes the
common factor exp(84β). Exterior elimination has a maximum observed neighbor
width of 4, an upper bound from the chosen elimination order, not a proof of
optimal treewidth. The largest spin table contains 32 entries.

For any observed boundary subset O, the contact spins shield the six remaining
interior spins, giving

\[
I(A:O)=I(S:O),\qquad H(A\mid O)=H(A\mid S)+H(S\mid O).
\]

Thus the information calculations use a 256 × 256 contact/boundary joint table,
while retaining the entropy and individual predictions for all fourteen
interior spins. For a specific observed outcome, the entropy of the six inner
spins must be averaged with the posterior contact distribution. The explorer
and CLI perform that reweighting explicitly.

## Temperature response

Both models were evaluated at 48 temperatures in 0 ≤ β ≤ 3, with all 256
sensor subsets at each temperature: **24,576 model/temperature/subset cases**.

| Clean model β | Interior entropy H(A), bits | I(A:B), bits | H(A\|B), bits | Mean spin prediction accuracy |
|---:|---:|---:|---:|---:|
| 0 | 14.000 | 0.000 | 14.000 | 50.00% |
| 0.30 | 12.809 | 0.497 | 12.312 | 60.76% |
| 0.64 | 6.784 | 1.549 | 5.236 | 86.85% |
| 1.00 | 1.575 | 1.029 | 0.546 | 99.40% |
| 3.00 | 1.000007 | 1.000000 | 0.0000067 | 99.999998% |

All conditional entropies and accuracies in this table are averaged over
boundary outcomes. Accuracy is the mean Bayes success probability for each
individual spin. It differs from reconstructing all fourteen spins correctly
at once, whose averaged probability at β = 0.64 is only 40.895%.

The largest sampled information peak was refined numerically to β =
0.612034912, I(A:B) = 1.570738751 bits for clean bonds. For the reversed bond,
the corresponding peak is β = 0.639253283 and 1.615752449 bits. These are local
numerical refinements of the largest grid peaks, without certified global
extremum enclosures.

The clean model has two exact ground states, related by global spin reversal.
In the zero-temperature limit, the boundary identifies their common orientation:
I(A:B) tends to 1 bit, and H(A|B) tends to zero. At β = 0 the interior consists
of fourteen independent fair spins, so boundary observations reveal nothing.
The intermediate peak sits between those two limits. See
`Boundary_Temperature_Response.png` and its vector SVG copy.

## Where to put the sensors

At β = 0.64 in the clean model, boundary spins 24 and 40 tie for the best
individual sensor, each revealing approximately 0.532757 bits. The best pair
is **6 + 22**, revealing 0.824711 bits. It excludes both best singles.
A greedy method beginning with a best single selects 24 + 40 instead, giving
0.818424 bits. Its shortfall is small but measurable, about 0.006287 bits.

| Sensor budget | Largest information, bits | Greedy information, bits |
|---:|---:|---:|
| 1 | 0.532757 | 0.532757 |
| 2 | 0.824711 | 0.818424 |
| 3 | 1.011487 | 0.985982 |
| 4 | 1.156227 | 1.134949 |
| 6 | 1.379393 | 1.362768 |
| 8 | 1.548670 | 1.548670 |

These are exhaustive comparisons over the finite set of sensor masks, using
numerical entropies and a 1e-9-bit tie tolerance. They do not certify that
floating-point near-ties are exact mathematical equalities. Selecting sensors
in advance also differs from choosing the next sensor after seeing an outcome;
the explorer offers a one-step conditional recommendation for the latter.
See `Boundary_Sensor_Placement.png`.

## A surprising boundary need not reduce uncertainty

At β = 0.64, observing all eight boundary spins as +1 occurs with probability
21.5816%. It leaves 3.687180 bits of interior uncertainty and gives 93.9434%
mean individual-spin accuracy.

The alternating boundary pattern +−+−+−+−, in the declared vertex order,
occurs with probability only 0.0167957%. It leaves **8.313285 bits**, exceeding
the unobserved interior entropy of 6.784461 bits. That rare outcome selects
a more disordered conditional interior distribution. There is no contradiction:
conditioning reduces entropy on average, not necessarily for each outcome.
The average after observing all of B is 5.235792 bits.

## The reversed bond separates information from predictability

| Quantity at β = 0.64 | Clean bonds | Reversed boundary bond |
|---|---:|---:|
| H(A), bits | 6.784461 | 7.853468 |
| I(A:B), bits | 1.548670 | 1.615740 |
| H(A\|B), bits | 5.235792 | 6.237728 |
| Mean spin accuracy, averaged | 86.8489% | 81.7929% |
| Full-configuration accuracy, averaged | 40.8946% | 26.6870% |

The defect increases information revealed while also increasing the amount
left unknown. It raises the interior's initial entropy by more than it raises
the boundary information. If the true model contains the defect but we use
clean-model posterior predictions, mean spin accuracy falls further to
77.9731%. The expected conditional divergence from the correct posterior is
0.903627 bits. These quantities evaluate model misspecification, without
estimating the bond from data. See `Boundary_Outcomes_and_Defect.png`.

## Verification and limits

The exact coefficient checks include every interior fiber count, all 256
exterior fiber counts, spin reversal, six alternate exterior elimination
orders, reconstruction of the whole-graph density of states in both models,
and agreement with independent integer partition calculations at edge-weight
ratios 1, 2, 4, and 8. The graph matches the supplied MTFT source exactly.

A separate 16,384 × 256 interior/boundary calculation at β = 0.64 checks
interior entropy, spin marginals, and five sensor masks in both models against
the compressed route, agreeing within 1e-11. All 24,576 sensor cases obey the
information bounds, contact-entropy floor, and monotonicity of averaged
information and optimal prediction accuracy when an observation is added.

An independent six-spin model verifies conditional independence by an exact
integer cross-product identity. Adding an A–C shortcut breaks that identity,
as intended, and produces conditional mutual information of 0.063839 bits.
This negative control checks that the calculation detects a failed separator.

The JavaScript engine passed 1,064 scalar comparisons across forty fixtures,
including hot and cold limits, mixed observations, both bond models, and the
one-step sensor advice; its largest absolute discrepancy was 1.38 × 10^-13.
Three CLI examples and outcome spin-reversal checks passed. The exported
scientific figures were visually inspected. Interactive browser layout and
theme checks could not be completed because local previews were blocked by
the workspace browser URL policy.

The finite integer counts and their identities are exact. Entropies, posterior
probabilities, numerical peaks, and sensor rankings are diagnostic numerical
results, not interval-certified quantities. We used no Monte Carlo estimates.
This first experiment establishes what can be learned for one fixed graph
separator and two specified bond models. A future boundary study could vary
the region and separator size to test scaling; that experiment has not been
performed here.
