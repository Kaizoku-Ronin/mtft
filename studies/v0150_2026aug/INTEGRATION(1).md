# mtft v0.14.1 — arithmetic area geometry, and two corrections

One new study, one patched module, 13 new tests. The study came from an
outside engine's exploration against v0.14.0; it was audited here, its
central theorem was proved and tested out of sample, and two corrections
against the package fell out of the audit.

## What to copy

```
studies/arithmetic_area_geometry.py  -> studies/
studies/arithmetic_area_geometry_ledger.json -> studies/
src/curvature.py                     -> src/mtft/curvature.py   (REPLACES)
tests/test_area_geometry.py          -> tests/
```

Windows, from the repo root:

```
robocopy .\src .\..\mtft\src\mtft curvature.py
robocopy .\tests .\..\mtft\tests test_area_geometry.py
robocopy .\studies .\..\mtft\studies *
py -m pytest tests/test_area_geometry.py -q
py -m pytest -q
```

13 tests, about 2 seconds. The study itself runs in 5 seconds and gates
9 results.

## The result worth having

`curvature.py` previously said the `{1,2,3,4}` lock at K = 1/4 held
because "atom 4 is inert", with a Pp mechanism sketch. That is now a
theorem with a stated support criterion:

> A two-dimensional discrete exponential family whose support consists
> of any number (>= 2) of statistic points on one affine line, together
> with **exactly one** point off that line, with arbitrary positive base
> weights, has K = 1/4 identically on the whole parameter manifold.

The proof is four lines once the metric is diagonalized: normalize the
line to y = 0 and the lone point to (0,1), so Z = A(u) + e^v; with
q = v - log A(u) the mixed term cancels identically and the metric
becomes the warped form ds^2 = dy^2 + cos^2(y/2) drho^2, for which
K = -f''/f = 1/4. `{1,2,3,4}` is in the class because w_4 = (1/2)log 4
gives X_4 = 2 X_2 **exactly**.

The theorem was then tested on supports that did not motivate it —
`{1,2,4,8}` and `{1,2,4,16}` are in the class for the same dyadic reason
and both land on 1/4 to better than 1e-115. The hypothesis is sharp in
both directions and both directions are gated: zero off-line points makes
the two statistics affinely dependent so det g vanishes identically, and
two off-line points destroys the lock (`{1,2,4,8,16}` gives -11,
`{1,2,3,4,5}` gives 476). A new helper `rigidity_class(points)` reports
the off-line count so a support can be classified before it is trusted.

Alongside it: `det g` equals a weighted sum of squared arithmetic
triangle areas, and K is the normalized oriented-area chirality of those
triangles about the Gibbs centroid — checked against the independent
Brioschi route, which shares no step with the cumulant determinant.

## CORRECTION CC-02 (append-only) — the cold amplitude was wrong

`CURVATURE["cold_amplitude"]` was `0.270126465305424759517602`. It is
wrong in the 16th digit. The correct value has a closed form in
log 2, log 3, log 5:

```
c = (9 L5^2)/(25 L2 L3) [1 - (9/5) L5/L3 + (4/5) L5/L2]
  = 0.27012646530542495706433719670365...
```

exposed as `curvature.cold_amplitude()`.

**Diagnosis, which is the reason to trust the correction.** The old value
was extracted from the six-atom model at beta = 200, and the six-atom
model still carries an atom-6 contamination of relative size
(5/6)^beta. At beta = 200 that is 1.46e-16; the observed error is
1.98e-16. The study reproduces the retracted value from the six-atom
model to the digit, which is what proves the diagnosis rather than
merely asserting it. The retracted value is kept in the dictionary as
`cold_amplitude_RETRACTED_CC02` and a test asserts it is still
reproducible.

**The physics behind it also changes.** v0.14.0 said the deep-cold
geometry was "the geometry of the first six integers". The cold core is
actually `{1,2,3,5}`: the rate 6/5 is 36/30, the squared metric triangle
(1,2,3) against the leading curvature triangle (2,3,5). The intermediate
candidate (2,3,4), with the smaller product 24, would give a faster
(3/2)^beta mode — but its factor cross(X_2, X_4) vanishes by the same
collinearity that produces the rigidity lock. Atom 4 shields the cold
regime, and perturbing w_4 by delta restores the suppressed mode with
amplitude 3 delta / log 3.

## CORRECTION CC-03 (append-only) — precision could be lowered silently

`finite_atom_curvature` took a fixed floor of 0.7*beta + 45 digits and
honoured a caller's lower `dps`. Both halves were unsafe. The required
precision is **support-dependent, not a function of beta**: when the
support contains a collinear triple the leading metric triangle
vanishes, det g drops to the next triangle product, and the cancellation
deepens far past any beta-based estimate. For `{1,2,4,8}` at beta = 80
the old floor of 101 digits returned about **-2e+19 instead of 1/4**,
with no error raised.

`finite_atom_curvature` is now adaptive: it recomputes at doubled
precision until two successive evaluations agree to `rtol` (default
1e-20), raising `ValueError` if it cannot stabilise. `dps` is a starting
precision that may only be raised. The exact call that used to fail
silently is a regression test.

Note this is a behaviour change, not just a bug fix: callers who passed
`dps` to make things *faster* will now get correctness instead of speed.
Nothing in the existing suite depends on the old behaviour.

## Prime boundary — a statement to fix before it propagates

The study states and gates the boundary identity in its corrected form:

```
w_n >= (log n)/n,  with equality iff n = 1 OR n is prime
```

The n = 1 case is not decorative: X_1 = (0,0) lies **on** the curve
y = x e^{-x} and 1 is not prime, so the unqualified "iff n prime" is
false. Any write-up of the prime-boundary picture should carry the
qualifier.

## Version bump

Three-way guard as usual: `pyproject.toml`, `src/mtft/__init__.py`,
`CITATION.cff`. Suggested `0.14.1`, since this is a correction plus a
study rather than new module surface.

```
## 0.14.1
Added: studies/arithmetic_area_geometry.py (9 gates) — Fisher volume as
squared arithmetic triangle area, curvature as oriented-area chirality,
the K = 1/4 rigidity theorem with out-of-sample confirmation, the
arithmetic area zeta A(s) = zeta(s)^3 det g(s,0) with fifth-order
Hagedorn wall. curvature.cold_amplitude(), curvature.rigidity_class().
Corrections: CC-02, the cold amplitude was contaminated by atom 6 at the
beta = 200 extraction point and is replaced by a closed form in
log 2, log 3, log 5; the cold core is {1,2,3,5}, not the first six
integers. CC-03, finite_atom_curvature is now adaptive in precision and
a dps argument may only raise it.
```

## Provenance

The exploration was proposed by an outside engine reading v0.14.0 and is
credited as such in the study docstring. What the audit added: the proof
of the rigidity theorem, the out-of-sample tests and both sharpness
controls, the CC-02 diagnosis with its (5/6)^beta error model, the CC-03
failure and adaptive fix, and the n = 1 correction to the prime boundary.
The generalized-variance/simplex-volume identity is classical and the
study says so; the novelty question is confined to the arithmetic
specialization, and no literature claim is made.
