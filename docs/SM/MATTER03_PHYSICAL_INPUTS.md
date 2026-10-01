# Matter-03: which physical inputs are still missing?

Matter-03 separates classical Higgs-branch existence, selection of a small
light spectrum, and conversion of holomorphic couplings into physical masses.
Passing one gate does not pass the other two. The underlying construction is
the specified blown-up surface B3 of Matter-02; the older two-real-dimensional
C3X breathing model supplies neither its metric nor its mode norms.

## Higgs vacuum and Standard Model metastability

A Higgs vacuum is a background expectation value at a stationary point of a
specified potential. The words *Higgs vacuum* do not mean *false vacuum*.
Metastability means that a locally stable state can decay to another state of
lower energy. It requires the effective potential, its applicable energy
range, and a tunneling calculation. A classical flat family of supersymmetric
backgrounds does not establish any of these features.

Extrapolating the Standard Model to high scales with representative central
inputs can favor a metastable electroweak vacuum. That conclusion depends on
top-mass extraction, the strong coupling, matching, running, and assumptions
about new physics. It is not an experimental determination that our universe
occupies a false vacuum. A universal measured lifetime such as 10^65 years
cannot be inferred from the short Gemini summary.

Primary references checked on 2026-10-01:

- [Vacuum Stability in the Standard Model and Beyond,
  arXiv:2401.08811v4](https://arxiv.org/html/2401.08811v4), discussion of input
  uncertainties and the still unresolved stability determination.
- [Renormalization-group analysis of the SM: Loops, uncertainties, and vacuum
  stability, Phys. Rev. D 113, 036018 (2026)](https://journals.aps.org/prd/abstract/10.1103/lgnl-2nxv),
  theoretical uncertainties in running and matching.
- [Stephen Martin, A Supersymmetry Primer, hep-ph/9709356v7](https://arxiv.org/pdf/hep-ph/9709356),
  the standard superpotential, kinetic, auxiliary-field and Higgs-potential
  framework. This reference does not identify B3 with the MSSM.

## Field-coordinate mass convention

For a lower-index four-dimensional superpotential tensor use

\[
 W=\kappa\sum_{abc}Y_{abc}L_aE_bH_c,
 \qquad K_{\rm quad}=\sum_{X=L,E,H}X^\dagger G_X X.
\]

Every supplied G must be positive Hermitian. Define raw fields X=C_X X_can
with C_X^dagger G_X C_X=I. With a unit canonical Higgs direction n and neutral
expectation value (v_H/sqrt(2))n, the mass matrix is

\[
 M_{\rm can}=\frac{\kappa v_H}{\sqrt2}
 C_L^T\Big[\sum_cY_c(C_Hn)_c\Big]C_E.
\]

The transpose on the holomorphic left slot is ordinary, not conjugate.
Nonnegative singular values of M are tree-level Dirac masses in this
specified effective theory. Invertible kinetic normalization preserves rank
for nonzero kappa; it does not supply a mass scale or remove unwanted modes.
When changing the raw Higgs basis the canonical vector n must also be
transported, so that the same physical expectation value is compared.

The historical `surface.hym.normalise_yukawa` uses a section-basis tensor with
a dual/target Higgs index. It is unchanged. The new `mtft.matter.kinetics` API
declares the lower-index field convention explicitly.

## Required data

| Gate | What is needed | Current result |
|---|---|---|
| Classical Higgs background | A holomorphic canonical-twisted pair solving F/D equations | A paired branch exists conditionally on the stated compact bulk parent |
| Higgs scale | A derived potential or other specified mechanism selecting an expectation value | Common scaling remains a modulus |
| Light families | The deformed zero-mode complex and a mechanism lifting excess modes | The present restricted pair leaves far too many charged classes |
| Kinetic normalization | Smooth surface and Hermitian background metrics; global harmonic representatives; parent kinetic coefficients | None of the B3 charged mode Gram matrices has been computed |
| Full coupling | Complete tensor and overall superpotential coefficient kappa | Matter-02 supplies a rank-three residue slice, not the complete operator |
| Photon | Higgs and axion/Stueckelberg kernels with a common physical dictionary | Q is preserved by the classical pair; the axion gate remains open |

The split Matter-02 slice uses dimensions 18, 44, 10, so its raw metric problem
would require 18-by-18, 44-by-44, and 10-by-10 matrices. Once the background is
deformed, physical light modes must be computed in that deformed complex and
given their own metrics. Normalizing a convenient split minor does not select
the surviving light families. Historical curve Petersson, curve HYM and M2
data are not B3 surface mode metrics.

## Reusable numerical API

The following is an INPUT example in arbitrary units. It illustrates how
different positive metrics change a known rank-three coupling. It is not a
lepton mass prediction or a selected family sector.

```python
import math
import numpy as np
from mtft.matter.collinear import yukawa_witness
from mtft.matter.kinetics import canonical_mass_matrix, mass_singular_values

Y = np.array(yukawa_witness().residue_matrix, dtype=complex)[:, :, None]
M = canonical_mass_matrix(
    Y, np.diag([1, 4, 25]), np.diag([1, 9, 16]), np.array([[1]]),
    np.array([1]), vev=math.sqrt(2), kappa=1,
)
print(mass_singular_values(M))  # [1, 1/8, 1/15], arbitrary input units
```

No observed mass, Standard Model vev, or parent normalization is silently
inserted by this API. Kinetic calculations are finite numerical diagnostics;
the exact basis-covariance identity is proved in the accompanying derivation.

The next family-selection investigation should change or augment the matter
background and recompute its full charged complex before fitting masses.
Additional interactions, flat twists, extensions, defects or a different
parent are possible inputs to test, not assumed repairs. Any proposal must
retain an admissible Higgs vacuum and photon while explicitly accounting for
the modes it lifts. Only then is the expensive surface-metric solve aimed at
an identified low-energy sector.
