# Matter-01: a working U(1) reference and an audited arithmetic bridge

1 October 2026 · unreleased candidate on the supplied MTFT v0.34.0 source

Matter-01 adds a massive tree-level electron–positron annihilation engine, an
offline interactive lab, and an explicit audit of the proposed MTFT charged
sector. It calculates amplitudes and event cross sections for
\(e^-e^+\rightarrow\gamma\gamma\), with independently checked conventions.
It also establishes conditional normalization and gauge–breathing identities.
The electron mass and electromagnetic coupling are **supplied reference
inputs**; the reference QED action assumes a massless photon. These are not
outputs of the MTFT vacuum.

The most consequential new MTFT result is a selection warning supported by
exact arithmetic: the background with the lowest archived net-index extras
score need not have the smallest cohomology spaces. A simpler collinear
control merits the next explicit mode calculation.

## What is implemented

| Component | Result | Evidence class |
|---|---|---|
| `mtft.matter.qed` | On-shell massive spinors, both exchange amplitudes, polarization sums, CM/boosted kinematics, differential and integrated event cross sections | REFERENCE / DIAGNOSTIC |
| `mtft.matter.oracle` | Independent invariant and CM formulas, analytic total and target-rest total | REFERENCE / DIAGNOSTIC |
| `mtft.matter.normalization` | Same-charge kinetic/vertex normalization, conditional electroweak matching and breathing gauge prefactor | EXACT conditional identities; numerical evaluation DIAGNOSTIC |
| B3 vacuum audit | A frozen archived witness, slopes, net indices, charges, anomaly checks, RR lower bounds and collinear comparison | EXACT conditional arithmetic; physical completion OPEN |
| Offline lab | Kinematics, two Feynman terms, interference/gauge evidence, angular distributions, phase illustration, breathing coefficient and completion ledger | Visualization of the stated reference and conditional formulas |
| Atlas ledger | Links the reference electron–photon vertex to this benchmark while retaining every MTFT mapping as `unmapped` | Evidence bookkeeping |

There is no real-time fermion/photon field evolution, event generator,
positronium calculation, proton/neutron model, or particle condensation in
this milestone. The lab's phase illustration explains a U(1) transformation;
it is not a generated quantum state. Feynman diagrams encode amplitude terms,
not measured microscopic trajectories.

## The QED benchmark

The reference action is ordinary QED,

\[
\mathcal L=-\frac14F_{\mu\nu}F^{\mu\nu}
 +\bar\psi(i\gamma^\mu D_\mu-m_e)\psi,
\qquad e^2=4\pi\alpha.
\]

Conventions are \((+,-,-,-)\), natural units, an average over the four incoming
spin combinations, and a sum over four outgoing physical polarization
combinations. The electron and positron have opposite charges and equal
mass. The two electron-exchange diagrams must be added before squaring.
Replacing a photon polarization by its momentum cancels the **full**
amplitude; either individual diagram can have a nonzero contraction.

For \(\beta=\sqrt{1-4m_e^2/s}\), the event cross section is

\[
\sigma=\frac{2\pi\alpha^2}{s\beta}
\left[\frac{3-\beta^4}{2\beta}
\log\frac{1+\beta}{1-\beta}-(2-\beta^2)\right].
\]

Full-sphere differential event distributions include the identical-photon
factor \(1/2!\). A photon-inclusive histogram has twice the event integral.
An absolute rate additionally needs an incident luminosity: \(R=\mathscr L\sigma\).
The retained mass regulates the collinear regions; a massless fixed-angle
formula cannot be integrated there as a finite physical total.

Frozen PDG 2026 reference inputs are
\(m_e=0.00051099895069\,\mathrm{GeV}\),
\(\alpha^{-1}=137.035999177\), and
\(1\,\mathrm{GeV}^{-2}=389379.3721\,\mathrm{nb}\).
CLI overrides remain supplied inputs and carry separate provenance.
The calculation uses a fixed coupling and Born amplitudes. Radiative
corrections, running, Coulomb enhancement, positronium and material effects
are outside its scope. The formal Born threshold limit is not a complete
physical threshold model.

See [the derivation](QED_DERIVATION.md) for the spin-averaged invariant formula,
event-counting derivation, limits and primary references.

## Conditional geometry–matter bridge

Assume a constant, unbroken photon zero mode, a positive charged-mode Gram
matrix \(K\), and equal charge \(q\) throughout a block. With
\(C^\dagger K C=I\), its vertex normalizes to

\[
C^\dagger(gqK)C=gqI.
\]

Thus changing a same-charge mode basis does not change its electric charge.
This identity does not identify an electron, derive a Yukawa mass, choose
the physical photon, or fix the parent coupling. In an unwarped reduction,
\(g_4=g_{4+d}/\sqrt{\mathrm{Vol}}\), with compatible dimensions and both
quantities supplied by a selected higher-dimensional action.

For uniform Einstein-frame breathing,

\[
ds_{4+d}^2=e^{-d\sigma}ds_{4,E}^2+e^{2\sigma}ds_{\rm int}^2,
\qquad \chi=M_P\sqrt{d(d+2)/2}\,\sigma,
\]

the external gauge term has

\[
\mathcal L_{\rm gauge}=-\frac14Z(\chi)F^2,
\qquad Z=e^{\sqrt{2d/(d+2)}\,\chi/M_P}.
\]

For C3X, \(d=2\) and \(Z=e^{\chi/M_P}\). For a uniform \(d=4\)
reduction, the exponent is \(2/\sqrt3\); this is a different ansatz.
In a **constant** breathing background,
\(g_{\rm eff}/g_0=Z^{-1/2}\). For varying \(\chi\), retain the
\(ZF^2\) action; a field-dependent canonical rescaling introduces
derivative terms. Matter-01 does not evolve this coupled system.

A null plane wave can have \(F^2=0\) while carrying positive stress-energy.
Its vanishing direct \(\chi F^2\) source does not remove gravitational
backreaction. Superpositions and standing waves can have nonzero \(F^2\).
The helper reports canonical Maxwell energy with \(Z=1\); unrescaled energy
and flux in the breathing action acquire a factor \(Z\).

## The selected-vacuum audit changes the next step

The archived B3 witness is D-flat at \((A_E,A_X,\epsilon)=(8,3,1)\)
and has three net family indices. Its recorded extras score is 31, a sum of
different index-bookkeeping units. It is not a physical particle count or a
global optimum. Matter-01 checks the witness and a small normal slice, not
the entire previous enumeration.

Riemann–Roch gives
\(h^1(L)\geq\max(0,-\chi(L))\), independently of whether \(h^0\)
or \(h^2\) vanishes. This exposes a substantial difference:

| Hom sector | Score-31 baseline, both orientations | Collinear control, score 36 |
|---|---:|---:|
| Quark-doublet line c/L | \(h^1\ge222,225\) | \(h^1\ge3,6\) |
| Higgs-type line L/B | \(h^1\ge897,891\) | \(h^1\ge21,15\) |

These are cohomology lower bounds, not exact dimensions or physical
massless-particle counts after couplings. They show why net indices alone
cannot select a minimal spectrum. The simpler collinear control is an
appropriate next calculation, but its viability is also OPEN.

The hypercharge assignment closes the lepton/Higgs bundle triangles and
assigns the expected charges. It does not establish nonzero normalized
Yukawas. Likewise \(\mathrm{Tr}(Yc_1)=0\) does not certify a massless photon.
The axion/Stückelberg mass matrix, gauge kinetic mixing and Higgs mass matrix
must have the required common electromagnetic null direction. This remains
uncomputed. The B3 construction and six-dimensional C3X breathing benchmark
are separate; no common parent has been established.

See [the complete audit](VACUUM_AUDIT.md) for all witnesses, exact bounds and
eleven OPEN completion gates.

## Verification and reproducibility

The focused Python gate passes 84 tests, and the offline JavaScript gate
passes 54 checks. The recorded gate also contains 45 exact vacuum checks, 13 conditional algebraic
reduction identities, eight direct numerical breathing samples, and 30
sampled mode-basis changes. QED has 12 structural numerical gates and 14
independent gates. The oracle compares 126 spinor evaluations and 42
independently boosted configurations with scalar formulas; the maximum
relative discrepancies are approximately \(1.2\times10^{-14}\) and
\(3.7\times10^{-14}\). The integrated engine total agrees with the analytic
oracle to approximately \(1.7\times10^{-14}\) on the declared grid.

These numerical residuals are **finite-sample diagnostics**, not global
certificates. JavaScript is independently compared with 25 Python angular
points, and its event integral is checked against the analytic total. Mock
DOM/canvas verification exercises controls and rendering calls; actual
browser layout was not verified in the recorded environment. The full MTFT
release suite and hours-long archived studies are not part of this focused
gate. Exact test counts, timings and dependency versions are in
`results/run_summary.json` and `results/focused_tests.log`.

From the candidate repository root:

```bash
python -m pip install -r studies/matter01/requirements.txt
python studies/matter01/run_all.py
```

To use the public CLI after installing the candidate (`pip install -e .`):

```bash
python -m mtft.matter benchmark --sqrt-s-gev 0.01 --cos-theta 0.37 --output benchmark.json
python -m mtft.matter report --output qed_report.json
```

The prebuilt lab opens locally without installation or network resources.

## Next milestone

The next MTFT-specific gate should construct actual holomorphic bundles at
a specified blow-up point for the simpler collinear control, compute the
lepton/Higgs cohomology and kinetic Grams, and evaluate a nonzero normalized
electron Yukawa. In parallel, select the axion and gauge kinetic action and
check the photon kernel. Those outputs can replace the supplied coefficients
of the working U(1) engine.

For computational coverage, proceed from electromagnetic U(1) to the joint
electroweak \(SU(2)_L\times U(1)_Y\) sector, then add \(SU(3)_c\).
These are interacting factors, not successive group inclusions.
Protons and neutrons require nonperturbative QCD bound-state calculations;
adding gauge vertices alone will not construct them. Oloid/triangle scalar
identities and the earlier dynamics remain useful diagnostics, but neither
supplies those missing fermion or QCD results. No Feigenbaum cascade was
established in the previous AdS scan.

## Primary references

- [PDG 2026 physical constants](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-phys-constants.pdf).
- [Kaplunovsky, annihilation lecture notes](https://web2.ph.utexas.edu/~vadim/Classes/2024f-qft/annihilation.pdf).
- [Geant4 positron–electron annihilation reference](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/electron_incident/annil.html).
- [Tong, quantum field theory notes](https://davidtong.org/pdfs/teaching/quantum-field-theory/qft.pdf).
