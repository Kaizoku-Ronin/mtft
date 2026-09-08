# Where MTFT stands relative to the Standard Model

Assessment of the attached MTFT 0.26.0 release · 7 September 2026

**MTFT now supports a concrete Standard Model compatibility program. The checked implementation establishes substantial arithmetic geometry, field-theory models, and conditional parameter formulas. It does not yet establish their identification with the complete observed Standard Model.** This is a useful, testable position: the next work can connect fields, interactions, and observables through one specified action.

The assessment uses the attached source, the preceding Hecke-sector experiment, the six new images, and a bounded read-only numerical audit. It is not an audit of every paper, a complete package test run, or a new precision fit. The source archive and existing experimental artifacts were not modified. Findings below are appended records.

## What the available tools establish

| Layer | Present evidence | What remains for Standard Model identification |
|---|---|---|
| Arithmetic geometry | Exact Manin/cycle data, intersection form, Hecke projectors, periods and numerical Hodge structure; all 19 frozen gates passed again. | An explicit map from this internal structure to physical fields and observables. A homology dimension is not a particle count. |
| Dynamics | Hamiltonian evolution with metric-aware checks, exact symmetry selection rules, and controlled mixing experiments. | Physical field content, units, kinetic normalization, and interaction identification. The preceding leakage percentages are dimensionless model diagnostics. |
| Field theory on the surface | A finite Ising model and gauge-theory tools on a real two-dimensional modular surface, including standard U(1), SU(2), and SU(3) Yang–Mills character sums. | A demonstrated relation to the local, relativistic theory in 3+1-dimensional spacetime. The group choice in these routines is an input. |
| Four-dimensional lattice prototype | A separate four-dimensional SU(N) lattice implementation with Wilson and arithmetic Polyakov-loop terms. | Validation of group membership and updates, controlled continuum/volume limits, matter content, and matching to observables. This is a prototype to audit, not an absent capability. |
| Mass and coupling relations | Arithmetic formulas, a dimensional electroweak anchor, a calibrated Hosotani potential, conditional mass-ratio relations, and a particle reference catalog. | One consistent normalization and scale prescription, an input/output dependency ledger, and independent predictions from a common field theory. |
| Full Standard Model structure | Standard Lie-algebra utilities and particle representation labels are available. | A derived assignment of left/right-handed fields, hypercharges, generations, Yukawa matrices, anomaly cancellation, and interaction vertices. Representation strings alone are not such a construction. |

The physical comparison concerns more than masses. The Standard Model specifies gauge representations, chiral matter, Higgs/Yukawa interactions, and the relation between mass and interaction bases. The electron's left-handed component belongs to a weak doublet, while its right-handed component is a weak singlet. These are testable structural requirements; see [Tong, Electroweak Interactions](https://davidtong.org/pdfs/teaching/standard-model/standardmodel5.pdf).

It is reasonable to supply ordinary 3+1-dimensional spacetime as an explicit modeling assumption and investigate whether MTFT fixes the internal operators or effective couplings. Deriving spacetime and gravity is not a prerequisite for this conditional comparison. It must remain clear which structure is assumed and which is derived. Neutrino mass terms would also require identifying the chosen extension of the minimal renormalizable Standard Model.

## What “we have particle masses” currently means

The source contains several different routes, which must be tracked separately.

1. **An arithmetic electroweak formula.** In constants.py, SM.m_W is computed from the arithmetic alpha and sin²(theta_W) formulas together with v = 246.22 GeV. It is not copied from PDG.m_W. SM.m_Z follows from the prescribed angle. HIGGS.m_H = v gamma/(2 Omega) is a separate arithmetic formula with a dimensional anchor.
2. **A calibrated potential.** HosotaniMTFT takes sin²(theta_0) = 3/13 and chooses its potential coefficients so that the curvature reproduces HIGGS.m_H. Its gauge-mass method reuses SM.m_W. Reproducing those numbers in the potential is a consistency construction, not a second independent prediction of them.
3. **A reference particle catalog.** particles.py explicitly supplies measured/reference masses, charges, spin labels, representation strings, and kappa values. Its compute_kappa method inverts an observed mass. Reconstructing the mass from that inferred kappa is an identity.
4. **Conditional relations.** predict_tau_mass uses two measured inputs, electron and muon masses, the Koide constraint, and the larger-root choice. It is a conditional tau relation, not an independent calculation of all three lepton masses.

The read-only dependency check reproduced the default Hosotani outputs:

| Quantity | Default output | Evidence class |
|---|---:|---|
| m_W | 77.605406 GeV | Numerical evaluation of the implemented phenomenological formula |
| m_Z | 88.483777 GeV | Same |
| m_H | 125.296414 GeV | Arithmetic target reproduced by potential calibration |
| m_KK | 323.097071 GeV | Model radius converted to an inverse length/energy scale |

These values are not new experimental measurements or a completed renormalized spectrum. In particular, electroweak comparisons must distinguish low-energy alpha, running couplings, pole masses, and the renormalization convention for the weak angle. A disagreement in a bare/default calculation cannot be assigned a meaningful precision significance without that prescription. The [2026 PDG electroweak review](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-standard-model.pdf) details the relevant schemes and radiative relations.

Replacing the constants.PDG object with a sentinel that refuses every lookup did not affect this mass pipeline. Conversely, doubling only the Hosotani instance's v_ew does not double its W/Z/H masses: it still reads the global SM and HIGGS objects. Doubling the instance VEV and both global anchors does double the masses. The API therefore does not expose one unified scale control, even though the default values coincide.

## Concrete findings from this audit

**SM-AUDIT-01 — a factor-of-four Higgs quartic discrepancy.** Use the explicit convention

\[
V(H)=-\mu^2 H^\dagger H+\lambda(H^\dagger H)^2,
\qquad m_H^2=2\lambda v^2.
\]

HosotaniMTFT.higgs_self_coupling returns m_H²/(2v²) = **0.1294793867**. The public HIGGS.lambda_quartic property instead returns **0.5179175467**, exactly four times as much. The package's prediction test computes gamma²/(8 Omega²) directly, so that test can pass while missing the inconsistent public property. This is an implementation/normalization issue with a clear standard-convention remedy; it is not evidence against the exact Hecke geometry. Any intentional alternate convention would need an explicit conversion and consistent API naming.

**SM-AUDIT-02 — the legacy random SU(N) initializer does not generally return SU(N).** random_su_n divides only one row of Q by det(Q)^(1/N). Its determinant consequently becomes det(Q)^(1−1/N), which need not be 1. All eight tested SU(3) samples failed the determinant-one tolerance of 10⁻¹⁰, with errors from 0.126 to 1.401, while unitarity errors stayed below 7.5 × 10⁻¹⁶. The separate near-identity sampler passed the determinant and unitarity checks on the same eight seeds. The exact correction is to divide one row by det(Q), or the entire matrix by an appropriate Nth root. This audit records the defect without changing the release. Hot starts using the affected function should be corrected before they support physics claims.

**SM-AUDIT-03 — two W/Z ratio prescriptions require reconciliation.** The mass routine uses sqrt(10/13) = **0.8770580193**, whereas GAUGE.W_Z_ratio returns 1/(2 Omega) = **0.8816114172**, a relative difference of **0.519167%**. They cannot both be the same ratio in the same approximation and convention. A running-angle versus pole-mass interpretation could be investigated through an explicit matching calculation; it cannot simply be asserted after observing the discrepancy.

**SM-AUDIT-04 — dependent outputs are not independent evidence.** The mass/potential/catalog distinctions above were verified from source and numerical dependency checks. In particular, a close Higgs number does not independently establish the Higgs self-coupling, its couplings to fermions, or the full scalar dynamics. Those observables must follow from one normalized action.

No package files were changed. These four findings provide a concrete starting list for a subsequent correction and matching pass. The numerical record and source hashes are in audit_results.json.

## Where the Hecke results help

The preceding experiment established an exact oldspace–quartic channel compatible with both Atkin–Lehner involutions. That can constrain a future physical interaction matrix once the fields and operators are identified. It is not yet a weak charged-current vertex or a mixing angle measured in particle physics.

The good-prime Hecke operators commute. Their eigenlabels alone therefore cannot be the noncommuting generators of su(2) or su(3). Additional operators may act on multiplicity spaces or enlarged modules; their commutators, representations, kinetic terms, and compatibility with Hecke structure must be tested. The existence of standard Gell-Mann matrix utilities in the package supplies computational machinery, not a derivation of the physical gauge group from the modular curve.

The packaged AF09-style doubled-space census also reports vanishing internal one-forms for the tested involutive twists. That specific construction does not generate a nontrivial internal gauge/Higgs fluctuation sector. This is a scoped finding about that route, not a no-go theorem for MTFT. It was inspected in the release record here, not rerun as part of this bounded audit. The 26-dimensional homology and 52-dimensional doubling should not be interpreted as counts of physical particles.

## Turning the six images into an experiment sequence

| Image | Productive question | A meaningful pass condition |
|---|---|---|
| Bosons | Which modes are vector, scalar, or fermionic fields? | Lorentz representations, statistics, gauge transformations, and physical polarizations are explicit; a photon candidate has two transverse helicities. A graviton belongs to a separate gravity extension. |
| Speed of light / Maxwell | Does the proposed electromagnetic sector have a common causal propagation law? | A gauge-consistent wave operator, correct physical modes, and controlled long-wavelength dispersion. A numerical eigenmode on the internal curve is not automatically an electromagnetic wave. |
| Lorentz force | Does the same charge appearing in the interaction govern motion? | The coupling yields q(E + v × B) in the classical limit, without independently fitting a second charge. The image's qvB sin(theta) is only the magnetic-force magnitude. |
| Bar magnets | Does the same theory recover dipole response? | Magnetic moment, dipole field, torque, and response follow from the specified action/state. Electron spin response should precede a macroscopic ferromagnet simulation. |
| Bohr atom | Do the same charge and mass yield an atomic spectrum? | Leading hydrogenic levels and later spin/fine-structure corrections follow from the model's Hamiltonian. Inserting fitted values into the textbook formula alone is a baseline, not new dynamical evidence. |
| Phase-conjugate resonator | Can boundary driving and coherent waves be modeled consistently? | First specify Maxwell boundary conditions and momentum/energy flow. The displayed F = ma = 2c Delta-nu is dimensionally inconsistent as a force law: c Delta-nu has units of acceleration. It should not be adopted as an established law. |

The polarization and local-gauge requirements are explained in [Tong, Symmetries](https://davidtong.org/pdfs/teaching/standard-model/standardmodel1.pdf). For atomic benchmarks, [NIST's discussion of the fine-structure constant](https://physics.nist.gov/cuu/Constants/alpha.html) links the leading atomic energy scale to alpha and the electron mass. Accurate comparisons additionally involve reduced mass, relativistic effects, and radiative/nuclear corrections.

## The next physical target

The most focused target is an **electron–photon effective theory**. In natural units, a reference form is

\[
\mathcal L=-\frac14F_{\mu\nu}F^{\mu\nu}
+\bar\psi(i\gamma^\mu D_\mu-m_e)\psi
+\sum_{d>4,i}\frac{C_i^{(d)}}{\Lambda^{d-4}}\mathcal O_i^{(d)}.
\]

This is a proposed matching target, not an action already derived from MTFT by this audit. The work is to identify the photon, spinor, charge operator, mass term, and any higher-order coefficients in the MTFT construction, then normalize them consistently. If spacetime or the Dirac field is supplied externally, record it as an assumption.

A practical sequence is: correct and gate the two implementation defects; reconcile the mass/coupling conventions and input ledger; define a candidate matter-plus-gauge action; test its symmetries and physical modes; derive charge response and magnetic moment; compute atomic and scattering observables using the same fixed inputs. A later chiral electroweak extension needs explicit left/right-handed multiplets and anomaly cancellation before detailed mass fitting.

This approach can determine whether MTFT supplies a compatible internal geometry, relations among Standard Model parameters, or additional effective interactions constrained by experiment. It gives the existing computational work a precise physical role without assuming the answer in advance.
