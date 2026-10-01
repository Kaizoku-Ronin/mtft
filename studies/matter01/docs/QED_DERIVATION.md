# Matter-01: independent massive QED reference

Status: **REFERENCE / DIAGNOSTIC**, not an MTFT prediction. This calculation is
the independent analytic oracle for the explicit-spinor engine. It contains no
imports of that engine's gamma matrices, spinors, amplitudes, or kinematics.
Both calculations implement established Born QED, so agreement validates the
implementation rather than deriving QED or its parameters from arithmetic.

## Inputs and conventions

Use natural units $\hbar=c=1$, metric $(+,-,-,-)$, and $e^2=4\pi\alpha$.
Frozen reproducibility inputs are $m=0.00051099895069$ GeV and
$\alpha=1/137.035999177$. These PDG 2026 numerical inputs are frozen for
reproducibility, with no subsequent silent update. Results in GeV$^{-2}$ convert to nb by multiplying by
389379.3721. The mass and electromagnetic coupling are supplied independently
of MTFT geometry, Higgs overlaps, and the AdS breathing action.

Incoming electron and positron momenta are $p_-$ and $p_+$; outgoing photons
are $k_1,k_2$. Define

$$s=(p_-+p_+)^2,\quad t=(p_--k_1)^2,\quad u=(p_--k_2)^2,$$

with $s+t+u=2m^2$. In the center-of-momentum frame,

$$\beta=\sqrt{1-4m^2/s},\quad c=\cos\theta,\quad
T=t-m^2=-\frac{s}{2}(1-\beta c),\quad
U=u-m^2=-\frac{s}{2}(1+\beta c).$$

Here $\theta$ is the angle from the incoming electron to photon 1.

## Amplitude and independent scalar reduction

The two electron-exchange diagrams have the common tree amplitude, up to a
physically irrelevant overall sign,

$$\mathcal M=e^2\bar v(p_+)\left[
\slashed\epsilon_2^*\frac{\slashed p_--\slashed k_1+m}{T}\slashed\epsilon_1^*
+\slashed\epsilon_1^*\frac{\slashed p_--\slashed k_2+m}{U}\slashed\epsilon_2^*
\right]u(p_-).$$

Each diagram separately has a nonzero longitudinal contraction; they cancel
in the sum. Therefore testing one diagram's Ward contraction as zero is an
incorrect gauge test. The full amplitude must vanish when either photon
polarization is replaced by its momentum.

Using the spin completeness relations and contracting the full, gauge-invariant
amplitude over photon polarizations gives the spin-average scalar result

$$\overline{|\mathcal M|^2}=2e^4\left[
\frac{T}{U}+\frac{U}{T}+
\frac{4m^2s}{TU}-4m^4\left(\frac1T+\frac1U\right)^2\right].$$

This averages over four initial spin combinations and sums over four final
physical polarization combinations. It is symmetric in $T,U$. Substitution of
the CM kinematics gives a separately evaluated expression

$$\overline{|\mathcal M|^2}=4e^4 A(\beta,c),\qquad
A=\frac{3-2\beta^2+\beta^2c^2}{1-\beta^2c^2}
-\frac{2(1-\beta^2)^2}{(1-\beta^2c^2)^2}.$$

These expressions retain the electron mass. The massless approximation fails
near threshold and does not supply the physical collinear cutoff in the total.

## Event counting and exact total

For the full sphere $-1\le c\le1$, $0\le\phi<2\pi$, the final photons are
identical, so the phase space must contain $1/2!$. The event differential rate is

$$\frac{d\sigma_{\rm event}}{d\Omega}
=\frac{\overline{|\mathcal M|^2}}{128\pi^2s\beta}
=\frac{\alpha^2}{2s\beta}A(\beta,c).$$

$$\frac{d\sigma_{\rm event}}{dc}
=\frac{\overline{|\mathcal M|^2}}{64\pi s\beta}.$$

Equivalently, remove the $1/2!$ and integrate only one hemisphere. A
single-photon inclusive distribution over the full sphere counts two photons
per event, so its integral is $2\sigma_{\rm event}$. These conventions must
not be mixed.

To integrate, let $L=\log[(1+\beta)/(1-\beta)]$. The elementary integrals are

$$I_1=\int_{-1}^{1}\frac{dc}{1-\beta^2c^2}=\frac{L}{\beta},\qquad
I_2=\int_{-1}^{1}\frac{dc}{(1-\beta^2c^2)^2}
=\frac1{1-\beta^2}+\frac{L}{2\beta}.$$

Since $\beta^2c^2/(1-\beta^2c^2)=1/(1-\beta^2c^2)-1$,

$$\int_{-1}^{1}A\,dc
=\frac{3-\beta^4}{\beta}L-2(2-\beta^2).$$

Consequently,

$$\boxed{\sigma_{\rm event}
=\frac{2\pi\alpha^2}{s\beta}
\left[\frac{3-\beta^4}{2\beta}L-(2-\beta^2)\right].}$$

No luminosity is supplied by this calculation. For specified incident
luminosity $\mathscr L$, the Born event rate is $R=\mathscr L\sigma$; the
photon production rate for this exclusive two-photon process is $2R$.

## Independent rest-target total and limiting checks

Geant4 documents the target-rest Heitler cross section per electron. Evaluate
that formula using

$$\gamma_{\rm lab}=\frac{s}{2m^2}-1
=\frac{1+\beta^2}{1-\beta^2},\qquad r_e=\alpha/m.$$

It then equals the CM event total above. The oracle evaluates the two forms
without sharing a total-rate helper, and numerically integrates the differential
formula as a third total-rate route.

The threshold limits are

$$\overline{|\mathcal M|^2}\longrightarrow4e^4,\qquad
\sigma\beta\longrightarrow\frac{\pi\alpha^2}{2m^2},\qquad
\sigma v_{\text{Møller}}\longrightarrow\frac{\pi\alpha^2}{m^2}.$$

In CM, $v_{\text{Møller}}=2\beta$. It is a flux velocity and can exceed 1;
it must not be interpreted as one particle's physical speed. The relativistic
relative speed is $2\beta/(1+\beta^2)$ and agrees at leading threshold order.
At $\beta=0$ the cross section is singular because the incident flux is zero;
the finite annihilation rate coefficient has the last limit instead.

At large $s/m^2$,

$$\sigma\sim\frac{2\pi\alpha^2}{s}[\log(s/m^2)-1],\qquad
\frac{d\sigma}{d\Omega}\sim\frac{\alpha^2}{2s}\frac{1+c^2}{1-c^2}$$

at fixed angles away from the beam axis. Retaining finite $m$ is necessary
before integrating through the forward/backward regions.

## Scope and numerical protocol

The primary verification grid is $\beta=0.05,0.2,0.5,0.8,0.95,0.99$ and
$c=-0.999,-0.95,-0.5,0,0.37,0.8,0.99$. A tight quadrature compares the
analytic CM total, the target-rest total, and the explicit-spinor engine total.
Additional limit checks exercise stable $\beta\to0$ arithmetic and the
high-energy asymptotic trend. Gates are finite-grid numerical diagnostics,
not a global interval-arithmetic certification.

Born QED neglects radiative corrections, Coulomb enhancement, bound-state
formation, material effects, weak interactions, and running of the imported
coupling. Near threshold $\alpha/\beta$ enhancement means the mathematically
correct Born limit is not a complete physical threshold description. No
positronium production or particle condensation is simulated.

## Primary sources

- Particle Data Group, 2026 *Physical Constants*: supplied electron mass,
  fine-structure constant, and unit conversion.
  https://pdg.lbl.gov/2026/reviews/rpp2026-rev-phys-constants.pdf

- Vadim Kaplunovsky, *Annihilation*, UT Austin QFT notes, equations 5, 12–14
  and the scalar trace calculation:
  https://web2.ph.utexas.edu/~vadim/Classes/2024f-qft/annihilation.pdf
- Geant4 Physics Reference Manual, *Positron–Electron Annihilation*, target-rest
  total cross section per electron (Heitler):
  https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/electron_incident/annil.html
- A. Afanasev et al., *Quasi-Monoenergetic Photon Source Based on
  Electron-Positron In-Flight Annihilation*, IPAC2010, MOPEA043. Supports
  the invariant and identical-photon counting convention:
  https://proceedings.jacow.org/IPAC10/papers/mopea043.pdf

The derivation, CM reduction, elementary integration, and numerical oracle
are independently written for this milestone; source prose is not reproduced.
