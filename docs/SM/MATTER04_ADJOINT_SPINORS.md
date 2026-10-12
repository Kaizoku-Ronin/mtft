# Matter-04 — the adjoint question, answered through spinors (2026-10-02)

Module `mtft.matter.spin_twists`; tests `tests/test_matter04_spin_twists.py`. Classes as marked.

## The problem (EXACT within the BHV dictionary of Matter-02/03)
The topological twist turns the 8D gaugino into (0,p)-forms and the 8D scalar into a (2,0)-form, so every unbroken gauge generator carries
H^{0,*}(S, O) = (1, 14, 13): a gaugino and 27 adjoint chiral multiplets on X0(143) x 143a1 (h^{0,1} = 13 + 1, h^{0,2} = 13 x 1), unchanged by any
blow-up (birational invariants). For SU(3) with six flavours, N = 1 running gives b0 = 9 - 81 - 6 = -78: no asymptotic freedom while they are light.

## The spinor lever (EXACT)
With the R-symmetry bundle fixed, the fermions' square root of K can be changed by any 2-torsion line bundle: H^1(S, Z2) = J0(143)[2] x E[2] = Z2^28.
Fermions then live in Omega^{0,*}(zeta x eta); bosons, being tensors, do not see zeta or eta. Adjoint fermion zero modes = H^*(C, zeta) x H^*(E, eta):
- zeta = eta = O (BHV, supersymmetric): (1, 14, 13) — 28 per generator;
- zeta != O, eta = O: (0, 12, 12) — 24 (a nontrivial flat bundle on a genus-13 curve keeps g - 1 = 12);
- eta != O: **(0, 0, 0)** — on an elliptic curve a nontrivial flat bundle has no cohomology.
Only the genus-one factor removes every adjoint fermion, through any of its three **even spin structures**: the nontrivial 2-torsion points of 143a1,
roots of 4x^3 - 4x^2 - 4x - 7 (y = -1/2), irreducible over Q, discriminant 16 Delta = -2^4 * 11 * 13^2, Galois group S3, quadratic resolvent
Q(sqrt(-11)). This is a Scherk–Schwarz twist: the would-be preserved supercharge is antiperiodic around a cycle of 143a1.

## What survives
- Family chirality (EXACT): eta is torsion, so c1(L).(2R - K~) = e - 24b is unchanged; the three families from the exceptional flux remain.
- Charged fermions (EXACT): every block with E-degree b = 0 keeps only its exceptional jets, h^1 = (9q^2 + 3q)/2, h^0 = h^2 = 0 — quarks, leptons and
  higgsinos alike lose the genus-13 part. Matter-02 background, |Q| = 1, split level: 279 -> 153 per sign. The bosonic Higgs keeps every section (unlike a
  gauge Wilson line, the spin structure does not act on bosons).
- Bosonic adjoints (EXACT at tree level): 1 vector + 14 Wilson-line + 13 Phi complex scalars per generator, now without supersymmetric protection.
- One loop (STRUCTURAL, flat-circle form only): Scherk–Schwarz along a cycle is a thermal boundary condition; with matched spectra the Wilson-line potential
  keeps only odd windings, V = -2 sum_{n odd} B_n cos(n theta), minimised at theta = 0 with positive curvature. The remaining 26 bosonic adjoints are
  expected to get positive thermal-type masses — the calculation on the blow-up is OPEN. Below that scale QCD runs with b0 = 11 - 2*6/3 = 7 > 0.

## Flux consequence (EXACT given the twist)
Under the twist each unordered stack pair with exceptional difference e contributes e^2 fermion zero modes (|e| net, (e^2 - |e|)/2 vector-like pairs). Among
D-flat b = 0 assignments at eps/A_E = 1/3 with exact Standard-Model net indices, the total is minimised at 288 (phi_A1 in {0, 3}); the Matter-02 collinear
choice gives 1152. At the minimum the classes coincide pairwise (A1 with c, A2 with B): gauge Wilson lines on 143a1 must separate them, or the gauge group
enhances to U(4) x U(2) x U(2), of Pati–Salam type. Neither branch is tested yet.

## Costs and open gates
Supersymmetry is broken at the compactification scale by construction; the D-flat/BPS statements become tree-level approximations. The Higgs also receives a
positive one-loop term, so electroweak breaking must come from the tree-level slope (Fayet–Iliopoulos) term — the competition sets the scale (OPEN). The
index-forced jet pairs need couplings to become massive (OPEN). The full one-loop spectrum on the blow-up, including the 13 Phi and 13 curve Wilson-line
adjoints, is the next calculation.
