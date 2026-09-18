# KK-TOWER-02 — mesh certification of the towers; Yukawa tensor from FEM eigenmodes; KK decay overlaps (2026-09-14)

## Certification (three meshes: h = 0.25, 0.2, 0.15)
   family tower first excitation above the zero-mode level:  0.0175 / 0.0227 / 0.0232  -> Delta lambda_1 = 0.023 (two finer meshes agree to 2%)
   generic degree-15 control (h = 0.2): 0.066                -> the light state is a property of the M1 bundle S0(sum P)
   S0 first excitation: 0.052 / 0.054 (stable)
   Higgs gap above the 18-fold Landau level: 0.590 / 0.552 / 0.544 -> 0.55 +- 0.03
   The "19th Higgs mode" of KK-TOWER-01 (+0.05 at h = 0.2 with nev = 28) is absent at h = 0.25, h = 0.15 and at h = 0.2
   with nev = 20: an eigensolver/mesh artefact — RETRACTED.
   Landau-level offsets converge: family +1.1% / +0.7% / +0.4%; Higgs +0.9% / +0.6% / +0.3%.

## Yukawa tensor from FEM eigenmodes (EXACT construction, DIAGNOSTIC numbers)
Kinetic-orthonormal eigenmodes need no Gram matrix: Y_ijk = int psi_i psi_j conj(psi^H_k) dA_c.  The Higgs connection is
exactly twice the family connection on every edge (K(2 sum P) = S0(sum P)^2), so the integrand is gauge-invariant nodewise;
the tensor is symmetric to 1e-17 and the zero-mode Gram is the identity to 1e-12.
   M1 up ratios, FEM pipeline (5/50/95%):   m2/m3 0.329 / 0.606 / 0.874    m1/m3 0.065 / 0.255 / 0.533
   CC-26-corrected modular-form pipeline:   m2/m3 0.342 / 0.615 / 0.868    m1/m3 0.074 / 0.258 / 0.530
   retracted (transposed contraction):      m2/m3 ~ 1.7e-4                 m1/m3 ~ 1.4e-8
Two pipelines sharing only the curve agree to 1–3% and exclude the retracted values by 6–8 orders.  This is the
independent confirmation of CC-26 and of the anarchic leading order.  Gate: `test_fem_yukawa_reproduces_corrected_normalisation`.

## KK decay overlaps (DIAGNOSTIC)
   KK1 (Delta lambda 0.023) -> f_i + H_k: rms overlap 0.020 vs zero-mode Yukawa rms 0.022 (ratio 0.92)
   KK2 (Delta lambda 0.059):               rms overlap 0.017
The light vector-like family state couples to fermion + Higgs at full Yukawa strength: its decay is unsuppressed, and its
mass is M_1 = sqrt(0.023) ~ 0.15 in curvature units (the compact metric has K = -1, area 48 pi).  Whether that is a
collider prediction or a GUT-scale state depends on the curvature radius, which nothing in the present class fixes.

## Next (dependency order)
Overlaps of KK modes with zero modes are now available from `magnetic.m1_fem_yukawa`; the one-loop Coleman–Weinberg
potential on the Higgs moduli space and the family wavefunction renormalisation are the next two computations.
