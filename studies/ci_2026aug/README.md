# ci_2026aug study bundle — the canonical-ideal arc (2026-08-16)

Wave artifacts from the four-session canonical-ideal arc (CI main /
AL-adapted / CI-AB / CI-C / CI-D), staged for v0.16.0 integration, plus
the auditor's correction and verification artifacts.

**Byte-preservation policy.** Every wave artifact listed under "Wave files"
is staged exactly as received — including statements the wave itself later
corrected (the C3 class-number miss). Corrections live in the append-only
`CC-0X_*.md` files and in the auditor artifacts; nothing upstream was
rewritten.

## The arc in one paragraph

X0(143) in its canonical embedding into P^12, degree 24, dim I_2 = 55
(Petri), quadrics generate; Atkin–Lehner adapted basis with W_11, W_13
simultaneously diagonal; I_2 graded into classes (26, 5, 4, 20); descent
of K_X to four eigen-line-bundles on E = X0(143)* = 143a1 of degrees
(0, 6, 5, 1); L_-- = O_E(P) with P = (4,-7) a generator of E(Q);
Jac(X0(143)/W_143) ~ 143a1 x 11a1; h(-143) = h(-572) = 10 (CC-09);
AL eigenvalues constant on Galois orbits (CC-08).

## Wave files (byte-preserved)

| file | content |
|---|---|
| `PREREG_canonical_ideal.md` | CI main pre-registration (P1–P9) |
| `PREREG_CI_AB.md` | CI-A/CI-B pre-registration |
| `PREREG_CI_C.md` | CI-C pre-registration (incl. C3, the honest miss) |
| `PREREG_CI_D.md` | CI-D pre-registration (D1–D3) |
| `X0_143_canonical_ideal_REPORT.md` | CI main report: dim I_2 = 55, Petri, generation |
| `X0_143_AL_adapted_ideal_REPORT.md` | adapted basis, grading, sparsity negative, projection table |
| `X0_143_CI_AB_REPORT.md` | descent structure, ten bundle-rank tests, CI-A/CI-B resolution |
| `X0_143_CI_C_REPORT.md` | ramification, (0,2,10), E = 143a1, CI-A forced, CC-09 proposed |
| `X0_143_CI_D_REPORT.md` | P = (4,-7), C2 equation, Route 2 (32/32), descent table |
| `X0_143_S2_qexpansions.txt` | PARI mfbasis f1..f13, q^0..q^140 |
| `X0_143_AL_adapted_basis.txt` | B (13x13), det B = -1078272 |
| `X0_143_AL_adapted_qexpansions.txt` | adapted basis e1..e13, q^0..q^140, integral |
| `X0_143_I2_quadric_basis.txt` | I_2 basis, 91 monomials x 55 quadrics (x-coords) |
| `X0_143_I2_by_AL_sector.txt` | I_2 split by AL class: 38x26, 11x5, 11x4, 31x20 (y-coords, LLL) |
| `X0_143_I2_sample_quadrics.txt` | sample quadrics |
| `X0_143_CI_A_quadric.txt` | Q*: a = -2439613813 = -7^2·13·1957^2 |
| `X0_143_canonical_ideal_cert.json` | wave certificate (updated with grading + projection table) |

SHA-1 manifest: see `sha1_manifest.txt` (generated at staging).

## Auditor artifacts

| file | content |
|---|---|
| `CC-08_AL_orbit_purity.md` | AG Pr 3.7.5 correction: AL constant on Galois orbits |
| `CC-09_class_number_143.md` | AG Pr 7.8.1 correction: h(-143) = h(-572) = 10 |
| `ci_verify_kimi.py` | self-contained exact replay of the whole arc (stdlib only) |
| `ci_verify_kimi.json` | its output |

## Audit verdict (2026-08-18, exact unless noted)

- QE = F·B on all 1833 entries; det B = -1078272.
- dim I_2 = 55 exact; kernel = span of the shipped basis; P8 rank 395 at
  two primes (Route A, prior session).
- e1 = 72·f_143a1 (140 coefficients); e8/e13 = g(q) ± 13·g(q^13),
  g = 11a1 = eta^2 eta(11.)^2 — exact.
- W_Q rebuilt independently (Cremona endpoint route on `mtft.hecke`,
  calibrated against the model's star involution): involutions, Klein
  group, [W_Q, T_p] = [W_Q, iota*] = 0, traces (2,-2,-18), eigenspaces
  (14,12)/(12,14)/(4,22), genera 7/6/2, X* = 1, block purity
  q4 = (-,+), q6 = (+,-) — **CC-08 confirmed computationally**.
- Ten bundle-rank tests 10/10; H^0(2K) grading 12/6/7/11; projection
  table 9/9 including all three honest excesses (newspace 3, f2+f3+old 2,
  f1+f2 1); CI-B ghost table 4/5/7/7; deficiency (0,0,3,0).
- Sector blocks: dims 26/5/4/20, support confined, max|coeff|
  55/4684/10008/72, all 55 quadrics residual exactly 0.
- CI-A: Q* residual 0; a = -7^2·13·1957^2 with 1957 = disc(g4); decoy
  20/20; a invariant under 6/6 random unimodular changes; "forced, not
  f2-specific" stands.
- CI-D: curve arithmetic exact over Q(i) (Q1+Q2 = 2P = -2G = (2,0);
  slope -3; h0(P) = 52; torsion trivial by gcd of point counts);
  disc -52 fixed-point identity exact; #Fix = (0,4,20) by both routes;
  cusp action (1 11)(13 143) / (1 13)(11 143), one free orbit;
  **Route 2 replayed 32/32 with both sides counted independently**.
  h(P) corroborated against the LMFDB regulator (Cert, numerical).
- Modular degree 4 / X* = 143a1 upgraded to an exact deduction:
  Hom(J0(143), E) = Z (certified charpoly factorization), deg pi = 4,
  deg phi_0 = 1 excluded by genus, so pi = ±phi_0 and deg phi_0 = 4.

## Honest miss recorded (wave's own)

C3 (h(-572) = 13) **failed** in the wave's CI-C session; tracing it
produced CC-09. Recorded here per the append-only policy.

## Internal → repo correction-label mapping

The wave's "CC-08" and "CC-09" proposals are adopted with the same ids
in the repo ledger (`src/mtft/hodge_polarization.py` docstring). No
renumbering.
