# Wave C — model records and the gate battery (2026-09-18); Gate 4 decision: route 2

Decision (R. Tano, 2026-09-18): the MTFT-native target is ONE parent with THREE internal modes; C3X's consistency
machinery (anomaly polynomial, charge lattice, discrete anomaly, tensor lattice, bordism) becomes the set of gates any
such parent must pass.  C3X remains a separate record (three input copies); it is a benchmark, not the target.

## Modules (`mtft.research`)
| Module | Content | Independent check in the tests |
|---|---|---|
| `parents` | immutable records M1, ONE_PARENT_INDEX3 (control), U8_ADJOINT (control), C3X, with conventions and multiplicity origin | immutability; records distinct |
| `discrete_anomalies` | Spin x Z_n fermion test (Hsieh), M1 Z3 ledger with generator g = (0,0,1,0,1), generator certificate | S3 = 24 recomputed from the anomaly polynomial on the generator charges; residue 12 mod 18 (FAIL) |
| `tensor_gs` | native-scalar flux gate (K m, k_Delta m), flux transgression (K_hat, SNF (1,3)), factorization check, C3X I8, integral lattice gate | K m = (12,12); k_Delta m = 0; SNF by Smith form; I8 = X4 Y4 symbolically; Omega even unimodular of signature (1,1) |
| `compactification` | 6D/7D restricted radius potentials, unwarped product background, C3X AdS control, radion mass, Planck reduction | lambda4, k reproduced by the direct Einstein solve; ell4^2 < 3R^2; control radion mass^2 = 8 |
| `bordism` | recorded Sq^2 matrices (F2) between H^4, H^6, H^8 of BG; ranks, composite, kernel = image, E3^{6,1} = 0 | F2 rank recomputation: (2, 7), kernel dim 2 = image dim 2 |
| `pipeline` | `m1_gate_report`, `c3x_gate_report`, `route_2_requirements` — every gate returns pass/None/False WITH a witness and a status; no "viable" boolean | reports carry witnesses; M1's failing gates listed |

## Where M1 stands against the route-2 gates (today)
PASS: 4D local anomaly (Y, B-L, phase exact); charge lattice SNF (1,1); native-scalar flux (k_Delta m = 0); three
internal modes (h^0 = 3, pure).  FAIL / OPEN: Spin x Z3 fermion test (residue 12 mod 18, fermion-only); 6D anomaly
lift (irreducible p2, tensor signature (2,1) — AXG-02/03); 6D Yukawa chirality (no 6D chiralities assigned; same-
chirality bilinear obstruction); classical radius (runaway); bare Majorana forbidden (B-L must break).
That list is the route-2 research programme: a parent whose 6D chirality assignment, tensor sector and discrete
remnants pass the gates C3X passes, while keeping the three families as internal zero modes.

## Where C3X stands (benchmark)
PASS: local factorization; integral even unimodular lattice; ordinary spin bordism; one pure mode per copy.
OPEN: global GS construction; families (input); background (AdS control, no scale separation); scale.
