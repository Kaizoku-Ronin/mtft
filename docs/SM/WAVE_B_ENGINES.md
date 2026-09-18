# Wave B — exact engines extracted from the handoff studies (2026-09-18)

Pure functions, no side effects, exact arithmetic; each test uses an independent construction of the expected answer.

| Module | Engine | Independent check in the test |
|---|---|---|
| `surface.spin_circle` | Seifert presentation, Fox-derivative cochains, rank-one cohomology with inverse-character duality, fiber holonomy, connection metric, 7D radius potential | fiber-acyclicity THEOREM vs the Fox computation (degrees +-3, +-6 -> [0,0,0,0]); Gysin controls [1,26,26,1] and [0,24,24,0]; Kaluza–Klein curvature -2/R^2 - r^2/(8R^4); volume 96 pi^2 R^2 r |
| `surface.hopf_geometry` | half-form pencil, Hopf pullback, line degrees | Riemann–Hurwitz: 24 = 12(-2) + 48; degrees +-12 = +-deg S0 |
| `surface.crt_dessin` | CRT projective line (78, 66), S/T permutations, cycles, genus formula | widths 1, 11, 13, 143 from T-cycles; 4 - 84 + 56 = 2 - 2g; agreement with the package's Manin-symbol permutations |
| `surface.hodge_blocks` | coprime block projectors P12/P14 (Bezout), quaternionic-multiplicity gate | charpoly x^2 (x+2)^4 g^2 s^2; ranks 12/14, idempotent, complementary (slow); (2,2,2,2) obstructed vs (4,4,4,4) |
| `research.anomalies` | cubic/gravity/mixed anomaly polynomials, polarised traces, mixed-anomaly matrix, primitive shift matrix, kernel | AXG-01 ideal factorization P = r1 Q1 + r2 Q2 and the mixed factorizations, reproduced from the charge sums; kernel = span{phase, Y, B-L} |
| `research.charge_lattices` | shift kernel, Smith remnants, integer dressing (Smith-form solve), kinetic-normalised vector masses | H_u H_d dressable, nu^c nu^c not; SNF (1,1); masses g^2 f^2 (4 +- 2 sqrt 2) |
| `research.mode_operators` | index vs cohomology, purity certificate, Bochner bound, 6D chirality gate | O(P) -> (2,1) vs O(P+Q-R) -> (1,0); M1 family h^0 = 3 pure; bound 1/4 = the magnetic-mesh Landau value |

One correction made on the way: my first transcription of the connection curvature used f = e/(2(2g-2)) instead of
f = 2 pi e / area = e/((2g-2) R^2); the Kaluza–Klein formula R_P = R_X - (r^2/2) f^2 then reproduces SC7-01's
-2/R^2 - r^2/(8R^4) — caught by the test, not by inspection.
