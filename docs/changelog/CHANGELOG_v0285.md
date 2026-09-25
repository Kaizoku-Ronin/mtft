# mtft v0.28.5 — WL-01: Wilson-line lifting verified against the Hodge metric (2026-09-13)

`surface.spectral.wilson_spectrum(N, theta, q)`: scalar Laplacian on the cusped Y0(N) twisted by the
flat U(1) connection with cycle coordinates theta (phases e^{i q zeta(e)}, zeta = B theta; flat at
cusps).  `wilson_line_mass_check`: ground state vs q^2 theta^T G theta / (56 pi - 4/Y0) with the
frozen Hodge metric — ratio 1.005 ± 0.001 across directions on X0(143), q^2 scaling 8.98/9.
Physics (overlay): the O-block vector-like pair of KK-01 is lifted by any nontrivial Wilson line
(cohomologically exact), with mass^2 given by the Hodge norm of the Wilson line; AL-preserving lines
live on the 143a1 block torus.  Protected chiral modes are not lifted.  Tests: 1 fast, 1 slow.
Base: live 0.28.3 (+ unpushed 0.28.4).  Pin four-way.
