# CC-29 and R2C-04 — the corrected obstruction theorem and where a parent can live (2026-09-19)

## CC-29 (my error in R2C-03, found by re-deriving with R2C-01's conventions)
The "escape route (b)" flux (3, 0, 0, 0, -3) is not family-preserving: d_ca = +3 where R2C-01 needs -3, so u^c, d^c and L flip into u, d, L-bar.
Solving the six family targets gives (m_c, m_L, m_a, m_b, m_d) = (m_L + 3, m_L, m_L + 6, m_L + 6, m_L + 3): with no U(1)_L, (3, 0, 6, 6, 3),
i.e. M1 shifted by a constant, so every sector degree — the Higgs -6 included — is unchanged.  Removing U(1)_L cures only the f_L term;
the u^c and d^c sectors leave -f_a/2 and -f_b/2.  The test I had written asserted my own wrong expectation; it now asserts the solved flux.
Lesson: a test must check an independent derivation, not the author's claim.

## The theorem, corrected and stronger (EXACT, Hom-type content)
Every U(1)_x whose stack pairs with the colour stack in a family sector carries an irreducible U(1)_x–SU(3)^3 (2-form x 6-form) term sourced by
that family sector alone.  Hypercharge needs U(1)_a, U(1)_b, U(1)_d.  Therefore no Hom-type unitary-stack flux model with Standard-Model
bifundamental families is a consistent six-dimensional gauge theory.  Escape routes: non-Hom sectors (N_i, N_j) (three vector-like coloured
pairs per cancelled U(1)), or a parent whose simple factor has no cubic Casimir.

## R2C-04: admissible parents and the arithmetic families (EXACT representation facts)
- No cubic Casimir: SU(2), SO(N != 6), Sp, G2, F4, E6, E7, E8.  No independent quartic Casimir (6D gauge anomaly factorizes): SU(2), SU(3), G2, F4, E6, E7, E8.
- One complex block of U(1) charge q with flux degree d gives q d net families (index of S0 (x) L^q).
- E7 -> E6 x U(1): 133 = 78 + 1 + 27_{+1} + 27bar_{-1}; with L = O(P1+P2+P3) (degree 3, pure): THREE 27's of E6, each 16 + 10 + 1 of SO(10) —
  a family, a Higgs 10 and a singlet from one arithmetic bundle.  Record M3 (research record).
- E6 -> SO(10) x U(1): 16_{-3}: three 16's need a degree-1 flux O(P) (S0(3P) pure).
- E8: the family blocks carry weights of a traceless factor (SU(3), SU(4)); net families vanish on a curve.
- Gravitational anomaly: a pure gaugino theory needs dim(adj) = 0 mod 4 for tensor cancellation (E6, E7 fail; E8 passes); the (1,0) SUSY
  condition H - V + 29T = 273 gives H = 377 for E7 with one tensor.  OPEN for M3: the Higgs 10's mass (the 27 block's vector modes are
  tachyonic at the compactification scale, m^2 = -3/8), the gauge-anomaly spectrum with those hypers, vacuum, scale.

## Gravity in this class (EXACT reduction facts)
Kaluza–Klein reduction over a curve is a scalar–tensor theory with the area modulus as Brans–Dicke scalar, omega_BD = -(n-1)/n = -1/2: a
massless radion violates Cassini (omega > 4e4), so radius stabilisation is required by solar-system gravity.  1/g_4^2 = A/g_6^2: the internal
area plays the role of the vacuum permittivity of the gauge sector.  `compactification.radion_brans_dicke`.


## Addendum (2026-09-23, v0.33.0): CC-33
The pure-gaugino p2 rule becomes dim(adj) = 0 mod 28: E6, E7 and E8 all fail (the E8 survivor of R2C-04 is withdrawn).  `unified_parent.gaugino_p2_integrality`.
