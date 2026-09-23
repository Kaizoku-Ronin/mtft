# R2C-05 — both polarisations of a flux block; the recombination scale (2026-09-19)

## CC-30 (my error in the M3 record)
The R2C-02 mass formula was derived for the (0,1)-polarisation of a NEGATIVE-degree block (M1's Higgs, d = -6).  Applied blindly to M3's
degree +3 block it gave -3/8, which violates the block symmetry (the (1,0)-polarisation of L is the conjugate of the (0,1)-polarisation of L^-1).
Enforcing that symmetry fixes the sign of the magnetic-moment term for each polarisation.

## Polarisation theorem (EXACT; mesh-verified at two degrees)
A gauge block with line bundle of degree d has two internal polarisations whose Landau levels sit at m^2 = +-d/(2g-2) = +-d/24 in
curvature units: the tachyonic one (Bochner on K (x) L_neg^-1, degree 24 + |d|) with multiplicity h^1(L_neg) = |d| + g - 1, and the massive
one (Bochner on K (x) L_neg, degree 24 - |d|) with multiplicity 12 - |d| + h^0(O(|d| points)).
   M1 (|d| = 6):  18 modes at -1/4 (measured -0.243) and 7 at +1/4 (measured +0.261; six tight, the h^0(O(2 sum P)) section shifted, then a gap)
   M3 (|d| = 3):  15 modes at -1/8 (measured -0.120, gap 0.36) and 10 at +1/8 (measured +0.137, gap 0.22) — h = 0.25 mesh
`vector_higgs.block_spectrum`; slow mesh test `test_route2_r2c05.py`.

## Consequence (EXACT structure, negative for the vector-Higgs hope)
Every flux-charged gauge block is tachyonic at the compactification scale, m^2 = -|d|/24 with the curvature radius as the unit: the split
flux background of any vector-Higgs parent is not a vacuum but recombines at M_KK.  M1's Higgs blocks (d = -6) and M3's 27 block (d = 3,
which contains the 10_H) are both flux-charged.  A light Higgs (m_H << M_KK) requires a flux-NEUTRAL Higgs block (d = 0: massless moduli with
no leading-order potential); E7 -> E6 x U(1) has no such 10.  Together with R2C-02 (the vector Higgs is the only option for M1's families,
y_t/g_4 ~ 0.6) this locates the electroweak hierarchy problem of the class exactly: the leading-order Higgs is either tachyonic at M_KK
(flux-charged) or an exact modulus (flux-neutral).  A TeV compactification scale does not rescue the charged case (m_H ~ M_KK/2).
Next: parents with a flux-neutral Higgs block (d = 0) and families in charged blocks — the block structure of E7/E8 decompositions with
TWO U(1) fluxes, searched with the same index arithmetic — and, for the neutral block, the mechanism that gives its moduli a small potential.


## Addendum (2026-09-23, v0.33.0): the h^0 values in the massive multiplicities
M3's 10 = 9 + h^0(O(3 CM points)) is proved (h^0 = 1: X0(143) has gonality >= 4, `hecke.gonality_lower_bound`).  M1's 7 = 6 + h^0(O(2 sum P)) uses
h^0(O(2 sum P)) = 1, which remains a register value supported by the mesh; an exact route is proposed in `V0330_COMPENDIUM_NOTES.md` item 1.
Later the same day: h^0(O(2 sum P)) = 1 is proved (`canonical.gates.gate_petri_w13_quotient`: the genus-6 quotient X0(143)/W13 has gonality >= 4 by Petri,
and h^0(O_X(2 sum P)) = h^0(O_Y(Q1+Q2+Q3))); M1's 7 = 6 + 1 is exact.  `V0330_COMPENDIUM_NOTES.md` item 1b.
