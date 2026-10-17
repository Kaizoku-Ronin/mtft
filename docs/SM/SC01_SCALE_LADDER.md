# SC-01 — the scale ladder of the twisted compactification (2026-10-17)

Module `mtft.matter.scales`; tests `tests/test_sc01_scales.py`. Addresses puzzle III.5 ("assign M_KK and the Scherk–Schwarz scale"), which
Matter-04 and Matter-05 left open. The question has two parts: the ratios between the scales, and the absolute scale. The first is answered
here; the second is bounded on one side and waits on the volume (I.6).

## Two numbers

Everything on X0(143) x 143a1 with the even twist is set by the Kaluza–Klein scale M_KK = 2 pi / sqrt(A_E) and the dimensionless 8D coupling at
that scale, lambda = g_8^2 / A_E^2 = g_4^2 r, r = A_X / A_E.

## 1. The ladder (EXACT in PA-01's product approximation, under one convention)

In units of M_KK:

| scale | value | origin |
|---|---|---|
| messenger (Scherk–Schwarz) scale M_SS | 1/(2 sqrt(Im tau)) = 0.4943 | PA-01's M sqrt(A_E) = 3.106, now by a second route (momentum space, PG-01 2026-10-16) |
| messenger-pair threshold 2 M_SS | 0.9886 | equals the next boson level exactly |
| the saddles' gap | |tau|/(2 sqrt(Im tau)) = 0.5629 | both towers, both saddles (PG-01) |
| Wilson-line masses | 0.2175 sqrt(lambda), 0.4039 sqrt(lambda) | Matter-05's closed form, sqrt(96 T_eta / pi^6); ratio 1.857 = PG-6's metric curvatures |
| Phi-adjoint masses | 0.1324 sqrt(lambda) | Matter-05's closed form, sqrt(8 S_eta / pi^6) |

Convention C1: Matter-05's closed forms with n = 8 real bosonic degrees of freedom per charged direction and sum_alpha |alpha(e)|^2 = 8
(Tr e^2 = 1/2). The repository leaves the group factor open; a different choice rescales lambda by an O(1) number and nothing else.

## 2. The measured background (NUMERICAL)

FH-3 puts r = 230–280. GC-01 puts the Pati–Salam matching at 8.9e11–2.1e14 GeV; read as a lower bound on M_KK (R1: the breaking is a 4D
condensate below the KK scale). With g_4^2 = 4 pi alpha_3(M_KK) = 0.31–0.37 at those scales:

| | r = 230 | r = 280 |
|---|---|---|
| lambda | 72–85 | 88–103 |
| Wilson-line masses / M_KK | 1.84–2.00, 3.43–3.72 | 2.03–2.21, 3.78–4.10 |
| Phi masses / M_KK | 1.12–1.22 | 1.24–1.35 |

Every one-loop "modulus" is heavier than the KK scale. The light-modulus picture behind PA-01's straight path, PG-01's cosmological run and
the "ringing jet vacuum" question does not hold on the measured background: the Wilson line is a KK-scale state, not a light field.

## 3. Dissipation (EXACT given 2 and PG-01's thresholds; the width an ESTIMATE)

PG-01 (2026-10-16) found that at the jet vacuum the massless charged scalars couple to the Wilson line only quadratically, so the first
linear channel is a messenger pair at 2 M_SS = 0.9886 M_KK, reached when lambda > 20.7 (light mode) or 6.0 (heavy), while broad parametric
resonance needs lambda < 1.3. With lambda = 72–104 both modes sit above the pair threshold: the jet vacuum damps its Wilson line by linear
decay into messenger pairs; broad resonance and moduli trapping are excluded. Width: Gamma / omega = (N g_4^2 / 8 pi) beta^3 with beta = 0.84–0.89
and N the number of Dirac messengers at threshold (25 pairs per matter entry times their multiplicity, 25 <= N <= 300) gives 0.2–3. At best the
light mode is a broad resonance; at N in the hundreds it is not a particle.

## 4. The expansion parameters (STRUCTURAL)

Matter-05's closed forms are leading order in A_E / R_C^2 = 56 pi / r (the curve's Gauss–Bonnet area over r), which is 0.63–0.77 on the
measured background: not small. The 8D loop parameter at the KK scale, g_8^2 M_KK^4 / (4 pi)^4 = lambda / 16, is 4.5–6.5. The numbers of §2
are therefore O(1) uncertain, and the compactification is not weakly coupled at its own KK scale. What would fix the masses is the full heat
kernel on X0(143) — the sum over closed geodesics of Gamma_0(143), whose shortest has length 2 arccosh 2 = 2.634 (Matter-05) — in place of the
small-t expansion.

## 5. What fixes M_KK itself (STRUCTURAL)

A lower bound from GC-01's matching window; an upper bound would need the volume stabilised (I.6), which the programme does not have (PV-01's
runaway, MS-6's missing Casimir minimum, PG-4). The absolute scale stays unassigned; every ratio is fixed.

## Open

The group factor in C1; the geodesic-sum heat kernel for the one-loop masses at r = 230–280; the absolute M_KK (I.6). A reading R2 — M_KK equal
to the matching scale rather than above it — changes g_4^2 by less than 20% and none of the conclusions.

Sources: Matter-05 (MATTER05_ONE_LOOP_ADJOINTS.md), PA-01, GC-01, FH-01, PG-01 (section 6 of PG01_ENERGY_TWO_CLOCKS.md and the 2026-10-16 R&D notes).
