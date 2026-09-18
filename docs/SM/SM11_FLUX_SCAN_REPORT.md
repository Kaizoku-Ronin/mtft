# SM-11 — flux-divisor scan and the origin of the hierarchies (2026-09-14)

## Flux-divisor scan (untwisted family spaces, h = 0.2)
For every pure degree-3 divisor D supported on the 24 AL-fixed points (W13 triples, W143 triples, mixed), the
family Gram spectrum/max is (7.8e-5, m, 1) with the SAME smallest ratio 7.8e-5 to 1–3%, while the middle
eigenvalue m ranges over 6.8e-4 … 3.9e-2 (a factor 500).  One tested W143 triple, (0,1,2), is impure (rank 32:
u takes one sign on all three points) and is excluded.  Tool: `hym.family_gram_for_divisor`.
Interpretation (exact): every untwisted family space contains the two sections of the spin structure itself,
1 and u (pole-free); the universal 7.8e-5 is ‖1‖²/‖u‖² in the HYM metric of S0 — a number of the spin structure,
not of the flux.  The flux divisor moves only the third (pole-carrying) section.  Hence in any untwisted sector
the first-generation suppression cannot be tuned by the choice of CM points.

## Twisted sectors: h^0(S0 (x) t) = 0 for the quadratic, cubic and sextic twists (exact, 12 conditions of rank 12 on
the 12-dimensional weight-2 character spaces), so the twisted family spaces contain NO pole-free sections.
Two-mesh certification (h = 0.25 vs 0.2) of the twisted Gram spectra, agreement 1–2%:
   untwisted Q     (7.7e-5, 6.9e-4, 1)
   quadratic d^c   (1.9e-4, 1.1e-3, 1)
   cubic L         (8.1e-3, 1.1e-1, 1)
   sextic e^c      (2.6e-5, 3.0e-4, 1)
The hierarchy strength is a certified function of the twist character.  The cubic sector's mildness is real; the
quadratic and sextic sectors' strong hierarchies are real and NOT explained by pole-free sections — an open
structural question (the sections' growth toward the cusps 0 and 1/11, where the metric weight is concentrated,
is the obvious suspect).

## Consequence for model building
- Realistic first-generation quark masses (m_d/m_b ~ 1e-3) are not reachable in M1 by changing the flux divisor;
  they require the quark sectors to carry cubic-type twists (spread 1e2 → m1/m3 ~ 1e-5) or a different spin
  structure / parent.  A consistent all-cubic assignment exists: twists (L,c,a,b,d) = (0,2,4,4,2) in Z/6 make
  every bifundamental cubic-twisted and every Yukawa triangle character-neutral.  Its predicted pattern:
  m2/m3 ~ 1e-2, m1/m3 ~ 1e-5 in all sectors — right for the up quarks' first generation and the lepton second,
  low by 30–100 for the down/lepton first generation.  Computing it is the next run of the pipeline.
- The exact, direction-independent statements now certified: (i) 7.8e-5 is a spin-structure invariant;
  (ii) twist character ⇒ hierarchy class; (iii) the down–lepton transposition theorem for quadratic twists.
