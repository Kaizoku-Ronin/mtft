# SM-12 — twist assignments: the all-cubic model and the Georgi–Jarlskog variant (2026-09-14; DIAGNOSTIC)

## Constraint structure (exact)
Stack twists t_x in Z/6 (classes: 0 trivial, 1 sextic chi_56, 2 cubic chi_133, 3 quadratic chi_13, 4 cubic-bar chi_100,
5 sextic-bar chi_23).  Bifundamental (x, y-bar) carries t_x − t_y; a Yukawa triangle is character-neutral automatically;
the shared H_d forces  t_Q + t_dc = t_L + t_ec  (mod 6).  Family hierarchy class by twist (certified, SM-11):
1e4 for {0, 3, 1/5}, 1e2 for {2, 4}.

## All-cubic model: twists (L, c, a, b, d) = (0, 2, 4, 4, 2) — every sector in class 2 (chi_133), Higgs in class 4 (chi_100)
chi_100 Higgs space from S_6(143, chi_100): 50 of 68 conditions, dim 18.  Tensor residual 1.3e-8
(chi_133^2 = chi_100 verified by the product test).  One tensor for all five sectors; up and down differ only in the
Higgs direction; down = lepton IDENTICALLY (same tensor, same H_d) — Georgi–Jarlskog fails as an identity.
   m2/m3 (5, 50, 95%) = 1.8e-2 / 7.2e-2 / 3.4e-1        m1/m3 = 6.4e-4 / 2.8e-3 / 1.4e-2
These ranges bracket the observed down (2e-2, 1e-3) and lepton (6e-2, 3e-4) ratios but are far too mild for the up
sector (m_u/m_t = 6e-6): the observed up hierarchy lies BETWEEN the two classes the twists generate (1e-8 and 1e-3).

## Georgi–Jarlskog variant: (Q, d^c) = (2, 2), (L, e^c) = (4, 0), shared chi_100 H_d
First attempt built the chi_100 family space by conjugating the chi_133 q-series; the product test rejected it
(1.7e-1) because conjugation reflects the CM points — a different flux divisor.  With the chi_100 space built from its
own Atkin–Lehner data (37 of 40, dim 3) the lepton tensor passes at 3.4e-8.  Results (same H_d directions):
   down:    m_s/m_b 1.8e-2 / 7.2e-2 / 3.4e-1      m_d/m_b 6.4e-4 / 2.8e-3 / 1.4e-2
   lepton:  m_mu/m_tau 4.4e-4 / 1.8e-3 / 1.8e-2   m_e/m_tau 8.5e-7 / 4.3e-6 / 2.6e-5
   (m_mu/m_tau)/(m_s/m_b) = 0.003 / 0.025 / 0.32 ;  (m_e/m_tau)/(m_d/m_b) = 1.6e-4 / 1.5e-3 / 1.5e-2   [observed 3, 1/3]
The mixed lepton pair splits leptons from down quarks in the WRONG direction (leptons lighter), by factors 40 and 700 at
the median.  In SM-10 the (cubic L, sextic e^c) pair split them in the right direction but by too much (65, 290).

## Reading
The twist classes control hierarchies by discrete jumps of ~1e2–1e4; observation asks for splittings of order 3.
Within the five-stack, three-W13-point, S0-based class on X0(143), the down–lepton relation is either an identity
(equal twists) or off by one or two orders in either direction (mixed twists).  Three exact ingredients would be
needed to change that: a different flux divisor for the lepton stacks (moves the middle family only, SM-11), a
different spin structure (changes the universal 7.8e-5), or lifting the Higgs-direction moduli (a potential).
The direction-independent statements stand; the specific Georgi–Jarlskog factor is not in reach of this class.
Status: DIAGNOSTIC (single mesh for the twisted Grams of this wave; Higgs Gram regularised).  Data: cubic_model.npz,
c100_family.npz, gj_variant_proper.npz (session files); the chi_100 weight-6 and weight-4 AL data are dumped for freezing.
