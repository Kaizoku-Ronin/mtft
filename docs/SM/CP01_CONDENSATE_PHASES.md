# CP-01 — the condensate's jet alignment: spin-1 phases and a Georgi–Jarlskog structure (2026-10-03)

Module `mtft.matter.condensate_phases`; tests `tests/test_cp01_condensate_phases.py`; command `python -m mtft matter phases`.

1. **Spin 1 (EXACT).** The blown-up point is locally U(2)-symmetric; its degree-2 jets x^2, x w, w^2 form the spin-1 representation. An SU(2)-invariant quartic has two
   invariants, |c|^4 and |theta|^2 with theta = c_0^2 - 2 c_+ c_-, hence two phases: coherent (a perfect square (a x + b w)^2) and polar (orthogonal lines).
2. **Tree level prefers polar (EXACT).** The D-term quartic is the fourth moment of the profile: E|s|^4 = (1/5)|c|^4 - (1/15)|theta|^2, normalised from 1.8 (coherent, maximum)
   to 1.2 (polar, minimum).
3. **Textures (EXACT within the Hankel class).** m2/m3 = |D2 - D1^2|, m1/m3 = |D2|^3/m2. Polar along x w: rank one; polar along x^2 - w^2: no hierarchy; tilted coherent:
   D2 = D1^2/2, the quark-like hierarchy. The down-quark fit has rho = D2/D1^2 = 0.59, near coherent — so the observed vacuum needs the spin-dependent quartic to change
   sign relative to tree level, a definite target for the one-loop effective potential.
4. **A Georgi–Jarlskog structure (STRUCTURAL; hypothesis).** If leptons see the charge-2 coupling with the opposite sign: m_tau = m_b, and m_mu/m_s = F, m_e/m_d = 1/F with
   F = (1 + rho)/(1 - rho), so m_mu m_e = m_s m_d exactly (data 0.90); rho = 1/2 gives 3 and 1/3. From the down fit: m_mu/m_s = 3.57 (data 3.06), m_e/m_d = 0.282 (0.295).
   Spotted after fitting, so a hypothesis to derive; a reflection 1 - 2 P_lepton through the condensate's direction on lepton lines would supply the sign.
5. **[Retracted by CC-40] Where each sector sits (MODEL).** |theta|/|c|^2 for the equal-overlap fits: up 0.006, down 0.033 (coherent), leptons 0.926 (polar). Coherent means D2 = +D1^2/2, polar-leaning
   D2 < 0: the sign flip of item 4 is the same distinction. Leptons — the sector the condensate points along — follow the tree-level (polar) phase; what makes quarks coherent is open.

## CP-02 (2026-10-03): zero flips, parity, and the sharpened one-loop target

6. **Zero flips (EXACT; CC-40).** The texture sees the profile only through |p|^2 on the unit circle; reflecting a zero through the circle leaves it unchanged. The lepton
   fit (1, 3.26, -0.31) and its flip (1, 0.205, -0.028) give identical masses, so item 5 is an artefact and is retracted. The invariants are D1 and rho; two mass ratios
   fix rho only up to two branches per sector: up 0.453 or -4.80, down 0.577 or 3.73, leptons 0.302 or -0.765. CP-4's sign flip is one branch choice.
7. **Parity (EXACT in the model).** In the Landau problem — twist anisotropy against the tree-level quartic — the ground state is x^2 (kappa >= 2) or an x^2/w^2 mixture,
   always with no x w component. Then D1 = 0 and the light masses are geometric, m1 m3 = m2^2.
8. **Data (EXACT given the texture class).** m1 m3/m2^2 = 0.57 (up), 2.55 (down), 0.081 (leptons): not geometric. The families need the charge-1 component — w -> -w broken
   in the condensate's alignment.
9. **The one-loop target (STRUCTURAL).** In the ratio-3 background the condensate entry is gamma-neutral; parity-odd terms in its potential come only from loops of
   gamma-charged fields, and they set D1.

## CP-03 (2026-10-03): the parity-odd term acts as a transverse field

10. **Zeeman structure (EXACT).** The field strength at the point splits into a trace and a spin-1 vector; by Wigner–Eckart the jets feel only <c|B.J|c>. Curve/torus
    flux differences give the longitudinal (P-even) part; the mixed gamma-flux gives the transverse (P-odd) part.
11. **Coherent vacuum (EXACT).** The minimum of -B.J is the coherent state along B — a tilted perfect square — so rho = 1/2 exactly for every tilt (D1 = 0.202 at 0.29 rad).
12. **Not enough (EXACT within the Hankel class).** rho = 1/2 gives m1 m3 = m2^2, exactly as D1 = 0 does: parity breaking is necessary (CP-8) but not sufficient. The data need
    rho - 1/2 = -0.047 (up) and +0.077 (down); these must come from unequal overlaps X_qq', which the equal-overlap texture sets to 1.
13. **The quartic does not help (EXACT).** The coherent vacuum stays exact up to lambda = 1/2, then rho runs away (-5.3 at 0.55); small departures need fine tuning.
14. **Size (STRUCTURAL).** Tree-level for a gamma-charged condensate (EPS-01's estimate is this tilt); one-loop for the gamma-neutral ratio-3 condensate (CP-9), needing a large trace.
