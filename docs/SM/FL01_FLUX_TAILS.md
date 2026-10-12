# FL-01 — which tail sets which angle, and what the data demand (2026-10-05)

Module `mtft.matter.flux_tails`; tests `tests/test_fl01_flux_tails.py`.

1. **Which tail sets which angle (EXACT).** In Y_ij = int Q_i d^c_j H each off-diagonal entry is dominated by the degree-0 field's tail reaching the jet point: Y_23 ~ T(d^c),
   Y_32 ~ T(Q). The CKM angle theta_L is set by the right-handed block's tail; theta_R by Q's, shared by up and down. TO-01 item 4 had the direction backwards (CC-51).
2. **What the data demand (EXACT given the induced texture).** m_c/m_t = theta_L^u theta_R, m_s/m_b = theta_L^d theta_R, V_cb = |theta_L^d - e^{i phi} theta_L^u|: the tails must stand
   T(Q) : T(d^c) : T(u^c) = 40 : 5.3 : 1 (aligned) to 86 : 5.3 : 1 (opposite phases). Q longest, u^c shortest.
3. **Curve flux cannot supply it (ESTIMATE).** Tails fall like exp(-d sqrt(2 pi a / A_X)) beyond the magnetic length; for a separation along the curve the ratio is independent of
   A_X and reaches at most 13.8 under a_Q + a_u = 24 — and stability's a_Q > a_u points it the wrong way.
4. **The remaining source (STRUCTURAL).** The condensate lives in the right-handed block; Pati–Salam-scale effects act on u^c and d^c but not Q, and its SU(2)_R D-term splits them.

Open: right-handed tail ranges with the condensate's D-terms; the sign of the u^c–d^c splitting.
