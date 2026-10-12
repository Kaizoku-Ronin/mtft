# PB-01 — messenger-exchange Yukawas: the holomorphic route closes, the non-holomorphic one opens (2026-10-02)

Module `mtft.matter.messenger_yukawa`; tests `tests/test_pb01_messenger_yukawa.py`.

1. **Holomorphic exchange gives nothing (EXACT).** Heavy-fermion exchange between superpotential vertices yields Q u^c H times holomorphic vevs. One condensate
   insertion is never SU(4)-invariant (4 x 4bar x 4 has no singlet); Sigma^n needs n = 0 mod 4 and an epsilon on at least three identical bosons, which forces
   antisymmetry in three SU(2)_R doublet indices — zero; Sigma Sigmabar is invariant, but the D-term vacuum has one orientation only (the D-terms fix the
   difference of the two vevs, the positive one-loop masses set the conjugate one to zero; checked analytically and by scan). With YK-01 this closes every
   tree-level superpotential source of jet-family masses.
2. **The non-holomorphic route (STRUCTURAL).** The one-orientation vacuum leaves nonabelian D-terms <D^a> ~ x_* (T^a)_{44}, x_* = (g_U^2 xi - m^2)/(g_U^2 + C g^2).
   They mix the light Higgs with every heavy bidoublet boson X: theta_X = g^2 <D^a> c_X / M_X^2. Heavy bosons are not closed forms, so jet–jet–X couplings are not
   cup products: Y_ij = sum_X theta_X int psi_{Q,i} ^ psi_{u,j} ^ X — diagonal entries from the heavy modes' values at cusp i, off-diagonal from the tails of the
   harmonic jet representatives, size y ~ g^2 y_X v_R^2 / M_X^2 (an order-one top needs v_R near the first heavy bidoublet scale).
3. **Arithmetic skeleton (NUMERICAL).** Cusp values of weight-2 forms are AL eigenvalues times a_1, so the Bergman kernel of S_2(Gamma_0(143)) at the cusps is a
   Klein-group circulant with sector weights N_{++} = 2.485, N_{+-} = 4.438, N_{-+} = 4.502, N_{--} = 1.066 (sizes 1, 6, 5, 1; Petersson Gram from the frozen periods).
   Eigenvalues 4 N_chi on four cusps; on any three cusps 6.44, 13.16, 17.88 (0.36 : 0.74 : 1): no large hierarchy here — it must come from theta_X and the Higgs
   alignment.
OPEN: the first heavy bidoublet modes on X0(143) and their cusp values; the overlaps c_X; harmonic jet representatives and their tails (CKM).
