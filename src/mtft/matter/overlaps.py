"""XQ-01 (2026-10-03): the overlaps X_qq' — family normalisation, the curve factor from X_0(143), and what they do to the texture.

1. Factorisation (EXACT in the PB-01 messenger model).  The messengers are curve copies: a holomorphic differential on X_0(143) times a massive torus mode.  The
   Yukawa diagram is H -> (D-term vertex, X-bar) -> <X X-bar> -> (family vertex, X), so the sum over curve factors is the Bergman kernel of X_0(143) and
   X_qq' = nu_q nu_q' x (Bergman jets of X_0(143) at the cusp) x (torus-propagator jets) x (condensate components).  All four cusps are Atkin–Lehner images of
   infinity, and Atkin–Lehner preserves the Bergman kernel, so the curve factor is the same at whichever cusp is blown up.
2. Family normalisation (EXACT).  The jets x^2, x w, w^2 have norms in ratio 2 : 1 : 2 on the exceptional curve, so canonically normalised families carry
   nu = (1, sqrt 2, 1) and the texture is nu H nu, not H.  Coherent and parity-even condensates then give m1 m3/m2^2 = 1/4, not 1 (CC-41).  The data sit at
   0.57, 2.55, 0.081 — still excluded (CP-8 survives) — and the rho targets become 0.568 (up), 0.684 (down), 0.408 (leptons).
3. The curve factor (NUMERICAL).  B_nm = a_n^dagger (Im Omega)^-1 a_m from the period-normalised q-expansions (n, m = 1..6, the Taylor orders n-1, m-1 in the
   cusp coordinate q).  Relative to B_11: diagonal 1, 2.094, 2.926, 4.160, 4.922, 5.208; off-diagonal correlations |B_nm| / sqrt(B_nn B_mm) all below 0.32
   (B_42 is 0.0018 of its diagonal scale).  On a disc the Bergman kernel is exactly diagonal in powers of q; the off-diagonal entries measure how much the
   global curve breaks the cusp's rotation.  So the curve factor nearly conserves the cusp charge.
4. Consequence under charge-conserving pairing (MODEL, conditional).  If the family vertex takes curve order 4 - n and the condensate vertex order n (each
   rotation charge conserved), the charge-n entry carries B_{5-n, n+1}: 0.316, 0.0054, 2.926 for n = 0, 1, 2.  The charge-1 channel is then suppressed about
   180-fold against its neighbours, the texture behaves as parity-even (m1 m3/m2^2 -> 1/4) whatever the tilt, and the data exclude it.  Either the family
   vertex pairs differently — which the harmonic representatives of the jet classes decide (CK-01 open item) — or the torus factor must supply
   t1^2/(t0 t2) of order 3e4, or the charge-1 entries arise elsewhere.
OPEN: harmonic jet representatives (the family-vertex pairing); torus-propagator jets; the condensate vertex beyond leading order."""
import math
import numpy as np

NU = np.array([1.0, math.sqrt(2.0), 1.0])


def cusp_bergman_jets(N=6, nmax=60, dps=30):
    """B_nm / B_11 for n, m = 1..N (real symmetric; imaginary parts vanish to working precision)."""
    from mtft.periods.channels import bergman_bilinear
    B = np.array([[complex(bergman_bilinear(n, m, nmax=nmax, dps=dps)) for m in range(1, N + 1)] for n in range(1, N + 1)])
    assert np.max(np.abs(B.imag)) < 1e-20
    return B.real / B.real[0, 0]

def correlations(B):
    d = np.sqrt(np.diag(B)); return B / np.outer(d, d)

def antidiagonal_weights(B):
    """Curve weights B_{5-n, n+1} for the charge-n entries (n = 0..4) under charge-conserving pairing."""
    return np.array([B[4 - n, n] for n in range(5)])

def texture_masses(a, b, normalised=True):
    H = np.array([[1.0, a, b], [a, b, 0.0], [b, 0.0, 0.0]])
    Y = np.diag(NU) @ H @ np.diag(NU) if normalised else H
    m = np.sort(np.linalg.svd(Y, compute_uv=False))[::-1]
    return m / m[0]

def geometric_ratio(a, b, normalised=True):
    m = texture_masses(a, b, normalised); return m[2] / m[1] ** 2

def data_rho(normalised=True):
    """First-branch rho from m1 m3/m2^2 = c (rho/(1-rho))^3, with c = 1/4 (normalised families) or 1 (plain Hankel)."""
    from mtft.matter.leptons import UP_MZ, DOWN_MZ, LEPTON_POLE
    c = 0.25 if normalised else 1.0; out = {}
    for k, m in (("up", UP_MZ), ("down", DOWN_MZ), ("leptons", LEPTON_POLE)):
        r = (m[0] * m[2] / m[1] ** 2 / c) ** (1 / 3); out[k] = r / (1 + r)
    return out

def report():
    B = cusp_bergman_jets(); C = correlations(B); w = antidiagonal_weights(B)
    off = np.abs(C - np.diag(np.diag(C)))
    return {"diagonal": np.diag(B), "max_offdiag_correlation": float(off.max()), "B42_correlation": float(abs(C[3, 1])),
            "charge_weights": w[:3], "charge1_suppression": float(math.sqrt(w[0] * w[2]) / abs(w[1])),
            "coherent_ratio_normalised": geometric_ratio(0.02, 0.0002, True), "coherent_ratio_plain": geometric_ratio(0.02, 0.0002, False),
            "rho_normalised": data_rho(True), "rho_plain": data_rho(False)}
