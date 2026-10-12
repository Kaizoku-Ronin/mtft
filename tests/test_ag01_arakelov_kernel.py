"""AG-01: the canonical Green's-function kernel and the pre-registered test."""
import math
from mtft.matter import arakelov_kernel as K


def test_construction_and_validation():
    assert K.multiplier_orders() == (-6.0, 6.0, 6.0, -6.0)
    assert abs(K.VALIDATION["ford_area"] / K.VALIDATION["target_area"] - 1) < 1e-5 and K.VALIDATION["gram_rel_err"] < 1e-5
    assert abs(K.VALIDATION["mu_cusp"] - 0.02434) < 1e-5 and all(K.spread(h) < 0.01 for h in K.KERNEL)
    assert 0.72 < K.ratio_to_flat("tau/2 preimage") < 0.73 and 0.85 < K.ratio_to_flat("Fricke point") < 0.86 and 0.59 < K.ratio_to_flat("cusp 0") < 0.60


def test_preregistered_outcome():
    lo, hi = K.r_range("Fricke point"); assert 165 < lo < 169 and 304 < hi < 308 and K.passes("Fricke point")
    assert not K.passes("tau/2 preimage") and K.touches("tau/2 preimage") and 68 < K.r_from_vcb("tau/2 preimage", 0.041) < 70
    assert all(not K.passes(h) and not K.touches(h) for h in ("cusp 0", "cusp 1/11", "cusp 1/13"))
    assert [h for h in K.ARITHMETIC if K.passes(h)] == ["Fricke point"]
