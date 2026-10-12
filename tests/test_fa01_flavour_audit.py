"""FA-01: the flavour audit."""
from mtft.matter import flavour_audit as F


def test_images_close_the_r_gap():
    assert 26 < F.r_single("over 1/2", 0.041) < 27.5 and 27 < F.r_single("over tau/2", 0.041) < 28.5
    a = [F.r_images("over 1/2", v) for v in F.VCB]; b = [F.r_images("over tau/2", v) for v in F.VCB]; f = [F.r_images("Fricke", v) for v in F.VCB]
    assert 110 < a[0] < 120 and 60 < a[2] < 65 and 125 < b[0] < 135 and 69 < b[2] < 72 and 190 < f[2] < 200 and 350 < f[0] < 365
    assert F.overlap((a[2], a[0])) and F.overlap((b[2], b[0])) and F.overlap((f[2], f[0]))
    assert not F.overlap((F.r_images("opposite cusp", 0.055), F.r_images("opposite cusp", 0.030)))
    assert abs(F.vcb_images(2.5, "1/2") / F.vcb_2d(2.5) - 1) < 1e-3


def test_curvature_and_arithmetic():
    assert 12.27 < F.curvature_radii(1.0) < 12.29 and 2.6 < F.curvature_radii(0.214) < 2.65 and 1.3 < F.curvature_radii(0.107) < 1.32
    assert F.cusps_rational() == (4, True) and F.no_rational_isogeny() and not F.no_rational_isogeny(163)
    assert len(F.SCOREBOARD) == 8 and sum(1 for s in F.SCOREBOARD if s[1].startswith("failed")) == 2
