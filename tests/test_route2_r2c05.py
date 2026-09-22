import sympy as sp, pytest
from mtft.research import vector_higgs as VH, unified_parent as UP

def test_block_spectrum_polarisation_theorem():
    b6 = VH.block_spectrum(-6); assert b6["tachyon_m2"] == -sp.Rational(1, 4) and b6["tachyon_multiplicity"] == 18 and b6["massive_m2"] == sp.Rational(1, 4) and b6["massive_multiplicity_M1_M3"] == 7
    b3 = VH.block_spectrum(3); assert b3["tachyon_m2"] == -sp.Rational(1, 8) and b3["tachyon_multiplicity"] == 15 and b3["massive_multiplicity_M1_M3"] == 10 and not b3["flux_neutral"]
    assert VH.block_spectrum(0)["flux_neutral"] and VH.block_spectrum(0)["tachyon_m2"] == 0
    # Riemann–Roch cross-check of the multiplicities: h^1(L_neg) = h^0(K L_neg^-1) = (24 + |d|) - 12 ; massive: (24 - |d|) - 12 + h^0(O(|d| pts))
    for d, h0 in ((6, 1), (3, 1)):
        assert VH.block_spectrum(d)["tachyon_multiplicity"] == 24 + d - 12 and VH.block_spectrum(d)["massive_multiplicity_M1_M3"] == 24 - d - 12 + h0
    assert "-1/8" in UP.M3_RECORD["higgs"]

@pytest.mark.slow
def test_m3_block_on_mesh():
    r = VH.m3_block_check(0.3, 5); assert abs(r["tachyon_m2"] + 0.125) < 0.03 and abs(r["massive_m2"] - 0.125) < 0.04 and r["tachyon_gap"] > 0.2 and r["massive_gap"] > 0.05   # measured 0.35 and 0.09 (h = 0.3), 0.36 and 0.22 (h = 0.25)
