import numpy as np
from mtft.surface import rrspace as RR

def test_cubic_higgs_spaces_and_mixed_tensors():
    cm = RR.cm_classes()
    for a in (100, 133):
        H = RR.cubic_higgs_space(a, cm); assert H["rank"] == 50 and H["nullspace"].shape == (68, 18)
    sp = RR.cubic_family_space_100(cm); assert sp["rank"] == 37 and sp["nullspace"].shape == (40, 3)
    up = RR.mixed_yukawa(133, 0, 133, cm); assert up["residual"] < 1e-6            # Q cubic x u^c untwisted -> chi_133 Higgs
    dn = RR.mixed_yukawa(133, 133, 100, cm); assert dn["residual"] < 1e-6          # all-cubic down/lepton -> chi_100 Higgs
    bad = RR.mixed_yukawa(133, 0, 100, cm); assert bad["residual"] > 1e-2          # character-violating pair is rejected by the product test
