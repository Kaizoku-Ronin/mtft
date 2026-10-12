"""Matter-05: one-loop masses of the bosonic adjoints under the Scherk–Schwarz twist, and the Pati–Salam first look."""
import math
import numpy as np
from mtft.matter import one_loop_adjoints as OL, pati_salam as PSM, spin_twists as ST


def test_lattice_sums_split_one_plus_two():
    tab = OL.spin_structure_table(N=200)
    assert abs(tab[(0, 1)]["S"] - 2.107219) < 1e-5 and abs(tab[(1, 0)]["S"] - 3.238510) < 1e-5
    assert abs(tab[(1, 0)]["S"] - tab[(1, 1)]["S"]) < 1e-9                                      # the complex-conjugate pair agrees
    assert np.allclose(tab[(1, 0)]["T_eigenvalues"], tab[(1, 1)]["T_eigenvalues"], atol=1e-9)
    assert np.allclose(tab[(0, 1)]["T_eigenvalues"], (0.473605, 1.633614), atol=1e-5)
    assert all(min(v["T_eigenvalues"]) > 0 for v in tab.values())                                # positive-definite torus Wilson-line masses
    S1, _ = OL.lattice_sums((0, 1), N=150); S2, _ = OL.lattice_sums((0, 1), N=300); assert abs(S1 - S2) < 1e-7   # converged


def test_sign_theorem_ingredients():
    for t in (0.003, 0.03, 0.3, 3):                                                             # method of images: twisted trace below untwisted
        assert all(OL.circle_twisted_heat_trace(t, a) <= OL.circle_twisted_heat_trace(t, 0) + 1e-15 for a in (0.05, 0.2, 0.5))
    assert all(OL.torus_P(t) > 0 for t in (0.02, 0.1, 1, 10)) and abs(OL.torus_P(50) - 1) < 1e-6   # P(t) -> 1 (boson zero mode, none for fermions)
    e0, e1 = OL.minmax_shift(lambda x: 2 * (1 + math.cos(2 * math.pi * x)) ** 2); assert np.all(e1 >= e0 - 1e-9)   # positive potential raises every level
    assert OL.mass_coefficients()["K_theta"] / OL.mass_coefficients()["K_Phi"] == 12


def test_systole_and_windows():
    s = OL.systole_X0_143()
    (a, b), (c, d) = s["example"]
    assert s["trace"] == 4 and a * d - b * c == 1 and c % 143 == 0 and a + d == 4 and abs(s["length"] - 2 * math.acosh(2)) < 1e-12
    assert abs(s["area"] - 56 * math.pi) < 1e-9
    assert OL.qcd_b0_windows() == {"below_all_adjoints": 7.0, "curve_Wilson_lines_light_only": -6.0, "all_27_complex_light": -20.0}


def test_pati_salam_first_look():
    fl = PSM.first_look()
    assert fl["group_phi"] == {"4": 0, "2L": -3, "2R": 3} and fl["total"] == 288 == ST.ss_fermion_total(ST.MINIMAL_E, ST.RANKS)["fermion_zero_modes"]
    reps = {r["rep"]: r for r in fl["blocks"]}
    assert reps["(4,2,1)"]["net"] == 24 and reps["(4,2,1)"]["families"] == 3 and reps["(4bar,1,2)"]["families"] == 3
    assert reps["(1,2,2)"]["weyl_total"] == 144
    assert set(fl["right_handed"].values()) == {"d^c", "u^c", "e^c", "nu^c"}
