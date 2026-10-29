"""CC-71 (2026-10-29): the Tano metric of a diatomic gas divides S°(X2) by the mass of X2, not by the mass of one atom.

Second route: the entropies of H2, N2 and O2 counted from first principles (Sackur–Tetrode, the quantum rotor with nuclear-spin statistics, the
harmonic oscillator, the spin multiplicity) reproduce the stored S° per mole of molecules, which is what makes the division by 2M the right one.
H2 needs Pauli: even J pair with the one antisymmetric nuclear-spin state and odd J with the three symmetric ones; 16O has spin 0, so the O2
ground state 3Sigma_g^- keeps only odd rotational levels."""
import math

import pytest

from mtft import tano_metric as TM

R, KB, H, NA, T, P, C2 = 8.314462618, 1.380649e-23, 6.62607015e-34, 6.02214076e23, 298.15, 1e5, 1.438776877  # SI; C2 = hc/k in cm K
H2_LEVELS = [0.0, 118.487, 354.373, 705.519, 1168.798, 1740.189, 2414.769, 3187.57, 4051.73, 5001.97]  # H2, v = 0, J = 0..9, cm^-1
DIATOMIC = ("H", "N", "O")


def s_trans(m_molar):
    m = m_molar / 1000 / NA
    return R * (math.log((2 * math.pi * m * KB * T / H ** 2) ** 1.5 * KB * T / P) + 2.5)


def s_levels(levels):
    z = sum(g * math.exp(-C2 * e / T) for e, g in levels)
    u = sum(g * C2 * e / T * math.exp(-C2 * e / T) for e, g in levels) / z
    return R * (math.log(z) + u)


def s_vib(we):
    x = C2 * we / T
    return R * (x / math.expm1(x) - math.log(-math.expm1(-x)))


def counted(symbol):
    """S° per mole of X2 by counting states (the 'practical' convention: nuclear-spin entropy left out)."""
    el = TM.get_element(symbol)
    if symbol == "H":
        rot = s_levels([(e, (2 * j + 1) * (3 if j % 2 else 1)) for j, e in enumerate(H2_LEVELS)]) - R * math.log(4)
        return s_trans(2 * el.molar_mass) + rot
    if symbol == "N":   # 14N: I = 1, bosons; even J take the 6 symmetric spin states, odd J the 3 antisymmetric ones
        rot = s_levels([(1.99824 * j * (j + 1), (2 * j + 1) * (3 if j % 2 else 6)) for j in range(120)]) - R * math.log(9)
        return s_trans(2 * el.molar_mass) + rot + s_vib(2358.57)
    rot = s_levels([(1.43768 * j * (j + 1), 2 * j + 1) for j in range(1, 160, 2)])   # 16O2: only odd J
    return s_trans(2 * el.molar_mass) + rot + s_vib(1580.19) + R * math.log(3)       # R ln 3: the electron-spin triplet


def test_diatomic_rows_are_per_molecule():
    assert {s for s, e in TM.ELEMENTS.items() if e.atoms == 2} == set(DIATOMIC)
    assert all(e.atoms == 1 for s, e in TM.ELEMENTS.items() if s not in DIATOMIC)


@pytest.mark.parametrize("symbol", DIATOMIC)
def test_tano_metric_divides_by_the_mass_of_the_molecule(symbol):
    el = TM.get_element(symbol)
    assert el.tano_metric == pytest.approx(el.molar_entropy / (2 * el.molar_mass), rel=1e-15)
    assert TM.TANO_METRIC_RETRACTED_CC71[symbol] == pytest.approx(2 * el.tano_metric, rel=1e-15)   # the retracted value was twice too large


def test_monatomic_rows_unchanged():
    for s, e in TM.ELEMENTS.items():
        if s not in DIATOMIC:
            assert e.tano_metric == e.molar_entropy / e.molar_mass


@pytest.mark.parametrize("symbol, tol", [("H", 0.01), ("N", 0.1), ("O", 0.05)])
def test_counted_entropy_is_the_stored_entropy_per_molecule(symbol, tol):
    el = TM.get_element(symbol)
    assert counted(symbol) == pytest.approx(el.molar_entropy, abs=tol)                  # per mole of X2: 130.68, 191.61, 205.14
    assert abs(counted(symbol) / 2 - el.molar_entropy) > 50                              # non-vacuity: per mole of atoms it would be far off


def test_pauli_in_the_hydrogen_row():
    no_pauli = s_trans(2 * TM.get_element("H").molar_mass) + s_levels([(e, 2 * j + 1) for j, e in enumerate(H2_LEVELS)])
    gap = no_pauli - counted("H")
    assert 5.6 < gap < R * math.log(2)                                                   # 5.71: close to R ln 2 = 5.76


def test_seebeck_diagnostic_uses_the_same_normalisation():
    for s in ("Au", "Ni", "O"):
        el = TM.get_element(s)
        psi = TM.seebeck_diagnostic(s, 20.0)["psi"]
        m_kg = el.molar_mass * 1e-3
        assert psi == pytest.approx(20e-6 * m_kg / (el.molar_entropy / (el.atoms * m_kg) * TM.E_CHARGE_SI), rel=1e-12)


def test_default_screen_moves_with_the_correction():
    rows = {(r["light"], r["heavy"]): r for r in TM.materials_screening()}
    la = rows[("H", "La")]
    assert la["delta_T"] == pytest.approx(TM.get_element("H").tano_metric - TM.get_element("La").tano_metric, rel=1e-12)
    assert la["lambda_eff"] == pytest.approx(1.4441, abs=1e-4) and la["T_c_K"] == pytest.approx(750.5, abs=0.1)   # 930.1 K before CC-71
    assert rows[("O", "Ba")]["T_c_K"] == pytest.approx(109.8, abs=0.1)                                             # 125.6 K before CC-71
    assert rows[("B", "Mg")]["T_c_K"] == pytest.approx(117.9, abs=0.1)                                             # no diatomic: unchanged
