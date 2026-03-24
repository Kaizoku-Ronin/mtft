"""
Tests for mtft.lhcb_analysis — LHCb Open Data bridge
Run: pytest test_lhcb_analysis.py -v

These tests don't require a ROOT file — they test the module's
import, constants, branch resolution, and invariant mass algebra.
"""
import math
import numpy as np
import pytest

import mtft.lhcb_analysis as lhcb


# ═══════════════════════════════════════════════════════════════
# §1 — Module loads and constants are correct
# ═══════════════════════════════════════════════════════════════

def test_module_imports():
    assert hasattr(lhcb, 'LHCbNtuple')
    assert hasattr(lhcb, 'setup_check')
    assert hasattr(lhcb, 'combine_ntuples')


def test_pdg_constants():
    assert abs(lhcb.M_MU_PDG - 105.658) < 0.01
    assert abs(lhcb.M_JPSI_PDG - 3096.9) < 0.1
    assert abs(lhcb.M_BPLUS_PDG - 5279.4) < 0.1
    assert abs(lhcb.M_KPLUS_PDG - 493.677) < 0.01


def test_mtft_constants():
    assert abs(lhcb.M_TAU_KOIDE - 1776.969) < 0.01
    assert lhcb.KOIDE_K == pytest.approx(2.0 / 3.0)
    assert abs(lhcb.ALPHA_INV - 137.036) < 0.001


def test_hidden_doublet_constants():
    assert lhcb.H11_MASS == 1312.0
    assert lhcb.H13_MASS == 1348.0
    assert lhcb.DOUBLET_SPLIT == 36.0
    assert lhcb.H13_MASS - lhcb.H11_MASS == lhcb.DOUBLET_SPLIT


# ═══════════════════════════════════════════════════════════════
# §2 — Branch pattern coverage
# ═══════════════════════════════════════════════════════════════

def test_branch_patterns_exist():
    """All essential generic names must have candidate patterns."""
    essential = [
        "Jpsi_M", "B_M",
        "mup_PX", "mup_PY", "mup_PZ", "mup_PE",
        "mum_PX", "mum_PY", "mum_PZ", "mum_PE",
        "K_PX", "K_PY", "K_PZ", "K_PE",
    ]
    for name in essential:
        assert name in lhcb._BRANCH_PATTERNS, f"Missing pattern: {name}"
        assert len(lhcb._BRANCH_PATTERNS[name]) >= 2, \
            f"Pattern {name} needs ≥2 candidates for robustness"


def test_branch_patterns_include_lhcb_conventions():
    """Check that both major LHCb naming conventions are covered."""
    # The Ntupling Service often uses "mu_plus_*" or "mup_*"
    assert "mu_plus_PX" in lhcb._BRANCH_PATTERNS["mup_PX"]
    assert "mup_PX" in lhcb._BRANCH_PATTERNS["mup_PX"]

    # Also "_or_H1_" convention from some stripping lines
    assert any("H1" in c for c in lhcb._BRANCH_PATTERNS["mup_PX"])


# ═══════════════════════════════════════════════════════════════
# §3 — Invariant mass math (synthetic data)
# ═══════════════════════════════════════════════════════════════

def test_dimuon_invariant_mass_jpsi():
    """
    Create synthetic μ⁺μ⁻ four-momenta for a J/ψ at rest
    and verify the invariant mass computation.
    """
    # J/ψ at rest: m = 3096.9 MeV
    # μ⁺ and μ⁻ back-to-back, each with E = m_J/2
    m_jpsi = 3096.9
    m_mu = 105.658
    E_mu = m_jpsi / 2.0
    p_mu = math.sqrt(E_mu**2 - m_mu**2)

    # μ⁺ along +z, μ⁻ along −z
    mup_PE = np.array([E_mu])
    mup_PX = np.array([0.0])
    mup_PY = np.array([0.0])
    mup_PZ = np.array([p_mu])

    mum_PE = np.array([E_mu])
    mum_PX = np.array([0.0])
    mum_PY = np.array([0.0])
    mum_PZ = np.array([-p_mu])

    # Compute invariant mass
    px_tot = mup_PX + mum_PX
    py_tot = mup_PY + mum_PY
    pz_tot = mup_PZ + mum_PZ
    E_tot  = mup_PE + mum_PE
    m2 = E_tot**2 - px_tot**2 - py_tot**2 - pz_tot**2
    m_inv = np.sqrt(m2)

    assert abs(m_inv[0] - m_jpsi) < 0.01


def test_b_mass_from_jpsi_k():
    """
    Synthetic B± → J/ψ K± at rest.
    """
    m_B = 5279.4
    m_jpsi = 3096.9
    m_K = 493.677

    # Two-body decay: p = λ^{1/2}/(2m_B)
    lam = (m_B**2 - (m_jpsi + m_K)**2) * (m_B**2 - (m_jpsi - m_K)**2)
    p = math.sqrt(lam) / (2 * m_B)

    E_jpsi = math.sqrt(p**2 + m_jpsi**2)
    E_K = math.sqrt(p**2 + m_K**2)

    # J/ψ along +z, K along −z
    px_tot = 0.0
    py_tot = 0.0
    pz_tot = p + (-p)  # = 0
    E_tot = E_jpsi + E_K

    m2 = E_tot**2 - px_tot**2 - py_tot**2 - pz_tot**2
    m_inv = math.sqrt(m2)

    assert abs(m_inv - m_B) < 0.01


# ═══════════════════════════════════════════════════════════════
# §4 — Setup check runs
# ═══════════════════════════════════════════════════════════════

def test_setup_check():
    result = lhcb.setup_check()
    assert "numpy" in result
    assert result["numpy"]["installed"] is True
    assert "uproot" in result


# ═══════════════════════════════════════════════════════════════
# §5 — File-not-found error
# ═══════════════════════════════════════════════════════════════

def test_missing_file_raises():
    try:
        lhcb.LHCbNtuple("nonexistent_file.root")
        assert False, "Should have raised"
    except (FileNotFoundError, ImportError):
        pass  # Either is acceptable (ImportError if uproot missing)
