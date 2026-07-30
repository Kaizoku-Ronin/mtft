"""test_spectral_stage1_3.py — primon-gas spectral reconstruction (v0.10.0 stages 1–4).

Wraps each module's own certification suite as pytest, per the
Integration Plan §6 tiers:

- ``ledger.verify()``      — 46 checks, the certified constants and the
  RELATIONS identities (fast, < 1 s of ledger time; the mu checks build
  one N=400 f64 chain.internal)
- ``chain.selftest()``     — rung-4 kernel/gaps/dressed-B/limits (fast)
- ``expansion.selftest()`` — monomial extraction and closed forms (fast)
- ``coupled.selftest()``   — 15 checks: Bloch/Kesten measures, moments,
  tau_c scaling, binding thresholds (fast, ~2 s)
- ``ep.selftest()``        — 62 checks incl. the mp-tier winding
  certificates and the resolved S3-1 relation (slow, ~75 s)

Every asserted number comes from ``ledger.py`` or the closed forms —
never a literal (Integration Plan §6; Add. BI.F2), via the ``_L``/``_LF``
guards that raise KeyError on any unregistered name.  Two literals
remain in ``coupled.selftest`` pending its switch to ``_L`` (BN-F1;
the numbers are registered as ``tau_c_star`` / ``V_b_tree``, Add. BN §7).
``richardson()`` has no test surface here pending BI.F1 — the audited
Neville fix is verified (Add. BN §4); the re-export decision is the
author's.
"""
import pytest

from mtft import ledger, chain, expansion, coupled, ep


def test_ledger_verify():
    assert ledger.verify(verbose=False)


def test_chain_selftest():
    assert chain.selftest(verbose=False)


def test_expansion_selftest():
    assert expansion.selftest(verbose=False)


def test_coupled_selftest():
    assert coupled.selftest(verbose=False)


@pytest.mark.slow
def test_ep_selftest():
    assert ep.selftest(verbose=False)
