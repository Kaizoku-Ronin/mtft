"""
tests/test_music.py
===================
First test suite for mtft.music (new in v0.7.0): canonical constants,
scale construction, deterministic composition, and WAV rendering.
"""
import wave

import pytest

from mtft.constants import CriticalDepths
from mtft.music import (
    F_BASE,
    Y_C,
    ModularScale,
    MonsterComposer,
    VacuumSonifier,
    spectral_fingerprint,
)


# ── canonical constants ──────────────────────────────────────────

def test_yc_is_canonical():
    assert Y_C == CriticalDepths.y_conf == 0.18174


def test_f_base_is_the_level():
    assert F_BASE == 143.0


# ── scales ───────────────────────────────────────────────────────

def test_supersingular_scale():
    ss = ModularScale().supersingular_scale()
    # the 15 supersingular primes of the Monster
    assert sorted(ss) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
    assert ss[2] == pytest.approx(160.061181, rel=1e-9)
    assert all(f > 0 for f in ss.values())


def test_koide_triad_is_augmented():
    t = ModularScale().koide_triad()
    assert t[0] == pytest.approx(F_BASE)
    assert t[1] / t[0] == pytest.approx(2 ** (1 / 3))
    assert t[2] / t[0] == pytest.approx(2 ** (2 / 3))
    assert t == pytest.approx((143.0, 180.16871013, 226.99835043), rel=1e-9)


def test_spectral_fingerprint_deterministic():
    fp1 = spectral_fingerprint(n_max=50)
    fp2 = spectral_fingerprint(n_max=50)
    assert set(fp1) == {"brightness", "consonance", "spectrum", "top_13"}
    assert fp1["brightness"] == pytest.approx(fp2["brightness"])
    assert fp1["consonance"] == pytest.approx(fp2["consonance"])
    assert list(fp1["top_13"]) == list(fp2["top_13"])


# ── composition ──────────────────────────────────────────────────

def test_generate_phrase_deterministic():
    mc = MonsterComposer()
    p1 = mc.generate_phrase(0, "electron")
    p2 = mc.generate_phrase(0, "electron")
    key = lambda ph: [(n.frequency, n.duration, n.velocity) for n in ph.notes]
    assert key(p1) == key(p2)
    assert len(p1.notes) == 13  # genus-length phrase
    assert all(n.voice == "electron" for n in p1.notes)


# ── WAV rendering ────────────────────────────────────────────────

def test_render_depth_tone_writes_valid_wav(tmp_path):
    out = tmp_path / "tone.wav"
    VacuumSonifier().render_depth_tone(str(out), y=Y_C, duration=0.05)
    assert out.exists()
    with wave.open(str(out), "rb") as w:
        assert w.getnchannels() == 1
        assert w.getsampwidth() == 2
        assert w.getframerate() == 44100
        n = w.getnframes()
        assert n > 0
        frames = w.readframes(n)
    assert any(b != 0 for b in frames)  # nonzero audio
