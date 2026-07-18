"""
Hearing the Vacuum — MTFT Sonification
=======================================

The v0.7.0 music module maps MTFT structures to sound: modular scales
built from the supersingular primes and the Koide circle, vacuum tones
at the canonical confinement depth, and a deterministic composer driven
by MonsterHash digests.

  pip install mtft
  python examples/15_sonification.py

Writes three small WAV files to the current directory.
"""

from mtft.constants import CriticalDepths
from mtft.music import (
    ModularScale,
    MonsterComposer,
    VacuumSonifier,
    spectral_fingerprint,
)

# ── 1. Modular scales ────────────────────────────────────────
print("MODULAR SCALES  (root = 143 Hz, the level)")
print("=" * 55)
ms = ModularScale()

ss = ms.supersingular_scale()
firsts = {p: ss[p] for p in sorted(ss)[:5]}
print("  Supersingular scale (15 Monster primes), first five:")
for p, f in firsts.items():
    print(f"    p={p:>2d} -> {f:8.2f} Hz")

t = ms.koide_triad()
print(f"  Koide triad (120 deg spacing): "
      f"{t[0]:.1f}, {t[1]:.1f}, {t[2]:.1f} Hz — the augmented chord")

fi = ms.feigenbaum_intervals(n_levels=4)
print(f"  Feigenbaum cascade intervals: {[f'{f:.1f}' for f in fi]}")
print()

# ── 2. Spectral fingerprint of the vacuum ────────────────────
print("VACUUM SPECTRAL FINGERPRINT at y_conf = 0.18174")
print("=" * 55)
fp = spectral_fingerprint(n_max=143)
print(f"  brightness (harmonic centroid): {fp['brightness']:.2f}")
print(f"  consonance score:               {fp['consonance']:.4f}")
top = ", ".join(f"n={n}" for n, _amp in fp["top_13"][:5])
print(f"  top-13 partials (by weight):    {top}, ...")
print()

# ── 3. Render audio ──────────────────────────────────────────
print("RENDERING WAV FILES")
print("=" * 55)
son = VacuumSonifier()

son.render_depth_tone("mtft_vacuum_tone.wav",
                      y=CriticalDepths.y_conf, duration=1.5)
print("  mtft_vacuum_tone.wav      — the vacuum at confinement depth")

son.render_stiffness_rhythm("mtft_stiffness_rhythm.wav", duration=2.0)
print("  mtft_stiffness_rhythm.wav — mu_3(y) as amplitude modulation")

MonsterComposer().compose_to_wav("mtft_monster_composition.wav", bars=2)
print("  mtft_monster_composition.wav — 2 bars of MonsterHash counterpoint")
print()
print("Deterministic by construction: same digests, same music.")
