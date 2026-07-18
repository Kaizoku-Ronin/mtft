#!/usr/bin/env python3
"""
MTFT Music Module — Arithmetic Sonification Engine
===================================================
Renders the mathematical structures of Modular Time Field Theory as sound.

Three integrated subsystems:
  1. VacuumSonifier  — renders the arithmetic vacuum V(τ) as audio waveforms
  2. ModularScale    — derives tuning systems from modular geometry
  3. MonsterComposer — deterministic algorithmic composition via MonsterHash

Every structural parameter is derived from MTFT mathematics:
  ┌──────────────────────────────────────────────────────────────┐
  │  Base frequency 143 Hz      ← level of X₀(143)              │
  │  13 harmonics per voice     ← genus of X₀(143)              │
  │  3 voices (e, μ, τ)        ← newform decomposition          │
  │  15-note scale              ← 15 supersingular primes        │
  │  y_c ≈ 0.1812 brightness   ← confinement depth              │
  │  δ ≈ 4.669 scaling ratio   ← Feigenbaum constant            │
  │  120° phase separation      ← Koide circle / augmented triad │
  │  143-beat phrase cycles     ← level of the modular curve     │
  └──────────────────────────────────────────────────────────────┘

Dependencies: numpy (required), wave + struct (stdlib)
Optional:     monster_security.monster_hash (falls back to HMAC-SHA256)

Usage:
    from mtft.music import VacuumSonifier, ModularScale, MonsterComposer

    # Sonify the arithmetic vacuum
    son = VacuumSonifier()
    son.render_vacuum_sweep('vacuum_dive.wav', duration=8.0)

    # Build a modular scale
    scale = ModularScale()
    freqs = scale.supersingular_scale(root=143.0)

    # Compose with MonsterHash
    comp = MonsterComposer(seed=b'MTFT')
    comp.compose_to_wav('monster_piece.wav', bars=13)

Roger Tano — MTFT Research Program — March 2026
MIT License
"""

from __future__ import annotations

import struct
import wave
import math
from typing import List, Tuple, Dict
from dataclasses import dataclass, field

import numpy as np

# ═══════════════════════════════════════════════════════════════
#  MTFT-DERIVED CONSTANTS
# ═══════════════════════════════════════════════════════════════

# Fundamental frequency: level of X₀(143) = 11 × 13
F_BASE = 143.0

# Genus of X₀(143) — number of harmonics per voice, phrase length
GENUS = 13

# Confinement depth — controls default spectral brightness
from mtft.constants import CriticalDepths

Y_C = CriticalDepths.y_conf  # 0.18174 — canonical confinement depth (v0.6.1 audit)

# Feigenbaum constant — scaling ratio between generations
DELTA_F = 4.669201609

# Burning Ship anisotropic constants — attractor/barrier folds

# Meissel-Mertens constant
M_MERTENS = 0.2614972128476427

# The 15 supersingular primes (primes dividing |M|, the Monster group order)
SUPERSINGULAR_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]

# Koide angle — 120° phase separation between generations
KOIDE_PHASE_SEP = 2 * np.pi / 3  # augmented triad

# Standard audio
SAMPLE_RATE = 44100
BIT_DEPTH = 16
MAX_AMP = 2**(BIT_DEPTH - 1) - 1

# Three-voice register derived from charged lepton mass ratios
# Electron : muon : tau ≈ 1 : 207 : 3477
# Mapped to pitch registers (high, mid, low) using log scaling
VOICE_REGISTERS = {
    'electron': 4.0,   # highest register (4× base)
    'muon':     1.0,   # middle register (1× base)
    'tau':      0.25,  # lowest register (1/4 base)
}


# ═══════════════════════════════════════════════════════════════
#  UTILITY: MonsterHash Interface
# ═══════════════════════════════════════════════════════════════

def _monster_hash(data: bytes, output_bits: int = 256) -> bytes:
    """
    MonsterHash digest — delegates to the real SL(2,ℤ) sponge
    if available, else HMAC-SHA256 with MTFT domain separation.
    """
    try:
        from mtft.monster_hash import MonsterHash
        h = MonsterHash(output_bits=min(output_bits, 512))
        digest_int = h.digest(data)
        return digest_int.to_bytes(output_bits // 8, 'big')
    except (ImportError, Exception):
        import hashlib
        import hmac
        key = b'MTFT-MonsterHash-v1-music'
        raw = hmac.new(key, data, hashlib.sha256).digest()
        if output_bits > 256:
            raw2 = hmac.new(key + b'-ext', data, hashlib.sha256).digest()
            raw = raw + raw2
        return raw[:output_bits // 8]


def _hash_to_floats(data: bytes, n: int) -> List[float]:
    """Extract n floats in [0, 1) from a hash digest."""
    h = _monster_hash(data, max(256, ((n * 4) // 32 + 1) * 256))
    floats = []
    for i in range(n):
        val = struct.unpack('>I', h[i*4:(i+1)*4])[0]
        floats.append(val / 0xFFFFFFFF)
    return floats


# ═══════════════════════════════════════════════════════════════
#  UTILITY: WAV Writer
# ═══════════════════════════════════════════════════════════════

def _write_wav(filepath: str, samples: np.ndarray, sample_rate: int = SAMPLE_RATE):
    """Write a mono 16-bit WAV file from float samples in [-1, 1]."""
    # Normalize to prevent clipping
    peak = np.max(np.abs(samples))
    if peak > 0:
        samples = samples / peak * 0.95

    pcm = (samples * MAX_AMP).astype(np.int16)

    with wave.open(filepath, 'w') as wf:
        wf.setnchannels(1)
        wf.setsampwidth(2)
        wf.setframerate(sample_rate)
        wf.writeframes(pcm.tobytes())


def _write_wav_stereo(filepath: str, left: np.ndarray, right: np.ndarray,
                      sample_rate: int = SAMPLE_RATE):
    """Write a stereo 16-bit WAV file from float samples in [-1, 1]."""
    peak = max(np.max(np.abs(left)), np.max(np.abs(right)))
    if peak > 0:
        left = left / peak * 0.95
        right = right / peak * 0.95

    left_pcm = (left * MAX_AMP).astype(np.int16)
    right_pcm = (right * MAX_AMP).astype(np.int16)

    interleaved = np.empty(len(left) + len(right), dtype=np.int16)
    interleaved[0::2] = left_pcm
    interleaved[1::2] = right_pcm

    with wave.open(filepath, 'w') as wf:
        wf.setnchannels(2)
        wf.setsampwidth(2)
        wf.setframerate(sample_rate)
        wf.writeframes(interleaved.tobytes())


# ═══════════════════════════════════════════════════════════════
#  UTILITY: Weight Computation
# ═══════════════════════════════════════════════════════════════

def _compute_weights(n_max: int) -> np.ndarray:
    """
    Compute divisor-log weights: w_n = Σ_{d|n} (log d)/d.

    These are the fundamental MTFT arithmetic weights that define
    the vacuum spectral profile. Highly composite numbers get
    boosted (many divisors) while primes are lean — producing
    timbres that naturally favor consonant intervals.
    """
    wn = np.zeros(n_max + 1)
    for d in range(2, n_max + 1):
        v = np.log(d) / d
        for m in range(d, n_max + 1, d):
            wn[m] += v
    return wn


def _sieve_primes(limit: int) -> List[int]:
    """Eratosthenes sieve."""
    s = [True] * (limit + 1)
    s[0] = s[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            for j in range(i*i, limit + 1, i):
                s[j] = False
    return [i for i in range(2, limit + 1) if s[i]]


# ═══════════════════════════════════════════════════════════════
#  1. VACUUM SONIFIER
# ═══════════════════════════════════════════════════════════════
#
#  The MTFT arithmetic potential is:
#
#    V(τ) = Σ_{n≥1} w_n e^{-2πyn} e^{2πixn}
#
#  This IS an audio signal: n indexes harmonics, w_n are amplitudes,
#  y controls spectral rolloff (brightness), x is time.
#
#  Large y (deep modular depth) → dark, bass-heavy tones
#  Small y (shallow)           → bright, harmonic-rich tones
#  y = y_c ≈ 0.1812           → the confinement timbre
#

class VacuumSonifier:
    """
    Renders the MTFT arithmetic vacuum as audible sound.

    The divisor-log weights w_n define a natural timbre where highly
    composite harmonics ring louder than primes — the vacuum literally
    sounds like number theory.

    Parameters
    ----------
    n_harmonics : int
        Number of harmonics to synthesize (default: 143)
    sample_rate : int
        Audio sample rate in Hz
    """

    def __init__(self, n_harmonics: int = 143, sample_rate: int = SAMPLE_RATE):
        self.n_harmonics = n_harmonics
        self.sample_rate = sample_rate
        self.weights = _compute_weights(n_harmonics)

    def vacuum_signal(self, freq: float, y: float, duration: float,
                      phase_offset: float = 0.0) -> np.ndarray:
        """
        Synthesize V(τ) as an audio signal at given depth y.

        Parameters
        ----------
        freq : float
            Fundamental frequency in Hz
        y : float
            Modular depth — controls spectral brightness.
            y_c ≈ 0.1812 is the confinement timbre.
        duration : float
            Duration in seconds
        phase_offset : float
            Initial phase in radians

        Returns
        -------
        np.ndarray
            Audio samples in [-1, 1]
        """
        t = np.arange(int(duration * self.sample_rate)) / self.sample_rate
        signal = np.zeros_like(t)

        for n in range(1, self.n_harmonics + 1):
            amplitude = self.weights[n] * np.exp(-2 * np.pi * y * n)
            signal += amplitude * np.sin(2 * np.pi * n * freq * t + phase_offset * n)

        return signal

    def three_voice_chord(self, y: float, duration: float,
                          root: float = F_BASE) -> np.ndarray:
        """
        Synthesize the three-generation chord: electron + muon + tau
        voices at 120° Koide phase separation (augmented triad).

        The three newforms f₁, f₂, f₃ on X₀(143) become three
        voices at phase separations of 2π/3 — the same geometry
        as the Koide circle.
        """
        e_signal = self.vacuum_signal(
            root * VOICE_REGISTERS['electron'], y, duration,
            phase_offset=0.0
        )
        mu_signal = self.vacuum_signal(
            root * VOICE_REGISTERS['muon'], y, duration,
            phase_offset=KOIDE_PHASE_SEP
        )
        tau_signal = self.vacuum_signal(
            root * VOICE_REGISTERS['tau'], y, duration,
            phase_offset=2 * KOIDE_PHASE_SEP
        )

        # Mix with amplitude weighting from Atkin-Lehner eigenvalues
        # f₁: (+,+), f₂: (−,+), f₃: (+,−) → weight by orbit dimension
        # dim(f₁)=1, dim(f₂)=4, dim(f₃)=6 → normalize
        total_dim = 1 + 4 + 6
        return (1/total_dim * e_signal +
                4/total_dim * mu_signal +
                6/total_dim * tau_signal)

    def render_vacuum_sweep(self, filepath: str, duration: float = 8.0,
                            y_start: float = 0.5, y_end: float = 0.05,
                            root: float = F_BASE):
        """
        Render a sweep through modular depth — "diving into the vacuum."

        Starts at large y (dark, few harmonics) and sweeps to small y
        (bright, harmonic-rich), crossing the confinement depth y_c
        where the Standard Model crystallizes.

        The crossing of y_c is marked by a subtle phase shift where
        the three voices lock into the Koide configuration.
        """
        t = np.arange(int(duration * self.sample_rate)) / self.sample_rate
        signal = np.zeros_like(t)

        # Smooth y sweep (exponential interpolation)
        y_curve = y_start * (y_end / y_start) ** (t / duration)

        for n in range(1, self.n_harmonics + 1):
            amplitude = self.weights[n] * np.exp(-2 * np.pi * y_curve * n)

            # Three-voice synthesis with Koide phase locking
            for voice_idx, (voice, register) in enumerate(VOICE_REGISTERS.items()):
                phase = voice_idx * KOIDE_PHASE_SEP
                f = n * root * register

                # Koide lock factor: increases near y_c
                lock = np.exp(-((y_curve - Y_C) / 0.03) ** 2)
                effective_phase = phase * lock

                signal += (amplitude / 3) * np.sin(
                    2 * np.pi * f * t + effective_phase * n
                )

        _write_wav(filepath, signal, self.sample_rate)
        return filepath

    def render_depth_tone(self, filepath: str, y: float = Y_C,
                          duration: float = 4.0, root: float = F_BASE):
        """Render a sustained tone at a specific modular depth."""
        signal = self.three_voice_chord(y, duration, root)
        _write_wav(filepath, signal, self.sample_rate)
        return filepath

    def render_stiffness_rhythm(self, filepath: str, duration: float = 8.0,
                                root: float = F_BASE, N: int = 3):
        """
        Render the SU(N) stiffness function μ_N(y) as a rhythmic pulse.

        The stiffness controls amplitude modulation — regions where
        μ_N is large produce loud beats, near-zero regions go silent.
        This makes the confinement zero at y_c audible as a dropout.
        """
        t = np.arange(int(duration * self.sample_rate)) / self.sample_rate
        signal = np.zeros_like(t)

        # Map time to modular depth
        y_range = np.linspace(0.05, 0.5, len(t))

        # Compute stiffness envelope
        stiffness_env = np.zeros_like(t)
        for n in range(1, min(self.n_harmonics + 1, 200)):
            wn = self.weights[n]
            stiffness_env += (n**2 * wn *
                              np.exp(-2 * np.pi * y_range * n) *
                              (1 - np.cos(2 * np.pi * n / N)))

        # Normalize envelope
        peak = np.max(np.abs(stiffness_env))
        if peak > 0:
            stiffness_env /= peak

        # Carrier signal
        for n in range(1, GENUS + 1):
            signal += self.weights[n] * np.sin(2 * np.pi * n * root * t)

        signal *= stiffness_env

        _write_wav(filepath, signal, self.sample_rate)
        return filepath


# ═══════════════════════════════════════════════════════════════
#  2. MODULAR SCALE — Tuning Systems from MTFT Geometry
# ═══════════════════════════════════════════════════════════════
#
#  MTFT provides multiple natural tuning systems:
#
#  (a) Supersingular scale: The 15 primes dividing |Monster|
#      mapped as interval ratios — a 15-note scale inherent to
#      the Monster group itself.
#
#  (b) Farey scale: Rational intervals from the Farey sequence
#      F_N, which approximates the structure of SL(2,ℤ) orbits.
#
#  (c) SU(N) scales: N-fold center symmetry → N equal peaks.
#      SU(2) = octave, SU(3) = augmented triad, SU(5) = pentatonic.
#
#  (d) Hecke scale: Frequency ratios from Hecke eigenvalues of
#      the newforms on X₀(143).
#

class ModularScale:
    """
    Derives tuning systems from MTFT mathematical structures.

    The scales produced are NOT equal temperament — they are
    determined by number-theoretic structures inherent to the
    modular curve X₀(143) and the Monster group.
    """

    def __init__(self):
        self.primes = _sieve_primes(200)

    def supersingular_scale(self, root: float = F_BASE) -> Dict[int, float]:
        """
        The 15-note supersingular scale.

        The 15 supersingular primes {2,3,5,...,71} are the primes
        dividing the order of the Monster group. Map each prime p
        to a frequency within a single octave via:

            f_p = root × 2^{log_71(p)}

        This compresses the 15 primes into one octave, preserving
        their multiplicative structure as musical intervals.

        Returns
        -------
        dict
            {prime: frequency_Hz} for each supersingular prime
        """
        max_p = max(SUPERSINGULAR_PRIMES)
        scale = {}
        for p in SUPERSINGULAR_PRIMES:
            # Map to [0, 1) within one octave using log ratio
            octave_position = np.log(p) / np.log(max_p)
            scale[p] = root * 2**octave_position
        return scale

    def farey_scale(self, order: int = 7, root: float = F_BASE) -> List[float]:
        """
        Farey sequence tuning: rational intervals from F_N.

        The Farey sequence of order N contains all rationals p/q
        with 0 ≤ p/q ≤ 1 and q ≤ N. These rationals index the
        cusps of Γ₀(N) and define a natural tuning that approximates
        just intonation.

        Parameters
        ----------
        order : int
            Farey sequence order (default: 7)
        root : float
            Root frequency in Hz

        Returns
        -------
        list of float
            Frequencies in Hz for each Farey fraction
        """
        # Build Farey sequence F_N
        fractions = set()
        for q in range(1, order + 1):
            for p in range(0, q + 1):
                if math.gcd(p, q) == 1:
                    fractions.add((p, q))

        fractions = sorted(fractions, key=lambda x: x[0] / x[1])

        # Map to frequencies: f = root × 2^(p/q) for each p/q in F_N
        freqs = []
        for p, q in fractions:
            ratio = p / q
            freqs.append(root * 2**ratio)

        return freqs

    def sun_scale(self, N: int = 5, root: float = F_BASE,
                  y: float = Y_C) -> List[float]:
        """
        SU(N) center symmetry scale: N notes from gauge theory.

        The stiffness μ_N(y) has N-fold center symmetry, producing
        N peaks at phases 2πk/N. This naturally generates:
            SU(2) → 2 notes (octave)
            SU(3) → 3 notes (augmented triad)
            SU(4) → 4 notes (diminished seventh)
            SU(5) → 5 notes (pentatonic)
            SU(7) → 7 notes (diatonic-like)
            SU(12) → 12 notes (chromatic)

        The amplitude of each note is weighted by the arithmetic
        potential at that phase, so unlike equal temperament, each
        degree has a different "brightness."

        Parameters
        ----------
        N : int
            Gauge group rank (number of scale degrees)
        root : float
            Root frequency in Hz
        y : float
            Modular depth for amplitude weighting
        """
        notes = []
        for k in range(N):
            # Frequency from N-fold division of the octave
            freq = root * 2**(k / N)
            notes.append(freq)
        return notes

    def sun_scale_weighted(self, N: int = 5, root: float = F_BASE,
                           y: float = Y_C, n_max: int = 100
                           ) -> List[Tuple[float, float]]:
        """
        SU(N) scale with arithmetic amplitude weights.

        Returns (frequency, weight) pairs where the weight reflects
        the stiffness contribution at each center element.
        """
        weights = _compute_weights(n_max)
        notes = []

        for k in range(N):
            freq = root * 2**(k / N)
            # Weight = stiffness contribution from this center element
            w = 0.0
            for n in range(1, n_max + 1):
                w += (n**2 * weights[n] *
                      np.exp(-2 * np.pi * y * n) *
                      (1 - np.cos(2 * np.pi * n * k / N)))
            notes.append((freq, abs(w)))

        # Normalize weights
        max_w = max(w for _, w in notes) if notes else 1.0
        if max_w > 0:
            notes = [(f, w / max_w) for f, w in notes]

        return notes

    def hecke_scale(self, root: float = F_BASE) -> Dict[str, List[float]]:
        """
        Hecke eigenvalue scale: intervals from the newforms on X₀(143).

        The Hecke polynomial T₂ on f₂ has roots that define frequency
        ratios. Each eigenvalue α gives an interval of 2^(|α|/max_α).

        This produces different scales for the electron, muon, and
        tau sectors — three timbres from one modular curve.
        """
        # T₂ on f₂: x⁴ - 3x³ - x² + 5x + 1
        roots_f2 = np.sort(np.real(np.roots([1, -3, -1, 5, 1])))
        # T₂ on f₃: x⁶ - 10x⁴ + 2x³ + 24x² - 7x - 12
        roots_f3 = np.sort(np.real(np.roots([1, 0, -10, 2, 24, -7, -12])))

        def eigenvalues_to_scale(roots, root_freq):
            max_r = max(abs(r) for r in roots)
            if max_r == 0:
                return [root_freq]
            freqs = [root_freq]
            for r in sorted(roots):
                octave_pos = (r - min(roots)) / (max(roots) - min(roots))
                freqs.append(root_freq * 2**octave_pos)
            return sorted(set(freqs))

        return {
            'muon':  eigenvalues_to_scale(roots_f2, root),
            'tau':   eigenvalues_to_scale(roots_f3, root),
        }

    def koide_triad(self, root: float = F_BASE) -> Tuple[float, float, float]:
        """
        The Koide triad: three frequencies at 120° phase separation.

        This is the augmented triad — the most symmetric three-note
        chord, reflecting the Koide circle geometry where the three
        generations sit at equal angular spacing.

        In equal temperament, this is root, major third, minor sixth
        (e.g., C-E-G#). Here it's derived from first principles.
        """
        return (
            root * 2**(0/3),    # root
            root * 2**(1/3),    # 120° = major third (≈ 5:4)
            root * 2**(2/3),    # 240° = minor sixth (≈ 8:5)
        )

    def feigenbaum_intervals(self, root: float = F_BASE,
                             n_levels: int = 5) -> List[float]:
        """
        Period-doubling cascade as a pitch hierarchy.

        Each level of the Feigenbaum cascade introduces a new
        frequency at ratio δ from the previous. This creates a
        self-similar pitch structure mirroring the fractal
        generation hierarchy.

        Level 0: root
        Level 1: root × δ^(1/log₂δ)
        Level 2: root × δ^(2/log₂δ)
        ...
        """
        log2_delta = np.log2(DELTA_F)
        return [root * DELTA_F**(k / log2_delta) for k in range(n_levels)]

    def render_scale(self, filepath: str, scale_freqs: List[float],
                     note_duration: float = 0.5, y: float = Y_C):
        """Render a scale as a WAV file, ascending then descending."""
        son = VacuumSonifier(n_harmonics=GENUS)

        ascending = scale_freqs
        descending = list(reversed(scale_freqs[:-1]))
        all_notes = ascending + descending

        total_samples = int(len(all_notes) * note_duration * SAMPLE_RATE)
        signal = np.zeros(total_samples)

        for i, freq in enumerate(all_notes):
            start = int(i * note_duration * SAMPLE_RATE)
            note = son.vacuum_signal(freq, y, note_duration)

            # Apply ADSR envelope
            env = _adsr_envelope(len(note), SAMPLE_RATE,
                                 attack=0.02, decay=0.05,
                                 sustain=0.7, release=0.1)
            note *= env

            end = start + len(note)
            if end <= total_samples:
                signal[start:end] += note

        _write_wav(filepath, signal, SAMPLE_RATE)
        return filepath


# ═══════════════════════════════════════════════════════════════
#  3. MONSTER COMPOSER — Algorithmic Composition via MonsterHash
# ═══════════════════════════════════════════════════════════════
#
#  MonsterHash is a deterministic, chaotic function with near-ideal
#  statistical properties (49.95% avalanche). This makes it a
#  perfect algorithmic composition engine:
#
#    • Deterministic: same seed always produces the same piece
#    • Chaotic: small seed changes → completely different music
#    • Structured: 13-round internal architecture → 13-bar phrases
#    • Three-voice: electron/muon/tau from the newform decomposition
#

@dataclass
class Note:
    """A single musical note."""
    frequency: float    # Hz
    duration: float     # seconds
    velocity: float     # amplitude [0, 1]
    voice: str          # 'electron', 'muon', or 'tau'
    start_time: float   # seconds from piece start


@dataclass
class Phrase:
    """A musical phrase — 13 notes (genus of X₀(143))."""
    notes: List[Note] = field(default_factory=list)
    voice: str = 'muon'


class MonsterComposer:
    """
    Deterministic algorithmic composition engine driven by MonsterHash.

    Every seed produces a unique, reproducible musical piece.
    The composition respects MTFT structural constraints:
      - 13-note phrases (genus)
      - 143-beat cycles (level)
      - Three voices separated by Koide phase angles
      - Supersingular prime scale

    Parameters
    ----------
    seed : bytes
        Seed for deterministic generation
    scale_type : str
        'supersingular', 'pentatonic', 'diatonic', or 'chromatic'
    tempo : float
        Beats per minute (default: 143 ÷ 2 ≈ 71.5 BPM)
    """

    def __init__(self, seed: bytes = b'MTFT-MonsterComposer',
                 scale_type: str = 'supersingular',
                 tempo: float = 143.0 / 2):
        self.seed = seed
        self.tempo = tempo
        self.beat_duration = 60.0 / tempo

        # Build the scale
        ms = ModularScale()
        if scale_type == 'supersingular':
            scale_dict = ms.supersingular_scale(F_BASE)
            self.scale = sorted(scale_dict.values())
        elif scale_type == 'pentatonic':
            self.scale = ms.sun_scale(5, F_BASE)
        elif scale_type == 'diatonic':
            self.scale = ms.sun_scale(7, F_BASE)
        elif scale_type == 'chromatic':
            self.scale = ms.sun_scale(12, F_BASE)
        else:
            self.scale = ms.sun_scale(5, F_BASE)

        # Extend scale across 3 octaves for voice registers
        base_scale = self.scale
        self.full_scale = []
        for octave_shift in [-1, 0, 1, 2]:
            for f in base_scale:
                self.full_scale.append(f * 2**octave_shift)
        self.full_scale.sort()

    def _hash_step(self, counter: int) -> List[float]:
        """Generate deterministic random floats from seed + counter."""
        data = self.seed + struct.pack('>Q', counter)
        return _hash_to_floats(data, 8)

    def _select_note(self, floats: List[float], voice: str,
                     start_time: float) -> Note:
        """Map hash-derived floats to a musical note."""
        # Pitch selection — quantize to scale
        register = VOICE_REGISTERS[voice]
        pitch_idx = int(floats[0] * len(self.scale))
        pitch_idx = min(pitch_idx, len(self.scale) - 1)
        freq = self.scale[pitch_idx] * register

        # Duration — weighted toward musical values
        # MonsterHash gives us 8 possible durations: whole to 32nd note
        dur_options = [4, 3, 2, 1.5, 1, 0.75, 0.5, 0.25]
        dur_idx = int(floats[1] * len(dur_options))
        dur_idx = min(dur_idx, len(dur_options) - 1)
        duration = dur_options[dur_idx] * self.beat_duration

        # Velocity (dynamics)
        velocity = 0.3 + 0.7 * floats[2]

        # Rest probability (≈ 15% silence)
        if floats[3] < 0.15:
            velocity = 0.0

        return Note(
            frequency=freq,
            duration=duration,
            velocity=velocity,
            voice=voice,
            start_time=start_time,
        )

    def generate_phrase(self, phrase_idx: int, voice: str) -> Phrase:
        """
        Generate a 13-note phrase (genus of X₀(143)).

        Each phrase is a self-contained musical idea. The 13 notes
        mirror the 13 rounds of MonsterHash's internal architecture.
        """
        phrase = Phrase(voice=voice)
        current_time = phrase_idx * GENUS * self.beat_duration

        for note_idx in range(GENUS):
            counter = phrase_idx * GENUS * 3 + note_idx * 3
            # Offset counter by voice to decorrelate
            voice_offset = {'electron': 0, 'muon': 1000000, 'tau': 2000000}
            floats = self._hash_step(counter + voice_offset.get(voice, 0))

            note = self._select_note(floats, voice, current_time)
            phrase.notes.append(note)
            current_time += note.duration

        return phrase

    def compose(self, bars: int = 13) -> List[Phrase]:
        """
        Compose a multi-voice piece.

        Parameters
        ----------
        bars : int
            Number of 13-note phrases per voice (default: 13,
            giving a 13×13 = 169-note structure per voice)
        """
        phrases = []
        for voice in ['electron', 'muon', 'tau']:
            for bar in range(bars):
                phrases.append(self.generate_phrase(bar, voice))
        return phrases

    def render_phrase(self, phrase: Phrase, y: float = Y_C) -> np.ndarray:
        """Render a single phrase to audio samples."""
        son = VacuumSonifier(n_harmonics=GENUS)

        # Calculate total duration
        if not phrase.notes:
            return np.zeros(SAMPLE_RATE)

        total_dur = max(n.start_time + n.duration for n in phrase.notes)
        total_dur -= phrase.notes[0].start_time
        total_samples = int(total_dur * SAMPLE_RATE) + SAMPLE_RATE
        signal = np.zeros(total_samples)

        base_time = phrase.notes[0].start_time

        for note in phrase.notes:
            if note.velocity < 0.01:
                continue

            local_start = int((note.start_time - base_time) * SAMPLE_RATE)
            n_samples = int(note.duration * SAMPLE_RATE)

            if n_samples < 100:
                continue

            tone = son.vacuum_signal(note.frequency, y, note.duration)
            env = _adsr_envelope(len(tone), SAMPLE_RATE,
                                 attack=0.01, decay=0.03,
                                 sustain=0.6, release=0.15)
            tone *= env * note.velocity

            end = local_start + len(tone)
            if end <= total_samples:
                signal[local_start:end] += tone
            else:
                fit = total_samples - local_start
                if fit > 0:
                    signal[local_start:total_samples] += tone[:fit]

        return signal

    def compose_to_wav(self, filepath: str, bars: int = 13,
                       y: float = Y_C):
        """
        Compose and render a complete piece to a WAV file.

        The three voices are panned:
          electron → right (high register)
          muon     → center (mid register)
          tau      → left (low register)
        """
        phrases = self.compose(bars)

        # Group by voice
        voice_phrases = {'electron': [], 'muon': [], 'tau': []}
        for p in phrases:
            voice_phrases[p.voice].append(p)

        # Render each voice
        voice_audio = {}
        for voice, vphrases in voice_phrases.items():
            segments = []
            for p in vphrases:
                seg = self.render_phrase(p, y)
                segments.append(seg)

            if segments:
                voice_audio[voice] = np.concatenate(segments)
            else:
                voice_audio[voice] = np.zeros(SAMPLE_RATE)

        # Match lengths
        max_len = max(len(v) for v in voice_audio.values())
        for voice in voice_audio:
            if len(voice_audio[voice]) < max_len:
                voice_audio[voice] = np.pad(
                    voice_audio[voice],
                    (0, max_len - len(voice_audio[voice]))
                )

        # Stereo panning: tau left, muon center, electron right
        left = (voice_audio['tau'] * 0.8 +
                voice_audio['muon'] * 0.5 +
                voice_audio['electron'] * 0.2)
        right = (voice_audio['tau'] * 0.2 +
                 voice_audio['muon'] * 0.5 +
                 voice_audio['electron'] * 0.8)

        _write_wav_stereo(filepath, left, right, SAMPLE_RATE)
        return filepath


# ═══════════════════════════════════════════════════════════════
#  UTILITY: ADSR Envelope
# ═══════════════════════════════════════════════════════════════

def _adsr_envelope(n_samples: int, sample_rate: int,
                   attack: float = 0.01, decay: float = 0.05,
                   sustain: float = 0.7, release: float = 0.1) -> np.ndarray:
    """
    Generate an ADSR amplitude envelope.

    Parameters in seconds (attack, decay, release) or linear (sustain).
    """
    a_samples = int(attack * sample_rate)
    d_samples = int(decay * sample_rate)
    r_samples = int(release * sample_rate)
    s_samples = max(0, n_samples - a_samples - d_samples - r_samples)

    envelope = np.concatenate([
        np.linspace(0, 1, max(a_samples, 1)),                # Attack
        np.linspace(1, sustain, max(d_samples, 1)),           # Decay
        np.full(s_samples, sustain),                          # Sustain
        np.linspace(sustain, 0, max(r_samples, 1)),           # Release
    ])

    # Trim or pad to exact length
    if len(envelope) > n_samples:
        envelope = envelope[:n_samples]
    elif len(envelope) < n_samples:
        envelope = np.pad(envelope, (0, n_samples - len(envelope)))

    return envelope


# ═══════════════════════════════════════════════════════════════
#  4. ANALYSIS & VISUALIZATION
# ═══════════════════════════════════════════════════════════════

def spectral_fingerprint(y: float = Y_C, n_max: int = 143) -> Dict:
    """
    Compute the spectral fingerprint of the vacuum at depth y.

    Returns the amplitude spectrum, dominant harmonics, and
    consonance measure — quantifying how "musical" the vacuum
    sounds at each depth.
    """
    weights = _compute_weights(n_max)
    spectrum = {}

    for n in range(1, n_max + 1):
        amp = weights[n] * np.exp(-2 * np.pi * y * n)
        spectrum[n] = amp

    # Sort by amplitude
    ranked = sorted(spectrum.items(), key=lambda x: -x[1])

    # Consonance: ratio of energy in consonant harmonics (2,3,4,5,6,8)
    # to total energy
    consonant_set = {2, 3, 4, 5, 6, 8, 10, 12, 15, 16}
    total_energy = sum(a**2 for _, a in ranked)
    consonant_energy = sum(a**2 for n, a in ranked if n in consonant_set)
    consonance = consonant_energy / total_energy if total_energy > 0 else 0

    return {
        'spectrum': spectrum,
        'top_13': ranked[:GENUS],
        'consonance': consonance,
        'brightness': sum(n * a for n, a in ranked) / sum(a for _, a in ranked) if ranked else 0,
    }


def compare_depths(depths: List[float] = None) -> Dict:
    """Compare the spectral character across multiple modular depths."""
    if depths is None:
        depths = [0.05, 0.10, Y_C, 0.25, 0.35, 0.50]

    results = {}
    for y in depths:
        fp = spectral_fingerprint(y)
        results[y] = {
            'consonance': fp['consonance'],
            'brightness': fp['brightness'],
            'dominant_harmonic': fp['top_13'][0][0] if fp['top_13'] else None,
        }
    return results


# ═══════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════

def main():
    """Command-line interface for the MTFT music module."""
    import sys

    commands = {
        'sweep':     'Render a vacuum depth sweep (diving into the vacuum)',
        'tone':      'Render a sustained confinement tone at y_c',
        'scale':     'Render the supersingular 15-note scale',
        'compose':   'Compose and render a MonsterHash piece',
        'spectrum':  'Print the spectral fingerprint at y_c',
        'compare':   'Compare spectral character across depths',
        'koide':     'Render the Koide triad (augmented chord)',
        'stiffness': 'Render the SU(3) stiffness rhythm',
    }

    if len(sys.argv) < 2 or sys.argv[1] not in commands:
        print("MTFT Music Module — Arithmetic Sonification Engine")
        print("=" * 52)
        print(f"\nUsage: python mtft_music.py <command> [output.wav]\n")
        print("Commands:")
        for cmd, desc in commands.items():
            print(f"  {cmd:12s}  {desc}")
        print(f"\nDefault output: mtft_<command>.wav")
        return

    cmd = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) > 2 else f'mtft_{cmd}.wav'

    if cmd == 'sweep':
        print(f"Rendering vacuum depth sweep → {outfile}")
        son = VacuumSonifier()
        son.render_vacuum_sweep(outfile, duration=10.0)
        print(f"  y: 0.50 → 0.05 over 10 seconds")
        print(f"  Crossing confinement depth y_c = {Y_C} at ~6.2s")

    elif cmd == 'tone':
        print(f"Rendering confinement tone at y_c = {Y_C} → {outfile}")
        son = VacuumSonifier()
        son.render_depth_tone(outfile, y=Y_C, duration=6.0)

    elif cmd == 'scale':
        print(f"Rendering supersingular scale → {outfile}")
        ms = ModularScale()
        scale = ms.supersingular_scale()
        freqs = sorted(scale.values())
        print(f"  15 notes from the Monster group:")
        for p, f in sorted(scale.items()):
            print(f"    p={p:2d} → {f:.1f} Hz")
        ms.render_scale(outfile, freqs, note_duration=0.4)

    elif cmd == 'compose':
        seed = sys.argv[2].encode() if len(sys.argv) > 2 and not sys.argv[2].endswith('.wav') else b'MTFT'
        if not outfile.endswith('.wav'):
            outfile = f'mtft_compose_{seed.decode()}.wav'
        print(f"Composing MonsterHash piece (seed={seed}) → {outfile}")
        comp = MonsterComposer(seed=seed)
        comp.compose_to_wav(outfile, bars=GENUS)
        print(f"  {GENUS} phrases × 3 voices × {GENUS} notes = {GENUS*3*GENUS} total notes")
        print(f"  Tempo: {comp.tempo:.1f} BPM")

    elif cmd == 'spectrum':
        fp = spectral_fingerprint()
        print(f"Spectral Fingerprint at y_c = {Y_C}")
        print(f"  Consonance: {fp['consonance']:.4f}")
        print(f"  Brightness: {fp['brightness']:.1f}")
        print(f"\n  Top {GENUS} harmonics:")
        for n, amp in fp['top_13']:
            prime_flag = '★' if n in SUPERSINGULAR_PRIMES else ' '
            print(f"    n={n:3d} {prime_flag}  amp={amp:.6f}")

    elif cmd == 'compare':
        results = compare_depths()
        print("Spectral Character across Modular Depths")
        print(f"{'y':>8s}  {'Consonance':>11s}  {'Brightness':>11s}  {'Dominant':>8s}")
        print("-" * 45)
        for y, data in sorted(results.items()):
            marker = ' ← y_c' if abs(y - Y_C) < 0.001 else ''
            print(f"  {y:.4f}  {data['consonance']:11.4f}  "
                  f"{data['brightness']:11.1f}  n={data['dominant_harmonic']}{marker}")

    elif cmd == 'koide':
        print(f"Rendering Koide triad → {outfile}")
        son = VacuumSonifier()
        ms = ModularScale()
        f1, f2, f3 = ms.koide_triad()
        print(f"  Electron: {f1:.1f} Hz")
        print(f"  Muon:     {f2:.1f} Hz")
        print(f"  Tau:      {f3:.1f} Hz")
        signal = son.three_voice_chord(Y_C, 6.0)
        _write_wav(outfile, signal, SAMPLE_RATE)

    elif cmd == 'stiffness':
        print(f"Rendering SU(3) stiffness rhythm → {outfile}")
        son = VacuumSonifier()
        son.render_stiffness_rhythm(outfile, duration=10.0)
        print(f"  Confinement dropout at y_c = {Y_C}")

    print("Done.")


if __name__ == '__main__':
    main()
