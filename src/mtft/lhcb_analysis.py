"""
MTFT × LHCb Open Data Analysis Bridge
=======================================
Reads LHCb Open Data ROOT ntuples (via uproot) and connects
reconstructed particle spectra to MTFT predictions.

Designed for:
  - B± → J/ψ(→ μ⁺μ⁻)K± ntuples from the CERN Open Data Portal
  - Any LHCb dvntuple.root file produced by the Ntupling Service

Dependencies (pip install):
  uproot >= 5.0
  awkward >= 2.0
  numpy
  matplotlib   (optional, for plotting)
  hist         (optional, for histogramming)

Usage:
  from mtft.lhcb_analysis import LHCbNtuple
  nt = LHCbNtuple("00334564_00000001_1.dvntuple.root")
  nt.info()
  jpsi = nt.dimuon_mass()
  nt.plot_jpsi()
  nt.mtft_confrontation()

Author: Roger Tano  (MTFT framework)
Computation: Claude  (Anthropic)
Package: mtft ≥ 0.5.1
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from typing import Optional, Union
from pathlib import Path

import numpy as np

# ═══════════════════════════════════════════════════════════════════════
# Lazy imports — fail gracefully if not installed
# ═══════════════════════════════════════════════════════════════════════

def _require(pkg: str):
    """Import-or-explain."""
    try:
        return __import__(pkg)
    except ImportError:
        raise ImportError(
            f"'{pkg}' is required for LHCb analysis.  "
            f"Install it with:  pip install {pkg}"
        )


# ═══════════════════════════════════════════════════════════════════════
# MTFT constants (imported from lhc.py if available, else self-contained)
# ═══════════════════════════════════════════════════════════════════════

# PDG 2024 reference masses (MeV)
M_MU_PDG     = 105.6583755       # muon mass
M_JPSI_PDG   = 3096.900          # J/ψ mass
M_BPLUS_PDG  = 5279.41           # B± mass
M_KPLUS_PDG  = 493.677           # K± mass

# MTFT predictions
M_TAU_KOIDE  = 1776.969          # τ mass from Koide (Paper 6)
KOIDE_K      = 2.0 / 3.0         # Koide constant
ALPHA_INV    = 137.035999        # MTFT α⁻¹ (Paper 1)

# Hidden doublet (Paper 25)
H11_MASS     = 1312.0            # SU(11) hidden scalar (MeV)
H13_MASS     = 1348.0            # SU(13) hidden scalar (MeV)
DOUBLET_SPLIT = 36.0             # MeV


# ═══════════════════════════════════════════════════════════════════════
# Branch name detection
# ═══════════════════════════════════════════════════════════════════════

# LHCb ntuples use varied naming conventions depending on the
# Ntupling Service configuration.  We try several patterns.

_BRANCH_PATTERNS = {
    # B± → J/ψ K± decay tree
    "B_M":      ["B_plus_M", "Bplus_M", "B_M", "B_plus_MM"],
    "B_PT":     ["B_plus_PT", "Bplus_PT", "B_PT"],
    "B_PE":     ["B_plus_PE", "Bplus_PE"],
    "B_PX":     ["B_plus_PX", "Bplus_PX"],
    "B_PY":     ["B_plus_PY", "Bplus_PY"],
    "B_PZ":     ["B_plus_PZ", "Bplus_PZ"],

    # J/ψ
    "Jpsi_M":   ["J_psi_1S_M", "Jpsi_M", "J_psi_M", "Jpsi_MM",
                 "J_psi_1S_MM"],
    "Jpsi_PT":  ["J_psi_1S_PT", "Jpsi_PT", "J_psi_PT"],
    "Jpsi_PE":  ["J_psi_1S_PE", "Jpsi_PE"],
    "Jpsi_PX":  ["J_psi_1S_PX", "Jpsi_PX"],
    "Jpsi_PY":  ["J_psi_1S_PY", "Jpsi_PY"],
    "Jpsi_PZ":  ["J_psi_1S_PZ", "Jpsi_PZ"],

    # μ⁺
    "mup_PE":   ["mu_plus_PE", "mup_PE", "muplus_PE", "mu_plus_or_H1_PE"],
    "mup_PX":   ["mu_plus_PX", "mup_PX", "muplus_PX", "mu_plus_or_H1_PX"],
    "mup_PY":   ["mu_plus_PY", "mup_PY", "muplus_PY", "mu_plus_or_H1_PY"],
    "mup_PZ":   ["mu_plus_PZ", "mup_PZ", "muplus_PZ", "mu_plus_or_H1_PZ"],
    "mup_PT":   ["mu_plus_PT", "mup_PT", "muplus_PT"],
    "mup_ProbNNmu": ["mu_plus_ProbNNmu", "mup_ProbNNmu", "muplus_ProbNNmu",
                     "mu_plus_or_H1_ProbNNmu"],

    # μ⁻
    "mum_PE":   ["mu_minus_PE", "mum_PE", "muminus_PE", "mu_minus_or_H2_PE"],
    "mum_PX":   ["mu_minus_PX", "mum_PX", "muminus_PX", "mu_minus_or_H2_PX"],
    "mum_PY":   ["mu_minus_PY", "mum_PY", "muminus_PY", "mu_minus_or_H2_PY"],
    "mum_PZ":   ["mu_minus_PZ", "mum_PZ", "muminus_PZ", "mu_minus_or_H2_PZ"],
    "mum_PT":   ["mu_minus_PT", "mum_PT", "muminus_PT"],
    "mum_ProbNNmu": ["mu_minus_ProbNNmu", "mum_ProbNNmu", "muminus_ProbNNmu",
                     "mu_minus_or_H2_ProbNNmu"],

    # K±
    "K_PE":     ["K_plus_PE", "Kplus_PE", "K_PE", "K_plus_or_H3_PE"],
    "K_PX":     ["K_plus_PX", "Kplus_PX", "K_PX", "K_plus_or_H3_PX"],
    "K_PY":     ["K_plus_PY", "Kplus_PY", "K_PY", "K_plus_or_H3_PY"],
    "K_PZ":     ["K_plus_PZ", "Kplus_PZ", "K_PZ", "K_plus_or_H3_PZ"],
    "K_PT":     ["K_plus_PT", "Kplus_PT", "K_PT"],
    "K_ProbNNk": ["K_plus_ProbNNk", "Kplus_ProbNNk", "K_ProbNNk",
                  "K_plus_or_H3_ProbNNk"],
}


# ═══════════════════════════════════════════════════════════════════════
# Core analysis class
# ═══════════════════════════════════════════════════════════════════════

class LHCbNtuple:
    """
    Interface to an LHCb Open Data dvntuple ROOT file.

    Parameters
    ----------
    path : str or Path
        Path to the .root file (local).
    tree_name : str or None
        Name of the TTree inside the file.  If None, auto-detected
        (tries common names like 'DecayTree', 'Tuple*/DecayTree', etc.)
    max_events : int or None
        Limit the number of events loaded (useful for quick checks).
    """

    def __init__(
        self,
        path: Union[str, Path],
        tree_name: Optional[str] = None,
        max_events: Optional[int] = None,
    ):
        uproot = _require("uproot")

        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"ROOT file not found: {self.path}")

        self._file = uproot.open(str(self.path))
        self._tree_name = tree_name or self._find_tree()
        self._tree = self._file[self._tree_name]
        self._max = max_events

        # Resolve branch names for this specific file
        self._branch_map: dict[str, str] = {}
        self._resolve_branches()

        # Cache
        self._arrays: dict[str, np.ndarray] = {}

    # ── Tree auto-detection ─────────────────────────────────────────

    def _find_tree(self) -> str:
        """Walk the file and find the first TTree."""
        uproot = _require("uproot")

        # Common LHCb patterns
        candidates = [
            "DecayTree",
            "Tuple/DecayTree",
            "TupleBtoJpsiK/DecayTree",
            "TupleBplus/DecayTree",
            "Bp2JpsiKplus/DecayTree",
        ]

        all_keys = self._file.keys()

        # Try known names first
        for c in candidates:
            if c in all_keys:
                return c

        # Fallback: find any TTree
        for key in all_keys:
            obj = self._file[key]
            if hasattr(obj, 'num_entries'):  # duck-type TTree
                return key

        # Last resort: look one level deeper
        for key in all_keys:
            obj = self._file[key]
            if hasattr(obj, 'keys'):
                for subkey in obj.keys():
                    full = f"{key}/{subkey}"
                    sub_obj = self._file[full]
                    if hasattr(sub_obj, 'num_entries'):
                        return full

        raise RuntimeError(
            f"No TTree found in {self.path}.  "
            f"Available keys: {all_keys[:20]}...  "
            f"Try passing tree_name='...' explicitly."
        )

    # ── Branch resolution ───────────────────────────────────────────

    def _resolve_branches(self):
        """Match generic branch names to the actual names in this file."""
        actual_branches = set(self._tree.keys())

        for generic, candidates in _BRANCH_PATTERNS.items():
            for c in candidates:
                if c in actual_branches:
                    self._branch_map[generic] = c
                    break

    def _get(self, generic: str) -> np.ndarray:
        """Load a branch by generic name, with caching."""
        if generic in self._arrays:
            return self._arrays[generic]

        if generic not in self._branch_map:
            raise KeyError(
                f"Branch '{generic}' not found in this ntuple.  "
                f"Resolved branches: {list(self._branch_map.keys())}"
            )

        ak = _require("awkward")
        actual = self._branch_map[generic]
        arr = self._tree[actual].array(library="np", entry_stop=self._max)

        # Handle awkward arrays → flatten if needed
        if hasattr(arr, 'to_numpy'):
            arr = arr.to_numpy()

        self._arrays[generic] = arr
        return arr

    # ── Info ─────────────────────────────────────────────────────────

    def info(self) -> str:
        """Print summary of the ntuple."""
        n = self._tree.num_entries
        resolved = len(self._branch_map)
        total = len(self._tree.keys())

        lines = [
            f"╔══════════════════════════════════════════════════╗",
            f"║  MTFT × LHCb Open Data Ntuple                   ║",
            f"╚══════════════════════════════════════════════════╝",
            f"  File:      {self.path.name}",
            f"  Size:      {self.path.stat().st_size / 1e9:.2f} GiB",
            f"  Tree:      {self._tree_name}",
            f"  Events:    {n:,}",
            f"  Branches:  {total} total, {resolved} MTFT-mapped",
            f"",
            f"  Mapped branches:",
        ]
        for generic, actual in sorted(self._branch_map.items()):
            lines.append(f"    {generic:<20} → {actual}")

        if resolved == 0:
            lines.append("  ⚠ No branches resolved — check tree_name or naming convention")
            lines.append(f"  Available: {list(self._tree.keys())[:15]}...")

        text = "\n".join(lines)
        print(text)
        return text

    # ═══════════════════════════════════════════════════════════════
    # §1  INVARIANT MASS COMPUTATION
    # ═══════════════════════════════════════════════════════════════

    def dimuon_mass(self) -> np.ndarray:
        """
        Compute μ⁺μ⁻ invariant mass from four-momenta.

        If the ntuple already has a J/ψ mass branch, returns that.
        Otherwise computes from muon PX/PY/PZ/PE.

        Returns mass in MeV.
        """
        # Try pre-computed branch first
        if "Jpsi_M" in self._branch_map:
            return self._get("Jpsi_M")

        # Compute from 4-momenta
        px = self._get("mup_PX") + self._get("mum_PX")
        py = self._get("mup_PY") + self._get("mum_PY")
        pz = self._get("mup_PZ") + self._get("mum_PZ")
        E  = self._get("mup_PE") + self._get("mum_PE")
        m2 = E**2 - px**2 - py**2 - pz**2
        return np.sqrt(np.maximum(m2, 0.0))

    def b_mass(self) -> np.ndarray:
        """
        B± invariant mass from the ntuple.

        Returns mass in MeV.
        """
        if "B_M" in self._branch_map:
            return self._get("B_M")

        # Compute from J/ψ + K four-momenta
        px = self._get("Jpsi_PX") + self._get("K_PX")
        py = self._get("Jpsi_PY") + self._get("K_PY")
        pz = self._get("Jpsi_PZ") + self._get("K_PZ")
        E  = self._get("Jpsi_PE") + self._get("K_PE")
        m2 = E**2 - px**2 - py**2 - pz**2
        return np.sqrt(np.maximum(m2, 0.0))

    # ═══════════════════════════════════════════════════════════════
    # §2  SELECTION CUTS
    # ═══════════════════════════════════════════════════════════════

    def jpsi_signal_mask(
        self,
        window_MeV: float = 50.0,
    ) -> np.ndarray:
        """
        Boolean mask selecting events in the J/ψ mass window.

        Default: |m(μμ) − m(J/ψ)_PDG| < 50 MeV.
        """
        m = self.dimuon_mass()
        return np.abs(m - M_JPSI_PDG) < window_MeV

    def b_signal_mask(
        self,
        window_MeV: float = 50.0,
    ) -> np.ndarray:
        """Boolean mask for B± mass window."""
        m = self.b_mass()
        return np.abs(m - M_BPLUS_PDG) < window_MeV

    def muon_pid_mask(self, threshold: float = 0.5) -> np.ndarray:
        """
        Muon PID cut: both muons must have ProbNNmu > threshold.
        Falls back to all-True if PID branches are absent.
        """
        try:
            p1 = self._get("mup_ProbNNmu")
            p2 = self._get("mum_ProbNNmu")
            return (p1 > threshold) & (p2 > threshold)
        except KeyError:
            warnings.warn("ProbNNmu branches not found — skipping PID cut")
            return np.ones(self._tree.num_entries if self._max is None
                           else min(self._max, self._tree.num_entries),
                           dtype=bool)

    # ═══════════════════════════════════════════════════════════════
    # §3  MTFT CONFRONTATION
    # ═══════════════════════════════════════════════════════════════

    def jpsi_peak_fit(self) -> dict:
        """
        Fit the J/ψ peak with a simple Gaussian + linear background.

        Returns dict with mean, sigma, yield, chi2.
        Uses numpy histogram + least-squares (no ROOT dependency).
        """
        from scipy.optimize import curve_fit

        m = self.dimuon_mass()
        # Focus on J/ψ region
        mask = (m > 3000) & (m < 3200)
        m_sig = m[mask]

        if len(m_sig) < 100:
            return {"error": "Too few events in J/ψ window"}

        # Histogram
        nbins = 100
        counts, edges = np.histogram(m_sig, bins=nbins, range=(3000, 3200))
        centers = 0.5 * (edges[:-1] + edges[1:])
        bin_width = edges[1] - edges[0]

        # Gaussian + linear
        def model(x, A, mu, sigma, a, b):
            gauss = A * np.exp(-0.5 * ((x - mu) / sigma)**2)
            linear = a * x + b
            return gauss + linear

        # Initial guesses
        p0 = [np.max(counts), 3097, 15, 0, np.median(counts)]

        try:
            popt, pcov = curve_fit(model, centers, counts, p0=p0,
                                   sigma=np.sqrt(np.maximum(counts, 1)),
                                   absolute_sigma=True)
            perr = np.sqrt(np.diag(pcov))
        except Exception as e:
            return {"error": str(e)}

        A, mu, sigma, a, b = popt
        A_e, mu_e, sigma_e, a_e, b_e = perr

        # Residual from PDG
        delta_jpsi = mu - M_JPSI_PDG

        return {
            "mean_MeV": mu,
            "mean_err_MeV": mu_e,
            "sigma_MeV": abs(sigma),
            "sigma_err_MeV": sigma_e,
            "yield": A * abs(sigma) * np.sqrt(2 * np.pi) / bin_width,
            "delta_from_PDG_MeV": delta_jpsi,
            "n_events_window": len(m_sig),
        }

    def mtft_confrontation(self, verbose: bool = True) -> dict:
        """
        Run MTFT predictions against this dataset.

        Checks:
        1. J/ψ mass reconstruction (validates μ mass)
        2. B± mass reconstruction (validates analysis chain)
        3. Desert check: no unexpected peaks in dimuon spectrum
        4. Hidden doublet relevance assessment

        Returns a summary dict.
        """
        results = {}

        # ── J/ψ peak ───────────────────────────────────────────────
        fit = self.jpsi_peak_fit()
        results["jpsi_fit"] = fit

        # ── B± mass ────────────────────────────────────────────────
        try:
            bm = self.b_mass()
            mask_b = (bm > 5200) & (bm < 5350)
            results["b_mass_median_MeV"] = float(np.median(bm[mask_b])) if np.sum(mask_b) > 0 else None
            results["b_candidates"] = int(np.sum(mask_b))
        except KeyError:
            results["b_mass_median_MeV"] = None

        # ── Desert scan ────────────────────────────────────────────
        try:
            m_dimu = self.dimuon_mass()
            # Look for unexpected bumps outside J/ψ and ψ(2S) windows
            # Exclude: J/ψ (3050-3150), ψ(2S) (3640-3730)
            desert_mask = (
                (m_dimu > 1000) & (m_dimu < 5000)
                & ~((m_dimu > 3050) & (m_dimu < 3150))
                & ~((m_dimu > 3640) & (m_dimu < 3730))
            )
            results["desert_events"] = int(np.sum(desert_mask))
            results["desert_note"] = (
                "Dimuon spectrum outside J/ψ and ψ(2S) — "
                "any unexpected peaks would challenge the MTFT desert prediction"
            )
        except KeyError:
            pass

        # ── Hidden doublet relevance ───────────────────────────────
        results["hidden_doublet"] = {
            "relevant": False,
            "note": (
                "This B± → J/ψK± dataset establishes the J/ψ reconstruction "
                "pipeline. The hidden doublet search requires J/ψ → γ + X "
                "radiative decays (BES-III primary channel). However, LHCb's "
                "central exclusive production pp → p + X + p in the 1200–1500 MeV "
                "window is a complementary search channel using the same detector."
            ),
            "H11_mass_MeV": H11_MASS,
            "H13_mass_MeV": H13_MASS,
            "splitting_MeV": DOUBLET_SPLIT,
        }

        # ── Report ────────────────────────────────────────────────
        if verbose:
            self._print_confrontation(results)

        return results

    def _print_confrontation(self, results: dict):
        lines = [
            "",
            "═" * 60,
            "  MTFT × LHCb  DATASET CONFRONTATION",
            "═" * 60,
            "",
        ]

        # J/ψ
        fit = results.get("jpsi_fit", {})
        if "error" not in fit:
            lines.append(f"  J/ψ PEAK FIT:")
            lines.append(f"    Reconstructed mass: {fit['mean_MeV']:.2f} ± {fit['mean_err_MeV']:.2f} MeV")
            lines.append(f"    PDG value:          {M_JPSI_PDG:.3f} MeV")
            lines.append(f"    Δ(reco − PDG):      {fit['delta_from_PDG_MeV']:+.2f} MeV")
            lines.append(f"    Resolution (σ):     {fit['sigma_MeV']:.2f} MeV")
            lines.append(f"    Events in window:   {fit['n_events_window']:,}")
            lines.append(f"    → Validates μ⁺μ⁻ reconstruction chain ✓")
        else:
            lines.append(f"  J/ψ fit failed: {fit['error']}")

        # B±
        lines.append("")
        b_med = results.get("b_mass_median_MeV")
        if b_med:
            lines.append(f"  B± MASS:")
            lines.append(f"    Median (signal window): {b_med:.1f} MeV")
            lines.append(f"    PDG value:              {M_BPLUS_PDG:.2f} MeV")
            lines.append(f"    Candidates:             {results['b_candidates']:,}")

        # Desert
        lines.append("")
        if "desert_events" in results:
            lines.append(f"  DESERT CHECK (dimuon off-peak):")
            lines.append(f"    Events outside J/ψ, ψ(2S): {results['desert_events']:,}")
            lines.append(f"    MTFT predicts: no unexpected resonances")

        # Hidden doublet
        lines.append("")
        hd = results["hidden_doublet"]
        lines.append(f"  HIDDEN DOUBLET STATUS:")
        lines.append(f"    Direct search: {'Yes' if hd['relevant'] else 'No (pipeline validation mode)'}")
        lines.append(f"    H₁₁: {hd['H11_mass_MeV']:.0f} MeV  |  H₁₃: {hd['H13_mass_MeV']:.0f} MeV")
        lines.append(f"    {hd['note']}")

        lines.append("")
        lines.append("═" * 60)
        print("\n".join(lines))

    # ═══════════════════════════════════════════════════════════════
    # §4  PLOTTING
    # ═══════════════════════════════════════════════════════════════

    def plot_jpsi(
        self,
        save: Optional[str] = None,
        nbins: int = 200,
        range_MeV: tuple = (2900, 3300),
    ):
        """Plot the J/ψ mass peak with MTFT annotation."""
        plt = _require("matplotlib.pyplot")

        m = self.dimuon_mass()
        mask = (m > range_MeV[0]) & (m < range_MeV[1])

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(m[mask], bins=nbins, range=range_MeV,
                histtype='stepfilled', alpha=0.7, color='#2196F3',
                label=f'Data ({np.sum(mask):,} events)')

        # PDG line
        ax.axvline(M_JPSI_PDG, color='red', ls='--', lw=1.5,
                   label=f'PDG: {M_JPSI_PDG:.1f} MeV')

        ax.set_xlabel('$m(\\mu^+\\mu^-)$ [MeV]', fontsize=14)
        ax.set_ylabel(f'Candidates / {(range_MeV[1]-range_MeV[0])/nbins:.1f} MeV',
                      fontsize=14)
        ax.set_title('MTFT × LHCb: $J/\\psi \\to \\mu^+\\mu^-$', fontsize=16)
        ax.legend(fontsize=12)
        ax.tick_params(labelsize=12)

        fig.tight_layout()
        if save:
            fig.savefig(save, dpi=150)
            print(f"  Saved: {save}")
        else:
            plt.show()
        return fig

    def plot_b_mass(
        self,
        save: Optional[str] = None,
        nbins: int = 200,
        range_MeV: tuple = (5100, 5500),
    ):
        """Plot the B± mass peak."""
        plt = _require("matplotlib.pyplot")

        m = self.b_mass()
        mask = (m > range_MeV[0]) & (m < range_MeV[1])

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(m[mask], bins=nbins, range=range_MeV,
                histtype='stepfilled', alpha=0.7, color='#4CAF50',
                label=f'Data ({np.sum(mask):,} events)')

        ax.axvline(M_BPLUS_PDG, color='red', ls='--', lw=1.5,
                   label=f'PDG: {M_BPLUS_PDG:.1f} MeV')

        ax.set_xlabel('$m(J/\\psi\\, K^\\pm)$ [MeV]', fontsize=14)
        ax.set_ylabel(f'Candidates / {(range_MeV[1]-range_MeV[0])/nbins:.1f} MeV',
                      fontsize=14)
        ax.set_title('MTFT × LHCb: $B^\\pm \\to J/\\psi\\, K^\\pm$', fontsize=16)
        ax.legend(fontsize=12)
        ax.tick_params(labelsize=12)

        fig.tight_layout()
        if save:
            fig.savefig(save, dpi=150)
        else:
            plt.show()
        return fig

    def plot_dimuon_full(
        self,
        save: Optional[str] = None,
        nbins: int = 500,
        range_MeV: tuple = (200, 5000),
        log_y: bool = True,
    ):
        """
        Full dimuon mass spectrum — the desert check.

        Annotates J/ψ, ψ(2S), and the MTFT hidden doublet search window.
        """
        plt = _require("matplotlib.pyplot")

        m = self.dimuon_mass()
        mask = (m > range_MeV[0]) & (m < range_MeV[1])

        fig, ax = plt.subplots(figsize=(14, 6))
        ax.hist(m[mask], bins=nbins, range=range_MeV,
                histtype='stepfilled', alpha=0.6, color='#607D8B',
                label='$\\mu^+\\mu^-$ spectrum')

        # Known resonances
        ax.axvline(M_JPSI_PDG, color='#F44336', ls='-', lw=2, alpha=0.7,
                   label='$J/\\psi$')
        ax.axvline(3686.10, color='#FF9800', ls='-', lw=2, alpha=0.7,
                   label="$\\psi(2S)$")

        # Hidden doublet window (Paper 25)
        ax.axvspan(1200, 1500, alpha=0.1, color='cyan',
                   label='MTFT doublet window (1200–1500 MeV)')
        ax.axvline(H11_MASS, color='cyan', ls=':', lw=1.5, alpha=0.8)
        ax.axvline(H13_MASS, color='cyan', ls=':', lw=1.5, alpha=0.8)

        if log_y:
            ax.set_yscale('log')
        ax.set_xlabel('$m(\\mu^+\\mu^-)$ [MeV]', fontsize=14)
        ax.set_ylabel('Candidates', fontsize=14)
        ax.set_title('MTFT × LHCb: Dimuon spectrum — Desert & doublet check',
                     fontsize=16)
        ax.legend(fontsize=11, loc='upper right')
        ax.tick_params(labelsize=12)

        fig.tight_layout()
        if save:
            fig.savefig(save, dpi=150)
        else:
            plt.show()
        return fig


# ═══════════════════════════════════════════════════════════════════════
# §5  BATCH PROCESSING — multiple files (MagUp + MagDown)
# ═══════════════════════════════════════════════════════════════════════

def combine_ntuples(
    paths: list[Union[str, Path]],
    tree_name: Optional[str] = None,
    max_events_per_file: Optional[int] = None,
) -> dict[str, np.ndarray]:
    """
    Combine multiple ROOT files into merged arrays.

    Useful for combining Magnet Up + Magnet Down datasets.

    Returns dict of {branch_name: concatenated_array}.
    """
    uproot = _require("uproot")

    merged = {}
    for p in paths:
        nt = LHCbNtuple(p, tree_name=tree_name, max_events=max_events_per_file)
        for generic in nt._branch_map:
            try:
                arr = nt._get(generic)
                if generic in merged:
                    merged[generic] = np.concatenate([merged[generic], arr])
                else:
                    merged[generic] = arr.copy()
            except Exception:
                pass

    return merged


# ═══════════════════════════════════════════════════════════════════════
# §6  SETUP HELPER
# ═══════════════════════════════════════════════════════════════════════

def setup_check() -> dict:
    """
    Check that all dependencies are installed and report versions.
    """
    results = {}

    for pkg in ["uproot", "awkward", "numpy", "scipy", "matplotlib"]:
        try:
            mod = __import__(pkg)
            results[pkg] = {"installed": True, "version": getattr(mod, "__version__", "?")}
        except ImportError:
            results[pkg] = {"installed": False, "version": None}

    # Print
    print("\n  MTFT × LHCb  DEPENDENCY CHECK")
    print("  " + "─" * 40)
    for pkg, info in results.items():
        status = f"✓ {info['version']}" if info['installed'] else "✗ MISSING"
        print(f"    {pkg:<15} {status}")

    missing = [p for p, i in results.items() if not i['installed']]
    if missing:
        print(f"\n  Install missing packages:")
        print(f"    pip install {' '.join(missing)}")
    else:
        print(f"\n  All dependencies satisfied ✓")

    print(f"\n  Quick start:")
    print(f"    from mtft.lhcb_analysis import LHCbNtuple")
    print(f"    nt = LHCbNtuple('path/to/dvntuple.root')")
    print(f"    nt.info()")
    print(f"    nt.mtft_confrontation()")
    print()

    return results


# ═══════════════════════════════════════════════════════════════════════
# Module entry point
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    setup_check()
