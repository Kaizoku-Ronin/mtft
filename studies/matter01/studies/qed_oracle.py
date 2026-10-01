"""Compatibility re-export of the independent analytic package oracle.

No engine gamma, spinor, amplitude, or kinematic helper is imported.
Run verify_qed_oracle.py for the reproducible independent comparison.
"""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT / "repo_candidate" if (ROOT / "repo_candidate").exists() else ROOT.parents[1]
sys.path.insert(0, str(REPO / "src"))
from mtft.matter.oracle import *  # noqa: F401,F403
