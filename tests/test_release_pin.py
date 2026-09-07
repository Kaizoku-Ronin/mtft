"""Four-way release version consistency: pyproject.toml, __init__, CITATION.cff, tests/_release_pin.py."""
import re
from pathlib import Path

import mtft
from _release_pin import RELEASE_VERSION

ROOT = Path(__file__).resolve().parents[1]


def _pyproject_version() -> str:
    return re.search(r'^version\s*=\s*"([^"]+)"', (ROOT / "pyproject.toml").read_text(), re.M).group(1)


def _citation_version() -> str:
    return re.search(r"^version:\s*(\S+)", (ROOT / "CITATION.cff").read_text(), re.M).group(1)


def test_release_version_four_way():
    assert mtft.__version__ == RELEASE_VERSION
    assert _pyproject_version() == RELEASE_VERSION
    assert _citation_version() == RELEASE_VERSION
