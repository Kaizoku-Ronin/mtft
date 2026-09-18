"""Single source of truth for the release version pin used by the version-gate tests.

A release bump touches exactly one test line (below) in addition to the three
package-side locations (pyproject.toml, src/mtft/__init__.py, CITATION.cff);
``test_release_pin.py`` verifies all four agree.  Introduced in v0.26.1 after
four consecutive releases shipped with stale per-file pins (Kimi).
"""
RELEASE_VERSION = "0.32.0"
