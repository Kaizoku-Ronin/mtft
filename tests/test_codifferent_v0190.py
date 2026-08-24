"""Canonical Codifferent Theorem at instance: pure-Python re-verification."""
from mtft import codifferent as CD


def test_trace_identification_f2():
    assert CD.verify_orbit("f2")


def test_trace_identification_f3():
    assert CD.verify_orbit("f3")


def test_orbit_indices():
    assert CD.orbit_indices("f2") == {"index_OK": 576,
                                      "index_codiff": 576 * 1957}
    assert CD.orbit_indices("f3") == {"index_OK": 2304,
                                      "index_codiff": 2304 * 194616205}
