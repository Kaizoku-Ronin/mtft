"""KN-1: genus characters and the cusp windings of the trace-4 lines (the systole of X₀(143))."""
import pytest

from mtft.lines import genus as G
from mtft.lines import rademacher as R


def _euler(a, p):
    v = pow(a % p, (p - 1) // 2, p)
    return 0 if a % p == 0 else (1 if v == 1 else -1)


def test_kronecker_against_euler_and_values():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 101, 143 // 11):
        for a in range(-30, 31):
            assert G.kronecker(a, p) == _euler(a, p)
    assert G.kronecker(-3, 2) == -1 and G.kronecker(5, 2) == -1 and G.kronecker(7, 2) == 1 and G.kronecker(4, 2) == 0
    assert G.jacobi(2, 15) == 1 and G.jacobi(7, 15) == -1


def test_prime_discriminants_and_genus_vectors():
    assert G.prime_discriminants(12) == (-3, -4)
    assert G.prime_discriminants(21) == (-3, -7)
    assert G.prime_discriminants(-24) == (-3, 8)
    assert G.prime_discriminants(-4) == (-4,) and G.prime_discriminants(5) == (5,) and G.prime_discriminants(8) == (8,)
    assert G.genus_vector(12, 11) == (-1, -1) and G.genus_vector(12, 13) == (1, 1) and G.genus_vector(12, 143) == (-1, -1)
    for p in (11, 13, 23, 37, 47, 59, 61, 71, 73):          # split primes: both genus characters agree
        a, b = G.genus_vector(12, p)
        assert a == b == G.chi_minus3(p)


def test_genus_law_against_dedekind_sums_at_split_primes():
    tot = 0
    for p in (11, 13, 23, 37, 47, 59, 61, 71, 73, 83, 97):
        res = G.check_trace4(p, kmax=4, amax=300)
        assert res["elements"] > 0 and res["disagreements"] == 0
        assert res["ratios"][p] == (1 if p % 3 == 1 else -1)
        tot += res["elements"]
    assert tot > 300


def test_the_systole_of_143():
    res = G.check_trace4(143, kmax=4, amax=600)
    assert res["elements"] == 62 and res["disagreements"] == 0
    assert (res["ratios"][1], res["ratios"][11], res["ratios"][13], res["ratios"][143]) == (1, -1, 1, -1)
    assert G.trace4_windings(1) == (0, 1, 0) and G.trace4_windings(-1) == (0, -1, 0)
    seen = set()
    for g in G.trace4_elements(143, 4, 600):
        psi = R.Psi(g)
        assert R.windings(g) == G.trace4_windings(psi)
        seen.add(psi)
    assert seen == {1, -1}


def test_no_trace_four_elements_where_three_is_not_a_square():
    for p in (5, 7, 17, 19, 29, 31, 41, 43):
        assert G.trace4_elements(p, 4, 300) == []


def test_raise_paths():
    with pytest.raises(ValueError):
        G.prime_discriminants(48)
    with pytest.raises(ValueError):
        G.genus_vector(12, 3)
    with pytest.raises(ValueError):
        G.trace4_ratio(6)
    with pytest.raises(ValueError):
        G.trace4_ratio(0)
    with pytest.raises(ValueError):
        G.jacobi(3, 4)
    with pytest.raises(ValueError):
        G.kronecker(3, 0)
    with pytest.raises(ValueError):
        G.trace4_elements(0)
    with pytest.raises(ArithmeticError):
        G.trace4_windings(1, units={"u": {1: 1, 11: -1}}, order=("u",))
