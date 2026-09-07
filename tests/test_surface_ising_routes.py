"""v0.26.1: three-route Ising sum rule (elimination / spin-structure Pfaffians / brute force)."""
import math
import pytest

from mtft.surface import ising as I

ASTRA_143_A_AT_TC = 14.751840885840   # independent engine (Astra, 2026-09-05) and this package agree


@pytest.mark.parametrize("N", [6, 11, 15, 35, 55])
def test_sum_rule_three_routes(N):
    r = I.sum_rule_gate(N)
    assert r["status"] == "PASS", [g for g in r["gates"] if g["status"] != "PASS"]


def test_density_of_states_143_exact_integers():
    D = I.density_of_states(143)
    assert sum(D) == 1 << 56 and D[0] == 2 and D[75] == 2 and all(x == 0 for x in D[76:])
    A = I.even_subgraph_polynomial(D, 56, 84)
    assert A[0] == 1 and A[1] == 1 and sum(A) == 1 << 29 and A[:11] == [1, 1, 0, 4, 4, 0, 10, 18, 28, 76, 124]
    At = sum(a * (1 / math.sqrt(3)) ** j for j, a in enumerate(A))
    assert abs(At - ASTRA_143_A_AT_TC) < 1e-9


def test_genus_zero_and_one_share_a_graph():
    # N=6 (genus 0) and N=11 (genus 1) have the same dual graph: scalar thermodynamics cannot see genus.
    assert I.density_of_states(6) == I.density_of_states(11) == [2, 2, 2, 6, 4, 0, 0]
