"""ZT-01: the twists of X₀(143) — the records' encoding, the cup form by route A, 143a1's pullbacks, the transfer and the flux
charge (records FF-14, FF-17, FF-18)."""
import random

import pytest

from mtft.lines import twists as T

# FF-17's cup form on the 26 basis twists, row i as a bitmask (record ff17_numbers.json, "form C (rows as 26-bit masks)").
FF17_C = [17826836, 29443092, 16778259, 58757808, 38406671, 253704, 896, 515912, 58757344, 29530360, 30054567, 58758832, 16781224,
          16384, 8036514, 58740648, 7995554, 29443760, 128, 29443088, 12664833, 4276224, 62604818, 22712074, 13278991, 4229400]


def _random_twists(n, seed=143):
    rng = random.Random(seed)
    return [rng.getrandbits(26) for _ in range(n)]


def test_the_records_integers():
    assert T.dim() == 26 and len(T.basis()) == 26
    assert T.W() == (0, 6665056, 14025081, 11776537)
    assert T.pullback(1, 0) == T.NAMED["φ*η"]
    assert T.act("W13", T.NAMED["minimum"]) == T.NAMED["its W13 image"]
    assert T.act("W11", T.NAMED["minimum"]) == T.NAMED["minimum"]


def test_encoding_round_trips_and_evaluates():
    G = T.graph()
    for z in _random_twists(200):
        x = T.to_character(z)
        assert T.to_twist(x) == z
        for f in range(len(G.faces)):                       # even on every cusp loop
            assert sum(c for j, c in enumerate(G.face_cycle(f)) if (x >> j) & 1) % 2 == 0
    from mtft.lines import walls as WL
    rng = random.Random(11)
    for L, ell, cls in WL.primitive_classes(6.0)[:60]:
        z = rng.getrandbits(26)
        direct = sum(c for j, c in enumerate(cls) if (T.to_character(z) >> j) & 1) % 2
        assert T.value(z, cls) == direct
    odd = next(1 << j for j in range(G.b1) if any(G.face_cycle(f)[j] % 2 for f in range(len(G.faces))))
    with pytest.raises(ValueError):
        T.to_twist(odd)


def test_the_involutions_on_twists():
    for i in range(26):
        e = 1 << i
        for g in ("W11", "W13", "W143", "mirror"):
            assert T.act(g, T.act(g, e)) == e
        assert T.act("W143", e) == T.act("W11", T.act("W13", e))
    assert T.act("mirror", T.pullback(1, 0)) == T.pullback(1, 0)
    assert T.act("mirror", T.pullback(0, 1)) == T.pullback(1, 1)


def test_cup_form_is_ff17s_and_symplectic():
    assert list(T.cup_form()) == FF17_C
    chk = T.check_cup_form()
    assert chk["alternating"] and chk["rank"] == 26 and chk["face words rotated give the same form"]
    assert chk["invariant under"] == {"W11": True, "W13": True, "W143": True, "mirror": True}
    for z in _random_twists(50, seed=7):
        assert T.form(z, z) == 0
        assert T.dual_class(T.wall_twist(z)) == z and T.wall_twist(T.dual_class(z)) == z
    ys = _random_twists(20, seed=8)
    assert all(T.intersection(a, b) == T.intersection(b, a) for a in ys for b in ys)


def test_transfer_and_flux_charge():
    assert T.check_transfer() == {"N(ζ) = φ*(q₂, q₁) on the basis": True, "W143 = W11 W13 on twists": True, "W isotropic": True}
    assert all(T.flux_charge(z) == (0, 0) for z in T.NAMED.values())                      # FF-18: q of the named twists
    assert T.flux_charge(14841521) == (0, 1) and T.flux_charge(15021403) == (1, 0) and T.flux_charge(57549835) == (1, 1)  # FF-18's orbit starts
    star = {(0, 0): 0, (1, 0): T.pullback(1, 0), (0, 1): T.pullback(0, 1), (1, 1): T.pullback(1, 1)}
    for z in _random_twists(100, seed=9):
        q = T.flux_charge(z)
        assert T.transfer(z) == star[(q[1], q[0])]
        if T.act("W11", z) == z or T.act("W13", z) == z or T.act("W143", z) == z:
            assert q == (0, 0)
    m = T.mirror_on_flux_classes()                          # FF-18: (0,1) fixed, (1,0) and (1,1) swapped
    assert m["exact on the basis"] and m["map"] == {(0, 1): (0, 1), (1, 0): (1, 1), (1, 1): (1, 0)}


def test_f2_helpers():
    assert T.rank_f2([0b011, 0b110, 0b101]) == 2 and T.parity(0b1011) == 1
    assert T.inverse_rows((0b01, 0b11)) == (0b01, 0b11)
    with pytest.raises(ValueError):
        T.inverse_rows((0b01, 0b01))
