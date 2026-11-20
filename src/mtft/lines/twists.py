r"""mtft.lines.twists — the twists of X₀(143) mod 2: the records' encoding, the cup form, 143a1's pullbacks, the flux charge.

EXACT; ZT-01, 2026-11-20.

A twist of the fermions along the curve is a spin structure relative to a reference: a character
ζ ∈ H¹(X₀(N), ℤ/2), one sign for each class of closed curve.  Bosons do not see it (MATTER-04).

Encoding (the records' 26-bit integers, FF-14 to FF-20).  The punctured curve's classes are ℤ²⁹ on the 29
generators of ``darts`` (the non-tree edges), so a character of it is 29 bits.  The closed curve's characters
are those even on every cusp loop (the face cycles): 26 dimensions at N = 143.  ``basis`` is FF-14's basis of
them, found from the face cycles by elimination; a twist is the integer whose bit i is the coefficient of
basis[i].  ``class_bits`` turns a class into 26 bits Y with ζ(γ) = parity(ζ & Y).

Layers:

  basis, to_twist, to_character, class_bits, value   the encoding
  act                                                 W_Q and the mirror on twists, from the class matrices of ``involutions``
  cup, cup_form, form                                 the cup product mod 2 by route A (FF-17): cone each cusp face from a
                                                      centre, halve each edge, order centres < midpoints < vertices and apply
                                                      Alexander–Whitney, giving
                                                      (α ∪ β)[X̄] = Σ_faces Σ_{j<i} α(e_j) β(e_i) + Σ_e α(e) β(e)
  pullback, W                                         143a1's twists on the curve: φ*(a, b)(γ) = a n₁ + b n₂ mod 2, with
                                                      (n₁, n₂) the class's 143a1 coordinates (``periods``); φ*η = φ*(1, 0)
  transfer, flux_charge                               N(ζ) = ζ + W₁₁ζ + W₁₃ζ + W₁₄₃ζ = (1 + W₁₁)(1 + W₁₃)ζ = φ*φ_!ζ (FF-18);
                                                      q(ζ) = (ζ ∪ φ*(1,0), ζ ∪ φ*(0,1)), kept by every change of the curve's
                                                      shape that keeps the γ-flux, and zero on every twist that an
                                                      Atkin–Lehner involution fixes

The cup form is the Poincaré duality of the closed curve mod 2: the class dual to a twist δ has bits Cδ
(``dual_class``), the intersection of two classes Y₁, Y₂ is (C⁻¹Y₁)·Y₂ (``class_form``, ``intersection``), and a
ζ-changing wall along a curve of class Y changes ζ by C⁻¹Y (``wall_twist``).

These twists are the ζ of MATTER-04 (``mtft.matter.spin_twists``): 2-torsion line bundles on X₀(143), J₀(143)[2] = (ℤ/2)²⁶.
The encoding works at any level ``darts`` handles; the involutions, the pullbacks, the transfer and the flux charge are
at N = 143.

Records: FF-14 (the encoding and the pullbacks), FF-17 (route A), FF-18 (the transfer and the flux charge), in the
MTFT Python project under claude/experiments/.  What this is not: the twist energy, or which twist a vacuum
takes (FF-15; not in the package).
"""
from __future__ import annotations

from functools import lru_cache

from . import darts as _d
from . import involutions as _inv
from . import periods as _per

# The records' names for some twists (MTFT Python project, FF-14 to FF-20).  "minimum" is the energy minimum of
# FF-15's one-loop landscape at r = 230; W₁₁ fixes it and W₁₃ sends it to "its W13 image".
NAMED = {"O": 0, "φ*η": 6665056, "φ*(0,1)": 14025081, "φ*(1,1)": 11776537, "minimum": 45796554, "its W13 image": 48201642,
         "lowest symmetric": 36661717, "FF-14 weakest": 394080, "φ1·φA": 61104341, "φA·φB": 39774856, "φ1·φB": 33196637}


def parity(x: int) -> int:
    return bin(x).count("1") & 1


def apply_rows(rows, x: int) -> int:
    """y = M x over F2, with M given as row bitmasks."""
    y = 0
    for i, r in enumerate(rows):
        if parity(r & x):
            y |= 1 << i
    return y


@lru_cache(maxsize=4)
def graph(N: int = 143):
    return _d.graph(N)


def _mask(row) -> int:
    return sum(((int(x) & 1) << j) for j, x in enumerate(row))


@lru_cache(maxsize=4)
def basis(N: int = 143) -> tuple:
    """FF-14's basis of the characters (29 bits) that are even on every cusp loop."""
    G = graph(N)
    piv, R = [], []
    for row in (G.face_cycle(f) for f in range(len(G.faces))):
        v = _mask(row)
        for p, rv in zip(piv, R):
            if (v >> p) & 1:
                v ^= rv
        if v:
            p = (v & -v).bit_length() - 1
            for k in range(len(R)):
                if (R[k] >> p) & 1:
                    R[k] ^= v
            piv.append(p)
            R.append(v)
    out = []
    for f in (j for j in range(G.b1) if j not in piv):
        x = 1 << f
        for p, rv in zip(piv, R):
            if (rv >> f) & 1:
                x |= 1 << p
        out.append(x)
    return tuple(out)


def dim(N: int = 143) -> int:
    return len(basis(N))


@lru_cache(maxsize=4)
def _pivots(N: int):
    piv = {}
    for i, v in enumerate(basis(N)):
        comb = 1 << i
        for p in sorted(piv, reverse=True):
            if (v >> p) & 1:
                v ^= piv[p][0]
                comb ^= piv[p][1]
        piv[v.bit_length() - 1] = (v, comb)
    return tuple(sorted(piv.items(), reverse=True))


def to_twist(x29: int, N: int = 143) -> int:
    """The twist (bits on ``basis``) of a 29-bit character; ValueError if it is odd on a cusp loop."""
    comb = 0
    for p, (v, c) in _pivots(N):
        if (x29 >> p) & 1:
            x29 ^= v
            comb ^= c
    if x29:
        raise ValueError("not a character of the closed curve")
    return comb


def to_character(z: int, N: int = 143) -> int:
    """The 29-bit character of a twist."""
    x = 0
    for i, v in enumerate(basis(N)):
        if (z >> i) & 1:
            x ^= v
    return x


def class_bits(cls, N: int = 143) -> int:
    """The 26 bits Y of a class of ℤ²⁹ (Y_i = basis[i] · cls mod 2), so that ζ(γ) = parity(ζ & Y)."""
    c = _mask(cls)
    return sum(parity(v & c) << i for i, v in enumerate(basis(N)))


def value(z: int, cls, N: int = 143) -> int:
    """ζ(γ) ∈ {0, 1}: 1 where the twist makes the fermions antiperiodic round the class."""
    return parity(z & class_bits(cls, N))


# ------------------------------------------------------------------ the involutions on twists
@lru_cache(maxsize=8)
def matrix(name: str) -> tuple:
    """Class matrices of W11, W13, W143 and the mirror (row e = the class of the image of generator e)."""
    if name == "mirror":
        return _inv.mirror_class_matrix()
    return _inv.al_matrix(int(name[1:]))


def act(name: str, z: int, N: int = 143) -> int:
    """The twist ζ∘g for g = W11, W13, W143 or the mirror."""
    M = matrix(name)
    x = to_character(z, N)
    out = 0
    for e, row in enumerate(M):
        s = 0
        for j, m in enumerate(row):
            if m & 1 and (x >> j) & 1:
                s ^= 1
        if s:
            out |= 1 << e
    return to_twist(out, N)


# ------------------------------------------------------------------ the cup form by route A (FF-17)
@lru_cache(maxsize=4)
def _face_words(N: int):
    G = graph(N)
    return tuple(tuple(G.eof[y] for y in o) for o in G.faces)


@lru_cache(maxsize=4)
def _nontree(N: int):
    return {k: j for j, k in enumerate(graph(N).nontree)}


def cup(x29: int, y29: int, N: int = 143, face_words=None) -> int:
    """(α ∪ β)[X̄] mod 2 for 29-bit characters, through the cusp faces."""
    NT = _nontree(N)
    s = parity(x29 & y29)
    for w in (face_words if face_words is not None else _face_words(N)):
        pre = 0
        for e in w:
            j = NT.get(e)
            if j is None:
                continue
            if pre and (y29 >> j) & 1:
                s ^= 1
            if (x29 >> j) & 1:
                pre ^= 1
    return s


@lru_cache(maxsize=4)
def cup_form(N: int = 143) -> tuple:
    """The 26 × 26 cup form on twists, as row bitmasks (row i, bit j = basis[i] ∪ basis[j])."""
    B = basis(N)
    return tuple(sum(cup(B[i], B[j], N) << j for j in range(len(B))) for i in range(len(B)))


def form(z1: int, z2: int, N: int = 143) -> int:
    return parity(apply_rows(cup_form(N), z1) & z2)


def dual_class(z: int, N: int = 143) -> int:
    """The 26 class bits of the Poincaré dual of a twist: a wall along it changes ζ by z."""
    return apply_rows(cup_form(N), z)


def inverse_rows(rows) -> tuple:
    """The inverse over F2 of a square matrix given as row bitmasks; ValueError if it is singular."""
    n = len(rows)
    aug = [(rows[i], 1 << i) for i in range(n)]
    for c in range(n):
        r = next((i for i in range(c, n) if (aug[i][0] >> c) & 1), None)
        if r is None:
            raise ValueError("singular over F2")
        aug[c], aug[r] = aug[r], aug[c]
        for i in range(n):
            if i != c and (aug[i][0] >> c) & 1:
                aug[i] = (aug[i][0] ^ aug[c][0], aug[i][1] ^ aug[c][1])
    return tuple(a[1] for a in aug)


@lru_cache(maxsize=4)
def class_form(N: int = 143) -> tuple:
    """C⁻¹: the intersection form mod 2 on class bits."""
    return inverse_rows(cup_form(N))


def intersection(y1: int, y2: int, N: int = 143) -> int:
    """The intersection number mod 2 of two classes given by their 26 bits."""
    return parity(apply_rows(class_form(N), y1) & y2)


def wall_twist(y: int, N: int = 143) -> int:
    """The twist that a wall along a curve of class Y adds: C⁻¹Y."""
    return apply_rows(class_form(N), y)


def rank_f2(vecs) -> int:
    piv = {}
    for v in vecs:
        v = int(v)
        while v:
            p = v.bit_length() - 1
            if p in piv:
                v ^= piv[p]
            else:
                piv[p] = v
                break
    return len(piv)


def check_cup_form(N: int = 143) -> dict:
    """Alternating, rank, invariance under the involutions, and independence of where each face word starts."""
    C = cup_form(N)
    n = len(C)
    alternating = all(((C[i] >> j) & 1) == ((C[j] >> i) & 1) for i in range(n) for j in range(n)) and all(not (C[i] >> i) & 1 for i in range(n))
    inv = {}
    if N == 143:
        for name in ("W11", "W13", "W143", "mirror"):
            imgs = [act(name, 1 << i) for i in range(n)]
            inv[name] = all(form(imgs[i], imgs[j]) == ((C[i] >> j) & 1) for i in range(n) for j in range(i, n))
    rotated = tuple(w[len(w) // 3:] + w[:len(w) // 3] for w in _face_words(N))
    B = basis(N)
    same = all(cup(B[i], B[j], N, rotated) == ((C[i] >> j) & 1) for i in range(n) for j in range(n))
    return {"alternating": alternating, "rank": rank_f2(C), "invariant under": inv, "face words rotated give the same form": same}


# ------------------------------------------------------------------ 143a1's pullbacks, the transfer and the flux charge (FF-14, FF-18)
@lru_cache(maxsize=1)
def ell_matrix() -> tuple:
    """The 143a1 coordinates (n₁, n₂) of the 29 generators, as two rows: φ_* on ℤ²⁹."""
    b1 = graph(143).b1
    cols = [_per.coordinates(tuple(1 if k == e else 0 for k in range(b1)))["ell"] for e in range(b1)]
    return (tuple(int(c[0]) for c in cols), tuple(int(c[1]) for c in cols))


@lru_cache(maxsize=4)
def pullback(a: int, b: int) -> int:
    """φ*(a, b): the twist γ ↦ a n₁(γ) + b n₂(γ) mod 2, with (n₁, n₂) the class's 143a1 coordinates."""
    n1, n2 = ell_matrix()
    return to_twist(sum(((a * u + b * v) % 2) << e for e, (u, v) in enumerate(zip(n1, n2))))


def W() -> tuple:
    """The four pullbacks: O, φ*(1,0) = φ*η, φ*(0,1), φ*(1,1)."""
    return (0, pullback(1, 0), pullback(0, 1), pullback(1, 1))


def transfer(z: int) -> int:
    """N(ζ) = ζ + W₁₁ζ + W₁₃ζ + W₁₄₃ζ, which equals φ*(q₂, q₁)."""
    return z ^ act("W11", z) ^ act("W13", z) ^ act("W143", z)


def flux_charge(z: int) -> tuple:
    """q(ζ) = (ζ ∪ φ*(1,0), ζ ∪ φ*(0,1)) ∈ F2²."""
    return (form(z, pullback(1, 0)), form(z, pullback(0, 1)))


def check_transfer() -> dict:
    """The transfer identity on all basis vectors (linear, so exact), and W₁₄₃ = W₁₁W₁₃ on them."""
    star = {(0, 0): 0, (1, 0): pullback(1, 0), (0, 1): pullback(0, 1), (1, 1): pullback(1, 1)}
    n = dim(143)
    ident = all(transfer(1 << i) == star[(flux_charge(1 << i)[1], flux_charge(1 << i)[0])] for i in range(n))
    law = all(act("W143", 1 << i) == act("W11", act("W13", 1 << i)) for i in range(n))
    isotropic = all(form(a, b) == 0 for a in W() for b in W())
    return {"N(ζ) = φ*(q₂, q₁) on the basis": ident, "W143 = W11 W13 on twists": law, "W isotropic": isotropic}


def mirror_on_flux_classes() -> dict:
    """How the mirror moves the flux classes: q ↦ q(mirror ζ), read off the basis (both sides linear)."""
    n = dim(143)
    img = {}
    for k in ((1, 0), (0, 1)):
        src = next(i for i in range(n) if flux_charge(1 << i) == k)
        img[k] = flux_charge(act("mirror", 1 << src))
    def A(k):
        return ((k[0] * img[(1, 0)][0] + k[1] * img[(0, 1)][0]) % 2, (k[0] * img[(1, 0)][1] + k[1] * img[(0, 1)][1]) % 2)
    exact = all(flux_charge(act("mirror", 1 << i)) == A(flux_charge(1 << i)) for i in range(n))
    return {"map": {k: A(k) for k in ((0, 1), (1, 0), (1, 1))}, "exact on the basis": exact}
