r"""mtft.plane.quantities — the lattice of physical quantities, its units and what lives at a site (EXACT; v1.1.0, 2026-11-16).

A quantity has a dimension M^a L^b T^c I^d, a point v = (a, b, c, d) of ℤ⁴ (mass, length, time,
current).  Products add and quotients subtract, so the quantities of mechanics and electromagnetism
sit on a lattice, and a system of natural units is a sublattice that is declared to be numbers:

    c = (0, 1, −1, 0),   ħ = (1, 2, −1, 0),   ε₀ = (−1, −3, 4, 2),   G = (−1, 3, −2, 0).

Setting c = ħ = ε₀ = 1 leaves ℤ⁴/⟨c, ħ, ε₀⟩ ≅ ℤ ⊕ ℤ/2 (Smith invariants 1, 1, 2): every quantity is
read as a *mass dimension* Δ(v) = a − b − c + d and a *charge parity* d mod 2, and ``natural`` is
that map; its kernel is exactly the span of the three units.  Adding G gives the Planck sublattice,
of index 4 in ℤ⁴ with Smith invariants 1, 1, 2, 2: in Planck units four classes remain, by the
parities of the mass dimension and of the charge.  Both statements are computed here through
``mtft.integral_lattice`` (``lattice``), not asserted.

``planck_decompose(v)`` writes v = Δ·(mass) + ε·(charge) + n_c c + n_ħ ħ + n_ε ε₀ with integer
coefficients.  ``site(v)`` is the lab's prospector: the catalogued quantities at v, the products
and quotients of two catalogued quantities that land there, and the quantities in the same class
of natural units.  The catalogue is the lab's: Sphereland 25's list (``CATALOGUE``), an extension
(``EXTENSION``) and the lab's own dials and readouts (``LAB``).

A plane sees none of this: a coordinate on the upper half-plane is a ratio of two quantities at
one site, and a change of units multiplies every quantity at a site by one factor, so it moves no
ratio (``rescale``).

Provenance: the Modular Plane lab's core (version 8.3), Room 7.
"""
from __future__ import annotations

from functools import lru_cache

AXES = ("M", "L", "T", "I")
UNITS = {"c": (0, 1, -1, 0), "hbar": (1, 2, -1, 0), "eps0": (-1, -3, 4, 2), "G": (-1, 3, -2, 0)}
"""The dimensions of the speed of light, Planck's constant, the permittivity of the vacuum and Newton's constant."""

MASS, CHARGE = (1, 0, 0, 0), (0, 0, 1, 1)

CATALOGUE = (
    ("absity", (0, 1, 2, 0), "ladder", "B"), ("absement", (0, 1, 1, 0), "ladder", "A"), ("position", (0, 1, 0, 0), "ladder", "x"),
    ("velocity", (0, 1, -1, 0), "ladder", "v"), ("acceleration", (0, 1, -2, 0), "ladder", "a"), ("jerk", (0, 1, -3, 0), "ladder", "j"),
    ("snap", (0, 1, -4, 0), "ladder", "s"), ("time", (0, 0, 1, 0), "mech", "t"), ("frequency · decay constant", (0, 0, -1, 0), "mech", "ν, λ"),
    ("area · cross section", (0, 2, 0, 0), "mech", "σ"), ("volume", (0, 3, 0, 0), "mech", "Vol"), ("mass", (1, 0, 0, 0), "mech", "m"),
    ("momentum · impulse", (1, 1, -1, 0), "mech", "p"), ("force", (1, 1, -2, 0), "mech", "F"), ("energy · work · torque", (1, 2, -2, 0), "mech", "E"),
    ("action · angular momentum · ħ", (1, 2, -1, 0), "mech", "S, ħ"), ("power", (1, 2, -3, 0), "mech", "P"),
    ("moment of inertia", (1, 2, 0, 0), "mech", "I_m"), ("density", (1, -3, 0, 0), "mech", "ρ"),
    ("pressure · momentum flux · energy density", (1, -1, -2, 0), "mech", "Π"), ("energy flux (Poynting)", (1, 0, -3, 0), "mech", "S_P"),
    ("mass flux", (1, -2, -1, 0), "mech", "j_m"), ("gravitational constant", (-1, 3, -2, 0), "const", "G"), ("charge", (0, 0, 1, 1), "em", "q"),
    ("current", (0, 0, 0, 1), "em", "I"), ("current density", (0, -2, 0, 1), "em", "J"), ("electric field", (1, 1, -3, -1), "em", "E-field"),
    ("magnetic field", (1, 0, -2, -1), "em", "B-field"), ("potential φ", (1, 2, -3, -1), "em", "φ"), ("magnetic flux", (1, 2, -2, -1), "em", "Φ_B"),
    ("electric flux", (1, 3, -3, -1), "em", "Φ_E"), ("capacitance", (-1, -2, 4, 2), "em", "C"), ("inductance", (1, 2, -2, -2), "em", "L"),
)
"""Sphereland 25's catalogue: (name, dimension, family, symbol)."""

EXTENSION = (
    ("wavenumber", (0, -1, 0, 0), "mech", "k"), ("linear mass density", (1, -1, 0, 0), "mech", "λ_m"),
    ("surface tension · spring constant", (1, 0, -2, 0), "mech", "γ_s"), ("dynamic viscosity", (1, -1, -1, 0), "mech", "μ_v"),
    ("kinematic viscosity · diffusivity", (0, 2, -1, 0), "mech", "ν_k"), ("power per volume", (1, -1, -3, 0), "mech", "P/V"),
    ("intensity · power per area", (1, 0, -3, 0), "mech", "I_w"), ("specific energy · gravitational potential", (0, 2, -2, 0), "mech", "E/m"),
    ("specific volume", (-1, 3, 0, 0), "mech", "V/m"), ("surface density", (1, -2, 0, 0), "mech", "σ_m"),
    ("speed of light c", (0, 1, -1, 0), "const", "c"), ("Planck constant h", (1, 2, -1, 0), "const", "h"),
    ("permittivity ε₀", (-1, -3, 4, 2), "const", "ε₀"), ("permeability μ₀", (1, 1, -2, -2), "const", "μ₀"),
    ("resistance · impedance", (1, 2, -3, -2), "em", "R"), ("conductance", (-1, -2, 3, 2), "em", "G_e"), ("resistivity", (1, 3, -3, -2), "em", "ρ_e"),
    ("conductivity", (-1, -3, 3, 2), "em", "σ_e"), ("charge density", (0, -3, 1, 1), "em", "ρ_q"),
    ("electric displacement · polarization", (0, -2, 1, 1), "em", "D"), ("magnetizing field H · magnetization", (0, -1, 0, 1), "em", "H"),
    ("electric dipole moment", (0, 1, 1, 1), "em", "p_e"), ("magnetic dipole moment", (0, 2, 0, 1), "em", "μ_m"),
    ("vector potential", (1, 1, -2, -1), "em", "A_μ"), ("magnetic flux density gradient", (1, -1, -2, -1), "em", "∇B"),
    ("elementary charge e", (0, 0, 1, 1), "const", "e"),
)
"""The lab's extension of the catalogue."""

LAB = (
    ("κ, the sign dial", (0, 0, -2, 0), "lab", "κ"), ("α, the dilation dial", (0, 0, -1, 0), "lab", "α"), ("t, the time dial", (0, 0, 1, 0), "lab", "t"),
    ("the entry M₀₁ of M(t)", (0, 0, 1, 0), "lab", "M₀₁"), ("the entry M₁₀ of M(t)", (0, 0, -1, 0), "lab", "M₁₀"),
    ("ω, a mode's frequency", (0, 0, -1, 0), "lab", "ω"), ("⟨x²⟩ of a mode", (0, 2, 0, 0), "lab", "σ₁₁"), ("⟨p²⟩ of a mode", (2, 2, -2, 0), "lab", "σ₂₂"),
    ("⟨xp⟩, and √det σ", (1, 2, -1, 0), "lab", "σ₁₂"), ("τ of a mode's state, (σ₁₂ + i√det σ)/σ₂₂", (-1, 0, 1, 0), "lab", "τ_σ"),
    ("mω of a mode", (1, 0, -1, 0), "lab", "mω"), ("β, an inverse temperature (k_B = 1: a time with ħ = 1)", (0, 0, 1, 0), "lab", "β"),
    ("the particle number N, d, ψ, D, the entries of a motion", (0, 0, 0, 0), "lab", "N, d, ψ"),
)
"""The lab's own quantities: its dials and what its instruments read."""


def _site(v) -> tuple:
    v = tuple(int(x) for x in v)
    if len(v) != 4:
        raise ValueError("a site has four exponents: mass, length, time, current")
    return v


def mass_dimension(v) -> int:
    """Δ(v) = a − b − c + d: the power of mass a quantity becomes when c = ħ = ε₀ = 1."""
    a, b, c, d = _site(v)
    return a - b - c + d


def charge_parity(v) -> int:
    """d mod 2: whether a quantity carries an odd power of charge."""
    return _site(v)[3] % 2


def natural(v) -> tuple:
    """(mass dimension, charge parity): the class of a site in ℤ⁴/⟨c, ħ, ε₀⟩ ≅ ℤ ⊕ ℤ/2."""
    return mass_dimension(v), charge_parity(v)


def planck_class(v) -> tuple:
    """(mass dimension mod 2, charge parity): the class of a site in ℤ⁴/⟨c, ħ, ε₀, G⟩ ≅ ℤ/2 ⊕ ℤ/2."""
    return mass_dimension(v) % 2, charge_parity(v)


_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def dim_str(v) -> str:
    """A site as text, 'M L² T⁻¹' for (1, 2, −1, 0); '1' for the origin."""
    parts = [s + ("" if e == 1 else str(e).translate(_SUP)) for s, e in zip(AXES, _site(v)) if e]
    return " ".join(parts) or "1"


def planck_decompose(v) -> dict:
    """v = Δ·mass + ε·charge + n_c c + n_ħ ħ + n_ε ε₀, with integers: {'Delta', 'eps', 'n_c', 'n_hbar', 'n_eps0'}."""
    v = _site(v)
    D, eps = mass_dimension(v), charge_parity(v)
    w = [v[0] - D, v[1], v[2] - eps, v[3] - eps]
    ne, rem = divmod(w[3], 2)
    if rem:
        raise ArithmeticError("the charge exponent left over is odd")
    w2 = [x - ne * u for x, u in zip(w, UNITS["eps0"])]
    nh = w2[0]
    nc = w2[1] - 2 * nh
    back = tuple(D * m + eps * q + nc * c + nh * h + ne * e for m, q, c, h, e in zip(MASS, CHARGE, UNITS["c"], UNITS["hbar"], UNITS["eps0"]))
    if back != v:
        raise ArithmeticError("the decomposition does not rebuild the site")
    return {"Delta": D, "eps": eps, "n_c": nc, "n_hbar": nh, "n_eps0": ne}


@lru_cache(maxsize=1)
def lattice() -> dict:
    """The two unit sublattices of ℤ⁴ by exact integer linear algebra: Smith invariants, ranks, the quotients and the index."""
    from mtft.integral_lattice import smith_invariants
    nat = [list(UNITS[k]) for k in ("c", "hbar", "eps0")]
    pl = nat + [list(UNITS["G"])]
    s_nat = [int(x) for x in smith_invariants(nat)]
    s_pl = [int(x) for x in smith_invariants(pl)]
    index = 1
    for x in s_pl:
        index *= x
    return {"natural": {"smith": s_nat, "rank": len(s_nat), "quotient": "Z + Z/2" if s_nat == [1, 1, 2] else None},
            "planck": {"smith": s_pl, "rank": len(s_pl), "index": index, "quotient": "Z/2 + Z/2" if s_pl == [1, 1, 2, 2] else None}}


def everything() -> tuple:
    """The whole catalogue: Sphereland 25's list, the extension and the lab's own."""
    return CATALOGUE + EXTENSION + LAB


def names_at(v, lab: bool = True) -> list:
    """The catalogued quantities at a site."""
    v = _site(v)
    pool = CATALOGUE + EXTENSION + (LAB if lab else ())
    return [q for q in pool if q[1] == v]


def site(v, limit: int = 12) -> dict:
    """The prospector's answer for a site: what lives there, what pairs build it, its class and its Planck decomposition."""
    v = _site(v)
    here = names_at(v, lab=False)
    products = []
    for i, a in enumerate(CATALOGUE):
        for j in range(i, len(CATALOGUE)):
            b = CATALOGUE[j]
            if tuple(x + y for x, y in zip(a[1], b[1])) == v:
                products.append(f"{a[3]} · {b[3]}")
            if tuple(x - y for x, y in zip(a[1], b[1])) == v:
                products.append(f"{a[3]} / {b[3]}")
            if i != j and tuple(y - x for x, y in zip(a[1], b[1])) == v:
                products.append(f"{b[3]} / {a[3]}")
    seen, uniq = set(), []
    for p in products:
        if p not in seen:
            seen.add(p)
            uniq.append(p)
    cls = natural(v)
    mates = [q for q in CATALOGUE + EXTENSION if natural(q[1]) == cls and q[1] != v]
    return {"site": v, "dimension": dim_str(v), "here": here, "lab": [q for q in LAB if q[1] == v], "occupied": bool(here),
            "products": uniq[:limit], "product_count": len(uniq), "mates": mates, "Delta": cls[0], "eps": cls[1],
            "planck": planck_decompose(v)}


def rescale(v, factors: dict) -> float:
    """The factor a quantity at site v takes when the units of mass, length, time and current are multiplied by the given
    factors ({'M': …, 'L': …, 'T': …, 'I': …}; missing ones are 1).  Two quantities at one site take the same factor."""
    out = 1.0
    for s, e in zip(AXES, _site(v)):
        out *= float(factors.get(s, 1.0)) ** (-e)
    return out
