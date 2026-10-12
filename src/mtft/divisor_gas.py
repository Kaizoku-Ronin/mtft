"""divisor_gas.py — the weights as the energy of the primon gas (PG-01, 2026-10-12).

In the primon gas (Julia 1990; Bost–Connes 1995) the primes are the particles, a number n = prod p^k is the multi-particle state with occupation
numbers k, its energy is log n, and energy adds when states multiply: log(mn) = log m + log n — the "loggy laws" are the gas's additivity.
Unique factorisation gives every energy level multiplicity one, so the partition function of the full gas is Z(beta) = sum n^-beta = zeta(beta),
with the Hagedorn wall at beta = 1 (the pole).  The programme's weight

    w_n = sum_{d|n} (log d)/d = sum over the sub-states d of n of (energy of d) x (Boltzmann factor d^-beta) at beta = 1

is the energy of n's divisor gas at the Hagedorn temperature: w_n(beta) = sum_{d|n} log d . d^-beta = -d/dbeta sigma_{-beta}(n), and w_n = w_n(1).
Its Dirichlet series is sum w_n(beta) n^-s = -zeta(s) zeta'(s + beta), which at beta = 1 is the alphabet -zeta(s) zeta'(s+1) = Z_D(s) of
`marked_gas` (the Dirichlet ensemble, a Gibbs state with Hamiltonian log n and a KMS condition).

Two energies sit side by side (EXACT):
  * the Dirichlet ensemble's internal energy U_D(beta) = -d/dbeta log Z_D = <log n>_D, with heat capacity C = beta^2 Var(log n);
  * the unmarked gas's U(beta) = -zeta'(beta)/zeta(beta) = sum Lambda(n) n^-beta, the von Mangoldt series, whose poles are the zeros of zeta and the
    pole beta = 1.
Exactly, U_D(beta) = U(beta) - zeta''(beta+1)/zeta'(beta+1): the marked energy keeps every pole of U at the zeros of zeta (residue -1) and adds poles
at beta = s' - 1 for every zero s' of zeta' — marking does not move the energy off the zeros, it adds singularities (CC-65 corrected the first
statement of this wave, which said the opposite).  What is RH-blind is the partition function and the summatory function: Z_D(s) = -zeta(s) zeta'(s+1)
is entire except at s = 1 and vanishes at the zeros of zeta, so sum_{n<=x} w_n = -zeta'(2) x - (1/4) log^2 x - (1/2)(ln 2 pi) log x + O(1) carries
no term from the zeros.

Also here (EXACT): the mean weight -zeta'(2) = (pi^2/6)(12 ln A - gamma - ln 2pi) with A the Glaisher–Kinkelin constant (the 12 multiplying ln A is
-1/zeta(-1) = 2 pi^2/zeta(2), the functional equation at s = 2; the prefactor is zeta(2)); the alphabet at the negative integers,
W(-k) = zeta(-k) . (-zeta'(1-k)) = zeta(-k) (ln A_{k-1} - B_k H_{k-1}/k) with A_k the Bendersky–Glaisher constants (A_0 = sqrt(2 pi), A_1 = A) — the
Bendersky identity is the definition of A_{k-1} rearranged, and the closed forms at odd k are the functional equation zeta'(-2m) = (-1)^m (2m)!
zeta(2m+1)/(2 (2 pi)^{2m}) — which reads PV-02's W(-1) = -(ln 2pi)/24 as (the regularised sum of all integers) x (the regularised log of their product);
and the Euler characteristic: chi(SL(2,Z)) = 1/4 + 1/6 - 1/2 = -1/12 is the rational (Wall–Serre) Euler characteristic of the amalgam Z/4 *_{Z/2} Z/6,
equal to zeta(-1) by the K = Q case of Harder's Gauss–Bonnet theorem chi(SL_2(O_K)) = zeta_K(-1); the genus formula
g - 1 = -mu chi(SL(2,Z)) - nu_2/4 - nu_3/3 - c/2 is orbifold Gauss–Bonnet (Area(X(1))/4 pi = 1/12), 168/12 - 2 = 12 at level 143; the q^{1/24} of
eta is the zeta-regularised (1/2) sum n = zeta(-1)/2.  That |SL(2,Z)^ab| = 12 = 1/|chi(SL(2,Z))| is a coincidence of size, not a map (AG-D5).

The Bost–Connes facts quoted in `bost_connes()` are standard and not recomputed: for 0 < beta <= 1 a unique KMS state (the symmetry Gal(Q^ab/Q) = Zhat^*
unbroken), for beta > 1 the extremal KMS states are parametrised by that group (symmetry broken), and at beta -> infinity the ground states
realise class field theory.  No map from beta to a physical temperature is claimed (AG-D5).
Tests: tests/test_pg01_energy.py.  Report: docs/SM/PG01_ENERGY_TWO_CLOCKS.md.

Roger Tano — MTFT Research Program — October 2026
"""
from __future__ import annotations

import math
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp

GLAISHER_A = mp.mpf("1.28242712910062263687534256886979172776768892732500119206374002174040630885882646112973649195820237439")


# ------------------------------------------------------------------------------------------------------------------------------ the divisor gas of n
def divisor_energy(n, beta=1.0):
    """w_n(beta) = sum_{d|n} log d . d^-beta: the energy numerator of n's divisor gas; w_n = divisor_energy(n, 1) (EXACT)."""
    return sum(math.log(d) * d ** (-beta) for d in sp.divisors(int(n)))


def divisor_partition(n, beta=1.0):
    """sigma_{-beta}(n) = sum_{d|n} d^-beta: the partition function of n's divisor gas."""
    return sum(d ** (-beta) for d in sp.divisors(int(n)))


def divisor_gas(n, beta=1.0):
    """Partition function, energy numerator, mean energy and energy variance of the divisor gas of n at inverse temperature beta."""
    ds = sp.divisors(int(n)); Z = sum(d ** (-beta) for d in ds); E1 = sum(math.log(d) * d ** (-beta) for d in ds)
    E2 = sum(math.log(d) ** 2 * d ** (-beta) for d in ds)
    return {"n": int(n), "beta": beta, "Z": Z, "energy_numerator": E1, "mean_energy": E1 / Z, "variance": E2 / Z - (E1 / Z) ** 2,
            "is_weight": abs(beta - 1.0) < 1e-15, "divisors": len(ds)}


def weight_is_energy(N=2000):
    """w_n = -d/dbeta sigma_{-beta}(n) at beta = 1 against the package's sieve, for n <= N (EXACT; checked numerically)."""
    with mp.workdps(30):
        return _weight_is_energy(N)


def _weight_is_energy(N):
    from mtft.marked_gas import weights_sieve
    w = weights_sieve(N)
    worst = 0.0
    for n in range(2, N + 1):
        f = lambda b: sum(mp.mpf(d) ** (-b) for d in sp.divisors(n))
        worst = max(worst, abs(float(-mp.diff(f, 1)) - w[n]))
    return {"N": N, "max_abs_difference": worst, "agree": worst < 1e-9}


# ------------------------------------------------------------------------------------------------------------------------------ the two gases
def dirichlet_ensemble(beta):
    """Z_D(beta) = -zeta(beta) zeta'(beta+1), U_D = <log n> = -d log Z_D/dbeta, C = beta^2 Var(log n) (EXACT via mpmath)."""
    with mp.workdps(30):
        return _dirichlet_ensemble(beta)


def _dirichlet_ensemble(beta):
    b = mp.mpf(beta)
    Z = lambda x: -mp.zeta(x) * mp.zeta(x + 1, derivative=1)
    U = -mp.diff(lambda x: mp.log(Z(x)), b)
    var = -mp.diff(lambda x: -mp.diff(lambda y: mp.log(Z(y)), x), b)
    return {"beta": float(b), "Z": float(Z(b)), "mean_log_n": float(U), "variance_log_n": float(var), "heat_capacity": float(b * b * var),
            "rh_sensitive": True, "decomposition": "U_D = -zeta'/zeta(beta) - zeta''(beta+1)/zeta'(beta+1)"}


def unmarked_gas(beta):
    """The full primon gas: Z = zeta(beta), U = -zeta'/zeta = sum Lambda(n) n^-beta (von Mangoldt; poles at the zeros), C = beta^2 (log zeta)''."""
    with mp.workdps(30):
        return _unmarked_gas(beta)


def _unmarked_gas(beta):
    b = mp.mpf(beta)
    U = -mp.zeta(b, derivative=1) / mp.zeta(b)
    var = mp.diff(lambda x: mp.log(mp.zeta(x)), b, 2)
    return {"beta": float(b), "Z": float(mp.zeta(b)), "mean_energy": float(U), "variance": float(var), "heat_capacity": float(b * b * var),
            "von_mangoldt_series": "sum Lambda(n) n^-beta", "rh_sensitive": True}


def marked_energy_decomposition(beta=2.0):
    """U_D(beta) = U(beta) - zeta''(beta+1)/zeta'(beta+1) exactly: the marked energy contains the von Mangoldt series as a summand (EXACT)."""
    with mp.workdps(30):
        b = mp.mpf(beta)
        UD = -mp.diff(lambda x: mp.log(-mp.zeta(x) * mp.zeta(x + 1, derivative=1)), b)
        U = -mp.zeta(b, derivative=1) / mp.zeta(b); extra = -mp.zeta(b + 1, derivative=2) / mp.zeta(b + 1, derivative=1)
        return {"beta": float(b), "U_D": float(UD), "U": float(U), "extra": float(extra), "agree": abs(UD - U - extra) < mp.mpf(10) ** -20}


def residue_at_first_zero():
    """(beta - rho_1) U_D(beta) -> -1 at the first zero rho_1 of zeta: the marked energy keeps the pole of the unmarked one (EXACT; checked numerically)."""
    with mp.workdps(30):
        rho = mp.zetazero(1)
        f = lambda b: (b - rho) * (-mp.diff(lambda x: mp.log(-mp.zeta(x) * mp.zeta(x + 1, derivative=1)), b))
        val = f(rho + mp.mpf("1e-8"))
        return {"rho_1": complex(rho), "residue": complex(val), "is_minus_one": abs(val + 1) < 1e-6}


def von_mangoldt_check(beta=2.0, pmax=20000):
    """sum_{p^k} log p . p^{-k beta} against -zeta'(beta)/zeta(beta) (the tail beyond pmax is the only difference)."""
    with mp.workdps(30):
        return _von_mangoldt_check(beta, pmax)


def _von_mangoldt_check(beta, pmax):
    s = mp.mpf(0)
    for p in sp.primerange(2, pmax):
        s += mp.log(p) * sum(mp.mpf(p) ** (-k * beta) for k in range(1, 60))
    exact = -mp.zeta(beta, derivative=1) / mp.zeta(beta)
    return {"partial": float(s), "exact": float(exact), "difference": float(exact - s)}


def hagedorn_wall():
    """Both partition functions have their first singularity at beta = 1: zeta's pole; Z_D's pole from the zeta(beta) factor (zeta'(2) is finite)."""
    with mp.workdps(30):
        res = mp.limit(lambda b: (b - 1) * (-mp.zeta(b) * mp.zeta(b + 1, derivative=1)), 1)
        return {"beta_wall": 1.0, "residue_Z_D": float(res), "residue_equals_minus_zeta_prime_2": abs(float(res + mp.zeta(2, derivative=1))) < 1e-12,
                "reading": "the weights are evaluated at the wall: w_n = w_n(beta = 1)"}


# ------------------------------------------------------------------------------------------------------------------------------ Glaisher identities
def mean_weight_glaisher():
    """-zeta'(2) = (pi^2/6)(12 ln A - gamma - ln 2 pi): the density of the weights (sum_{n<=x} w_n ~ -zeta'(2) x) in Glaisher's constant (EXACT)."""
    with mp.workdps(40):
        lhs = -mp.zeta(2, derivative=1); rhs = (mp.pi ** 2 / 6) * (12 * mp.log(GLAISHER_A) - mp.euler - mp.log(2 * mp.pi))
        return {"minus_zeta_prime_2": float(lhs), "glaisher_form": float(rhs), "agree": abs(lhs - rhs) < mp.mpf(10) ** -30}


def partial_sum_density(x=100000):
    """sum_{n<=x} w_n / x against -zeta'(2): the error is -log^2 x/(4x) + O(log x/x), from the triple pole of Z_D(s) x^s/s at s = 0
    (sum_{n<=x} w_n = -zeta'(2) x - (1/4) log^2 x - (1/2)(ln 2 pi) log x + O(1))."""
    S = 0.0
    for d in range(1, x + 1):
        S += math.log(d) / d * (x // d)
    lim = float(-mp.zeta(2, derivative=1)); lx = math.log(x)
    return {"x": x, "density": S / x, "limit": lim, "error": S / x - lim, "predicted_error": -(lx * lx / 4 + math.log(2 * math.pi) * lx / 2) / x}


def bendersky_constant_log(k):
    """ln A_k = B_{k+1} H_k/(k+1) - zeta'(-k) (A_0 = sqrt(2 pi), A_1 = Glaisher's A; B_1 = 0 convention, H_0 = 0)."""
    k = int(k)
    B = mp.mpf(0) if k + 1 == 1 else mp.bernoulli(k + 1)
    H = mp.harmonic(k) if k > 0 else mp.mpf(0)
    return B * H / (k + 1) - mp.zeta(-k, derivative=1)


def alphabet_negative_values(kmax=7):
    """W(-k) = -zeta(-k) zeta'(1-k) = zeta(-k) (ln A_{k-1} - B_k H_{k-1}/k): the alphabet at the negative integers as Bernoulli x log Bendersky.
    W(-1) = -(ln 2 pi)/24 (PV-02's identity), W(-2m) = 0, W(-3) = zeta(3)/(480 pi^2), W(-5) = zeta(5)/(336 pi^4)."""
    with mp.workdps(40):
        return _alphabet_negative_values(kmax)


def _alphabet_negative_values(kmax):
    rows = []
    for k in range(1, kmax + 1):
        s = -k
        W = -mp.zeta(s) * mp.zeta(s + 1, derivative=1)
        B = mp.mpf(0) if k == 1 else mp.bernoulli(k); H = mp.harmonic(k - 1) if k > 1 else mp.mpf(0)
        pred = mp.zeta(s) * (bendersky_constant_log(k - 1) - B * H / k)
        rows.append({"k": k, "W": float(W), "bendersky_form": float(pred), "agree": abs(W - pred) < mp.mpf(10) ** -20, "ln_A_km1": float(bendersky_constant_log(k - 1))})
    checks = {"W(-1) = -(ln 2pi)/24": abs(rows[0]["W"] + float(mp.log(2 * mp.pi) / 24)) < 1e-15,
              "W(-3) = zeta(3)/(480 pi^2)": abs(rows[2]["W"] - float(mp.zeta(3) / (480 * mp.pi ** 2))) < 1e-15,
              "W(-5) = zeta(5)/(336 pi^4)": abs(rows[4]["W"] - float(mp.zeta(5) / (336 * mp.pi ** 4))) < 1e-15,
              "W(-2) = W(-4) = 0": rows[1]["W"] == 0 and rows[3]["W"] == 0,
              "A at the derivative: W'(-2) = (zeta(3)/4pi^2)(1/12 - ln A)": abs(float(mp.diff(lambda x: -mp.zeta(x) * mp.zeta(x + 1, derivative=1), -2))
                                                                           - float(mp.zeta(3) / (4 * mp.pi ** 2) * (mp.mpf(1) / 12 - mp.log(GLAISHER_A)))) < 1e-15}
    return {"rows": rows, "checks": checks, "reading": "W(-k) = (regularised sum of n^k) x (regularised sum of n^(k-1) log n): the convolution 1 * (log d)/d at s = -k"}


# ------------------------------------------------------------------------------------------------------------------------------ the clock is zeta(-1)
def euler_characteristic_is_zeta_minus_one(N=143):
    """chi(SL(2,Z)) = 1/4 + 1/6 - 1/2 = -1/12, the rational Euler characteristic of the amalgam Z/4 *_{Z/2} Z/6, equal to zeta(-1) (Harder, K = Q);
    eta's q^{1/24} is the zeta-regularised (1/2) sum n = zeta(-1)/2; the genus formula of X0(N) as orbifold Gauss–Bonnet in terms of chi(SL(2,Z)),
    evaluated at N (EXACT).  The match |SL(2,Z)^ab| = 12 = 1/|chi| is a coincidence of size (AG-D5)."""
    from mtft.matter.star_quotients import census
    chi = Fr(1, 4) + Fr(1, 6) - Fr(1, 2)
    row = next((r for r in census() if r["N"] == N), None)
    out = {"chi_SL2Z": chi, "zeta_minus_1": Fr(-1, 12), "equal": chi == Fr(-1, 12), "chi_PSL2Z": Fr(1, 2) + Fr(1, 3) - 1,
           "eta_multiplier_exponent": Fr(1, 24), "equals_minus_zeta_minus_1_over_2": Fr(1, 24) == -Fr(-1, 12) / 2,
           "abelianisation_match": "|SL(2,Z)^ab| = 12 = 1/|chi|: a coincidence of size, not a map (AG-D5)"}
    if row is not None:
        mu, c, g = row["index"], row["cusps"], row["genus"]
        nu2, nu3 = _elliptic_points(N)
        out.update({"N": N, "index": mu, "cusps": c, "genus": g, "g_minus_1_from_zeta": -mu * Fr(-1, 12) - Fr(nu2, 4) - Fr(nu3, 3) - Fr(c, 2),
                    "matches": -mu * Fr(-1, 12) - Fr(nu2, 4) - Fr(nu3, 3) - Fr(c, 2) == g - 1})
    return out


def _elliptic_points(N):
    """nu_2 = prod (1 + (-4/p)), nu_3 = prod (1 + (-3/p)) over the primes of squarefree N (zero if 4 | N or 9 | N)."""
    nu2 = 0 if N % 4 == 0 else math.prod(1 + int(sp.jacobi_symbol(-4, p)) if p != 2 else 1 for p in sp.primefactors(N))
    nu3 = 0 if N % 9 == 0 else math.prod(1 + int(sp.jacobi_symbol(-3, p)) if p != 3 else 1 for p in sp.primefactors(N))
    return nu2, nu3


def bost_connes():
    """Standard facts (Bost–Connes 1995; uniqueness for beta <= 1 also Neshveyev 2002), quoted not recomputed."""
    return {"partition_function": "zeta(beta)", "wall": 1.0,
            "hot_phase": "0 < beta <= 1: a unique KMS_beta state; the symmetry Gal(Q^ab/Q), isomorphic to Zhat^* by Kronecker–Weber, is unbroken",
            "cold_phase": "beta > 1: Zhat^* acts freely and transitively on the extremal KMS_beta states; the symmetry is spontaneously broken",
            "ground_states": "beta -> infinity: the extremal states evaluated on the arithmetic subalgebra take values in Q^ab, and the Zhat^* action on "
                             "states intertwines the Galois action on those values — the isomorphism Gal(Q^ab/Q) = Zhat^* is an input, not an output",
            "marked_gas": "mtft.marked_gas is its own Gibbs construction on l^2(N >= 2) whose construction note deforms the Bost–Connes symmetry twist in "
                          "the IR; it is not a KMS state of the Bost–Connes algebra and nothing about that algebra is claimed here",
            "weights_at_the_wall": "w_n = w_n(beta) at beta = 1, the parameter value of the Bost–Connes critical point; the divisor gas of one n is finite "
                                   "and has no transition, and nothing follows from the coincidence of parameter values",
            "physical_temperature": "no map from beta to a physical temperature is claimed (AG-D5)"}


def report():
    return {"divisor_gas_143": divisor_gas(143), "divisor_gas_144": divisor_gas(144), "weight_is_energy": weight_is_energy(500),
            "ensemble": {b: dirichlet_ensemble(b) for b in (1.05, 1.5, 2.0, 3.0, 5.0)}, "unmarked": {b: unmarked_gas(b) for b in (1.05, 1.5, 2.0, 3.0, 5.0)},
            "von_mangoldt": von_mangoldt_check(), "wall": hagedorn_wall(), "mean_weight": mean_weight_glaisher(), "density": partial_sum_density(20000),
            "negative_values": alphabet_negative_values(), "clock": euler_characteristic_is_zeta_minus_one(), "decomposition": marked_energy_decomposition(),
            "first_zero": residue_at_first_zero(), "bost_connes": bost_connes()}


def text():
    r = report()
    L = ["PG-01 — the weights as the energy of the primon gas (2026-10-12)", ""]
    for n in (143, 144):
        g = r[f"divisor_gas_{n}"]
        L.append(f"divisor gas of {n}: {g['divisors']} divisors, Z = {g['Z']:.4f}, energy numerator w_{n} = {g['energy_numerator']:.4f}, mean energy {g['mean_energy']:.4f}")
    L.append(f"w_n = -d/dbeta sigma_(-beta)(n) at beta = 1 for n <= {r['weight_is_energy']['N']}: max difference {r['weight_is_energy']['max_abs_difference']:.1e}")
    L.append("beta   Z_D       <log n>_D   Var_D     | zeta(beta)  -zeta'/zeta (unmarked)")
    for b in (1.05, 1.5, 2.0, 3.0, 5.0):
        e, u = r["ensemble"][b], r["unmarked"][b]
        L.append(f"{b:<6} {e['Z']:<9.4f} {e['mean_log_n']:<11.4f} {e['variance_log_n']:<9.4f} | {u['Z']:<11.4f} {u['mean_energy']:.4f}")
    L.append(f"Hagedorn wall at beta = 1: residue of Z_D = {r['wall']['residue_Z_D']:.6f} = -zeta'(2): {r['wall']['residue_equals_minus_zeta_prime_2']}")
    L.append(f"mean weight -zeta'(2) = {r['mean_weight']['minus_zeta_prime_2']:.9f} = (pi^2/6)(12 ln A - gamma - ln 2pi): {r['mean_weight']['agree']}; density at x = {r['density']['x']}: {r['density']['density']:.5f}")
    L.append("alphabet at negative integers: " + "; ".join(f"W({-row['k']}) = {row['W']:.3e}" for row in r["negative_values"]["rows"]) + f"; checks {all(r['negative_values']['checks'].values())}")
    c = r["clock"]
    L.append(f"Euler characteristic: chi(SL2Z) = {c['chi_SL2Z']} = zeta(-1): {c['equal']}; genus formula at {c.get('N')}: g - 1 = {c.get('g_minus_1_from_zeta')} (matches {c.get('matches')})")
    L.append(f"marked energy: U_D = U - zeta''(b+1)/zeta'(b+1) at beta = 2: {r['decomposition']['agree']}; residue of U_D at the first zero: {r['first_zero']['residue'].real:.6f} (is -1: {r['first_zero']['is_minus_one']})")
    return "\n".join(L)


if __name__ == "__main__":
    print(text())
