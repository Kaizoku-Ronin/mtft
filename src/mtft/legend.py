#!/usr/bin/env python3
"""
The Legend — a map key to the arithmetic territory
====================================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

mtft computes artifacts of the number line — structures older than the
universe whose physical shadows are the constants of nature. A package
like that needs more than help(): it needs a LEGEND in the cartographic
sense. Every entry here answers four questions ordinary documentation
cannot:

  WHAT   — the nature of the object (identity / ensemble / diagnostic /
           engine / certificate / standard / constant / definition)
  WHERE  — its tier, and its Arithmetica Generale primitive signature:
           which of the five cognitive acts it exercises —
              I ITERATE (*)   II DIVIDE (/)   III ASSEMBLE (S)
              IV EXTRACT (^)  V CURVE (d)
  HOW TRUE — the corpus epistemic tag (Df definition, Pp proposition,
           Pr proved, Conj conjecture, Heur heuristic, Cert machine
           certificate) and an exactness class:
              EXACT        closed-form identity, machine-precision
              CERTIFIED(e) numeric, sealed to stated tolerance
              DIAGNOSTIC   monotone bounds / truncations / scans
              PHENO        confronts measured data
  WHENCE — upstream derivation links. Because the pipeline has zero
           free parameters, every chain terminates in the integers.
           `trace` walks it. No other scientific package can ship this
           honestly: everyone else's chains end in a fitted number.

Interfaces:
    python -m mtft.legend            the map (tiers, glyphs, tags)
    python -m mtft.legend card NAME  one entry, full detail
    python -m mtft.legend trace NAME derivation chain down to N
    python -m mtft.legend status     epistemic audit of the surface
    python -m mtft.legend status Pr  filter by tag (or EXACT, etc.)
    python -m mtft.legend search T   substring search

    >>> import mtft; mtft.what("dirichlet_curvature")

Roger Tano — MTFT Research Program — July 2026
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ── glyphs & tags ─────────────────────────────────────────────

PRIMITIVE_GLYPHS = {"I": "*", "II": "/", "III": "S", "IV": "^", "V": "d"}
PRIMITIVE_NAMES = {"I": "ITERATE", "II": "DIVIDE", "III": "ASSEMBLE",
                   "IV": "EXTRACT", "V": "CURVE"}

TAGS = ("Df", "Pp", "Pr", "Conj", "Heur", "Cert")
EXACTNESS = ("EXACT", "CERTIFIED", "DIAGNOSTIC", "PHENO", "GIVEN")

_ANSI = {"Df": "\033[36m", "Pp": "\033[33m", "Pr": "\033[32m",
         "Conj": "\033[35m", "Heur": "\033[90m", "Cert": "\033[92m",
         "reset": "\033[0m", "dim": "\033[2m", "bold": "\033[1m"}


def _color_on() -> bool:
    return sys.stdout.isatty() and not os.environ.get("NO_COLOR")


def _c(tag: str, s: str) -> str:
    if not _color_on():
        return s
    return f"{_ANSI.get(tag, '')}{s}{_ANSI['reset']}"


# ── the entry ─────────────────────────────────────────────────

@dataclass
class LegendEntry:
    name: str
    tier: str                    # "0".."11" or "GIVENS"
    kind: str                    # constant/definition/identity/ensemble/...
    primitives: Tuple[str, ...]  # subset of ("I","II","III","IV","V")
    tag: str                     # Df/Pp/Pr/Conj/Heur/Cert
    exactness: str               # EXACT / CERTIFIED(...) / DIAGNOSTIC / PHENO / GIVEN
    nature: str                  # one line: what it IS
    example: str = ""            # how to call it
    upstream: Tuple[str, ...] = ()   # names of parents in the derivation
    ref: str = ""                # paper / test anchor

    @property
    def glyphs(self) -> str:
        return "".join(PRIMITIVE_GLYPHS[p] for p in self.primitives)


REGISTRY: Dict[str, LegendEntry] = {}


def _reg(e: LegendEntry):
    REGISTRY[e.name] = e


# ═══════════════════════════════════════════════════════════════
#  THE GIVENS — older than the universe
# ═══════════════════════════════════════════════════════════════

_reg(LegendEntry("integers", "GIVENS", "constant", (), "Df", "GIVEN",
                 "N. The number line. Every chain below ends here.",
                 ref="Peano 1889; AG Ch.0"))
_reg(LegendEntry("N_143", "GIVENS", "constant", ("II",), "Df", "GIVEN",
                 "143 = 11 x 13 — the level. Semiprime; two odd primes.",
                 upstream=("integers",)))
_reg(LegendEntry("X0_143", "GIVENS", "constant", ("V",), "Df", "GIVEN",
                 "The modular curve X0(143): genus 13, four cusps.",
                 example="from mtft.modular_curve import X0",
                 upstream=("N_143",)))
_reg(LegendEntry("genus_13", "GIVENS", "constant", ("V",), "Pr", "EXACT",
                 "genus(X0(143)) = 13. The geometric backbone number.",
                 upstream=("X0_143",), ref="LMFDB 143.2"))
_reg(LegendEntry("index_168", "GIVENS", "constant", ("II",), "Pr", "EXACT",
                 "[SL2(Z):Gamma0(143)] = 168 = |PSL(2,7)| = genus^2 - 1.",
                 upstream=("X0_143",)))
_reg(LegendEntry("dim_S2_11", "GIVENS", "constant", ("III",), "Pr", "EXACT",
                 "dim S2^new(Gamma0(143)) = 11; Galois orbits [1,4,6].",
                 upstream=("X0_143",), ref="LMFDB"))
_reg(LegendEntry("monster_order", "GIVENS", "constant", ("III",), "Pr", "EXACT",
                 "|M| ~ 8.08e53 — Monster group order; dim V-natural = 196883.",
                 upstream=("integers",), ref="Griess 1982; moonshine"))
_reg(LegendEntry("w_n", "GIVENS", "definition", ("I", "II", "III"), "Df", "EXACT",
                 "w_n = sum_{d|n} (log d)/d — the holonomy weights. THE input.",
                 example="mtft.arithmetic.weights(n_max)",
                 upstream=("integers",)))
_reg(LegendEntry("tano_weight_lattice", "GIVENS", "definition",
                 ("I", "II", "III"), "Df", "EXACT",
                 "Two-species lattice w^(a,b): skeleton (Lambda, zero-tuned) "
                 "vs bulk (log, zero-blind); unifies stiffness, bulk, and Z_D "
                 "as lattice points.",
                 upstream=("w_n",), ref="L1/L2 draft Df 0 (two-engine, Jul 2026)"))

# ═══════════════════════════════════════════════════════════════
#  MARQUEE ENTRIES BY TIER (curated; introspection covers the rest)
# ═══════════════════════════════════════════════════════════════

_reg(LegendEntry("alpha_inverse", "3", "identity", ("I", "II", "III", "IV", "V"),
                 "Pr", "CERTIFIED(3.5ppm)",
                 "alpha^-1 = ln|M| + genus - 1/11 + O(1e-4) = 137.0355...",
                 example="python -m mtft verify",
                 upstream=("monster_order", "genus_13", "dim_S2_11"),
                 ref="Paper 27; corrected ln|M| (audit B2)"))
_reg(LegendEntry("X0_143_engine", "1", "engine", ("III", "IV"), "Pr", "EXACT",
                 "ModularCurve engine for X0(143): genus, cusps, Hecke "
                 "spectrum, homology; newform orbits [1,4,6] LMFDB-anchored "
                 "in x0_143.py.",
                 example="mtft.X0(143).summary()",
                 upstream=("X0_143",)))
_reg(LegendEntry("HosotaniMTFT", "2", "engine", ("I", "III", "V"), "Pp", "PHENO",
                 "Gauge-Higgs unification: Higgs as A_tau holonomy; "
                 "center-projected effective potential.",
                 example="mtft.HosotaniMTFT().find_vacuum()",
                 upstream=("w_n",)))
_reg(LegendEntry("mu_stiffness", "5c", "identity", ("I", "II", "III", "V"),
                 "Pr", "EXACT",
                 "mu_N(y) = sum n^2 w_n e^{-2pi y n}(1 - cos 2pi n/N): "
                 "the Laplace-ensemble mass-gap object.",
                 example="mtft.tower / viz.make_hero",
                 upstream=("w_n",)))
_reg(LegendEntry("filtered_moment_identity", "5d", "identity",
                 ("I", "III", "IV", "V"), "Pr", "EXACT",
                 "mu_N(y) = (1/4pi^2)[T''(y) - Re T''(y - i/N)] — mass gap "
                 "as twisted-minus-untwisted curvature. Machine precision.",
                 example="mtft.filtered_moment_identity(0.18174, N=3)",
                 upstream=("mu_stiffness",), ref="v0.7.1; N=3 = SU(3) center"))
_reg(LegendEntry("y_c", "5c", "constant", ("V",), "Pr", "CERTIFIED",
                 "Confinement depth y_conf = 0.18174 (canonical; "
                 "CriticalDepths). Bulk zero of the stiffness landscape.",
                 upstream=("mu_stiffness",)))
_reg(LegendEntry("corrected_rh_diagnostic", "5d", "diagnostic",
                 ("I", "III", "IV", "V"), "Pr", "DIAGNOSTIC",
                 "Th 1: RH <=> limsup |(kappa^L - kappa_Main)(2pi y)^{-3/2}| "
                 "< inf. Supersedes the false kappa>=0 propositions.",
                 example="mtft.corrected_rh_diagnostic(...)",
                 upstream=("tano_weight_lattice",),
                 ref="L1/L2 draft Th 1; conspiracy/crossover scans"))
_reg(LegendEntry("Z_D_closed_form", "9", "identity", ("III", "IV"), "Pr", "EXACT",
                 "Z_D(beta) = sum w_n n^-beta = -zeta(beta) zeta'(beta+1). "
                 "The Dirichlet ensemble has a closed form.",
                 example="mtft.dirichlet_ensemble(beta)",
                 upstream=("w_n",), ref="arithmetic_wick; audit-verified"))
_reg(LegendEntry("dirichlet_curvature", "5d", "identity", ("IV", "V"),
                 "Pr", "EXACT",
                 "g_D(beta) = d2 log zeta(beta) + d2 log(-zeta'(beta+1)) — "
                 "Fisher-Rao curvature splits: prime piece + zeta' piece.",
                 example="mtft.dirichlet_curvature(3.0)",
                 upstream=("Z_D_closed_form",),
                 ref="Speiser (1935) lives in the zeta' piece"))
_reg(LegendEntry("marked_gas", "5d", "ensemble", ("IV", "V"), "Pr", "EXACT",
                 "Marked primon gas: Gibbs state rho_n = (log n) n^-(beta+1) / "
                 "(-zeta'(beta+1)) on l2(N>=2); H = log Q; spectrum "
                 "E_n = (beta+1) log n - log log n; prime-shift isometries.",
                 example="mtft.marked_gas.gates()",
                 upstream=("Z_D_closed_form",),
                 ref="Marked Primon Gas note v0.1.1 (July 2026)"))
_reg(LegendEntry("kms_flow", "5d", "identity", ("IV", "V"), "Pr",
                 "CERTIFIED(1e-12)",
                 "Modular flow alpha_t = Ad e^{itK} satisfies KMS at t+i "
                 "termwise: rho_n e^{-Delta E_n} = rho_{pn}; the wrong-sign "
                 "control (t-i) fails by 4.574, as it must.",
                 example="mtft.marked_gas.kms_check(p=2, beta=2.0, t=1.0)",
                 upstream=("marked_gas",),
                 ref="note v0.1.1 sec.5; BC twist vanishes in the UV"))
_reg(LegendEntry("cold_gas_amplitude", "5d", "identity", ("III", "IV"), "Pr",
                 "EXACT",
                 "Cold gas: A_n ~ B n^alpha / Gamma(alpha+1) with "
                 "alpha = -zeta'(2), B = e^{-zeta''(2)} (gamma-cancellation); "
                 "certified 0.14040 vs closed form 0.14027492... (0.09%) at "
                 "N=1e5, converging from above.",
                 example="mtft.marked_gas.cold_gas_report(100_000)",
                 upstream=("Z_D_closed_form",),
                 ref="note v0.1.1 sec.6; Karamata Tauberian"))
_reg(LegendEntry("spectral_edge_soft", "5d", "diagnostic", ("IV", "V"), "Pr",
                 "DIAGNOSTIC",
                 "Edge softness at (beta+1) log p: mass(gap<eps) = "
                 "[M^-beta (log M/beta + 1/beta^2) + O(log M M^-beta-1)] / "
                 "(-zeta'(beta+1)), M = exp(log p/(e^eps-1)); per-level "
                 "convention pinned.",
                 example="mtft.marked_gas.edge_mass(p=2, beta=2.0, eps=0.1)",
                 upstream=("marked_gas",),
                 ref="note v0.1.1; audit Addendum U.4 (the auditor's law)"))
_reg(LegendEntry("hadamard_zetaprime_check", "5d", "identity", ("III", "IV", "V"),
                 "Pr", "CERTIFIED(1e-5)",
                 "d2 log(-zeta'(s)) = 2/(s-1)^2 - sum_{rho'}(s-rho')^-2 over "
                 "the certified 19-zero census; use s in [3,10] "
                 "(conditioning).",
                 example="mtft.hadamard_zetaprime_check(3.0)",
                 upstream=("dirichlet_curvature",), ref="v0.7.1 Addendum I"))
_reg(LegendEntry("li_lambda", "10", "identity", ("I", "III", "IV", "V"),
                 "Pr", "EXACT",
                 "Li coefficients lambda_n: Taylor family of log xi at s=1; "
                 "three cross-certified methods; lambda_1 = 1+gamma/2-ln(4pi)/2.",
                 example="mtft.li_criterion_report(12)",
                 upstream=("integers",),
                 ref="Bombieri-Lagarias caveat attached to every report"))
_reg(LegendEntry("li_lambda_zero_sum", "10", "diagnostic", ("I", "III"),
                 "Pr", "DIAGNOSTIC",
                 "Truncated zero-sum: monotone lower bounds; tail is "
                 "density-dominated (S(T) ~ 0.1% at gamma~236).",
                 example="mtft.li_lambda_zero_sum(3, 100)",
                 upstream=("li_lambda",)))
_reg(LegendEntry("falsify_engine", "5b", "engine", ("I", "IV"), "Pp", "PHENO",
                 "23 pre-registered zero-parameter predictions; ppm ledger; "
                 "no cherry-picking.",
                 example="mtft.falsify.honest_report()",
                 upstream=("alpha_inverse", "mu_stiffness")))
_reg(LegendEntry("monster_hash", "7", "engine", ("I", "II", "III", "IV"),
                 "Heur", "DIAGNOSTIC",
                 "SL(2,Z)-sponge hash; 49.95% avalanche; research-grade "
                 "(collision-resistance gap documented).",
                 example="from mtft.monster_hash import MonsterHash",
                 upstream=("X0_143",)))
_reg(LegendEntry("decompose_turing_machine", "9", "definition",
                 ("I", "II", "III", "IV", "V"), "Df", "EXACT",
                 "A Turing machine as a five-primitive AG object; "
                 "computation generated, not imported.",
                 example="mtft.decompose_turing_machine()",
                 upstream=("integers",), ref="arithmetic_machine v0.7.0"))
_reg(LegendEntry("lhcb_bridge", "8", "engine", ("IV",), "Pp", "PHENO",
                 "ROOT-ntuple confrontation via uproot; Run 1 open data "
                 "(DOI-cited) in scope.",
                 example="from mtft.lhcb_analysis import ...",
                 upstream=("falsify_engine",)))

# ── Tier 11: Certificates & Standards (July 2026) ─────────────

_reg(LegendEntry("jc_counterexample", "11", "certificate",
                 ("I", "II", "III", "IV"), "Cert", "EXACT",
                 "The degree-7 Jacobian Conjecture counterexample "
                 "(Alpoge/Claude Fable 5, Jul 20 2026): det DF = -2, "
                 "degree 3, S3 monodromy, missed curve. Self-certifying; "
                 "re-derived dependency-free on every call. JC false n>=3.",
                 example="mtft.jc_verify_all(verbose=True)",
                 upstream=("integers",),
                 ref="scripts/jc/; 7*6*4=168 flagged AG-D5, dismissed"))
_reg(LegendEntry("binned_log_slope", "11", "standard", ("I", "IV"),
                 "Df", "EXACT",
                 "A.7 slope estimator: terminal-bin guard (min_bin) kills "
                 "the 0.035 leverage bias.",
                 example="mtft.binned_log_slope(ys, vals)",
                 upstream=("integers",), ref="estimator_standards; L1/L2 A.7"))
_reg(LegendEntry("stride_resonance_check", "11", "standard", ("I", "II"),
                 "Df", "EXACT",
                 "Flags sliding-window strides near-resonant with a known "
                 "oscillation (gamma_3 x 6-decade: 54.994 cycles, the trap).",
                 example="mtft.stride_resonance_check(25.0109, 6.0)",
                 upstream=("integers",)))

TIER_TITLES = {
    "GIVENS": "THE GIVENS — older than the universe",
    "0": "Constants & Arithmetic", "1": "Modular Geometry & Forms",
    "2": "Gauge-Higgs", "3": "Phenomenology", "4": "Lattice & Fermions",
    "5": "Bridges & Dark Sector", "5b": "Falsifiability",
    "5c": "Multi-N Tower", "5d": "Riemann / Speiser-Hadamard",
    "5e": "LHC Materials", "6": "Quantum", "7": "Cryptography",
    "8": "LHC Confrontation", "9": "Computation & Wick",
    "10": "Critical Ensemble (Li)", "11": "Certificates & Standards",
}


# ═══════════════════════════════════════════════════════════════
#  VIEWS
# ═══════════════════════════════════════════════════════════════

def _entry_line(e: LegendEntry) -> str:
    tag = _c(e.tag, f"[{e.tag:>4s}]")
    gl = f"{e.glyphs:<5s}"
    ex = e.exactness.split("(")[0]
    return f"  {tag} {gl} {e.name:<28s} {ex:<10s} {e.nature}"


def legend_map() -> str:
    """The map: all registered entries grouped by tier."""
    tiers: Dict[str, List[LegendEntry]] = {}
    for e in REGISTRY.values():
        tiers.setdefault(e.tier, []).append(e)
    order = ["GIVENS"] + sorted((t for t in tiers if t != "GIVENS"),
                                key=lambda s: (len(s), s))
    out = [_c("bold", "MTFT LEGEND — map key to the arithmetic territory"),
           _c("dim", "  glyphs: * ITERATE  / DIVIDE  S ASSEMBLE  "
                     "^ EXTRACT  d CURVE"),
           _c("dim", "  tags:   Df Pp Pr Conj Heur Cert   ·   "
                     "trace NAME follows any chain to N"), ""]
    for t in order:
        out.append(_c("bold", f"─── Tier {t}: {TIER_TITLES.get(t, '')} ───"))
        for e in sorted(tiers[t], key=lambda x: x.name):
            out.append(_entry_line(e))
        out.append("")
    return "\n".join(out)


def card(name: str) -> str:
    """Full detail for one entry (falls back to live introspection)."""
    e = REGISTRY.get(name)
    if e is None:
        return _introspect_card(name)
    lines = [
        _c("bold", f"{e.name}") + f"   {_c(e.tag, '[' + e.tag + ']')} "
        f"{e.exactness}   tier {e.tier} · {e.kind}",
        f"  primitives : {e.glyphs}  "
        + " ".join(PRIMITIVE_NAMES[p] for p in e.primitives),
        f"  nature     : {e.nature}",
    ]
    if e.example:
        lines.append(f"  usage      : {e.example}")
    if e.upstream:
        lines.append(f"  upstream   : {', '.join(e.upstream)}")
    if e.ref:
        lines.append(f"  ref        : {e.ref}")
    return "\n".join(lines)


def _introspect_card(name: str) -> str:
    """Unregistered names: live signature + docstring first line."""
    try:
        import mtft as _m
        obj = getattr(_m, name)
    except Exception:
        return f"  '{name}' not in registry and not an mtft export."
    import inspect
    try:
        sig = str(inspect.signature(obj))
    except (TypeError, ValueError):
        sig = ""
    doc = (inspect.getdoc(obj) or "").split("\n")[0]
    return (f"{name}{sig}   [UNREGISTERED — introspected]\n"
            f"  {doc}\n"
            f"  (curated legend entry pending; epistemic tag unknown)")


def trace(name: str, _depth: int = 0,
          _path: Tuple[str, ...] = ()) -> str:
    """Walk the derivation chain down to the integers. Branches may
    share ancestors (the chains CONVERGE on N — that is the theorem);
    only a true cycle along the current path is cut."""
    e = REGISTRY.get(name)
    pad = "   " * _depth + ("└─ " if _depth else "")
    if e is None:
        return f"{pad}{name}  [unregistered]"
    line = (f"{pad}{_c(e.tag, e.name)}  "
            f"{_c('dim', '[' + e.tag + ', ' + e.exactness + ']')}")
    out = [line]
    if name in _path:
        out.append("   " * (_depth + 1) + "└─ (cycle — cut)")
        return "\n".join(out)
    for up in e.upstream:
        out.append(trace(up, _depth + 1, _path + (name,)))
    if not e.upstream and name == "integers":
        out.append("   " * (_depth + 1)
                   + _c("dim", "└─ end of chain. It was always the integers."))
    return "\n".join(out)


def status(filter_tag: Optional[str] = None) -> str:
    """Epistemic audit of the registered surface."""
    entries = list(REGISTRY.values())
    if filter_tag:
        entries = [e for e in entries
                   if e.tag == filter_tag
                   or e.exactness.split("(")[0] == filter_tag]
        hdr = f"entries with tag/exactness = {filter_tag}"
    else:
        hdr = "epistemic audit of the registered surface"
    out = [_c("bold", f"LEGEND STATUS — {hdr}"), ""]
    if filter_tag:
        for e in sorted(entries, key=lambda x: (x.tier, x.name)):
            out.append(_entry_line(e))
    else:
        by_tag: Dict[str, int] = {}
        by_ex: Dict[str, int] = {}
        for e in entries:
            by_tag[e.tag] = by_tag.get(e.tag, 0) + 1
            by_ex[e.exactness.split("(")[0]] = \
                by_ex.get(e.exactness.split("(")[0], 0) + 1
        out.append("  by tag:      "
                   + "  ".join(f"{_c(t, t)}={by_tag.get(t, 0)}" for t in TAGS))
        out.append("  by exactness: "
                   + "  ".join(f"{x}={by_ex.get(x, 0)}" for x in EXACTNESS))
        out.append(f"  total registered: {len(entries)} "
                   f"(unregistered exports fall back to introspection)")
    return "\n".join(out)


def search(term: str) -> str:
    t = term.lower()
    hits = [e for e in REGISTRY.values()
            if t in e.name.lower() or t in e.nature.lower()]
    if not hits:
        return f"  no legend entries matching '{term}'"
    return "\n".join(_entry_line(e)
                     for e in sorted(hits, key=lambda x: x.name))


def what(name: str) -> None:
    """REPL-speed: print the card."""
    print(card(name))


def legend() -> None:
    """REPL-speed: print the map."""
    print(legend_map())



# ── Tier 0b ancestry registrations (combinatorial, v0.13.0) ──
from mtft.combinatorial import ANCESTRY_LEGEND as _ANCESTRY_LEGEND
for _d in _ANCESTRY_LEGEND:
    if _d["name"] not in REGISTRY:
        _reg(LegendEntry(**_d))

# ── Tier 12 promotion registrations (v0.14.0) ─────────────────
# Constants nominated by INTEGRATION v0.14.0 (auditor-composed entries;
# values certified in mtft.moments / mtft.curvature / mtft.eisenstein).
PROMOTION_LEGEND: Tuple[Dict[str, object], ...] = (
    dict(name="tano_second_moment_closed_form", tier="0", kind="constant",
         primitives=("I", "IV"), tag="Pr", exactness="EXACT",
         nature="<w^2>_beta = zeta''(v)C0 − 2 zeta'(v)C1 + zeta(v)C2, "
                "v = beta+2, C-blocks in zeta(beta+1), zeta(2beta+2). "
                "chi_w = T − zeta'(beta+1)^2; Cov(log n, w) = zeta''(beta+1). "
                "Cold values T(1) = 1.70276979154901697001, "
                "chi_w(1) = 0.82377306237833093427.",
         example="mtft.moments.weight_second_moment(1)",
         upstream=("w_n",), ref="studies/promotion_2026aug/w2_susceptibility.py"),
    dict(name="hessian_cancellation_theorem", tier="0", kind="identity",
         primitives=("IV",), tag="Pr", exactness="EXACT",
         nature="For a Hessian metric the Brioschi second-derivative block "
                "is −1/2 k4 + k4 − 1/2 k4 = 0: the fourth cumulant "
                "contributes NOTHING to the curvature of any exponential "
                "family. Gaussian convention lock K = −1/2.",
         example="mtft.curvature.gaussian_family_curvature()",
         upstream=("tano_second_moment_closed_form",),
         ref="studies/promotion_2026aug/curvature_tano_manifold.py (11 gates)"),
    dict(name="hagedorn_slope_A", tier="0", kind="constant",
         primitives=("IV",), tag="Pr", exactness="CERTIFIED(1e-28)",
         nature="dK/dbeta at the Hagedorn edge: A = (zeta''(2) kappa3_cold − "
                "kappa_wwl_cold chi_cold) / (2 chi_cold^2) = "
                "0.423657463797093480081718158187. The wall is approached "
                "FLAT, linearly from positive curvature.",
         example="mtft.curvature.hagedorn_slope()",
         upstream=("hessian_cancellation_theorem",),
         ref="studies/promotion_2026aug/curvature_tano_manifold.py G5"),
    dict(name="curvature_profile_milestones", tier="0", kind="constant",
         primitives=("IV",), tag="Pr", exactness="CERTIFIED(1e-9)",
         nature="The (beta, lambda=0) profile: positive dome, summit "
                "beta* = 4.593591164956, K* = 1.19569598199193852905; "
                "flat temperature beta_0 = 8.8565170425 (sign change).",
         example="mtft.curvature.gaussian_curvature(4.593591164956)",
         upstream=("hessian_cancellation_theorem",),
         ref="studies/promotion_2026aug/curvature_tano_manifold.py G6"),
    dict(name="cold_dive_law", tier="0", kind="constant",
         primitives=("IV",), tag="Pr", exactness="CERTIFIED(1e-15)",
         nature="Cold tail K ~ −c (6/5)^beta, c = 0.270126465305424759517602 "
                "[RETRACTED by CC-04: six-atom beta=200 extraction carried "
                "atom-6 contamination (5/6)^200 = 1.46e-16; correct closed "
                "form in cold_amplitude_closed_form_CC04]; "
                "rate 6/5 = Boltzmann ratio of the first new prime (5) to "
                "the first mixed composite (6). Rigidity locks: K = 1/4 "
                "identically on {1,2,3} and {1,2,3,4} (atom 4 inert); "
                "atom 5 flips the sign. Deep-cold geometry = geometry of "
                "the first six integers.",
         example="mtft.curvature.finite_atom_curvature(60, (1,2,3,4,5,6))",
         upstream=("hessian_cancellation_theorem",),
         ref="studies/promotion_2026aug/curvature_tano_manifold.py G7"),
    dict(name="congruence_primes_X0143", tier="1", kind="constant",
         primitives=("II", "III"), tag="Pr", exactness="CERTIFIED(Sturm for 5; census E2-sampled)",
         nature="Eisenstein congruence norm-moduli of the four Hecke "
                "blocks of X0(143): 143a1 → 1 (no congruence), 11a1 ghost "
                "→ 5 (= Mazur numerator((11−1)/12), Sturm-certified), "
                "f2 quartic → 7, f3 sextic → 12 = 2^2·3.",
         example="mtft.eisenstein.congruence_census()",
         upstream=("X0_143_engine",),
         ref="studies/promotion_2026aug/x0143_hecke_particles.py; "
             "eisenstein_congruences.py (pending drop)"),
)
for _d in PROMOTION_LEGEND:
    if _d["name"] not in REGISTRY:
        _reg(LegendEntry(**_d))

# ── Tier 13 v0.15.0 registrations (certification wave) ─────────
V0150_LEGEND: Tuple[Dict[str, object], ...] = (
    dict(name="wseries_shift_chain_CC02", tier="0", kind="identity",
         primitives=("I",), tag="Pr", exactness="EXACT",
         nature="sum w_n n^-s = F(s+1) = −zeta(s)·zeta'(s+1) (CC-02; "
                "Paper 1 Prop 1.5's −zeta'(s+1) is the summand series, "
                "AG Pr 4.1.4 mistranscribes). Three-route E2: sieve "
                "7.1e-15, s=3 numeric 5.2e-12, F-consistency.",
         example="mtft.weil.w_series_check(3)",
         upstream=("w_n",),
         ref="studies/v0150_2026aug/CC-02_wseries_shift_chain.md"),
    dict(name="curvature_rigidity_theorem", tier="0", kind="identity",
         primitives=("IV",), tag="Pr", exactness="EXACT",
         nature="A 2-dim discrete exponential family supported on one "
                "affine line plus exactly one off-line point has K = 1/4 "
                "identically (warped form ds^2 = dy^2 + cos^2(y/2) drho^2). "
                "{1,2,3,4} is in the class: X_4 = 2 X_2 exactly. Out-of-"
                "sample {1,2,4,8}, {1,2,4,16} to 1e-115; sharpness gated "
                "both directions. Classify with rigidity_class().",
         example="mtft.curvature.rigidity_class([(0,0),(1,0),(0,1)])  # -> 1",
         upstream=("hessian_cancellation_theorem",),
         ref="studies/v0150_2026aug/arithmetic_area_geometry.py A1-A6"),
    dict(name="cold_amplitude_closed_form_CC04", tier="0", kind="constant",
         primitives=("IV",), tag="Pr", exactness="EXACT",
         nature="c = (9 L5^2)/(25 L2 L3) [1 − (9/5) L5/L3 + (4/5) L5/L2] "
                "= 0.27012646530542495706433719670365 (CC-04; supersedes "
                "cold_dive_law's retracted 16th digit). Cold core "
                "{1,2,3,5}, rate 6/5 = 36/30.",
         example="mtft.curvature.cold_amplitude()",
         upstream=("cold_dive_law",),
         ref="studies/v0150_2026aug/CC-04_cold_amplitude.md"),
    dict(name="weil_compressed_form_W1", tier="5d", kind="identity",
         primitives=("I", "V"), tag="Cert", exactness="E2(3.472e-09)",
         nature="Gabor-compressed Weil form: prime side vs zero side of "
                "the band-limited kernel agree to 3.472e-09 (W1-P1, "
                "pre-registered <= 1e-5); C/N = 0.73398 in [0.729, 0.739] "
                "(W1-P2, anchor 0.734). Lemma 3.2 re-proved independently. "
                "Fourth column of the Three Ensembles table.",
         example="python -m pytest tests/test_weil.py",
         upstream=("hadamard_zetaprime_check",),
         ref="studies/v0150_2026aug/W1_weil_compression_study.md"),
    dict(name="c9b_period_sign_adjudication", tier="1", kind="identity",
         primitives=("II", "III"), tag="Cert", exactness="EXACT + E2(1e-7)",
         nature="Re lambda_1 = −1/2 EXACTLY for the 143a1 period lattice "
                "(Route A: exact rational Manin symbols). Route B "
                "(q-expansion slash-integrals, own a_n sieve, no PARI "
                "mfsymboleval) reproduces v6 per_277 to 1e-7; {1/38,2/77} "
                "= −per_11 to 4e-15 by Manin integrality. v6 sign "
                "CORRECT; Paper 33 v2 archive needs sign-only correction.",
         example="studies/v0150_2026aug/c9b_exact_symbol.py",
         upstream=("X0_143_engine",),
         ref="studies/v0150_2026aug/c9b_exact_symbol.json; c9b_routeB.json"),
)
for _d in V0150_LEGEND:
    if _d["name"] not in REGISTRY:
        _reg(LegendEntry(**_d))

# ═══════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        print(legend_map())
    elif args[0] == "card" and len(args) > 1:
        print(card(args[1]))
    elif args[0] == "trace" and len(args) > 1:
        print(trace(args[1]))
    elif args[0] == "status":
        print(status(args[1] if len(args) > 1 else None))
    elif args[0] == "search" and len(args) > 1:
        print(search(" ".join(args[1:])))
    else:
        print("usage: python -m mtft.legend "
              "[card NAME | trace NAME | status [TAG] | search TERM]")
