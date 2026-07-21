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
