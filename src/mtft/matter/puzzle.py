"""mtft.matter.puzzle — the programme map: every piece of the Standard-Model puzzle on X0(143) x 143a1, its status, and its own puzzles.

A piece is one question the programme asks.  It owns the ledger waves that address it — every ledger claim belongs to exactly one piece: by its
wave, by an explicit override, or (for a correction) by the piece that owns the module that filed it — and it lists the reports in docs/SM that
address it, including those written before the ledger existed.  Each piece carries its own sub-puzzles, and depends only on pieces it uses.

Piece status        stands   its results hold and nothing open threatens them
                    partial  results hold, but sub-puzzles that matter are open
                    open     the piece's question is unanswered (failures may be recorded)
                    closed   investigated to a documented no-go: an honest negative
Sub-puzzle status   open, answered, failed (a decisive or pre-registered test failed), moot (made irrelevant by a closure elsewhere)

Also kept here: the pre-registered tests with their outcomes, the honest negatives, the chronology of the waves, and the disposition of every report
in docs/SM.  This module is bookkeeping over the ledger and the reports; it computes nothing.  tests/test_puzzle_map.py checks it against the ledger,
the modules and the reports; docs/SM/MATTER_PUZZLE_MAP.md is generated from it (python -m mtft matter puzzle --markdown)."""
from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

ORDER = ("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X")      # pieces are numbered in dependency order
PIECE_STATUS = ("stands", "partial", "open", "closed")
SUB_STATUS = ("open", "answered", "failed", "moot")
OUTCOMES = ("pass", "fail", "no hit", "mixed", "pending")


@dataclass(frozen=True)
class Sub:
    """A piece's own puzzle."""
    id: str
    question: str
    status: str
    answer: str = ""                 # answered / failed / moot: what settled it (or why it no longer matters)
    claims: tuple = ()               # ledger ids that settle or frame it
    docs: tuple = ()                 # reports in docs/SM that frame it
    settle: str = ""                 # open: what would settle it
    prereg: str = ""                 # id of the pre-registered test, if one was run on it
    reopen: bool = False             # open sub-puzzle of a closed piece: what would reopen the piece


@dataclass(frozen=True)
class Piece:
    id: str
    name: str
    question: str
    status: str
    verdict: str
    key_claims: tuple
    waves: tuple
    depends_on: tuple
    docs: tuple
    subs: tuple
    survives: tuple = ()             # results that stand inside a closed or open piece


# ------------------------------------------------------------------------------------------------------------------------------------ the pieces
PIECES = (
    Piece("I", "The parent action",
          "Which higher-dimensional theory, compactified on X0(143) x 143a1, is the parent, and is it consistent?",
          "open",
          "Recorded in the registers, not the ledger: no Hom-type unitary-stack flux model with Standard-Model bifundamental families is a consistent "
          "six-dimensional gauge theory (CC-29/R2C-04); the eight-dimensional anomaly gate admits the first chiral, quintic-free contents that land the "
          "Standard Model on X0(143) x 143a1 as net indices (PV-01); supersymmetry of the parent is undeclared (PV-01, PV-02). GF-01 put the three "
          "readings of the question side by side: with the supersymmetric twist the U(8) parent admits the exact Standard Model on its difference lines "
          "(PV-02's conditions, GF-6), while SO(16) and E8 through U(8) are closed — the sum line e_c + e_L doubles SO(16)'s quark-doublet count, E8's "
          "three quark doublets come with three chiral exotic doublets, and in both the lepton rows conflict (GF-11); PV-01's non-supersymmetric U(8) "
          "stands. Among these columns, declaring supersymmetry chooses between the first and the last.",
          (), (), (),
          ("PV01_PARENT_VACUUM_REGISTER.md", "PV02_BLOWUP_SUSY_REGISTER.md", "PARENT_ACTION_REQUIREMENTS.md", "CC29_R2C04_REGISTER.md",
           "R2C03_ANOMALY_LIFT_REGISTER.md", "CC33_TENSOR_COEFFICIENT.md", "WAVE_C_GATES.md", "V0340_PARENT_REVIEW_CORRECTIONS.md",
           "MATTER01_VACUUM_AUDIT.md", "CC35_REFLECTION_POSITIVITY.md", "CC34_DARK_ENERGY_FORMULA_AUDIT.md", "GF01_FORCES_FIRST.md",
           "VE01_VACUUM_ENERGY.md"),
          (Sub("I.1", "Is the parent supersymmetric?", "open", claims=("GF-6", "GF-11"),
               docs=("PV01_PARENT_VACUUM_REGISTER.md", "PV02_BLOWUP_SUSY_REGISTER.md", "GF01_FORCES_FIRST.md"),
               settle="Declare it: it decides whether the U(8) class needs a supersymmetric completion or only its two Green–Schwarz fields (PV-01). "
                      "The choice is now between two columns, since the third is closed: a supersymmetric U(8) with the particles on its difference "
                      "lines (GF-6), or PV-01's non-supersymmetric U(8) with them on the sum lines and the dots; a supersymmetric SO(16) or E8 reached "
                      "through U(8) with the stack hypercharge does not give the exact Standard Model (GF-11)."),
           Sub("I.2", "Global anomalies and the Green–Schwarz construction for the extra U(1)s", "open", claims=("PS1-2",),
               docs=("PV01_PARENT_VACUUM_REGISTER.md", "PARENT_ACTION_REQUIREMENTS.md"),
               settle="A global Green–Schwarz construction for the anomalous stack U(1)s, with the global anomalies checked."),
           Sub("I.3", "Non-split stable bundles with the B1 Chern classes", "open", docs=("PV01_PARENT_VACUUM_REGISTER.md",),
               settle="Decide the mu-stability of a non-split bundle with these classes (PV-01 leaves it undecided)."),
           Sub("I.4", "Strong CP after the lift", "open", claims=("TR-1", "TR-8", "CC-63"), prereg="T6",
               docs=("PV01_PARENT_VACUUM_REGISTER.md", "PARENT_ACTION_REQUIREMENTS.md", "TR01_TIME_REVERSAL.md"),
               settle="The lift must be phase-aligned, or an axion exists (parent-action ledger). Time reversal fixes the arithmetic form of the "
                      "test: the lifted strong coupling must sit on a mirror of the modular tiling whose class matches the quark-mass phase, theta-bar = "
                      "theta + arg det(M_u M_d) = 0; the rule is sealed as PR-S2 (TR-8; PR-S1 was superseded the same day, CC-63) and waits for "
                      "a candidate."),
           Sub("I.5", "Reflection positivity of the lattice action (OS2)", "answered", docs=("CC35_REFLECTION_POSITIVITY.md", "CC69_REPAIRED_ACTION.md"),
               answer="For the programme's action since CC-69 (2026-10-17) — the arithmetic tower on Sym^n characters on plaquettes, CC-35's repair "
                      "adopted — yes, by Osterwalder–Seiler: the weight is a positive-type class function, so link-reflection positivity, a positive "
                      "transfer matrix and a strong-coupling mass gap follow, and centre symmetry is restored. For the original Polyakov power-sum "
                      "placement it stays unproven (CC-35). What remains open is what is open for every lattice action: the continuum limit with a "
                      "uniform gap."),
           Sub("I.6", "A vacuum energy computed from the parent action", "open", claims=("VE-1", "VE-2", "VE-3", "VE-6"),
               docs=("CC34_DARK_ENERGY_FORMULA_AUDIT.md", "PG01_ENERGY_TWO_CLOCKS.md", "VE01_VACUUM_ENERGY.md"),
               settle="The dark-energy formula stays DIAGNOSTIC until it has a route from an action (CC-34). The classical flux and curvature terms "
                      "have no stationary point in the volume (PV-01; the atlas of PG-4 records the runaway), so the route must also stabilise it. "
                      "What the programme's own terms give (VE-01, 2026-10-17): the classical curvature term is positive, 10^105-10^110 rho_Lambda on the "
                      "measured background; the one-loop energy is minus a square, -c sum_odd K(w) |Tr U^w|^2 (VE-1), -3.7 to -4.5 M_KK^4 at the jet "
                      "vacuum (VE-6), and it vanishes exactly at the Bose-Fermi-balanced Z/2 point Tr W = 0 - PA-01's matter-untwisted background (VE-2) - "
                      "which is the global maximum of the Wilson-line potential (VE-3). A vanishing one-loop energy at a vacuum therefore needs the "
                      "Z/2 twist realised without continuous torus Wilson lines (a Green-Schwarz mass for the anomalous T handles two of the fourteen "
                      "unstable directions, I.2; Atkin-Lehner parity none that matter)."),
           Sub("I.7", "What would a big-bang-to-matter run of this vacuum need?", "open", claims=("PG-5",), docs=("PG01_ENERGY_TWO_CLOCKS.md",),
               settle="The gate of PG-01, in order: a declared parent (I.1, now a two-way choice); a vacuum energy from the parent action with a "
                      "stationary point in the volume (I.6; GS-1's fixed Kaehler moduli) — today the classical terms run away and the one-loop Casimir "
                      "has no minimum; the couplings as functions of the moduli at the stabilised point (GC-01, VI.1); a cosmological background to "
                      "supply H(t), since the package's Friedmann module is pheno-tier; then standard physics at confinement and nucleosynthesis, used "
                      "rather than re-derived. What runs today is the shape modulus and the volume's runaway (PG-5, DIAGNOSTIC).")),
          ),
    Piece("II", "The arena",
          "What do the curve, the torus and their product fix before any physics is chosen?",
          "stands",
          "Exact or certified to tolerance: X0(143)(Q) is the four cusps (FA-4); 143a1's lattice has only the rotations +-1 (LB-1) and its even spin "
          "structures have Galois group S3 (M04-2); the canonical metric's distances (SP-1) and Gauss–Bonnet's 12.28 curvature radii per unit length "
          "(FA-3); the tail kernel is a conformal invariant, computed at every home (AG-1, AG-3); the Fricke point maps to 0.222 omega_1 on 143a1, not a "
          "2-torsion point (FH-1); a symmetric bundle of torus degree N has floor((N-1)/2) forced odd modes (OM-1). Among the 105 supersingular "
          "semiprime levels, 13 fold onto a rank-1 elliptic curve (SQ-1); 143 is the only twin-prime one (SQ-3). The thirteen quotients are the "
          "curves of modular degree 4 in their classes, 65a2 among them (TR-3); 143a1's lattice is rhombic, on the theta = pi mirror of time "
          "reversal, like the oloid and triangle mirrors (TR-2), and every quotient's L-function is odd under t -> -t (TR-4). The tile's own clock "
          "is Z/12, the abelianisation of SL(2,Z) (S of order 4, ST of order 6), and its two shortest loops cross the mirrors only at CM points: the "
          "golden loop through i, the twelve loop through rho (GF-12).",
          ("FA-4", "LB-1", "M04-2", "SP-1", "FA-3", "AG-1", "AG-3", "FH-1", "OM-1", "SQ-1", "SQ-3", "TR-3", "TR-2", "TR-4", "GF-12"),
          ("SP-01", "AG-01", "SQ-01", "TR-01"), (),
          ("R2C07_PRODUCT_SURFACE_REGISTER.md", "V0330_COMPENDIUM_NOTES.md"),
          (Sub("II.1", "A second, independent route to the canonical kernel (the E2 rule)", "open", claims=("AG-3",),
               settle="Compute the Arakelov Green's function between the cusps infinity and 0 with the same Ford-domain machinery and compare with the "
                      "Abbes–Ullmo closed form (N squarefree, prime to 6). AG-3's two constructions share the integration, so they are one route."),
           Sub("II.2", "Which points are arithmetically distinguished?", "answered", claims=("FA-4", "FH-1"),
               answer="The four rational cusps (Faltings; Mazur–Kenku) and CM points such as the Fricke point i/sqrt(143)."),
           Sub("II.3", "Does the Atkin–Lehner involution pin the Fricke point's torus image?", "answered", claims=("FH-1", "NU2-1"),
               answer="No: W_143-eigenvalue +1 and root number -1 make phi o W = phi, and the image is 0.222 omega_1."),
           Sub("II.4", "How does curvature change the flat-curve tail?", "answered", claims=("AG-1", "AG-3", "FA-3"),
               answer="The kernel is conformal; at the homes it is 0.60–0.85 of the flat 1/l^2, and Gauss–Bonnet puts them 1.3–2.6 curvature radii away."),
           Sub("II.5", "The jet-to-jet tail beyond the two-point kernel", "open", claims=("AG-3",),
               settle="Differentiate the same Green's function to the jets' orders; needed only if a redesign reuses tails."),
           Sub("II.6", "How many odd modes does topology force under the torus flip?", "answered", claims=("OM-1", "OM-2"),
               answer="floor((N-1)/2) for torus degree N, so none for the design's degree-0 blocks."),
           Sub("II.7", "Is X0(143) the only supersingular semiprime level that folds onto an elliptic curve?", "answered",
               claims=("SQ-1", "SQ-2", "SQ-6"),
               answer="No: 13 of the 105 levels fold onto a rank-1 elliptic curve and 21 onto a sphere; with three or four primes 11 more do up to "
                      "700, each with 8 or 16 cusps."),
           Sub("II.8", "What singles out 143 among them?", "answered", claims=("SQ-3", "SQ-4"),
               answer="Arithmetic only: twin primes, consecutive primes with new dimension p, and the orbit pattern (1, 4, 6); 123 = 3 x 41 shares "
                      "the index, genus, cusps and quotient type."),
           Sub("II.9", "Do the sphere levels follow the moonshine pattern?", "answered", claims=("SQ-5",),
               answer="Yes: every sphere level in the scan is a Monster element order; 143 is not one."),
           Sub("II.10", "Which curves are the thirteen star quotients, and how do they sit under time reversal?", "answered",
               claims=("TR-3", "TR-4", "TR-5", "TR-1"),
               answer="The curve of modular degree 4 in each class (65a2, otherwise the optimal curve); nine have one real circle and four two; "
                      "all are odd under t -> -t. On X(1) the mirror classes theta = 0 and theta = pi never mix; the level-143 real structure "
                      "is finer and still open (TR-5)."),
           Sub("II.11", "What is the clock of the tile, and where does the golden ratio live?", "answered", claims=("GF-12", "CC-64"),
               answer="SL(2,Z)^ab = Z/12, generated by S (order 4, S^2 = -1) and ST (order 6), the character of eta^2; at 143 the twelve is g - 1 = "
                      "(11 + 13)/2. The golden ratio lives at level 5 (X(5) is the icosahedron, twelve cusps at (0, +-1, +-phi)) and on X(1)'s "
                      "systole (trace 3, length 4 log phi); the twelve loop (trace 4, length 2 log(2 + sqrt3)) is the next, with two self-crossings, one "
                      "on each mirror class. Mirror crossings of both are CM points; the golden loop is reversed by the half-turn S, the twelve "
                      "loop only by a mirror.")),
          ),
    Piece("III", "The twist and the families",
          "How does the twist on 143a1 leave three chiral families, and what can they couple to?",
          "stands",
          "An even spin structure removes every adjoint fermion zero mode (M04-1) and the bosonic adjoints get positive one-loop masses (M05-1); the "
          "untwisting Wilson line is a one-loop maximum, so the jet vacuum is stable (PA1-1); full Atkin–Lehner symmetry forbids three families (AL1-1); "
          "one unit of torus flux costs 24 chiral families (TP-1). The families are jets at a blown-up point (CK1-2, with the caveat that the third "
          "cannot be a jet).",
          ("M04-1", "M05-1", "PA1-1", "AL1-1", "TP-1"),
          ("Matter-04", "Matter-05", "PA-01", "AL-01", "SC-01"), ("I", "II"),
          ("SM01_PREREGISTRATION.md", "SM01_REPORT.md", "R2C01_CHIRALITY_REGISTER.md", "R2C08_SURFACE_CHIRALITY_REGISTER.md", "MATTER02_REPORT.md",
           "COHOMOLOGY_DERIVATION.md", "CHARGED_COMPLEX_DERIVATION.md"),
          (Sub("III.1", "Does the twist remove the adjoint fermions?", "answered", claims=("M04-1", "M04-3"),
               answer="Yes: every adjoint fermion zero mode goes; charged fermions of b = 0 blocks are exceptional jets only."),
           Sub("III.2", "Is the jet vacuum stable against untwisting?", "answered", claims=("PA1-1", "PA1-2"),
               answer="Yes: the untwisting Wilson line is a one-loop maximum; the curve copies become messengers."),
           Sub("III.3", "Can Atkin–Lehner symmetry organise the three families?", "answered", claims=("AL1-1", "CC-37"),
               answer="No: full Atkin–Lehner symmetry forbids three families, and one family per cusp is retired for phenomenology."),
           Sub("III.4", "What does torus flux cost?", "answered", claims=("TP-1", "CC-45"),
               answer="24 chiral families per unit (2g - 2 for genus 13); CC-45 retracts the design that assumed one."),
           Sub("III.5", "The physical scales of the compactification and the twist", "open", claims=("M05-2", "SC-1", "SC-2", "SC-3", "SC-4", "SC-5"),
               docs=("SC01_SCALE_LADDER.md",),
               settle="The ratios are fixed (SC-1): M_SS = 0.4943 M_KK, the pair threshold 0.9886 M_KK, and the Wilson-line and Phi masses as "
                      "0.2175, 0.4039 and 0.1324 times sqrt(lambda) M_KK, lambda = g_4^2 r. On the measured background (r = 230-280, M_KK above GC-01's "
                      "matching window) lambda = 72-104 and every one-loop modulus is heavier than the KK scale (SC-2), so the Wilson line decays "
                      "linearly into messenger pairs (SC-3); the closed forms are O(1) uncertain there (SC-4). What stays open is the absolute M_KK: "
                      "bounded below by GC-01, bounded above only once the volume is stabilised (I.6, SC-5)."),
           Sub("III.6", "The excess charged modes of Matter-03's Higgs branch", "open", docs=("CHARGED_COMPLEX_DERIVATION.md",),
               settle="Additional interactions, a changed background or different parent data to lift them (Matter-03).")),
          ),
    Piece("V", "The stack and its breaking",
          "Which gauge group does the stack give, and how is Pati–Salam broken to the Standard Model?",
          "partial",
          "The stack gives the Pati–Salam branching with nonabelian anomalies cancelled and two anomalous U(1)s that need Green–Schwarz masses (PS1-1, "
          "PS1-2); curve degrees alone give D-flat breaking (FS-1, CC-38); the condensate's D-term vacuum is v^2 = (8/13) xi in one orientation (CO-1); "
          "stability needs xi_4L >= (11/13) xi_R4, which orders a_Q > a_u (VS-2, VS-3); vector-like pairs get no tree-level mass (VM-1, VM-2). The "
          "FI-driven condensate needs fixed Kaehler moduli (GS-1). Read as eight dots, the stack's forces are the groupings and its particles the "
          "lines between groups (GF-1); hypercharge follows from anomalies, the normalisation of Q and a neutral Higgs component, by two routes "
          "(GF-4), and the Z6 of the Standard Model acts trivially on every line (GF-10, sealed); with the supersymmetric twist only the U(8) column "
          "of parents admits the exact Standard Model on its lines — SO(16) and E8 through U(8) are closed (GF-6, GF-11).",
          ("PS1-1", "PS1-2", "CO-1", "VS-2", "VS-3", "VM-1", "GS-1", "GF-1", "GF-4", "GF-11"),
          ("PS-01", "toolkit", "FS-01", "CP-01", "CP-02", "CP-03", "CO-01", "VS-01", "GS-01", "VM-01", "GF-01"), ("III", "IV"),
          (),
          (Sub("V.1", "Which phase does the condensate take, coherent or polar?", "open", claims=("CP-2", "CP-11", "CC-40", "LEP-6"),
               settle="The one-loop sign of the spin-dependent quartic; the equal-overlap fits cannot decide (CC-40, LEP-6)."),
           Sub("V.2", "Is the condensate's orientation unique?", "answered", claims=("CO-1",),
               answer="Yes: one orientation, v^2 = (8/13) xi; a conjugate condensate is lifted by any positive mass."),
           Sub("V.3", "Is the breaking vacuum stable?", "answered", claims=("VS-2", "VS-3", "CC-49"),
               answer="Only if xi_4L >= (11/13) xi_R4, which orders the curve degrees a_Q > a_u."),
           Sub("V.4", "How do the vector-like pairs get their masses?", "open", claims=("VM-1", "VM-2", "CO-2", "CO-4", "GS-3"),
               settle="The one-loop Dirac masses with Scherk–Schwarz gauginos and D-term-split partners; gauge loops alone give m_R/m_L = 0.90."),
           Sub("V.5", "Do the anomalous U(1)s get consistent Green–Schwarz masses?", "open", claims=("PS1-2", "GS-2"),
               docs=("PV01_PARENT_VACUUM_REGISTER.md",), settle="The global Green–Schwarz construction of the parent (I.2)."),
           Sub("V.6", "Does Pati–Salam breaking need the gamma direction?", "answered", claims=("CC-38", "FS-2"),
               answer="No: curve degrees suffice, and a global Higgs with a single D-flat condensate is compatible."),
           Sub("V.7", "Is charge quantised, and what fixes hypercharge?", "answered", claims=("GF-10", "GF-4", "CC-64"), prereg="T7",
               answer="The Z6 of [SU(3) x SU(2) x U(1)]/Z6 acts trivially on every line (sealed PR-S5; PR-S4 superseded the same day, CC-64; PASS), "
                      "so the matter allows G/Z6 or any of its covers and free particles carry whole charges, "
                      "Q = (2/3) triality mod 1. Anomalies alone leave <T_3R, B-L>; Y(Q) = 1/6 and a neutral Higgs component fix Y = T_3R + (B-L)/2 "
                      "up to the u <-> d relabelling, and the Pati–Salam route agrees (E2)."),
           Sub("V.8", "Which parent's lines can carry three families?", "answered", claims=("GF-6", "GF-11"), docs=("GF01_FORCES_FIRST.md",),
               answer="With the supersymmetric twist the chirality of a block is linear in the stack twist indices. U(8) with difference lines: yes, "
                      "PV-02's three conditions with six forced doublet pairs. SO(16) or E8 through U(8): not the exact Standard Model — SO(16)'s "
                      "quark-doublet count is 2 phi_c, E8's three quark doublets come with three chiral exotic doublets, and in both the lepton rows "
                      "force net(L) + net(e^c) = 0; inconsistent at every hypercharge shift. PV-01's non-supersymmetric U(8) is the other open column.")),
          ),
    Piece("IV", "Size and home",
          "What fixes the area ratio r = A_X/A_E, and where on the curve does the top live?",
          "partial",
          "Nothing fixes r: for genus >= 2 the one-loop Casimir energy has no minimum (MS-6), nothing in the model stabilises the Kaehler moduli (MS-7), "
          "and no canonical prescription hits the pre-registered window (CS-3). r is measured: four windows meet at 230–280 with the top at the Fricke "
          "point (FH-3), the only arithmetic home that passes the pre-registered kernel test (AG-4) — modest evidence (AG-5). Why the level is 143 is "
          "itself open (IV.6): the arithmetic singles it out without a physical reason, and a selection rule is pre-registered (SQ-8); the "
          "trials counter records seventeen properties tried, four singling 143 out, all found after the fact (PR-3, GF-9). Read as a dynamical "
          "system, the shape r sits in a cosh well with its floor at the self-dual point (PG-1) while the volume has no critical point at all (PG-4); "
          "the atlas's saddles are the two passes of the Wilson-line torus (PG-6), the half-periods at which BL-01's B–L line also sits (BL-3; "
          "CC-66 corrected the first filing's 'no saddle' on 2026-10-14), and the spiral type of a cosmological run is set by m/H, which nothing "
          "fixes (PG-5; CC-65 corrected the wave's first filing the same day, and CC-67 (2026-10-16) supplied the hypotheses its capture argument needs).",
          ("MS-6", "MS-7", "CS-3", "AG-4", "AG-5", "FH-3", "SQ-8", "PG-1", "PG-4", "PG-6"),
          ("MS-01", "MS-02", "CS-01", "FH-01", "PR-01", "PG-01", "VE-01"), ("II", "III"),
          (),
          (Sub("IV.1", "A mechanism for r = A_X/A_E", "open", claims=("MS-6", "MS-7", "CS-1", "CS-3"), prereg="T1",
               settle="Flux quantisation of the Kaehler class or a non-perturbative potential; nothing in the present model supplies one."),
           Sub("IV.2", "What stabilises the Kaehler moduli at nonzero FI?", "open", claims=("GS-1", "MS-7"),
               settle="A stabilisation compatible with the FI-driven condensate; with dynamical moduli the condensate relaxes away (GS-1)."),
           Sub("IV.3", "Where does the top live?", "answered", claims=("AG-4", "AG-5", "FA-4"), prereg="T3",
               answer="At the Fricke point i/sqrt(143): the only arithmetic home whose kernel puts r from V_cb inside the Cabibbo window; one pass in "
                      "four after a post-hoc flag, so modest evidence."),
           Sub("IV.4", "Do independent measurements of r agree?", "answered", claims=("FH-3", "CC-58"),
               answer="Yes at r = 230–280; the V_us window is order-of-magnitude only and the windows are broad."),
           Sub("IV.5", "Why did V_cb and V_us disagree on r by a factor 3.3?", "answered", claims=("FA-2", "CC-55", "CC-56"),
               answer="Single-image tails; image sums close the gap and the curved kernel refines it."),
           Sub("IV.6", "Why this level: what physical rule, if any, selects 143 — for instance by preferring twin primes?", "open",
               claims=("SQ-3", "SQ-4", "SQ-8", "PR-3", "GF-9"), docs=("PR01_PREREGISTRATION.md",),
               settle="A rule sealed in mtft.matter.prereg before it is run on the 13 torus levels (SQ-8, PR-1): it must select 143 and reject the "
                      "other twelve, 123 (same index, genus and cusps) and 77 (consecutive primes) included. Seventeen properties have been tried "
                      "after the fact and four single 143 out (PR-3; X17, the twelve-loop systole, selects 118, 141, 142 and 143, GF-9); none counts."),
           Sub("IV.7", "What kind of critical point is each modulus of the vacuum?", "answered", claims=("PG-1", "PG-4", "PG-5", "CC-65", "PG-6", "CC-66", "CC-67"), docs=("PG01_ENERGY_TWO_CLOCKS.md",),
               answer="Shape: a sink in a cosh well, floor at r* = sqrt(19/11), curvature 1 in floor units (PG-1). Volume: no finite critical point, a "
                      "runaway, given PV-01's positive coefficients (PG-4). Condensate: a degenerate line lifted to a sink by a mass term with "
                      "0 < m^2 < xi (eigenvalue 4 m^2). Wilson line: a 2D sink at the jet vacuum, a 2D source at the untwisted point and, on the "
                      "whole torus, two saddles at the other half-periods (at V = 0 exactly in the real spin structure), exact in PA-01's product "
                      "approximation; its masses follow the curvatures in 143a1's kinetic metric (PG-6). Adjoints: sinks; A_X at one loop: no "
                      "critical point. Spirals only with time, when m > 3H/2; capture of the shape follows from the energy argument when H is "
                      "bounded, of bounded variation and not integrable — constant H and the radiation and matter eras — not for every H > 0, and "
                      "the toy run only illustrates it (PG-5, CC-67). CC-65 corrected the first filing (hypotheses, the run-length-dependent count, and the arithmetic "
                      "half's RH sentence) the same day; CC-66 (2026-10-14) corrected its 'no saddle', which held only along the straight path "
                      "and conflicted with BL-3; CC-67 (2026-10-16) added the hypotheses of the capture argument.")),
          ),
    Piece("VI", "Couplings, the Pati–Salam scale and b–tau",
          "Do the gauge couplings and the bottom–tau relation agree with a Pati–Salam scale this geometry allows?",
          "open",
          "Not as things stand. Flux designs quantise R to 1/9 ... 9 (two-point) or 1/5, 1, 5, never the data's 3.678 (GC-3); the allowed matching "
          "scales 8.9e11–2.1e14 GeV leave y_b/y_tau = 0.66–0.71, so b–tau needs a 30–35 per cent threshold (GC-4); thresholds need about 17 light "
          "pairs (TH-2, TH-3). Couplings fix V/kappa, never the shape A_X/A_E (GC-1). On the force torus the couplings are lengths and the matching "
          "is Pythagoras (GF-5). GC-02 (2026-10-17) separates the two scales GC-01 had identified: matched at the curve's compactification scale "
          "M_C above M_PS, with Pati–Salam running between, R rises and the designs R = 5 and 9 are met at M_C = 1.7e7 and 2.8e8 GeV with M_PS at "
          "the b–tau crossing — no threshold correction needed on those designs (GC-6, GC-7).",
          ("GC-1", "GC-3", "GC-4", "TH-2", "TH-3", "GF-5", "GC-6"),
          ("GC-01", "TH-01", "GC-02"), ("V",),
          ("PHOTON_KERNEL_REPORT.md", "KKT01_MAGNETIC_TOWERS_REPORT.md"),
          (Sub("VI.1", "Can thresholds reconcile R with the allowed matching scales?", "open", claims=("TH-1", "TH-2", "TH-3", "CA-5"),
               settle="A mass source for about 17 right-handed pairs (gamma-charge 2); with charge 1 they sit near 3 TeV."),
           Sub("VI.2", "The Kaluza–Klein towers' beta coefficients", "open", claims=("TH-4",),
               settle="The curve's spectral invariant against the torus's DKL factor 0.3527: the only place the shape enters the couplings."),
           Sub("VI.3", "The b–tau normalisation at matching", "open", claims=("GC-4", "QL-1", "GC-7"),
               settle="Under GC-01's single-scale matching, a 30–35 per cent threshold correction at 8.9e11–2.1e14 GeV (the data want y_tau/y_b = "
                      "1.40–1.53 there); under GC-02's two-scale matching with M_PS at the b–tau crossing, none — the question becomes which design "
                      "(R = 5 or 9) the background realises and what the thresholds at M_C do (VI.2)."),
           Sub("VI.4", "Is the Pati–Salam scale where y_b = y_tau (1.1e6 GeV)?", "open", claims=("CC-48", "GC-4", "GC-6", "GC-7"),
               docs=("GC02_TWO_SCALES.md",),
               settle="Not if the couplings are matched at the Pati–Salam scale itself (GC-01: 8.9e11–2.1e14 GeV). It can be if the flux splitting is "
                      "matched at the curve's compactification scale M_C above M_PS (GC-6): the designs R = 5 and 9 then meet it at M_C = 1.7e7 and "
                      "2.8e8 GeV with M_PS = 1.1e6 GeV (GC-7). To settle: the background's design, the thresholds at M_C, and the K_L -> mu e bound "
                      "on the leptoquarks, which sits near 1e6 GeV."),
           Sub("VI.5", "Does an SU(2)_R-breaking Wilson-line split move the matching?", "open", claims=("GT-4", "QL-4"),
               settle="Its effect on GC-01's matching; its lepton purpose is closed (QL-4)."),
           Sub("VI.6", "Are the couplings the shape of the force torus?", "answered", claims=("GF-5", "GC-3"),
               answer="Yes: lengths on the torus are 1/g, the photon is the diagonal with 1/e^2 = 1/g^2 + 1/g'^2 and the weak angle its tilt; "
                      "GC-01's matching is the same theorem on the Pati–Salam torus, and R = 3.678 is a ratio of cycle-length differences the "
                      "flux designs quantise away from the data (GF-5, GC-3).")),
          ),
    Piece("VII", "Higgs, vertex and the top",
          "Which Higgs couples the families, and can any configuration give an order-one top Yukawa?",
          "partial",
          "An A-type Higgs gives no tree-level Yukawa with A-type families (HP-3, CC-43); holomorphic messengers are forbidden at the family vertex "
          "(HJ-1); the coupling Higgs is a torus-constant Phi-type pullback (HP-4, TP-4); the character audit leaves a rank-2 texture (YA-3); lattice "
          "breaking gives jet Yukawas <= 0.01 and no jet can host the top (LB-3, TF-1, LB-4). A constant-jet top at a second point is marginal (TP-2), "
          "and a bridged top fails its pre-registered test (YB-1).",
          ("HP-3", "HJ-1", "HP-4", "YA-3", "LB-3", "TF-1", "TP-2", "YB-1"),
          ("YK-01", "PB-01", "XQ-01", "HJ-01", "HP-01", "YA-01", "LB-01", "TF-01", "TP-01", "YB-01"), ("III", "V"),
          ("R2C02_VECTOR_HIGGS_REGISTER.md", "R2C05_POLARISATION_REGISTER.md", "R2C06_TELESCOPING_REGISTER.md", "CW01_ONE_LOOP_POTENTIAL.md",
           "HIGGS_VACUUM_DERIVATION.md", "MATTER03_REPORT.md", "MATTER03_PHYSICAL_INPUTS.md", "YUKAWA_DERIVATION.md", "SM02_STEP5_REPORT.md",
           "SM03_UP_YUKAWA_REPORT.md", "KKT02_CERTIFICATION_AND_FEM_YUKAWA.md", "KINETIC_NORMALIZATION_REPORT.md", "CC36_JET_YUKAWA_RETRACTION.md"),
          (Sub("VII.1", "Which Phi-type Higgs combination is light and takes the vev?", "open", claims=("HP-4", "TP-4"),
               settle="The Higgs potential on the Phi-type space: thirteen cusp forms, or one constant when a_Q + a_u = 24."),
           Sub("VII.2", "An order-one top Yukawa", "open", claims=("TF-1", "TF-2", "LB-4", "TP-2", "YB-1"), prereg="T5",
               settle="A channel that is neither a jet (TF-1) nor a bridge (YB-1); the constant-jet top at a second point is only marginal (TP-2)."),
           Sub("VII.3", "Can families couple through holomorphic messengers?", "answered", claims=("HJ-1", "PB1-1", "CC-42"),
               answer="No: the vertex charge is never zero, and the SU(4) invariants vanish."),
           Sub("VII.4", "Does the ratio-3 background carry a Yukawa?", "answered", claims=("HP-2", "HP-3", "CC-43"),
               answer="No: its Higgs is A-type, and A-type x A-type x A-type vanishes."),
           Sub("VII.5", "The heavy bidoublet spectrum behind Higgs–bidoublet mixing", "open", claims=("PB1-2", "CC-44"),
               settle="The first C-excited modes on X0(143) and their overlaps; CC-44 caveats the leading amplitude assumed."),
           Sub("VII.6", "Is there a leading-order electroweak hierarchy?", "answered",
               docs=("R2C06_TELESCOPING_REGISTER.md", "R2C07_PRODUCT_SURFACE_REGISTER.md"),
               answer="Not on a curve (chirality and Higgs mass are one number); on the product surface the electroweak scale is a Kaehler-modulus "
                      "deviation, so its size waits on IV."),
           Sub("VII.7", "Does one-loop physics select the Higgs direction?", "open", docs=("CW01_ONE_LOOP_POTENTIAL.md",),
               settle="CW-01's computation failed; a continuum-limit renormalisation at fixed cutoff is needed.")),
          ),
    Piece("IX", "Quark flavour",
          "Can the two-point design — two jet families at a blown-up cusp, the top at a second point — produce the CKM matrix and the quark mass ratios?",
          "closed",
          "No. The condensate reaches only u^c, e^c, nu^c, L and H_u (DS-1), so the down sector is parity-symmetric at first order and V_ub has a floor "
          "(DS-2); with tails capped at 13.8 the design predicts |V_ub| >= 0.0084, 2.3 times the data (DS-4). No odd source can reach the down sector, "
          "by topology (PA-1, OM-2), and the bridged top fails (YB-1). The design is closed (YB-2). GF-01 made the down-sector wall parent-level within "
          "the supersymmetric menus: no Standard-Model-singlet line of U(8) reaches d^c at first order, and the sum lines and sheets that do belong "
          "to SO(16) and E8, which do not give the exact Standard Model (GF-13, GF-11).",
          ("DS-1", "DS-2", "DS-4", "PA-1", "OM-2", "YB-2", "GF-13"),
          ("CK-01", "EPS-01", "RT-01", "TA-01", "DD-01", "CA-01", "TO-01", "FL-01", "CD-01", "DP-01", "FA-01", "LV-01", "TL-01", "DS-01", "PAR-01",
           "OM-01"),
          ("V", "IV", "VII"),
          ("SM04_TEXTURE_REPORT.md", "SM08_HYM_REPORT.md", "SM09_CERTIFICATION_REPORT.md", "SM11_FLUX_SCAN_REPORT.md", "SM13_MODEL_M2_REPORT.md",
           "SM15_W13_DEGENERACY_REPORT.md", "CC26_NORMALISATION_RETRACTION.md", "R2C09_THETA_TORUS_REGISTER.md"),
          (Sub("IX.1", "Does the torus parity allow the observed mixing?", "answered", claims=("DD-1", "DD-2", "DD-3"),
               answer="Only broken, with the odd family first: unbroken, the odd family decouples; broken by eps_P ~ 0.225 it gives the Wolfenstein pattern."),
           Sub("IX.2", "V_cb from the sectors' own mass ratios", "failed", claims=("TO-2",), prereg="T2",
               answer="0.078–0.199 for every phase, against the window 0.030–0.055."),
           Sub("IX.3", "Which tail sets which angle, and can curve flux make the asymmetry?", "failed", claims=("FL-1", "FL-2", "FL-3", "CC-51"),
               answer="The CKM angle is set by the right-handed block's tail; the data need T(Q):T(d^c):T(u^c) = 40:5.3:1 and curve flux reaches 13.8."),
           Sub("IX.4", "Can up-quark dressing supply the asymmetry?", "moot", claims=("CD-3", "DP-3", "YB-2"),
               answer="z = 0.25 needs a strong x-jet component and a gaugino tower near M_PS/2; the design is closed."),
           Sub("IX.5", "V_ub's coefficient c", "failed", claims=("TL-3", "TL-4", "DS-2"), prereg="T4",
               answer="c is the down quarks' share of V_cb; with the programme's tails and dressing c >= 0.565, and any parity-symmetric down sector "
                      "with theta_R <= 45 degrees has c >= 0.47."),
           Sub("IX.6", "An odd source in the down sector", "answered", claims=("PA-1", "OM-2", "OM-3"),
               answer="None: no odd source reaches Q, d^c or H_d, and topology forbids one on the B stack short of 72 extra families."),
           Sub("IX.7", "Where does the Cabibbo angle come from?", "answered", claims=("DS-1", "TL-2"),
               answer="From the up sector: the down sector has no first-order parity breaking, and a common gamma tilt cancels."),
           Sub("IX.8", "The CP phase gamma", "moot", claims=("LV-3", "YB-2"),
               answer="Real arms predict 17 or 141 degrees against 65.7; a twist between the arms fits; the design is closed."),
           Sub("IX.9", "What would reopen quark flavour?", "open", claims=("TP-1", "OM-3", "YB-2", "SQ-7"), reopen=True,
               settle="A geometry where the torus varies over the curve (an elliptic fibration rather than a product) changes the 24-family price of "
                      "torus flux and with it the odd-mode count; scope it before building. Another level would not: the exact obstructions do not "
                      "use it (SQ-7)."),
           Sub("IX.10", "Can any parent's condensate reach the down sector at first order?", "answered", claims=("GF-13", "GF-11", "DS-1"), prereg="T8",
               answer="Not in the supersymmetric menus. The sealed touch rule (PR-S3) finds no Standard-Model-singlet line of U(8) that reaches d^c "
                      "or H_d; the sum lines and sheets that do belong to SO(16) and E8, which do not give the exact Standard Model (GF-11). PV-01's "
                      "non-supersymmetric U(8) has the sum lines as fermions only and was not scanned. Within the supersymmetric menus the wall is "
                      "parent-level, not a design accident (GF-13).")),
          ("DD-1", "DD-2", "LV-1", "TL-2", "TL-3", "DS-2", "FH-2", "FL-1"),
          ),
    Piece("X", "Lepton flavour",
          "What makes the lepton masses differ from the down-quark masses in the pattern the data show?",
          "closed",
          "The data ask lepton/quark factors 3.09 and 0.30 relative to the third family — Georgi–Jarlskog's (3, 1/3, 1) within 3 and 10 per cent "
          "(QL-2) — that is, lepton tails 1.76 times quark tails with local entries untouched (CI-3). Kinetic dressing only lowers Yukawas (QL-1), the "
          "condensate channel's lever arm is too short (CI-4), and a B–L Wilson line gives the pattern or V_cb, never both, at a one-loop saddle (BL-2, "
          "BL-3): no mechanism in the current matter sector (BL-5).",
          ("QL-2", "CI-3", "QL-1", "CI-4", "BL-2", "BL-3", "BL-5"),
          ("LEP-01", "GT-01", "QL-01", "CI-01", "BL-01"), ("VI", "VII", "IX"),
          ("SM05_TWIST_STATUS.md", "SM06_TWIST_REPORT.md", "SM07_CUBIC_REPORT.md", "SM10_DOWN_LEPTON_REPORT.md", "SM12_TWIST_ASSIGNMENTS_REPORT.md"),
          (Sub("X.1", "Kinetic dressing by the condensate", "failed", claims=("QL-1", "GT-3", "CC-52"),
               answer="Dressing only lowers Yukawas: z_e = 3.09 forces y_tau/y_b <= 0.32 against 1.40–1.53."),
           Sub("X.2", "A condensate insertion at the vertex", "failed", claims=("CI-2", "CI-4", "CC-53"),
               answer="Tail-carried entries are SU(4)-universal, and the condensate channel's lever arm is 2.4–3.0 against 5.9–7.7 needed."),
           Sub("X.3", "A B–L Wilson line", "failed", claims=("BL-2", "BL-3", "CC-54"),
               answer="Lepton pattern or V_cb, never both, and the vacuum is a one-loop saddle."),
           Sub("X.4", "Does the equal-overlap texture force the charge-1 jet for leptons?", "answered", claims=("LEP-6", "CC-61"),
               answer="No: every sector has several exact fits, and the leptons also fit with the charge-0 jet dominant."),
           Sub("X.5", "Koide's Q = 2/3", "answered", claims=("LEP-5",),
               answer="Recorded as numerology: the texture does not single out 2/3."),
           Sub("X.6", "Is the predicted relation f1 f2 = f3^2 realised?", "open", claims=("CI-3",), reopen=True,
               settle="Data give 0.919 against 1, inside m_d/m_s's 8 per cent; any future mechanism must make lepton tails 1.76 times quark tails."),
           Sub("X.7", "Did the curve-only models reach Georgi–Jarlskog?", "answered",
               docs=("SM06_TWIST_REPORT.md", "SM07_CUBIC_REPORT.md", "SM12_TWIST_ASSIGNMENTS_REPORT.md"),
               answer="No: quadratic twists cannot cure it, the cubic twist lifts the down–lepton degeneracy, and the specific factor is out of reach of "
                      "that class.")),
          ("QL-2", "CI-3", "LEP-2", "LEP-6"),
          ),
    Piece("VIII", "Neutrinos",
          "How do the neutrinos get their small masses here?",
          "open",
          "Type-I seesaw needs M_R(3) = 2.9e14 GeV, excluded; an inverse seesaw needs mu(3) = 4.3 MeV, mu/M = 3.8e-9 (NU1-3). Gamma-localised bosons "
          "and fermions sit on disjoint sets (NU1-1). In the ratio-3 vacuum the jet-sector neutrinos come out 10^9 too heavy (RT-3); mu/M ~ "
          "exp(-beta sqrt r) (NU2-3, a caveat).",
          ("NU1-1", "NU1-3", "RT-3", "NU2-3"),
          ("NU-01", "NU-02"), ("V", "IV"),
          ("SM14_HIGGS_CKM_GAUGE_REPORT.md",),
          (Sub("VIII.1", "The B–L-conserving Dirac pairing M", "open", claims=("NU1-3", "NU2-3"),
               settle="M has no singlet-scalar vev at tree level; whether it is local to N decides mu/M ~ exp(-beta sqrt r)."),
           Sub("VIII.2", "Majorana or Dirac, and the PMNS structure", "open", claims=("NU1-3",), docs=("SM14_HIGGS_CKM_GAUGE_REPORT.md",),
               settle="Y_nu = Y_u has small mixing; large PMNS angles need another source."),
           Sub("VIII.3", "Can the jet sector supply light neutrinos?", "failed", claims=("RT-3",),
               answer="Ten to the ninth too heavy in the ratio-3 vacuum."),
           Sub("VIII.4", "A neutrino mechanism in the two-point design", "moot", claims=("YB-2",),
               answer="The two-point design is closed for quark flavour.")),
          ),
)

# claims whose subject is another piece's than their wave's (the wave decides everything else)
CLAIM_PIECE = {
    "M04-2": "II", "M05-3": "II", "NU2-1": "II", "NU2-2": "II", "LB-1": "II", "LB-2": "II", "TA-1": "II", "GC-5": "II", "FA-3": "II", "FA-4": "II",
    "FH-1": "II", "OM-1": "II",
    "TP-1": "III", "CK1-2": "III", "CC-37": "III",
    "SP-2": "IV", "SP-4": "IV", "AG-4": "IV", "AG-5": "IV", "FA-2": "IV", "BL-4": "IV", "CA-1": "IV", "CA-2": "IV", "CA-3": "IV",
    "CC-55": "IV", "CC-56": "IV", "CC-68": "IV",
    "PS1-4": "VI", "CA-5": "VI",
    "FH-2": "IX", "FH-4": "IX", "YB-2": "IX", "CC-39": "IX", "CC-47": "IX", "CC-60": "IX",
    "CD-4": "X", "DP-4": "X",
    "RT-3": "VIII",
    "SQ-7": "IX", "SQ-8": "IV",
    "GF-2": "IX", "GF-5": "VI", "GF-8": "II", "GF-9": "IV", "GF-12": "II", "GF-13": "IX",
}

# ------------------------------------------------------------------------------------------------------------------------------------ the record
WAVE_DATES = {
    "Matter-04": "2026-10-02", "Matter-05": "2026-10-02", "PS-01": "2026-10-02", "NU-01": "2026-10-02", "NU-02": "2026-10-02", "AL-01": "2026-10-02",
    "YK-01": "2026-10-02", "PA-01": "2026-10-02", "PB-01": "2026-10-02", "CK-01": "2026-10-02", "EPS-01": "2026-10-02", "toolkit": "2026-10-02",
    "FS-01": "2026-10-02", "RT-01": "2026-10-02", "LEP-01": "2026-10-02", "CP-01": "2026-10-03", "CP-02": "2026-10-03", "CP-03": "2026-10-03",
    "XQ-01": "2026-10-03", "HJ-01": "2026-10-04", "HP-01": "2026-10-04", "YA-01": "2026-10-04", "LB-01": "2026-10-04", "TF-01": "2026-10-05",
    "TP-01": "2026-10-05", "TA-01": "2026-10-05", "SP-01": "2026-10-05", "GC-01": "2026-10-05", "TH-01": "2026-10-05", "VM-01": "2026-10-05",
    "CO-01": "2026-10-05", "VS-01": "2026-10-05", "GS-01": "2026-10-05", "MS-01": "2026-10-05", "MS-02": "2026-10-05", "CS-01": "2026-10-05",
    "DD-01": "2026-10-05", "CA-01": "2026-10-05", "TO-01": "2026-10-05", "FL-01": "2026-10-05", "CD-01": "2026-10-05", "DP-01": "2026-10-05",
    "GT-01": "2026-10-05", "QL-01": "2026-10-06", "CI-01": "2026-10-06", "BL-01": "2026-10-06", "FA-01": "2026-10-06", "AG-01": "2026-10-07",
    "FH-01": "2026-10-07", "LV-01": "2026-10-07", "TL-01": "2026-10-07", "DS-01": "2026-10-07", "PAR-01": "2026-10-08", "OM-01": "2026-10-08",
    "YB-01": "2026-10-08", "SQ-01": "2026-10-08", "TR-01": "2026-10-10", "PR-01": "2026-10-10", "GF-01": "2026-10-11", "PG-01": "2026-10-12", "SC-01": "2026-10-17", "VE-01": "2026-10-17", "GC-02": "2026-10-17",
}

WAVE_REPORT = {
    "Matter-04": "MATTER04_ADJOINT_SPINORS.md", "Matter-05": "MATTER05_ONE_LOOP_ADJOINTS.md", "PS-01": "PS01_PATI_SALAM_BREAKING.md",
    "NU-01": "NU01_NEUTRINO_GATE.md", "NU-02": "NU02_NEUTRINO_MU.md", "AL-01": "AL01_CUSP_FAMILIES.md", "YK-01": "YK01_YUKAWA_JETS.md",
    "PA-01": "PA01_PATH_A.md", "PB-01": "PB01_MESSENGER_YUKAWA.md", "CK-01": "CK01_CKM_CUSPS.md", "EPS-01": "EPS01_EPSILON.md",
    "toolkit": "MATTER_PROGRAM_MAP.md", "FS-01": "FS01_FLUX_SCAN.md", "RT-01": "RT01_RATIO_THREE.md", "LEP-01": "LEP01_LEPTONS.md",
    "CP-01": "CP01_CONDENSATE_PHASES.md", "CP-02": "CP01_CONDENSATE_PHASES.md", "CP-03": "CP01_CONDENSATE_PHASES.md", "XQ-01": "XQ01_OVERLAPS.md",
    "HJ-01": "HJ01_JET_VERTEX.md", "HP-01": "HP01_HIGGS_PROFILE.md", "YA-01": "YA01_YUKAWA_AUDIT.md", "LB-01": "LB01_LATTICE_BREAKING.md",
    "TF-01": "TF01_TOP_FAMILY.md", "TP-01": "TP01_TOP_POINT.md", "TA-01": "TA01_TWO_POINT_AUDIT.md", "SP-01": "SP01_SPECIAL_POINTS.md",
    "GC-01": "GC01_GAUGE_COUPLINGS.md", "TH-01": "TH01_THRESHOLDS.md", "VM-01": "VM01_VECTOR_LIKE_MASSES.md", "CO-01": "CO01_CONJUGATE_ORIENTATION.md",
    "VS-01": "VS01_VACUUM_STABILITY.md", "GS-01": "GS01_GREEN_SCHWARZ.md", "MS-01": "MS01_MODULI_CASIMIR.md", "MS-02": "MS01_MODULI_CASIMIR.md",
    "CS-01": "CS01_CANONICAL_SHAPE.md", "DD-01": "DD01_DOUBLE_DUTY.md", "CA-01": "CA01_CABIBBO.md", "TO-01": "TO01_THIRD_OBSERVABLE.md",
    "FL-01": "FL01_FLUX_TAILS.md", "CD-01": "CD01_CONDENSATE_DRESSING.md", "DP-01": "DP01_DRESSING_PROFILES.md", "GT-01": "GT01_GAUGINO_TOWERS.md",
    "QL-01": "QL01_QUARK_LEPTON_SPLIT.md", "CI-01": "CI01_CONDENSATE_INSERTION.md", "BL-01": "BL01_BL_WILSON_LINE.md",
    "FA-01": "FA01_FLAVOUR_AUDIT.md", "AG-01": "AG01_ARAKELOV_KERNEL.md", "FH-01": "FH01_FRICKE_HOME.md", "LV-01": "LV01_LEVER.md",
    "TL-01": "TL01_TILT.md", "DS-01": "DS01_DOWN_PARITY.md", "PAR-01": "PA01_PARITY_AUDIT.md", "OM-01": "OM01_YB01_EXPERIMENTS.md",
    "YB-01": "OM01_YB01_EXPERIMENTS.md", "SQ-01": "SQ01_STAR_QUOTIENTS.md", "TR-01": "TR01_TIME_REVERSAL.md", "PR-01": "PR01_PREREGISTRATION.md",
    "GF-01": "GF01_FORCES_FIRST.md", "PG-01": "PG01_ENERGY_TWO_CLOCKS.md",
 "SC-01": "SC01_SCALE_LADDER.md", "VE-01": "VE01_VACUUM_ENERGY.md", "GC-02": "GC02_TWO_SCALES.md",
}

PREREG = (
    {"id": "T0", "wave": "SM-01", "registered": "2026-09-13", "piece": "III", "test": "Three chiral families from four or five stacks on X0(143)",
     "criterion": "The list of solutions or the exact negative, stated before the search", "outcome": "mixed",
     "result": "Four stacks: no model gives three chiral families; five stacks: 30 solutions, all with m_L = -3", "claims": (),
     "docs": ("SM01_PREREGISTRATION.md", "SM01_REPORT.md")},
    {"id": "T1", "wave": "CS-01", "registered": "2026-10-05", "piece": "IV", "test": "A canonical prescription for A_X/A_E",
     "criterion": "A_X/A_E inside V_cb's window [18.1, 36.4], four families of prescriptions", "outcome": "no hit",
     "result": "4, 1, and 44.3 or 22.2 (straddling); the genus-type counts 24–28 have no derivation", "claims": ("CS-1", "CS-2", "CS-3"), "docs": ()},
    {"id": "T2", "wave": "TO-01", "registered": "2026-10-05", "piece": "IX", "test": "V_cb from the sectors' own mass ratios, no fitted parameter",
     "criterion": "V_cb in [0.030, 0.055]", "outcome": "fail", "result": "0.078–0.199 for every phase: a factor 1.9", "claims": ("TO-2",), "docs": ()},
    {"id": "T3", "wave": "AG-01", "registered": "2026-10-06", "piece": "IV", "test": "r from V_cb at the arithmetic homes, with the canonical kernel",
     "criterion": "r (V_cb 0.030–0.055) inside the Cabibbo window 91–896", "outcome": "pass",
     "result": "Fricke point 167–306 passes; cusps 0, 1/11, 1/13 fail; the tau/2 preimage (not arithmetic) touches the edge. The Fricke point had "
               "been flagged post hoc, so the evidence is modest (AG-5)", "claims": ("AG-4", "AG-5"), "docs": ()},
    {"id": "T4", "wave": "TL-01", "registered": "2026-10-07", "piece": "IX", "test": "V_ub's coefficient c from the gamma tilt",
     "criterion": "c in [0.35, 0.45]; second target gamma = 65.7 +- 3 degrees", "outcome": "fail",
     "result": "c >= 0.565 with DP-01's dressing; c >= 0.47 for any parity-symmetric down sector with theta_R <= 45 degrees; gamma never computed "
               "because the design closed", "claims": ("TL-4", "DS-2"), "docs": ()},
    {"id": "T5", "wave": "YB-01", "registered": "2026-10-08", "piece": "VII", "test": "The bridged top Yukawa",
     "criterion": "y_t >= 0.5", "outcome": "fail", "result": "0.29–0.32 with a flat tail, 0.17–0.27 with curvature", "claims": ("YB-1",), "docs": ()},
    {"id": "T6", "wave": "TR-01", "registered": "2026-10-10", "piece": "I", "test": "Strong CP after the lift (sealed as PR-S2, SHA-256 in mtft.matter.prereg; PR-S1 superseded the same day, CC-63)",
     "criterion": "theta-bar = theta + arg det(M_u M_d) = 0 mod 2 pi at a CM point of the modular tiling: a theta = 0 mirror with arg det = 0, "
                  "or a theta = pi mirror with arg det = pi; theta-bar = pi fails; no mirror fails unless an axion exists", "outcome": "pending",
     "result": "No candidate yet: the lift (PV-01) is open", "claims": ("TR-8", "PR-1", "CC-63"), "docs": ("TR01_TIME_REVERSAL.md", "PR01_PREREGISTRATION.md")},
    {"id": "T7", "wave": "GF-01", "registered": "2026-10-11", "piece": "V", "test": "Charge quantisation: the global form contains Z6 (sealed as PR-S4)",
     "criterion": "t/3 + d/2 + Y integer on every Standard-Model line of the stack hypercharge assignment", "outcome": "pass",
     "result": "PASS on all eight lines; the rule holds for every odd T_3R coefficient, so it is necessary, not what fixes hypercharge (GF-10; "
               "sealed as PR-S4, re-sealed as PR-S5 the same day after review, CC-64)", "claims": ("GF-10", "CC-64"), "docs": ("GF01_FORCES_FIRST.md",)},
    {"id": "T8", "wave": "GF-01", "registered": "2026-10-11", "piece": "IX", "test": "The touch rule on the three block menus (sealed as PR-S3)",
     "criterion": "a Standard-Model-singlet block of the menu with nonzero dot product against the d^c or H_d line", "outcome": "mixed",
     "result": "U(8) FAIL (only the nu^c line, DS-1 reproduced); SO(16) and E8 PASS through sum lines and spinor sheets — in columns that do not "
               "give the exact Standard Model (GF-11); the non-supersymmetric U(8) was not scanned", "claims": ("GF-13",), "docs": ("GF01_FORCES_FIRST.md",)},
)

NEGATIVES = (
    {"id": "N1", "piece": "IX", "what": "Curve-only Yukawa textures are anarchic: the kinetic normalisation generates no hierarchy", "claims": (),
     "docs": ("CC26_NORMALISATION_RETRACTION.md",)},
    {"id": "N2", "piece": "I", "what": "No Hom-type unitary-stack flux model with Standard-Model bifundamental families is a consistent 6D gauge theory",
     "claims": (), "docs": ("CC29_R2C04_REGISTER.md",)},
    {"id": "N3", "piece": "VII", "what": "On a curve, chirality and Higgs mass are one number: no leading-order electroweak hierarchy", "claims": (),
     "docs": ("R2C06_TELESCOPING_REGISTER.md",)},
    {"id": "N4", "piece": "I", "what": "Reflection positivity of the lattice action: retracted from proven to open", "claims": (),
     "docs": ("CC35_REFLECTION_POSITIVITY.md",)},
    {"id": "N5", "piece": "VII", "what": "One-loop lifting of the Higgs direction: the computation failed", "claims": (), "docs": ("CW01_ONE_LOOP_POTENTIAL.md",)},
    {"id": "N6", "piece": "X", "what": "The Georgi–Jarlskog factor is out of reach of the curve-only twist class", "claims": (),
     "docs": ("SM12_TWIST_ASSIGNMENTS_REPORT.md",)},
    {"id": "N7", "piece": "VII", "what": "The ratio-3 background loses its coupling: an A-type Higgs carries no tree-level Yukawa", "claims": ("HP-3", "CC-43"),
     "docs": ()},
    {"id": "N8", "piece": "VII", "what": "Holomorphic messengers are forbidden at the family vertex", "claims": ("HJ-1", "CC-42"), "docs": ()},
    {"id": "N9", "piece": "VII", "what": "No jet can host the top", "claims": ("TF-1", "LB-4"), "docs": ()},
    {"id": "N10", "piece": "IV", "what": "No one-loop Casimir minimum in A_X for genus >= 2", "claims": ("MS-6",), "docs": ()},
    {"id": "N11", "piece": "IV", "what": "No canonical prescription for A_X/A_E hits the pre-registered window", "claims": ("CS-3",), "docs": ()},
    {"id": "N12", "piece": "IX", "what": "V_cb from the mass ratios misses the pre-registered window by a factor 1.9", "claims": ("TO-2",), "docs": ()},
    {"id": "N13", "piece": "VI", "what": "Flux designs quantise R; none reaches the data's 3.678", "claims": ("GC-3",), "docs": ()},
    {"id": "N14", "piece": "VI", "what": "b–tau needs a 30–35 per cent threshold at the allowed matching scales", "claims": ("GC-4",), "docs": ()},
    {"id": "N15", "piece": "X", "what": "The lepton pattern has no mechanism in the current matter sector", "claims": ("BL-5", "QL-1", "CI-4", "BL-2"),
     "docs": ()},
    {"id": "N16", "piece": "IX", "what": "A parity-symmetric down sector puts a floor under V_ub, 2.3 times the data with capped tails",
     "claims": ("DS-2", "DS-4"), "docs": ()},
    {"id": "N17", "piece": "IX", "what": "No odd source can reach the down sector: excluded by topology", "claims": ("OM-2",), "docs": ()},
    {"id": "N18", "piece": "VII", "what": "The bridged top fails y_t >= 0.5", "claims": ("YB-1",), "docs": ()},
    {"id": "N19", "piece": "IX", "what": "The two-point design does not carry the quark flavour structure", "claims": ("YB-2",), "docs": ()},
    {"id": "N20", "piece": "II", "what": "X0(143) is not the only supersingular semiprime level that folds onto a rank-1 elliptic curve: thirteen do",
     "claims": ("SQ-1", "SQ-2"), "docs": ()},
    {"id": "N22", "piece": "IV", "what": "The classical moduli flow has no finite critical point in the volume: the Einstein-frame terms all fall with Vol "
                                           "and the vacuum decompactifies; at one loop the Casimir energy has no minimum in A_X either", "claims": ("PG-4", "MS-6"), "docs": ()},
    {"id": "N21", "piece": "I", "what": "A supersymmetric SO(16) or E8 parent reached through U(8) with the stack hypercharge does not give the exact "
                                          "Standard Model: SO(16)'s quark doublets come in an even number, E8's three come with three chiral exotic doublets, "
                                          "and in both the lepton rows force net(L) + net(e^c) = 0", "claims": ("GF-11",), "docs": ()},
)

REFERENCE_DOCS = ("MATTER01_REPORT.md", "MATTER01_QED_DERIVATION.md", "INTERACTION_ATLAS.md")
AUDIT_TRAIL = ("REVIEW_V0311_RESPONSE.md", "RELEASE_NOTES_v0320_FOR_AUDITOR.md", "WAVE_B_ENGINES.md", "HANDOFF_RECONCILIATION.md",
               "OLOID_HODGE_RECONCILIATION.md", "RECONCILIATION_V0340_REVIEW_MATTER03.md", "MATTER_PROGRAMME_INTEGRATION_NOTES.md")
GENERATED_DOCS = ("MATTER_LEDGER.md", "MATTER_PUZZLE_MAP.md")

_ID = re.compile(r"\b(?:[A-Z]+[0-9]*-[0-9]+|M0[45]-[0-9])\b")


# ------------------------------------------------------------------------------------------------------------------------------------ lookups
def _ledger():
    from mtft.matter import ledger as L
    return L.E

def ordered():
    """The pieces in their numbering (dependency) order."""
    return sorted(PIECES, key=lambda p: ORDER.index(p.id))

def pieces():
    return {p.id: p for p in ordered()}

def subs():
    return {s.id: (p.id, s) for p in PIECES for s in p.subs}

def _wave_piece():
    return {w: p.id for p in PIECES for w in p.waves}

def _module_piece():
    """The piece that owns each module, through the waves of its non-correction claims."""
    wp, out = _wave_piece(), {}
    for e in _ledger():
        if e["wave"] != "correction" and e["wave"] in wp:
            out.setdefault(e["module"], set()).add(wp[e["wave"]])
    return out

def piece_of(entry):
    """The piece that owns a ledger entry: an explicit override, its wave, or (a correction) the piece of the module that filed it."""
    if entry["id"] in CLAIM_PIECE:
        return CLAIM_PIECE[entry["id"]]
    if entry["wave"] != "correction":
        return _wave_piece()[entry["wave"]]
    owners = _module_piece()[entry["module"]]
    if len(owners) != 1:
        raise ValueError(f"{entry['id']}: module {entry['module']} spans pieces {sorted(owners)}; add an override")
    return next(iter(owners))

def assignment():
    return {e["id"]: piece_of(e) for e in _ledger()}

def claims_of(pid):
    return [e for e in _ledger() if piece_of(e) == pid]

def modules_of(pid):
    return sorted({e["module"] for e in claims_of(pid)})

def layers():
    """Layer of each piece: 0 with no dependencies, otherwise one more than its deepest dependency."""
    P, out = pieces(), {}

    def layer(pid, seen=()):
        if pid in seen:
            raise ValueError(f"dependency cycle through {pid}")
        if pid not in out:
            out[pid] = 0 if not P[pid].depends_on else 1 + max(layer(d, seen + (pid,)) for d in P[pid].depends_on)
        return out[pid]

    for pid in P:
        layer(pid)
    return out

def chronology():
    from collections import Counter
    n = Counter(e["wave"] for e in _ledger())
    wp = _wave_piece()
    rows = [(WAVE_DATES[w], w, wp[w], n[w], WAVE_REPORT[w]) for w in WAVE_DATES]
    return sorted(rows, key=lambda r: (r[0], list(WAVE_DATES).index(r[1])))

def cited_ids():
    """Every ledger id cited by the map (verdicts, sub-puzzles, pre-registered tests, negatives, survivors), with where it is cited."""
    out = {}
    for p in PIECES:
        for c in p.key_claims:
            out.setdefault(c, []).append(f"{p.id} verdict")
        for c in p.survives:
            out.setdefault(c, []).append(f"{p.id} survives")
        for s in p.subs:
            for c in s.claims:
                out.setdefault(c, []).append(s.id)
    for t in PREREG:
        for c in t["claims"]:
            out.setdefault(c, []).append(t["id"])
    for n in NEGATIVES:
        for c in n["claims"]:
            out.setdefault(c, []).append(n["id"])
    return out

def all_docs():
    """Every report the map names, with its role."""
    out = {}
    for p in PIECES:
        for d in p.docs:
            out.setdefault(d, set()).add(f"{p.id} docs")
        for s in p.subs:
            for d in s.docs:
                out.setdefault(d, set()).add(s.id)
    for w, d in WAVE_REPORT.items():
        out.setdefault(d, set()).add(f"wave {w}")
    for t in PREREG:
        for d in t["docs"]:
            out.setdefault(d, set()).add(t["id"])
    for n in NEGATIVES:
        for d in n["docs"]:
            out.setdefault(d, set()).add(n["id"])
    for d in REFERENCE_DOCS:
        out.setdefault(d, set()).add("reference")
    for d in AUDIT_TRAIL:
        out.setdefault(d, set()).add("audit trail")
    for d in GENERATED_DOCS:
        out.setdefault(d, set()).add("generated")
    return out

def docs_dir():
    return Path(__file__).resolve().parents[3] / "docs" / "SM"


# ------------------------------------------------------------------------------------------------------------------------------------ checks
def check():
    """Integrity of the map against the ledger (and, in a source tree, the reports). Returns {check: bool}."""
    E = _ledger(); ids = {e["id"]: e for e in E}; P = pieces()
    out = {}
    try:
        a = assignment(); out["every_claim_assigned"] = len(a) == len(E) and all(v in P for v in a.values())
    except (KeyError, ValueError):
        out["every_claim_assigned"] = False
    waves = [w for p in PIECES for w in p.waves]
    out["waves_owned_once"] = len(waves) == len(set(waves)) and set(waves) == {e["wave"] for e in E if e["wave"] != "correction"}
    out["overrides_exist"] = set(CLAIM_PIECE) <= set(ids) and set(CLAIM_PIECE.values()) <= set(P)
    out["citations_exist"] = set(cited_ids()) <= set(ids)
    out["verdicts_cite_live_claims"] = all(ids[c]["status"] in ("active", "caveat") for p in PIECES for c in p.key_claims + p.survives)
    settled = [s for p in PIECES for s in p.subs if s.status in ("answered", "failed", "moot")]
    out["settled_subs_have_support"] = all(s.answer and (any(ids[c]["status"] in ("active", "caveat") for c in s.claims) or s.docs) for s in settled)
    out["open_subs_say_what_settles"] = all(s.settle for p in PIECES for s in p.subs if s.status == "open")
    out["statuses_valid"] = all(p.status in PIECE_STATUS for p in PIECES) and all(s.status in SUB_STATUS for p in PIECES for s in p.subs)
    out["numbered_in_order"] = [p.id for p in ordered()] == list(ORDER) and \
        all(ORDER.index(d) < ORDER.index(p.id) for p in PIECES for d in p.depends_on)
    out["sub_ids_well_formed"] = all(s.id.startswith(p.id + ".") for p in PIECES for s in p.subs) and len(subs()) == sum(len(p.subs) for p in PIECES)
    try:
        L = layers(); out["dependencies_acyclic"] = all(L[d] < L[p.id] for p in PIECES for d in p.depends_on)
    except (KeyError, ValueError):
        out["dependencies_acyclic"] = False
    out["status_coherent"] = all(_coherent(p) for p in PIECES)
    out["prereg_valid"] = all(t["outcome"] in OUTCOMES and t["piece"] in P for t in PREREG) and \
        all(s.prereg in {t["id"] for t in PREREG} for p in PIECES for s in p.subs if s.prereg)
    out["negatives_valid"] = all(n["piece"] in P and (n["claims"] or n["docs"]) for n in NEGATIVES)
    out["chronology_complete"] = set(WAVE_DATES) == set(WAVE_REPORT) == set(waves)
    d = docs_dir()
    if d.is_dir():
        on_disk = {f.name for f in d.glob("*.md")}
        named = set(all_docs())
        out["docs_exist"] = named <= on_disk
        out["every_doc_placed"] = on_disk <= named
    return out

def _coherent(p):
    st = [s.status for s in p.subs]
    if not st:
        return False
    if p.status == "stands":
        return "failed" not in st
    if p.status == "closed":
        return all(s.status != "open" or s.reopen for s in p.subs) and any(s in ("failed", "answered") for s in st)
    if p.status in ("partial", "open"):
        return "open" in st
    return False


# ------------------------------------------------------------------------------------------------------------------------------------ views
def summary():
    from collections import Counter
    a = assignment(); E = {e["id"]: e for e in _ledger()}
    rows = {}
    for p in ordered():
        c = Counter(E[i]["status"] for i, q in a.items() if q == p.id)
        s = Counter(x.status for x in p.subs)
        rows[p.id] = {"name": p.name, "status": p.status, "claims": sum(c.values()), "claim_status": dict(c), "subs": dict(s),
                      "waves": len(p.waves), "modules": len(modules_of(p.id)), "docs": len(p.docs)}
    return rows

def to_json():
    E = {e["id"]: e for e in _ledger()}; L = layers(); a = assignment()
    out = {"pieces": [], "prereg": [dict(t, claims=list(t["claims"]), docs=list(t["docs"])) for t in PREREG],
           "negatives": [dict(n, claims=list(n["claims"]), docs=list(n["docs"])) for n in NEGATIVES],
           "retracted": [{"id": e["id"], "piece": a[e["id"]], "statement": e["statement"], "note": e["note"]} for e in E.values() if e["status"] == "retracted"],
           "chronology": [{"date": d, "wave": w, "piece": p, "claims": n, "report": r} for d, w, p, n, r in chronology()],
           "ledger_size": len(E)}
    for p in ordered():
        d = asdict(p)
        d["layer"] = L[p.id]
        d["modules"] = modules_of(p.id)
        d["claims"] = [{"id": e["id"], "wave": e["wave"], "class": e["class"], "status": e["status"]} for e in claims_of(p.id)]
        d["key_claim_text"] = {c: E[c]["statement"] for c in p.key_claims}
        out["pieces"].append(d)
    return out

def to_dot():
    colour = {"stands": "#2e7d32", "partial": "#f9a825", "open": "#1565c0", "closed": "#757575"}
    lines = ["digraph puzzle {", "  rankdir=LR;", '  node [shape=box, style="rounded,filled", fontname="Helvetica", fontcolor=white];']
    for p in PIECES:
        lines.append(f'  {p.id} [label="{p.id} {p.name}\\n{p.status}", fillcolor="{colour[p.status]}"];')
    for p in PIECES:
        for d in p.depends_on:
            lines.append(f"  {d} -> {p.id};")
    lines.append("}")
    return "\n".join(lines) + "\n"

def to_markdown():
    E = {e["id"]: e for e in _ledger()}; L = layers(); S = summary()
    n_open = sum(1 for p in PIECES for s in p.subs if s.status == "open")
    out = ["# The matter programme — puzzle map", "",
           "Generated by `python -m mtft matter puzzle --markdown` from `mtft.matter.puzzle`; do not edit by hand. Each piece is one question the "
           "programme asks; it owns the ledger claims that address it and carries its own sub-puzzles. Statuses: stands, partial, open, closed "
           "(an honest negative); sub-puzzles are open, answered, failed or moot.", "",
           f"{len(PIECES)} pieces, {sum(len(p.subs) for p in PIECES)} sub-puzzles ({n_open} open), {len(E)} ledger claims, {len(PREREG)} pre-registered "
           f"tests, {len(NEGATIVES)} recorded negatives.", "",
           "## The pieces in dependency order", "",
           "| Layer | Piece | Status | Claims | Open sub-puzzles | Depends on |", "|---|---|---|---|---|---|"]
    for p in sorted(ordered(), key=lambda p: (L[p.id], ORDER.index(p.id))):
        n = sum(1 for s in p.subs if s.status == "open")
        out.append(f"| {L[p.id]} | {p.id} {p.name} | {p.status} | {S[p.id]['claims']} | {n} of {len(p.subs)} | {', '.join(_deps(p)) or '—'} |")
    for p in sorted(ordered(), key=lambda p: (L[p.id], ORDER.index(p.id))):
        out += ["", f"## {p.id} — {p.name} ({p.status})", "", f"**Question.** {p.question}", "", f"**Verdict.** {p.verdict}", ""]
        if p.waves:
            out.append(f"**Waves.** {', '.join(f'{w} ({WAVE_REPORT[w]})' for w in p.waves)}.")
            out.append(f"**Modules.** {', '.join('`' + m + '`' for m in modules_of(p.id))}.")
            c = S[p.id]["claim_status"]
            out.append(f"**Claims.** {S[p.id]['claims']} ({', '.join(f'{v} {k}' for k, v in sorted(c.items()))}).")
        if p.docs:
            out.append(f"**Reports.** {', '.join(p.docs)}.")
        if p.survives:
            out.append(f"**What survives.** {', '.join(p.survives)}.")
        out += ["", "| Sub-puzzle | Status | Answer, or what would settle it | Claims |", "|---|---|---|---|"]
        for s in p.subs:
            text = s.answer if s.status != "open" else s.settle
            tag = " (pre-registered " + s.prereg + ")" if s.prereg else ""
            tag += " — would reopen the piece" if s.reopen else ""
            refs = ", ".join(list(s.claims) + list(s.docs))
            out.append(f"| {s.id} {s.question} | {s.status}{tag} | {text} | {refs} |")
    out += ["", "## Pre-registered tests", "", "| Test | Wave | Registered | Criterion | Outcome | Result |", "|---|---|---|---|---|---|"]
    for t in PREREG:
        out.append(f"| {t['id']} {t['test']} | {t['wave']} | {t['registered']} | {t['criterion']} | {t['outcome']} | {t['result']} |")
    out += ["", "## Honest negatives", "", "| | Piece | What was closed | Record |", "|---|---|---|---|"]
    for n in NEGATIVES:
        out.append(f"| {n['id']} | {n['piece']} | {n['what']} | {', '.join(list(n['claims']) + list(n['docs']))} |")
    retracted = [e for e in E.values() if e["status"] == "retracted"]
    out += ["", f"Retracted claims ({len(retracted)}): " + "; ".join(f"{e['id']} ({e['note'].split(':')[0]})" for e in retracted) + ".", ""]
    out += ["## Chronology of the waves", "", "| Date | Wave | Piece | Claims | Report |", "|---|---|---|---|---|"]
    for d, w, p, n, r in chronology():
        out.append(f"| {d} | {w} | {p} | {n} | {r} |")
    out += ["", "## Other reports", "", f"Reference: {', '.join(REFERENCE_DOCS)}.", f"Audit trail: {', '.join(AUDIT_TRAIL)}.",
            f"Generated: {', '.join(GENERATED_DOCS)}.", ""]
    return "\n".join(out)

def _deps(p):
    return sorted(p.depends_on, key=ORDER.index)

def text(pid=None, only_open=False):
    """Plain-text view for the command line."""
    L = layers(); S = summary(); lines = []
    for p in sorted(ordered(), key=lambda p: (L[p.id], ORDER.index(p.id))):
        if pid and p.id != pid:
            continue
        lines.append(f"{p.id} [{p.status}] {p.name} — layer {L[p.id]}, {S[p.id]['claims']} claims; depends on {', '.join(_deps(p)) or 'nothing'}")
        if pid:
            lines.append(f"  Q: {p.question}")
            lines.append(f"  Verdict: {p.verdict}")
        for s in p.subs:
            if only_open and s.status != "open":
                continue
            lines.append(f"    {s.id} [{s.status}] {s.question}" + (f" — {s.settle if s.status == 'open' else s.answer}" if pid or only_open else ""))
    return "\n".join(lines)

def json_text():
    return json.dumps(to_json(), indent=2, ensure_ascii=False)
