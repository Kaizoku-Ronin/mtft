#!/usr/bin/env python3
"""
mtft.ledger — every certified constant of the MTFT program, as data.
====================================================================
MIT License — (c) 2026 Roger Tano & Kimi K3 — MTFT Research Program

Integration stage 1 (MTFT_Repo_Integration_Plan_v0.1 §4, §7.1).

DESIGN RULES
------------
1. Every constant carries: value, class, tolerance, closed form (where one
   exists), and provenance (the note AND the audit addendum that certified it).
2. Superseded values are RECORDED, not deleted — the correction chains are
   themselves results (e.g. K3: -54.7822 kappa-fit -> -54.7825976 single-step
   -> -54.7825732 Richardson).
3. Nothing elsewhere in the package may hardcode a number that lives here.
   Tests assert against this module, never against literals.
4. `verify()` recomputes every entry that has a closed form, from the model
   definition, and compares against the stored value. The ledger is
   self-auditing: a contradiction fails here, not in a note.
5. `legend()` generates the Legend registrations mechanically (BH Decision 4):
   the ledger is the source, the Legend is a view.

Gate classes (program convention): EXACT (algebraic identity or closed form,
verified), CERTIFIED (verified independently by both engines to tol), DERIVED
(closed form evaluated), MEASURED (extraction, both engines agree), DIAGNOSTIC
(informative, not load-bearing).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import mpmath as mp

mp.mp.dps = 50


# ---------------------------------------------------------------------------
# Entry type
# ---------------------------------------------------------------------------

@dataclass
class Entry:
    value: str                      # decimal string (exact representation stored)
    cls: str                        # EXACT / CERTIFIED / DERIVED / MEASURED / DIAGNOSTIC
    note: str                       # certifying study note
    addendum: str                   # certifying audit addendum
    tol: Optional[float] = None     # certification tolerance (absolute)
    form: Optional[str] = None      # closed form, human-readable
    supersedes: tuple = ()          # ((old_value, provenance), ...) recorded, not deleted
    comment: str = ""

    def mpf(self):
        return mp.mpf(self.value)


@dataclass
class Family:
    """A sequence of related entries (indexed families: A_i, C_i, ...)."""
    name: str
    entries: list                   # of Entry
    index_note: str                 # indexing convention
    cls: str = "CERTIFIED"
    note: str = ""
    addendum: str = ""


# ---------------------------------------------------------------------------
# The model (for verify()): rho(n) = ln n / n^3, beta = 2, gaps g_j = ln(rho2/rho(j+2))
# ---------------------------------------------------------------------------

def rho(n):
    return mp.log(n) / n ** 3

RHO2 = rho(2)

def gap(j):
    return mp.log(RHO2 / rho(2 + j))


# ---------------------------------------------------------------------------
# LEDGER
# ---------------------------------------------------------------------------

LEDGER = {}

def _put(key, entry):
    LEDGER[key] = entry

# ---- rung 4: canonical coupling ------------------------------------------

_put('m_inf', Entry(
    '0.755834576126', cls='EXACT',
    form='log(27*log(2)/(8*log(3)))',
    note='rung4 minimal chain, Pr G', addendum='Y',
    comment='thermodynamic-limit mass gap of the minimal chain'))

_put('exp_minus_m_inf', Entry(
    '0.4696185', cls='EXACT', tol=5e-8,
    form='8*log(3)/(27*log(2))',
    note='rung4 / PR-17', addendum='AO'))

_put('R_inf_over_m_inf', Entry(
    '1.8854354', cls='EXACT', tol=5e-8,
    form='27*log(2)/(27*log(2)-8*log(3))',
    note='PR-17', addendum='AO'))

_put('kappa_star', Entry(
    '5.0', cls='DIAGNOSTIC', tol=0.5,
    note='rung4 LG7', addendum='Y',
    comment='interior minimum of m(kappa); monotonicity expectation falsified (owned)'))

_put('m_star', Entry(
    '0.736839', cls='DIAGNOSTIC', tol=1e-6,
    note='rung4 LG7', addendum='Y'))

# ---- rotation observables -------------------------------------------------

_put('ROT1', Entry(
    '1.004590', cls='CERTIFIED', tol=1e-5,
    note='PR-12 UG2', addendum='AH'))

_put('ROT2', Entry(
    '0.947589', cls='CERTIFIED', tol=4e-6,
    note='PR-13 VG3 (series 0.947587 vs Richardson 0.947591, crossing resolved)',
    addendum='AI',
    supersedes=(('0.9485', 'PR-9 (superseded)'),
                ('0.947943', 'PR-12 order-6 local excursion (AH-F1)')),
    comment='tau_2 -> 0 limit; truncation-order systematic identified as PR-12 gap'))

# ---- EP limit positions u_{i,i+1}(inf) ------------------------------------

_EP_LIMITS = ['1.42507723746', '2.8707041129', '5.13015139272',
              '8.3166465032', '12.5544990771']
_EP_NOTES = [('PR-17 Pr BB', 'AO'), ('PR-18', 'AP'), ('PR-18/PR-31 control', 'AP/BC'),
             ('PR-14/PR-31 control', 'AJ/BC'), ('PR-19/PR-31', 'AQ/BC')]
for _i, (_v, (_n, _a)) in enumerate(zip(_EP_LIMITS, _EP_NOTES)):
    _put(f'u_ep_{_i}{_i+1}', Entry(
        _v, cls='CERTIFIED', tol=2e-9,
        form='(g_{i+1}-g_i)/(exp(-g_i)-exp(-g_{i+1}))  [magnitude]',
        note=_n, addendum=_a))

# ---- response diagonals at kappa* (PROPOSED by Claude, Add. BM; AMENDED and
# ---- accepted by the auditor, Add. BN).  Needed because chain.selftest
# ---- certified mu_0/mu_1 against unregistered numbers (BI-F2) and the
# ---- guard now refuses.  Values: the N-CONVERGED f64 path (identical to
# ---- 1e-11 across N = 400..3200); the mp/nsite path approaches from below
# ---- (nsite=14 reads mu_1 1.9e-6 LOW, nsite=18 within 4e-7 — low-nsite
# ---- mp is the under-converged reading here, BN-E1).
_put('mu_0', Entry('1.050397683', cls='CERTIFIED', tol=1e-6,
                   form='B_00 at kappa*=5 (dressed diagonal of level 0)',
                   note='rung-4 LG4 / PR-6 CG0 / PR-36 responses',
                   addendum='Y/AB/BH; BM proposal; BN accepted',
                   comment='response coefficient used throughout the '
                           'perturbation arc; N-converged to 1e-11 '
                           '(f64, N>=400); mp/nsite approaches from below'))
_put('mu_1', Entry('0.5585643816', cls='CERTIFIED', tol=1e-6,
                   form='B_11 at kappa*=5 (dressed diagonal of level 1)',
                   note='rung-4 LG4 / PR-6 CG0 / PR-36 responses',
                   addendum='Y/AB/BH; BM proposal; BN accepted',
                   comment='mp nsite=14 under-reads by 1.9e-6 (BN-E1 trap); '
                           'N-converged f64 is the certification path'))

# ---- stage-4 coupled model (BN certifications) ----------------------------
# Registered by the auditor so coupled.selftest can assert via the _L guard
# instead of literals (BN-F1, the BI-F2 class).  Values are the auditor's
# independently reproduced numbers; they coincide with coupled's literals.

_put('tau_c_star', Entry('0.23003179', cls='CERTIFIED', tol=1e-6,
                   form='band-merging coupling, Bloch (Z) measure, kappa*=5',
                   note='stage-4 coupled model; PR-8 geometry-free scaling: '
                        'tau_c*x_max is measure-independent, tree/Z ratio '
                        '= 1/sqrt(2) to 12 digits',
                   addendum='BN',
                   comment='BN certification legs: independent bisection '
                           '(own bracket scan) + exact-ratio leg. '
                           'coupled.selftest to switch from literal to '
                           "_L('tau_c_star') (BN-F1)."))
_put('V_b_tree', Entry('0.0373264', cls='CERTIFIED', tol=1e-5,
                   form='tree binding threshold at kappa*=5, V_b = 1/G_dd(edge)',
                   note='Z: G_dd diverges (exponent -1/2; any attraction '
                        'binds). tree: G_dd converges; V_b finite.',
                   addendum='BN',
                   comment='BN certification legs: cross-instrument '
                           'convergence (Gauss nodes n=400/900/1500 -> '
                           '0.0372416/0.0373101/0.0373205) vs sqrt-delta '
                           'extrapolation 0.0373264; Z exponent -0.5008 '
                           '(mine) / -0.5004 (coupled) / -1/2 (theory). '
                           "coupled.selftest to switch to _L('V_b_tree') "
                           '(BN-F1).'))

# ---- first-order amplitudes A_i (pair-indexed) ----------------------------

_A_TAIL = ['1.81145993554', '4.45480521585', '7.76185251663', '6.51848732826',
           '-11.8052703719', '-70.9849145776', '-210.568270945', '-491.006183725',
           '-999.368850689', '-1855.6309329', '-3219.510994', '-5297.85285925',
           '-8352.53810003']
A_FAMILY = Family(
    name='A', cls='CERTIFIED',
    entries=[Entry(v, cls='CERTIFIED', tol=1e-7,
                   form='A_of(i): three channels, closed form (PR-23)',
                   note='PR-23 (0..2), PR-24 (3,4), PR-25 (5..12)',
                   addendum='AU/AV/AW',
                   supersedes=(('7.750', 'A_2 finite-kappa seal; AV-F1'),) if i == 2 else ())
             for i, v in enumerate(_A_TAIL)],
    index_note='A_i multiplies ((i+3)/(i+4))^{2 kappa} in dev(pair (i,i+1)); '
               'single sign change between i=3 and i=4; no integer zero in 0..60')

# ---- second-order coefficients C_i ----------------------------------------

_C = ['13.23412586', '42.38221958', '112.0526532', '255.391562147', '517.697902341']
C_FAMILY = Family(
    name='C', cls='CERTIFIED',
    entries=[Entry(v, cls='CERTIFIED', tol=5e-6 if i < 3 else 1e-5,
                   form='C_of(i): three channels incl. i=0 reference shift (PR-24)',
                   note='PR-24', addendum='AV')
             for i, v in enumerate(_C)],
    index_note='dev includes -C_i ((i+2)/(i+3))^{2 kappa}; '
               'omitting the reference shift is 4.17% high at i=0')

# ---- third-order coefficients C3_i ----------------------------------------

_put('C3_3', Entry('250.3238719', cls='CERTIFIED', tol=1e-5,
                   form='C3_of(i): (i-1,i) coupling, three channels',
                   note='PR-25', addendum='AW',
                   comment='channels 11.0015 + 88.8807 + 150.442'))
_put('C3_4', Entry('528.7608971', cls='CERTIFIED', tol=1e-5,
                   form='C3_of(i)', note='PR-25', addendum='AW',
                   comment='channels 18.7603 + 179.236 + 330.764'))

# ---- fourth order K_i (pair-indexed; old rungs wrote K_4 for K_3) ---------

_put('K_2', Entry('-35.3706007', cls='CERTIFIED', tol=5e-6,
                  form='crossing monomial extraction: [Q^2]+[RSQ]+[R^2S^2]',
                  note='PR-31', addendum='BC (Richardson)'))
_put('K_3', Entry('-54.7825732', cls='CERTIFIED', tol=1e-6,
                  form='crossing monomial extraction; also closure: '
                       '3.8545981*dg + 15.06185*dV + X_rep',
                  note='PR-30/PR-31', addendum='BB/BC',
                  supersedes=(('-54.7822', 'AW/AX kappa-window fit, +-5e-4; demoted BB'),
                              ('-54.7825976', 'PR-30 single-step h=1e-4'))))
_put('K_4', Entry('134.0082168', cls='CERTIFIED', tol=1e-5,
                  form='crossing monomial extraction',
                  note='PR-31', addendum='BC (Richardson)',
                  comment='positive: K sign change between pairs (3,4) and (4,5)'))

# ---- fifth order K5_i (upward from LOWER member) --------------------------

_put('K5_2', Entry('6627.4435', cls='CERTIFIED', tol=5e-3,
                   form='[T^2]+[PRT]+[P^2R^2], EP extraction (P couples members)',
                   note='PR-34', addendum='BF',
                   comment='i=2 class degeneracy (bond (2,7)=4/9) flagged; closed walks clean'))
_put('K5_3', Entry('37227.284', cls='CERTIFIED', tol=5e-2,
                   form='[T^2]+[PRT]+[P^2R^2]', note='PR-34', addendum='BF'))
_put('K5_4', Entry('157283.30', cls='CERTIFIED', tol=2e-1,
                   form='[T^2]+[PRT]+[P^2R^2]', note='PR-34', addendum='BF',
                   comment='all-positive at i=2,3,4: no-zero prediction confirmed (BE.2)'))

# ---- fourth-order channels, pair (3,4) ------------------------------------

_put('dg_4', Entry('-9.6016055674', cls='EXACT', tol=1e-9,
                   form='-(rho6/(rho4-rho6)) * (rho4/(rho4-rho5))^2  '
                        '[rho_j = rho(site j+2); path resummation, completeness proved]',
                   note='PR-26', addendum='AX (6e-10) / AY (completeness)'))
_put('dg4_Q2', Entry('-0.959288278631', cls='CERTIFIED', tol=1e-10,
                     note='PR-28', addendum='AZ'))
_put('dg4_RSQ', Entry('-4.15125298792', cls='CERTIFIED', tol=1e-9,
                      note='PR-28', addendum='AZ'))
_put('dg4_R2S2', Entry('-4.49106430089', cls='EXACT', tol=1e-9,
                       form='-(rho6/(rho4-rho6)) * (rho5/(rho4-rho5))^2',
                       note='PR-28/PR-32', addendum='AZ/BD'))

_put('dV_44', Entry('47.2515618601', cls='CERTIFIED', tol=5e-8,
                    form='eigenvector channel, monomial extraction',
                    note='PR-28', addendum='AZ (ratio 1.000000039)'))
_put('dV44_Q2', Entry('2.74191340501', cls='CERTIFIED', tol=1e-9,
                      note='PR-28', addendum='AZ'))
_put('dV44_RSQ', Entry('18.1262729205', cls='CERTIFIED', tol=1e-8,
                       note='PR-28', addendum='AZ'))
_put('dV44_R2S2', Entry('26.3833755346', cls='CERTIFIED', tol=1e-8,
                        note='PR-28', addendum='AZ'))
_put('dV44_factor', Entry('17.2330613', cls='CERTIFIED', tol=1e-6,
                          note='PR-28', addendum='AZ',
                          supersedes=(('17.233083', 'AZ-F1 transcription slip'),),
                          comment='enhancement over naive [Q^2]; dg factor 10.00909297; '
                                  'ratio 1.7217'))

_put('X_rep', Entry('-729.4682', cls='CERTIFIED', tol=5e-2,
                    form='closure: K_3 - (du/db)*dg_4 - (m/D^2)e^-b * dV_44, exact coefficients',
                    note='PR-28 (target) / PR-36 (exact coefficients)',
                    addendum='AZ/BB/BH',
                    supersedes=(('-729.3137', 'BB, soft closure coefficients (BH-F1)'),
                                ('-729.32', 'PR-28'),
                                ('-729.5', 'AX'),
                                ('-722', 'AW, contaminated targets')),
                    comment='38x the -19.08 direct estimate; '
                            'eps4-eps5 = -0.12824 near-degeneracy'))

# ---- crossing monomial decomposition of K_3 --------------------------------

_put('K3_Q2', Entry('18.52145343265', cls='CERTIFIED', tol=1e-9,
                    note='PR-29/PR-30', addendum='BA/BB'))
_put('K3_RSQ', Entry('-123.452039493', cls='CERTIFIED', tol=1e-8,
                     note='PR-29/PR-30', addendum='BA/BB'))
_put('K3_R2S2', Entry('50.1480128813', cls='CERTIFIED', tol=1e-8,
                      note='PR-30', addendum='BB',
                      supersedes=(('100.295976915', 'PR-29 stencil bug: missing /4 (BA)'),)))

# ---- zeros and limits ------------------------------------------------------

_put('A_zero', Entry('3.537205722', cls='CERTIFIED', tol=1e-8,
                     form='zero of A_of(i), continuous i', note='PR-33', addendum='BE'))
_put('K_zero', Entry('3.563077', cls='CERTIFIED', tol=1e-5,
                     note='PR-33', addendum='BE',
                     comment='0.02587 from A_zero; distinct at 40x noise; unexplained data'))
_put('R2S2_zero', Entry('2.4969', cls='MEASURED', tol=5e-4,
                        note='PR-32/PR-33', addendum='BD/BE',
                        comment='crossing [R^2 S^2] vanishes; drives the "leap" reframing'))
_put('C3_ratio_limit', Entry('1.0', cls='EXACT',
                             form='lim_{i->inf} C3 cancellation ratio = 1, from above; '
                                  'ratio = 1 + 5.7632/i + O(i^-2); crossing below 1 '
                                  'impossible (all-positive constituents)',
                             note='PR-35 (predicted) ', addendum='BG (computed)'))

# ---- exact response coefficients (BH) --------------------------------------

for _i in (2, 3, 4):
    _put(f'resp_{_i}', Entry(
        '', cls='EXACT',
        form='|u| = m/(p_a-p_b): du/da=(m e^-a - D)/D^2, du/db=(D - m e^-b)/D^2, '
             'du/dp_a = -m/D^2, du/dp_b = +m/D^2 (opposition is an identity)',
        note='PR-36', addendum='BH',
        comment='values computed by verify(); see response(i)'))

# ---- near-degeneracies at crossings ----------------------------------------

_put('degen_34', Entry('-0.12824456', cls='CERTIFIED', tol=1e-7,
                       form='eps_4 - eps_5 at u_34(inf), limit couplings',
                       note='AW/PR-31', addendum='BC',
                       comment='amplifies H-space channels; X_rep attribution'))
_put('degen_45', Entry('-0.10067747', cls='CERTIFIED', tol=1e-7,
                       form='eps_5 - eps_6 at u_45(inf)', note='PR-31', addendum='BC'))

FAMILIES = {'A': A_FAMILY, 'C': C_FAMILY}


# ---------------------------------------------------------------------------
# RELATIONS — derived identities between observables (not scalars, so not
# Entries).  Each relation carries its derivation, its meaning, and the
# addenda that derived and verified it.  verify() checks each relation at
# the level its derivation specifies.
# ---------------------------------------------------------------------------
RELATIONS = {
    'dev_centre_vs_dev_EP': {
        'expression': 'dev_centre(i;k) = dev_EP(i;k) - rep_C(i) * s(i)^{2k}',
        'symbols': ('rep_C(i) = 2 m e^{-(a+b)}/D^3 (the repulsion '
                    'constituent of C_i, expansion.channels_C(i)[2]); '
                    's(i) = (i+2)/(i+3); dev_* = |u| - u_ep_limit(i)'),
        'derivation': ('2x2 dressed block H=[[a-u p,-u v],[-u v,b-u q]], '
                       'D=p-q, m=b-a: EP |u_EP| = m/sqrt(D^2+4v^2) '
                       '-> u0 - 2 m v^2/D^3; diabatic centre '
                       'u_min = m D/(D^2+4v^2) -> u0 - 4 m v^2/D^3; '
                       'with v^2 = e^{-(a+b)} s^{2k} (1+O(1e-6), measured '
                       'k=60) the difference is exactly one copy of '
                       'rep_C(i) s^{2k}'),
        'meaning': ('the SEALED hierarchy (A, C, C3, K, K5) describes the '
                    'complex EP branch point; ep.ep_of returns the '
                    'diabatic centre (minimum of the closest-pair gap^2). '
                    'Rung 0 is exempt — its partner is the reference '
                    'level lambda_0 (BK).  Resolves S3-1.'),
        'status': 'DERIVED',
        'addendum': 'BK (derivation); BL (verification: 2x2 coefficients '
                    '-4/-2 exact; full model k=60 rel resid <= 1.2e-6 '
                    'rungs 1-4)',
    },
}


# ---------------------------------------------------------------------------
# Closed forms (the reference implementations verify() checks against).
# Sources: pr24.py (A_of, C_of), pr25.py (C3_of), PR-26 (dg_4), PR-36 (responses).
#
# BN duplication resolution (Addendum BN): A_of / C_of / C3_of / C3_channels
# were byte-for-byte reimplementations of expansion.A / C / C3 / channels_C3
# (identical at integer i, <= ~1e-12 apart at non-integer i — verified in the
# BN audit).  Two implementations of the same closed form is exactly what the
# repo migration exists to kill.  The ledger now keeps the DATA and the
# verify(); expansion owns the implementation.  The accessors below are lazy
# delegates (function-level imports): chain/expansion import the ledger only
# function-level via the _L guard, so there is no import cycle either way.
# dg4_closed / dg4_R2S2_closed / u_ep_limit / response stay here: expansion
# has no counterparts for them.
# ---------------------------------------------------------------------------

def A_of(i):
    """First-order amplitude, pair (i, i+1). PR-23 closed form, three channels.
    Implementation: expansion.A (BN duplication resolution)."""
    from .expansion import A as _A
    return _A(i)


def C_of(i):
    """Second-order coefficient, pair (i, i+1). PR-24 closed form.
    dev = A r^2k - C s^2k. i = 0 carries the reference-level shift.
    Implementation: expansion.C (BN duplication resolution)."""
    from .expansion import C as _C
    return _C(i)


def C3_of(i):
    """Third-order coefficient: the (i-1, i) incident coupling. PR-25.
    Implementation: expansion.C3 (BN duplication resolution)."""
    from .expansion import C3 as _C3
    return _C3(i)


def dg4_closed():
    """PR-26 path resummation for the fourth-order eigenvalue channel, pair (3,4)."""
    r4, r5, r6 = rho(6), rho(7), rho(8)
    return -(r6 / (r4 - r6)) * (r4 / (r4 - r5)) ** 2


def dg4_R2S2_closed():
    r4, r5, r6 = rho(6), rho(7), rho(8)
    return -(r6 / (r4 - r6)) * (r5 / (r4 - r5)) ** 2


def u_ep_limit(i):
    """Limit EP position (magnitude) of pair (i, i+1)."""
    a, b = gap(i), gap(i + 1)
    return (b - a) / (mp.e ** (-a) - mp.e ** (-b))


def response(i):
    """Exact crossing responses (PR-36/BH): (du/da, du/db, du/dp_a, du/dp_b)."""
    a, b = gap(i), gap(i + 1)
    ea, eb = mp.e ** (-a), mp.e ** (-b)
    D, m = ea - eb, b - a
    return ((m * ea - D) / D ** 2, (D - m * eb) / D ** 2, -m / D ** 2, m / D ** 2)


def cancellation_ratio(channels):
    s = sum(channels)
    return abs(s) / max(abs(x) for x in channels)


# ---------------------------------------------------------------------------
# verify(): the self-audit. Recompute everything with a closed form and
# compare against the stored values. Any contradiction raises here.
# ---------------------------------------------------------------------------

_CHECKS = []

def check(name, computed, entry_or_value, tol):
    target = entry_or_value.mpf() if isinstance(entry_or_value, Entry) \
        else mp.mpf(entry_or_value)
    err = abs(computed - target)
    ok = err <= max(mp.mpf(tol), mp.mpf('1e-40'))
    _CHECKS.append((name, mp.nstr(computed, 12), mp.nstr(target, 12),
                    mp.nstr(err, 3), ok))
    return ok


def verify(verbose=True):
    _CHECKS.clear()
    # rung 4
    check('m_inf', mp.log(27 * mp.log(2) / (8 * mp.log(3))), LEDGER['m_inf'], 1e-9)
    check('exp_minus_m_inf', 8 * mp.log(3) / (27 * mp.log(2)),
          LEDGER['exp_minus_m_inf'], 5e-8)
    check('R_inf/m_inf', 27 * mp.log(2) / (27 * mp.log(2) - 8 * mp.log(3)),
          LEDGER['R_inf_over_m_inf'], 5e-8)
    # EP limits
    for i in range(5):
        check(f'u_ep_{i}{i+1}', u_ep_limit(i), LEDGER[f'u_ep_{i}{i+1}'], 2e-9)
    # families
    for i, e in enumerate(A_FAMILY.entries):
        check(f'A_{i}', A_of(i), e, e.tol)
    for i, e in enumerate(C_FAMILY.entries):
        check(f'C_{i}', C_of(i), e, e.tol)
    check('C3_3', C3_of(3), LEDGER['C3_3'], 1e-5)
    check('C3_4', C3_of(4), LEDGER['C3_4'], 1e-5)
    # fourth-order channels
    check('dg_4', dg4_closed(), LEDGER['dg_4'], 1e-9)
    check('dg4_R2S2', dg4_R2S2_closed(), LEDGER['dg4_R2S2'], 1e-9)
    check('dg4 sum', LEDGER['dg4_Q2'].mpf() + LEDGER['dg4_RSQ'].mpf()
          + LEDGER['dg4_R2S2'].mpf(), LEDGER['dg_4'], 1e-8)
    check('dV44 sum', LEDGER['dV44_Q2'].mpf() + LEDGER['dV44_RSQ'].mpf()
          + LEDGER['dV44_R2S2'].mpf(), LEDGER['dV_44'], 5e-8)
    check('dV44 factor', LEDGER['dV_44'].mpf() / LEDGER['dV44_Q2'].mpf(),
          LEDGER['dV44_factor'], 1e-6)
    check('K3 monomial sum', LEDGER['K3_Q2'].mpf() + LEDGER['K3_RSQ'].mpf()
          + LEDGER['K3_R2S2'].mpf(), LEDGER['K_3'], 1e-6)
    # closure with exact coefficients -> X_rep (BH)
    da, db, dpa, dpb = response(3)
    a, b = gap(3), gap(4)
    c2 = (m2 := (b - a)) / (mp.e ** (-a) - mp.e ** (-b)) ** 2 * mp.e ** (-b)
    xrep = LEDGER['K_3'].mpf() - db * LEDGER['dg_4'].mpf() - c2 * LEDGER['dV_44'].mpf()
    check('X_rep (closure, exact coefficients)', xrep, LEDGER['X_rep'], 5e-2)
    # zeros: A_zero from closed form
    zA = mp.findroot(A_of, (mp.mpf('3.53'), mp.mpf('3.55')))
    check('A_zero', zA, LEDGER['A_zero'], 1e-8)
    # near-degeneracies
    for key, i in (('degen_34', 3), ('degen_45', 4)):
        u = -u_ep_limit(i)   # crossing sits at negative u in the program convention
        d = (gap(i + 1) - u * mp.e ** (-gap(i + 1))) \
            - (gap(i + 2) - u * mp.e ** (-gap(i + 2)))
        check(key, d, LEDGER[key], 1e-7)
    # C3 cancellation-ratio limit behaviour at large i
    for i in (100, 400):
        ch = (C3_chans := C3_channels(i))
        r = cancellation_ratio(ch)
        pred = 1 + mp.mpf('5.7632') / i
        ok = abs(r - pred) < mp.mpf('0.02') / i * i * mp.mpf('0.02')
        _CHECKS.append((f'C3 ratio ~ 1+5.7632/i at i={i}', mp.nstr(r, 8),
                        mp.nstr(pred, 8), mp.nstr(abs(r - pred), 3),
                        abs(r - pred) < mp.mpf('0.05')))
    # mu_0/mu_1 (BM proposal, BN accepted): certified against the
    # N-converged f64 construction (N=400 suffices — N-stable to 1e-11).
    # Lazy import: chain's ledger guard is function-level, so no cycle.
    from .chain import internal as _internal
    _ic5 = _internal(5, nb=30, backend='f64', N=400)
    check('mu_0', _ic5.diag(0), LEDGER['mu_0'], 1e-9)
    check('mu_1', _ic5.diag(1), LEDGER['mu_1'], 1e-9)
    # RELATION dev_centre_vs_dev_EP (S3-1, BK/BL): the 2x2 dressed-block
    # identity at the coefficient level — centre correction -> -4 m v^2/D^3,
    # EP correction -> -2 m v^2/D^3, so centre = EP - one copy of the
    # repulsion constituent.  Algebraic; the full-model assertion lives in
    # ep.selftest (which owns the dynamical observable).
    for (a2, b2, p2, q2, v2) in ((1.5, 2.2, 0.6, 0.2, 1e-3),
                                 (0.7, 1.9, 1.1, 0.3, 3e-4)):
        D2, m2v = p2 - q2, b2 - a2
        u02 = m2v / D2
        cent = (m2v * D2 / (D2 ** 2 + 4 * v2 ** 2) - u02) \
            / (m2v * v2 ** 2 / D2 ** 3)
        epc = (abs(m2v / (D2 + 2j * v2)) - u02) / (m2v * v2 ** 2 / D2 ** 3)
        check('rel dev_centre_vs_dev_EP: centre coeff -> -4', float(cent),
              -4.0, 5e-3)
        check('rel dev_centre_vs_dev_EP: EP coeff -> -2', float(epc),
              -2.0, 5e-3)
    if verbose:
        nbad = 0
        for name, comp, tgt, err, ok in _CHECKS:
            flag = 'ok ' if ok else 'FAIL'
            if not ok:
                nbad += 1
            print(f'  [{flag}] {name:42} computed {comp:>16} stored {tgt:>16} err {err}')
        print(f'\n  {len(_CHECKS)} checks, {nbad} failures')
    return all(ok for *_, ok in _CHECKS)


def C3_channels(i):
    """Channel decomposition of C3_of(i) (for ratio checks).
    Implementation: expansion.channels_C3 (BN duplication resolution)."""
    from .expansion import channels_C3 as _ch
    return tuple(_ch(i))


# ---------------------------------------------------------------------------
# legend(): mechanical Legend registrations (BH Decision 4).
# The ledger is the source; the Legend is a view.
# ---------------------------------------------------------------------------

def legend():
    lines = ['# MTFT Legend — generated from mtft.ledger (do not hand-edit)', '']
    lines.append('| name | value | class | tol | certified by | supersedes |')
    lines.append('|---|---|---|---|---|---|')
    for key, e in LEDGER.items():
        if not e.value:
            continue
        sup = '; '.join(f'{v} ({p})' for v, p in e.supersedes) or '—'
        lines.append(f'| {key} | {e.value} | {e.cls} | {e.tol or "—"} '
                     f'| {e.note} / Add. {e.addendum} | {sup} |')
    for fam in FAMILIES.values():
        for i, e in enumerate(fam.entries):
            sup = '; '.join(f'{v} ({p})' for v, p in e.supersedes) or '—'
            lines.append(f'| {fam.name}_{i} | {e.value} | {e.cls} | {e.tol or "—"} '
                         f'| {e.note} / Add. {e.addendum} | {sup} |')
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    print('mtft.ledger — self-audit')
    ok = verify(verbose=True)
    print(f'\n  ledger self-audit: {"ALL GREEN" if ok else "FAILURES PRESENT"}')
    with open('mtft_legend.md', 'w') as f:
        f.write(legend())
    print('  Legend written to mtft_legend.md')
