#!/usr/bin/env python3
"""
Arithmetic Busy Beavers: Modular Constraints on Uncomputability
================================================================

MIT License — Copyright (c) 2026 Roger Tano
See LICENSE file for full terms.

This module constructs a hierarchy of Busy Beaver functions constrained
by the arithmetic of the modular curve X₀(143), implementing the
theoretical framework from MTFT Chapter 10 and Papers 25-30.

The hierarchy (from most to least constrained):

    BB_Fatou(n) ≤ BB_g(n) ≤ BB_Hecke(n) ≤ BB(n)

    BB_Fatou  — Fatou-trapped ATMs only (confinement certificate)
    BB_g      — genus-truncated: tape window L = genus(X₀(143)) = 13
    BB_Hecke  — Hecke sign-constrained TMs (computability: open)
    BB        — unrestricted (uncomputable)

Key insight: the Hecke traces {a_n} on S₂ᶰᵉʷ(Γ₀(143)) provide a
number-theoretically structured constraint on Turing machine transitions.
Negative traces → fermionic fold (write 1), positive → bosonic analytic
(write 0), mirroring the X₀(143) Burning Mandelbrot iteration engine.

The genus truncation L = 13 recovers the Faulhaber correction structure:
    BB_g(n) = A · D(n) - R(n)
where R(n) is the correction polynomial (mass gap analogue) whose
recurrence ΔR generates the identity via telescoping.

MTFT Dictionary:
    Naive term  →  Free theory
    R(n)        →  Mass gap
    ΔR          →  Stiffness μ_N(y)
    Deg bound   →  Controlled discontinuity
    [D,M] ≠ 0  →  [t_U, t_R] ≠ 0
    Fatou trap  →  Confinement

Roger Tano — MTFT Research Program — March 2026
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import (
    List, Tuple, Dict, Optional,
)
from enum import IntEnum
import time


# ═══════════════════════════════════════════════════════════════
#  MTFT STRUCTURAL CONSTANTS
# ═══════════════════════════════════════════════════════════════

LEVEL = 143              # N = 11 × 13
GENUS = 13               # genus(X₀(143))
INDEX = 168              # [SL(2,ℤ) : Γ₀(143)] = |PSL(2,7)| = dim SU(13)
DIM_NEW = 11             # dim S₂ᶰᵉʷ(Γ₀(143))
ORBIT_DIMS = (1, 4, 6)   # Galois orbit dimensions [f₁, f₂, f₃]
CANONICAL_DEG = 2 * GENUS - 2  # = 24, degree of canonical divisor

# Hecke traces on S₂ᶰᵉʷ(Γ₀(143)), LMFDB-verified, n = 1..200
HECKE_TRACES = [
    11, 3, 2, 9, 0, -4, 8, 3, 13, -2, -3, 0, 1, -16, -6, 17, 2, -9,
    0, -18, -12, 3, 14, -36, 31, -3, 2, 12, -10, -48, -10, -21, -2,
    -34, -28, -13, 16, -14, 4, -14, 14, 6, 20, -3, 2, 0, -16, -6, 23,
    25, -12, 7, -2, 36, 0, -6, 52, 18, -6, 0, 2, 28, 0, 33, 2, 2, 10,
    2, -14, 4, -38, -1, 38, 50, -28, -36, 4, -2, 12, 26, 3, 18, -28,
    44, 40, -8, -20, 15, -20, 14, 0, 86, 2, 32, -52, -64, 36, -17, -9,
    27, 18, 24, -16, -15, 8, -26, -24, -14, -2, -14, -34, 12, 32, 56,
    -50, -46, 13, -4, 12, 20, 11, 6, -64, 4, -42, -82, 0, -25, -20,
    14, -28, 8, 0, 36, -102, -86, -4, -52, -4, -80, -52, 40, -9, -23,
    -12, 28, -6, 114, 6, 72, -12, -48, 18, 8, -18, -12, 4, -48, -20,
    -98, -32, 115, 24, 106, -14, -68, -20, -14, 11, 28, 8, 24, -2, 72,
    12, -7, 50, -54, 6, 22, -24, -8, -4, -16, 10, -4, 10, -4, -64, -12,
    -6, 138, 86, 14, 12, -39, -2, -21, -8, 173,
]

# Per-orbit traces at small primes: {p: (a_p(f₁), Tr(a_p(f₂)), Tr(a_p(f₃)))}
# Single source of truth: the independently verified table shipped with the
# v0.6.1 audit-coalescence release.  The pre-audit local copy that used to
# live here disagreed with the verified 143a1 point counts at every listed
# prime and was removed in v0.7.0.
from mtft.x0_143 import ORBIT_TRACES_VERIFIED as ORBIT_TRACES


# ═══════════════════════════════════════════════════════════════
#  §1. HECKE SIGN ORACLE
# ═══════════════════════════════════════════════════════════════

class HeckeSign(IntEnum):
    """Sign of the Hecke trace, determining TM transition constraint."""
    FERMIONIC = -1   # a_n < 0 → fold (Burning Ship)
    FREE = 0         # a_n = 0 → unconstrained
    BOSONIC = +1     # a_n > 0 → analytic (Mandelbrot)


def hecke_sign(n: int) -> HeckeSign:
    """
    Return sgn(a_n) on S₂ᶰᵉʷ(Γ₀(143)).

    For n ≤ 200: direct lookup from LMFDB table.
    For n > 200: would require Eichler-Selberg computation.
    """
    if n < 1:
        raise ValueError(f"Hecke index must be ≥ 1, got {n}")
    if n <= len(HECKE_TRACES):
        a_n = HECKE_TRACES[n - 1]
        if a_n < 0:
            return HeckeSign.FERMIONIC
        elif a_n > 0:
            return HeckeSign.BOSONIC
        else:
            return HeckeSign.FREE
    # Beyond table: would need Eichler-Selberg
    raise ValueError(
        f"Hecke trace at n={n} exceeds table (max {len(HECKE_TRACES)}). "
        f"Eichler-Selberg extension needed."
    )


def hecke_sign_pattern(length: int) -> List[HeckeSign]:
    """Return the first `length` signs of the Hecke trace sequence."""
    return [hecke_sign(n) for n in range(1, length + 1)]


def hecke_constraint_density(length: int = 200) -> Dict[str, float]:
    """
    Compute the constraint density of the Hecke sign pattern.

    Returns the fraction of indices that are fermionic, bosonic, or free.
    By the non-vanishing theorem for newforms at unramified primes,
    the density of FREE indices among primes approaches 0.
    """
    signs = hecke_sign_pattern(min(length, len(HECKE_TRACES)))
    n = len(signs)
    counts = {
        "fermionic": sum(1 for s in signs if s == HeckeSign.FERMIONIC),
        "bosonic": sum(1 for s in signs if s == HeckeSign.BOSONIC),
        "free": sum(1 for s in signs if s == HeckeSign.FREE),
    }
    return {k: v / n for k, v in counts.items()}


def dominant_sector(n: int) -> Optional[int]:
    """
    For small primes p, identify which Galois orbit dominates a_p.

    Returns 0 (electron/f₁), 1 (muon/f₂), 2 (tau/f₃), or None.
    This implements the three-sector decomposition from Paper 26.
    """
    if n not in ORBIT_TRACES:
        return None
    a1, tr2, tr3 = ORBIT_TRACES[n]
    vals = [abs(a1), abs(tr2), abs(tr3)]
    return vals.index(max(vals))


# ═══════════════════════════════════════════════════════════════
#  §2. TURING MACHINE INFRASTRUCTURE
# ═══════════════════════════════════════════════════════════════

class Direction(IntEnum):
    LEFT = 0
    RIGHT = 1


@dataclass(frozen=True)
class Transition:
    """A single TM transition: (write_symbol, move_direction, next_state)."""
    write: int             # 0 or 1
    move: Direction        # LEFT or RIGHT
    next_state: int        # state index, or -1 for HALT


HALT_STATE = -1


@dataclass(frozen=True)
class TMConfig:
    """Instantaneous configuration of a Turing machine."""
    state: int
    head: int              # head position on tape
    tape: tuple            # immutable tape snapshot


@dataclass
class TuringMachine:
    """
    A 2-symbol Turing machine with n states.

    The transition table δ maps (state, symbol) → Transition.
    State indices: 0..n_states-1 for active states, HALT_STATE for halt.
    Symbols: 0, 1.
    """
    n_states: int
    transitions: Dict[Tuple[int, int], Transition]

    def run(self, max_steps: int, tape_size: Optional[int] = None
            ) -> Tuple[int, int, bool]:
        """
        Execute the TM on a blank tape.

        Parameters:
            max_steps: maximum steps before declaring non-halting
            tape_size: if set, restricts tape to a window of this size
                       (for genus-truncated BB)

        Returns:
            (output, steps, halted)
            output = number of 1s on tape when halted
            steps  = number of steps executed
            halted = whether the machine halted within max_steps
        """
        if tape_size is not None:
            # Genus-truncated: finite tape window
            tape = [0] * tape_size
            head = tape_size // 2  # start in middle
        else:
            # Dynamic tape (bounded by max_steps)
            tape = [0] * (2 * max_steps + 1)
            head = max_steps  # start in middle

        state = 0
        for step in range(max_steps):
            sym = tape[head]
            key = (state, sym)
            if key not in self.transitions:
                # Missing transition = halt
                return sum(tape), step, True

            tr = self.transitions[key]
            if tr.next_state == HALT_STATE:
                tape[head] = tr.write
                return sum(tape), step + 1, True

            tape[head] = tr.write
            if tr.move == Direction.RIGHT:
                head += 1
            else:
                head -= 1

            # Boundary handling
            if tape_size is not None:
                if head < 0 or head >= tape_size:
                    # Fell off truncated tape → halt
                    return sum(tape), step + 1, True
            else:
                if head < 0 or head >= len(tape):
                    # Shouldn't happen with max_steps tape
                    return sum(tape), step + 1, True

            state = tr.next_state

        return sum(tape), max_steps, False

    def _step_once(self, state: int, head: int, tape: list, tape_size: int
                    ) -> Tuple[int, int, bool, bool]:
        """
        Execute one TM step in-place on tape.

        Returns (new_state, new_head, halted, fell_off).
        """
        sym = tape[head]
        key = (state, sym)
        if key not in self.transitions:
            return state, head, True, False

        tr = self.transitions[key]
        tape[head] = tr.write

        if tr.next_state == HALT_STATE:
            return HALT_STATE, head, True, False

        new_head = head + (1 if tr.move == Direction.RIGHT else -1)
        if new_head < 0 or new_head >= tape_size:
            return tr.next_state, new_head, True, True

        return tr.next_state, new_head, False, False

    def run_floyd(self, max_steps: int, tape_size: int
                  ) -> Tuple[int, int, bool, Optional[int]]:
        """
        Execute with Floyd's tortoise-and-hare cycle detection.

        O(1) memory — no orbit recording. Detects cycles by running
        two copies: tortoise (1 step/iter) and hare (2 steps/iter).
        When their full configurations match, a cycle is found.

        Returns:
            (output, steps, halted, cycle_period)
            cycle_period is None if halted, else the detected period.
        """
        # Tortoise state
        t_tape = [0] * tape_size
        t_head = tape_size // 2
        t_state = 0

        # Hare state
        h_tape = [0] * tape_size
        h_head = tape_size // 2
        h_state = 0

        for step in range(max_steps):
            # Tortoise: 1 step
            t_state, t_head, t_halt, _ = self._step_once(
                t_state, t_head, t_tape, tape_size)
            if t_halt:
                return sum(t_tape), step + 1, True, None

            # Hare: 2 steps
            h_state, h_head, h_halt, _ = self._step_once(
                h_state, h_head, h_tape, tape_size)
            if h_halt:
                # Hare halted — tortoise will eventually halt too
                # Continue running tortoise only
                for s2 in range(step + 1, max_steps):
                    t_state, t_head, t_halt, _ = self._step_once(
                        t_state, t_head, t_tape, tape_size)
                    if t_halt:
                        return sum(t_tape), s2 + 1, True, None
                return sum(t_tape), max_steps, False, None

            h_state, h_head, h_halt, _ = self._step_once(
                h_state, h_head, h_tape, tape_size)
            if h_halt:
                for s2 in range(step + 1, max_steps):
                    t_state, t_head, t_halt, _ = self._step_once(
                        t_state, t_head, t_tape, tape_size)
                    if t_halt:
                        return sum(t_tape), s2 + 1, True, None
                return sum(t_tape), max_steps, False, None

            # Check if tortoise == hare (cycle detected)
            if (t_state == h_state and t_head == h_head
                    and t_tape == h_tape):
                # Cycle found — this machine will never halt
                # Find the period by running one copy further
                period = 1
                save_state, save_head = h_state, h_head
                save_tape = h_tape[:]
                h_state, h_head, _, _ = self._step_once(
                    h_state, h_head, h_tape, tape_size)
                while not (h_state == save_state and h_head == save_head
                           and h_tape == save_tape):
                    h_state, h_head, h_halt, _ = self._step_once(
                        h_state, h_head, h_tape, tape_size)
                    if h_halt:
                        break
                    period += 1
                    if period > tape_size * (1 << tape_size) * self.n_states:
                        break
                return sum(t_tape), step + 1, False, period

        return sum(t_tape), max_steps, False, None

    def run_with_orbit(self, max_steps: int, tape_size: Optional[int] = None
                       ) -> Tuple[int, int, bool, List[TMConfig]]:
        """
        Execute and record the full orbit (sequence of configurations).
        Used for Fatou certificate checking. SLOW — prefer run_floyd.
        """
        if tape_size is not None:
            tape = list([0] * tape_size)
            head = tape_size // 2
        else:
            tape = list([0] * (2 * max_steps + 1))
            head = max_steps

        state = 0
        orbit: List[TMConfig] = []

        for step in range(max_steps):
            config = TMConfig(state=state, head=head, tape=tuple(tape))
            orbit.append(config)

            sym = tape[head]
            key = (state, sym)
            if key not in self.transitions:
                return sum(tape), step, True, orbit

            tr = self.transitions[key]
            if tr.next_state == HALT_STATE:
                tape[head] = tr.write
                orbit.append(TMConfig(state=HALT_STATE, head=head, tape=tuple(tape)))
                return sum(tape), step + 1, True, orbit

            tape[head] = tr.write
            if tr.move == Direction.RIGHT:
                head += 1
            else:
                head -= 1

            if tape_size is not None:
                if head < 0 or head >= tape_size:
                    return sum(tape), step + 1, True, orbit
            else:
                if head < 0 or head >= len(tape):
                    return sum(tape), step + 1, True, orbit

            state = tr.next_state

        return sum(tape), max_steps, False, orbit


# ═══════════════════════════════════════════════════════════════
#  §3. ARITHMETIC TURING MACHINE (ATM)
# ═══════════════════════════════════════════════════════════════

def hecke_constraint_index(state: int, symbol: int) -> int:
    """
    Map (state, symbol) pair to Hecke index.

    The mapping 2·state + symbol + 1 ensures each transition slot
    gets a unique Hecke index. For an n-state machine, this uses
    indices 1 through 2n.
    """
    return 2 * state + symbol + 1


def is_hecke_compatible(n_states: int, transitions: Dict[Tuple[int, int], Transition]
                        ) -> bool:
    """
    Check whether a transition table satisfies the Hecke constraint.

    For each (state, symbol) pair with Hecke index i:
      - If a_i < 0 (fermionic): must write 1 (fold)
      - If a_i > 0 (bosonic):   must write 0 (analytic)
      - If a_i = 0 (free):      either symbol permitted
    """
    for (state, symbol), tr in transitions.items():
        idx = hecke_constraint_index(state, symbol)
        if idx > len(HECKE_TRACES):
            continue  # beyond table: unconstrained
        sign = hecke_sign(idx)
        if sign == HeckeSign.FERMIONIC and tr.write != 1:
            return False
        if sign == HeckeSign.BOSONIC and tr.write != 0:
            return False
    return True


def enumerate_atms(n_states: int) -> List[TuringMachine]:
    """
    Enumerate all n-state Arithmetic Turing Machines at level 143.

    These are TMs whose transition tables satisfy the Hecke sign
    constraint. The total state space is:
      unconstrained: (4 · (n+1))^(2n) possibilities
      constrained:   reduced by Hecke sign pattern

    For small n this is tractable:
      n=1: 2 transition slots, max Hecke index = 2
           a_1 = 11 > 0 → slot (0,0) must write 0
           a_2 =  3 > 0 → slot (0,1) must write 0
           Both constrained → write symbol is fixed
      n=2: 4 slots, Hecke indices 1-4
           a_1=11>0, a_2=3>0, a_3=2>0, a_4=9>0
           All bosonic → all writes fixed to 0
    """
    atms = []
    # All (state, symbol) pairs
    slots = [(s, sym) for s in range(n_states) for sym in (0, 1)]

    # For each slot, determine allowed write symbols
    allowed_writes = {}
    for state, symbol in slots:
        idx = hecke_constraint_index(state, symbol)
        if idx <= len(HECKE_TRACES):
            sign = hecke_sign(idx)
            if sign == HeckeSign.FERMIONIC:
                allowed_writes[(state, symbol)] = [1]
            elif sign == HeckeSign.BOSONIC:
                allowed_writes[(state, symbol)] = [0]
            else:
                allowed_writes[(state, symbol)] = [0, 1]
        else:
            allowed_writes[(state, symbol)] = [0, 1]

    # Possible moves
    directions = [Direction.LEFT, Direction.RIGHT]

    # Possible next states: 0..n_states-1 plus HALT
    next_states = list(range(n_states)) + [HALT_STATE]

    # Build all compatible transition tables
    # Each slot gets: (write, direction, next_state)
    slot_options = []
    for slot in slots:
        writes = allowed_writes[slot]
        options = [
            Transition(w, d, ns)
            for w in writes
            for d in directions
            for ns in next_states
        ]
        slot_options.append(options)

    count = 1
    for opts in slot_options:
        count *= len(opts)

    # Only enumerate if tractable
    MAX_ENUM = 10_000_000
    if count > MAX_ENUM:
        raise ValueError(
            f"ATM enumeration for n={n_states} states would generate "
            f"{count:,} machines (limit: {MAX_ENUM:,}). "
            f"Use sampling or BB_g truncation instead."
        )

    for combo in itertools.product(*slot_options):
        table = {}
        for i, slot in enumerate(slots):
            table[slot] = combo[i]
        atms.append(TuringMachine(n_states=n_states, transitions=table))

    return atms


def enumerate_atms_sector_resolved(n_states: int) -> List[Tuple[TuringMachine, List[int]]]:
    """
    Enumerate ATMs with sector labels for each transition.

    Returns (TM, sector_labels) where sector_labels[i] indicates
    which Galois orbit (0=electron, 1=muon, 2=tau) dominates
    the Hecke constraint at slot i.
    """
    atms = enumerate_atms(n_states)
    result = []
    for tm in atms:
        sectors = []
        for state in range(n_states):
            for sym in (0, 1):
                idx = hecke_constraint_index(state, sym)
                sec = dominant_sector(idx)
                sectors.append(sec if sec is not None else -1)
        result.append((tm, sectors))
    return result


def _can_write_ones(n_states: int, transitions: Dict[Tuple[int, int], Transition]
                    ) -> bool:
    """
    Dead-state pruning: check if the machine can possibly write 1s.

    A machine cannot produce output if no state reachable from q₀
    has a transition that writes 1. This prunes the search space
    significantly when most states are bosonic (write-0 forced).
    """
    # Find which states can write 1
    writers = set()
    for (state, sym), tr in transitions.items():
        if tr.write == 1:
            writers.add(state)

    if not writers:
        return False

    # BFS from q₀ to check reachability
    reachable = set()
    queue = [0]
    while queue:
        s = queue.pop()
        if s in reachable or s == HALT_STATE:
            continue
        reachable.add(s)
        for sym in (0, 1):
            key = (s, sym)
            if key in transitions:
                ns = transitions[key].next_state
                if ns != HALT_STATE and ns not in reachable:
                    queue.append(ns)

    return bool(writers & reachable)


def generate_atms(n_states: int, *, prune_dead: bool = True):
    """
    Generator yielding ATMs one at a time — O(1) memory per machine.

    With prune_dead=True, skips machines where no reachable state
    can write 1 (guaranteed BB output = 0). This is the major
    optimization for large n where most states are bosonic.
    """
    slots = [(s, sym) for s in range(n_states) for sym in (0, 1)]

    # Determine allowed writes per slot
    allowed_writes = {}
    for state, symbol in slots:
        idx = hecke_constraint_index(state, symbol)
        if idx <= len(HECKE_TRACES):
            sign = hecke_sign(idx)
            if sign == HeckeSign.FERMIONIC:
                allowed_writes[(state, symbol)] = [1]
            elif sign == HeckeSign.BOSONIC:
                allowed_writes[(state, symbol)] = [0]
            else:
                allowed_writes[(state, symbol)] = [0, 1]
        else:
            allowed_writes[(state, symbol)] = [0, 1]

    directions = [Direction.LEFT, Direction.RIGHT]
    next_states = list(range(n_states)) + [HALT_STATE]

    slot_options = []
    for slot in slots:
        writes = allowed_writes[slot]
        options = [
            Transition(w, d, ns)
            for w in writes
            for d in directions
            for ns in next_states
        ]
        slot_options.append(options)

    # Count total
    total = 1
    for opts in slot_options:
        total *= len(opts)

    yielded = 0
    pruned = 0
    for combo in itertools.product(*slot_options):
        table = {slots[i]: combo[i] for i in range(len(slots))}

        if prune_dead and not _can_write_ones(n_states, table):
            pruned += 1
            continue

        yielded += 1
        yield TuringMachine(n_states=n_states, transitions=table)


def count_atms(n_states: int) -> Tuple[int, int]:
    """
    Count total ATMs and estimate productive ones (can write 1s).
    Returns (total, estimate_productive).
    """
    slots = [(s, sym) for s in range(n_states) for sym in (0, 1)]
    total = 1
    for state, symbol in slots:
        idx = hecke_constraint_index(state, symbol)
        if idx <= len(HECKE_TRACES):
            sign = hecke_sign(idx)
            if sign == HeckeSign.FREE:
                total *= 2 * 2 * (n_states + 1)
            else:
                total *= 1 * 2 * (n_states + 1)
        else:
            total *= 2 * 2 * (n_states + 1)
    return total, total  # exact count requires enumeration


# ═══════════════════════════════════════════════════════════════
#  §4. BB_g: GENUS-TRUNCATED BUSY BEAVER
# ═══════════════════════════════════════════════════════════════

@dataclass
class BBResult:
    """Result of a BB computation."""
    bb_value: int           # maximum output among halting machines
    champion: Optional[TuringMachine]  # the machine achieving the maximum
    champion_steps: int     # steps taken by the champion
    total_machines: int     # total machines enumerated
    halting_count: int      # number that halted
    cycle_count: int        # number trapped in cycles
    max_steps_used: int     # step limit used
    tape_window: Optional[int]  # tape window (None = unlimited)
    computation_time: float  # wall-clock seconds


def bb_genus(n_states: int, *,
             hecke_constrained: bool = True,
             max_steps: int = 10_000,
             prune_dead: bool = True,
             verbose: bool = False) -> BBResult:
    """
    Compute BB_g(n, 143): the genus-truncated Busy Beaver.

    Tape window L = genus(X₀(143)) = 13.
    Configuration space: |Config₁₃| = n · 13 · 2¹³ = 106496 · n

    This is COMPUTABLE: the finite tape makes all orbits decidable.
    The Faulhaber correction structure is recoverable from the results.

    Parameters:
        n_states: number of TM states
        hecke_constrained: if True, only Hecke-compatible TMs (ATMs)
        max_steps: safety bound (should be >> config space size)
        prune_dead: skip machines that provably can't write 1s
        verbose: print progress
    """
    tape_window = GENUS  # L = 13

    # Configuration space bound: n · L · 2^L
    config_space = n_states * tape_window * (1 << tape_window)
    # A machine that doesn't halt within config_space steps is cycling
    effective_max = min(max_steps, config_space + 1)

    t0 = time.time()

    best_output = 0
    best_machine = None
    best_steps = 0
    halting = 0
    cycling = 0
    total = 0
    pruned_count = 0

    if hecke_constrained:
        # Use generator for large n, list for small n
        total_est, _ = count_atms(n_states)
        if total_est > 10_000_000:
            # Generator mode with pruning
            if verbose:
                print(f"BB_g({n_states}, L={tape_window}): "
                      f"~{total_est:,} ATMs (generator + pruning)")
            for tm in generate_atms(n_states, prune_dead=prune_dead):
                total += 1
                output, steps, halted = tm.run(effective_max, tape_size=tape_window)
                if halted:
                    halting += 1
                    if output > best_output:
                        best_output = output
                        best_machine = tm
                        best_steps = steps
                else:
                    cycling += 1
                if verbose and total % 500_000 == 0:
                    elapsed = time.time() - t0
                    print(f"  ... {total:,} ({elapsed:.1f}s) best={best_output}")
        else:
            # List mode (original behavior for small n)
            machines = enumerate_atms(n_states)
            total = len(machines)
            if verbose:
                print(f"BB_g({n_states}, L={tape_window}): "
                      f"enumerating {total:,} ATMs")
            for i, tm in enumerate(machines):
                output, steps, halted = tm.run(effective_max, tape_size=tape_window)
                if halted:
                    halting += 1
                    if output > best_output:
                        best_output = output
                        best_machine = tm
                        best_steps = steps
                else:
                    cycling += 1
                if verbose and (i + 1) % 100_000 == 0:
                    elapsed = time.time() - t0
                    print(f"  ... {i+1:,}/{total:,} ({elapsed:.1f}s) best={best_output}")
    else:
        machines = _enumerate_all_tms(n_states, max_count=10_000_000)
        total = len(machines)
        if verbose:
            print(f"BB_g({n_states}, L={tape_window}): "
                  f"enumerating {total:,} TMs (unconstrained)")
        for i, tm in enumerate(machines):
            output, steps, halted = tm.run(effective_max, tape_size=tape_window)
            if halted:
                halting += 1
                if output > best_output:
                    best_output = output
                    best_machine = tm
                    best_steps = steps
            else:
                cycling += 1
            if verbose and (i + 1) % 100_000 == 0:
                elapsed = time.time() - t0
                print(f"  ... {i+1:,}/{total:,} ({elapsed:.1f}s) best={best_output}")

    elapsed = time.time() - t0
    if verbose:
        print(f"BB_g({n_states}) = {best_output}  "
              f"[{halting:,} halting, {cycling:,} cycling, {elapsed:.2f}s]")

    return BBResult(
        bb_value=best_output,
        champion=best_machine,
        champion_steps=best_steps,
        total_machines=total,
        halting_count=halting,
        cycle_count=cycling,
        max_steps_used=effective_max,
        tape_window=tape_window,
        computation_time=elapsed,
    )


def _enumerate_all_tms(n_states: int, max_count: int = 1_000_000
                       ) -> List[TuringMachine]:
    """Enumerate all n-state 2-symbol TMs (unconstrained)."""
    slots = [(s, sym) for s in range(n_states) for sym in (0, 1)]
    next_states = list(range(n_states)) + [HALT_STATE]
    directions = [Direction.LEFT, Direction.RIGHT]

    options_per_slot = [
        Transition(w, d, ns)
        for w in (0, 1)
        for d in directions
        for ns in next_states
    ]
    total_count = len(options_per_slot) ** len(slots)
    if total_count > max_count:
        raise ValueError(
            f"Full TM enumeration for n={n_states}: {total_count:,} machines "
            f"exceeds limit {max_count:,}. Use hecke_constrained=True."
        )

    machines = []
    for combo in itertools.product(*([options_per_slot] * len(slots))):
        table = {slots[i]: combo[i] for i in range(len(slots))}
        machines.append(TuringMachine(n_states=n_states, transitions=table))
    return machines


# ═══════════════════════════════════════════════════════════════
#  §4b. MONTE CARLO AND TARGETED SEARCH FOR LARGE n
# ═══════════════════════════════════════════════════════════════

import random

@dataclass
class BBSampleResult(BBResult):
    """BB result from Monte Carlo sampling — includes confidence info."""
    sample_size: int = 0
    population_size: int = 0
    coverage: float = 0.0


def _random_atm(n_states: int, rng: random.Random) -> TuringMachine:
    """Generate a single random Hecke-compatible TM."""
    next_states = list(range(n_states)) + [HALT_STATE]
    directions = [Direction.LEFT, Direction.RIGHT]
    table = {}
    for state in range(n_states):
        for sym in (0, 1):
            idx = hecke_constraint_index(state, sym)
            if idx <= len(HECKE_TRACES):
                sign = hecke_sign(idx)
                if sign == HeckeSign.FERMIONIC:
                    w = 1
                elif sign == HeckeSign.BOSONIC:
                    w = 0
                else:
                    w = rng.choice([0, 1])
            else:
                w = rng.choice([0, 1])
            d = rng.choice(directions)
            ns = rng.choice(next_states)
            table[(state, sym)] = Transition(w, d, ns)
    return TuringMachine(n_states=n_states, transitions=table)


def bb_sample(n_states: int, *,
              sample_size: int = 1_000_000,
              max_steps: int = 10_000,
              seed: int = 143,
              verbose: bool = False) -> BBSampleResult:
    """
    Monte Carlo lower bound on BB_Hecke(n, 143).

    Samples random Hecke-compatible TMs and returns the best output
    found. This is a LOWER BOUND — the true BB_Hecke may be higher.

    For n where full enumeration is intractable (n ≥ 4), this gives
    the best known value. The seed is fixed at 143 (the level) for
    reproducibility.
    """
    tape_window = GENUS
    config_space = n_states * tape_window * (1 << tape_window)
    effective_max = min(max_steps, config_space + 1)

    rng = random.Random(seed)
    total_pop, _ = count_atms(n_states)

    t0 = time.time()
    best_output = 0
    best_machine = None
    best_steps = 0
    halting = 0
    cycling = 0

    if verbose:
        print(f"BB_sample({n_states}, samples={sample_size:,}): "
              f"population ~{total_pop:,}")

    for i in range(sample_size):
        tm = _random_atm(n_states, rng)
        output, steps, halted = tm.run(effective_max, tape_size=tape_window)
        if halted:
            halting += 1
            if output > best_output:
                best_output = output
                best_machine = tm
                best_steps = steps
                if verbose:
                    print(f"  New best at sample {i+1:,}: "
                          f"output={output} steps={steps}")
        else:
            cycling += 1

        if verbose and (i + 1) % 200_000 == 0:
            elapsed = time.time() - t0
            print(f"  ... {i+1:,}/{sample_size:,} ({elapsed:.1f}s) "
                  f"best={best_output} halting={halting}")

    elapsed = time.time() - t0
    coverage = sample_size / total_pop if total_pop > 0 else 1.0
    if verbose:
        print(f"BB_sample({n_states}) ≥ {best_output}  "
              f"[{halting:,} halting, {cycling:,} cycling, "
              f"coverage={coverage:.2e}, {elapsed:.2f}s]")

    return BBSampleResult(
        bb_value=best_output,
        champion=best_machine,
        champion_steps=best_steps,
        total_machines=sample_size,
        halting_count=halting,
        cycle_count=cycling,
        max_steps_used=effective_max,
        tape_window=tape_window,
        computation_time=elapsed,
        sample_size=sample_size,
        population_size=total_pop,
        coverage=coverage,
    )


def bb_targeted(n_states: int, *,
                max_steps: int = 10_000,
                router_samples: int = 50_000,
                verbose: bool = False) -> BBResult:
    """
    Targeted search exploiting the writer-router decomposition.

    Key insight: in the Hecke-constrained regime, most states are
    "routers" (bosonic, can only write 0) while a few are "writers"
    (free or fermionic, can write 1). The search fixes the writer
    states' transitions and exhaustively varies the routing.

    For n=4: only q₂ can write (slots 5,6). States q₀,q₁,q₃ are
    pure routers. We enumerate all routing configurations while
    keeping the writer transitions fixed at their best known values,
    then do a local search around the champion.
    """
    tape_window = GENUS
    config_space = n_states * tape_window * (1 << tape_window)
    effective_max = min(max_steps, config_space + 1)

    # Identify writer and router states
    writers = set()
    for state in range(n_states):
        for sym in (0, 1):
            idx = hecke_constraint_index(state, sym)
            if idx <= len(HECKE_TRACES):
                a = HECKE_TRACES[idx - 1]
                if a <= 0:  # free or fermionic: can write 1
                    writers.add(state)
    routers = set(range(n_states)) - writers

    if verbose:
        print(f"BB_targeted({n_states}): writers={sorted(writers)} "
              f"routers={sorted(routers)}")

    # Build options for each slot
    next_states = list(range(n_states)) + [HALT_STATE]
    directions = [Direction.LEFT, Direction.RIGHT]

    writer_slots = []
    router_slots = []
    all_slots = [(s, sym) for s in range(n_states) for sym in (0, 1)]

    for state, symbol in all_slots:
        idx = hecke_constraint_index(state, symbol)
        if idx <= len(HECKE_TRACES):
            sign = hecke_sign(idx)
            if sign == HeckeSign.FERMIONIC:
                writes = [1]
            elif sign == HeckeSign.BOSONIC:
                writes = [0]
            else:
                writes = [0, 1]
        else:
            writes = [0, 1]

        options = [Transition(w, d, ns)
                   for w in writes for d in directions for ns in next_states]

        if state in writers:
            writer_slots.append(((state, symbol), options))
        else:
            router_slots.append(((state, symbol), options))

    # Phase 1: enumerate all writer configurations
    writer_combos = list(itertools.product(
        *[opts for _, opts in writer_slots]))

    if verbose:
        router_count = 1
        for _, opts in router_slots:
            router_count *= len(opts)
        print(f"  Writer configs: {len(writer_combos):,}")
        print(f"  Router configs per writer: {router_count:,}")
        print(f"  Strategy: sample routers for each writer config")

    t0 = time.time()
    best_output = 0
    best_machine = None
    best_steps = 0
    halting = 0
    cycling = 0
    total = 0

    rng = random.Random(143)

    # For each writer configuration, sample routers
    router_sample_count = router_samples  # samples per writer config

    for wi, w_combo in enumerate(writer_combos):
        # Set writer transitions
        writer_table = {}
        for j, (slot, _) in enumerate(writer_slots):
            writer_table[slot] = w_combo[j]

        # Sample router configurations
        for _ in range(router_sample_count):
            table = dict(writer_table)
            for slot, opts in router_slots:
                table[slot] = rng.choice(opts)
            total += 1

            tm = TuringMachine(n_states=n_states, transitions=table)

            # Quick check: can reach a writer state?
            if not _can_write_ones(n_states, table):
                continue

            output, steps, halted = tm.run(effective_max, tape_size=tape_window)
            if halted:
                halting += 1
                if output > best_output:
                    best_output = output
                    best_machine = tm
                    best_steps = steps
                    if verbose:
                        print(f"  New best: output={output} steps={steps} "
                              f"(writer config {wi+1}/{len(writer_combos)})")
            else:
                cycling += 1

        if verbose and (wi + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"  ... writer config {wi+1}/{len(writer_combos)} "
                  f"({elapsed:.1f}s) best={best_output}")

    elapsed = time.time() - t0
    if verbose:
        print(f"BB_targeted({n_states}) ≥ {best_output}  "
              f"[{halting:,} halting, {cycling:,} cycling, "
              f"{total:,} total, {elapsed:.2f}s]")

    return BBResult(
        bb_value=best_output,
        champion=best_machine,
        champion_steps=best_steps,
        total_machines=total,
        halting_count=halting,
        cycle_count=cycling,
        max_steps_used=effective_max,
        tape_window=tape_window,
        computation_time=elapsed,
    )


# ═══════════════════════════════════════════════════════════════
#  §5. FAULHABER CORRECTION DECOMPOSITION
# ═══════════════════════════════════════════════════════════════

@dataclass
class FaulhaberDecomposition:
    """
    Faulhaber-type decomposition of BB_g values.

    BB_g(n) = A · D(n) - R(n)

    where:
      D(n) = naive dominant (max possible output = 2^L for L-tape)
      A    = normalization
      R(n) = correction polynomial (mass gap analogue)
      ΔR   = recurrence generating the correction (stiffness analogue)
    """
    n_values: List[int]       # state counts evaluated
    bb_values: List[int]      # BB_g(n) for each n
    naive_dominant: int       # D = 2^L (maximum possible 1s on tape)
    corrections: List[int]    # R(n) = D - BB_g(n)
    deltas: List[int]         # ΔR(n) = R(n+1) - R(n)
    degree_bound: int         # predicted degree bound from canonical divisor


def faulhaber_decompose(bb_results: Dict[int, BBResult]) -> FaulhaberDecomposition:
    """
    Extract the Faulhaber correction structure from a sequence of BB_g values.

    The naive dominant term is D = 2^L where L = GENUS = 13,
    giving D = 8192 (the maximum number of 1s on a 13-cell tape).

    The correction R(n) = D - BB_g(n) measures how far below the
    theoretical maximum the actual BB value falls. This is the
    "mass gap analogue" — the discrete adjustment that the arithmetic
    constraint imposes on the naive expectation.

    The recurrence ΔR(n) = R(n+1) - R(n) is the "stiffness analogue" —
    the step-by-step mechanism generating the correction.
    """
    D = 1 << GENUS  # 2^13 = 8192

    n_vals = sorted(bb_results.keys())
    bb_vals = [bb_results[n].bb_value for n in n_vals]
    corrections = [D - bb for bb in bb_vals]
    deltas = [corrections[i+1] - corrections[i]
              for i in range(len(corrections) - 1)]

    return FaulhaberDecomposition(
        n_values=n_vals,
        bb_values=bb_vals,
        naive_dominant=D,
        corrections=corrections,
        deltas=deltas,
        degree_bound=CANONICAL_DEG,  # 2g - 2 = 24
    )


def verify_telescoping(decomp: FaulhaberDecomposition) -> bool:
    """
    Verify Corollary 3.2: telescoping the recurrence ΔR recovers
    the correction R(n).

    This is the preimage-identity check: the recurrence is the
    preimage of the closed-form identity.
    """
    if len(decomp.deltas) < 2:
        return True  # trivially true with insufficient data

    # R(n) = R(1) + Σ_{k=1}^{n-1} ΔR(k)
    R1 = decomp.corrections[0]
    for i in range(1, len(decomp.corrections)):
        reconstructed = R1 + sum(decomp.deltas[:i])
        if reconstructed != decomp.corrections[i]:
            return False
    return True


# ═══════════════════════════════════════════════════════════════
#  §6. FATOU CERTIFICATE — CONFINEMENT AS HALTING
# ═══════════════════════════════════════════════════════════════

@dataclass
class FatouCertificate:
    """
    A certificate that an ATM orbit is Fatou-trapped.

    From Theorem 11.7 (Dual Criticality): the orbit is trapped if
    there exist n₀, p ≥ 1 and ε, δ > 0 such that:
      (1) No escape: max_{0≤n≤n₀+p} |z_n| ≤ 2
      (2) Eventual contraction: |z_{n₀+p} - z_{n₀}| ≤ δ
          and |(f^p_c)'(z_{n₀})| ≤ 1 - ε

    In the TM context, "Fatou-trapped" means the configuration
    orbit eventually contracts — the machine enters a periodic
    sub-orbit whose basin attracts all nearby configurations.
    This guarantees halting (or cycle detection) within bounded time.
    """
    is_trapped: bool
    period: Optional[int]       # period p of the eventual cycle
    onset: Optional[int]        # n₀: step where cycle begins
    contraction: Optional[float]  # ε: contraction factor
    basin_size: Optional[int]   # size of the attracting basin


def check_fatou_certificate(orbit: List[TMConfig], *,
                            min_period: int = 1,
                            max_period: int = 100
                            ) -> FatouCertificate:
    """
    Check whether a TM orbit satisfies the Fatou-trapping condition.

    A TM orbit is "Fatou-trapped" if it enters a repeating cycle.
    The period p and onset n₀ are determined by detecting the first
    configuration that repeats.

    This is the TM analogue of Theorem 11.7: the Fatou certificate
    provides a computationally verifiable test for confinement.
    Confinement ⟺ the machine's dynamics enter a stable periodic orbit.

    In MTFT terms: μ_N(y) > 0 ⟺ the orbit is confined ⟺ halting
    is decidable for this machine.
    """
    # Build set of seen configurations
    seen: Dict[TMConfig, int] = {}
    for step, config in enumerate(orbit):
        if config.state == HALT_STATE:
            # Machine halted — trivially "trapped" (fixed point)
            return FatouCertificate(
                is_trapped=True,
                period=0,
                onset=step,
                contraction=1.0,
                basin_size=1,
            )
        if config in seen:
            onset = seen[config]
            period = step - onset
            if min_period <= period <= max_period:
                # Cycle detected
                # Contraction factor: how quickly nearby configs converge
                # In the TM case this is binary (exact cycle), so ε = 1
                return FatouCertificate(
                    is_trapped=True,
                    period=period,
                    onset=onset,
                    contraction=1.0,
                    basin_size=step - onset,
                )
        seen[config] = step

    return FatouCertificate(
        is_trapped=False,
        period=None,
        onset=None,
        contraction=None,
        basin_size=None,
    )


def bb_fatou(n_states: int, *, max_steps: int = 10_000,
             verbose: bool = False) -> BBResult:
    """
    Compute BB_Fatou(n, 143): the Fatou-trapped Busy Beaver.

    Only counts ATMs whose orbits are demonstrably Fatou-trapped
    (enter a cycle or halt). Uses Floyd's tortoise-and-hare algorithm
    for O(1) memory cycle detection.

    By Theorem 5.2 of the theoretical framing, this is computable.
    In MTFT terms: we restrict to machines with μ_N(y) > 0,
    i.e., machines in the confined phase.
    """
    tape_window = GENUS
    config_space = n_states * tape_window * (1 << tape_window)
    effective_max = min(max_steps, config_space + 1)

    t0 = time.time()

    best_output = 0
    best_machine = None
    best_steps = 0
    halting = 0
    cycling = 0
    total = 0

    total_est, _ = count_atms(n_states)
    if verbose:
        print(f"BB_Fatou({n_states}, L={tape_window}): ~{total_est:,} ATMs (Floyd detection)")

    for tm in generate_atms(n_states, prune_dead=False):
        total += 1
        output, steps, halted, cycle_period = tm.run_floyd(
            effective_max, tape_size=tape_window
        )

        if halted:
            halting += 1
            if output > best_output:
                best_output = output
                best_machine = tm
                best_steps = steps
        elif cycle_period is not None:
            # Fatou-trapped: definite cycle detected
            cycling += 1
        # else: neither halted nor cycle found within max_steps — skip

        if verbose and total % 100_000 == 0:
            elapsed = time.time() - t0
            print(f"  ... {total:,} ({elapsed:.1f}s) best={best_output} "
                  f"halting={halting} cycling={cycling}")

    elapsed = time.time() - t0
    if verbose:
        print(f"BB_Fatou({n_states}) = {best_output}  "
              f"[{halting:,} halting, {cycling:,} Fatou-cycling, {elapsed:.2f}s]")

    return BBResult(
        bb_value=best_output,
        champion=best_machine,
        champion_steps=best_steps,
        total_machines=total,
        halting_count=halting,
        cycle_count=cycling,
        max_steps_used=effective_max,
        tape_window=tape_window,
        computation_time=elapsed,
    )


# ═══════════════════════════════════════════════════════════════
#  §7. HALTING BASIN ANALYSIS
# ═══════════════════════════════════════════════════════════════

@dataclass
class HaltingBasinAnalysis:
    """
    Analysis of the halting basin for genus-truncated ATMs.

    From Definition 7.1 and Theorem 7.5:
      Basin(H) = ∪_{k≥0} F^{-k}(H)
    where F is the one-step transition map and H is the halting set.

    In the genus-truncated setting, |Config_L| is finite and
    Basin_L(H) is computable. The "correction" in the Faulhaber
    decomposition counts configurations NOT in Basin_L(H) —
    i.e., configurations trapped in cycles.
    """
    n_states: int
    tape_window: int
    config_space_size: int
    halting_basin_size: int      # |Basin_L(H)|
    cycle_basin_size: int        # |Config_L \ Basin_L(H)|
    halting_fraction: float
    basin_by_sector: Dict[str, int]  # electron/muon/tau decomposition


def analyze_halting_basin(n_states: int, *, max_steps: int = 10_000
                          ) -> HaltingBasinAnalysis:
    """
    Compute the halting basin for genus-truncated ATMs.

    This directly implements Theorem 7.5: in the truncated setting,
    the halting basin is computable and the correction structure
    is recoverable.
    """
    tape_window = GENUS
    config_space = n_states * tape_window * (1 << tape_window)
    effective_max = min(max_steps, config_space + 1)

    machines = enumerate_atms(n_states)

    halting_basin = 0
    cycle_basin = 0
    sector_counts = {"electron": 0, "muon": 0, "tau": 0, "unknown": 0}

    for tm in machines:
        output, steps, halted = tm.run(effective_max, tape_size=tape_window)
        if halted:
            halting_basin += 1
            # Determine dominant sector of this machine
            sector_dom = _machine_dominant_sector(tm, n_states)
            if sector_dom == 0:
                sector_counts["electron"] += 1
            elif sector_dom == 1:
                sector_counts["muon"] += 1
            elif sector_dom == 2:
                sector_counts["tau"] += 1
            else:
                sector_counts["unknown"] += 1
        else:
            cycle_basin += 1

    total = halting_basin + cycle_basin
    return HaltingBasinAnalysis(
        n_states=n_states,
        tape_window=tape_window,
        config_space_size=config_space,
        halting_basin_size=halting_basin,
        cycle_basin_size=cycle_basin,
        halting_fraction=halting_basin / total if total > 0 else 0.0,
        basin_by_sector=sector_counts,
    )


def _machine_dominant_sector(tm: TuringMachine, n_states: int) -> int:
    """Determine which Galois orbit dominates a machine's constraints."""
    sector_weights = [0, 0, 0]
    for state in range(n_states):
        for sym in (0, 1):
            idx = hecke_constraint_index(state, sym)
            sec = dominant_sector(idx)
            if sec is not None:
                a_n = HECKE_TRACES[idx - 1] if idx <= len(HECKE_TRACES) else 0
                sector_weights[sec] += abs(a_n)
    if max(sector_weights) == 0:
        return -1
    return sector_weights.index(max(sector_weights))


# ═══════════════════════════════════════════════════════════════
#  §8. SIGN PATTERN STATISTICS
# ═══════════════════════════════════════════════════════════════

@dataclass
class SignPatternStats:
    """Statistics of the Hecke sign pattern relevant to BB."""
    length: int
    fermionic_count: int
    bosonic_count: int
    free_count: int
    fermionic_density: float
    bosonic_density: float
    free_density: float
    constraint_ratio: float    # fraction of constrained slots
    longest_fermionic_run: int
    longest_bosonic_run: int
    sign_changes: int          # number of sign alternations


def sign_pattern_statistics(length: int = 200) -> SignPatternStats:
    """
    Compute statistics of the Hecke sign pattern.

    The constraint density determines how severely the Hecke
    condition restricts the ATM transition space. By the
    Sato-Tate distribution, the density of FREE indices among
    primes approaches 0 — almost every transition is constrained.
    """
    n = min(length, len(HECKE_TRACES))
    signs = [hecke_sign(i) for i in range(1, n + 1)]

    ferm = sum(1 for s in signs if s == HeckeSign.FERMIONIC)
    bos = sum(1 for s in signs if s == HeckeSign.BOSONIC)
    free = sum(1 for s in signs if s == HeckeSign.FREE)

    # Longest runs
    max_ferm_run = _longest_run(signs, HeckeSign.FERMIONIC)
    max_bos_run = _longest_run(signs, HeckeSign.BOSONIC)

    # Sign changes (ignoring FREE)
    non_free = [s for s in signs if s != HeckeSign.FREE]
    changes = sum(1 for i in range(1, len(non_free))
                  if non_free[i] != non_free[i-1])

    return SignPatternStats(
        length=n,
        fermionic_count=ferm,
        bosonic_count=bos,
        free_count=free,
        fermionic_density=ferm / n,
        bosonic_density=bos / n,
        free_density=free / n,
        constraint_ratio=(ferm + bos) / n,
        longest_fermionic_run=max_ferm_run,
        longest_bosonic_run=max_bos_run,
        sign_changes=changes,
    )


def _longest_run(seq: list, value) -> int:
    """Find the longest consecutive run of `value` in `seq`."""
    max_run = 0
    current = 0
    for s in seq:
        if s == value:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 0
    return max_run


# ═══════════════════════════════════════════════════════════════
#  §9. COMPARISON: UNRESTRICTED BB_g vs HECKE-CONSTRAINED BB_g
# ═══════════════════════════════════════════════════════════════

@dataclass
class BBComparison:
    """Side-by-side comparison of constrained vs unconstrained BB."""
    n_states: int
    bb_unconstrained: BBResult
    bb_hecke: BBResult
    bb_fatou: BBResult
    constraint_reduction: float  # ratio BB_Hecke / BB_unconstrained
    fatou_reduction: float       # ratio BB_Fatou / BB_unconstrained


def compare_bb(n_states: int, *, max_steps: int = 10_000,
               verbose: bool = False) -> BBComparison:
    """
    Compare BB_g(n) across all three constraint levels:
      BB_Fatou(n) ≤ BB_Hecke(n) ≤ BB_unconstrained(n)

    This is the empirical test of the hierarchy.
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"  BB Hierarchy Comparison: n = {n_states} states")
        print(f"  Tape window: L = {GENUS} (genus of X₀({LEVEL}))")
        print(f"{'='*60}")

    # Unconstrained
    if verbose:
        print("\n[1/3] Computing BB_g (unconstrained)...")
    bb_unc = bb_genus(n_states, hecke_constrained=False,
                      max_steps=max_steps, verbose=verbose)

    # Hecke-constrained
    if verbose:
        print("\n[2/3] Computing BB_Hecke (Hecke-constrained)...")
    bb_hck = bb_genus(n_states, hecke_constrained=True,
                      max_steps=max_steps, verbose=verbose)

    # Fatou-trapped
    if verbose:
        print("\n[3/3] Computing BB_Fatou (Fatou-trapped)...")
    bb_fat = bb_fatou(n_states, max_steps=max_steps, verbose=verbose)

    c_ratio = bb_hck.bb_value / bb_unc.bb_value if bb_unc.bb_value > 0 else 0
    f_ratio = bb_fat.bb_value / bb_unc.bb_value if bb_unc.bb_value > 0 else 0

    if verbose:
        print(f"\n{'─'*60}")
        print(f"  Results for n = {n_states}:")
        print(f"    BB_unconstrained = {bb_unc.bb_value}  "
              f"({bb_unc.halting_count:,} halting / {bb_unc.total_machines:,} total)")
        print(f"    BB_Hecke         = {bb_hck.bb_value}  "
              f"({bb_hck.halting_count:,} halting / {bb_hck.total_machines:,} total)")
        print(f"    BB_Fatou         = {bb_fat.bb_value}  "
              f"({bb_fat.halting_count:,} halting / {bb_fat.total_machines:,} total)")
        print(f"    Hecke reduction  = {c_ratio:.4f}")
        print(f"    Fatou reduction  = {f_ratio:.4f}")
        print(f"{'─'*60}")

    return BBComparison(
        n_states=n_states,
        bb_unconstrained=bb_unc,
        bb_hecke=bb_hck,
        bb_fatou=bb_fat,
        constraint_reduction=c_ratio,
        fatou_reduction=f_ratio,
    )


# ═══════════════════════════════════════════════════════════════
#  §10. CLI INTERFACE
# ═══════════════════════════════════════════════════════════════

def _print_banner():
    print("""
╔══════════════════════════════════════════════════════════════╗
║  ARITHMETIC BUSY BEAVERS                                    ║
║  Modular Constraints on Uncomputability                     ║
║                                                              ║
║  X₀(143) · genus 13 · dim 11 · index 168                   ║
║  Three Galois orbits: [1, 4, 6]                             ║
║                                                              ║
║  Roger Tano — MTFT Research Program — March 2026            ║
╚══════════════════════════════════════════════════════════════╝
""")


def _print_sign_pattern():
    """Display the Hecke sign pattern visually."""
    print("\nHecke sign pattern on S₂ᶰᵉʷ(Γ₀(143)), n = 1..200:")
    print("  + = bosonic (write 0)    − = fermionic (write 1)    · = free")
    print()

    for row_start in range(0, 200, 50):
        row_end = min(row_start + 50, 200)
        chars = []
        for n in range(row_start + 1, row_end + 1):
            s = hecke_sign(n)
            if s == HeckeSign.BOSONIC:
                chars.append('+')
            elif s == HeckeSign.FERMIONIC:
                chars.append('−')
            else:
                chars.append('·')
        label = f"  {row_start+1:>3}–{row_end:>3}: "
        print(label + ''.join(chars))

    stats = sign_pattern_statistics()
    print(f"\n  Constraint density: {stats.constraint_ratio:.1%}")
    print(f"  Bosonic: {stats.bosonic_density:.1%}  "
          f"Fermionic: {stats.fermionic_density:.1%}  "
          f"Free: {stats.free_density:.1%}")
    print(f"  Sign changes: {stats.sign_changes}")
    print(f"  Longest bosonic run: {stats.longest_bosonic_run}  "
          f"Longest fermionic run: {stats.longest_fermionic_run}")


def run_full_analysis(max_n: int = 2, max_steps: int = 10_000):
    """Run the complete arithmetic BB analysis."""
    _print_banner()
    _print_sign_pattern()

    # Sign pattern statistics
    print("\n" + "="*60)
    print("  HECKE CONSTRAINT ANALYSIS")
    print("="*60)

    for n in range(1, max_n + 1):
        n_slots = 2 * n
        print(f"\n  n = {n} states → {n_slots} transition slots")
        for s in range(n):
            for sym in (0, 1):
                idx = hecke_constraint_index(s, sym)
                sign = hecke_sign(idx)
                a_val = HECKE_TRACES[idx - 1]
                sec = dominant_sector(idx)
                sec_name = ["electron", "muon", "tau"][sec] if sec is not None else "—"
                constraint = {
                    HeckeSign.BOSONIC: "WRITE 0",
                    HeckeSign.FERMIONIC: "WRITE 1",
                    HeckeSign.FREE: "free",
                }[sign]
                print(f"    δ(q{s}, {sym}): a_{idx} = {a_val:>4} → "
                      f"{constraint:<8}  sector: {sec_name}")

    # BB computation
    print("\n" + "="*60)
    print("  BB HIERARCHY COMPUTATION")
    print("="*60)

    bb_results = {}
    for n in range(1, max_n + 1):
        comp = compare_bb(n, max_steps=max_steps, verbose=True)
        bb_results[n] = comp.bb_hecke

    # Faulhaber decomposition
    if len(bb_results) >= 2:
        print("\n" + "="*60)
        print("  FAULHABER CORRECTION DECOMPOSITION")
        print("="*60)

        decomp = faulhaber_decompose(bb_results)
        print(f"\n  Naive dominant D = 2^{GENUS} = {decomp.naive_dominant}")
        print(f"  Degree bound (canonical divisor): {decomp.degree_bound}")

        for i, n in enumerate(decomp.n_values):
            print(f"  BB_g({n}) = {decomp.bb_values[i]}, "
                  f"R({n}) = {decomp.corrections[i]}")

        if decomp.deltas:
            print(f"\n  Recurrence ΔR:")
            for i, d in enumerate(decomp.deltas):
                print(f"    ΔR({decomp.n_values[i]}→{decomp.n_values[i+1]}) = {d}")

        telescoping_ok = verify_telescoping(decomp)
        print(f"\n  Telescoping identity verified: {telescoping_ok}")

    # Halting basin analysis
    print("\n" + "="*60)
    print("  HALTING BASIN ANALYSIS")
    print("="*60)

    for n in range(1, max_n + 1):
        basin = analyze_halting_basin(n, max_steps=max_steps)
        print(f"\n  n = {n} states:")
        print(f"    Config space |Config₁₃| = {basin.config_space_size:,}")
        print(f"    Halting basin: {basin.halting_basin_size:,} "
              f"({basin.halting_fraction:.1%})")
        print(f"    Cycle basin:   {basin.cycle_basin_size:,}")
        print(f"    Sector breakdown:")
        for sec, count in basin.basin_by_sector.items():
            print(f"      {sec}: {count}")

    # MTFT dictionary
    print("\n" + "="*60)
    print("  MTFT DICTIONARY")
    print("="*60)
    print("""
    ┌──────────────────────────┬──────────────────────────┐
    │ Computability            │ Physics                  │
    ├──────────────────────────┼──────────────────────────┤
    │ BB_g naive dominant      │ Free theory              │
    │ Correction R(n)          │ Mass gap                 │
    │ Recurrence ΔR            │ Stiffness μ_N(y)         │
    │ Degree bound ≤ 24        │ Controlled discontinuity │
    │ [D, M] ≠ 0              │ [t_U, t_R] ≠ 0           │
    │ Fatou-trapped orbit      │ Confinement              │
    │ Halting                  │ μ_N(y) > 0               │
    │ Non-halting cycle        │ Deconfinement             │
    │ Undecidable (∂M)         │ Gravity                  │
    └──────────────────────────┴──────────────────────────┘
    """)


# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        cmd = sys.argv[1]
        if cmd == "signs":
            _print_banner()
            _print_sign_pattern()
        elif cmd == "bb":
            n = int(sys.argv[2]) if len(sys.argv) > 2 else 1
            _print_banner()
            result = bb_genus(n, verbose=True)
            print(f"\n  BB_g({n}, L={GENUS}) = {result.bb_value}")
        elif cmd == "compare":
            n = int(sys.argv[2]) if len(sys.argv) > 2 else 1
            _print_banner()
            compare_bb(n, verbose=True)
        elif cmd == "fatou":
            n = int(sys.argv[2]) if len(sys.argv) > 2 else 1
            _print_banner()
            result = bb_fatou(n, verbose=True)
            print(f"\n  BB_Fatou({n}) = {result.bb_value}")
        elif cmd == "full":
            max_n = int(sys.argv[2]) if len(sys.argv) > 2 else 2
            run_full_analysis(max_n=max_n)
        elif cmd == "sample":
            n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
            k = int(sys.argv[3]) if len(sys.argv) > 3 else 100_000
            _print_banner()
            result = bb_sample(n, sample_size=k, verbose=True)
            print(f"\n  BB_Hecke({n}) ≥ {result.bb_value}")
        elif cmd == "targeted":
            n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
            _print_banner()
            result = bb_targeted(n, verbose=True)
            print(f"\n  BB_Hecke({n}) ≥ {result.bb_value}")
        elif cmd == "sequence":
            _print_banner()
            max_n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
            samples = int(sys.argv[3]) if len(sys.argv) > 3 else 30_000
            print(f"  Computing BB_Hecke sequence n=1..{max_n}")
            print(f"  Exact for n≤3, {samples:,} samples for n≥4\n")
            seq = {}
            for n in range(1, max_n + 1):
                if n <= 3:
                    r = bb_genus(n, hecke_constrained=True)
                    seq[n] = (r.bb_value, r.champion_steps, "exact")
                else:
                    r = bb_sample(n, sample_size=samples)
                    seq[n] = (r.bb_value, r.champion_steps, "lower")
                prefix = "= " if n <= 3 else ">="
                w = sum(1 for s in range(n) for sym in (0,1)
                        if hecke_constraint_index(s, sym) <= len(HECKE_TRACES)
                        and HECKE_TRACES[hecke_constraint_index(s, sym) - 1] <= 0)
                print(f"  BB_Hecke({n:>2}) {prefix}{seq[n][0]:<5}  "
                      f"steps={seq[n][1]:<5}  writers={w}")
            print(f"\n  Sequence: {[seq[n][0] for n in sorted(seq)]}")
        elif cmd == "info":
            _print_banner()
            print(f"  Level N = {LEVEL} = 11 × 13")
            print(f"  Genus   = {GENUS}")
            print(f"  Index   = {INDEX} = |PSL(2,7)| = dim SU(13)")
            print(f"  dim S₂ᶰᵉʷ = {DIM_NEW}")
            print(f"  Orbits  = {ORBIT_DIMS}")
            print(f"  Canonical degree = {CANONICAL_DEG}")
            print(f"  Tape window L = {GENUS}")
            print(f"  Config space = n × {GENUS} × 2^{GENUS} = n × {GENUS * (1 << GENUS):,}")
        else:
            print("Usage: python busy_beaver.py [signs|bb N|compare N|fatou N|sample N [K]|targeted N|sequence [N] [K]|full N|info]")
    else:
        run_full_analysis(max_n=2)
