# RT-0 — the symmetry ceiling of the canonical ideal

**Run against live v0.17.0, 2026-08-18. Script:
`rt0_symmetry_ceiling.py` (imports `mtft.canonical`, no other inputs).**

## The question, made decidable

Sol proposes pre-registering an "SM representation gate suite" and seeing
whether su(3) + su(2) + u(1) emerges from v0.17.0. Before spending a wave
on that, compute the ceiling: the Lie algebra of the stabiliser of I_2,

    stab(I_2) = { M in gl_13 : D_M(Q) in I_2 for all Q in I_2 }

where D_M is the derivation from x -> x + eps M x. Any algebraic group
acting linearly on P^12 and preserving the canonical curve has its Lie
algebra inside this. Scalars always qualify, so dim >= 1 unconditionally.

## Result

| ideal | dim I_2 | dim stab(I_2), p = 2147483647 | p = 1000003 |
|---|---|---|---|
| rational normal curve, deg 12 (**control**) | 66 | **4** | **4** |
| **X0(143) canonical ideal** | 55 | **1** | **1** |

The control is the degree-12 rational normal curve, whose automorphisms
are PGL_2; there the answer must be dim sl_2 + scalars = 4, and it is, at
both primes and from a manifestly integral +-1 basis. So the code detects
continuous symmetry when it exists.

**For X0(143) the answer is 1: scalars only. Nothing continuous acts.**
For reference, dim(su(3) + su(2) + u(1)) = 12.

## Why this was always going to be 1

For a non-hyperelliptic curve the canonical embedding is intrinsic, so the
projective linear maps preserving the canonical image are exactly Aut(X),
which is finite for g >= 2 (Hurwitz: |Aut| <= 84(g-1) = 1008; for X0(N)
with N squarefree it is expected to be just the Atkin-Lehner group,
order 4 here). The computation certifies the infinitesimal version from
the shipped data rather than resting on the citation.

## What this closes, and what it does not

**Closed.** No gauge algebra can act on the 55 quadrics. An "SM
representation gate suite" pre-registered against v0.17.0 fails by
construction, not by evidence — it would burn a pre-registration and
produce a negative that means nothing about MTFT.

**Not closed.** This says nothing about whether continuous symmetry exists
elsewhere in the program. It says the *canonical ideal* is the wrong place
to look, which is worth knowing before the search rather than after.

## Where continuous structure actually lives — RT-1, to pre-register

The Mumford-Tate group of J_0(143). From the certified isogeny blocks
(1 + 4 + 6 + 2 = 13, with the level-11 factor appearing twice), the
predicted derived algebra over Qbar is

    sl_2^(+)12   plus a one-dimensional torus,   dim 37,

indexed by the Galois orbits: 1 (143a1) + 4 (f2) + 6 (f3) + 1 (11a1).
Falsifiable: extra endomorphisms — CM or real multiplication beyond the
Hecke field — on any factor would shrink it. Weight-2 trivial nebentypus
forces totally real Hecke fields (CC-01), so CM is already excluded, but
this should be verified rather than assumed.

If that prediction holds, the honest statement about the Standard Model is
sharper than "no demonstrated color sector": the Mumford-Tate group
supplies **u(1) and abundant su(2), and no su(3) at all**. Every factor is
type A_1. A GL_2-type abelian variety cannot produce A_2. That is a
structural obstruction, not a gap in the search.

## The constructive route — RT-2

If A_2 is wanted, the principled move is rank-3 automorphic data, not more
GL_2 curve geometry. The canonical functoriality is the symmetric square
lift Sym^2 : GL_2 -> GL_3 (Gelbart-Jacquet), cuspidal unless the form is
dihedral.

The arc just completed is the geometric shadow of exactly that lift: the
map Sym^2 H^0(K) -> H^0(2K) whose kernel is the 55 quadrics *is*
symmetric square on the geometric side. So Sym^2 f1, Sym^2 f2, Sym^2 f3
are the natural next objects — computable conductors, functional
equations, cuspidality, and L-values that connect to the Petersson
machinery already in the package.

This does not claim su(3) = colour. It asks whether a rank-3 object exists
at all, which is a prerequisite, and it is decidable.

## One of Sol's questions, answerable now

Sol asks whether the three ghost-only dimensions "carry a coherent
representation" rather than the number 3 being read directly.

They cannot, from the Atkin-Lehner group. (Z/2)^2 is abelian, so **all its
irreducible representations are one-dimensional**. The 3 is a multiplicity
of the (-,+) character — from CI-B, the part of H^0(M_-+) (degree 7)
reachable only through the degree-7 channel L_+- (x) L_--, the other
channel L_++ (x) L_-+ having degree 5 and the newspace supplying only 4 of
its 5 dimensions. Three copies of a one-dimensional character. There is no
representation-theoretic content in the number 3 here, and no reading of
it as generations survives contact with the character table.

## Proposed RT series

| item | content | status |
|---|---|---|
| RT-0 | symmetry ceiling of I_2 = 1 (scalars) | **done, above** |
| RT-1 | Mumford-Tate group of J_0(143); predict sl_2^12 + torus | pre-register |
| RT-2 | Sym^2 lifts to GL_3: conductors, cuspidality, L-functions | pre-register |
| RT-3 | anomaly/hypercharge gate — only meaningful after RT-2 | blocked on RT-2 |

RT-3 is deliberately last. Anomaly cancellation is a test on charges
assigned to representations; with no candidate gauge algebra there is
nothing for it to test, and running it early would produce a
meaningless pass or a meaningless fail.
