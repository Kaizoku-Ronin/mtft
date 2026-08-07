#!/usr/bin/env python3
"""
mtft.coset_reps — PSL(2,p) accounting for the X0(143) coset layer.
==================================================================
MIT License — (c) 2026 Roger Tano — MTFT Research Program

BIJECTION [EXACT]: split torus of PSL(2,p) = (Z/p)*/{±1}; its characters are
the EVEN Dirichlet characters mod p = the SU(p) gauge channels (mtft.lchannels)
= the principal-series parameters.  Discrete series (nonsplit torus) are
GAUGE-INVISIBLE; for p = 3 mod 4 the quadratic character lives there.

C[P1(F11) x P1(F13)] = (1+St11) x (1+St13):  168 = 1 + 11 + 13 + 143,
and 143 = dim(St11 x St13): the MTFT level is the Steinberg (x) Steinberg
block of the coset function space.
"""
PSL_IRREPS = {
 13: [1, 13, 14, 14, 7, 7, 12, 12, 12],   # principal block: 1+St,14,14,7+7
 11: [1, 11, 12, 12, 10, 10, 5, 5],       # principal block: 1+St,12,12
}
GAUGE_VISIBLE = {13: [1, 13, 14, 14, 7, 7], 11: [1, 11, 12, 12]}

def order_psl2(p): return p*(p*p - 1)//2

def sum_squares_check(p):
    return sum(d*d for d in PSL_IRREPS[p]) == order_psl2(p)

def torus_char_count(p):
    """= number of even Dirichlet characters mod p."""
    return (p - 1)//2

def stage_decomposition():
    return {'1': 1, 'St11': 11, 'St13': 13, 'St11xSt13': 143, 'total': 168}
