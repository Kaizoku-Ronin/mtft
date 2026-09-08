"""Exact Bernoulli residues and prime regularity, using integer arithmetic."""
from math import isqrt
import numpy as np


def prime_sieve(limit):
    """Boolean primality for integers 0 <= n < limit."""
    flags = np.ones(limit, dtype=bool)
    flags[:2] = False
    for p in range(2, isqrt(limit - 1) + 1):
        if flags[p]:
            flags[p * p::p] = False
    return flags


def factorial_tables(p):
    facts = np.ones(p - 1, dtype=np.int64)
    for k in range(1, p - 1):
        facts[k] = int(facts[k - 1]) * k % p
    invfacts = np.ones(p - 1, dtype=np.int64)
    invfacts[-1] = pow(int(facts[-1]), -1, p)
    for k in range(p - 2, 0, -1):
        invfacts[k - 1] = int(invfacts[k]) * k % p
    return facts, invfacts


def bernoulli_residues(p):
    """[B_2, B_4, ..., B_(p-3)] modulo an odd prime p.

    Newton series inversion: b <- b * (2 - A*b), doubling precision.
    No floating-point operation enters the computation.
    """
    p = int(p)
    if p < 3 or p % 2 == 0:
        raise ValueError('p must be an odd prime')
    size = p - 2
    if size * (p - 1) ** 2 >= np.iinfo(np.int64).max:
        raise OverflowError('integer convolution could overflow int64')
    facts, invfacts = factorial_tables(p)
    a = invfacts[1:p - 1]
    b = np.ones(1, dtype=np.int64)
    while len(b) < size:
        m = min(2 * len(b), size)
        ab = np.convolve(a[:m], b)[:m] % p
        correction = (-ab) % p
        correction[0] = (int(correction[0]) + 2) % p
        b = np.convolve(b, correction)[:m] % p
    return (b * facts[:size] % p)[2::2]


def verify_series_product(p, even_residues):
    """Check the full inverse identity, reconstructed from the saved residues.

    B_0=1, B_1=-1/2, and odd B_k=0 for k>1. A has constant term 1,
    so a verified inverse is unique in F_p[t]/(t^(p-2)).
    """
    p = int(p)
    size = p - 2
    if len(even_residues) != (p - 3) // 2:
        return False
    if size * (p - 1) ** 2 >= np.iinfo(np.int64).max:
        raise OverflowError('integer convolution could overflow int64')
    _, invfacts = factorial_tables(p)
    b = np.zeros(size, dtype=np.int64)
    b[0] = 1
    if size > 1:
        b[1] = -pow(2, -1, p) % p
    b[2::2] = np.asarray(even_residues, dtype=np.int64) * invfacts[2:size:2] % p
    product = np.convolve(invfacts[1:p - 1], b)[:size] % p
    return bool(product[0] == 1 and np.all(product[1:] == 0))


def power_sum_residue(p, k):
    """Independent modular-power route to B_k mod p in the admissible range."""
    p, k = int(p), int(k)
    if not (2 <= k <= p - 3 and k % 2 == 0):
        raise ValueError('need even 2 <= k <= p-3')
    remainder = sum(pow(a, k, p * p) for a in range(1, p)) % (p * p)
    if remainder % p:
        raise ArithmeticError('power sum is not divisible by p')
    return remainder // p


def witness_indices(residues):
    return (2 + 2 * np.flatnonzero(np.asarray(residues) == 0)).tolist()
