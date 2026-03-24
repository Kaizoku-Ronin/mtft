"""
MonsterHash — SL(2,Z)-Sponge Hash Function
=============================================

A cryptographic hash function built from MTFT structures that fixes
the ArithmeticHash avalanche problem (6.8% → target ≥50%).

Architecture: Sponge construction over SL(2,Z)
    - State: 2×2 integer matrix in SL(2,Z) (det = 1)
    - Absorb: input blocks drive S^a · T^b multiplications
    - Permute: multi-round SL(2,Z) mixing + Burning Ship nonlinearity
    - Squeeze: extract output from matrix entries mod large primes

Why this works:
    1. SL(2,Z) multiplication is NON-COMMUTATIVE → order matters → avalanche
    2. Matrix products grow entries exponentially (hyperbolic geometry)
    3. Burning Ship fold provides TRUE nonlinearity (|x| operation)
    4. Multiple rounds ensure full diffusion

The old ArithmeticHash failed because:
    - Single 13-bit rotation (not enough mixing)
    - w_n is smooth (nearby n → similar weights)
    - XOR + add is nearly linear

⚠️  Research primitive. Requires independent audit before production use.

Roger Tano — MTFT Research Program — March 2026
"""

from __future__ import annotations

import struct
import math
from typing import Optional

import numpy as np


# ═══════════════════════════════════════════════════════════════
#  SL(2,Z) Matrix Arithmetic (arbitrary precision integers)
# ═══════════════════════════════════════════════════════════════

def _mat_mul(A: tuple, B: tuple) -> tuple:
    """
    Multiply two 2×2 matrices represented as (a,b,c,d) tuples.
    A = [[a,b],[c,d]], entries are Python ints (arbitrary precision).
    """
    a1, b1, c1, d1 = A
    a2, b2, c2, d2 = B
    return (
        a1*a2 + b1*c2,
        a1*b2 + b1*d2,
        c1*a2 + d1*c2,
        c1*b2 + d1*d2,
    )

def _mat_mod(M: tuple, mod: int) -> tuple:
    """Reduce matrix entries modulo mod."""
    return (M[0] % mod, M[1] % mod, M[2] % mod, M[3] % mod)

# SL(2,Z) generators
_S = (0, -1, 1, 0)     # τ → -1/τ
_T = (1, 1, 0, 1)      # τ → τ+1
_ST = _mat_mul(_S, _T)  # τ → -1/(τ+1)
_TS = _mat_mul(_T, _S)  # τ → τ-1/τ = (τ²-1)/τ
_I = (1, 0, 0, 1)       # Identity


def _s_power(n: int) -> tuple:
    """S^n in SL(2,Z). S has order 4: S^4 = I."""
    n = n % 4
    if n == 0: return _I
    if n == 1: return _S
    if n == 2: return ((-1, 0, 0, -1))  # S² = -I
    return (0, 1, -1, 0)  # S³ = -S^T


def _t_power(n: int) -> tuple:
    """T^n in SL(2,Z). T^n = [[1,n],[0,1]]."""
    return (1, n, 0, 1)


# ═══════════════════════════════════════════════════════════════
#  Burning Ship Nonlinear Fold
# ═══════════════════════════════════════════════════════════════

def _burning_ship_fold(x: int, y: int, mod: int) -> tuple:
    """
    Apply Burning Ship-style nonlinearity to two integers.
    
    The key insight: the absolute value |x| operation is what makes
    the Burning Ship non-analytic. In integer arithmetic, this becomes
    conditional negation — a genuine nonlinear operation that breaks
    any linear/affine structure in the state.
    
    Maps (x, y) → (|x|² - |y|² + c₁, 2|x||y| + c₂) mod p
    where c₁, c₂ are derived from the inputs.
    """
    ax = abs(x % mod - mod // 2)  # Center and take abs → fold
    ay = abs(y % mod - mod // 2)
    
    # Burning Ship iteration step (integer version)
    new_x = (ax * ax - ay * ay + x) % mod
    new_y = (2 * ax * ay + y) % mod
    
    return new_x, new_y


# ═══════════════════════════════════════════════════════════════
#  MTFT Round Constants (from arithmetic weights and j-function)
# ═══════════════════════════════════════════════════════════════

def _compute_weight(n: int) -> float:
    """w_n = Σ_{d|n} log(d)/d — MTFT arithmetic weight."""
    s = 0.0
    for d in range(1, n + 1):
        if n % d == 0:
            s += math.log(d) / d
    return s


# Precompute round constants from j-function coefficients and weights
# j(τ) = q⁻¹ + 744 + 196884q + 21493760q² + ...
_J_COEFFS = [
    744, 196884, 21493760, 864299970, 20245856256,
    333202640600, 4252023300096, 44656994071935,
    401490886656000, 3176440229784420, 22567393309593600,
    146211911499519294, 874313719685775360,
]

# Round constants: j-function coefficients XORed with weight-derived integers
_ROUND_CONSTANTS = []
for i, jc in enumerate(_J_COEFFS):
    w = _compute_weight(i + 2)
    w_int = int(w * 2**48) & 0xFFFFFFFFFFFFFFFF
    _ROUND_CONSTANTS.append(jc ^ w_int)

# Additional constants from Monster primes and MTFT structural numbers
_MONSTER_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
_STRUCTURAL = [
    143,      # level = 11 × 13
    168,      # index = |PSL(2,7)|
    196883,   # dim(V♮)
    1093,     # first Wieferich prime
    3511,     # second Wieferich prime
    78,       # arithmetic kernel = gcd(1092, 3510)
    744,      # j-function constant term = 24 × 31
]


# ═══════════════════════════════════════════════════════════════
#  MonsterHash Core
# ═══════════════════════════════════════════════════════════════

# Working modulus: a large prime close to 2^64 for efficiency
# Using 2^61 - 1 (Mersenne prime M61) for fast reduction
_MODULUS = (1 << 61) - 1

# Output modulus for final extraction
_OUTPUT_MOD = (1 << 64) - 59  # largest 64-bit prime


class MonsterHash:
    """
    SL(2,Z)-Sponge Hash Function.
    
    Construction:
        1. Initialize state as S·T in SL(2,Z) (the standard cusp form generator)
        2. For each 64-bit input block:
           a. Decompose block into (a, b) pair
           b. Multiply state by S^a · T^b (absorb)
           c. Apply Burning Ship fold to matrix entries (nonlinear permute)
           d. Reduce entries mod working modulus
        3. Run ROUNDS additional mixing permutations (squeeze prep)
        4. Extract output from matrix entries
    
    The non-commutativity of SL(2,Z) ensures that reordering input blocks
    changes the output completely (avalanche). The Burning Ship fold
    ensures that nearby values diverge exponentially (chaos).
    
    Parameters:
        output_bits: 256 or 512 (default 512)
        rounds: mixing rounds per block (default 13 = genus of X₀(143))
    """
    
    ROUNDS = 13  # genus(X₀(143))
    
    def __init__(self, output_bits: int = 512):
        if output_bits not in (256, 512):
            raise ValueError("output_bits must be 256 or 512")
        self.output_bits = output_bits
        self._output_words = output_bits // 64
    
    def _init_state(self) -> tuple:
        """Initial state: S·T matrix, the generator of the cusp."""
        # Start with [[0,-1],[1,1]] = S·T, reduced mod MODULUS
        return _mat_mod(_ST, _MODULUS)
    
    def _permute(self, state: tuple, round_key: int) -> tuple:
        """
        One round of the permutation.
        
        Combines:
        1. SL(2,Z) multiplication by round-key-derived matrix
        2. Burning Ship nonlinear fold on entries
        3. Cross-mixing of entries
        """
        a, b, c, d = state
        
        # Step 1: Derive SL(2,Z) action from round key
        s_exp = (round_key & 0x3) + 1        # 1-4
        t_exp = ((round_key >> 2) & 0xFFFF) + 1  # 1-65536
        action = _mat_mul(_s_power(s_exp), _t_power(t_exp))
        state = _mat_mul(state, action)
        state = _mat_mod(state, _MODULUS)
        a, b, c, d = state
        
        # Step 2: Burning Ship fold on (a,b) and (c,d) pairs
        a, b = _burning_ship_fold(a, b, _MODULUS)
        c, d = _burning_ship_fold(c, d, _MODULUS)
        
        # Step 3: Cross-mix (like a Feistel round)
        a = (a + c * 196883) % _MODULUS  # Monster dimension
        b = (b ^ (d * 744 % _MODULUS)) % _MODULUS  # j-constant
        c = (c + a * 143) % _MODULUS     # MTFT level
        d = (d ^ (b * 168 % _MODULUS)) % _MODULUS  # MTFT index
        
        # Step 4: Rotate bits within entries for diffusion
        shift = (round_key >> 18) % 53 + 5  # 5-57 bit rotation
        a = ((a << shift) | (a >> (61 - shift))) & ((1 << 61) - 1)
        d = ((d << (61 - shift)) | (d >> shift)) & ((1 << 61) - 1)
        
        return (a % _MODULUS, b % _MODULUS, c % _MODULUS, d % _MODULUS)
    
    def _absorb_block(self, state: tuple, block: int) -> tuple:
        """Absorb one 64-bit block into the state."""
        # Split block into two 32-bit halves
        hi = (block >> 32) & 0xFFFFFFFF
        lo = block & 0xFFFFFFFF
        
        # XOR block data into state entries (injection)
        a, b, c, d = state
        a = (a ^ (hi * 2654435761)) % _MODULUS  # Knuth multiplicative
        b = (b ^ (lo * 2246822519)) % _MODULUS
        c = (c ^ (hi * 3266489917)) % _MODULUS
        d = (d ^ (lo * 668265263)) % _MODULUS
        state = (a, b, c, d)
        
        # Run ROUNDS permutation rounds with round constants
        for r in range(self.ROUNDS):
            rc = _ROUND_CONSTANTS[r % len(_ROUND_CONSTANTS)]
            # Mix in round index and block value for domain separation
            key = rc ^ (r * 0x9E3779B97F4A7C15) ^ block
            state = self._permute(state, key)
        
        return state
    
    def _squeeze(self, state: tuple) -> list:
        """
        Squeeze phase: extract output words from the state.
        
        Run additional permutation rounds, extracting 64-bit words
        from matrix entries after each pair of rounds.
        """
        output = []
        
        for i in range(self._output_words):
            # Two rounds of mixing per output word
            key1 = _STRUCTURAL[i % len(_STRUCTURAL)] ^ (i * 0xDEADBEEF)
            key2 = _STRUCTURAL[(i + 3) % len(_STRUCTURAL)] ^ ((i + 1) * 0xCAFEBABE)
            state = self._permute(state, key1)
            state = self._permute(state, key2)
            
            # Extract from all four entries, mixed
            a, b, c, d = state
            word = ((a ^ b) * 196883 + (c ^ d)) % _OUTPUT_MOD
            output.append(word & 0xFFFFFFFFFFFFFFFF)
        
        return output
    
    def _pad_and_block(self, data: bytes) -> list:
        """
        Pad input and split into 64-bit blocks.
        
        Uses MTFT-specific padding:
        - Append 0x8B (143 in hex... well, 0x8F, but we use 0x80 + len)
        - Append length as 64-bit big-endian
        - Pad to multiple of 8 bytes with Burning Ship marker 0xB5
        """
        # Standard sponge padding
        msg_len = len(data)
        padded = bytearray(data)
        padded.append(0x80)  # Standard padding bit
        
        # Pad to 8-byte alignment (leaving room for 8-byte length)
        while (len(padded) + 8) % 8 != 0:
            padded.append(0x00)
        
        # Append original length as 64-bit big-endian
        padded.extend(msg_len.to_bytes(8, 'big'))
        
        # Convert to 64-bit blocks
        blocks = []
        for i in range(0, len(padded), 8):
            blocks.append(int.from_bytes(padded[i:i+8], 'big'))
        
        return blocks
    
    def digest(self, data: bytes) -> int:
        """
        Compute MonsterHash digest.
        
        Returns an integer of self.output_bits bits.
        """
        state = self._init_state()
        blocks = self._pad_and_block(data)
        
        # Absorb phase
        for block in blocks:
            state = self._absorb_block(state, block)
        
        # Domain separation: absorb the output length
        state = self._absorb_block(state, self.output_bits)
        
        # Squeeze phase
        output_words = self._squeeze(state)
        
        # Combine into single integer
        result = 0
        for word in output_words:
            result = (result << 64) | (word & 0xFFFFFFFFFFFFFFFF)
        
        # Mask to exact output size
        mask = (1 << self.output_bits) - 1
        return result & mask
    
    def hexdigest(self, data: bytes) -> str:
        """Compute hash and return as hex string."""
        d = self.digest(data)
        return format(d, f'0{self.output_bits // 4}x')
    
    def verify_avalanche(self, data: bytes, n_trials: int = 100) -> dict:
        """
        Statistical avalanche test.
        
        Flips each bit position in the first byte and measures
        how many output bits change. Reports mean, min, max, and
        standard deviation of the avalanche ratio.
        
        A good hash should have mean ≈ 0.50, std < 0.05.
        """
        h1 = self.digest(data)
        
        ratios = []
        for bit_pos in range(min(n_trials, len(data) * 8)):
            modified = bytearray(data)
            byte_idx = bit_pos // 8
            bit_idx = bit_pos % 8
            if byte_idx < len(modified):
                modified[byte_idx] ^= (1 << bit_idx)
            h2 = self.digest(bytes(modified))
            
            xor = h1 ^ h2
            flipped = bin(xor).count('1')
            ratios.append(flipped / self.output_bits)
        
        import statistics
        return {
            "hash": format(h1, f'0{self.output_bits // 4}x')[:32] + "...",
            "output_bits": self.output_bits,
            "trials": len(ratios),
            "mean_avalanche": statistics.mean(ratios),
            "min_avalanche": min(ratios),
            "max_avalanche": max(ratios),
            "std_avalanche": statistics.stdev(ratios) if len(ratios) > 1 else 0,
            "ideal": 0.50,
            "pass": statistics.mean(ratios) > 0.40,
        }


# ═══════════════════════════════════════════════════════════════
#  Quick comparison test
# ═══════════════════════════════════════════════════════════════

def compare_hashes():
    """Compare old ArithmeticHash vs new MonsterHash avalanche."""
    from mtft.crypto import ArithmeticHash
    
    msg = b"Modular Time Field Theory - MonsterCoin whitepaper test"
    
    print("=" * 70)
    print("HASH AVALANCHE COMPARISON")
    print("=" * 70)
    
    # Old hash
    print("\n--- ArithmeticHash (OLD) ---")
    old = ArithmeticHash(output_bits=256)
    old_result = old.verify_avalanche(msg)
    print(f"  Hash:      {old_result['original_hash'][:32]}...")
    print(f"  Avalanche: {old_result['avalanche_ratio']:.1%}")
    print(f"  Status:    {'PASS' if old_result['avalanche_ratio'] > 0.4 else 'FAIL'}")
    
    # New hash (256-bit)
    print("\n--- MonsterHash 256-bit (NEW) ---")
    new256 = MonsterHash(output_bits=256)
    new_result = new256.verify_avalanche(msg, n_trials=64)
    print(f"  Hash:      {new_result['hash']}")
    print(f"  Mean:      {new_result['mean_avalanche']:.1%}")
    print(f"  Range:     [{new_result['min_avalanche']:.1%}, {new_result['max_avalanche']:.1%}]")
    print(f"  Std:       {new_result['std_avalanche']:.4f}")
    print(f"  Status:    {'PASS' if new_result['pass'] else 'FAIL'}")
    
    # New hash (512-bit)
    print("\n--- MonsterHash 512-bit (NEW) ---")
    new512 = MonsterHash(output_bits=512)
    new_result512 = new512.verify_avalanche(msg, n_trials=64)
    print(f"  Hash:      {new_result512['hash']}")
    print(f"  Mean:      {new_result512['mean_avalanche']:.1%}")
    print(f"  Range:     [{new_result512['min_avalanche']:.1%}, {new_result512['max_avalanche']:.1%}]")
    print(f"  Std:       {new_result512['std_avalanche']:.4f}")
    print(f"  Status:    {'PASS' if new_result512['pass'] else 'FAIL'}")
    
    # Collision resistance spot check
    print("\n--- Collision Spot Check (10,000 random inputs) ---")
    import os
    hashes = set()
    collisions = 0
    h = MonsterHash(512)
    for i in range(10000):
        data = os.urandom(32) + i.to_bytes(4, 'big')
        digest = h.digest(data)
        if digest in hashes:
            collisions += 1
        hashes.add(digest)
    print(f"  Unique hashes: {len(hashes)}/10000")
    print(f"  Collisions:    {collisions}")
    
    print("\n" + "=" * 70)
    return new_result


if __name__ == "__main__":
    compare_hashes()
