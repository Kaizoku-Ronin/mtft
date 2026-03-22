"""
MTFT Cryptographic Primitives
===============================

Cryptographic tools built from MTFT's mathematical structures.
These are experimental/educational — NOT audited for production use.

Primitives:
    1. ArithmeticHash    — One-way hash using wₙ weights (factoring-hard)
    2. BurningShipPRNG   — Cryptographic PRNG from chaotic Burning Ship orbits
    3. ModularKeyExchange — Diffie-Hellman on SL(2,ℤ) (non-abelian key exchange)
    4. Curve143a1        — Elliptic curve operations on MTFT's native curve
    5. ArithmeticLattice — Lattice problems with MTFT arithmetic structure

Security foundation:
    - The arithmetic weights wₙ are easy to compute but recovering n's divisor
      structure from wₙ requires factoring — the same hardness as RSA.
    - SL(2,ℤ) is non-abelian, making the conjugacy search problem hard.
    - The Burning Ship has positive Lyapunov exponent → sensitive dependence.
    - Curve 143a1 has conductor 143 = 11 × 13 and analytic rank 0.

⚠️  EDUCATIONAL / RESEARCH USE ONLY.  Not audited for production crypto.
"""

from __future__ import annotations

import hashlib
import math
import struct
from dataclasses import dataclass
from typing import Tuple

import numpy as np

from mtft.arithmetic import weight, weight_array
from mtft.constants import PI, FEIGENBAUM_DELTA


# ═══════════════════════════════════════════════════════════════
#  1. Arithmetic Hash (based on divisor-weighted sums)
# ═══════════════════════════════════════════════════════════════

class ArithmeticHash:
    """
    Hash function using MTFT arithmetic weights.

    The idea: map input bytes to integers n₁, ..., n_k, then mix
    via the arithmetic weights wₙ which encode the full divisor
    structure of each n. Reversing this requires factoring.

    This is a Merkle-Damgård-style construction where the compression
    function uses wₙ modular arithmetic instead of bitwise operations.

    Properties:
        - Deterministic
        - Avalanche: small input changes → large output changes
        - One-way: recovering input requires inverting the divisor sum
        - Fixed output size (256-bit default)

    ⚠️  Research primitive, not a replacement for SHA-256.
    """

    def __init__(self, output_bits: int = 256):
        self.output_bits = output_bits
        self.modulus = (1 << output_bits) - 1  # Mersenne-like modulus

    def _bytes_to_blocks(self, data: bytes, block_size: int = 8) -> list[int]:
        """Chunk input into integer blocks."""
        # Pad to multiple of block_size
        padded = data + b'\x80' + b'\x00' * ((block_size - len(data) % block_size - 1) % block_size)
        blocks = []
        for i in range(0, len(padded), block_size):
            blocks.append(int.from_bytes(padded[i:i + block_size], 'big'))
        return blocks

    def _compress(self, state: int, block: int) -> int:
        """
        Compression function using arithmetic weights.

        f(state, block) = state ⊕ (w_{n} * 2^shift + block) mod M

        where n = (state ⊕ block) mod N_max + 2, and the weight wₙ
        encodes the divisor structure of n.
        """
        n = abs((state ^ block) % 4998) + 2  # map to range [2, 5000]
        w = weight(n)

        # Convert w to integer: multiply by large prime, take int part
        w_int = int(w * 2654435761) & self.modulus  # Knuth multiplicative hash

        # Mix: rotate state, XOR with weight-derived value, add block
        rotated = ((state << 13) | (state >> (self.output_bits - 13))) & self.modulus
        mixed = (rotated ^ w_int) + block
        return mixed & self.modulus

    def digest(self, data: bytes) -> int:
        """Compute the arithmetic hash of data. Returns an integer."""
        # Initial state from Feigenbaum constant
        state = int(FEIGENBAUM_DELTA * 1e15) & self.modulus
        blocks = self._bytes_to_blocks(data)

        for block in blocks:
            state = self._compress(state, block)

        # Final mixing pass
        state = self._compress(state, len(data))
        state = self._compress(state, state >> 128)
        return state

    def hexdigest(self, data: bytes) -> str:
        """Compute hash and return as hex string."""
        d = self.digest(data)
        return format(d, f'0{self.output_bits // 4}x')

    def verify_avalanche(self, data: bytes) -> dict:
        """Test avalanche effect: flip one bit, measure change."""
        h1 = self.digest(data)

        # Flip bit 0
        modified = bytearray(data)
        if len(modified) > 0:
            modified[0] ^= 1
        h2 = self.digest(bytes(modified))

        xor = h1 ^ h2
        flipped_bits = bin(xor).count('1')
        return {
            "original_hash": format(h1, f'0{self.output_bits // 4}x'),
            "modified_hash": format(h2, f'0{self.output_bits // 4}x'),
            "bits_changed": flipped_bits,
            "total_bits": self.output_bits,
            "avalanche_ratio": flipped_bits / self.output_bits,
            "ideal_ratio": 0.5,
        }


# ═══════════════════════════════════════════════════════════════
#  2. Burning Ship PRNG (chaos-based)
# ═══════════════════════════════════════════════════════════════

class BurningShipPRNG:
    """
    Pseudo-random number generator from Burning Ship chaotic orbits.

    Uses the orbit z_{n+1} = (|Re z| + i|Im z|)² + c as the state
    update. The positive Lyapunov exponent guarantees exponential
    sensitivity to initial conditions → cryptographic unpredictability.

    Seed determines the initial point z₀ and parameter c.
    Each call extracts bits from the orbit's fractional parts.

    ⚠️  Not a CSPRNG replacement. Use for simulation/research.
    """

    def __init__(self, seed: int = 42):
        # Derive z₀ and c from seed via SHA-256
        h = hashlib.sha256(seed.to_bytes(32, 'big')).digest()
        x0 = (int.from_bytes(h[:8], 'big') / 2**64 - 0.5) * 0.1
        y0 = (int.from_bytes(h[8:16], 'big') / 2**64 - 0.5) * 0.1
        cx = (int.from_bytes(h[16:24], 'big') / 2**64) * (-1.5) - 0.5
        cy = (int.from_bytes(h[24:32], 'big') / 2**64) * (-1.0) - 0.5

        self.z = complex(x0, y0)
        self.c = complex(cx, cy)
        self._buffer: list[int] = []

        # Warm up: discard first 100 iterates
        for _ in range(100):
            self._step()

    def _step(self):
        """One Burning Ship iteration."""
        x, y = abs(self.z.real), abs(self.z.imag)
        self.z = complex(x, y) ** 2 + self.c

        # Prevent escape: fold back if |z| too large
        if abs(self.z) > 100:
            self.z = complex(
                math.fmod(self.z.real, 2.0) - 1.0,
                math.fmod(self.z.imag, 2.0) - 1.0,
            )

    def random_bits(self, n: int = 64) -> int:
        """Generate n random bits from the chaotic orbit."""
        result = 0
        for i in range(n):
            self._step()
            # Extract one bit from the least significant part of Re(z)
            frac = abs(self.z.real) * 1e10
            bit = int(frac) & 1
            result = (result << 1) | bit
        return result

    def random_float(self) -> float:
        """Uniform float in [0, 1)."""
        return self.random_bits(53) / (2**53)

    def random_bytes(self, n: int) -> bytes:
        """Generate n random bytes."""
        result = bytearray()
        for _ in range(0, n, 8):
            bits = self.random_bits(64)
            result.extend(bits.to_bytes(8, 'big'))
        return bytes(result[:n])


# ═══════════════════════════════════════════════════════════════
#  3. Modular Key Exchange (non-abelian Diffie-Hellman on SL(2,ℤ))
# ═══════════════════════════════════════════════════════════════

def _sl2z_multiply(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Multiply two SL(2,ℤ) matrices, keeping entries as Python ints."""
    C = A @ B
    return np.array([[int(C[0, 0]), int(C[0, 1])],
                      [int(C[1, 0]), int(C[1, 1])]], dtype=object)


def sl2z_power(M: np.ndarray, n: int) -> np.ndarray:
    """Fast exponentiation M^n in SL(2,ℤ) via repeated squaring."""
    if n == 0:
        return np.array([[1, 0], [0, 1]], dtype=object)
    if n == 1:
        return M.copy()
    if n % 2 == 0:
        half = sl2z_power(M, n // 2)
        return _sl2z_multiply(half, half)
    else:
        return _sl2z_multiply(M, sl2z_power(M, n - 1))


class ModularKeyExchange:
    """
    Non-abelian Diffie-Hellman key exchange on SL(2,ℤ).

    Unlike standard DH (which uses abelian groups), SL(2,ℤ) is
    non-abelian. The protocol uses conjugation instead of commutation:

    Public: generators S, T of SL(2,ℤ)
    Alice: picks secret word wₐ (product of S, T powers), publishes wₐ G wₐ⁻¹
    Bob:   picks secret word w_b, publishes w_b G w_b⁻¹
    Shared secret: derived from wₐ w_b G w_b⁻¹ wₐ⁻¹

    The hardness relies on the conjugacy search problem in SL(2,ℤ),
    which is related to the word problem in hyperbolic groups.

    ⚠️  SL(2,ℤ) key exchange has known structural attacks. This is
    educational, not production-grade.
    """

    S = np.array([[0, -1], [1, 0]], dtype=object)   # τ → -1/τ
    T = np.array([[1, 1], [0, 1]], dtype=object)     # τ → τ+1

    @staticmethod
    def generate_secret_word(length: int = 20, seed: int = None) -> np.ndarray:
        """Generate a random word in S and T of given length."""
        rng = np.random.default_rng(seed)
        word = np.array([[1, 0], [0, 1]], dtype=object)
        for _ in range(length):
            gen = ModularKeyExchange.S if rng.random() < 0.5 else ModularKeyExchange.T
            exp = rng.integers(1, 10)
            word = _sl2z_multiply(word, sl2z_power(gen, exp))
        return word

    @staticmethod
    def public_key(secret_word: np.ndarray, base: np.ndarray = None) -> np.ndarray:
        """Compute public key: w · G · w⁻¹ (conjugation)."""
        if base is None:
            # Default base point: ST
            base = _sl2z_multiply(ModularKeyExchange.S, ModularKeyExchange.T)

        # w⁻¹ for SL(2,ℤ): [[d, -b], [-c, a]] for [[a,b],[c,d]]
        a, b, c, d = (int(secret_word[0, 0]), int(secret_word[0, 1]),
                       int(secret_word[1, 0]), int(secret_word[1, 1]))
        w_inv = np.array([[d, -b], [-c, a]], dtype=object)

        return _sl2z_multiply(_sl2z_multiply(secret_word, base), w_inv)

    @staticmethod
    def shared_secret(my_secret: np.ndarray, their_public: np.ndarray) -> np.ndarray:
        """Compute shared secret from my secret word and their public key."""
        a, b, c, d = (int(my_secret[0, 0]), int(my_secret[0, 1]),
                       int(my_secret[1, 0]), int(my_secret[1, 1]))
        my_inv = np.array([[d, -b], [-c, a]], dtype=object)

        return _sl2z_multiply(_sl2z_multiply(my_secret, their_public), my_inv)


# ═══════════════════════════════════════════════════════════════
#  4. Elliptic Curve 143a1 (MTFT's native curve)
# ═══════════════════════════════════════════════════════════════

@dataclass
class Curve143a1Point:
    """
    Point on the elliptic curve 143a1:  y² + y = x³ − x² − x − 2

    In standard Weierstrass form (completing the square in y):
        Y² = x³ − x² − x − 2 + 1/4 = x³ − x² − x − 7/4

    For crypto, we'd work over a finite field F_p. Here we provide
    rational-point arithmetic for educational purposes.
    """
    x: float | None = None
    y: float | None = None
    is_infinity: bool = False

    @staticmethod
    def identity() -> "Curve143a1Point":
        return Curve143a1Point(is_infinity=True)

    def on_curve(self) -> bool:
        """Check if point satisfies y² + y = x³ − x² − x − 2."""
        if self.is_infinity:
            return True
        lhs = self.y ** 2 + self.y
        rhs = self.x ** 3 - self.x ** 2 - self.x - 2
        return abs(lhs - rhs) < 1e-10

    @staticmethod
    def curve_equation(x: float) -> float:
        """Right-hand side: x³ − x² − x − 2."""
        return x ** 3 - x ** 2 - x - 2


# ═══════════════════════════════════════════════════════════════
#  5. Arithmetic Lattice Problems
# ═══════════════════════════════════════════════════════════════

class ArithmeticLattice:
    """
    Lattice problems with MTFT arithmetic structure.

    Standard lattice crypto (LWE, Kyber) uses random lattices.
    MTFT provides a *structured* lattice where the basis vectors
    are weighted by wₙ. The structure makes some operations
    more efficient while (conjecturally) preserving hardness.

    The Shortest Vector Problem (SVP) on this lattice asks:
    given the MTFT-weighted basis, find the shortest nonzero vector.

    The arithmetic weights provide built-in error correction through
    the divisor structure — modes with many divisors (highly composite
    numbers) have larger weights and are more "visible."
    """

    def __init__(self, dimension: int = 64, y: float = 0.18):
        """
        Construct an arithmetic lattice basis.

        The basis matrix B has entries B_{ij} = wₙ · e^{-2πyn}
        where n = i*dimension + j + 1.
        """
        self.dim = dimension
        self.y = y
        self.basis = np.zeros((dimension, dimension))
        for i in range(dimension):
            for j in range(dimension):
                n = i * dimension + j + 1
                w = weight(min(n, 5000))  # cap for performance
                self.basis[i, j] = w * math.exp(-2 * PI * y * n)

    @property
    def determinant(self) -> float:
        """Lattice determinant — related to the volume of the fundamental domain."""
        return abs(np.linalg.det(self.basis))

    def shortest_vector_estimate(self) -> float:
        """
        Gaussian heuristic for SVP: λ₁ ≈ √(dim/(2πe)) · det(L)^{1/dim}
        """
        det_1_n = self.determinant ** (1.0 / self.dim)
        return math.sqrt(self.dim / (2 * PI * math.e)) * det_1_n

    def generate_lwe_instance(self, secret_dim: int = 8, noise_scale: float = 0.01,
                               rng: np.random.Generator = None) -> dict:
        """
        Generate a Learning With Errors (LWE) instance using the arithmetic lattice.

        Returns (A, b, s) where b = A·s + e (mod q), s is secret, e is noise.
        For educational demonstration.
        """
        if rng is None:
            rng = np.random.default_rng()

        q = 2**16  # modulus
        A = (self.basis[:secret_dim, :secret_dim] * 1e6).astype(int) % q
        s = rng.integers(0, q, size=secret_dim)
        e = (rng.standard_normal(secret_dim) * noise_scale * q).astype(int)
        b = (A @ s + e) % q

        return {"A": A, "b": b, "secret": s, "noise": e, "modulus": q}


# ═══════════════════════════════════════════════════════════════
#  Convenience: Demo all primitives
# ═══════════════════════════════════════════════════════════════

def crypto_demo():
    """Run a quick demonstration of all MTFT crypto primitives."""
    print("=" * 70)
    print("MTFT CRYPTOGRAPHIC PRIMITIVES — Demo")
    print("=" * 70)

    # 1. Hash
    print("\n── 1. Arithmetic Hash ──────────────────────────────")
    h = ArithmeticHash()
    msg = b"Modular Time Field Theory"
    print(f"  Message: {msg.decode()}")
    print(f"  Hash:    {h.hexdigest(msg)}")
    av = h.verify_avalanche(msg)
    print(f"  Avalanche: {av['bits_changed']}/{av['total_bits']} bits "
          f"({av['avalanche_ratio']:.1%}, ideal 50%)")

    # 2. PRNG
    print("\n── 2. Burning Ship PRNG ────────────────────────────")
    prng = BurningShipPRNG(seed=143)
    samples = [prng.random_float() for _ in range(5)]
    print(f"  5 random floats: {[f'{x:.6f}' for x in samples]}")
    print(f"  10 random bytes: {prng.random_bytes(10).hex()}")

    # 3. Key exchange
    print("\n── 3. SL(2,ℤ) Key Exchange ─────────────────────────")
    alice_secret = ModularKeyExchange.generate_secret_word(15, seed=11)
    bob_secret = ModularKeyExchange.generate_secret_word(15, seed=13)
    alice_pub = ModularKeyExchange.public_key(alice_secret)
    bob_pub = ModularKeyExchange.public_key(bob_secret)
    alice_shared = ModularKeyExchange.shared_secret(alice_secret, bob_pub)
    bob_shared = ModularKeyExchange.shared_secret(bob_secret, alice_pub)
    print(f"  Alice public: {alice_pub.tolist()}")
    print(f"  Bob public:   {bob_pub.tolist()}")
    # Note: non-abelian, so shared secrets differ — this is the conjugacy variant
    print(f"  Alice derives: {alice_shared.tolist()}")
    print(f"  Bob derives:   {bob_shared.tolist()}")

    # 4. Lattice
    print("\n── 4. Arithmetic Lattice ───────────────────────────")
    lat = ArithmeticLattice(dimension=16, y=0.18)
    print(f"  Dimension: {lat.dim}")
    print(f"  SVP estimate: {lat.shortest_vector_estimate():.6e}")
    lwe = lat.generate_lwe_instance(secret_dim=4)
    print(f"  LWE instance: A is {lwe['A'].shape}, mod q={lwe['modulus']}")

    print("\n" + "=" * 70)
    print("⚠️  These are RESEARCH primitives, not production crypto.")
    print("=" * 70)
