"""
MTFT Cryptographic Primitives
==============================

Experimental crypto tools built on MTFT's mathematical structures:
  - ArithmeticHash:    one-way hash from divisor-weighted w_n
  - BurningShipPRNG:   chaotic PRNG from fractal orbits
  - ModularKeyExchange: Diffie-Hellman on non-abelian SL(2,Z)
  - ArithmeticLattice:  LWE instances with structured hardness

⚠️  RESEARCH / EDUCATIONAL USE ONLY.  Not audited for production.

  pip install mtft
  python examples/09_cryptography.py
"""

from mtft.crypto import (
    ArithmeticHash, BurningShipPRNG,
    ModularKeyExchange, ArithmeticLattice,
)

# ── 1. Arithmetic hash ───────────────────────────────────────
print("ARITHMETIC HASH (divisor-weighted one-way function)")
print("=" * 55)
h = ArithmeticHash(output_bits=256)

messages = [
    b"Modular Time Field Theory",
    b"Roger Tano",
    b"X0(143)",
    b"pip install mtft",
    b"",
]
for msg in messages:
    digest = h.hexdigest(msg)
    label = msg.decode() if msg else "(empty)"
    print(f"  {label:30s} -> {digest[:40]}...")

# Avalanche test
print(f"\n  Avalanche effect (1 bit flip):")
for msg in [b"MTFT", b"Hello World", b"SU(3) confinement"]:
    av = h.verify_avalanche(msg)
    print(f"    '{msg.decode()}': {av['bits_changed']}/{av['total_bits']} bits "
          f"({av['avalanche_ratio']:.1%})")

print(f"\n  Hardness: inverting w_n requires factoring n -> RSA-equivalent")

# ── 2. Burning Ship PRNG ─────────────────────────────────────
print(f"\nBURNING SHIP PRNG (chaotic orbits)")
print("-" * 55)
for seed in [42, 143, 1093]:
    prng = BurningShipPRNG(seed=seed)
    floats = [prng.random_float() for _ in range(5)]
    formatted = [f"{x:.4f}" for x in floats]
    print(f"  seed={seed:>5d}: {formatted}")

prng = BurningShipPRNG(seed=143)
print(f"\n  Random bytes: {prng.random_bytes(16).hex()}")
print(f"  Random bits:  {bin(prng.random_bits(32))}")

# Uniformity check
vals = [BurningShipPRNG(seed=i).random_float() for i in range(1000)]
mean = sum(vals) / len(vals)
print(f"\n  Uniformity check (1000 seeds):")
print(f"    Mean = {mean:.4f} (ideal: 0.5)")
print(f"    Min  = {min(vals):.4f}, Max = {max(vals):.4f}")

# ── 3. SL(2,Z) key exchange ──────────────────────────────────
print(f"\nSL(2,Z) KEY EXCHANGE (non-abelian Diffie-Hellman)")
print("=" * 55)

# Alice and Bob generate secret words
alice_secret = ModularKeyExchange.generate_secret_word(15, seed=11)
bob_secret = ModularKeyExchange.generate_secret_word(15, seed=13)

# Compute public keys (conjugation: w * G * w^-1)
alice_pub = ModularKeyExchange.public_key(alice_secret)
bob_pub = ModularKeyExchange.public_key(bob_secret)

print(f"  Alice's public key: {alice_pub.flatten().tolist()}")
print(f"  Bob's public key:   {bob_pub.flatten().tolist()}")

# Verify det = 1 (still in SL(2,Z))
alice_det = int(alice_pub[0,0]*alice_pub[1,1] - alice_pub[0,1]*alice_pub[1,0])
bob_det = int(bob_pub[0,0]*bob_pub[1,1] - bob_pub[0,1]*bob_pub[1,0])
print(f"  det(Alice_pub) = {alice_det}  (must be 1)")
print(f"  det(Bob_pub)   = {bob_det}  (must be 1)")

# Shared secrets (non-abelian: conjugacy variant)
alice_shared = ModularKeyExchange.shared_secret(alice_secret, bob_pub)
bob_shared = ModularKeyExchange.shared_secret(bob_secret, alice_pub)
print(f"  Alice derives: {alice_shared.flatten().tolist()}")
print(f"  Bob derives:   {bob_shared.flatten().tolist()}")
print(f"\n  Hardness: conjugacy search problem in SL(2,Z)")

# ── 4. Arithmetic lattice ────────────────────────────────────
print(f"\nARITHMETIC LATTICE PROBLEMS (structured LWE)")
print("-" * 55)
for dim in [8, 16, 32]:
    lat = ArithmeticLattice(dimension=dim, y=0.18)
    svp = lat.shortest_vector_estimate()
    print(f"  dim={dim:3d}:  SVP estimate = {svp:.4e}")

print(f"\n  LWE instance generation:")
lat = ArithmeticLattice(dimension=16, y=0.18)
lwe = lat.generate_lwe_instance(secret_dim=4)
print(f"    A shape:  {lwe['A'].shape}")
print(f"    b shape:  {lwe['b'].shape}")
print(f"    modulus:  q = {lwe['modulus']}")
print(f"    secret:   s = {lwe['secret'][:4]}...")

# Verify: b = A*s + e (mod q)
import numpy as np
b_check = (lwe['A'] @ lwe['secret'] + lwe['noise']) % lwe['modulus']
print(f"    A*s + e mod q = {b_check}")
print(f"    b             = {lwe['b']}")
print(f"    Match: {np.array_equal(b_check, lwe['b'])}")
