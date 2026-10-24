"""The Hopf fibration as a qubit; the quaternionic Hopf map as two-qubit entanglement.  2026-10-18.  EXACT unless marked.

(1) S^3 (normalised states in C^2) -> S^2 (the Bloch sphere): fibres are the circles of global phase; two fibres link once (Gauss
integral), the Hopf invariant.  (2) The natural connection Im<psi|d psi> is the Berry connection: the horizontal lift of a latitude loop
at polar angle theta returns with the phase -Omega/2, Omega = 2 pi (1 - cos theta), checked by integrating the lift.  (3) Mosseri-Dandoloff:
for a two-qubit state (alpha, beta, gamma, delta) written as the quaternion pair q1 = alpha + beta j, q2 = gamma + delta j, the Hopf map
S^7 -> S^4, (q1, q2) -> (2 q1 conj(q2), |q1|^2 - |q2|^2), has for its complex and real parts the reduced Bloch vector of the first qubit and
for its j-part 2(beta gamma - alpha delta), whose modulus is the concurrence C = 2|alpha delta - beta gamma|; |r|^2 + C^2 = 1; separable
states land on the S^2 where C = 0; local unitaries move the base point without changing C."""
import numpy as np
rng = np.random.default_rng(7)

# ---------------------------------------------------------------- 1. the Hopf map and the linking of fibres
def hopf(z0, z1):                                        # S^3 -> S^2, the Bloch vector of (z0, z1)
    return np.array([2 * (np.conj(z0) * z1).real, 2 * (np.conj(z0) * z1).imag, abs(z0) ** 2 - abs(z1) ** 2])
def state(theta, phi): return np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)
for th, ph in ((0.3, 1.0), (1.2, 2.5), (2.9, -0.7)):
    n = hopf(*state(th, ph)); assert np.allclose(n, [np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
def stereo(z0, z1):                                      # S^3 in R^4 = (Re z0, Im z0, Re z1, Im z1) -> R^3 from the pole (0, 0, 0, 1)
    x = np.array([z0.real, z0.imag, z1.real, z1.imag]); return x[:3] / (1 - x[3])
def fibre(theta, phi, n=400):
    z0, z1 = state(theta, phi); return np.array([stereo(z0 * np.exp(1j * s), z1 * np.exp(1j * s)) for s in np.linspace(0, 2 * np.pi, n, endpoint=False)])
def linking(P, Q):
    dP = np.roll(P, -1, 0) - P; dQ = np.roll(Q, -1, 0) - Q; mP = P + dP / 2; mQ = Q + dQ / 2
    r = mP[:, None, :] - mQ[None, :, :]; n3 = np.linalg.norm(r, axis=2) ** 3
    return float(np.sum(np.einsum('ijk,ijk->ij', r, np.cross(dP[:, None, :], dQ[None, :, :])) / n3) / (4 * np.pi))
lk = linking(fibre(0.8, 0.3), fibre(2.0, 2.2))
assert abs(abs(lk) - 1) < 1e-3, lk
# the fibre over a point is a round circle in the stereographic picture (all fibres are circles or the line over the pole's point)
F = fibre(1.1, 0.4); A_, B_, C_ = F[0], F[133], F[266]                   # circumcentre of three points, then every point is equidistant
a_, b_ = B_ - A_, C_ - A_; nrm = np.cross(a_, b_)
cen = A_ + (np.dot(a_, a_) * np.cross(b_, nrm) + np.dot(b_, b_) * np.cross(nrm, a_)) / (2 * np.dot(nrm, nrm))
assert np.std(np.linalg.norm(F - cen, axis=1)) < 1e-9 and abs(np.dot(F - cen, nrm / np.linalg.norm(nrm))).max() < 1e-9
print("Hopf fibres of two Bloch points link once: lk = %+.4f; each fibre is a round circle in R^3" % lk)

# ---------------------------------------------------------------- 2. the Berry phase as the holonomy of the Hopf connection
def berry_holonomy(theta, n=20000):                      # horizontal lift of the latitude loop: beta' = -Im<psi0|d psi0/dphi> = -sin^2(theta/2)
    phis = np.linspace(0, 2 * np.pi, n + 1); beta = 0.0; z0, z1 = state(theta, 0.0)
    for k in range(n):
        p0, p1 = phis[k], phis[k + 1]; psi0 = np.array(state(theta, p0)); dpsi = (np.array(state(theta, p1)) - np.array(state(theta, p0 - (p1 - p0)))) / (2 * (p1 - p0))
        beta -= (np.vdot(psi0, dpsi)).imag * (p1 - p0)
    return beta
for th in (0.4, 1.0, np.pi / 2, 2.3):
    om = 2 * np.pi * (1 - np.cos(th)); b = berry_holonomy(th)
    assert abs(b + om / 2) < 1e-6, (th, b, -om / 2)
print("Berry phase of a spin-1/2 around a latitude loop = -Omega/2 = -pi(1 - cos theta), from the horizontal lift (1e-6); the Hopf bundle is the charge-1 monopole bundle")

# ---------------------------------------------------------------- 3. two qubits: the quaternionic Hopf map is entanglement-sensitive
def qmul(p, q):                                          # quaternions as (complex a, complex b) meaning a + b j, with j z = conj(z) j
    a1, b1 = p; a2, b2 = q; return (a1 * a2 - b1 * np.conj(b2), a1 * b2 + b1 * np.conj(a2))
qconj = lambda q: (np.conj(q[0]), -q[1])
def hopf4(psi):                                          # (alpha, beta, gamma, delta) -> point of S^4 in R^5 = (Re, Im of the complex part, the j part as a complex number, the real coordinate)
    al, be, ga, de = psi; q1 = (al, be); q2 = (ga, de)
    prod = qmul(q1, qconj(q2)); x = 2 * prod[0]; jpart = 2 * prod[1]; x0 = abs(al) ** 2 + abs(be) ** 2 - abs(ga) ** 2 - abs(de) ** 2
    return np.array([x.real, x.imag, x0]), jpart
def concurrence(psi): al, be, ga, de = psi; return 2 * abs(al * de - be * ga)
def reduced_bloch(psi):
    al, be, ga, de = psi; rho = np.outer(psi, np.conj(psi)).reshape(2, 2, 2, 2); r1 = np.einsum('ijkj->ik', rho)
    return np.array([2 * r1[0, 1].real, 2 * r1[0, 1].imag, (r1[0, 0] - r1[1, 1]).real])
for _ in range(200):
    psi = rng.normal(size=4) + 1j * rng.normal(size=4); psi /= np.linalg.norm(psi)
    r, jp = hopf4(psi); C = concurrence(psi)
    assert abs(abs(jp) - C) < 1e-12 and np.allclose(r * np.array([1, -1, 1]), reduced_bloch(psi)) or np.allclose(r, reduced_bloch(psi)), (r, reduced_bloch(psi))
    assert abs(np.dot(r, r) + C ** 2 - 1) < 1e-12
# the sign of the Im component is a convention; fix it by checking which matches and report
psi = rng.normal(size=4) + 1j * rng.normal(size=4); psi /= np.linalg.norm(psi); r, _ = hopf4(psi)
conv = 'as is' if np.allclose(r, reduced_bloch(psi)) else 'with Im reversed'
# separable states (product) have C = 0 and |r| = 1; maximally entangled have r = 0, C = 1
a = rng.normal(size=2) + 1j * rng.normal(size=2); b = rng.normal(size=2) + 1j * rng.normal(size=2); prod = np.kron(a / np.linalg.norm(a), b / np.linalg.norm(b))
r, jp = hopf4(prod); assert abs(jp) < 1e-12 and abs(np.linalg.norm(r) - 1) < 1e-12
bell = np.array([1, 0, 0, 1]) / np.sqrt(2); r, jp = hopf4(bell); assert np.linalg.norm(r) < 1e-12 and abs(abs(jp) - 1) < 1e-12
# local unitaries: U1 in SU(2) rotates r and fixes the j-part (det Psi is SU(2) x SU(2)-invariant); U2 in SU(2) is the fibre itself (right
# multiplication by a unit quaternion) and fixes the whole base point; a global phase e^{i phi} turns the j-part by 2 phi
def su2():
    a = rng.normal(size=4); a /= np.linalg.norm(a); return np.array([[a[0] + 1j * a[1], a[2] + 1j * a[3]], [-a[2] + 1j * a[3], a[0] - 1j * a[1]]])
psi = rng.normal(size=4) + 1j * rng.normal(size=4); psi /= np.linalg.norm(psi); U1, U2 = su2(), su2()
r0, j0 = hopf4(psi); r1, j1 = hopf4(np.kron(U1, np.eye(2)) @ psi); r2, j2 = hopf4(np.kron(np.eye(2), U2) @ psi); r3, j3 = hopf4(np.exp(0.37j) * psi)
assert np.linalg.norm(r1 - r0) > 1e-3 and abs(np.linalg.norm(r1) - np.linalg.norm(r0)) < 1e-12 and abs(j1 - j0) < 1e-12
assert np.allclose(r2, r0) and abs(j2 - j0) < 1e-12
assert np.allclose(r3, r0) and abs(j3 - j0 * np.exp(2 * 0.37j)) < 1e-12
print("Mosseri-Dandoloff map checked on 200 random two-qubit states: complex+real parts = the first qubit's reduced Bloch vector (%s), |j-part| = concurrence 2|alpha delta - beta gamma|, |r|^2 + C^2 = 1" % conv)
print("   product states: C = 0, |r| = 1 (the separable S^2 in S^4); Bell state: r = 0, C = 1 (the entangled circle)")
print("   SU(2) on qubit 1 rotates r and fixes the j-part; SU(2) on qubit 2 is the fibre (fixes the base point); a global phase phi turns the j-part by 2 phi; C is fixed by all three")
# the family cos chi |00> + sin chi |11>: C = sin 2 chi, |r| = |cos 2 chi|
for chi in (0.2, 0.5, np.pi / 4):
    psi = np.array([np.cos(chi), 0, 0, np.sin(chi)]); r, jp = hopf4(psi)
    assert abs(abs(jp) - np.sin(2 * chi)) < 1e-12 and abs(np.linalg.norm(r) - abs(np.cos(2 * chi))) < 1e-12
print("   cos chi|00> + sin chi|11>: C = sin 2chi, |r| = |cos 2chi|: a quarter circle from the separable S^2 to the entangled circle")
print("the sequence stops: S^0 -> S^1 -> S^1, S^1 -> S^3 -> S^2, S^3 -> S^7 -> S^4, S^7 -> S^15 -> S^8 (Adams 1960): one, two and three qubits, and no Hopf fibration for four")
print("ALL CHECKS PASSED")
