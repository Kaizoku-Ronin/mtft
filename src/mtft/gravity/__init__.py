r"""mtft.gravity — knotted weak fields on flat space: the tube, the squared hopfion's lines of force, a frozen tetrad (v1.1.0).

  fields   MODEL    the labels (α, β), the (p, q) light fields F_pq, the squared hopfion W = n ⊗ n/d⁵, the null
                    flow that carries the labels, and the package's own integrator for lines of force
  knotted  MODEL    the tube round a core line: a rotation by 2πν a circuit, ν = √(2pq(p + q)); the profile ι(h);
                    the shear at the core line in closed form; helicity and energy of the (p, q) fields
  hopfion  MODEL    the squared hopfion's lines of force after t = 0 as an integrable system: the Hamiltonian, the
                    second invariant I₂, the winding W as a ratio of two periods, the two closed lines
  tetrad   MODEL    a tetrad e = dx + B for the squared hopfion that is frozen into its null flow; its field
                    strength; the helicities (2/3)π² and π²

What this is, and is not.
  * Every field here is a weak wave on flat Minkowski space: Maxwell's equations and the linearised vacuum Einstein
    equations.  Nothing is a solution of the full Einstein equations.
  * Newton's constant is not derived anywhere in MTFT (the gravity map's negative N-G1), and nothing here bears on it.
  * No link to X₀(143) is claimed.  These are tools for the gravity map's pieces III and IV.

Whose work is whose (a search of the literature, 16 November 2026).
  * The fields F_pq, their invariant Re(α^p β^q) and its nested tori: Kedia, Bialynicki-Birula, Peralta-Salas and Irvine,
    Phys. Rev. Lett. 111, 150404 (2013).  The hopfion itself is Rañada's.
  * The gravitational hopfion and its tendex and vortex lines as the light's lines turned by a local duality angle:
    Swearngin, Thompson, Wickes, Dalhuisen and Bouwmeester, arXiv:1302.1431; Thompson, Swearngin and Bouwmeester,
    J. Phys. A 47, 355205 (2014); Thompson, Wickes, Swearngin and Bouwmeester, J. Phys. A 48, 205202 (2015).
  * The integrable quartic potential that ``hopfion`` lands on is Garnier's system, classical; it separates in elliptic
    coordinates.
  * Frozen-in gravitational fields and their helicity e_a ∧ de^a (the Nieh–Yan 3-form of a teleparallel frame): Asenjo,
    Winkler and Comisso, arXiv:2609.25240.  A helicity of the gravitational hopfion in transverse-traceless gauge,
    which changes in time: Huang, Long, Zhao and Zhou, arXiv:2608.30154, section 6.4.
  * Derived in the lab's records MP09 and MP21 to MP23 (14 November 2026), tested there and here, unreviewed, and not
    found in those sources: the rotation number ν, the profile and the shear; the reduction of the later lines of
    force to the integrable system, with I₂ and the winding; the frozen tetrad and its helicity (2/3)π².  Not found
    is not the same as new.

Open: T-G3 (a three-dimensional resistive example) has not been attempted; H-G7 (``knotted.OPEN``) is untested.
"""
from __future__ import annotations

from . import fields, hopfion, knotted, tetrad

OPEN = {"T-G3": "a three-dimensional example with resistivity: dH/dt = -2 int E.B to tolerance, and the linking lost per unit of the rate; not attempted",
        "H-G7": knotted.OPEN["H-G7"]}
NOT_DERIVED = {"N-G1": "Newton's constant is not derived"}

__all__ = ["NOT_DERIVED", "OPEN", "fields", "hopfion", "knotted", "tetrad"]
