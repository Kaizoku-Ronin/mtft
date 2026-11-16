# studies/gravity_2026nov — the lab's records behind `mtft.gravity`

`mtft.gravity` (added 16 November 2026) is a port of four records of the Modular Plane lab. This folder keeps them as they were
written: the plan of each, with its predictions registered before the run and its result and faults after; the script; the numbers.

| folder | record | what it holds |
|---|---|---|
| `lab/mp09` | MP09, 10 November 2026 (T-G5) | the matrix of the field lines round a core line: a rotation by 2πν, ν = √(2pq(p + q)); the profile ι(h). `mp81_verify.py` needs `mp8_verify_part1.py` and `mp8_numbers.json` beside it |
| `lab/mp21` | MP21, 14 November 2026 (T-G7) | the shear at the core line in closed form, by three routes |
| `lab/mp22` | MP22, 14 November 2026 (T-G8) | the squared hopfion's lines of force after t = 0 as an integrable system: H, I₂, the winding, the two closed lines |
| `lab/mp23` | MP23, 14 November 2026 (T-G1, T-G2) | a tetrad frozen into the null flow, and the helicities (2/3)π² and π² |

The package's tests read `mp81_numbers.json`, `mp21_numbers.json`, `mp22_numbers.json` and `mp23_numbers.json` from here and compare the
package with them.

## How the package differs from these scripts

- **No scipy.** The scripts integrate every line with scipy, which is an optional extra of the package. The package integrates with a
  routine of its own (`mtft.gravity.fields.integrate`, an extrapolation of the modified midpoint rule), so each agreement with these
  numbers is also an agreement between two integrators. Along the lines of force the second invariant is kept to 10⁻¹³ of its range
  here, against 2.7·10⁻¹⁰ in MP22.
- **Exact derivatives.** MP23 took the tetrad's derivatives by nested differences at 40 digits. The package carries every quantity with
  its first and second derivatives, so the check runs in double precision, at any number of digits, or over the Gaussian rationals,
  where the differences are exactly zero.
- **A third route to the shear.** MP21's third route was the lab page's own function, which is not shipped. The package's is the period
  of the cross-section's Hamiltonian system, integrated as a differential equation.
- **Not ported.** MP09's return map from the variational equation in space; MP21's lines in space; MP22's re-reading of the older
  records' numbers (P5 a and c) and its torus of closed lines followed in space (P7; the package checks only that W = 65/64 there);
  MP23's search for the tetrad by least squares and its kernel run. They stay results of these scripts.

## Two things this build found in the records

- MP09 writes its line field on the label sphere "up to a positive factor". With B = Im F_pq the factor is negative. Nothing depends on it.
- `mp23_numbers.json` holds the two helicities on the label sphere with a minus sign, and the script's test took absolute values. The
  order of coordinates it integrates in, (η, ξ₁, ξ₂), is the orientation opposite to the one space gives the sphere. The package uses
  (ξ₁, η, ξ₂), and its sphere and its slices agree in sign.

The search of the literature that the release plan asked for is in `docs/gravity/KNOTTED_FIELDS.md`, section 1.
