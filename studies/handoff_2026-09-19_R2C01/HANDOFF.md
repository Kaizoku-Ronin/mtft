# R2C-01 integration handoff

## Finding to register

Suggested identifier: R2C-01. Status: EXACT, CONDITIONAL ON THE DECLARED CLASS.

For the M1 stack group, fluxes and ten-sector field content, preserving the six
named three-family indices forces all six associated 6D Weyl chiralities to be
positive in the stated convention. All four local nonderivative elementary-
scalar Yukawa contractions require opposite signs within the corresponding
fermion pairs. The two requirements have no simultaneous assignment.

Exhaustion: 1,024 assignments; 16 retain the six family indices; four retain the
whole sector ledger; 64 allow four scalar bilinears; zero overlap. The output
also keeps the index-zero modes and the extra doublet sectors visible.

This extends the existing AXG-03 scalar-bilinear witness by exhaustively closing
the previously unassigned M1 chirality choice. It does not withdraw the pure
three-mode arithmetic construction or make the unassigned parent into an
inconsistent completed theory.

## Integration contract

An appropriate new research module could be named research.chirality. Suggested
operations and required inputs are below; they are interface proposals, not a
claim that these APIs have already been installed in MTFT.

| Operation | Required inputs | Output and witness |
|---|---|---|
| signed_index | Oriented representation, 6D chirality, internal index | Signed 4D index; preserve orientation explicitly |
| enumerate_chiralities | Sector list, permitted signs, target ledger, interaction list | Assignment table, separate gate outcomes and conflicting sector pair |
| bilinear_selection | Two Weyl chiralities, Clifford parity, contraction type | Vanishing/nonvanishing algebraic contraction; no overlap claim |
| spectrum_anomaly | Explicit representations and 6D multiplicities/signs | Complete fermionic I8, including p2 and color-cubic terms |
| flux_pushforward | Polynomial and integral flux vector | I6, checked against the signed 4D representation ledger |
| primitive_norm_witness | Global group, allowed cocharacter, cancellation assumptions | Exact quartic coefficient and required tensor norm |

A parent record should carry both:

1. The number of six-dimensional field copies.
2. The internal (h0,h1) for each field's actual bundle and operator.

Never implement family degeneracy by multiplying the bulk I8 by h0−h1.
Use the pushforward to compute the reduced anomaly. Changing chirality must
trigger recomputation of the 4D ledger; current unsigned-flux helper functions
cannot be reused silently with changed signs.

## Essential regression witnesses

- M1 signs all positive: original signed indices and anomaly ledger, but four
  elementary-scalar Lorentz contractions vanish.
- Flipping ca at fixed M1 flux: up-antiquark target index changes from −3 to +3.
  The scalar contraction with cL opens, but the intended family is reversed.
- Flipping cd at fixed M1 flux: net 4D chiral ledger unchanged, but bulk I8 and
  its gravitational/color-cubic coefficients change. Its explicit (2,2)
  internal kernel remains visible.
- Every sector's Chern-character I8 agrees with a direct sum of fourth powers
  of Cartan weights, including mixed gravitational and gauge terms.
- Every formal stack derivative of a sector I8 equals its oriented I6 times
  the corresponding fundamental charge.
- All same-chirality scalar contractions vanish for both Dirac and transpose
  contractions; odd-Clifford/vector contractions have rank four.
- Changed-flux family equations have exactly four branches after fixing the
  common flux shift. Their Higgs degrees become zero; the original 18-mode
  Higgs data must not be inherited.
- The six-sector changed-flux controls require primitive norms ±8/3 or ±4/3.
- Of the 64 ten-sector changed-flux completions, eight cancel the fermionic
  pure gravitational term, but all eight fail the displayed integral norm.
- C3X's factorization is recovered from a charge-spectrum calculation. A test
  that merely multiplies the two factors and compares with their own product
  does not establish this independent route.

Use the explicit logical derivations and direct-weight route as the expected
answers. Do not save opaque numeric outputs as unexplained truth.

## Scope of the additional anomaly result

The fixed-flux family-preserving set has signed bulk representation dimension
at least eight. Its fermionic p2 coefficient therefore never vanishes.
All 1,024 sign assignments also have a nonzero independent color-cubic
coefficient for the full five-stack gauge background.

This excludes cancellation by ordinary characteristic 4-form products alone.
It does not exclude new charged matter, a declared 0-form/6-form cancellation
mechanism, additional anomaly inflow or a different global gauge construction.
If tensors or a supergravity spectrum are added, include their own anomaly
polynomials before claiming a total gravitational cancellation.

The primitive norm tests assume an integral tensor-charge lattice and allowed
backgrounds along the stated cocharacter. They do not settle arbitrary axionic
or higher-group extensions. Record those assumptions in the gate result.

## Next experiment

Study the physical origin of the existing Dolbeault Higgs space while retaining
M1's flux and family bundles. The report provides four gauge- and Lorentz-
invariant vector interaction structures with equal 6D family chiralities.

A proposed parent must specify the vector/form representation, reality
conditions, kinetic action, interaction, bundle connection, and complete
anomaly-cancelling field content. Derive its internal quadratic operator and
verify whether its physical kernel is the desired Dolbeault space. Then compute
the induced overlap map and the spectrum of unwanted modes.

The odd-Clifford matrix check does not establish a healthy charged-vector
theory. Do not revive the rejected U(8) parent or identify a Dolbeault class
with a massless particle solely because the contraction is allowed.

Gravity, stabilization and absolute masses remain uncomputed by R2C-01.

## Provenance and review

The attached v0.32.0 archive is the actual baseline. Its SHA-256 differs from the
hash quoted in the conversation; both are recorded in PROVENANCE.json. This
handoff contains source snapshots and standalone calculations. It does not
modify or publish the MTFT release.

The checked run has 293 exact checks, zero failures. The release-wide test suite
was not rerun for this study.
