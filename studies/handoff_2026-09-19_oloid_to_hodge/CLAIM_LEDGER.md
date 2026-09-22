# Claim and evidence ledger

This is a navigation aid to the preserved reports, not a replacement for their hypotheses and proofs. EXACT denotes a stated algebraic identity or theorem consequence; COMPUTED denotes finite calculations with their stated certification conditions; DIAGNOSTIC denotes numerical evidence; OPEN denotes a construction or identification not established here.

Study paths are beneath `studies/`:

- **G**: `MTFT_Graded_Moonshine_Investigation/INVESTIGATION.md`
- **T**: `MTFT_Triangle_Feigenbaum_Hodge_Investigation/REPORT.md`
- **H**: `MTFT_Hodge_Dynamics_Transport/REPORT.md`
- **O**: `notes/OLOID_RETURN.md` at the handoff root.

## Graded geometry and the Monster

| ID | Result | Evidence and boundary |
|---|---|---|
| G01 | For X0(143), genus 13 and canonical Hilbert series (1+11t+11t²+t³)/(1−t)². | EXACT from the stated curve data and Riemann–Roch; G. The grading is tensor degree, not time or Fourier exponent. |
| G02 | Canonical product dimensions through degrees 1–5 are 13,36,60,84,108, with exact AL sectors. | COMPUTED over two finite fields, plus matching geometric upper bounds and Sturm precision. Conditional on the frozen series being the stated forms with the stated labels; G and `graded_results.json`. |
| G03 | Canonical, log-canonical, and ordinary cusp-form spaces differ. | EXACT; for example dim R2=36 while dim S4=40. See G. |
| G04 | No nontrivial grading-preserving Monster action on either of the two specified generated section rings. | EXACT obstruction using generator dimensions; independent geometric argument in G. It does not exclude every construction using additional states. |
| G05 | Substitution by T2 and the Leibniz lift of T2 both fail to preserve the canonical ideal. | EXACT nonzero q-coefficient witnesses; all 55 shipped basis relations fail for each proposed lift. G and `operator_results.json`. |
| G06 | The intrinsic bracket {f(dz)^n,g(dz)^m}=(nfg′−mf′g)(dz)^(n+m+1) preserves the canonical ring. | EXACT coordinate/cusp arguments and Poisson interpretation; finite checks support implementation. Ordinary Poisson bracket, not automatically a super-bracket. |
| G07 | The first bracket map wedge²R1 -> R3 has rank 47, kernel 31, cokernel 13. | COMPUTED exact finite-field lower bounds meet sector upper bounds; same fixture assumptions as G02. G and `bracket_results.json`. No particle count follows. |
| G08 | A counting ensemble with number operator N has partition function H_R(exp(−β)). | EXACT for the defined ensemble. A physical energy scale is additional input; G. |

## Triangle symmetry, iteration, and Hodge structures

| ID | Result | Evidence and boundary |
|---|---|---|
| T01 | T(n)=n(n+1)/2 is invariant under n -> −1−n; (2n+1)²=8T(n)+1. | EXACT, present in reviewed source; T §1. The separate three-square forbidden-indicator scaling is a different statement. |
| T02 | F_r(n)=2rT(n) is conjugate to the logistic family; z=r(n+1/2) gives z -> z²+c with c=r(2−r)/4. | EXACT after choosing the evolution and r; T §2. The unscaled T alone does not select chaos. |
| T03 | C_(m,c): y²=f_c composed m times(x) has generic genus 2^(m−1)−1 and degree-two tower maps. | EXACT on smooth fibers; T §3. Iteration index m is not period-doubling generation. |
| T04 | Rational H1 splits into inherited and new Prym pieces; the new part has an explicit quotient-curve model. | EXACT with standard double-cover theory; T §3. Isogeny/rational splitting is not an integral product identity or automatic principal polarization. |
| T05 | Critical orbit returns control the discriminant; a primitive period L first gives a node in C_L. | EXACT resultant and local arguments; T §4. At the Feigenbaum accumulation parameter every finite tower member remains smooth. |
| T06 | C2 and its free-involution quotient are degree-two isogenous, yielding the standard X0(2) pair of j-functions. | EXACT explicit quotient; T §5. No identification with X0(143). |
| T07 | The nodal branch integral has coefficient 1/(2√2) multiplying log(1/ε). | Exact asymptotic from the normalized elliptic integral, with high-precision DIAGNOSTIC checks; T §6. |
| T08 | The principal superstable cascade through period 1024 gives δ10≈4.66919515603 and abs(α10)≈2.50290757151. | DIAGNOSTIC, 60/90-digit calculations; T §7. No new universality theorem or physical-constant derivation. |

## Transport and proposed dynamics

| ID | Result | Evidence and boundary |
|---|---|---|
| H01 | The elliptic 2×2 period connection is given explicitly and verified by polynomial identities. | EXACT; H §1, `root_checks/elliptic_exact.py`. |
| H02 | Compact nonsingular interval transport has zero asymptotic Lyapunov exponents, even with chaotic base motion. | EXACT telescoping/boundedness proof under the stated path and norm hypotheses; H §2. It is not a no-go theorem for all Hodge dynamics. |
| H03 | Specified loop words give bounded, linear, and exponential growth; U²S has spectral radius 2+√3. | EXACT trace calculation with DIAGNOSTIC ODE support; H §3. The chosen periodic loop protocol does not require chaos or select a physical rate. |
| H04 | The genus-three connection is exactly A_E direct-sum A_D for inherited elliptic and new genus-two sectors. | EXACT direct reduction and flat projectors; H §4 and `geometry/TRANSPORT_NOTE.md`. |
| H05 | Around the primitive period-three degeneration, inherited monodromy is identity and the new part is rank-one unipotent. | EXACT local geometry; DIAGNOSTIC matrix integration. H §4 and `monodromy/NEW_COMPONENT_NOTE.md`. Integral normalizations across isogenies need care. |
| H06 | Pullback by a finite degree-d holomorphic map multiplies Hodge norm squared by d. | EXACT for the stated H1 setting; H §5. Degree-normalized pullback is isometric on its image. |
| H07 | The explicit Fricke lift intertwines the connection; the unnormalized round trip is 2I with compatible branches. | EXACT formulas and rational checks; H §5. Other square-root choices can add a central minus sign. No expanding normalized round trip. |
| H08 | Smooth period samples can inherit δ from preselected Feigenbaum parameters. | EXACT Taylor argument under the derivative hypotheses; H §6. Observing that ratio is not an independent derivation. |

## Return to the oloid

| ID | Result | Evidence and boundary |
|---|---|---|
| O01 | Oloid contact angles satisfy cos(u)=−cos(t)/(1+cos(t)). | Published geometric result, equation (6) in Dirnböck–Stachel; O. |
| O02 | The coordinate n=−1/(q+2) conjugates w(q)=−q/(1+q) to n -> −1−n. | EXACT rational substitution, derived in the final reflection; O. Real parameter statement, not integer-lattice identification. |
| O03 | The same rational formula occurs as the elliptic Fricke parameter involution. | EXACT comparison with T06/H07. No lift from physical oloid motion to the curve family has been supplied. |
| O04 | Ideal no-slip oloid constraints are integrable. | Published mechanical model in Kuleshov et al.; O. Holonomic constraints and transport holonomy are distinct notions. |
| O05 | Two solid oloids intersected in R3 cannot form S4. | Dimensional obstruction. If the intersection has interior, convexity gives a 3-ball with S2 boundary. A configuration space is a separate object. |

## Open or unsupported identifications

The studies do not establish any of the following:

1. Selection of a logistic parameter, a renormalization operator, or physical time by the reviewed MTFT data.
2. Feigenbaum δ as a local nodal-monodromy eigenvalue, or as expansion from degree-normalized holomorphic pullback.
3. A new exact identity relating the physical fine-structure constant to the Feigenbaum constants. The previously discussed numerical expressions remain numerical observations/conjectures.
4. A Monster vertex algebra or physical Monster symmetry derived from X0(143), 13!×11!, or the canonical ring studied here.
5. An identification of the quadratic tower's X0(2) correspondence with X0(143).
6. A four-dimensional gravitational action, spacetime metric, Standard Model, or quantum-gravity theory from these calculations.
7. Chaotic oloid motion or a Neimark–Sacker mechanism inferred merely from the order-two contact exchange.

The distinctions and negative results are part of the handoff, not optional qualifications to discard during continuation.
