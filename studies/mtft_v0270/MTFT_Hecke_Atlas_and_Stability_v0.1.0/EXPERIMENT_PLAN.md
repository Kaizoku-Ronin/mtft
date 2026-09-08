# Hecke atlas and sector stability — bounded experiment plan

2026-09-06. This plan is fixed after inspecting the defining algebra and numerical metric gates, before counting splitting types or simulating perturbed flows. This is an exploratory computational study, not a preregistered statistical trial.

## Inputs

MTFT 0.26.0's frozen X0(143) data, verified against its live cycle basis, supply integral T2, T3, W11, W13, the intersection form, and the numerical period-derived complex structure. The four rational Hecke sectors have real dimensions 2, 4, 8, 12. The quartic and sextic polynomials are the factors of T2, with the package's sign convention. Reuse the preceding exact prime/irregularity census for all primes p < 10201. Store compact copies of the inputs and hashes.

## A. Splitting atlas

Factor both polynomials modulo every prime in the census. Exclude primes dividing either polynomial discriminant from the joint unramified atlas; record ramification and level primes 11,13 separately. Verify each factor product and independently recover all unramified factor degrees through gcd(x^(p^k)-x,f). Compute exact Galois groups, and supply a second certificate using irreducible, (1,n-1), and transposition factor patterns. Establish the intersection of the two splitting fields algebraically, not by a frequency correlation. Enumerate all partition types and their exact permutation densities. Tabulate all 77 joint types, including zeros in this sample, and their product densities when the splitting fields are disjoint.

Join square-interval position, irregularity index, and base-2 Wieferich status to each prime. These are descriptive overlays. Do not infer infinitude of intersections with irregular primes or coverage of every square interval. Do not fit or select a statistical dependence claim from this finite atlas; asymptotic density predictions and finite counts are separate outputs.

## B. Sector stability

Construct rational CRT projectors from the squarefree minimal polynomial of T2. Verify idempotence, orthogonality, completeness, ranks, and exact Atkin–Lehner intersections. Work in coordinates orthonormal for the frozen Hodge metric G. Use the real 26-dimensional, J-commuting, symmetric energy model

    A0 = I + T2 / (4 ||T2||_op), A = A0 + epsilon V, F = J A, U(t) = exp(t F).

All norms refer to these orthonormal coordinates. Each V is symmetric, commutes with J, and has Frobenius norm 1. Thus ||V||_op <= 1 and A is positive for the planned epsilon values. This restricts the experiment to conservative linear mixing. It is a dimensionless model, with chosen parameters and no claim of physical frequencies or particle masses.

Use 16 paired random seeds 2026090600 through 2026090615; epsilon in [0, .05, .1, .2, .4]; and 161 equally spaced times in [0,80]. From the same Gaussian symmetric J-commuting matrix per seed construct five classes: (i) scalar offsets on the four Hecke sectors, (ii) arbitrary within-sector perturbations, (iii) the W11/W13 group average, (iv) only the oldspace–quartic cross blocks of that group average, and (v) unrestricted perturbations within the stated J-commuting symmetric class. Class (iv) has every diagonal Hecke block zero, so every block trace is balanced, while allowing the shared-character coupling.

Record commutators, off-block norms, positivity, frequency shifts, and transfer M[j,i](t) = ||Pj U(t) Pi||_F^2 / dim(i). Because U is orthogonal in this restricted model, columns sum to one. Report mean peak leakage over the fixed time grid, and time-averaged transfer at epsilon .4. Seed variation describes chosen perturbations, not uncertainty about an arithmetic theorem. Check representative flows independently using scipy.linalg.expm. Numerical tolerance 1e-9 for conservation/selection rules, with actual residuals retained; exact matrix identities use rational arithmetic.

## Image connections

Prove by induction that the matrix iteration Z_next = Z^2 + C preserves every Hecke block whenever Z0 and C do; if they belong to the same commutative Hecke algebra, every iterate belongs to it. Give an exact two-dimensional counterexample showing that zero diagonal block traces do not imply absence of mixing. The ham-sandwich theorem motivates a future spatial balance experiment; no spatial cut is being computed here and a balance identity is not a commutation identity.

## Deliverables and stopping point

Scripts, exact certificates, atlas, simulation arrays, validation, a report, two standalone figures, and an inspector command. Finish after the planned census, 400 flows (including baseline repeats), and meaningful verification. No expansion of the prime range or perturbation family based on an attractive pattern.
