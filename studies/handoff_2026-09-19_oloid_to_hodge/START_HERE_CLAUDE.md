# MTFT: from the oloid to Hodge transport

Consolidated handoff for Claude · Roger Tano · 19 September 2026

## Purpose and scope

Continue the investigation begun with the oloid and developed through triangular-number symmetry, quadratic dynamics, modular geometry, graded rings, and Hodge transport. This package collects all three completed study bundles available from this conversation, their programs and recorded outputs, the reviewed MTFT v0.32.0 source distribution, the supplied PDFs and images, and the final oloid follow-up.

The earlier 69-page document mentioned at the start of the conversation is outside this package: its contents were not among the available inputs to this detour. This is not an export of every MTFT conversation or the live GitHub repository.

## Read in this order

1. `CLAIM_LEDGER.md` — the main results, their evidence, and their limits.
2. `notes/OLOID_RETURN.md` — the latest connection, added after the three original bundles.
3. `studies/MTFT_Graded_Moonshine_Investigation/INVESTIGATION.md` — the canonical ring, exact Atkin–Lehner sectors, operator obstructions, and compatible bracket.
4. `studies/MTFT_Triangle_Feigenbaum_Hodge_Investigation/REPORT.md` — triangle symmetry, the chosen dynamical family, curve tower, critical degenerations, and level-two correspondence.
5. `studies/MTFT_Hodge_Dynamics_Transport/REPORT.md` — explicit connections, monodromy, and transport-growth obstructions.
6. `REPRODUCE.md` — commands and dependencies. Each study also retains its original README and manifest.

`SESSION_CONTEXT.md` records the broader questions and side investigations. It is a curated summary, not a verbatim conversation transcript. `notes/INTERMEDIATE_DERIVATION.md` preserves an earlier working derivation; the finished reports take precedence.

## The current mathematical position

For the compact modular curve X0(143), the first study computes a canonical graded ring and symmetry-resolved product dimensions. A direct grading-preserving Monster action on that particular ring is obstructed. Two proposed Hecke lifts fail to preserve its relations. An intrinsic Poisson/Rankin–Cohen bracket does survive; its first bracket map has rank 47 and cokernel dimension 13, under the frozen modular-form identifications.

Separately, the triangle reflection can be carried into a chosen quadratic iteration and an explicit tower of algebraic curves. The tower splits inherited and new Hodge components. Its first elliptic member has a degree-two isogeny described by X0(2). That modular curve has not been identified with the project's X0(143).

The transport continuation constructs exact elliptic and genus-two period connections. Parameter loops can produce nontrivial monodromy and, for specified repeated loop words, exponential linear transport. Chaotic motion confined to a compact nonsingular real interval gives zero asymptotic growth for the ordinary flat transport considered here. Degree-normalized holomorphic pullback is isometric. These are constraints on proposed mechanisms, not failed numerical searches.

The final return to the oloid identifies the same rational involution in two settings:

    w(q) = -q/(1+q),       w(w(q)) = q.

It relates oloid contact parameters in a published geometric analysis and is the Fricke parameter transformation of our elliptic family. The change n = -1/(q+2) conjugates it to the triangle reflection n -> -1-n. This establishes a shared algebraic operation. No map has yet been constructed that carries the oloid's rolling dynamics to the Hodge transport system.

## Keep these distinctions explicit

- An invariant does not select an evolution law. The choice to iterate F_r(n)=r n(n+1), and the choice of r, are additional inputs.
- Exact finite-field ranks need the matching geometric upper bounds and the stated upstream fixture assumptions to imply rational ranks.
- High-precision roots and floating-point ODE diagnostics are not interval certificates.
- A parameter-space connection and its gauge covariance do not yet specify a physical spacetime gauge theory.
- A shared involution, a numerical resemblance, and a physical identification are different levels of evidence.
- No physical fine-structure constant, gravitational action, Monster vertex algebra, or Feigenbaum renormalization operator on the Hodge tower was derived here.

## Suggested next audit

First inspect the explicit formulas, frozen inputs, and recorded outputs before expanding the claims. The concrete open target is a selected renormalization correspondence acting compatibly on parameters, curves, cohomology, and norms. It must survive the bounded-flat-transport and degree-normalization obstructions. For the oloid connection, specify the intended configuration space and an actual lift before assigning Hodge or physical meaning to the common parameter involution.

The compact canonical-ring bracket is another available direction, but a real phase space, a Hamiltonian, an energy scale, and a physical time would still need to be specified.

## Packaging provenance

The three original ZIPs are retained byte-for-byte in `original_archives/`; their contents are expanded in `studies/` for immediate reading. All 47 payload checksums in their original manifests were verified. Their standalone report copies matched the versions inside the archives. The source tarball has SHA-256:

    46ffc563a0cd81a7a82d975497bf64d82f8471eb424409e3d154532b5d8d6069

The consolidation does not modify MTFT or rerun its full test suite. Existing study calculations were preserved rather than rerun during packaging. Run `python3 checks/verify_bundle.py` before reproduction to verify the delivered file checksums. The attached infographics and Wikipedia PDFs are reference inputs; their presence is not an endorsement of every statement they contain.
