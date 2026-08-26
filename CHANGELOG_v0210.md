# mtft v0.21.0 — periods, true Hodge geometry, and the X0(143) physics bridge

This release promotes the long-standing `mtft.periods` debt into the wheel.  It
is deliberately a **physics bridge without a particle-fit layer**: the package
now computes the genuine genus-13 period/Hodge geometry of `X0(143)` on the
same 26-dimensional stage used by `mtft.hecke`, then exposes mathematical
charge and mixing diagnostics.  Turning those dimensionless objects into
masses or couplings still requires an explicit field theory and scale.

## New subpackage: `mtft.periods`

### `periods.core`
- Frozen v6 `Omega_cusp`, intersection matrices `Q,E`, exact symplectic change
  `S`, and the frozen `tau` record ship in the wheel.
- `omega_symplectic()` is **reconstructed at call time** as `Omega_cusp @ S`.
- `riemann_matrix()` recomputes `tau=A^-1 B` with mpmath; the v6 record replays
  to ~1e-49 at 60 dps.
- `hodge_complex_structure()` is pinned by `[I tau] J = i[I tau]` and
  `J^2=-I`; `hodge_metric()` returns the positive polarization metric `E J`.
- `charge_energy(n,m)` exposes the positive dimensionless quadratic form.  It
  is mathematics; a physical mass interpretation is explicitly PHENO unless a
  separate model supplies a scale.

### Append-only finding: v6 derived-field serialization defect
The old `x0143_period_data_v6.json` field `Omega_symplectic_13x26` contains a
GP-scientific-notation tokenization defect (e.g. a tiny `...E-59` exponent was
lost after a whitespace split), producing O(10) corrupt entries.  **Primitive
`Omega_cusp`, exact `S`, and frozen `tau` are not affected**: `Omega_cusp @ S`
reproduces `tau` to ~1e-49 and the first Riemann bilinear residual remains
~1e-50 at full precision.  The bad derived field stays shipped for provenance
but is never used by the API.  No CC number is assigned here.

### `periods.bridge` — exact basis blocker closed
The v6 period computation and `mtft.hecke` used different Manin bases.  The
bridge reconstructs the 29 period paths by continued fractions inside the
live promoted Manin model.  Two exact changes of basis result:
- `relative_basis_change()` has det +1;
- `cuspidal_basis_change()` has det +1 and satisfies `B_Hecke=B_period C`;
- `hecke_to_symplectic_change()` = `S^-1 C`, also integral unimodular.

This is the missing basis identification from the M8 program.  Analytically,
`diag(a_p) Omega_H = Omega_H T_p` replays at p=2,3 below 1e-45.

### `periods.forms` — no-PARI q-series and Bergman density
- Reconstructs all 13 raw weight-2 q-series through q^140 in the v6 row order:
  143a1, four f2 real embeddings, six f3 real embeddings, 11a(q), 11a(q^13).
- f2/f3 use the exact power-basis `a_n` shipped by `mtft.codifferent`; elliptic
  rows use exact point counts + Hecke recurrences.
- The reconstructed rows span the frozen canonical q-expansion space with
  max absolute replay residual ~1.4e-12 in the float E2 gate.
- `bergman_density(z)` evaluates the actual canonical analytic density from
  alpha-normalized differentials and `(Im tau)^-1`.
- `q_tail_bound()` supplies a conservative Deligne-style q^140 truncation
  bound; at y=1/sqrt(143) it is <3e-29.

### `periods.physics` — quantitative M8 with the true period-derived J
- `hodge_structure_hecke()` transports the genuine transcendental Hodge J to
  the promoted `mtft.hecke` basis.  It commutes with T_p for p=2,3,5,7,11,13
  and anticommutes with the star involution to ~1e-14.
- `graph_coupling()` promotes the M7 width/degree/distance graph channel using
  a frozen exact-rational harmonic embedding with call-time G and M rebuilds.
- `complex_linear_decomposition()` splits V relative to J into commuting and
  J-antilinear pieces.
- `cp_channel_report()` measures the J-antilinear fraction with the Hodge
  Hilbert-Schmidt metric.  For the canonical width potential:

      ||V_-||_H / ||V||_H = 0.06796942853668321
      power fraction       = 0.004619843215603286

  This is a **new mathematical diagnostic, not the historical M8b amplitude
  and not a Standard-Model CP observable**.  The structural result survives:
  `[V,J] != 0`, while the commutator is star-odd as predicted.
- `finite_charge_partition()` intentionally accepts an explicit finite charge
  set; it does not masquerade as the full 26-dimensional theta series.

## Gates and CLI
`python -m mtft.periods verify` runs six no-PARI gates: exact symplectic
integrality, period reconstruction, exact basis bridge, Hodge/Hecke
compatibility, native q-expansion span, and quantitative M8.

Other commands:
- `python -m mtft.periods tau`
- `python -m mtft.periods physics`
- `python -m mtft.periods bergman --x 0 --y 0.0836242`

## Packaging correction
`sympy>=1.10` is now declared in core dependencies.  v0.20's origami package
hard-imported SymPy while `pyproject.toml` omitted it.
