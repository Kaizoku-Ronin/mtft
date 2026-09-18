"""INT-10: immutable model records — M1, the rejected parent controls, and C3X — with the conventions each carries.

A record separates: global gauge group (and quotient), representations with multiplicities and 6D chiralities,
flux/cocharacter data, internal bundle and mode operator, the ORIGIN of the family multiplicity (internal zero modes vs
input copies), supplied couplings vs computed consequences, and the local/integral/global consistency status.
Records are data; gates live in `research.pipeline` and return witnesses.  C3X is a research candidate identifier,
not a replacement for M1 (Gate 4, route 2: one parent x three internal modes is the MTFT-native target)."""
from types import MappingProxyType as _frozen

def _rec(**kw): return _frozen(kw)

M1 = _rec(model_id="M1", kind="five-stack line-flux model on X0(143)", global_group="U(3)_c x U(2)_L x U(1)_a x U(1)_b x U(1)_d (direct product; independent line fluxes)",
    stacks=("c", "L", "a", "b", "d"), ranks=(3, 2, 1, 1, 1), degrees=(0, -3, 3, 3, 0),
    hypercharge=("1/6", "0", "-1/2", "1/2", "-1/2"), B_minus_L=("1/3", "0", "0", "0", "-1"),
    family_bundle="S0 (x) O(P1+P2+P3), degree 15, h^0 = 3, h^1 = 0 (purity: both signs of u on the divisor)",
    higgs_space="H^1(O(-2 sum P)) = H^0(K(2 sum P))^*, dimension 18 (Dolbeault classes, not elementary-scalar zero modes)",
    multiplicity_origin="internal zero modes (three sections of one bundle)", yukawas="section-product tensors; 6D interaction not specified (same-chirality bilinear obstruction, AXG-03)",
    six_d_chiralities="not assigned in the specification", supersymmetry="not declared", status=_rec(local_4d_anomaly="Y, B-L, phase anomaly-free (exact)", six_d_lift="obstructed (AXG-02/03)", vacuum="none", scale="none"))

ONE_PARENT_INDEX3_CONTROL = _rec(model_id="ONE_PARENT_INDEX3", kind="single parent with index-three flux, six sectors (AXG-04 §9)",
    result="primitive tensor-charge integrality witness: norm -8/3 (a charged repair still gives -4/3); constant-Higgs limit gives degenerate singular values",
    model_gate="FAIL", source="AXG-04 larger_parent_notes.md")

U8_DIRAC_ADJOINT_CONTROL = _rec(model_id="U8_ADJOINT", kind="U(8) Dirac-adjoint parent (AXG-04 §9)", flux_centralizer="U(4) x U(2) x U(2)",
    result="mirror content; split-bundle Yang–Mills Hessian unstable in the fixed-metric test: complex Morse index 312 (180 colored)",
    model_gate="FAIL (fixed-metric test; not a coupled-AdS stability result)", source="AXG-04 larger_parent_audit.py")

C3X = _rec(model_id="C3X", kind="explicit charged 6D EFT candidate (AXG-04)", global_group="SU(3) x SU(2) x U(1)_h x U(1)_X (direct product, h = 6Y)",
    fields=(("Q", 3, 2, 1, 1, "+"), ("U", 3, 1, 4, 1, "-"), ("D", 3, 1, -2, 1, "-"), ("L", 1, 2, -3, 1, "+"), ("E", 1, 1, -6, 1, "-"), ("N", 1, 1, 0, 1, "-"), ("H", 1, 2, 3, 0, "scalar")),
    copies=3, flux="one unit of X flux: L_X = O(P+Q-R) with u(P) != u(Q); (h^0, h^1)(S0 L_X) = (1, 0): one pure mode per copy",
    multiplicity_origin="three INPUT copies x one internal mode", yukawas="four free 3x3 parent matrices (opposite 6D chiralities allow the scalar bilinear)",
    anomaly="I8 = (W + 9 h^2)(3 C + W + eta - 27 h^2 - 6 x^2), integral factors; Omega = [[0,1],[1,0]]", flux_reduction="I6 = -12 x (W + 9 h^2); BF level -12; candidate Z12",
    bordism="Omega_7^Spin(BG) = 0 (ordinary category); global GS construction OPEN", background="unwarped product: AdS, no scale separation (classical control)",
    status=_rec(local_factorization="PASS", integral_lattice="PASS", bordism="PASS (necessary, not sufficient)", global_GS="OPEN", vacuum="AdS control only", flavor="input", scale="none"))

RECORDS = _frozen({"M1": M1, "ONE_PARENT_INDEX3": ONE_PARENT_INDEX3_CONTROL, "U8_ADJOINT": U8_DIRAC_ADJOINT_CONTROL, "C3X": C3X})
