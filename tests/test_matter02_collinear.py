"""Exact geometry, chirality conventions and global cup-product witnesses."""
from fractions import Fraction as F
from dataclasses import asdict
import json
import pytest
from mtft.matter.collinear import (
    cohomology,cusp_jet_certificate,yukawa_witness,report,_det,_rref,
)


@pytest.mark.parametrize("q",range(-15,16))
def test_rr_matches_archived_independent_intersection(q):
    from mtft.surface.blowup import chi as archived_chi
    row=cohomology(q)
    assert row.chi==archived_chi((q,0,-3*q))
    assert row.h0-row.h1+row.h2==row.chi
    assert min(row.h0,row.h1,row.h2)>=0


@pytest.mark.parametrize("q",[True,False,1.0,"1",None,F(1),[1]])
def test_exponents_strict_integer(q):
    with pytest.raises(TypeError):
        cohomology(q)


def test_exact_boundary_dimensions_and_zero_gauge_exception():
    assert (cohomology(-1).h0,cohomology(-1).h1,cohomology(-1).h2)==(0,15,12)
    assert (cohomology(1).h0,cohomology(1).h1,cohomology(1).h2)==(0,18,12)
    assert (cohomology(0).h0,cohomology(0).h1,cohomology(0).h2)==(1,14,13)
    assert cohomology(0).holo_chiral_dim is None
    assert cohomology(12).h2==1 and cohomology(13).h2==0
    assert cohomology(-6).h2==2 and cohomology(-7).h2==0
    assert cohomology(2).phi_dim==10


def test_conjugation_and_both_named_chirality_conventions():
    for q in range(1,20):
        forward,reverse=cohomology(q),cohomology(-q)
        assert forward.phi_dim==reverse.h2
        assert forward.holo_chiral_dim-reverse.holo_chiral_dim==forward.bhv_net==3*q
        assert forward.chi-reverse.chi==forward.archived_net==-3*q
        assert reverse.archived_net==-forward.archived_net
        assert reverse.bhv_net==-forward.bhv_net


def test_primary_and_adapted_cusp_basis_and_independent_matrix_change():
    from mtft.canonical import s2_qexpansions,adapted_qexpansions,adapted_basis
    original=s2_qexpansions(); adapted=adapted_qexpansions(); basis=adapted_basis()
    assert [[sum(row[i]*basis[i][j] for i in range(13)) for j in range(13)]
            for row in original]==adapted
    cert=cusp_jet_certificate()
    assert cert["passed"] and cert["gap_sequence"]==list(range(1,14))
    assert cert["full13_jet_determinant"]=="-18196056815164784640"
    assert cert["basis_change_determinant"]=="-1078272"
    assert all(row["rank_original"]==row["rank_adapted"]==row["vanishing_order"] for row in cert["jet_rows"])


def test_actual_canonical_combination_is_order3_and_a5_zero():
    from mtft.canonical import s2_qexpansions
    witness=yukawa_witness(); data=s2_qexpansions()
    coefficients=tuple(sum(F(row[i])*witness.canonical_basis_coefficients[i] for i in range(13)) for row in data)
    assert coefficients==witness.canonical_q_coefficients
    assert coefficients[1:6]==(0,0,0,1,0)
    assert dict(witness.phi_coefficients)[5]==1
    assert 6 not in dict(witness.phi_coefficients)


def test_independent_residue_by_formal_polynomial_multiplication():
    witness=yukawa_witness()
    phi=dict(witness.phi_coefficients)
    # Build each Laurent polynomial product before reading its z^-1 w^-1 term.
    residue=[]
    for a,b in witness.positive_jets:
        row=[]
        for i,j in witness.negative_principal_parts:
            product={(n+a-i,b-j):value for n,value in phi.items()}
            row.append(product.get((-1,-1),F(0)))
        residue.append(tuple(row))
    assert tuple(residue)==witness.residue_matrix==((1,0,0),(0,0,1),(0,1,0))
    assert _det(residue)==witness.determinant==-1
    assert witness.full_matrix_shape==(18,44)
    assert witness.higgs_space_dimension==10


def test_global_lift_and_target_jet_conditions():
    witness=yukawa_witness()
    assert all(a+b<3 and b>0 for a,b in witness.positive_jets)
    assert all(i+j<=9 and i>=1 and j>=2 for i,j in witness.negative_principal_parts)
    target=[(5+a,b) for a,b in witness.positive_jets]
    assert len(set(target))==3 and all(a+b<8 and b>0 for a,b in target)
    # Base global sections depend only on z, so these distinct w-dependent
    # target monomials survive its image. j>=2 principal parts have zero
    # residue against every w-independent base dual section.
    for i,j in witness.negative_principal_parts:
        base_dual={(n,0):F(n*n+1,3) for n in range(-3,25)}
        product={(n-i,w-j):value for (n,w),value in base_dual.items()}
        assert product.get((-1,-1),F(0))==0


def test_nonzero_u8_root_bracket_with_independent_matrix_units():
    def unit(i,j):
        return [[int(r==i and c==j) for c in range(8)] for r in range(8)]
    def mul(a,b):
        return [[sum(a[i][k]*b[k][j] for k in range(8)) for j in range(8)] for i in range(8)]
    a,b,c=unit(5,3),unit(7,5),unit(3,7)
    ab,ba=mul(a,b),mul(b,a)
    bracket=[[ab[i][j]-ba[i][j] for j in range(8)] for i in range(8)]
    product=mul(bracket,c)
    assert sum(product[i][i] for i in range(8))==yukawa_witness().root_bracket_trace==-1


@pytest.mark.parametrize("basis",[
    [[1,2,0],[0,1,3],[0,0,1]],
    [[0,1,0],[1,0,0],[0,0,-2]],
    [[F(1,3),0,1],[0,2,0],[1,1,0]],
])
def test_positive_basis_changes_preserve_rank3(basis):
    matrix=yukawa_witness().residue_matrix
    transformed=[[sum(F(basis[i][k])*matrix[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    assert _det(basis)!=0
    assert len(_rref(transformed)[1])==3
    assert _det(transformed)==_det(basis)*_det(matrix)


def test_report_serialization_and_physical_limits():
    output=report()
    json.dumps(output)
    assert output["higgs_vev_selected"] is False
    assert output["physical_mass_matrix"] is False
    assert output["yukawa_witness"]["determinant"]=="-1"
    assert output["cohomology"][5]["holo_chiral_dim"] is None
