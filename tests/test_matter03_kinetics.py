"""Field-coordinate normalization, same-vacuum covariance and invalid inputs."""
import math

import numpy as np
import pytest

from mtft.matter.kinetics import (
    canonical_frame, canonical_yukawa, canonical_mass_matrix, mass_singular_values,
)


def example():
    # Known upper-triangular canonical frames: all metrics are supplied INPUT.
    frames = (
        np.array([[1,(1+1j)/3,0],[0,1/2,1j/4],[0,0,2/3]],dtype=complex),
        np.array([[1/2,-1j/3,1/5],[0,2/3,(1-1j)/4],[0,0,3/2]],dtype=complex),
        np.array([[2/3,1j/5],[0,3/4]],dtype=complex),
    )
    inverse = tuple(np.linalg.inv(c) for c in frames)
    grams = tuple(b.conj().T@b for b in inverse)
    tensor = np.empty((3,3,2),dtype=complex)
    tensor[:,:,0] = [[1,0,0],[0,0,1],[0,1,0]]
    tensor[:,:,1] = [[.2+.1j,0,.3j],[.4,.5j,0],[0,-.2j,.1]]
    direction = np.array([1+.2j,-.3+.5j]); direction /= np.linalg.norm(direction)
    return tensor,grams,frames,direction


def test_known_complex_frame_and_holomorphic_component():
    tensor,grams,frames,_ = example()
    for gram,known in zip(grams,frames):
        actual=canonical_frame(gram)
        assert np.allclose(actual,known,atol=2e-14)
        assert np.allclose(actual.conj().T@gram@actual,np.eye(len(gram)),atol=2e-14)
    yn=canonical_yukawa(tensor,*grams,kappa=2j)
    # Higgs canonical coordinate 0 pulls back only to raw coordinate 0.
    expected=2j*(2/3)*frames[0].T@tensor[:,:,0]@frames[1]
    assert np.allclose(yn[:,:,0],expected,atol=2e-14)
    # det known exactly: bare normalized slice determinant -4/81.
    assert np.linalg.det(yn[:,:,0]) == pytest.approx((2j)**3*(-4/81),abs=2e-14)


@pytest.mark.parametrize("seed",range(8))
def test_all_field_basis_changes_keep_same_vacuum_masses(seed):
    tensor,grams,frames,n=example()
    matrix=canonical_mass_matrix(tensor,*grams,n,vev=3,kappa=.7+.2j)
    rng=np.random.default_rng(14303+seed)
    changes=tuple(2*np.eye(size)+.4*(rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))) for size in (3,3,2))
    grams2=tuple(b.conj().T@g@b for b,g in zip(changes,grams))
    tensor2=np.einsum("ai,bj,ck,abc->ijk",*changes,tensor)
    frames2=tuple(canonical_frame(g) for g in grams2)
    # Transport the raw expectation value directly, independently of Yukawa contraction.
    raw_h2=np.linalg.solve(changes[2],frames[2]@n)
    n2=np.linalg.solve(frames2[2],raw_h2)
    assert np.vdot(n2,n2).real == pytest.approx(1,abs=2e-13)
    matrix2=canonical_mass_matrix(tensor2,*grams2,n2,vev=3,kappa=.7+.2j)
    assert np.allclose(mass_singular_values(matrix2),mass_singular_values(matrix),rtol=4e-12,atol=2e-13)
    left_u=np.linalg.solve(frames[0],changes[0]@frames2[0])
    right_u=np.linalg.solve(frames[1],changes[1]@frames2[1])
    assert np.allclose(matrix2,left_u.T@matrix@right_u,rtol=4e-12,atol=2e-13)


def test_complex_holomorphic_slot_cannot_use_conjugate_transpose():
    tensor,grams,frames,n=example()
    raw_h=frames[2]@n
    raw_y=np.einsum("abc,c->ab",tensor,raw_h)
    correct=canonical_mass_matrix(tensor,*grams,n,vev=math.sqrt(2),kappa=1)
    wrong=frames[0].conj().T@raw_y@frames[1]
    difference=np.max(np.abs(mass_singular_values(correct)-mass_singular_values(wrong)))/mass_singular_values(correct)[0]
    assert difference > .09


def test_supplied_metrics_change_magnitude_and_ratios_without_changing_rank():
    from mtft.matter.collinear import yukawa_witness
    y=np.array(yukawa_witness().residue_matrix,dtype=complex)[:,:,None]
    mass=canonical_mass_matrix(y,np.diag([1,4,25]),np.diag([1,9,16]),[[1]],[1],vev=math.sqrt(2),kappa=2)
    assert np.allclose(mass_singular_values(mass),[2,1/4,2/15],atol=2e-14)
    assert np.linalg.matrix_rank(mass)==3
    changed_h=canonical_mass_matrix(y,np.diag([1,4,25]),np.diag([1,9,16]),[[4]],[1],vev=math.sqrt(2),kappa=2)
    assert np.allclose(changed_h,mass/2,atol=2e-14)


def test_zero_vev_and_zero_parent_coefficient_are_explicit_inputs():
    tensor,grams,_,n=example()
    assert np.array_equal(canonical_mass_matrix(tensor,*grams,n,vev=0,kappa=1),np.zeros((3,3)))
    assert np.array_equal(canonical_mass_matrix(tensor,*grams,n,vev=1,kappa=0),np.zeros((3,3)))
    with pytest.raises(TypeError):
        canonical_yukawa(tensor,*grams)
    with pytest.raises(TypeError):
        canonical_mass_matrix(tensor,*grams,n,kappa=1)
    with pytest.raises(TypeError):
        canonical_mass_matrix(tensor,*grams,n,vev=1)


@pytest.mark.parametrize("scale",[1e-300,1e300])
def test_finite_uniform_kinetic_scales_do_not_overflow_in_symmetrization(scale):
    gram=scale*np.eye(2)
    c=canonical_frame(gram)
    assert np.allclose(c.conj().T@gram@c,np.eye(2),atol=2e-14)


def test_rectangular_canonical_mass_matrix_has_declared_singular_values():
    matrix=np.array([[0,0,3j,0],[0,2,0,0]],dtype=complex)
    assert np.allclose(mass_singular_values(matrix),[3,2],atol=2e-14)
    assert np.array_equal(mass_singular_values(np.zeros((3,4))),np.zeros(3))


@pytest.mark.parametrize("gram",[
    [[1,2],[2,1]], [[1,0],[0,0]], [[1,1j],[0,1]],
    [[float("nan")]], [[complex(1,float("inf"))]], [], [1,2], [[1,2]],
    [[True]], [[1,False],[False,1]], [["1"]], [[0]],
])
def test_invalid_kinetic_forms_are_rejected(gram):
    with pytest.raises(ValueError): canonical_frame(gram)


@pytest.mark.parametrize("value",[-1,float("nan"),float("inf"),1j,True,"246"])
def test_invalid_vev_is_rejected(value):
    tensor,grams,_,n=example()
    with pytest.raises(ValueError): canonical_mass_matrix(tensor,*grams,n,vev=value,kappa=1)


@pytest.mark.parametrize("value",[float("nan"),float("inf"),complex(1,float("nan")),True,"1",[1]])
def test_invalid_parent_coefficient_is_rejected(value):
    tensor,grams,_,_=example()
    with pytest.raises(ValueError): canonical_yukawa(tensor,*grams,kappa=value)


@pytest.mark.parametrize("direction",[[0,0],[1,1],[1],[1,0,0],[[1,0]],[float("nan"),0],[True,False],["1","0"]])
def test_invalid_or_nonunit_higgs_direction_is_never_silently_normalized(direction):
    tensor,grams,_,_=example()
    with pytest.raises(ValueError): canonical_mass_matrix(tensor,*grams,direction,vev=1,kappa=1)


@pytest.mark.parametrize("tensor",[[],[[1]],np.zeros((0,3,2)),np.zeros((3,3)),np.full((3,3,2),np.nan)])
def test_invalid_tensor_shapes_or_values_are_rejected(tensor):
    _,grams,_,_=example()
    with pytest.raises(ValueError): canonical_yukawa(tensor,*grams,kappa=1)


def test_gram_sizes_must_match_tensor_slots():
    tensor,grams,_,_=example()
    for slot in range(3):
        invalid=list(grams); invalid[slot]=np.eye(len(grams[slot])+1)
        with pytest.raises(ValueError,match="dimensions"):
            canonical_yukawa(tensor,*invalid,kappa=1)


@pytest.mark.parametrize("matrix",[[],[1,2],[[float("inf")]],np.zeros((2,0)),[[True]],[["2"]]])
def test_invalid_mass_matrices_are_rejected(matrix):
    with pytest.raises(ValueError): mass_singular_values(matrix)


def test_finite_inputs_whose_products_overflow_are_rejected():
    with pytest.raises(ValueError,match="finite numerical domain"):
        canonical_yukawa([[[1e200]]],[[1]],[[1]],[[1]],kappa=1e200)
