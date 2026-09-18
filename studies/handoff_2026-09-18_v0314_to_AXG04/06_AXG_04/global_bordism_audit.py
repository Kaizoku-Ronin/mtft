#!/usr/bin/env python3
"""Exact low-degree spin-AHSS check for BSU3 x BSU2 x BU1 x BU1.

The mathematical input is the polynomial integral cohomology ring, the low
spin-bordism coefficients, and the standard spin-AHSS d2 (dual to Sq^2, with
mod-two reduction on the q=0 row). The finite-field computation is independent
of the anomaly-polynomial audit. It does not construct a differential GS/Wu
counterterm and does not compute a background category enlarged by B fields.
"""
from pathlib import Path
import itertools
import json

ROOT=Path(__file__).resolve().parent
generators=('h','u','C','L','D')
degrees=(2,2,4,4,6)
zero=(0,0,0,0,0)
checks=[]

def check(label,condition):
    assert condition,label
    checks.append(label)

def mono(h=0,u=0,C=0,L=0,D=0):
    return (h,u,C,L,D)

def degree(m):
    return sum(a*b for a,b in zip(m,degrees))

def name(m):
    out=[]
    for generator,power in zip(generators,m):
        if power: out.append(generator if power==1 else generator+'^'+str(power))
    return '*'.join(out) or '1'

def add_mod2(polynomial,monomial):
    if monomial in polynomial: polynomial.remove(monomial)
    else: polynomial.add(monomial)

def sq2(m):
    """Cartan rule; Sq1=0 since every generator has an integral lift."""
    answer=set()
    # Sq2 h=h^2, Sq2 u=u^2, Sq2 C=D; Sq2 L=Sq2 D=0.
    replacements={0:mono(h=2),1:mono(u=2),2:mono(D=1)}
    for i,replacement in replacements.items():
        if m[i]%2:
            target=list(m)
            target[i]-=1
            target=tuple(a+b for a,b in zip(target,replacement))
            add_mod2(answer,target)
    return answer

def sq2_poly(poly):
    out=set()
    for m in poly:
        for image in sq2(m): add_mod2(out,image)
    return out

def homogeneous_basis(d):
    bounds=[range(d//w+1) for w in degrees]
    return sorted([e for e in itertools.product(*bounds) if degree(e)==d],reverse=True)

B4=[mono(h=2),mono(h=1,u=1),mono(u=2),mono(C=1),mono(L=1)]
B6=[mono(h=3),mono(h=2,u=1),mono(h=1,u=2),mono(u=3),
    mono(h=1,C=1),mono(u=1,C=1),mono(h=1,L=1),mono(u=1,L=1),mono(D=1)]
B8=homogeneous_basis(8)
check('H4 basis is complete',set(B4)==set(homogeneous_basis(4)))
check('H6 basis is complete',set(B6)==set(homogeneous_basis(6)))
check('H4 dimension five',len(B4)==5)
check('H6 dimension nine',len(B6)==9)
check('H8 dimension sixteen',len(B8)==16)
check('no odd integral cohomology through degree seven',all(not homogeneous_basis(i) for i in [1,3,5,7]))
check('Sq2 h',sq2(mono(h=1))=={mono(h=2)})
check('Sq2 u',sq2(mono(u=1))=={mono(u=2)})
check('Sq2 color c2',sq2(mono(C=1))=={mono(D=1)})
check('Sq2 weak c2',not sq2(mono(L=1)))
check('Sq2 color c3',not sq2(mono(D=1)))
check('Cartan mixed Abelian image',sq2(mono(h=1,u=1))=={mono(h=2,u=1),mono(h=1,u=2)})
check('Sq2 squares vanish on relevant integral classes',all(not sq2_poly(sq2(m)) for m in B4+B6))

def operation_matrix(domain,codomain):
    images=[sq2(m) for m in domain]
    check('Sq2 images stay in expected homogeneous basis '+str(degree(domain[0])),
          all(image<=set(codomain) for image in images))
    return [[int(target in image) for image in images] for target in codomain]

def rref2(matrix):
    A=[row[:] for row in matrix]
    nrows=len(A); ncols=len(A[0]) if nrows else 0
    pivot_row=0; pivots=[]
    for col in range(ncols):
        found=next((row for row in range(pivot_row,nrows) if A[row][col]),None)
        if found is None: continue
        A[pivot_row],A[found]=A[found],A[pivot_row]
        for row in range(nrows):
            if row!=pivot_row and A[row][col]:
                A[row]=[a^b for a,b in zip(A[row],A[pivot_row])]
        pivots.append(col)
        pivot_row+=1
        if pivot_row==nrows: break
    return A,pivots

def apply(M,v):
    return tuple(sum(a*b for a,b in zip(row,v))%2 for row in M)

def transpose(M):
    return [list(row) for row in zip(*M)]

A=operation_matrix(B4,B6)
B=operation_matrix(B6,B8)
RA,pa=rref2(A)
RB,pb=rref2(B)
check('Sq2 H4 to H6 rank two',len(pa)==2)
check('Sq2 H6 to H8 rank seven',len(pb)==7)
check('cohomology square differential zero',all(not any(apply(B,col)) for col in transpose(A)))

kernel={v for v in itertools.product([0,1],repeat=9) if not any(apply(B,v))}
image={apply(A,v) for v in itertools.product([0,1],repeat=5)}
check('exhaustive nine-dimensional kernel has four vectors',len(kernel)==4)
check('exhaustive five-dimensional image has four vectors',len(image)==4)
check('exhaustive cohomology exactness kernel equals image',kernel==image)
v_mixed=tuple(int(m in {mono(h=2,u=1),mono(h=1,u=2)}) for m in B6)
v_D=tuple(int(m==mono(D=1)) for m in B6)
explicit_span={tuple(a*x^b*y for x,y in zip(v_mixed,v_D)) for a,b in itertools.product([0,1],repeat=2)}
check('explicit kernel basis is mixed cubic and color c3',kernel==explicit_span)

# E2_{6,1}: outgoing d2=A^T, incoming d2=B^T composed with integral mod2
# reduction. H8(BG;Z) is free with this polynomial basis, so reduction onto
# H8(BG;F2) is surjective.
outgoing=transpose(A)
incoming=transpose(B)
check('homology outgoing d2 rank two',len(rref2(outgoing)[1])==2)
check('homology incoming d2 rank seven',len(rref2(incoming)[1])==7)
check('homology differential square zero',all(not any(apply(outgoing,col)) for col in transpose(incoming)))
hom_kernel={v for v in itertools.product([0,1],repeat=9) if not any(apply(outgoing,v))}
# Independent exhaustive image uses the seven RREF-independent columns only.
incoming_pivots=rref2(incoming)[1]
incoming_basis=[transpose(incoming)[i] for i in incoming_pivots]
hom_image=set()
for coeffs in itertools.product([0,1],repeat=len(incoming_basis)):
    hom_image.add(tuple(sum(c*v[i] for c,v in zip(coeffs,incoming_basis))%2 for i in range(9)))
check('homological kernel has dimension seven',len(hom_kernel)==128)
check('homological image has dimension seven',len(hom_image)==128)
check('homological E3 term vanishes exactly',hom_kernel==hom_image)

# Spin bordism of a point in dimensions 0..7, and polynomial-ring homology,
# isolate the only possible total-degree-seven term.
spin={0:'Z',1:'Z2',2:'Z2',3:'0',4:'Z',5:'0',6:'0',7:'0'}
possible=[]
for qq in range(8):
    pp=7-qq
    if spin[qq]!='0' and homogeneous_basis(pp): possible.append([pp,qq])
check('unique E2 total-degree-seven term',[ [6,1] ]==possible)
check('spin-bordism coefficient at dimension seven is zero',spin[7]=='0')

result={
 'investigation':'AXG-04 Spin-BG degree-seven bordism check',
 'class':'Exact finite-field AHSS calculation conditional on stated standard topological inputs',
 'gauge_group':'SU(3) x SU(2) x U(1)_h x U(1)_X, direct product',
 'integral_cohomology':'Z[h,u,C,L,D], degrees (2,2,4,4,6), torsion-free; homology likewise free and even',
 'generator_meanings':{'h':'c1(U1_h)','u':'c1(U1_X)','C':'c2(SU3)','L':'c2(SU2)','D':'c3(SU3)'},
 'Sq1':'zero on all generators because they lift integrally',
 'Sq2_rules':{'h':'h^2','u':'u^2','C':'D','L':'0','D':'0'},
 'bases':{'H4':[name(v) for v in B4],'H6':[name(v) for v in B6],'H8':[name(v) for v in B8]},
 'Sq2_H4_to_H6':{'matrix':A,'rank':len(pa)},
 'Sq2_H6_to_H8':{'matrix':B,'rank':len(pb)},
 'kernel_H6_equals_image_H4':{'basis':['h^2*u+h*u^2','D'],'dimension':2,'exhaustive_vectors':len(kernel)},
 'spin_point_coefficients':spin,
 'AHSS':{'E2_definition':'E2_pq=Hp(BG;Omega_q^Spin(pt))',
  'only_total_degree_7_term':possible,'E2_6_1_dimension_F2':9,
  'd2_outgoing_rank':2,'d2_incoming_rank':7,'E3_6_1_dimension_F2':0,
  'reason_incoming_integral_reduction_surjective':'H8(BG;Z) is free in the displayed monomial basis',
  'higher_differentials':'No total-degree-seven associated-graded group survives E3; later pages cannot create one.',
  'extensions':'All degree-seven associated-graded pieces are zero; there is no extension ambiguity.'},
 'conclusion':'Omega_7^Spin(B[SU3 x SU2 x U1 x U1])=0',
 'physical_scope':[
  'For ordinary spin manifolds with bundles for exactly this direct-product group, no additional Spin-BG bordism torsion anomaly remains once cancellation of the full anomaly theory by globally defined GS terms is constructed.',
  'This calculation does not itself define the differential-cohomological Green-Schwarz/Wu counterterm or its coupling to dynamical tensors.',
  'If the background category includes a chosen B-field trivialization such as dH=X4, it may be a homotopy fiber rather than BG. Its bordism group must be calculated separately.',
  'The result does not apply automatically to quotient gauge groups, spin-charge identifications, defects, boundaries, or non-spin manifolds.',
  'This is not a UV-completion, positive-energy, scalar-potential, or vacuum theorem.'
 ],
 'references':[
  {'url':'https://arxiv.org/abs/2012.11622','sections':'Appendix B.1-B.3',
   'support':'Spin AHSS d2 dual to Sq2 with mod2 reduction, and SU2/SU3 cohomology and Sq2 c2=c3.'},
  {'url':'https://arxiv.org/abs/2012.11622','sections':'Section 3.3 and Appendix B.6',
   'support':'The B-field-trivialization background space can differ from BG and have nonzero degree-seven bordism.'}
 ],
 'checks':checks,'summary':{'passed':len(checks),'failed':0}}
(ROOT/'global_bordism_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result['summary']))
