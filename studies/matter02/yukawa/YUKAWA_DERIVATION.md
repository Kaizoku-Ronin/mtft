# Matter-02: a bulk lepton–Higgs Yukawa survives the blow-up

For the specified holomorphic background, a bulk lepton–Higgs tensor is **nonzero**, and one explicitly selected holomorphic Higgs direction has a lepton coupling matrix of **rank at least three**. This is an exact holomorphic statement. It does not identify the three observed light families, compute a physical electron Yukawa, select a Higgs vacuum, or supply Hermitian normalization.

The existence proof is a finite minor of the full tensor rather than a calculation of the full tensor. Its three rows lie in an 18-dimensional space and its three columns in a 44-dimensional space. The dimensions are established by the separate Matter-02 cohomology certificate. The six local classes used below are not assumed to be the physical light families.

## 1. The actual holomorphic background

Let

\[
C=X_0(143),\qquad E_0=143a1,\qquad
\pi:\widetilde S=\operatorname{Bl}_{(P,Q)}(C\times E_0)\longrightarrow C\times E_0,
\]

where \(P=\infty\) is the rational cusp and \(Q\) is the elliptic origin. Write \(E_{\rm exc}\) for the exceptional curve and choose

\[
M=\mathcal O_C(P),\qquad
\mathcal L=\pi^*\operatorname{pr}_C^*M\otimes\mathcal O_{\widetilde S}(-3E_{\rm exc}).
\]

There is no degree-zero torus twist. Stack powers are

| Stack | \(c\) | \(L\) | \(A_1\) | \(A_2\) | \(B\) |
|---|---:|---:|---:|---:|---:|
| Exponent of \(\mathcal L\) | 0 | 1 | 2 | −3 | −1 |
| Rank | 3 | 2 | 1 | 1 | 1 |

The gauge background is split and its adjoint scalar \(\Phi\) is zero. All appearances of \(\Phi\) in this certificate refer to fluctuations or candidate Higgs directions. This study uses the partially twisted 8D supersymmetric Yang–Mills action adopted conditionally by the existing B3 construction.

## 2. The twist changes the field dictionary

For an oriented Hom block of exponent \(q\), the relevant four-dimensional chiral fields arise from

\[
a\in H^1(\widetilde S,\mathcal L^q),\qquad
\varphi\in H^0(\widetilde S,K_{\widetilde S}\otimes\mathcal L^q).
\]

They are respectively a \((0,1)\)-form and a \((2,0)\)-form internally. The holomorphic action has

\[
W\propto\int_{\widetilde S}\operatorname{Tr}
\bigl(F^{0,2}\wedge\Phi^{2,0}\bigr),\qquad
W_3\propto\int_{\widetilde S}\operatorname{Tr}(a\wedge a\wedge\varphi).
\]

Thus the Yukawa is a cup product followed by Serre duality:

\[
H^1(\mathcal L^a)\otimes H^1(\mathcal L^b)
\longrightarrow H^2(\mathcal L^{a+b}),
\quad
H^2(\mathcal L^{a+b})\otimes
H^0(K_{\widetilde S}\mathcal L^{-a-b})\longrightarrow\mathbb C.
\]

Three \(a\) fields alone have internal antiholomorphic degree three and do not give this cubic on a complex surface. Two \(a\) fields and one \(\varphi\) have degree \((2,2)\), which does. The older untwisted product-surface Künneth-parity tests cannot be substituted for this twisted blow-up calculation.

In the holomorphic chirality convention of this dictionary, a charged block has \(n_R=h^1(\mathcal L^q)+h^0(K_{\widetilde S}\mathcal L^q)\), and \(n_R-n_{R^*}=\chi(\mathcal L^{-q})-\chi(\mathcal L^q)=3q\). The existing B3 ledger calls the opposite sign \(e-24b=-3q\) its net index. Comparing them requires an explicit common chirality/orientation convention. Neither this choice of sign nor the finite minor proves that the observed three left-handed families have been selected.

Charged \(\varphi\) modes are present. With \(K_{\widetilde S}=\pi^*K_{C\times E_0}\otimes\mathcal O(E_{\rm exc})\), pushforward gives

\[
H^0(K_{\widetilde S}\mathcal L^q)=
\begin{cases}
H^0(C,K_C(-(2q-1)P)),&q>0,\\
H^0(C,K_C(-|q|P)),&q<0,\\
H^0(C,K_C)\otimes H^0(E_0,K_{E_0}),&q=0.
\end{cases}
\]

For the occurring charged exponents \(|q|\le5\), these spaces are nonzero. Although \(H^0(\mathcal L^q)=0\) for every \(q\ne0\), that vanishing does not remove \(\varphi\) modes because their bundles include \(K_{\widetilde S}\). The del Pezzo/Hirzebruch bulk-Yukawa vanishing theorem does not apply: this surface has \(h^{2,0}=13\), and its canonical geometry has a different sign.

## 3. An electron-type gauge triangle

The \(A_1\) triangle is

\[
\ell:A_1\overline L\ (q=1,Y=-\tfrac12),\quad
e^c:B\overline{A_1}\ (q=-3,Y=1),\quad
H_d:L\overline B\ (q=2,Y=-\tfrac12).
\]

Both exponent and hypercharge close. Its gauge Lie bracket also survives:

\[
\operatorname{Tr}\left([E_{A_1,L},E_{B,A_1}]E_{L,B}\right)=-1.
\]

The diagonal background commutes, but these off-diagonal fluctuations do not. Background commutativity therefore does not kill their cubic. We place the two lepton fluctuations in \(H^1(\mathcal L)\) and \(H^1(\mathcal L^{-3})\), and the Higgs in \(H^0(K_{\widetilde S}\mathcal L^2)\).

## 4. The positive classes are genuine cohomology classes

Choose local coordinates \((z,w)\) at \((P,Q)\), with \(z=q\) the cusp coordinate and a nonzero elliptic differential represented by \(dw\). On the product,

\[
\pi_*\mathcal L=M\otimes\mathcal I_p^3,\qquad
R^1\pi_*\mathcal L=0.
\]

The exact sequence

\[
0\longrightarrow M\mathcal I_p^3\longrightarrow M
\longrightarrow M\otimes\mathcal O_{Z_3}\longrightarrow0
\]

supplies a connecting map from fat-point jets to \(H^1(\mathcal L)\). In the local frame of \(M=\mathcal O(P)\), its one global section is represented by \(z\), and all global sections are independent of \(w\). Consequently the three jets

\[
p_1=w,\qquad p_2=w^2,\qquad p_3=zw
\]

have total degree below three, survive the global-section image, and yield three linearly independent cohomology classes. They belong to the fat-point contribution created by the blow-up; the product-only argument misses them.

## 5. The negative principal parts lift globally

For \(\mathcal L^{-3}=\pi^*M^{-3}\otimes\mathcal O(9E_{\rm exc})\), the local \(R^1\pi_*\mathcal O(9E_{\rm exc})\) is represented by

\[
z^{-i}w^{-j},\qquad i,j\ge1,\quad i+j\le9.
\]

This basis follows from the two blow-up charts \((z,w)=(u,uv)\) and \((z,w)=(st,t)\). Relative to the \(u^{-9}\) frame, overlap monomials \(u^av^b\) survive the two chart subspaces precisely when \(a\ge0\) and \(a-9<b<0\). Setting \(j=-b\), \(i=9-a-j\) gives the displayed principal parts.

The Leray obstruction sends these local classes to \(H^2(C\times E_0,M^{-3})\). By duality, it is evaluated by residues against sections of \(K_C\otimes M^3\), pulled from \(C\), times \(dw\). Those sections have no \(w\) dependence. Every principal part with \(j\ge2\) therefore has zero obstruction and lifts to global \(H^1(\mathcal L^{-3})\).

Choose

\[
n_1=z^{-6}w^{-2},\qquad
n_2=z^{-7}w^{-2},\qquad
n_3=z^{-6}w^{-3}.
\]

Each has \(i+j\le9\) and \(j\ge2\). Nonuniqueness of the global lift adds product cohomology; it does not alter this pairing with the fat-point connecting classes. Equivalently, the principal parts define residue functionals on the target jet quotient and annihilate its entire global-section image.

## 6. A certified holomorphic Higgs direction

The archived full canonical bases identify a differential

\[
\omega=\left(z^3+0\,z^4+O(z^5)\right)dz
\in H^0(C,K_C(-3P)).
\]

For \(\omega=f(z)dz/z\), this is the exact linear system \(a_1=a_2=a_3=a_5=0\), \(a_4=1\) on the coefficients of \(f(z)\). Exact row ranks in both the original and Atkin–Lehner-adapted archived bases are \(1,2,\ldots,9\) through the first nine cusp jets. The certificate saves an explicit rational combination of the original 13 basis forms.

Relative to the \(M^2\) frame, the coefficient of the resulting \(\varphi_2\) is

\[
c(z)=z^2\,\omega/dz=z^5+0\,z^6+O(z^7).
\]

The coefficient normalization is holomorphic, not Hermitian. The overall choice of differential and residue convention cannot turn a nonzero minor into zero.

## 7. The rank-three minor

The natural residue pairing extracts the coefficient of \(z^{-1}w^{-1}\) in \(p_i n_j c(z)\,dz\wedge dw\). For a general allowed Higgs coefficient \(c(z)=v_5z^5+v_6z^6+\cdots\), it yields

\[
Y_{\rm minor}=\begin{pmatrix}
v_5&v_6&0\\
0&0&v_5\\
0&v_5&0
\end{pmatrix},\qquad
\det Y_{\rm minor}=-v_5^3.
\]

Our explicit direction has \(v_5=1,v_6=0\), so the determinant is \(-1\). Higher Taylor coefficients do not contribute to these nine residues. This proves a rank lower bound of three for the full holomorphic matrix for this Higgs direction. The nonzero group bracket supplies the remaining gauge factor.

A second way to verify nonvanishing uses only exact sequences. Multiplication by \(\varphi_2\) maps

\[
H^1(M\mathcal I_p^3)\longrightarrow
H^1(K_C M^3\mathcal I_p^8).
\]

It sends \(w^2\) to \(z^5w^2\) modulo \(\mathcal I_p^8\). This jet has total degree seven and a nonzero \(w\) power, so it survives every global section of \(K_CM^3\). Its connecting class is nonzero, and perfect Serre duality supplies a partner in \(H^1(\mathcal L^{-3})\). This independent formulation also verifies that the local residue is a global cup-product witness.

## 8. Neutral scalar modes and physical limitations

For \(q=0\), there are 13 holomorphic \(\varphi\) modes per available neutral adjoint component. A scalar proportional to the full identity matrix has zero Lie bracket and decouples from this cubic. A relative stack-diagonal scalar can couple to an opposite-orientation pair because \(\operatorname{Tr}([E_{ij},E_{ji}]\Phi)=\Phi_i-\Phi_j\). Neutrality under the background flux alone does not make every such mode central.

Invertible canonical normalization preserves matrix rank, but it changes physical coupling magnitudes and singular values. Those kinetic Gram matrices have not been computed here. A charged Higgs VEV must also satisfy the coupled F- and D-term equations; a nonzero holomorphic overlap does not select that vacuum or prove that all unwanted partners become heavy. A projection onto a proposed three-dimensional light-family sector must be separately justified.

Open gates are the full tensor, selected Higgs direction and VEV, a controlled light-family sector, Hermitian normalization, exotic lifting, the global gauge/axion action, and the physical massless-photon kernel. The QED engine from Matter-01 still uses supplied electron mass and fine-structure constant.

## Reproduce

```bash
python yukawa_certificate.py
```

The script uses exact rational arithmetic from the Python standard library. It writes `yukawa_exact.json` with basis hashes, an explicit differential, the residue matrix and **34 passing finite checks**. The mathematical global-lift and field-dictionary arguments are above; these checks verify their finite inputs and algebra and are not a computer formalization of sheaf cohomology.

## Sources and scope

- [Beasley, Heckman and Vafa, *GUTs and Exceptional Branes in F-theory—I*, arXiv:0802.3391](https://arxiv.org/abs/0802.3391), sections 3.3.2–3.3.3, particularly equations 3.44 and 3.54–3.55, supply the standard twisted field dictionary and bulk trilinear.
- The two canonical q-expansion tables are copied from the preserved MTFT v0.34.0 archive. Their SHA-256 hashes appear in the result.
- The blow-up, jet and residue calculations in this note are conditional derivations for the specified B3 background. They do not assert an F-theory compactification or establish a complete quantum-gravity model.
