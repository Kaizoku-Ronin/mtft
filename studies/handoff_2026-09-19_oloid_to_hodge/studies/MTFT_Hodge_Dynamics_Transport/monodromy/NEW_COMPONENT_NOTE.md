# A loop detected only by the new Hodge component

The genus-three curve \(C_{3,c}:y^2=f_c^{\circ3}(x)\) splits, at the
rational-cohomology/isogeny level, into its inherited elliptic component
\(C_{2,c}\) and the new genus-two component

\[
D_{2,c}:v^2=(u-c)P_2(u,c),\qquad P_2(u,c)=(u^2+c)^2+c.
\]

This experiment compares their period-transport matrices on the same
parameter loop, using the exact rational connections in the bundle.

Let \(a=-1.75487766624669\ldots\) be the real root of
\(H(c)=c^3+2c^2+c+1\). It is a primitive period-three critical parameter:

\[
p_3(c)=P_2(c,c)=cH(c),\quad p_3(a)=0,\quad p_1(a)p_2(a)\ne0.
\]

The based loop starts at \(c=-2\), travels along the real axis to
\(a-0.05\), circles \(a\) counterclockwise with radius 0.05, and retraces
the stem. It encloses none of the other singular parameters.

**The inherited component has exactly trivial local monodromy.** Its
connection has poles only at \(0,-1\), so it extends holomorphically
across this loop's interior. Flat transport on that simply connected
region is path independent.

**The new component has nontrivial rank-one unipotent monodromy.** At
\(c=a\), precisely two branch points collide: the branch \(u=c\) and
one simple root of \(P_2\). In local coordinates \(s=u-c\),

\[
v^2=s\bigl[p_3(c)+B(c)s+O(s^2)\bigr],\qquad
B(a)=\partial_uP_2(a,a)=4a^2(a+1)\ne0.
\]

Since \(p_3'(a)\ne0\), completing the square gives a nodal smoothing
parameter proportional to \((c-a)^2\). Its Picard–Lefschetz monodromy
is therefore a nonzero unipotent transvection power: \(M_D=I+N_D\),
\(N_D^2=0\), rank \(N_D=1\). The square in the smoothing parameter is
also visible in the \(H(c)^2\) discriminant factor.

The original \(C_3\) first acquires a node with a simple smoothing
parameter at this same critical return. Its relationship to \(D_2\)
uses a degree-two quotient/isogeny. One must not equate integral
vanishing-cycle normalizations across that isogeny from these
complex-frame matrices alone.

## Numerical result

SciPy DOP853 integrations were repeated at relative tolerances
\(10^{-8},10^{-10},10^{-13}\). At the finest tolerance:

| Diagnostic | Value |
|---|---:|
| Inherited \(\|M_E-I\|\) | \(7.88\times10^{-16}\) |
| New \(\|M_D-I\|\) | 19.1323750280 |
| Singular values of \(M_D-I\) | \(19.1324,\;4.29\times10^{-15},\;2.99\times10^{-15},\;3.76\times10^{-16}\) |
| \(\|(M_D-I)^2\|\) | \(2.46\times10^{-13}\) |
| Relative nilpotence defect \(\|N_D^2\|/\|N_D\|^2\) | \(6.72\times10^{-16}\) |
| \(|\operatorname{tr}M_D-4|\) | \(1.25\times10^{-14}\) |
| \(|\det M_D-1|\) | \(3.04\times10^{-14}\) |
| Matrix difference from the \(10^{-10}\) run | \(6.94\times10^{-11}\) |

The singular values provide a numerical rank diagnostic, not an exact
rank proof. The local geometric argument gives the exact classification.
No eigenvalue splitting near one is used as a chaos diagnostic.

This is an explicit sense in which an added Hodge component can detect
a dynamical critical-return parameter that the inherited component does
not detect. The loop and its timing remain externally chosen; this
does not determine a natural parameter evolution or a Feigenbaum
renormalization operator.

## Reproduction

```sh
python3 new_component_monodromy.py
```

This reads `../geometry/transport_results.json` and writes
`new_component_monodromy_results.json`. It requires NumPy and SciPy.
All numerical values are complex128 computations, without interval
certification.
