# Returning to the oloid

This note records the final reflection after the three completed studies. It adds an explicit comparison of rational involutions and preserves the limits of that comparison. It is not a new physical model or a novelty claim.

## The contact relation

An oloid is the convex hull of two equal circles in perpendicular planes, each passing through the center of the other. For unit circles, let t and u be the angular parameters of the endpoints of the same straight generator/contact segment on its boundary.

Dirnböck and Stachel, *The Development of the Oloid*, Journal for Geometry and Graphics 1 (1997), 105–118, equation (6), establish

    cos(u) = -cos(t)/(1+cos(t)).

Put q=cos(t) and q'=cos(u). The relation is

    q' = w(q) = -q/(1+q).

Applying it twice gives q. In the real contact geometry the allowed q values lie between −1/2 and 1, with endpoint limits needed where the chosen smooth chart degenerates. The map exchanges the two ends of this interval and fixes q=0. This relation describes two contact points of one oloid; it does not arise from intersecting two independently oriented oloids.

Primary source: https://www.heldermann-verlag.de/jgg/jgg01_05/jgg0113.pdf

## Exact conjugacy to the triangle reflection

Define a new real coordinate

    n(q) = -1/(q+2).

Then

    w(q)+2 = (q+2)/(q+1),
    n(w(q)) = -(q+1)/(q+2)
             = -1 + 1/(q+2)
             = -1 - n(q).

Thus w becomes the triangle reflection sigma(n)=−1−n. Consequently T(n)=n(n+1)/2 is invariant under this transformed operation. This is an exact identity of rational maps on their domains, extending projectively if desired.

The physical interval q in [−1/2,1] maps to n in [−2/3,−1/3]. These are real coordinates, not nonnegative integer labels of triangular numbers. No integer-lattice or three-square representation theorem has been realized by the physical oloid through this coordinate change.

## Comparison with the elliptic family

The triangle/Hodge study uses

    C_(2,c): y²=x⁴+2cx²+c²+c,   c != 0,-1,
    j(c)=64(4c+3)³/(c+1).

Its degree-two isogenous partner has

    j_partner(c)=64(3-c)³/(c+1)².

The exact substitution

    c' = -c/(1+c)

exchanges these two j-functions. It is the Fricke parameter involution of the X0(2) parametrization described in the completed studies. It has the same rational formula as the oloid contact exchange.

This comparison identifies a common order-two algebraic structure. Such fractional involutions commonly become reflections after a coordinate change, so the formula by itself is not a unique fingerprint of a unified geometry. In the oloid it acts on contact parameters; in the elliptic construction it acts on a parameter labelling curves and isogenies.

There is also a concrete domain issue: q=0 is an ordinary oloid contact value, whereas setting c=0 in the displayed quartic gives a singular model. Therefore simply relabelling q as c does not produce a globally smooth correspondence of the two systems. The explicit isogeny lift in the transport report additionally requires a square-root choice. Neither a lift from oloid configurations nor compatibility with the rolling equations has been constructed.

## Holonomic rolling and the original sphere idea

Kuleshov, Hubbard, Peterson, and Gede, *Motion of the Oloid-toy*, ENOC 2011, analyze ideal rolling without slip. In their model the no-slip constraints are integrable, so describing that system as essentially holonomic is supported. This statement concerns mechanical constraints. Holonomy of a connection concerns transport around loops; the two terms must not be identified merely because they sound similar.

Primary author-uploaded paper: https://www.researchgate.net/publication/266032307_Motion_of_the_Oloid-toy

The original proposal to intersect two oppositely oriented solid oloids and obtain S4 cannot hold literally in R3. S4 is a four-dimensional manifold. The intersection of two convex bodies remains convex; when it has nonempty interior it is homeomorphic to a three-ball, with boundary homeomorphic to S2. Lower-dimensional or empty intersections do not resolve the dimensional obstruction. A higher-dimensional configuration space of arrangements would be a different object and must be defined separately.

The physical boundary has trivial first cohomology, H1(S2)=0. The nontrivial elliptic and higher-genus Hodge structures in the studies therefore do not live directly on that boundary. Associated algebraic curves or parameter/configuration spaces could carry them, but the required construction remains open.

## What this adds to the handoff

The oloid provides a tangible occurrence of the same rational involution already present in our modular calculation, and an exact coordinate change relates it to triangle reflection. It does not supply Feigenbaum renormalization: repeating an involution only produces fixed points or two-cycles. A future dynamical claim would require a distinct evolution law, the relevant configuration space, and explicit compatible maps.
