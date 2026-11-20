# WL-01 — walls between twists: the wall distance, the systole band and one-way costs (2026-11-20)

Module `mtft.lines.walls`; tests `tests/test_wl01_walls.py`. Ported from the project records FF-14 (closed geodesics by class),
FF-17 (walls) and FF-20 (the wall metric), 18–20 November 2026, filed under `claude/experiments/` in the MTFT Python project.

A wall in 4D across which the twist changes is a vortex wrapping a curve on X₀(143) (FF-17). A curve of class Y changes ζ by
C⁻¹Y (`mtft.lines.twists`), so the least wall from ζ₁ to ζ₂ has length N(C(ζ₁ ⊕ ζ₂)). Here N(y) is the least total length of
closed geodesics whose classes add to y mod 2. d_W(ζ₁, ζ₂) = N(C(ζ₁ ⊕ ζ₂)) is a translation-invariant metric on the 2²⁶ twists.

What it holds:
- `level1_classes`, `lifts`, `primitive_classes`, `shortest_by_class`: every primitive oriented closed geodesic up to a length,
  with its class, and the shortest one in each class mod 2.
- `wall_ball`: every class within a radius, with N and one least wall. It is a Dijkstra search over the classes reached, and it
  is exact because a wall of length ≤ r uses geodesics ≤ r only. `three_sums` and `compare_routes` give the second route.
- `wall_distance`, `wall`: d_W and a least wall within a radius. None certifies that the distance exceeds the radius.
- `systoles`, `systole_band`, `check_systoles`: the four systoles (2 arccosh 2), and how many of them a twist makes antiperiodic.
- `thin_wall_exponent`, `bottleneck`: Coleman's B = 27π²σ⁴/(2ε³), and the least over routes of a route's hardest step,
  c(x, z) ≤ max(c(x, y), c(y, z)).

Status: EXACT within the radius. Gates and routes:
- The level-1 classes are found two ways up to trace 20: canonical block words, and all positive words by minimal rotation.
- FF-20's ball sizes on its grid up to 12: 1, 5, 7, 11, 33, 289, 1,020 and 5,559.
- The sum of at most three geodesics agrees on all 1,019 walls up to radius 10 and never beats the ball.
- FF-17's records: δ = 6533120 (φ*η to FF-14's weakest twist) has dual class 12664832 and a one-geodesic wall of 9.8974. The
  minimum's dual class is 40693698, and it has no wall within 10.
- Systoles: eight oriented geodesics in four classes, of rank 4 mod 2, with pairwise intersection 0. So the bands hold
  2²²·(1, 4, 6, 4, 1) twists, and a wall along a systole never changes a twist's band.
- Outside the suite: the package was compared class by class with FF-20's numba programme (h8 plus one long piece) on all 120,827
  classes within 16. The classes are the same and |ΔN| ≤ 4.8 × 10⁻⁷, which is the record's float32 precision. The radius-16 ball
  takes 63 s.

How many geodesics a least wall uses is not an invariant. Within 16, 64 classes have least walls of two sizes, equal in length
to 10⁻⁹: 24 classes have walls of one or two geodesics, and 40 have walls of two or three. Every piece of every tie is a
multiple of the systole length L = 2 arccosh 2: 5L = 2L + 3L (16 classes), 6L = 2L + 4L (8) and L + 5L = L + 2L + 3L (40).
Geodesics of length kL lie over level-1 classes of trace 2T_k(2) = 4, 14, 52, 194, 724, whose units are powers of 2 + √3. Of
the 368 geodesics of length 5L, 8 are fifth lifts of the systole's class and 360 lie over primitive classes of trace 724. By
number of geodesics, the package's walls are 60,218 / 47,619 / 12,288 / 691 / 10 and FF-20's are 60,202 / 47,595 / 12,328 /
691 / 10. The differences (+16, +24, −40) are what these ties allow. The distance is the invariant.

What this is not: the twist energy. The depth of each band therefore comes from the record (FF-20: each systole twisted adds about
0.23 of the minimum's depth). The wall between the two energy minima stays bracketed in [16, 19.12] (FF-20).

Correction (20 November 2026, before any push): the first version of this note said the tied geodesics were lifts of the
systole's level-1 class. Only 8 of them are; the rest lie over their own level-1 classes in orders of ℚ(√3).
