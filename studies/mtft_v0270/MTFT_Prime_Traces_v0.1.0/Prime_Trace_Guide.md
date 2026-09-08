# Prime traces: infinitude proved and infinitude open

MTFT companion experiment v0.1.0 · 6 September 2026

The trace separates the mathematical status of a family from membership of
an individual prime. Each recorded membership is decided by exact finite
arithmetic. Whether the family continues infinitely is a separate field in
`family_catalogue.json`.

## The two groups

| Group | Families used in this first trace | Mathematical basis |
|---|---|---|
| Infinitude proved | All primes; each class p ≡ r mod 30 for r = 1, 7, 11, 13, 17, 19, 23, 29 | Euclid's theorem and Dirichlet's theorem |
| Infinitude open | Twin-prime members; Sophie Germain primes; safe primes; Wieferich primes in base 2 | The respective infinitude problems remain open |

Dirichlet's theorem supplies infinitely many primes in every fixed arithmetic
progression whose residue and modulus are coprime. All eight residues above
satisfy that condition. Every prime greater than 5 lies on exactly one of
these rows. The primes 2, 3, and 5 are retained as explicit finite exceptions
to this row system. See [Sutherland's MIT notes, Theorem 17.1](https://math.mit.edu/classes/18.785/2015fa/LectureNotes17.pdf).

A twin-prime member is a prime p for which p−2 or p+2 is prime. We count each
prime once, even if it has two twin partners. Infinitely many bounded prime
gaps are known; that theorem does not establish infinitely many gaps of size
2. See [the Polymath retrospective](https://arxiv.org/abs/1409.8361).

A Sophie Germain prime p has 2p+1 prime. A safe prime s has (s−1)/2 prime.
Their infinitude questions are equivalent under s = 2p+1. Their counts at
the same upper bound differ: safe primes up to x correspond exactly to
Sophie Germain primes up to (x−1)/2. The underlying open problem is discussed
in [Wikstrom's thesis on Cunningham chains](https://www.csc.kth.se/utbildning/forskar/avhandlingar/doktor/2005_2006/WikstromDouglas.pdf).

For an odd prime p, being Wieferich in base 2 means
2^(p−1) ≡ 1 mod p². The only known examples are 1093 and 3511, and infinitude
is unresolved. The base is an essential part of the definition. Infinitude
of their non-Wieferich complement is also not known unconditionally, so we
do not automatically promote an unmarked complement to the proved group.
See [Graves and Weiss, Introduction](https://arxiv.org/html/2503.19144).

## Reading the web

`Prime_Trace_Web.png` shows two finite windows around the two known base-2
Wieferich primes. Each horizontal row represents a residue family whose
infinitude is proved. Every grey point is a prime; colored symbols add
memberships whose infinitude remains open. Vertical offsets keep overlapping
symbols visible and carry no numerical meaning. Alternating bands delimit
consecutive-square intervals. The two panels use their own horizontal scales.

The samples were chosen to show known Wieferich examples, so these windows
are not a random sample for estimating rarity. Full-range counts are in the
second figure and the data files.

At 1093, the row is 13 mod 30, and the labels are twin-prime member and
Wieferich in base 2. At 3511, the row is 1 mod 30 and the Wieferich label is
present, with no twin partner. This is how the open families appear within
the proved-infinite rows: a prime retains every applicable membership.

Infinitude of a containing row does not establish infinitude of a marked
subset or of an intersection. A missing symbol in an interval means that
the tested property has no member there; it is not evidence that the family
is finite. The labels are intentionally overlapping and must not be added
as though they partitioned all primes.

## Some of the pattern is forced by congruences

| Open family | Necessary residue classes modulo 30 | Exceptions handled separately |
|---|---|---|
| Twin-prime member | 1, 11, 13, 17, 19, 29 | p ≤ 7 |
| Sophie Germain | 11, 23, 29 | p ≤ 5 |
| Safe prime | 17, 23, 29 | p ≤ 11 |

These restrictions follow from divisibility by 2, 3, and 5. For example, a
Sophie Germain prime p > 5 must have both p and 2p+1 coprime to 30. Checking
the eight possible residues leaves 11, 23, and 29. Mapping those residues
through p ↦ 2p+1 gives the safe-prime restriction. For twin members p > 7,
residue 7 makes p−2 divisible by 5 and p+2 divisible by 3; residue 23 makes
p−2 divisible by 3 and p+2 divisible by 5.

These are necessary conditions, not sufficient ones. They explain visible
empty rows before any statistical interpretation. The finite data satisfy
all three restrictions exactly. No analogous extra modulo-30 restriction
is assumed for the Wieferich overlay.

## First finite census

The census contains all **78,650 primes below 1,002,001**, covering the first
**1,000 square intervals**, (n²,(n+1)²) for 1 ≤ n ≤ 1000. The endpoint is a
composite square. Partner primality is tested beyond each displayed interval
and beyond the census endpoint when the predicate requires it.

| Open-infinitude property | Members in this census |
|---|---:|
| Twin-prime member | 16,371 |
| Sophie Germain | 7,760 |
| Safe prime | 4,333 |
| Wieferich in base 2 | 2 |

The eight residue-row counts, in the declared order, are
9,827; 9,827; 9,828; 9,846; 9,831; 9,808; 9,856; 9,824. They sum to 78,647;
adding the three exceptional primes recovers the full count.

`Prime_Trace_Growth.png` plots finite cumulative counts and each open family's
fraction of all primes. A fraction can decline while the number of observed
members grows. Neither behavior by itself resolves infinitude. Zero counts
are omitted from the logarithmic panels; no curve is extrapolated beyond
the computed range.

Every square interval in this finite run is nonempty. This is a reproducible
small-range check, not a new record or a proof of Legendre's conjecture.
The conjecture requires a prime in every such interval. A proved-infinite
family can still leave individual square intervals empty. See
[Chamberland and Straub's 2026 survey](https://arxiv.org/html/2602.22502v1).

## Exact trace fields and verification

For each prime we store its integer value, square interval, exact fractional
position within that interval, residue modulo 30, four membership flags,
the previous prime, and the gap to the previous member of each marked
family. First occurrences use an explicit missing sentinel for that gap.

We also retain the Fermat quotient q₂(p) = (2^(p−1)−1)/p mod p, computed with
modular exponentiation. It is exactly zero for the Wieferich predicate and
is not applicable at p = 2. A theorem status of OPEN never means a finite
membership test was left unfinished.

The checks include independent trial division for all integers through
5,000; 209 independently checked prime/family rows; six known prime-count
anchors through one million; residue partition conservation at every square
endpoint; the Sophie-Germain/safe-prime cutoff correspondence at all 1,000
endpoints; and unreduced integer-power checks for both Wieferich positives.
Four CLI checks passed. Both scientific figures were rendered and visually
inspected.

Counts, membership, and rational positions are exact finite results. The
status of infinitude comes from the stated mathematics and literature, not
from these checks. At prime p the MTFT weight is log(p)/p; the additional
properties here provide the type diagnostics discussed in our experiment
design. This companion adds data and tracing tools without modifying MTFT
0.26.0.
