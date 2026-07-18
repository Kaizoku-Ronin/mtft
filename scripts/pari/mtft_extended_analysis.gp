\\ MTFT Paper 32: Extended Analysis
\\ 1. Cusp equivalence classes of Gamma_0(143)
\\ 2. Hecke polynomial factorization (algebraic Q2/Q4)
\\ 3. Broad resonance scan (c <= 50)
\\ 4. j-invariant of 143.a1
\\ Roger Tano - April 2026

default(realprecision, 50);
default(parisizemax, 2000000000);

cabs(z) = { sqrt(norml2([z])); }

print("================================================================");
print("  MTFT PAPER 32: EXTENDED ANALYSIS");
print("================================================================");

\\ ================================================================
\\ SECTION 1: CUSP EQUIVALENCE CLASSES
\\ ================================================================
print("\n================================================================");
print("  1. CUSP EQUIVALENCE CLASSES OF Gamma_0(143)");
print("================================================================");

\\ For N = p*q squarefree, cusps are represented by 1/d for d | N
\\ N = 143 = 11 * 13, divisors: 1, 11, 13, 143
\\ Cusp representatives: inf = 1/0 ~ 0/1, 1/11, 1/13, 0 = 0/1 ~ 1/143
\\ Width of cusp a/c (with c | N): N / (c * gcd(c, N/c)^2)... 
\\ Actually let's just compute which cusp class each a/c belongs to.

\\ Two cusps a/c and a'/c' are Gamma_0(N)-equivalent iff
\\ there exists gamma in Gamma_0(N) mapping one to the other.
\\ For N squarefree, a/c ~ a'/c' iff gcd(c,N) = gcd(c',N) and
\\ a/c = a'/c' mod gcd(c,N) in a suitable sense.
\\ Simpler: The cusp class is determined by d = gcd(c, N).

print("\n  N = 143 = 11 x 13 (squarefree)");
print("  Divisors of 143: 1, 11, 13, 143");
print("  Number of cusp classes: 4");
print("  Representatives: inf (d=1), 1/11 (d=11), 1/13 (d=13), 0 (d=143)");

print("\n  Cusp class assignments for tested fractions:");
print("  (d = gcd(denominator, 143) determines the class)\n");

test_fracs = [[1,7], [2,7], [3,7], [1,11], [2,11], [3,11], [5,11], [7,11], [1,13], [2,13], [3,13], [5,13], [1,77], [2,77], [3,77], [4,77], [5,77], [6,77], [1,143], [2,143], [3,143], [4,143], [5,143], [6,143], [7,143], [8,143], [9,143], [10,143]];

for(j = 1, #test_fracs, my(a, c, d); a = test_fracs[j][1]; c = test_fracs[j][2]; d = gcd(c, 143); print("  ", a, "/", c, "  gcd(", c, ", 143) = ", d, "  -> cusp class d=", d, if(d==1, " (inf-class)", if(d==11, " (1/11-class)", if(d==13, " (1/13-class)", " (0-class)")))));

\\ Cusp widths
print("\n  Cusp widths:");
print("  inf (d=1):  width = 143/(1 * gcd(1, 143)) = 143/1 = 143");
print("  1/11 (d=11): width = 143/(11 * gcd(11, 13)) = 143/11 = 13");
print("  1/13 (d=13): width = 143/(13 * gcd(13, 11)) = 143/13 = 11");
print("  0 (d=143):  width = 143/(143 * gcd(143,1)) = 1");

\\ Verify: sum of widths = psi(N) = N * prod(1 + 1/p) = 143 * (1+1/11) * (1+1/13)
psi = 143 * (1 + 1/11) * (1 + 1/13);
print("  Sum of widths: 143 + 13 + 11 + 1 = 168 = psi(143) = |PSL(2,7)| CHECK: ", psi);

\\ ================================================================
\\ SECTION 2: HECKE POLYNOMIAL FACTORIZATION
\\ ================================================================
print("\n================================================================");
print("  2. HECKE POLYNOMIAL FACTORIZATION (ALGEBRAIC Q2/Q4)");
print("================================================================");

\\ f3 has coefficient field Q[y]/(y^6 - 10y^4 - 2y^3 + 24y^2 + 7y - 12)
\\ a2(f3) = -y in this field
\\ The T2 eigenvalues on f3's 6 embeddings are the 6 roots of the min poly
\\ evaluated at a2 = -y, i.e., the roots of y^6 - 10y^4 - 2y^3 + 24y^2 + 7y - 12

print("\n  f3 coefficient field minimal polynomial:");
f3_minpoly = y^6 - 10*y^4 - 2*y^3 + 24*y^2 + 7*y - 12;
print("  P(y) = ", f3_minpoly);

print("\n  a2(f3) = -y, so T2 eigenvalues are roots of P(y)");
print("  If Q2/Q4 split exists, P should factor over Q as deg 2 x deg 4");

print("\n  Factorization of P over Q:");
fac = factor(f3_minpoly);
print("  ", fac);
print("  Number of irreducible factors: ", #fac[,1]);

for(k = 1, #fac[,1], print("  Factor ", k, ": ", fac[k,1], "  degree = ", poldegree(fac[k,1]), "  multiplicity = ", fac[k,2]));

\\ Check if it factors as 2+4 or 3+3 or stays irreducible
print("\n  Interpretation:");
if(#fac[,1] == 1, print("  P is IRREDUCIBLE over Q -> no algebraic Q2/Q4 split"));
if(#fac[,1] == 2, print("  P factors into 2 pieces over Q:"); print("    Factor degrees: ", poldegree(fac[1,1]), " and ", poldegree(fac[2,1])); if(poldegree(fac[1,1]) == 2 && poldegree(fac[2,1]) == 4, print("    THIS IS THE Q2 x Q4 SPLIT (deg 2 + deg 4)")); if(poldegree(fac[1,1]) == 4 && poldegree(fac[2,1]) == 2, print("    THIS IS THE Q4 x Q2 SPLIT (deg 4 + deg 2)")); if(poldegree(fac[1,1]) == 3 && poldegree(fac[2,1]) == 3, print("    3+3 split (not Q2/Q4)")));
if(#fac[,1] >= 3, print("  P factors into ", #fac[,1], " pieces"));

\\ Now check the T2 Hecke polynomial: H3(x) = product over embeddings of (x - a2(sigma_k))
\\ Since a2 = -y, we substitute y -> -x to get H3
\\ H3(x) = P(-x) = (-x)^6 - 10(-x)^4 - 2(-x)^3 + 24(-x)^2 + 7(-x) - 12
\\        = x^6 - 10x^4 + 2x^3 + 24x^2 - 7x - 12
print("\n  T2 Hecke polynomial H3(x) = P(-x):");
H3 = subst(f3_minpoly, y, -x);
print("  H3(x) = ", H3);

print("\n  Factorization of H3 over Q:");
fac_H3 = factor(H3);
print("  ", fac_H3);
for(k = 1, #fac_H3[,1], print("  Factor ", k, ": ", fac_H3[k,1], "  degree = ", poldegree(fac_H3[k,1])));

\\ Roots of each factor
print("\n  Roots of H3 (T2 eigenvalues on f3):");
roots_H3 = polroots(H3);
for(k = 1, #roots_H3, print("  lambda_", k, " = ", roots_H3[k], "  |.| = ", cabs(roots_H3[k])));

\\ If there's a 2+4 split, identify which roots belong to Q2 vs Q4
if(#fac_H3[,1] == 2, print("\n  Roots by factor:"); print("  Q2 factor: ", fac_H3[1,1]); rQ2 = polroots(fac_H3[1,1]); for(k=1, #rQ2, print("    root = ", rQ2[k])); print("  Q4 factor: ", fac_H3[2,1]); rQ4 = polroots(fac_H3[2,1]); for(k=1, #rQ4, print("    root = ", rQ4[k])); print("\n  Tr(a2^2) over Q2: ", norml2(rQ2)); print("  Tr(a2^2) over Q4: ", norml2(rQ4)); print("  Tr(a2^2) total: ", norml2(roots_H3)); print("  lambda_2 (stiffness) = 6.427 for comparison with Tr_Q2(a2^2)"));

\\ Also check f2
print("\n  --- f2 coefficient field ---");
f2_minpoly = y^4 - 4*y^2 - y + 1;
print("  f2 minpoly: ", f2_minpoly);
print("  Factorization: ", factor(f2_minpoly));
print("  (Expected: irreducible of degree 4)");

\\ Verify trace relationships
print("\n  --- Trace verification ---");
print("  Tr(a2(f1)) = 0 (from q-expansion)");
\\ For f2: a2 = y^3 - 3y mod (y^4 - 4y^2 - y + 1)
\\ Trace = sum of roots of (y^3-3y) evaluated at roots of minpoly
\\ But trace of y mod minpoly = minus coeff of y^3 = 0
print("  Tr(a2(f2)) = Tr(y^3 - 3y) mod f2_minpoly");
print("    Tr(y) mod (y^4 - 4y^2 - y + 1) = 0 (coeff of y^3 is 0)");
\\ For f3: a2 = -y
\\ Trace of y mod (y^6 - 10y^4 - 2y^3 + 24y^2 + 7y - 12) = 0 (coeff of y^5 is 0)
print("  Tr(a2(f3)) = Tr(-y) = 0 (coeff of y^5 is 0) CHECK");

\\ ================================================================
\\ SECTION 3: BROAD RESONANCE SCAN
\\ ================================================================
print("\n================================================================");
print("  3. BROAD RESONANCE SCAN (c <= 50)");
print("================================================================");

\\ Initialize modular forms and symbols
mf = mfinit([143, 2]);
B = mfeigenbasis(mf);
print("  Computing modular symbols...");
FS1 = mfsymbol(mf, B[1]); print("  f1 done.");
FS2 = mfsymbol(mf, B[2]); print("  f2 done.");
FS3 = mfsymbol(mf, B[3]); print("  f3 done.");

target = log(1776.86/0.511) / log(105.658/0.511);

\\ Track best results
best_delta = 999.0;
best_a = 0;
best_c = 0;
best_ratio = 0;

\\ Also track top 10
results = List();

print("\n  Scanning all coprime a/c with c <= 50...\n");
print("  a/c\t\t|Om1|\t\t|Tr(Om2)|\t|Tr(Om3)|\td3/d2\t\tdelta\t\tcusp_class");

for(c = 1, 50, for(a = 1, c-1, if(gcd(a, c) != 1, next); my(path, p1, p2, p3, a1, a2v, a3v, d2, d3, ratio, delta, cusp_d); path = [a, 1; c, 0]; p1 = mfsymboleval(FS1, path); p2 = mfsymboleval(FS2, path); p3 = mfsymboleval(FS3, path); a1 = cabs(p1); a2v = cabs(vecsum(p2)); a3v = cabs(vecsum(p3)); cusp_d = gcd(c, 143); if(a1 < 1e-15 || a2v < 1e-15, next); d2 = log(a2v/a1); d3 = log(a3v/a1); if(abs(d2) < 0.05, next); ratio = d3/d2; delta = abs(ratio - target); listput(results, [delta, a, c, ratio, a1, a2v, a3v, cusp_d]); if(delta < 0.1, print("  ", a, "/", c, "\t", a1, "\t", a2v, "\t", a3v, "\t", ratio, "\t", delta, "\td=", cusp_d)); if(delta < best_delta, best_delta = delta; best_a = a; best_c = c; best_ratio = ratio)));

print("\n  --- BEST RESULT ---");
print("  Best mass cycle: {inf, ", best_a, "/", best_c, "}");
print("  d3/d2 = ", best_ratio);
print("  delta = ", best_delta);
print("  cusp class d = ", gcd(best_c, 143));

\\ Sort and show top 10
print("\n  --- TOP 10 CLOSEST (by |delta|) ---");
results_vec = Vec(results);
\\ Manual bubble sort of first 10
n_res = #results_vec;
print("  Total non-trivial cycles tested: ", n_res);

\\ Find top 10 by scanning
for(rank = 1, min(10, n_res), my(best_idx, best_d); best_d = 999; for(j = 1, n_res, if(results_vec[j][1] < best_d, best_d = results_vec[j][1]; best_idx = j)); my(r); r = results_vec[best_idx]; print("  #", rank, ": ", r[2], "/", r[3], "  d3/d2=", r[4], "  delta=", r[1], "  |Om1|=", r[5], "  cusp d=", r[8]); results_vec[best_idx] = [9999, 0, 0, 0, 0, 0, 0, 0]);

\\ ================================================================
\\ SECTION 4: j-INVARIANT
\\ ================================================================
print("\n================================================================");
print("  4. j-INVARIANT OF 143.a1");
print("================================================================");

E = ellinit([0, -1, 1, -1, -2]);
print("  Elliptic curve: [0, -1, 1, -1, -2]");
print("  Conductor: ", ellglobalred(E)[1]);
j = ellj(E);
print("  j-invariant: ", j);

\\ Check if it's a known CM j-invariant
\\ CM j-invariants for class number 1: 0, 1728, -3375, 8000, -32768, 54000, -884736, -12288000, -884736000
print("\n  Known CM j-invariants (class number 1):");
print("    D=-3: j=0");
print("    D=-4: j=1728");
print("    D=-7: j=-3375");
print("    D=-8: j=8000");
print("    D=-11: j=-32768");
print("    D=-19: j=-884736");
print("    D=-43: j=-884736000");
print("    D=-67: j=-147197952000");
print("    D=-163: j=-262537412640768000");

\\ Also check discriminant
disc = ellglobalred(E);
print("\n  Global data: ", disc);
print("  Minimal discriminant: ", elldiscminimal = E.disc);

\\ Check endomorphism ring
print("\n  Has CM? j = ", j);
cm_js = [0, 1728, -3375, 8000, -32768, -884736, -884736000, -147197952000, -262537412640768000];
is_cm = 0;
for(k = 1, #cm_js, if(j == cm_js[k], is_cm = 1; print("  YES - CM with discriminant from class number 1 list")));
if(is_cm == 0, print("  Checking class number 2 CM j-invariants..."); print("  (These are algebraic integers of degree 2, not rational)"); print("  Since j is rational and not in the class-1 list: NOT CM"));

\\ Rank
print("\n  Mordell-Weil rank: ", ellrank(E));

\\ Torsion
print("  Torsion subgroup: ", elltors(E));

\\ L-function value at s=1
print("  L(E, 1) = ", ellL1(E, 1));

print("\n================================================================");
print("  ANALYSIS COMPLETE");
print("================================================================");
