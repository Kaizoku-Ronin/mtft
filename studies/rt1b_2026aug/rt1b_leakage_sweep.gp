\\ RT-1B — the cusp boundary layer of X0(143), corrected and extended.
\\
\\ SUPERSEDES the version staged 2026-08-18, which applied the T_p formula
\\ at p = 11, 13.  Those divide the level, where the operator is U_p with no
\\ p^(k-1) term.  The wrong formula gave ungraded rank 41 > dim S_4 = 40,
\\ which is impossible; that is how it was caught.  This script selects the
\\ operator by divisibility and reproduces the corrected numbers.
\\
\\   T_p(f)_n = a_{pn} + p^(k-1) a_{n/p}    for p not dividing 143
\\   U_p(f)_n = a_{pn}                      for p dividing 143
\\
\\ Weight k = 4, Sturm bound 56.  Needs q-expansions to q^(56 p), so the
\\ depth below caps the prime range; v0.17.0 ships q^140 and reaches p = 2
\\ only.  Everything here is exact rational linear algebra.

default(parisize, 2^32);

{
readB(path) =
  my(L = readstr(path), R = List(), p, s);
  for(i = 1, #L,
    s = L[i];
    if(#s == 0, next);
    if(Vec(s)[1] == "#", next);
    p = strsplit(s, ",");
    if(p[1] == "i", next);
    listput(R, vector(#p-1, k, eval(p[k+1]))));
  matrix(#R, #R[1], r, c, R[r][c]);
}

DEPTH = 2700;          \\ 56 * 47 = 2632, so good primes up to 47
STB   = 56;            \\ weight-4 Sturm bound for level 143

mf = mfinit([143,2],1);
A  = mfcoefs(mf, DEPTH);
\\ Read the adapted basis from the RELEASED package data, not a scratch file.
\\ (A previous revision read studies/Bmat.txt, which is not in the v0.17.0
\\ sdist -- caught in audit. This path IS in the release.)
DATA = "src/mtft/canonical/_data/X0_143_AL_adapted_basis.txt";
Bm = readB(DATA);
print("adapted basis read from ", DATA, "; det = ", matdet(Bm), " [-1078272]");
E  = A * Bm;
print("adapted q-expansions to q^", DEPTH, "; integral: ", E == liftall(E));

ser(k) = sum(n=0, DEPTH, E[n+1,k]*z^n) + O(z^(DEPTH+1));
S = vector(13, k, ser(k));
pairs = List(); for(i=1,13, for(j=i,13, listput(pairs,[i,j]))); pairs = Vec(pairs);
P = vector(91, c, S[pairs[c][1]] * S[pairs[c][2]]);
vecof(f) = vectorv(STB+1, n, polcoef(f, n-1));

\\ operator selected by divisibility -- the whole point of this revision
op(f, p) = if(143 % p == 0, \
  sum(n=0, STB, polcoef(f, n*p)*z^n) + O(z^(STB+1)), \
  sum(n=0, STB, (polcoef(f, n*p) + p^3*if(n%p==0, polcoef(f, n/p), 0))*z^n) + O(z^(STB+1)));

sec  = [1,2,2,2,2,2,2,3,3,3,3,3,4];
mult = [1,2,3,4; 2,1,4,3; 3,4,1,2; 4,3,2,1];
names = ["(+,+)","(+,-)","(-,+)","(-,-)"];
idxof(s) = select(z -> mult[sec[pairs[z][1]], sec[pairs[z][2]]] == s, vector(91,z,z));
Vof(idx) = matconcat(vector(#idx, k, vecof(P[idx[k]])));
Tof(idx, p) = matconcat(vector(#idx, k, vecof(op(P[idx[k]], p))));
leak(idx, p) = matrank(matconcat([Vof(idx), Tof(idx,p)])) - matrank(Vof(idx));

Vall = Vof(vector(91,z,z));
print("dim H^0(2K) = ", matrank(Vall), "   [36]   dim S_4 = 40, boundary B = 4");
print("");
print("GOOD PRIMES  (T_p commutes with every W_Q, so the AL grading is defined)");
print("   p    (+,+) (+,-) (-,+) (-,-)   total");
forprime(p = 2, 47, if(143 % p != 0, \
  r = vector(4, s, leak(idxof(s), p)); \
  print("  ", if(p<10," ",""), p, "      ", r[1], "     ", r[2], "     ", r[3], "     ", r[4], "       ", vecsum(r))));

print("");
print("BAD PRIMES  p | 143  (U_p; per-sector ranks are NOT defined because");
print("U_p does not commute with W_p -- report ungraded, plus the partial");
print("grading by the OTHER involution, which U_p does commute with)");
for(t = 1, 2, p = [11,13][t]; \
  print("  p = ", p, "   ungraded rank l_p = ", leak(vector(91,z,z), p)));

print("");
print("  partial grading of im l_11 by W_13  (sectors with W13 = +1 are");
print("  (+,+) and (-,+); with W13 = -1 are (+,-) and (-,-))");
gp13p = concat(idxof(1), idxof(3)); gp13m = concat(idxof(2), idxof(4));
print("    W13 = +1 block: dim ", matrank(Vof(gp13p)), "   leakage ", leak(gp13p, 11));
print("    W13 = -1 block: dim ", matrank(Vof(gp13m)), "   leakage ", leak(gp13m, 11));

print("  partial grading of im l_13 by W_11  (sectors with W11 = +1 are");
print("  (+,+) and (+,-); with W11 = -1 are (-,+) and (-,-))");
gp11p = concat(idxof(1), idxof(2)); gp11m = concat(idxof(3), idxof(4));
print("    W11 = +1 block: dim ", matrank(Vof(gp11p)), "   leakage ", leak(gp11p, 13));
print("    W11 = -1 block: dim ", matrank(Vof(gp11m)), "   leakage ", leak(gp11m, 13));
quit;
