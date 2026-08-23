\\ RT-1B.5 — Cusp Leakage Matrices at N = 143.
\\ Replaces rank statistics with the actual map, in the AL character basis.
default(parisize, 2^32);
N = 143; p = 11; q = 13;
mf2 = mfinit([N,2],1); g = mfdim(mf2);
mf4 = mfinit([N,4],1); d4 = mfdim(mf4);
idx = N*(1+1/p)*(1+1/q); STB = idx\3; D = STB + 60;
A2 = mfcoefs(mf2, D); A4 = mfcoefs(mf4, D);
SS = vector(g, k, sum(n=0, D, A2[n+1,k]*z^n) + O(z^(D+1)));
pr = List(); for(i=1,g, for(j=i,g, listput(pr, SS[i]*SS[j]))); pr = Vec(pr);
PM = matconcat(vector(#pr, k, vectorv(D+1, n, polcoef(pr[k], n-1))));
M4 = matrix(D+1, d4, r, c, A4[r,c]);
H0full = matsolve(mattranspose(M4)*M4, mattranspose(M4)*PM);
H0 = H0full;
H0 = matimage(H0);
print("dim S_4 = ", d4, "   dim H^0(2K) = ", matsize(H0)[2], "   [36]");
aw = mfatkininit(mf4, p); W11 = aw[2]/aw[3];
aw = mfatkininit(mf4, q); W13 = aw[2]/aw[3];
v = matrix(1, d4, r, c, A4[2,c]);          \\ a_1 functional
C = [v; v*W11; v*W13; v*W11*W13];          \\ c_1, c_11, c_13, c_143
print("all four cusp functionals vanish on H^0(2K): ", matconcat(C)*H0 == 0);
print("");
print("VERIFY Sol's identity  c_d(T_l F) = a_l(F | W_d):");
print("  equivalently  v . T_l = (a_l functional), since T_l commutes with W_d");
forprime(l = 2, 13, if(N % l != 0, \
  print("   l = ", l, " :  v.T_l == a_l row ?  ", \
        v * mfheckemat(mf4, l) == matrix(1, d4, r, c, A4[l+1,c]))));
print("   a_1(U_p F) = a_p(F) holds for the level primes too -- the identity");
print("   itself is NOT where good and bad primes differ. The difference is");
print("   that T_l commutes with W_d and U_p does not, so only for good l can");
print("   the identity be transported to the other cusps as a_l(F | W_d).");
forprime(l = 11, 13, if(N % l == 0, \
  print("   p = ", l, " :  v.U_p == a_p row ?  ", \
        v * mfheckemat(mf4, l) == matrix(1, d4, r, c, A4[l+1,c]))));
print("");
\\ character basis of C^4 indexed by cusps (1, 11, 13, 143)
U = [1,1,1,1; 1,1,-1,-1; 1,-1,1,-1; 1,-1,-1,1];   \\ rows: ++, +-, -+, --
nm = ["(+,+)","(+,-)","(-,+)","(-,-)"];
print("cusp-leakage image, expressed in the AL character basis");
print("   operator   rank   image described by");
{
report(A, lab) =
  my(Mm, img, k, co);
  Mm = matconcat([v*A; v*W11*A; v*W13*A; v*W11*W13*A]) * H0;
  img = matimage(Mm);
  k = matsize(img)[2];
  co = matsolve(mattranspose(U), img);   \\ coordinates in the character basis
  print("   ", lab, "      ", k);
  for(j = 1, k, print("        component ", j, " : ",
      concat(vector(4, t, Str(nm[t], "=", co[t,j], "  ")))));
}
forprime(l = 2, 7, if(N % l != 0, report(mfheckemat(mf4, l), Str("T_", l))));
report(mfheckemat(mf4, 11), "U_11");
report(mfheckemat(mf4, 13), "U_13");

print("");
print("HARD GATES on the RT-1B.5 relations (normalization: c_d(F) = a_1(F|W_d)");
print("with W_d from mfatkininit; vanishing SUPPORT and DIMENSIONS are");
print("normalization-invariant, the coefficient 1 in c_1 = c_11 is not)");
{
gate(A, lab, rk, zrows, extra) =
  my(M, n, sub);
  M = matconcat([v*A; v*W11*A; v*W13*A; v*W11*W13*A]) * H0;
  if(matrank(M) != rk, error(lab, ": rank ", matrank(M), " != ", rk));
  for(i = 1, #zrows,
    if(M[zrows[i],] != 0, error(lab, ": row d-index ", zrows[i], " not zero")));
  if(M == 0, error(lab, ": image is zero -- no witness"));
  if(extra && M[1,] - M[2,] != 0, error(lab, ": c_1 - c_11 not identically zero"));
  print("   ", lab, "  rank ", rk, "  zero rows ", zrows,
        if(extra, "  and c_1 = c_11", ""), "   PASS");
}
gate(mfheckemat(mf4, 11), "U_11", 2, [2,4], 0);
gate(mfheckemat(mf4, 13), "U_13", 1, [3,4], 1);
forprime(l = 2, 7, if(N % l != 0, gate(mfheckemat(mf4, l), Str("T_", l), 4, [], 0)));
quit;
