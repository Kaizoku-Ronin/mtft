\\ Second theorem target: the surviving pair for the LARGER level prime is
\\ always dependent.  What is the constant?  At N = 143 it was c_1 = c_11
\\ exactly.  Discovery levels only (N <= 160), p >= 5 (for p = 3 the whole
\\ U_q signature is zero, so there is no pair to relate).
\\ NOTE: the constant is tied to the mfatkininit normalization of W_d;
\\ comparing across levels in ONE normalization is what is meaningful here.
default(parisize, 2^32);
{
ratio(r1, r2) =
  my(j = 0, lam);
  for(k = 1, #r1, if(r1[k] != 0, j = k; break));
  if(j == 0, return("row1 is zero"));
  lam = r2[j]/r1[j];
  if(r2 - lam*r1 != 0, return("NOT proportional"));
  lam;
}
print("N     p   q   | U_q: surviving rows d=1, d=p    | U_p pair (if dependent)");
forprime(p = 5, 100, forprime(q = p+1, 160\p, N = p*q; if(N > 160, next); \
  mf2 = mfinit([N,2],1); g = mfdim(mf2); if(g < 3, next); \
  mf4 = mfinit([N,4],1); d4 = mfdim(mf4); \
  idx = N*(1+1/p)*(1+1/q); STB = idx\3; D = STB + 12; \
  A2 = mfcoefs(mf2, D); A4 = mfcoefs(mf4, D); \
  SS = vector(g, k, sum(n=0, D, A2[n+1,k]*z^n) + O(z^(D+1))); \
  pr = List(); for(i=1,g, for(j=i,g, listput(pr, SS[i]*SS[j]))); pr = Vec(pr); \
  PM = matconcat(vector(#pr, k, vectorv(D+1, n, polcoef(pr[k], n-1)))); \
  if(matrank(PM) != 3*g-3, next); \
  M4 = matrix(D+1, d4, r, c, A4[r,c]); \
  H0 = matsolve(mattranspose(M4)*M4, mattranspose(M4)*PM); \
  aw = mfatkininit(mf4, p); Wp = aw[2]/aw[3]; aw = mfatkininit(mf4, q); Wq = aw[2]/aw[3]; \
  v = matrix(1, d4, r, c, A4[2,c]); \
  C = matconcat([v; v*Wp; v*Wq; v*Wp*Wq]); \
  Mq = C * mfheckemat(mf4, q) * H0; \
  Mp = C * mfheckemat(mf4, p) * H0; \
  lq = ratio(Mq[1,], Mq[2,]); \
  rp = matrank(matrix(2, matsize(Mp)[2], r, c, if(r==1, Mp[1,c], Mp[3,c]))); \
  lp = if(rp == 1, Str("c_q = ", ratio(Mp[1,], Mp[3,]), " * c_1"), "independent (rank 2)"); \
  print(N, "  ", p, " ", q, "   | c_", p, " = ", lq, " * c_1", "   | ", lp)));
quit;
