default(parisize, 2^32);
nu2(N) = my(f = factor(N)[,1]~); prod(i=1, #f, 1 + kronecker(-1, f[i]));
nu3(N) = my(f = factor(N)[,1]~); prod(i=1, #f, 1 + kronecker(-3, f[i]));
{
signature(mf4, ell, H0, S42, dH, dS, nell, lab) =
  my(TH, tt, cc, re);
  TH = mfheckemat(mf4, ell) * H0;
  tt = matrank(matconcat([H0, TH])) - dH;
  cc = matrank(matconcat([S42, TH])) - dS;
  re = tt - cc;
  \\ GATE 4, as preregistered: every signature, both components.
  if(cc < 0 || cc > 4,
     error("GATE 4 VIOLATION ", lab, ": r_cusp = ", cc, " not in [0,4]"));
  if(re < 0 || re > nell,
     error("GATE 4 VIOLATION ", lab, ": r_ell = ", re, " not in [0,", nell, "]"));
  [cc, re];
}
{
goodprimes(N, howmany) =
  my(L = List(), ll = 2);
  while(#L < howmany, if(N % ll != 0, listput(L, ll)); ll = nextprime(ll+1));
  Vec(L);
}
XMAX = 160;
print("N    p   q   g  b_N (4,n2,n3)   good l   signatures            U_p     U_q");
nvalid = 0; ngood = 0; q1 = 1; q2 = 1; q3 = 1; q4 = 1; viol = 0; badq = List();
forprime(p = 3, 100, forprime(q = p+1, XMAX\p, N = p*q; if(N > XMAX, next); \
  mf2 = mfinit([N,2],1); g = mfdim(mf2); if(g < 3, next); \
  mf4 = mfinit([N,4],1); d4 = mfdim(mf4); bN = d4 - (3*g-3); \
  if(bN != 4 + nu2(N) + nu3(N), error("GATE 1 VIOLATION N=", N, ": b_N = ", bN, " != 4 + nu_2 + nu_3 = ", 4+nu2(N)+nu3(N), " -- halt and re-derive")); \
  idx = N*(1+1/p)*(1+1/q); STB = idx\3; D = STB + 12; \
  A2 = mfcoefs(mf2, D); A4 = mfcoefs(mf4, D); \
  SS = vector(g, k, sum(n=0, D, A2[n+1,k]*z^n) + O(z^(D+1))); \
  pr = List(); for(i=1,g, for(j=i,g, listput(pr, SS[i]*SS[j]))); pr = Vec(pr); \
  PM = matconcat(vector(#pr, k, vectorv(D+1, n, polcoef(pr[k], n-1)))); \
  d0 = matrank(PM); \
  if(d0 != 3*g-3, print(N, "  ", p, " ", q, "  g=", g, "  INVALID_FOR_CANONICAL_LEAKAGE (rank mu2 = ", d0, " < ", 3*g-3, ")"); next); \
  M4 = matrix(D+1, d4, r, c, A4[r,c]); \
  H0 = matsolve(mattranspose(M4)*M4, mattranspose(M4)*PM); \
  aw = mfatkininit(mf4, p); W1 = aw[2]/aw[3]; aw = mfatkininit(mf4, q); W2 = aw[2]/aw[3]; \
  S42 = matker(matconcat([v0 = matrix(1, d4, r, c, A4[2,c]); v0*W1; v0*W2; v0*W1*W2])); \
  ds = matsize(S42)[2]; \
  if(ds != d4 - 4, error("GATE 3 VIOLATION N=", N, ": dim S_4^(2) = ", ds, " != d4 - 4 = ", d4-4, " -- halt and re-derive")); \
  nvalid++; \
  gl = goodprimes(N, 3); \
  nell = nu2(N)+nu3(N); \
  gs = vector(3, t, signature(mf4, gl[t], H0, S42, 3*g-3, ds, nell, Str("N=",N," T_",gl[t]))); \
  for(t = 1, 3, ngood++; if(gs[t] != [4, nell], q1 = 0)); \
  sp = signature(mf4, p, H0, S42, 3*g-3, ds, nell, Str("N=",N," U_",p)); \
  sq = signature(mf4, q, H0, S42, 3*g-3, ds, nell, Str("N=",N," U_",q)); \
  listput(badq, [N, p, q, nu2(N)+nu3(N), sp, sq]); \
  if(sq[1] != 1, q2 = 0); \
  if(nu2(N)+nu3(N) > 0 && !(sp[2] > 0 && sq[2] == 0), q3 = 0); \
  if(sp[2] > 1 || sq[2] > 1, q4 = 0); \
  print(N, "  ", p, " ", q, "  ", g, "  ", bN, "  (4,", nu2(N), ",", nu3(N), ")  ", gl, " ", gs, "  ", sp, "  ", sq)));
print("");
print("valid levels ", nvalid, "   good (N,l) pairs ", ngood);
print("gates 1, 3, 4 are hard errors: reaching this line means all passed");
print("Q1 good primes saturate BOTH channels        : ", if(q1, "HOLDS", "FALSIFIED"));
print("Q2 r_cusp(U_q) = 1 for the larger level prime: ", if(q2, "HOLDS", "FALSIFIED"));
print("Q3 smaller level prime reaches elliptic      : ", if(q3, "HOLDS", "FALSIFIED"));
print("Q4 r_ell <= 1 for both level primes          : ", if(q4, "HOLDS", "FALSIFIED"));
quit;
