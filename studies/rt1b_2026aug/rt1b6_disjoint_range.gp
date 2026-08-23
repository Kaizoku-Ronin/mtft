\\ RT-1B.6 — FROZEN RUNNER for the disjoint holdout range.
\\ Tests ONLY H1, Q1b, Q3.  Q2 and Q4 are dead and are not computed here.
\\ Nothing outside 160 < N <= 400 is reported.  Gates 1, 3, 4 are hard errors.
\\ Do not edit after freezing; the discovery census lives in
\\ rt1b4_boundary_census.gp and is untouched.
default(parisize, 2^32);
NMIN = 160; NMAX = 400;
LBOUND = 20;      \\ Q1b: ALL good primes l <= LBOUND, an intrinsic test set

nu2(N) = my(f = factor(N)[,1]~); prod(i=1, #f, 1 + kronecker(-1, f[i]));
nu3(N) = my(f = factor(N)[,1]~); prod(i=1, #f, 1 + kronecker(-3, f[i]));
{
signature(mf4, ell, H0, S42, dH, dS, nell, lab) =
  my(TH, tt, cc, re);
  TH = mfheckemat(mf4, ell) * H0;
  tt = matrank(matconcat([H0, TH])) - dH;
  cc = matrank(matconcat([S42, TH])) - dS;
  re = tt - cc;
  if(cc < 0 || cc > 4, error("GATE4 ", lab, ": r_cusp = ", cc));
  if(re < 0 || re > nell, error("GATE4 ", lab, ": r_ell = ", re, " vs ", nell));
  [cc, re];
}
{
goodlist(N, L) =
  my(R = List());
  forprime(l = 2, L, if(N % l != 0, listput(R, l)));
  Vec(R);
}
print("RT-1B.6 holdout run.  ", NMIN, " < N <= ", NMAX, "   Q1b bound L = ", LBOUND);
print("N    p   q   g  b_N (4,n2,n3)   Q1b   U_p      U_q     H1(U_q invariant)");
nlev = 0; ngood = 0; h1n = 0; h1ok = 1; q1 = 1; q3 = 1; q3n = 0;
forprime(p = 3, 200, forprime(q = p+1, NMAX\p, N = p*q; \
  if(N <= NMIN || N > NMAX, next); \
  mf2 = mfinit([N,2],1); g = mfdim(mf2); if(g < 3, next); \
  mf4 = mfinit([N,4],1); d4 = mfdim(mf4); bN = d4 - (3*g-3); nell = nu2(N)+nu3(N); \
  if(bN != 4 + nell, error("GATE1 N=", N, ": b_N = ", bN, " != ", 4+nell)); \
  idx = N*(1+1/p)*(1+1/q); STB = idx\3; D = STB + 12; \
  A2 = mfcoefs(mf2, D); A4 = mfcoefs(mf4, D); \
  SS = vector(g, k, sum(n=0, D, A2[n+1,k]*z^n) + O(z^(D+1))); \
  pr = List(); for(i=1,g, for(j=i,g, listput(pr, SS[i]*SS[j]))); pr = Vec(pr); \
  PM = matconcat(vector(#pr, k, vectorv(D+1, n, polcoef(pr[k], n-1)))); \
  if(matrank(PM) != 3*g-3, print(N, "  ", p, " ", q, "   INVALID_FOR_CANONICAL_LEAKAGE (mu_2 not surjective)"); next); \
  M4 = matrix(D+1, d4, r, c, A4[r,c]); \
  H0 = matsolve(mattranspose(M4)*M4, mattranspose(M4)*PM); \
  aw = mfatkininit(mf4, p); Wp = aw[2]/aw[3]; aw = mfatkininit(mf4, q); Wq = aw[2]/aw[3]; \
  v0 = matrix(1, d4, r, c, A4[2,c]); \
  S42 = matker(matconcat([v0; v0*Wp; v0*Wq; v0*Wp*Wq])); ds = matsize(S42)[2]; \
  if(ds != d4 - 4, error("GATE3 N=", N, ": dim S4^(2) = ", ds, " != ", d4-4)); \
  nlev++; \
  gl = goodlist(N, LBOUND); \
  allsat = 1; \
  for(t = 1, #gl, ngood++; \
      if(signature(mf4, gl[t], H0, S42, 3*g-3, ds, nell, Str("N=",N," T_",gl[t])) != [4, nell], allsat = 0; q1 = 0)); \
  sp = signature(mf4, p, H0, S42, 3*g-3, ds, nell, Str("N=",N," U_",p)); \
  sq = signature(mf4, q, H0, S42, 3*g-3, ds, nell, Str("N=",N," U_",q)); \
  if(nell > 0, q3n++; if(!(sp[2] > 0 && sq[2] == 0), q3 = 0)); \
  h1 = -1; \
  if(p == 3, h1n++; h1 = (sq == [0,0]); if(!h1, h1ok = 0)); \
  h1s = if(h1 < 0, "n/a", if(h1, "HOLDS", "*** FAILS ***")); \
  print(N, "  ", p, " ", q, "  ", g, "  ", bN, "  (4,", nu2(N), ",", nu3(N), ")  ", \
        if(allsat, "SAT", "*FAIL*"), "  ", sp, "  ", sq, "   ", h1s)));
print("");
print("valid levels ", nlev, "   good (N,l) pairs ", ngood, "   N=3q levels ", h1n, "   elliptic levels ", q3n);
print("gates 1,3,4 are hard errors: reaching this line means all passed");
\\ Guard against vacuous truth: a hypothesis with no test cases is NOT a pass.
verdict(flag, n, lab) = print(lab, if(n == 0, "  NO TEST CASES -- vacuous, not a pass", if(flag, "  HOLDS", "  FALSIFIED")), "   (n = ", n, ")");
verdict(h1ok, h1n, "H1  U_q H^0(2K) inside H^0(2K), N = 3q  :");
verdict(q1, ngood, "Q1b every good l <= L gives (4,nu2+nu3) :");
verdict(q3, q3n, "Q3  smaller prime reaches elliptic       :");
quit;
