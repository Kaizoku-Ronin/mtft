\\ RT-1B cross-level test — self-auditing version.
\\
\\ SUPERSEDES the version staged 2026-08-18, which (a) printed the false
\\ header "boundary dim is 4 for every squarefree N = p*q", (b) computed the
\\ hyperelliptic failure but proceeded anyway, leaving the gate to human
\\ reading of the output, and (c) never checked rank l_p against the
\\ boundary dimension.  All three are enforced in code below.
\\
\\ TRUE statement:  dim B_N = dim S_4 - (3g-3) = nu_infinity + nu_2 + nu_3.
\\ For squarefree N = p*q, nu_infinity = 4, but nu_2 and nu_3 need not
\\ vanish -- N = 91 has nu_3 = 4 and boundary dimension 8, not 4.

default(parisize, 2^32);
nu2(N) = my(f = factor(N)[,1]~); prod(i=1, #f, 1 + kronecker(-1, f[i]));
nu3(N) = my(f = factor(N)[,1]~); prod(i=1, #f, 1 + kronecker(-3, f[i]));

LEV = [35, 39, 55, 65, 77, 85, 91, 95, 143];
PL  = [[2,3,5,7],[2,5,3,13],[2,3,5,11],[2,5,13],[2,7,11],[2,5,17],[2,3,5,7,13],[2,5,19],[2,3,11,13]];

for(w = 1, #LEV, N = LEV[w]; pl = PL[w]; \
  mf2 = mfinit([N,2],1); g = mfdim(mf2); \
  mf4 = mfinit([N,4],1); d4 = mfdim(mf4); \
  bN = d4 - (3*g-3); \
  if(bN != 4 + nu2(N) + nu3(N), print("N = ", N, "  *** boundary formula MISMATCH"); next); \
  fac = factor(N)[,1]~; idx = N * prod(i=1, #fac, 1 + 1/fac[i]); STB = idx \ 3; \
  DEPTH = STB * vecmax(pl) + 5; \
  A2 = mfcoefs(mf2, DEPTH); \
  S = vector(g, k, sum(n=0, DEPTH, A2[n+1,k]*z^n) + O(z^(DEPTH+1))); \
  pr = List(); for(i=1,g, for(j=i,g, listput(pr, S[i]*S[j]))); pr = Vec(pr); \
  V = matconcat(vector(#pr, k, vectorv(STB+1, n, polcoef(pr[k], n-1)))); \
  d0 = matrank(V); \
  print("N = ", N, "  genus ", g, "  dim S_4 ", d4, "  b_N ", bN, \
        "  (nu_inf 4, nu_2 ", nu2(N), ", nu_3 ", nu3(N), ")"); \
  if(d0 != 3*g-3, \
     print("   rank mu_2 = ", d0, " < 3g-3 = ", 3*g-3, \
           "   ==> INVALID_FOR_CANONICAL_LEAKAGE (hyperelliptic: canonical", \
           " multiplication is not surjective, so the space tested is not H^0(2K))"); \
     print(""); next); \
  for(t = 1, #pl, p = pl[t]; \
    TV = matconcat(vector(#pr, k, vectorv(STB+1, n, \
       if(N % p == 0, polcoef(pr[k], (n-1)*p), \
          polcoef(pr[k], (n-1)*p) + p^3*if((n-1)%p==0, polcoef(pr[k], (n-1)/p), 0))))); \
    lk = matrank(matconcat([V, TV])) - d0; \
    if(lk < 0 || lk > bN, \
       print("   *** GATE VIOLATION: rank l_", p, " = ", lk, " not in [0, ", bN, "]"); \
       print("   *** this is the check that would have caught the 41 > 40 error")); \
    print("     p = ", p, if(N%p==0, "  U_p", "  T_p"), "   rank l_p = ", lk, \
          if(N%p!=0 && lk == bN, "   (saturates)", ""))); \
  print(""));
quit;
