\\ Test Sol's theorem candidate on the DISCOVERY levels (N <= 160):
\\   im( cusp leakage of U_p )  contained in  span{ e_d : d | N, p does not divide d }
\\ For N = pq the right side is 2-dimensional, so this would give
\\ r_cusp(U_p) <= 2 with no census at all.
\\ Retrospective test on discovery data -- labelled as such, not out of sample.
default(parisize, 2^32);
XMAX = 160;
print("N     p  q  |  U_p: rows(d=p, d=pq) zero? surviving rank, dependent?  |  U_q: rows(d=q,d=pq) zero? rank, dep?");
ok = 1; nlev = 0;
forprime(p = 3, 100, forprime(q = p+1, XMAX\p, N = p*q; if(N > XMAX, next); \
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
  C = [v; v*Wp; v*Wq; v*Wp*Wq]; \
  nlev++; \
  Mp = matconcat(C) * mfheckemat(mf4, p) * H0; \
  Mq = matconcat(C) * mfheckemat(mf4, q) * H0; \
  zp = (Mp[2,] == 0) && (Mp[4,] == 0); \
  zq = (Mq[3,] == 0) && (Mq[4,] == 0); \
  rp = matrank(matrix(2, matsize(Mp)[2], r, c, if(r==1, Mp[1,c], Mp[3,c]))); \
  rq = matrank(matrix(2, matsize(Mq)[2], r, c, if(r==1, Mq[1,c], Mq[2,c]))); \
  if(!zp || !zq, ok = 0); \
  print(N, "  ", p, " ", q, "  |  U_", p, ": ", if(zp,"ZERO ok","*** NONZERO ***"), \
        "  rank ", rp, if(rp<2, "  (surviving rows DEPENDENT)", "  (independent)"), \
        "  |  U_", q, ": ", if(zq,"ZERO ok","*** NONZERO ***"), \
        "  rank ", rq, if(rq<2, "  (dep)", "  (indep)"))));
print("");
\\ Vacuity guard: zero valid levels is not a pass. The discovery census is
\\ known to contain 22 valid levels at XMAX = 160; anything else means the
\\ eligibility logic has changed and the run must not be read as support.
print("levels tested: ", nlev);
if(nlev == 0, error("NO TEST CASES -- vacuous, not a pass"));
if(nlev != 22, error("expected 22 valid levels at XMAX = 160, got ", nlev, " -- eligibility logic changed, halt"));
print("vanishing-support prediction holds at every level: ", if(ok, "YES", "NO"));
quit;
