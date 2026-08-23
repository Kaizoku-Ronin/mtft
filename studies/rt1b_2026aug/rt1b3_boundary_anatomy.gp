\\\ RT-1B.3 Boundary Anatomy. Splits B_N = S_4/H^0(2K) into cusp and elliptic
\\\\ channels via the filtration H^0(2K) < S_4^(2) < S_4, where S_4^(2) is cut
\\\\ by a_1(f|W) = 0 over the four Atkin-Lehner elements (valid because AL acts
\\\\ freely transitively on the four cusps of squarefree N = pq).
\\\\ Uses PARI mfheckemat/mfatkininit -- an INDEPENDENT operator implementation
\\\\ from the q-expansion convolutions in rt1b_leakage_sweep.gp.

default(parisize, 2^32);
nu2(N) = my(f = factor(N)[,1]~); prod(i=1, #f, 1 + kronecker(-1, f[i]));
nu3(N) = my(f = factor(N)[,1]~); prod(i=1, #f, 1 + kronecker(-3, f[i]));
LEV = [55, 65, 77, 85, 91, 95, 143];
PL  = [[2,5,11],[2,5,13],[2,7,11],[2,5,17],[2,7,13],[2,5,19],[2,11,13]];
for(w = 1, #LEV, N = LEV[w]; pl = PL[w]; fac = factor(N)[,1]~; \
  mf2 = mfinit([N,2],1); g = mfdim(mf2); \
  mf4 = mfinit([N,4],1); d4 = mfdim(mf4); b = d4 - (3*g-3); \
  print("N = ", N, "  genus ", g, "  dim S_4 = ", d4, "  dim H0(2K) = ", 3*g-3, "  boundary b_N = ", b); \
  print("   nu_inf = 4  nu_2 = ", nu2(N), "  nu_3 = ", nu3(N), "  predicted b_N = ", 4+nu2(N)+nu3(N), if(b == 4+nu2(N)+nu3(N), "   OK", "   MISMATCH")); \
  D = 3*d4 + 30; A2 = mfcoefs(mf2, D); A4 = mfcoefs(mf4, D); \
  np = g*(g+1)/2; prodm = matrix(D+1, np); kk = 0; \
  for(i = 1, g, for(j = i, g, kk++; for(n = 1, D+1, prodm[n,kk] = sum(a=1, n, A2[a,i]*A2[n+1-a,j])))); \
  M4 = matrix(D+1, d4, r, c, A4[r,c]); \
  H0 = matsolve(mattranspose(M4)*M4, mattranspose(M4)*prodm); \
  print("   H0(2K) in mf4 coords: rank ", matrank(H0), "   [", 3*g-3, "]"); \
  aw = mfatkininit(mf4, fac[1]); W1 = aw[2]/aw[3]; \
  aw = mfatkininit(mf4, fac[2]); W2 = aw[2]/aw[3]; \
  v = matrix(1, d4, r, c, A4[2,c]); \
  F = matconcat([v; v*W1; v*W2; v*W1*W2]); \
  S42 = matker(F); ds = matsize(S42)[2]; \
  print("   S_4^(2) dim ", ds, "   cusp channels ", d4-ds, "   elliptic channels ", ds-(3*g-3)); \
  print("      p     total   cusp   elliptic"); \
  for(t = 1, #pl, p = pl[t]; TH = mfheckemat(mf4, p) * H0; \
    tot = matrank(matconcat([H0, TH])) - (3*g-3); \
    cu = matrank(matconcat([S42, TH])) - ds; \
    if(tot > b, print("   *** VIOLATION rank ", tot, " > b_N ", b)); \
    print("     ", if(p<10," ",""), p, if(N%p==0, "  U", "  T"), "      ", tot, "       ", cu, "       ", tot-cu)); \
  print(""));
quit;
