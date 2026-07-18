\\ MTFT Paper 32: Exact Period Matrix of X0(143) - v4
\\ Uses norml2([z]) for absolute value (proven to work in v2)
\\ Roger Tano - April 2026

default(realprecision, 50);
default(parisizemax, 2000000000);

\\ Absolute value via norml2 (avoids t_SER issues)
cabs(z) = { sqrt(norml2([z])); }

\\ Phase in degrees via one-arg atan
myphase(z) = { my(r, i, ph); r = real(z); i = imag(z); ph = if(r == 0, if(i > 0, Pi/2, -Pi/2), atan(i/r) + if(r < 0, if(i >= 0, Pi, -Pi), 0)); 180/Pi * ph; }

print("================================================================");
print("  MTFT PAPER 32: PERIOD MATRIX OF X0(143) - v4");
print("================================================================");

\\ --- 1. Initialize ---
print("\n--- 1. Modular form space ---");
mf = mfinit([143, 2]);
print("  Level 143 = 11 x 13, Weight 2");
print("  Cuspidal dim = ", mfdim(mf, 0));
print("  New dim = ", mfdim(mf, 1));
B = mfeigenbasis(mf);
print("  Orbits: ", #B);
print("  f1 poly: ", B[1].mod, "  dim = ", poldegree(B[1].mod));
print("  f2 poly: ", B[2].mod, "  dim = ", poldegree(B[2].mod));
print("  f3 poly: ", B[3].mod, "  dim = ", poldegree(B[3].mod));

\\ --- 2. q-expansions ---
print("\n--- 2. q-expansion verification ---");
c1 = mfcoefs(B[1], 20); print("  f1: ", c1);
print("  Elliptic curve 143.a1:");
E = ellinit([0, -1, 1, -1, -2]);
print("    Conductor: ", ellglobalred(E)[1]);
print("    a2=", ellap(E,2), " a3=", ellap(E,3), " a5=", ellap(E,5), " a7=", ellap(E,7));

\\ --- 3. Modular symbols ---
print("\n--- 3. Modular symbols ---");
FS1 = mfsymbol(mf, B[1]); print("  f1 done.");
FS2 = mfsymbol(mf, B[2]); print("  f2 done.");
FS3 = mfsymbol(mf, B[3]); print("  f3 done.");

\\ --- 4. Petersson norms ---
print("\n--- 4. Petersson norms ---");
pet1 = mfpetersson(FS1); print("  <f1,f1> = ", pet1);

\\ --- 5. Basic periods ---
print("\n--- 5. Basic periods {0,inf} ---");
print("  f1: ", mfsymboleval(FS1, [0,1;1,0]));
print("  f2: ", mfsymboleval(FS2, [0,1;1,0]));
print("  f3: ", mfsymboleval(FS3, [0,1;1,0]));

\\ ================================================================
\\ MASS CYCLE {inf, 2/77}
\\ ================================================================
print("\n================================================================");
print("  MASS CYCLE: {inf, 2/77}");
print("================================================================");

path_mc = [2, 1; 77, 0];
target = log(1776.86/0.511) / log(105.658/0.511);
print("  Target d3/d2 = ", target);

\\ --- Evaluate periods ---
p1 = mfsymboleval(FS1, path_mc);
p2 = mfsymboleval(FS2, path_mc);
p3 = mfsymboleval(FS3, path_mc);

print("\n  f1 period: ", p1);
print("  type(f1) = ", type(p1));
print("  f2 period: ", p2);
print("  type(f2) = ", type(p2), "  #components = ", if(type(p2)=="t_VEC", #p2, 1));
print("  f3 period: ", p3);
print("  type(f3) = ", type(p3), "  #components = ", if(type(p3)=="t_VEC", #p3, 1));

\\ --- Absolute values ---
omega1 = cabs(p1);
print("\n  |Omega_1| = ", omega1);

\\ f2 per-embedding
print("\n  f2 per-embedding:");
for(k = 1, #p2, print("    sigma_", k, ": ", p2[k], "  |.| = ", cabs(p2[k]), "  phase = ", myphase(p2[k]), " deg"));

\\ f3 per-embedding
print("\n  f3 per-embedding:");
for(k = 1, #p3, print("    sigma_", k, ": ", p3[k], "  |.| = ", cabs(p3[k]), "  phase = ", myphase(p3[k]), " deg"));

\\ --- Norms ---
print("\n--- Norms ---");
tr2 = vecsum(p2); tr3 = vecsum(p3);
norm_tr2 = cabs(tr2);
norm_tr3 = cabs(tr3);
print("  Tr(Om2) = ", tr2);
print("  |Tr(Om2)| = ", norm_tr2);
print("  Tr(Om3) = ", tr3);
print("  |Tr(Om3)| = ", norm_tr3);

norm_l2_2 = sqrt(norml2(p2));
norm_l2_3 = sqrt(norml2(p3));
print("  L2(Om2) = ", norm_l2_2);
print("  L2(Om3) = ", norm_l2_3);

norm_l1_2 = sum(k=1, #p2, cabs(p2[k]));
norm_l1_3 = sum(k=1, #p3, cabs(p3[k]));
print("  L1(Om2) = ", norm_l1_2);
print("  L1(Om3) = ", norm_l1_3);

\\ --- Depth ratios ---
print("\n--- Depth ratios ---");

d2_tr = log(norm_tr2 / omega1);
d3_tr = log(norm_tr3 / omega1);
r_tr = d3_tr / d2_tr;
print("  Trace:  d3/d2 = ", r_tr, "  delta = ", r_tr - target);

d2_l2 = log(norm_l2_2 / omega1);
d3_l2 = log(norm_l2_3 / omega1);
r_l2 = d3_l2 / d2_l2;
print("  L2:     d3/d2 = ", r_l2, "  delta = ", r_l2 - target);

d2_l1 = log(norm_l1_2 / omega1);
d3_l1 = log(norm_l1_3 / omega1);
r_l1 = d3_l1 / d2_l1;
print("  L1:     d3/d2 = ", r_l1, "  delta = ", r_l1 - target);

\\ ================================================================
\\ Q2/Q4 HECKE FACTORIZATION
\\ ================================================================
print("\n================================================================");
print("  Q2/Q4 HECKE FACTORIZATION");
print("================================================================");

\\ Identify Q2 (phase near +105) vs Q4 (rest) by phase
print("\n  Phase structure of f3:");
for(k = 1, #p3, print("    sigma_", k, ": phase = ", myphase(p3[k]), " deg  |.| = ", cabs(p3[k])));

\\ Q2 = sigma_1, sigma_4 (phases near +105 deg)
\\ Q4 = sigma_2, sigma_3, sigma_5, sigma_6
q4_re = real(p3[2]) + real(p3[3]) + real(p3[5]) + real(p3[6]);
q4_im = imag(p3[2]) + imag(p3[3]) + imag(p3[5]) + imag(p3[6]);
abs_tr_q4 = sqrt(q4_re^2 + q4_im^2);
q4_l2 = sqrt(norml2([p3[2], p3[3], p3[5], p3[6]]));

q2_re = real(p3[1]) + real(p3[4]);
q2_im = imag(p3[1]) + imag(p3[4]);
abs_tr_q2 = sqrt(q2_re^2 + q2_im^2);
q2_l2 = sqrt(norml2([p3[1], p3[4]]));

print("\n  Q2 (barrier fold, dropped):");
print("    sigma_1 + sigma_4");
print("    |Tr_Q2| = ", abs_tr_q2);
print("    L2_Q2 = ", q2_l2);

print("\n  Q4 (attractor + reflected, kept):");
print("    sigma_2 + sigma_3 + sigma_5 + sigma_6");
print("    |Tr_Q4| = ", abs_tr_q4);
print("    L2_Q4 = ", q4_l2);

print("\n  --- Q4-restricted depth ratios ---");

d3_q4_l2 = log(q4_l2 / omega1);

r_q4_l2_tr = d3_q4_l2 / d2_tr;
print("  Q4_L2 / f2_Trace:  d3/d2 = ", r_q4_l2_tr, "  delta = ", r_q4_l2_tr - target);

r_q4_l2_l2 = d3_q4_l2 / d2_l2;
print("  Q4_L2 / f2_L2:     d3/d2 = ", r_q4_l2_l2, "  delta = ", r_q4_l2_l2 - target, "  <-- PAPER 32");

d3_q4_tr = log(abs_tr_q4 / omega1);
r_q4_tr_tr = d3_q4_tr / d2_tr;
print("  Q4_Tr / f2_Trace:  d3/d2 = ", r_q4_tr_tr, "  delta = ", r_q4_tr_tr - target);

print("\n  Dimension equalization: dim_eff(f3) = 4 = dim(f2) = 4");

\\ ================================================================
\\ STIFFNESS-PERIOD FORMULA
\\ ================================================================
print("\n================================================================");
print("  STIFFNESS-PERIOD FORMULA");
print("================================================================");

lam1 = 1.52655721; lam2 = 6.42733234; lam3 = 45.46811045;
m_e = 0.51100; m_mu = 105.658; m_tau = 1776.86;

ln_mu = log(m_mu / m_e);
ln_tau = log(m_tau / m_e);
dl2 = log(lam2/lam1);
dl3 = log(lam3/lam1);

dp2 = log(norm_tr2 / omega1);
dp3 = log(norm_tr3 / omega1);

print("  ln(m_mu/m_e) = ", ln_mu);
print("  ln(m_tau/m_e) = ", ln_tau);
print("  ln(lam2/lam1) = ", dl2);
print("  ln(lam3/lam1) = ", dl3);
print("  dp2 = ln(|Tr(Om2)|/|Om1|) = ", dp2);
print("  dp3 = ln(|Tr(Om3)|/|Om1|) = ", dp3);

det_A = dl2*dp3 - dl3*dp2;
alpha_sol = (ln_mu*dp3 - ln_tau*dp2) / det_A;
beta_sol = (dl2*ln_tau - dl3*ln_mu) / det_A;
print("  alpha = ", alpha_sol);
print("  beta = ", beta_sol);

pred_mu = m_e * (lam2/lam1)^alpha_sol * (norm_tr2/omega1)^beta_sol;
pred_tau = m_e * (lam3/lam1)^alpha_sol * (norm_tr3/omega1)^beta_sol;
print("  m_mu predicted = ", pred_mu, " MeV (obs: ", m_mu, ")");
print("  m_tau predicted = ", pred_tau, " MeV (obs: ", m_tau, ")");

\\ ================================================================
\\ MASS CYCLE SCAN (c = k*11, a = 2)
\\ ================================================================
print("\n================================================================");
print("  MASS CYCLE SCAN: a=2, c=k*11 for k=1..13");
print("================================================================");

scan_cycle(FS1, FS2, FS3, a, c, omega1, target) = { my(pp1, pp2, pp3, a1, a2, a3, dd2, dd3, rr); pp1 = mfsymboleval(FS1, [a,1;c,0]); pp2 = mfsymboleval(FS2, [a,1;c,0]); pp3 = mfsymboleval(FS3, [a,1;c,0]); a1 = cabs(pp1); a2 = cabs(vecsum(pp2)); a3 = cabs(vecsum(pp3)); if(a1 < 1e-10 || a2 < 1e-10, print("  ", a, "/", c, ": SKIP (near zero)"); return(0)); dd2 = log(a2/a1); dd3 = log(a3/a1); rr = dd3/dd2; print("  ", a, "/", c, "  |f1|=", a1, "  |Tr(f2)|=", a2, "  |Tr(f3)|=", a3, "  d3/d2=", rr, "  delta=", rr-target); rr; }

r = scan_cycle(FS1, FS2, FS3, 2, 11, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 22, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 33, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 44, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 55, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 66, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 77, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 88, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 99, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 110, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 121, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 132, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 2, 143, omega1, target);

\\ Also scan a=1 for comparison
print("\n  --- a=1 comparison ---");
r = scan_cycle(FS1, FS2, FS3, 1, 11, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 1, 77, omega1, target);
r = scan_cycle(FS1, FS2, FS3, 1, 143, omega1, target);

print("\n================================================================");
print("  COMPUTATION COMPLETE");
print("================================================================");
