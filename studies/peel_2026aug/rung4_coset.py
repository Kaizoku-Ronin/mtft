"""
PART A: Rung 4 mass gap in the corrected two-temperature ensemble language.

The Dirichlet ensemble factorizes EXACTLY (log = Lambda * 1):
   Z_D(beta) = sum w_n n^-beta = zeta(beta) * zeta(beta+1) * sum_m Lambda(m) m^-(beta+1)
             = [gas at beta] x [gas at beta+1] x [ONE marked excitation at beta+1]
The n^2-tilted (stiffness) ensemble is Z_D(s-2); its dominant pole: s = 3, residue T_inf.

CLAIM: the Rung-4 strong-coupling transfer spectrum is the single-marked-primon
Boltzmann spectrum at the stiffness condensation point beta_c = 3:
   lambda_m = Lambda(m) m^-3,   eps_m = -ln lambda_m = 3 ln m - ln ln p
   m_inf = eps_3 - eps_2 = ln( 27 ln2 / (8 ln3) )   -- exactly.

PART B: coset-layer representation bookkeeping for PSL(2,11) x PSL(2,13).
"""
from mpmath import mp, mpf, log, exp, pi, zeta, diff, nstr
mp.dps = 30

# ---------- A1: E2 the factorization ----------
Nw = 4000
lw = [mpf(0)]*(Nw+1)
for d in range(2, Nw+1):
    c = log(mpf(d))/mpf(d)
    for n in range(d, Nw+1, d):
        lw[n] += c
def lam_sum(s, M=200000):
    # sum Lambda(m) m^-s via -zeta'/zeta  (closed) -- use closed form as one route
    return -diff(zeta, mpf(s))/zeta(mpf(s))
print("="*74)
print("  A1. Two-temperature factorization  Z_D = zeta(b) zeta(b+1) [Lambda at b+1]")
print("="*74)
for b in [mpf('2.5'), mpf('3.5')]:
    direct = sum(lw[n]*mpf(n)**(-b) for n in range(1, Nw+1))
    route1 = -zeta(b)*diff(zeta, b+1)                      # -zeta(b) zeta'(b+1)
    route2 = zeta(b)*zeta(b+1)*lam_sum(b+1)                # 3-factor form
    d1 = float(-mp.log10(abs(direct-route1)/abs(route1)))
    d2 = float(-mp.log10(abs(route1-route2)/abs(route2)))
    print(f"  beta={float(b)}:  direct vs -zeta zeta' : {d1:.1f} digits"
          f"   |   -zeta zeta' vs 3-factor : {d2:.1f} digits")

# ---------- A2: the gap identity ----------
print("\n" + "="*74)
print("  A2. m_inf as marked-level spacing at beta_c = 3 (the T_inf pole)")
print("="*74)
LN2, LN3 = log(mpf(2)), log(mpf(3))
m_inf = log(27*LN2/(8*LN3))
alt   = 3*log(mpf(3)/2) - log(LN3/LN2)
print(f"  m_inf = ln(27 ln2 / 8 ln3) = {nstr(m_inf, 22)}")
print(f"  3 ln(3/2) - ln(log_2 3)    = {nstr(alt, 22)}   [identical: "
      f"{nstr(abs(m_inf-alt),3)}]")
lam2, lam3 = LN2/8, LN3/27
print(f"  lambda_2 = Lambda(2) 2^-3 = ln2/8  = {nstr(lam2,12)}")
print(f"  lambda_3 = Lambda(3) 3^-3 = ln3/27 = {nstr(lam3,12)}")
print(f"  ln(lambda_2/lambda_3) = {nstr(log(lam2/lam3),22)}  = m_inf  EXACTLY")
Rinf = m_inf/(1 - exp(-m_inf))
print(f"\n  R(inf) = m_inf/(1-e^-m_inf) = {nstr(Rinf, 15)}"
      f"   (ledger: 1.42507723746)")
print(f"  e^-m_inf = 8 ln3/(27 ln2)  = {nstr(exp(-m_inf),12)}"
      f"  = lambda_3/lambda_2 (transfer-eigenvalue ratio)")

# ---------- A3: the falsifiable ladder ----------
print("\n" + "="*74)
print("  A3. Predicted transfer spectrum at beta_c = 3:  eps_m = 3 ln m - ln Lambda(m)")
print("="*74)
marks = [2,3,4,5,7,8,9,11,13]
def Lam_of(m):
    for p in [2,3,5,7,11,13]:
        q = p
        while q <= m:
            if q == m: return log(mpf(p))
            q *= p
    return None
levels = sorted((3*log(mpf(m)) - log(Lam_of(m)), m) for m in marks)
print(f"  {'rank':>5} {'mark m':>7} {'eps_m':>14} {'gap to previous':>16}")
prev = None
for k,(e,m) in enumerate(levels):
    gap = "" if prev is None else nstr(e-prev, 8)
    print(f"  {k:>5} {m:>7} {nstr(e,10):>14} {gap:>16}")
    prev = e
print("\n  ORDERING PREDICTION: 2, 3, 5, 4, 7, 11, 9, 8, 13")
print("  -> p=11 sits BELOW the prime powers 9 and 8. Distinctive; auditable")
print("     against the chain spectrum in studies/ (PR-5.1 record).")

print("\n  Temperature structure of the gap:  Delta(sigma) = sigma ln(3/2) - ln(log_2 3)")
for s in [2, 3]:
    d = s*log(mpf(3)/2) - log(LN3/LN2)
    tag = "  <-- m_inf (Rung 4)" if s==3 else "  (marked sector AT the pole factorization)"
    print(f"    Delta({s}) = {nstr(d,12)}{tag}")
print(f"    identity: Delta(3) = Delta(2) + ln(3/2)   [{nstr(2*log(mpf(3)/2)-log(LN3/LN2)+log(mpf(3)/2),12)}]")
sig_star = log(LN3/LN2)/log(mpf(3)/2)
print(f"    gap closes at sigma* = ln(log_2 3)/ln(3/2) = {nstr(sig_star,12)}")
print(f"    [Prox note: ln(pi) = {nstr(log(pi),10)} sits "
      f"{nstr(100*abs(sig_star-log(pi))/log(pi),3)}% away -- index-only, not evidence]")

# ---------- B: coset layer bookkeeping ----------
print("\n" + "="*74)
print("  B. COSET LAYER: PSL(2,11) x PSL(2,13) representation accounting")
print("="*74)
def psl_data(q, dims, names):
    ssq = sum(d*d for d in dims)
    print(f"  PSL(2,{q}): |G| = {q*(q*q-1)//2}, irreps {dims}")
    print(f"     sum of squares = {ssq}  ({'OK' if ssq==q*(q*q-1)//2 else 'FAIL'})"
          f"   classes = {len(dims)}")
    for d,n in zip(dims,names): print(f"       dim {d:>3} : {n}")

psl_data(13, [1,13,14,14,7,7,12,12,12],
    ["trivial (chi_0, with Steinberg: the 1+St block = C[P1(F13)])",
     "Steinberg St_13",
     "principal series, cubic torus char pair  <-> our chi_4 channel",
     "principal series, sextic torus char pair <-> our chi_2 channel",
     "split principal (quadratic torus char)   <-> our chi_6 channel",
     "split principal (quadratic, 2nd half)",
     "discrete series (nonsplit torus) -- GAUGE-INVISIBLE",
     "discrete series -- GAUGE-INVISIBLE",
     "discrete series -- GAUGE-INVISIBLE"])
print()
psl_data(11, [1,11,12,12,10,10,5,5],
    ["trivial (with Steinberg: 1+St = C[P1(F11)])",
     "Steinberg St_11",
     "principal series, quintic torus char pair <-> SU(11) channel",
     "principal series, quintic torus char pair <-> SU(11) channel",
     "discrete series -- GAUGE-INVISIBLE",
     "discrete series -- GAUGE-INVISIBLE",
     "biplane rep (discrete, quadratic ODD char) -- GAUGE-INVISIBLE",
     "biplane rep -- GAUGE-INVISIBLE"])

print("\n  Torus-character <-> gauge-channel bijection:")
print("    split torus of PSL(2,p) = (Z/p)*/{+-1}; its characters = EVEN Dirichlet chars mod p")
print("    p=13: torus order 6  = 6 even chars  = our six SU(13) channels   OK")
print("    p=11: torus order 5  = 5 even chars  = SU(11): principal + two quintic pairs   OK")
print("    (11 = 3 mod 4: quadratic char is ODD -> lives in the DISCRETE series -> no gauge channel)")

print("\n  The 168-point stage decomposes:")
print("    C[P1(F11) x P1(F13)] = (1 + St11) x (1 + St13)")
print("    168 = 1 + 11 + 13 + 143      143 = dim(St11 x St13)")
print("    -> the MTFT level 143 is the Steinberg (x) Steinberg component of the coset layer.")
