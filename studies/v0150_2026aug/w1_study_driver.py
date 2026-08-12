"""W1 reproduction driver (module-based; supersedes the raw in-session script,
whose synthetic section carried correction W1-c2 -- see the study report).
Reproduces every number in report Sec. 3. Zero sweep: loads
zeros_gamma_T100.npy if present, else computes via mpmath (~70 s) and saves.
Run:  py w1_study_driver.py   (Windows) /  python3 ...  (elsewhere)"""
import os
import numpy as np
import mpmath as mp
import weil

T, LAM, ETA = 100.0, 1.0, 0.2
g = weil.gabor(T, LAM)
print(f"T={T}  l={g['l']:.4f}  X={g['X']:.4f}  d={g['d']}  h={g['h']:.4f}")

win = weil.Window(g["L"], ETA)
print(f"taper: a={win.a:.6f} (law {1-0.603*ETA:.4f})  b={win.b:.6f} (law {1-0.688*ETA:.4f})")

# Lemma 2.2 frame identity
kk = np.arange(-4000, 4000 + g["d"]); tk = T + g["h"] * kk
rng = np.random.default_rng(11); err = 0.0
for _ in range(5):
    ta, tb = rng.uniform(80, 220, 2)
    err = max(err, abs(np.sum(win.phihat(ta - tk) * win.phihat(tb - tk))
                       - g["L"] * win.Phi(ta - tb)) / (g["L"] * abs(win.Phi(ta - tb))))
print(f"Lemma 2.2 max rel err = {err:.2e}")

# zeros
cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "zeros_gamma_T100.npy")
if os.path.exists(cache):
    gam = np.load(cache)
else:
    mp.mp.dps = 15
    gam, n = [], 1
    while True:
        z = float(mp.zetazero(n).imag); gam.append(z)
        if z > 640: break
        n += 1
        if n % 200 == 0: np.save(cache, np.array(gam))
    gam = np.array(gam); np.save(cache, gam)
print(f"{gam.size} zeros, gamma_max={gam[-1]:.2f}")

# the two routes and the E2 identity
Gp, g, win = weil.G_prime(T, LAM, ETA, pad=250.0, dtau=0.005, win=win)
Gz, _ = weil.G_zero(gam, T, LAM, ETA, win=win)
rel = np.max(np.abs(Gp - Gz)) / np.max(np.abs(Gp))
print(f"E2 prime-vs-zero: max rel discrepancy = {rel:.3e}")
print(f"  tr {np.trace(Gp):.6f} / {np.trace(Gz):.6f}   trG^2 {np.sum(Gp**2):.4f} / {np.sum(Gz**2):.4f}")

in_I = (gam > T) & (gam <= 2 * T)
in_Ip = (gam > T - np.sqrt(T)) & (gam <= 2 * T + np.sqrt(T))
N_I, N_Ip = int(in_I.sum()), int(in_Ip.sum())
Uc = win.phihat(gam[in_Ip][:, None] - g["tau_k"][None, :], direct=True)
cert = weil.certificates(Gz, g, win, N_I, N_Ip, A_Ip=Uc.T @ Uc)
print(f"N(I)={N_I} (RvM {T/(2*np.pi)*g['ell1']:.2f})  N(I')={N_Ip}")
print(f"C/N = {cert['C_over_N']:.4f}  law = {cert['law']:.4f}  F(1)=0.75")
print(f"spectrum: min {cert['eigs'].min():.4f}  n_+ = {(cert['eigs']>0).sum()} of {g['d']}")
print(f"rank-trace cert = {cert['cert_rank_trace']:+.2f} ({cert['cert_rank_trace']/N_I:+.3f} N);  "
      f"CS cert = {cert['cert_cs']:+.2f} ({cert['cert_cs']/N_I:+.3f} N);  "
      f"||A||^2/tr = {cert['frob_over_tr']:.4f}")

# synthetic inertia, non-vacuous (p < d) configurations
ords = gam[in_Ip]
for label, o, deps, m in [
    ("20 pairs @0.30", ords[:20], np.full(20, 0.30), None),
    ("10 pairs @0.25 + 8 doubled on-line", ords[:18],
     np.concatenate([np.full(10, 0.25), np.zeros(8)]),
     np.concatenate([np.ones(10), 2 * np.ones(8)])),
]:
    Gs, cnt = weil.G_zero(o, T, LAM, ETA, depths=deps, mults=m, mirror=False, win=win)
    ev = np.linalg.eigvalsh(Gs / g["L"])
    npos, nneg = int((ev > 1e-8).sum()), int((ev < -1e-8).sum())
    Ah = Gs / (win.a * g["L"] ** 2)
    cs = 4 * np.trace(Ah) - 2 * cnt["N"] - np.sum(Ah * Ah)
    bound = cnt["s1"] + cnt["s2"] + cnt["p"]
    print(f"synthetic [{label}]: n+={npos}<={bound} {npos<=bound}; n-={nneg}<=p {nneg<=cnt['p']}; "
          f"cert={cs:+.1f}<=s1 {cs<=cnt['s1']+1e-9}")

# lemma audit, MT constant, CC-02
a = weil.audit_rank_trace(100_000)
print(f"Lemma 3.2: {a['trials']} instances, {a['violations']} violations, min gap {a['min_gap']:.3e}")
print(f"MT: c*_1 = {weil.mt_constant(1.0):.12f};  2-1/c* = {2-1/weil.mt_constant(1.0):.12f}")
r = weil.w_series_check(3, 300_000)
print(f"CC-02 @ s=3: |sum - (-zeta zeta')| = {r['diff_correct']:.2e}; "
      f"wrong forms miss by {r['diff_paper1']:.2e} / {r['diff_ag']:.2e}")
