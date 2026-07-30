# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
import math, numpy as np, importlib.util, os
spec = importlib.util.spec_from_file_location(
    "c", os.path.join(os.path.dirname(os.path.abspath(__file__)), "ep_census.py"))
c = importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
g, B = c.internal()

# S5-1: retained verbatim -- complex-u f64 gsq predates mtft.ep's
# real-axis diabatic-centre design; frozen so this study's record stands.
def gsq(u):
    ev = np.linalg.eigvals(np.diag(g) - u * B)
    ev = ev[np.argsort(ev.real)]
    d = ev[1:] - ev[:-1]
    return d[int(np.argmin(np.abs(d)))] ** 2

def newton(u0, iters=80, h=1e-7):
    u = complex(u0)
    for _ in range(iters):
        f = gsq(u); fp = (gsq(u+h) - gsq(u-h)) / (2*h)
        if fp == 0: break
        s = f/fp; u -= s
        if abs(s) < 1e-14: break
    return u

roots = []
for th in np.linspace(0, 2*math.pi, 90, endpoint=False):
    for r in np.linspace(0.15, 1.45, 40):
        u = newton(r*np.exp(1j*th))
        if abs(u) < 1.5 and abs(gsq(u)) < 1e-12 and all(abs(u-v) > 1e-6 for v in roots):
            roots.append(u)
roots.sort(key=abs)
print("EPs found inside |u|<1.5:", len(roots), "(census says 2)", flush=True)
for u in roots[:4]:
    # continue the colliding pair back to u=0
    ev = np.linalg.eigvals(np.diag(g) - u*B); ev = ev[np.argsort(ev.real)]
    d = ev[1:]-ev[:-1]; i = int(np.argmin(np.abs(d)))
    pair = np.array([ev[i], ev[i+1]])
    for t in np.linspace(1.0, 0.0, 500)[1:]:
        cur = np.linalg.eigvals(np.diag(g) - u*t*B)
        new, used = [], set()
        for p in pair:
            j = int(np.argmin([abs(p-cc) if k not in used else 9e9 for k, cc in enumerate(cur)]))
            used.add(j); new.append(cur[j])
        pair = np.array(new)
    idx = sorted(int(np.argmin(np.abs(g - p.real))) for p in pair)
    print(f"  u={u.real:+.5f}{u.imag:+.5f}i  |u|={abs(u):.5f}  levels {idx}", flush=True)
