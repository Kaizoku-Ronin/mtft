# --- stage 5 (Integration Plan v0.1): re-pointed to the mtft package. ---
# internal() delegates to mtft.chain.internal (certified spectrally
# identical: g to 9e-15, full H(u) spectrum to 2.5e-14; B differs only
# by gauge inside near-degenerate blocks).  Where a local gsq() remains
# it is retained VERBATIM: it evaluates complex u in the f64 backend,
# which mtft.ep deliberately does not offer (PR-20 floors) -- S5-1.
# ------------------------------------------------------------------------
import math, numpy as np, importlib.util, os
spec = importlib.util.spec_from_file_location("p14", os.path.join(os.path.dirname(os.path.abspath(__file__)), "pr14_ep.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
out = open("ks3.out", "w", buffering=1)
out.write("kappa   R=|u_c|      s*        R/s*     arg     nroots<6\n")
for kap in (3.0, 4.0, 5.0, 6.0, 8.0):
    g, Bm = m.internal(N=700, kappa=kap, nb=30)
    roots = []
    for th in np.linspace(1.4, 2.9, 16):          # local net around the
        for r in np.linspace(1.5, 5.0, 15):       # known EP sector
            u = m.newton_ep(r*np.exp(1j*th), g, Bm, iters=50)
            if (np.isfinite(u) and abs(u) < 6.0
                    and abs(m.gsq(u, g, Bm)) < 1e-11
                    and all(abs(u-v) > 1e-5 for v in roots)):
                roots.append(u)
    roots.sort(key=abs)
    ss = g[1] / (Bm[0, 0] + Bm[1, 1])
    if roots:
        u = roots[0]
        out.write(f"{kap:5.1f}  {abs(u):9.5f}  {ss:8.5f}  {abs(u)/ss:8.3f}"
                  f"  {math.atan2(u.imag,u.real):7.4f}  {len(roots)}\n")
    else:
        out.write(f"{kap:5.1f}  none found in |u|<6\n")
out.write("done\n"); out.close()
