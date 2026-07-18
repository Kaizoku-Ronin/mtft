import { useState, useEffect, useRef, useCallback } from "react";

// ═══════════════════════════════════════════════════════════════
// MTFT HYPERBOLIC TILING — Poincaré Disk Model
// The modular group SL(2,ℤ) tessellates the hyperbolic plane.
// Each tile is a copy of the fundamental domain F.
// Colored by the three Galois orbits of newforms on X₀(143).
// ═══════════════════════════════════════════════════════════════

const SIZE = 520;
const HALF = SIZE / 2;
const R_DISK = HALF - 20;

// Three newform orbit colors (MTFT palette)
const C_ELECTRON = [0, 235, 200];   // f₁ — moonshine cyan
const C_MUON     = [255, 200, 50];  // f₂ — amber gold
const C_TAU      = [220, 70, 120];  // f₃ — rose
const C_BG       = [8, 8, 14];      // void black
const C_DISK_BG  = [12, 12, 22];
const C_BORDER   = [0, 245, 212];

const lerp = (a, b, t) => a.map((v, i) => Math.round(v + (b[i] - v) * t));

function inverseCayley(u, v) {
  const dn = (1 - u) * (1 - u) + v * v;
  if (dn < 1e-12) return [0, 1000];
  return [-2 * v / dn, Math.max((1 - u * u - v * v) / dn, 1e-10)];
}

function cayley(x, y) {
  const dn = x * x + (y + 1) * (y + 1);
  return [(x * x + y * y - 1) / dn, -2 * x / dn];
}

function modularReduce(x, y) {
  let steps = 0, sCount = 0, tCount = 0;
  let c_acc = 0; // accumulate c component for Γ₀(143) index
  let ma = 1, mb = 0, mc = 0, md = 1;

  for (let i = 0; i < 200; i++) {
    const n = Math.round(x);
    if (n !== 0) {
      x -= n;
      ma -= n * mc; mb -= n * md;
      tCount += Math.abs(n);
      steps += Math.abs(n);
    }
    const norm = x * x + y * y;
    if (norm < 0.9999) {
      const nx = -x / norm, ny = y / norm;
      x = nx; y = ny;
      const ta = mc, tb = md;
      mc = -ma; md = -mb;
      ma = ta; mb = tb;
      sCount++;
      steps++;
    } else break;
  }
  const cMod143 = ((mc % 143) + 143) % 143;
  return { x, y, steps, sCount, tCount, cMod143, ma, mb, mc, md };
}

function farey(order) {
  const fracs = new Set();
  for (let q = 1; q <= order; q++)
    for (let p = 0; p <= q; p++)
      if (gcd(p, q) === 1) fracs.add(`${p}/${q}`);
  return [...fracs].map(s => {
    const [p, q] = s.split("/").map(Number);
    return { p, q };
  });
}

function gcd(a, b) { while (b) { [a, b] = [b, a % b]; } return a; }

function fordCirclePoints(p, q, nPts = 48) {
  const cx = p / q, cy = 1 / (2 * q * q), r = cy;
  const pts = [];
  for (let i = 0; i <= nPts; i++) {
    const th = Math.PI * i / nPts; // upper semicircle only
    const hx = cx + r * Math.cos(th);
    const hy = cy + r * Math.sin(th);
    if (hy > 0.001) {
      const [du, dv] = cayley(hx, hy);
      if (du * du + dv * dv < 0.999) pts.push([du, dv]);
    }
  }
  return pts;
}

export default function HyperbolicTiling() {
  const canvasRef = useRef(null);
  const [colorMode, setColorMode] = useState("orbits");
  const [showFord, setShowFord] = useState(false);
  const [showLabels, setShowLabels] = useState(true);
  const [rendering, setRendering] = useState(true);
  const [hoverInfo, setHoverInfo] = useState(null);
  const tileDataRef = useRef(null);

  const renderTiling = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const imgData = ctx.createImageData(SIZE, SIZE);
    const data = imgData.data;

    const tileInfo = new Int16Array(SIZE * SIZE * 3); // steps, sCount, cMod143

    for (let py = 0; py < SIZE; py++) {
      for (let px = 0; px < SIZE; px++) {
        const idx = (py * SIZE + px) * 4;
        const u = (px - HALF) / R_DISK;
        const v = (py - HALF) / R_DISK;
        const r2 = u * u + v * v;

        if (r2 >= 1.0) {
          data[idx] = C_BG[0]; data[idx+1] = C_BG[1]; data[idx+2] = C_BG[2]; data[idx+3] = 255;
          continue;
        }

        const [hx, hy] = inverseCayley(u, v);
        const res = modularReduce(hx, hy);
        const tIdx = (py * SIZE + px) * 3;
        tileInfo[tIdx] = res.steps;
        tileInfo[tIdx+1] = res.sCount;
        tileInfo[tIdx+2] = res.cMod143;

        let cr, cg, cb;
        const edgeDist = Math.min(
          Math.abs(Math.abs(res.x) - 0.5),
          Math.abs(res.x * res.x + res.y * res.y - 1.0)
        );
        const edgeFade = Math.min(1, edgeDist * 60);
        const depthFade = Math.min(1, res.y * 3);

        if (colorMode === "orbits") {
          const orbit = res.steps % 3;
          const base = orbit === 0 ? C_ELECTRON : orbit === 1 ? C_MUON : C_TAU;
          const bright = 0.3 + 0.7 * depthFade * edgeFade;
          cr = Math.round(base[0] * bright);
          cg = Math.round(base[1] * bright);
          cb = Math.round(base[2] * bright);
        } else if (colorMode === "checker") {
          const v2 = res.steps % 2 === 0 ? 0.85 : 0.2;
          const bright = v2 * edgeFade * depthFade;
          cr = Math.round(C_ELECTRON[0] * bright * 0.4);
          cg = Math.round(C_ELECTRON[1] * bright);
          cb = Math.round(C_ELECTRON[2] * bright * 0.9);
        } else if (colorMode === "depth") {
          const t = Math.min(1, (res.y - 0.866) / 3.0);
          const bright = edgeFade * depthFade;
          const c = lerp(C_TAU, C_ELECTRON, t);
          cr = Math.round(c[0] * bright);
          cg = Math.round(c[1] * bright);
          cb = Math.round(c[2] * bright);
        } else if (colorMode === "gamma143") {
          const hue = (res.cMod143 / 143) * 360;
          const [rr, gg, bb] = hslToRgb(hue, 0.8, 0.35 + 0.3 * edgeFade * depthFade);
          cr = rr; cg = gg; cb = bb;
        }

        // Darken near edge of Poincaré disk
        const rimDark = Math.pow(1 - r2, 0.15);
        data[idx]   = Math.min(255, Math.round(cr * rimDark));
        data[idx+1] = Math.min(255, Math.round(cg * rimDark));
        data[idx+2] = Math.min(255, Math.round(cb * rimDark));
        data[idx+3] = 255;
      }
    }

    ctx.putImageData(imgData, 0, 0);

    // Disk border
    ctx.strokeStyle = `rgba(${C_BORDER.join(",")}, 0.6)`;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(HALF, HALF, R_DISK, 0, 2 * Math.PI);
    ctx.stroke();

    // Ford circles
    if (showFord) {
      const fracs = farey(7);
      ctx.strokeStyle = "rgba(255,255,255,0.35)";
      ctx.lineWidth = 0.8;
      for (const { p, q } of fracs) {
        if (q > 7) continue;
        const pts = fordCirclePoints(p, q);
        if (pts.length < 3) continue;
        ctx.beginPath();
        ctx.moveTo(HALF + pts[0][0] * R_DISK, HALF + pts[0][1] * R_DISK);
        for (let i = 1; i < pts.length; i++)
          ctx.lineTo(HALF + pts[i][0] * R_DISK, HALF + pts[i][1] * R_DISK);
        ctx.stroke();
      }
    }

    // Labels
    if (showLabels) {
      ctx.font = "bold 9px 'JetBrains Mono', 'Courier New', monospace";
      const labels = [
        { tau: [0, 1], label: "τ = i", desc: "center" },
        { tau: [0, 1.732], label: "ρ = e^{iπ/3}", desc: "elliptic" },
        { tau: [0, 10], label: "cusp ∞", desc: "" },
      ];
      for (const lbl of labels) {
        const [du, dv] = cayley(lbl.tau[0], lbl.tau[1]);
        const sx = HALF + du * R_DISK, sy = HALF + dv * R_DISK;
        if (sx > 30 && sx < SIZE - 30 && sy > 30 && sy < SIZE - 30) {
          ctx.fillStyle = "rgba(255,255,255,0.9)";
          ctx.beginPath();
          ctx.arc(sx, sy, 3, 0, 2 * Math.PI);
          ctx.fill();
          ctx.fillStyle = "rgba(255,255,220,0.85)";
          ctx.fillText(lbl.label, sx + 6, sy - 3);
        }
      }
      // Fundamental domain outline
      ctx.strokeStyle = "rgba(255,255,255,0.25)";
      ctx.lineWidth = 1;
      ctx.setLineDash([3, 3]);
      const fdPts = [];
      for (let t = 0; t <= 1; t += 0.02) {
        const angle = Math.PI / 3 + t * Math.PI / 3;
        const [du, dv] = cayley(Math.cos(angle), Math.sin(angle));
        fdPts.push([HALF + du * R_DISK, HALF + dv * R_DISK]);
      }
      for (let t = 0; t <= 1; t += 0.02) {
        const hy = 1 + t * 8;
        const [du, dv] = cayley(-0.5, hy);
        fdPts.push([HALF + du * R_DISK, HALF + dv * R_DISK]);
      }
      if (fdPts.length > 2) {
        ctx.beginPath();
        ctx.moveTo(fdPts[0][0], fdPts[0][1]);
        for (const [fx, fy] of fdPts) ctx.lineTo(fx, fy);
        ctx.stroke();
      }
      ctx.setLineDash([]);
    }

    tileDataRef.current = tileInfo;
    setRendering(false);
  }, [colorMode, showFord, showLabels]);

  useEffect(() => {
    setRendering(true);
    const timer = setTimeout(renderTiling, 50);
    return () => clearTimeout(timer);
  }, [renderTiling]);

  const handleMouseMove = useCallback((e) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const px = Math.round((e.clientX - rect.left) * (SIZE / rect.width));
    const py = Math.round((e.clientY - rect.top) * (SIZE / rect.height));
    const u = (px - HALF) / R_DISK, v = (py - HALF) / R_DISK;
    if (u * u + v * v >= 0.98) { setHoverInfo(null); return; }
    const [hx, hy] = inverseCayley(u, v);
    const res = modularReduce(hx, hy);
    setHoverInfo({
      tau: `${hx.toFixed(3)} + ${hy.toFixed(3)}i`,
      reduced: `${res.x.toFixed(3)} + ${res.y.toFixed(3)}i`,
      steps: res.steps,
      orbit: res.steps % 3,
      coset: res.cMod143,
    });
  }, []);

  const modeLabel = { orbits: "Three Galois Orbits", checker: "Modular Checkerboard",
                      depth: "Modular Depth", gamma143: "Γ₀(143) Cosets" };

  return (
    <div style={{
      background: `rgb(${C_BG.join(",")})`,
      minHeight: "100vh", display: "flex", flexDirection: "column",
      alignItems: "center", fontFamily: "'JetBrains Mono', 'Courier New', monospace",
      color: "#d0d0d0", padding: "16px 8px",
    }}>
      {/* Title */}
      <div style={{ textAlign: "center", marginBottom: 12 }}>
        <h1 style={{
          fontSize: 18, fontWeight: 700, letterSpacing: 3,
          color: `rgb(${C_BORDER.join(",")})`, margin: 0,
          textShadow: `0 0 20px rgba(${C_BORDER.join(",")},0.4)`,
        }}>
          SL(2,ℤ) HYPERBOLIC TILING
        </h1>
        <div style={{ fontSize: 10, color: "#666", marginTop: 4, letterSpacing: 1 }}>
          POINCARÉ DISK MODEL — X₀(143) STRUCTURE — {modeLabel[colorMode]}
        </div>
      </div>

      <div style={{ display: "flex", gap: 16, alignItems: "flex-start", flexWrap: "wrap", justifyContent: "center" }}>
        {/* Canvas */}
        <div style={{ position: "relative" }}>
          <canvas
            ref={canvasRef}
            width={SIZE}
            height={SIZE}
            onMouseMove={handleMouseMove}
            onMouseLeave={() => setHoverInfo(null)}
            style={{
              borderRadius: "50%",
              boxShadow: `0 0 40px rgba(${C_BORDER.join(",")},0.15), 0 0 80px rgba(${C_BORDER.join(",")},0.05)`,
              cursor: "crosshair",
              maxWidth: "min(520px, 85vw)", height: "auto",
            }}
          />
          {rendering && (
            <div style={{
              position: "absolute", top: "50%", left: "50%", transform: "translate(-50%,-50%)",
              color: `rgb(${C_BORDER.join(",")})`, fontSize: 12, letterSpacing: 2,
            }}>
              TESSELLATING...
            </div>
          )}
        </div>

        {/* Controls */}
        <div style={{
          background: "rgba(255,255,255,0.03)", border: "1px solid rgba(255,255,255,0.08)",
          borderRadius: 8, padding: 16, width: 200, fontSize: 11,
        }}>
          <div style={{ color: `rgb(${C_BORDER.join(",")})`, fontSize: 9, letterSpacing: 2,
                        textTransform: "uppercase", marginBottom: 10 }}>
            Color Mode
          </div>
          {Object.entries(modeLabel).map(([k, v]) => (
            <label key={k} style={{ display: "flex", alignItems: "center", gap: 6,
                                    marginBottom: 5, cursor: "pointer" }}>
              <input type="radio" name="cm" value={k} checked={colorMode === k}
                     onChange={() => setColorMode(k)}
                     style={{ accentColor: `rgb(${C_BORDER.join(",")})` }} />
              <span style={{ color: colorMode === k ? `rgb(${C_BORDER.join(",")})` : "#888" }}>
                {v}
              </span>
            </label>
          ))}

          <div style={{ borderTop: "1px solid rgba(255,255,255,0.06)", margin: "12px 0" }} />

          <label style={{ display: "flex", alignItems: "center", gap: 6, cursor: "pointer", marginBottom: 5 }}>
            <input type="checkbox" checked={showFord}
                   onChange={() => setShowFord(!showFord)}
                   style={{ accentColor: `rgb(${C_BORDER.join(",")})` }} />
            <span style={{ color: "#aaa" }}>Ford Circles (F₇)</span>
          </label>
          <label style={{ display: "flex", alignItems: "center", gap: 6, cursor: "pointer" }}>
            <input type="checkbox" checked={showLabels}
                   onChange={() => setShowLabels(!showLabels)}
                   style={{ accentColor: `rgb(${C_BORDER.join(",")})` }} />
            <span style={{ color: "#aaa" }}>Labels & Domain</span>
          </label>

          <div style={{ borderTop: "1px solid rgba(255,255,255,0.06)", margin: "12px 0" }} />

          {/* Legend */}
          <div style={{ color: `rgb(${C_BORDER.join(",")})`, fontSize: 9, letterSpacing: 2,
                        textTransform: "uppercase", marginBottom: 8 }}>
            {colorMode === "orbits" ? "Galois Orbits" : colorMode === "gamma143" ? "168 Cosets" : "Legend"}
          </div>
          {colorMode === "orbits" && [
            { c: C_ELECTRON, l: "f₁ — Electron (dim 1)" },
            { c: C_MUON, l: "f₂ — Muon (dim 4)" },
            { c: C_TAU, l: "f₃ — Tau (dim 6)" },
          ].map(({ c, l }) => (
            <div key={l} style={{ display: "flex", alignItems: "center", gap: 6, marginBottom: 4 }}>
              <div style={{ width: 10, height: 10, borderRadius: 2,
                            background: `rgb(${c.join(",")})` }} />
              <span style={{ color: "#999", fontSize: 10 }}>{l}</span>
            </div>
          ))}
          {colorMode === "gamma143" && (
            <div style={{ color: "#999", fontSize: 10, lineHeight: 1.5 }}>
              168 cosets of Γ₀(143) in SL(2,ℤ)<br/>
              Index = N∏(1+1/p) = 168<br/>
              = |PSL(2,7)| = dim SU(13)
            </div>
          )}

          {/* Hover info */}
          {hoverInfo && (
            <>
              <div style={{ borderTop: "1px solid rgba(255,255,255,0.06)", margin: "12px 0" }} />
              <div style={{ fontSize: 9, color: `rgb(${C_BORDER.join(",")})`, letterSpacing: 2,
                            textTransform: "uppercase", marginBottom: 6 }}>
                Cursor
              </div>
              <div style={{ color: "#ccc", fontSize: 10, lineHeight: 1.7 }}>
                τ = {hoverInfo.tau}<br/>
                reduced: {hoverInfo.reduced}<br/>
                steps: {hoverInfo.steps}<br/>
                orbit: f{hoverInfo.orbit + 1}<br/>
                Γ₀(143) coset: {hoverInfo.coset}
              </div>
            </>
          )}
        </div>
      </div>

      {/* Footer */}
      <div style={{ marginTop: 16, textAlign: "center", fontSize: 9, color: "#444", letterSpacing: 1 }}>
        MTFT — Modular Time Field Theory — The fundamental domain F of SL(2,ℤ) tessellates the
        hyperbolic plane ℍ into infinitely many copies, rendered here in the Poincaré disk model.
      </div>
    </div>
  );
}

function hslToRgb(h, s, l) {
  h /= 360;
  const a = s * Math.min(l, 1 - l);
  const f = (n) => {
    const k = (n + h * 12) % 12;
    return l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1));
  };
  return [Math.round(f(0) * 255), Math.round(f(8) * 255), Math.round(f(4) * 255)];
}
