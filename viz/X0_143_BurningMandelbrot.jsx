import { useState, useEffect, useRef, useCallback } from "react";

// ════════════════════════════════════════════════════════════════
// MTFT X₀(143) BURNING MANDELBROT
// The vacuum fractal modulated by the arithmetic of level 143
// ════════════════════════════════════════════════════════════════

// Hecke traces Tr(T_n) for newspace 143.2.a from LMFDB (n=1..200)
const HECKE_TRACES = [
  11, 3, 2, 9, 0, -4, 8, 3, 13, -2, -3, 0, 1, -16, -6, 17, 2, -9,
  0, -18, -12, 3, 14, -36, 31, -3, 2, 12, -10, -48, -10, -21, -2,
  -34, -28, -13, 16, -14, 4, -14, 14, 6, 20, -3, 2, 0, -16, -6, 23,
  25, -12, 7, -2, 36, 0, -6, 52, 18, -6, 0, 2, 28, 0, 33, 2, 2, 10,
  2, -14, 4, -38, -1, 38, 50, -28, -36, 4, -2, 12, 26, 3, 18, -28,
  44, 40, -8, -20, 15, -20, 14, 0, 86, 2, 32, -52, -64, 36, -17, -9,
  27, 18, 24, -16, -15, 8, -26, -24, -14, -2, -14, -34, 12, 32, 56,
  -50, -46, 13, -4, 12, 20, 11, 6, -64, 4, -42, -82, 0, -25, -20,
  14, -28, 8, 0, 36, -102, -86, -4, -52, -4, -80, -52, 40, -9, -23,
  -12, 28, -6, 114, 6, 72, -12, -48, 18, 8, -18, -12, 4, -48, -20,
  -98, -32, 115, 24, 106, -14, -68, -20, -14, 11, 28, 8, 24, -2, 72,
  12, -7, 50, -54, 6, 22, -24, -8, -4, -16, 10, -4, 10, -4, -64, -12,
  -6, 138, 86, 14, 12, -39, -2, -21, -8, 173
];

// MTFT Constants
const FEIGENBAUM_DELTA = 4.669201609;
const FEIGENBAUM_ALPHA = 2.502907875;
const FEIGENBAUM_C = -1.401155189; // Feigenbaum point on real axis
const T_INF = 0.4687741272; // Vacuum torque flux -ζ'(2)/2

// Monster group data
const LN_MONSTER = 124.1264; // ln(|M|)
const ALPHA_INV_MONSTER = LN_MONSTER + 13 - 1/11; // ≈ 137.0355
const ALPHA_INV_FEIGENBAUM = 2 * Math.PI * FEIGENBAUM_DELTA**2 + 1/(4*FEIGENBAUM_DELTA); // leading terms ≈ 137.1

// Heegner number
const HEEGNER_163 = Math.exp(Math.PI * Math.sqrt(163));

// Key Mandelbrot anatomy points
const ANATOMY = {
  cardioid_cusp: { re: 0.25, im: 0, label: "U(1) Cusp" },
  period2: { re: -1.0, im: 0, label: "SU(2)" },
  period3_top: { re: -0.1226, im: 0.7449, label: "SU(3)" },
  feigenbaum: { re: FEIGENBAUM_C, im: 0, label: "δ = 4.669..." },
  tip: { re: -2.0, im: 0, label: "Tip c=-2" },
};

// Color palettes inspired by the two-vacuum architecture
const PALETTES = {
  bosonic: (t) => {
    // Mandelbrot sector: deep blue → cyan → white
    const r = Math.floor(255 * Math.pow(t, 3));
    const g = Math.floor(255 * Math.pow(t, 1.5));
    const b = Math.floor(255 * Math.pow(t, 0.7));
    return [r, g, b];
  },
  fermionic: (t) => {
    // Burning Ship sector: deep red → orange → gold
    const r = Math.floor(255 * Math.pow(t, 0.6));
    const g = Math.floor(255 * Math.pow(t, 1.8));
    const b = Math.floor(255 * Math.pow(t, 4));
    return [r, g, b];
  },
  modular: (t, fold_ratio) => {
    // Interpolate between bosonic and fermionic based on fold ratio
    const f = Math.max(0, Math.min(1, fold_ratio));
    const bos = PALETTES.bosonic(t);
    const fer = PALETTES.fermionic(t);
    return [
      Math.floor(bos[0] * (1 - f) + fer[0] * f),
      Math.floor(bos[1] * (1 - f) + fer[1] * f),
      Math.floor(bos[2] * (1 - f) + fer[2] * f),
    ];
  },
  // Three-sector coloring: R=dim6(tau), G=dim4(muon), B=dim1(electron)
  threeOrbit: (t, sector_weights) => {
    const [w1, w4, w6] = sector_weights;
    const norm = w1 + w4 + w6 || 1;
    const base = Math.pow(t, 0.8) * 255;
    return [
      Math.floor(base * (0.3 + 0.7 * w6 / norm)), // R: tau (dim-6)
      Math.floor(base * (0.3 + 0.7 * w4 / norm)), // G: muon (dim-4)
      Math.floor(base * (0.3 + 0.7 * w1 / norm)), // B: electron (dim-1)
    ];
  }
};

// ════════════════════════════════════════════════════════════════
// ITERATION ENGINES
// ════════════════════════════════════════════════════════════════

function iterateMandelbrot(cr, ci, maxIter) {
  let zr = 0, zi = 0;
  for (let n = 0; n < maxIter; n++) {
    const zr2 = zr * zr, zi2 = zi * zi;
    if (zr2 + zi2 > 4) return n / maxIter;
    zi = 2 * zr * zi + ci;
    zr = zr2 - zi2 + cr;
  }
  return -1; // inside set
}

function iterateBurningShip(cr, ci, maxIter) {
  let zr = 0, zi = 0;
  for (let n = 0; n < maxIter; n++) {
    const zr2 = zr * zr, zi2 = zi * zi;
    if (zr2 + zi2 > 4) return n / maxIter;
    zi = 2 * Math.abs(zr) * Math.abs(zi) + ci;
    zr = zr2 - zi2 + cr;
  }
  return -1;
}

// THE X₀(143) BURNING MANDELBROT
// At each iteration step n, the Hecke trace a(n) determines
// whether we fold (Burning Ship) or stay analytic (Mandelbrot).
// Negative traces → fermionic fold; positive → bosonic analytic.
function iterateX0_143(cr, ci, maxIter, mode = "hecke") {
  let zr = 0, zi = 0;
  let foldCount = 0;
  let sectorWeights = [0, 0, 0]; // [dim1, dim4, dim6]

  for (let n = 0; n < maxIter; n++) {
    const zr2 = zr * zr, zi2 = zi * zi;
    if (zr2 + zi2 > 4) {
      return {
        t: n / maxIter,
        foldRatio: foldCount / (n || 1),
        sectorWeights
      };
    }

    const traceIdx = n % HECKE_TRACES.length;
    const a_n = HECKE_TRACES[traceIdx];

    if (mode === "hecke") {
      // Hecke-modulated: negative trace → fold (fermionic)
      if (a_n < 0) {
        zi = 2 * Math.abs(zr) * Math.abs(zi) + ci;
        foldCount++;
        // Weight by which sector contributes most
        const absA = Math.abs(a_n);
        if (absA <= 11) sectorWeights[0] += absA; // dim-1 scale
        else if (absA <= 44) sectorWeights[1] += absA; // dim-4 scale
        else sectorWeights[2] += absA; // dim-6 scale
      } else {
        zi = 2 * zr * zi + ci;
        sectorWeights[0] += Math.max(0, a_n) * 0.1;
      }
    } else if (mode === "smooth") {
      // Smooth interpolation: mixing angle from trace
      const theta = Math.atan2(Math.max(0, -a_n), Math.max(1, a_n)) * 2 / Math.PI;
      const absZr = Math.abs(zr), absZi = Math.abs(zi);
      const analyticZi = 2 * zr * zi + ci;
      const foldedZi = 2 * absZr * absZi + ci;
      zi = analyticZi * (1 - theta) + foldedZi * theta;
      foldCount += theta;
      sectorWeights[0] += (1 - theta);
      sectorWeights[1] += theta * 0.6;
      sectorWeights[2] += theta * 0.4;
    } else if (mode === "168") {
      // Mod-168 iteration: fold every 168th step (ψ(143) = |PSL(2,7)|)
      if (n % 168 < 84 && a_n < 0) {
        zi = 2 * Math.abs(zr) * Math.abs(zi) + ci;
        foldCount++;
      } else {
        zi = 2 * zr * zi + ci;
      }
      sectorWeights[n % 3] += 1;
    }

    zr = zr2 - zi2 + cr;
  }

  return { t: -1, foldRatio: foldCount / maxIter, sectorWeights };
}

// ════════════════════════════════════════════════════════════════
// REACT COMPONENT
// ════════════════════════════════════════════════════════════════

export default function BurningMandelbrotX0143() {
  const canvasRef = useRef(null);
  const [fractalType, setFractalType] = useState("x0_143_hecke");
  const [colorMode, setColorMode] = useState("modular");
  const [maxIter, setMaxIter] = useState(128);
  const [showOverlay, setShowOverlay] = useState(true);
  const [rendering, setRendering] = useState(false);
  const [stats, setStats] = useState(null);

  // Viewport state
  const [view, setView] = useState({
    centerX: -0.5, centerY: 0, zoom: 1.5
  });

  const WIDTH = 720;
  const HEIGHT = 540;

  const render = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const imgData = ctx.createImageData(WIDTH, HEIGHT);
    const data = imgData.data;

    setRendering(true);
    let totalFolds = 0;
    let totalPixels = 0;
    let insideCount = 0;

    const aspect = WIDTH / HEIGHT;
    const xRange = view.zoom * aspect;
    const yRange = view.zoom;

    for (let py = 0; py < HEIGHT; py++) {
      for (let px = 0; px < WIDTH; px++) {
        const cr = view.centerX + (px / WIDTH - 0.5) * 2 * xRange;
        const ci = view.centerY + (py / HEIGHT - 0.5) * 2 * yRange;

        let r = 0, g = 0, b = 0;

        if (fractalType === "mandelbrot") {
          const t = iterateMandelbrot(cr, ci, maxIter);
          if (t < 0) { insideCount++; }
          else {
            const c = PALETTES.bosonic(t);
            r = c[0]; g = c[1]; b = c[2];
          }
        } else if (fractalType === "burning_ship") {
          const t = iterateBurningShip(cr, -ci, maxIter);
          if (t < 0) { insideCount++; }
          else {
            const c = PALETTES.fermionic(t);
            r = c[0]; g = c[1]; b = c[2];
          }
        } else {
          // X₀(143) variants
          const mode = fractalType === "x0_143_smooth" ? "smooth"
                     : fractalType === "x0_143_168" ? "168" : "hecke";
          const result = iterateX0_143(cr, ci, maxIter, mode);

          if (result.t < 0) {
            insideCount++;
            // Color the interior by fold ratio
            const fr = result.foldRatio;
            r = Math.floor(20 * fr);
            g = Math.floor(10 * (1 - fr));
            b = Math.floor(15);
          } else {
            totalFolds += result.foldRatio;
            totalPixels++;

            if (colorMode === "modular") {
              const c = PALETTES.modular(result.t, result.foldRatio);
              r = c[0]; g = c[1]; b = c[2];
            } else if (colorMode === "three_orbit") {
              const c = PALETTES.threeOrbit(result.t, result.sectorWeights);
              r = c[0]; g = c[1]; b = c[2];
            } else {
              // Classic escape time
              const c = PALETTES.bosonic(result.t);
              r = c[0]; g = c[1]; b = c[2];
            }
          }
        }

        const idx = (py * WIDTH + px) * 4;
        data[idx] = Math.min(255, r);
        data[idx + 1] = Math.min(255, g);
        data[idx + 2] = Math.min(255, b);
        data[idx + 3] = 255;
      }
    }

    ctx.putImageData(imgData, 0, 0);

    // Overlay annotations
    if (showOverlay && fractalType !== "burning_ship") {
      ctx.save();
      ctx.font = "bold 10px 'Courier New', monospace";

      Object.entries(ANATOMY).forEach(([key, pt]) => {
        const sx = ((pt.re - view.centerX) / (2 * xRange) + 0.5) * WIDTH;
        const sy = ((pt.im - view.centerY) / (2 * yRange) + 0.5) * HEIGHT;
        if (sx > 10 && sx < WIDTH - 10 && sy > 10 && sy < HEIGHT - 10) {
          ctx.fillStyle = "rgba(255,255,255,0.9)";
          ctx.strokeStyle = "rgba(255,200,50,0.8)";
          ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.arc(sx, sy, 4, 0, Math.PI * 2);
          ctx.fill();
          ctx.stroke();
          ctx.fillStyle = "rgba(255,255,200,0.95)";
          ctx.fillText(pt.label, sx + 7, sy - 4);
        }
      });
      ctx.restore();
    }

    setStats({
      avgFoldRatio: totalPixels ? (totalFolds / totalPixels).toFixed(4) : "N/A",
      insideFraction: ((insideCount / (WIDTH * HEIGHT)) * 100).toFixed(2),
      area: (insideCount / (WIDTH * HEIGHT) * 4 * xRange * yRange).toFixed(6),
    });
    setRendering(false);
  }, [fractalType, colorMode, maxIter, showOverlay, view]);

  useEffect(() => {
    const timer = setTimeout(render, 50);
    return () => clearTimeout(timer);
  }, [render]);

  // Zoom on click
  const handleCanvasClick = (e) => {
    const rect = canvasRef.current.getBoundingClientRect();
    const px = e.clientX - rect.left;
    const py = e.clientY - rect.top;
    const aspect = WIDTH / HEIGHT;
    const xRange = view.zoom * aspect;
    const yRange = view.zoom;
    const cr = view.centerX + (px / WIDTH - 0.5) * 2 * xRange;
    const ci = view.centerY + (py / HEIGHT - 0.5) * 2 * yRange;
    const zoomFactor = e.shiftKey ? 2 : 0.5;
    setView({ centerX: cr, centerY: ci, zoom: view.zoom * zoomFactor });
  };

  const resetView = () => setView({ centerX: -0.5, centerY: 0, zoom: 1.5 });

  return (
    <div style={{
      background: "#0a0a0f",
      color: "#e0ddd5",
      minHeight: "100vh",
      fontFamily: "'JetBrains Mono', 'Fira Code', 'Courier New', monospace",
      padding: "20px",
      boxSizing: "border-box",
    }}>
      {/* Header */}
      <div style={{
        borderBottom: "1px solid rgba(255,200,50,0.3)",
        paddingBottom: "12px",
        marginBottom: "16px",
      }}>
        <h1 style={{
          margin: 0,
          fontSize: "18px",
          fontWeight: 700,
          color: "#ffd866",
          letterSpacing: "2px",
          textTransform: "uppercase",
        }}>
          X₀(143) Burning Mandelbrot
        </h1>
        <p style={{
          margin: "4px 0 0",
          fontSize: "11px",
          color: "#8a8a7a",
          letterSpacing: "1px",
        }}>
          MTFT · Hecke-Modulated Vacuum Fractal · Bosonic ↔ Fermionic Phase Transition
        </p>
      </div>

      {/* Main layout */}
      <div style={{ display: "flex", gap: "16px", flexWrap: "wrap" }}>
        {/* Canvas */}
        <div style={{ flex: "1 1 720px" }}>
          <canvas
            ref={canvasRef}
            width={WIDTH}
            height={HEIGHT}
            onClick={handleCanvasClick}
            style={{
              width: "100%",
              maxWidth: WIDTH,
              border: "1px solid rgba(255,200,50,0.2)",
              cursor: "crosshair",
              imageRendering: "pixelated",
            }}
          />
          <div style={{
            display: "flex",
            justifyContent: "space-between",
            fontSize: "9px",
            color: "#666",
            marginTop: "4px",
            maxWidth: WIDTH,
          }}>
            <span>Click to zoom in · Shift+click to zoom out</span>
            <span>
              center: ({view.centerX.toFixed(6)}, {view.centerY.toFixed(6)}) · zoom: {view.zoom.toFixed(6)}
            </span>
          </div>
        </div>

        {/* Controls panel */}
        <div style={{
          flex: "0 0 260px",
          background: "rgba(255,255,255,0.03)",
          border: "1px solid rgba(255,200,50,0.15)",
          padding: "14px",
          fontSize: "11px",
        }}>
          {/* Fractal Type */}
          <div style={{ marginBottom: "14px" }}>
            <label style={{ color: "#ffd866", fontSize: "10px", letterSpacing: "1px", textTransform: "uppercase", display: "block", marginBottom: "6px" }}>
              Fractal Type
            </label>
            {[
              ["x0_143_hecke", "X₀(143) Hecke"],
              ["x0_143_smooth", "X₀(143) Smooth"],
              ["x0_143_168", "X₀(143) mod-168"],
              ["mandelbrot", "Mandelbrot (bosonic)"],
              ["burning_ship", "Burning Ship (fermionic)"],
            ].map(([val, label]) => (
              <div key={val} style={{ marginBottom: "3px" }}>
                <label style={{ cursor: "pointer", display: "flex", alignItems: "center", gap: "6px" }}>
                  <input
                    type="radio"
                    name="fractalType"
                    value={val}
                    checked={fractalType === val}
                    onChange={() => setFractalType(val)}
                    style={{ accentColor: "#ffd866" }}
                  />
                  <span style={{ color: fractalType === val ? "#ffd866" : "#999" }}>{label}</span>
                </label>
              </div>
            ))}
          </div>

          {/* Color Mode */}
          <div style={{ marginBottom: "14px" }}>
            <label style={{ color: "#ffd866", fontSize: "10px", letterSpacing: "1px", textTransform: "uppercase", display: "block", marginBottom: "6px" }}>
              Color Mode
            </label>
            {[
              ["modular", "Bosonic↔Fermionic"],
              ["three_orbit", "Three Galois Orbits"],
              ["classic", "Classic Escape"],
            ].map(([val, label]) => (
              <div key={val} style={{ marginBottom: "3px" }}>
                <label style={{ cursor: "pointer", display: "flex", alignItems: "center", gap: "6px" }}>
                  <input
                    type="radio"
                    name="colorMode"
                    value={val}
                    checked={colorMode === val}
                    onChange={() => setColorMode(val)}
                    style={{ accentColor: "#ffd866" }}
                  />
                  <span style={{ color: colorMode === val ? "#ffd866" : "#999" }}>{label}</span>
                </label>
              </div>
            ))}
          </div>

          {/* Max Iterations */}
          <div style={{ marginBottom: "14px" }}>
            <label style={{ color: "#ffd866", fontSize: "10px", letterSpacing: "1px", textTransform: "uppercase", display: "block", marginBottom: "6px" }}>
              Max Iterations: {maxIter}
            </label>
            <input
              type="range"
              min={32}
              max={512}
              step={16}
              value={maxIter}
              onChange={(e) => setMaxIter(Number(e.target.value))}
              style={{ width: "100%", accentColor: "#ffd866" }}
            />
          </div>

          {/* Overlay toggle */}
          <div style={{ marginBottom: "14px" }}>
            <label style={{ cursor: "pointer", display: "flex", alignItems: "center", gap: "6px" }}>
              <input
                type="checkbox"
                checked={showOverlay}
                onChange={() => setShowOverlay(!showOverlay)}
                style={{ accentColor: "#ffd866" }}
              />
              <span style={{ color: "#999" }}>Show MTFT anatomy</span>
            </label>
          </div>

          {/* Reset */}
          <button
            onClick={resetView}
            style={{
              background: "rgba(255,200,50,0.1)",
              border: "1px solid rgba(255,200,50,0.3)",
              color: "#ffd866",
              padding: "6px 12px",
              cursor: "pointer",
              fontSize: "10px",
              letterSpacing: "1px",
              width: "100%",
              marginBottom: "14px",
            }}
          >
            RESET VIEW
          </button>

          {/* Stats */}
          {stats && (
            <div style={{
              borderTop: "1px solid rgba(255,200,50,0.15)",
              paddingTop: "10px",
              fontSize: "10px",
              lineHeight: "1.8",
            }}>
              <div style={{ color: "#ffd866", fontSize: "9px", letterSpacing: "1px", marginBottom: "4px", textTransform: "uppercase" }}>
                Diagnostics
              </div>
              <div><span style={{ color: "#888" }}>Avg fold ratio:</span> <span style={{ color: "#ff9944" }}>{stats.avgFoldRatio}</span></div>
              <div><span style={{ color: "#888" }}>Inside fraction:</span> {stats.insideFraction}%</div>
              <div><span style={{ color: "#888" }}>Viewport area:</span> {stats.area}</div>
            </div>
          )}

          {/* MTFT Constants */}
          <div style={{
            borderTop: "1px solid rgba(255,200,50,0.15)",
            paddingTop: "10px",
            marginTop: "10px",
            fontSize: "9px",
            lineHeight: "1.9",
            color: "#777",
          }}>
            <div style={{ color: "#ffd866", fontSize: "9px", letterSpacing: "1px", marginBottom: "4px", textTransform: "uppercase" }}>
              MTFT Constants
            </div>
            <div>α⁻¹<sub>Monster</sub> = ln|M|+13−1/11 = <span style={{ color: "#ff9944" }}>{ALPHA_INV_MONSTER.toFixed(4)}</span></div>
            <div>α⁻¹<sub>Feigen.</sub> = 2πδ²+... ≈ <span style={{ color: "#6699ff" }}>{ALPHA_INV_FEIGENBAUM.toFixed(4)}</span></div>
            <div>δ = {FEIGENBAUM_DELTA}</div>
            <div>T∞ = −ζ′(2)/2 = {T_INF}</div>
            <div>ψ(143) = 168 = |PSL(2,7)|</div>
            <div>dim S<sub>new</sub> = 11 = [1+4+6]</div>
            <div>genus(X₀(143)) = 13</div>
            <div style={{ marginTop: "6px", color: "#998866", fontStyle: "italic" }}>
              e<sup>π√163</sup> ≈ 640320³ + 744
            </div>
          </div>
        </div>
      </div>

      {/* Theory box */}
      <div style={{
        marginTop: "16px",
        padding: "14px",
        background: "rgba(255,200,50,0.03)",
        border: "1px solid rgba(255,200,50,0.12)",
        fontSize: "10px",
        lineHeight: "1.7",
        color: "#999",
        maxWidth: "1000px",
      }}>
        <div style={{ color: "#ffd866", fontSize: "10px", letterSpacing: "1px", marginBottom: "8px", textTransform: "uppercase" }}>
          How This Works
        </div>
        <p style={{ margin: "0 0 6px" }}>
          The <strong style={{ color: "#e0ddd5" }}>X₀(143) Burning Mandelbrot</strong> replaces
          the standard z→z²+c iteration with a <em>Hecke-modulated</em> iteration. At each step n, the
          Hecke trace Tr(T<sub>n</sub>) of the level-143 newspace determines the dynamics:
        </p>
        <p style={{ margin: "0 0 6px" }}>
          <span style={{ color: "#6699ff" }}>Tr(T<sub>n</sub>) ≥ 0</span> → analytic (Mandelbrot/bosonic) step: z<sub>n+1</sub> = z<sub>n</sub>² + c
          <br/>
          <span style={{ color: "#ff6644" }}>Tr(T<sub>n</sub>) &lt; 0</span> → folded (Burning Ship/fermionic) step: z<sub>n+1</sub> = (|Re z<sub>n</sub>| + i|Im z<sub>n</sub>|)² + c
        </p>
        <p style={{ margin: "0 0 6px" }}>
          The three Galois orbits (dim 1 + 4 + 6 = 11) of the newforms on X₀(143) encode the
          electron, muon, and tau sectors. The fold ratio at each pixel measures the
          <em> local fermionic fraction</em> of the vacuum — how much the iteration "burns."
        </p>
        <p style={{ margin: "0" }}>
          The Monster enters through α⁻¹ = ln|M| + genus − 1/11 = {ALPHA_INV_MONSTER.toFixed(4)},
          and Ramanujan's e<sup>π√163</sup> connects via the j-function whose coefficient 196884 = 196883+1
          obeys the moonshine residue map: 196884 ≡ 6 (mod 11), 196883 ≡ 11 (mod 13).
        </p>
      </div>
    </div>
  );
}
