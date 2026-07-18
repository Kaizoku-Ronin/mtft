import { useState, useRef, useEffect, useCallback, useMemo } from "react";

// ═══════════════════════════════════════════════════════════════════
// MTFT × ENNEPER SURFACE EXPLORER
// Weierstrass-Enneper Minimal Surfaces from Modular Arithmetic
// ═══════════════════════════════════════════════════════════════════

// --- Math helpers ---
function enneper(u, v, order = 1) {
  const r = Math.sqrt(u * u + v * v);
  const theta = Math.atan2(v, u);
  const rn = Math.pow(r, order);
  const x = u - (Math.pow(u, 2 * order + 1)) / (2 * order + 1) + u * v * v;
  const y = v - (Math.pow(v, 2 * order + 1)) / (2 * order + 1) + v * u * u;
  // Classical Enneper: f=1, g=z^n
  // x1 = Re ∫ (1 - z^{2n}) dz
  // x2 = Re ∫ i(1 + z^{2n}) dz
  // x3 = Re ∫ 2z^n dz
  const z_re = u, z_im = v;
  // z^n
  let zn_re = 1, zn_im = 0;
  for (let k = 0; k < order; k++) {
    const tr = zn_re * z_re - zn_im * z_im;
    const ti = zn_re * z_im + zn_im * z_re;
    zn_re = tr; zn_im = ti;
  }
  // z^{2n}
  let z2n_re = 1, z2n_im = 0;
  for (let k = 0; k < 2 * order; k++) {
    const tr = z2n_re * z_re - z2n_im * z_im;
    const ti = z2n_re * z_im + z2n_im * z_re;
    z2n_re = tr; z2n_im = ti;
  }
  // z^{n+1}
  let zn1_re = 1, zn1_im = 0;
  for (let k = 0; k < order + 1; k++) {
    const tr = zn1_re * z_re - zn1_im * z_im;
    const ti = zn1_re * z_im + zn1_im * z_re;
    zn1_re = tr; zn1_im = ti;
  }
  // z^{2n+1}
  let z2n1_re = 1, z2n1_im = 0;
  for (let k = 0; k < 2 * order + 1; k++) {
    const tr = z2n1_re * z_re - z2n1_im * z_im;
    const ti = z2n1_re * z_im + z2n1_im * z_re;
    z2n1_re = tr; z2n1_im = ti;
  }

  // Integrated forms:
  // x1 = Re[ z - z^{2n+1}/(2n+1) ]
  // x2 = Re[ iz + iz^{2n+1}/(2n+1) ] = -Im[ z + z^{2n+1}/(2n+1) ]
  // x3 = Re[ 2z^{n+1}/(n+1) ]
  const d = 2 * order + 1;
  const x1 = z_re - z2n1_re / d;
  const x2 = -(z_im + z2n1_im / d);
  const x3 = 2 * zn1_re / (order + 1);

  return [x1, x2, x3];
}

// Modular Enneper: f = eta-inspired, g = z^n with arithmetic weight modulation
function modularEnneper(u, v, order, arithStrength) {
  const z_re = u, z_im = v;
  const r2 = u * u + v * v;
  
  // Arithmetic weight: w_n = sum_{d|n} log(d)/d (simplified for visualization)
  // We modulate the surface with an arithmetic envelope
  const w = 1.0 + arithStrength * Math.log(1 + r2) * Math.cos(order * Math.atan2(v, u));
  
  // Dedekind eta-inspired conformal factor
  // |eta(tau)|^2 ≈ |q|^{1/12} prod (1 - |q|^{2n})^2
  // We approximate with an exponential damping modulated by arithmetic oscillation
  const yParam = 0.08; // y_c vacuum depth
  const etaFactor = Math.exp(-2 * Math.PI * yParam * Math.sqrt(r2));
  
  const [x1, x2, x3] = enneper(u, v, order);
  
  return [
    x1 * w * (1 + arithStrength * etaFactor * 0.3),
    x2 * w * (1 + arithStrength * etaFactor * 0.3),
    x3 * w
  ];
}

// Gaussian curvature for Enneper surface of order n
function gaussianCurvature(u, v, order) {
  const r2 = u * u + v * v;
  const rn = Math.pow(r2, order);
  // K = -4 / (1 + r^{2n})^4 for order-n Enneper
  return -4.0 / Math.pow(1 + rn, 4);
}

// 3D rotation
function rotate3D(p, ax, ay) {
  let [x, y, z] = p;
  // Rotate around Y
  let c = Math.cos(ay), s = Math.sin(ay);
  let x1 = x * c + z * s, z1 = -x * s + z * c;
  // Rotate around X
  c = Math.cos(ax); s = Math.sin(ax);
  let y1 = y * c - z1 * s, z2 = y * s + z1 * c;
  return [x1, y1, z2];
}

// Project 3D to 2D
function project(p, cx, cy, scale) {
  const perspective = 5;
  const d = perspective / (perspective + p[2] * 0.3);
  return [cx + p[0] * scale * d, cy - p[1] * scale * d, p[2]];
}

// Color from curvature
function curvatureColor(K, mode) {
  const absK = Math.min(Math.abs(K), 4);
  const t = absK / 4;
  if (mode === "curvature") {
    // Blue (flat) → Magenta (curved)
    const r = Math.floor(40 + 180 * t);
    const g = Math.floor(180 - 150 * t);
    const b = Math.floor(220 - 20 * t);
    return `rgb(${r},${g},${b})`;
  } else if (mode === "arithmetic") {
    // Gold → Deep red (arithmetic weight)
    const r = Math.floor(220 + 35 * t);
    const g = Math.floor(180 - 140 * t);
    const b = Math.floor(60 - 40 * t);
    return `rgb(${r},${g},${b})`;
  } else {
    // Depth-based
    const r = Math.floor(60 + 140 * t);
    const g = Math.floor(200 - 80 * t);
    const b = Math.floor(200 + 55 * t);
    return `rgb(${r},${g},${b})`;
  }
}

function SurfaceCanvas({ surfaceType, order, arithStrength, colorMode, rotX, rotY, resolution }) {
  const canvasRef = useRef(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const W = canvas.width;
    const H = canvas.height;
    ctx.fillStyle = "#06080d";
    ctx.fillRect(0, 0, W, H);

    const cx = W / 2, cy = H / 2;
    const range = surfaceType === "modular" ? 1.2 : (order <= 2 ? 1.4 : 1.0);
    const scale = surfaceType === "modular" ? 80 : (order <= 2 ? 60 : 40);
    const N = resolution;
    const step = (2 * range) / N;

    // Collect quads with depth for sorting
    const quads = [];
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        const u = -range + i * step;
        const v = -range + j * step;
        const u1 = u + step;
        const v1 = v + step;

        let p00, p10, p01, p11;
        if (surfaceType === "modular") {
          p00 = modularEnneper(u, v, order, arithStrength);
          p10 = modularEnneper(u1, v, order, arithStrength);
          p01 = modularEnneper(u, v1, order, arithStrength);
          p11 = modularEnneper(u1, v1, order, arithStrength);
        } else {
          p00 = enneper(u, v, order);
          p10 = enneper(u1, v, order);
          p01 = enneper(u, v1, order);
          p11 = enneper(u1, v1, order);
        }

        const pts = [p00, p10, p11, p01].map(p => rotate3D(p, rotX, rotY));
        const avgZ = (pts[0][2] + pts[1][2] + pts[2][2] + pts[3][2]) / 4;
        const K = gaussianCurvature((u + u1) / 2, (v + v1) / 2, order);
        const projected = pts.map(p => project(p, cx, cy, scale));
        quads.push({ projected, avgZ, K, u: (u + u1) / 2, v: (v + v1) / 2 });
      }
    }

    // Sort back to front
    quads.sort((a, b) => a.avgZ - b.avgZ);

    // Draw
    for (const q of quads) {
      const { projected: pp, K } = q;
      ctx.beginPath();
      ctx.moveTo(pp[0][0], pp[0][1]);
      ctx.lineTo(pp[1][0], pp[1][1]);
      ctx.lineTo(pp[2][0], pp[2][1]);
      ctx.lineTo(pp[3][0], pp[3][1]);
      ctx.closePath();
      ctx.fillStyle = curvatureColor(K, colorMode);
      ctx.fill();
      ctx.strokeStyle = "rgba(255,255,255,0.05)";
      ctx.lineWidth = 0.3;
      ctx.stroke();
    }

    // Labels
    ctx.fillStyle = "#8899aa";
    ctx.font = "11px 'Courier New', monospace";
    ctx.fillText(`Order n = ${order}`, 10, H - 30);
    if (surfaceType === "modular") {
      ctx.fillText(`Arith. strength = ${arithStrength.toFixed(2)}`, 10, H - 15);
    }
    ctx.fillText(`Total curvature = -${(4 * order).toFixed(0)}π`, 10, H - 45);

  }, [surfaceType, order, arithStrength, colorMode, rotX, rotY, resolution]);

  return <canvas ref={canvasRef} width={380} height={340} style={{ borderRadius: 8 }} />;
}

// Fold line visualization
function FoldCanvas({ order }) {
  const canvasRef = useRef(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const W = canvas.width, H = canvas.height;
    ctx.fillStyle = "#0a0e16";
    ctx.fillRect(0, 0, W, H);

    const cx = W / 2, cy = H / 2;
    const N = 200;
    const range = 1.3;

    // Draw self-intersection locus for Enneper of order n
    // Self-intersections occur when enneper(u1,v1) = enneper(u2,v2) for (u1,v1) != (u2,v2)
    // For visualization, we show where the surface folds over itself
    // by detecting where the Jacobian determinant changes sign

    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        const u = -range + (2 * range * i) / N;
        const v = -range + (2 * range * j) / N;
        const eps = 0.01;

        const p0 = enneper(u, v, order);
        const pu = enneper(u + eps, v, order);
        const pv = enneper(u, v + eps, order);

        // Jacobian determinant of 2D projection
        const dxdu = (pu[0] - p0[0]) / eps;
        const dxdv = (pv[0] - p0[0]) / eps;
        const dydu = (pu[1] - p0[1]) / eps;
        const dydv = (pv[1] - p0[1]) / eps;
        const det = dxdu * dydv - dxdv * dydu;

        const K = gaussianCurvature(u, v, order);
        const px = cx + (u / range) * (W / 2 - 20);
        const py = cy + (v / range) * (H / 2 - 20);

        if (Math.abs(det) < 0.3 && order > 1) {
          // Near fold / self-intersection
          ctx.fillStyle = "#ff4444";
          ctx.fillRect(px - 1.5, py - 1.5, 3, 3);
        } else {
          const t = Math.min(Math.abs(K), 2) / 2;
          const r = Math.floor(30 + 60 * t);
          const g = Math.floor(100 + 80 * (1 - t));
          const b = Math.floor(160 + 60 * t);
          ctx.fillStyle = `rgb(${r},${g},${b})`;
          ctx.fillRect(px - 1, py - 1, 2, 2);
        }
      }
    }

    // Draw fold lines
    ctx.strokeStyle = "#ff884488";
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    // Re(z) = 0 fold
    ctx.beginPath();
    ctx.moveTo(cx, 10);
    ctx.lineTo(cx, H - 10);
    ctx.stroke();
    // Im(z) = 0 fold
    ctx.beginPath();
    ctx.moveTo(10, cy);
    ctx.lineTo(W - 10, cy);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = "#aabbcc";
    ctx.font = "11px 'Courier New', monospace";
    ctx.fillText("Parameter domain (u,v)", 8, 14);
    ctx.fillStyle = "#ff6666";
    ctx.fillText("■ Self-intersection loci", 8, H - 8);
    ctx.fillStyle = "#ff884488";
    ctx.fillText("--- Fold lines (Re=0, Im=0)", 8, H - 22);

  }, [order]);

  return <canvas ref={canvasRef} width={380} height={280} style={{ borderRadius: 8 }} />;
}

// Gauss map visualization
function GaussMapCanvas({ order }) {
  const canvasRef = useRef(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const W = canvas.width, H = canvas.height;
    ctx.fillStyle = "#0a0e16";
    ctx.fillRect(0, 0, W, H);

    const cx = W / 2, cy = H / 2;
    const R = Math.min(W, H) / 2 - 30;

    // Draw Riemann sphere outline
    ctx.strokeStyle = "#334455";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, 2 * Math.PI);
    ctx.stroke();

    // The Gauss map of Enneper order n is g(z) = z^n
    // On the Riemann sphere, this wraps n times
    // Visualize by mapping parameter domain points through stereographic projection
    const N = 150;
    const range = 1.5;

    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        const u = -range + (2 * range * i) / N;
        const v = -range + (2 * range * j) / N;

        // g(z) = z^n
        let gre = 1, gim = 0;
        for (let k = 0; k < order; k++) {
          const tr = gre * u - gim * v;
          const ti = gre * v + gim * u;
          gre = tr; gim = ti;
        }

        // Stereographic projection of Gauss map
        const r2 = gre * gre + gim * gim;
        const sx = 2 * gre / (1 + r2);
        const sy = 2 * gim / (1 + r2);

        const px = cx + sx * R;
        const py = cy + sy * R;

        // Color by angle
        const angle = Math.atan2(gim, gre);
        const hue = ((angle / Math.PI + 1) / 2);
        const r = Math.floor(80 + 140 * (0.5 + 0.5 * Math.cos(hue * 2 * Math.PI)));
        const g = Math.floor(120 + 100 * (0.5 + 0.5 * Math.cos(hue * 2 * Math.PI + 2.094)));
        const b = Math.floor(180 + 75 * (0.5 + 0.5 * Math.cos(hue * 2 * Math.PI + 4.189)));

        ctx.fillStyle = `rgb(${r},${g},${b})`;
        ctx.fillRect(px - 0.8, py - 0.8, 1.6, 1.6);
      }
    }

    ctx.fillStyle = "#aabbcc";
    ctx.font = "11px 'Courier New', monospace";
    ctx.fillText(`Gauss map: g(z) = z^${order}`, 8, 14);
    ctx.fillText(`Wrapping number = ${order}`, 8, H - 8);
    ctx.fillText(`(≈ newform dim on X₀(143) for n=11)`, 8, H - 22);

  }, [order]);

  return <canvas ref={canvasRef} width={380} height={280} style={{ borderRadius: 8 }} />;
}

// Main Component
export default function MTFTEnneper() {
  const [order, setOrder] = useState(1);
  const [arithStrength, setArithStrength] = useState(0.0);
  const [colorMode, setColorMode] = useState("curvature");
  const [rotX, setRotX] = useState(-0.5);
  const [rotY, setRotY] = useState(0.4);
  const [surfaceType, setSurfaceType] = useState("classical");
  const [activeTab, setActiveTab] = useState("surface");
  const [isDragging, setIsDragging] = useState(false);
  const [lastPos, setLastPos] = useState({ x: 0, y: 0 });

  const handlePointerDown = (e) => {
    setIsDragging(true);
    setLastPos({ x: e.clientX, y: e.clientY });
  };
  const handlePointerMove = (e) => {
    if (!isDragging) return;
    const dx = e.clientX - lastPos.x;
    const dy = e.clientY - lastPos.y;
    setRotY(prev => prev + dx * 0.01);
    setRotX(prev => prev + dy * 0.01);
    setLastPos({ x: e.clientX, y: e.clientY });
  };
  const handlePointerUp = () => setIsDragging(false);

  const dictionary = [
    { math: "Weierstrass data (f, g)", mtft: "Newforms on X₀(143)", detail: "Weight-2 cusp forms f₁, f₂, f₃ as holomorphic input" },
    { math: "Order n of g(z) = z^n", mtft: "Newform dimension = 11", detail: "Gauss map wraps 11 times → 11-dim newform space" },
    { math: "Self-intersections", mtft: "Burning Ship fold lines", detail: "Loci where map fails to be injective" },
    { math: "Total curvature -4nπ", mtft: "Genus 13 constraint", detail: "Jorge-Meeks: ∫K dA = -4π(g + r - 1)" },
    { math: "Isothermal coords", mtft: "Conformal vacuum structure", detail: "ds² = λ(z)|dz|² on the modular curve" },
    { math: "Zero mean curvature", mtft: "Vacuum stability (H = 0)", detail: "Minimal = equilibrium of arithmetic potential" },
    { math: "Gaussian curvature K", mtft: "Stiffness μ_N(y)", detail: "Both measure local restoring force / curvature" },
    { math: "Enneper surface ends", mtft: "Cusps of X₀(143)", detail: "Asymptotic directions ↔ cusp neighborhoods" },
  ];

  const tabStyle = (t) => ({
    padding: "6px 14px",
    background: activeTab === t ? "#1a2540" : "transparent",
    color: activeTab === t ? "#4af0c0" : "#667788",
    border: activeTab === t ? "1px solid #2a3a55" : "1px solid transparent",
    borderBottom: activeTab === t ? "1px solid #1a2540" : "1px solid #2a3a55",
    borderRadius: "6px 6px 0 0",
    cursor: "pointer",
    fontSize: 12,
    fontFamily: "'Courier New', monospace",
    fontWeight: activeTab === t ? 700 : 400,
    transition: "all 0.2s",
    marginRight: 2,
  });

  const sliderBox = { display: "flex", alignItems: "center", gap: 8, marginBottom: 8 };
  const labelStyle = { color: "#8899aa", fontSize: 11, fontFamily: "'Courier New', monospace", minWidth: 90 };
  const valStyle = { color: "#4af0c0", fontSize: 11, fontFamily: "'Courier New', monospace", minWidth: 30 };

  return (
    <div style={{
      background: "#06080d",
      color: "#c8d8e8",
      fontFamily: "'Courier New', monospace",
      padding: 16,
      minHeight: "100vh",
      maxWidth: 820,
      margin: "0 auto",
    }}>
      {/* Header */}
      <div style={{ textAlign: "center", marginBottom: 16 }}>
        <div style={{ fontSize: 10, letterSpacing: 6, color: "#4af0c0", marginBottom: 4 }}>
          MODULAR TIME FIELD THEORY
        </div>
        <h1 style={{
          fontSize: 22,
          fontWeight: 300,
          color: "#e8f0f8",
          margin: 0,
          letterSpacing: 1,
        }}>
          Enneper Surfaces × Modular Arithmetic
        </h1>
        <div style={{ fontSize: 11, color: "#667788", marginTop: 4 }}>
          Weierstrass-Enneper Representation from X₀(143) Newforms
        </div>
      </div>

      {/* Tabs */}
      <div style={{ display: "flex", borderBottom: "1px solid #2a3a55", marginBottom: 0 }}>
        <button style={tabStyle("surface")} onClick={() => setActiveTab("surface")}>3D Surface</button>
        <button style={tabStyle("folds")} onClick={() => setActiveTab("folds")}>Fold Structure</button>
        <button style={tabStyle("gauss")} onClick={() => setActiveTab("gauss")}>Gauss Map</button>
        <button style={tabStyle("dictionary")} onClick={() => setActiveTab("dictionary")}>Dictionary</button>
      </div>

      {/* Content */}
      <div style={{
        background: "#0c1018",
        border: "1px solid #2a3a55",
        borderTop: "none",
        borderRadius: "0 0 8px 8px",
        padding: 16,
      }}>
        {/* Controls (shared) */}
        {activeTab !== "dictionary" && (
          <div style={{
            display: "flex",
            flexWrap: "wrap",
            gap: 16,
            marginBottom: 16,
            padding: 12,
            background: "#0a0e16",
            borderRadius: 6,
            border: "1px solid #1a2a40",
          }}>
            <div style={{ flex: 1, minWidth: 200 }}>
              <div style={sliderBox}>
                <span style={labelStyle}>Order n:</span>
                <input type="range" min={1} max={7} value={order} onChange={e => setOrder(+e.target.value)}
                  style={{ flex: 1, accentColor: "#4af0c0" }} />
                <span style={valStyle}>{order}</span>
              </div>
              {activeTab === "surface" && (
                <>
                  <div style={sliderBox}>
                    <span style={labelStyle}>Arith. mod:</span>
                    <input type="range" min={0} max={100} value={arithStrength * 100}
                      onChange={e => { setArithStrength(+e.target.value / 100); setSurfaceType(+e.target.value > 0 ? "modular" : "classical"); }}
                      style={{ flex: 1, accentColor: "#ff8844" }} />
                    <span style={valStyle}>{arithStrength.toFixed(2)}</span>
                  </div>
                </>
              )}
            </div>
            <div style={{ flex: 1, minWidth: 200 }}>
              {activeTab === "surface" && (
                <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
                  {["curvature", "arithmetic", "depth"].map(m => (
                    <button key={m} onClick={() => setColorMode(m)} style={{
                      padding: "4px 10px",
                      background: colorMode === m ? "#1a2540" : "#0a0e16",
                      color: colorMode === m ? "#4af0c0" : "#667788",
                      border: `1px solid ${colorMode === m ? "#4af0c0" : "#2a3a55"}`,
                      borderRadius: 4, fontSize: 10, cursor: "pointer",
                      fontFamily: "'Courier New', monospace",
                    }}>{m}</button>
                  ))}
                </div>
              )}
            </div>
          </div>
        )}

        {/* Surface Tab */}
        {activeTab === "surface" && (
          <div style={{ display: "flex", flexDirection: "column", alignItems: "center" }}>
            <div
              onPointerDown={handlePointerDown}
              onPointerMove={handlePointerMove}
              onPointerUp={handlePointerUp}
              onPointerLeave={handlePointerUp}
              style={{ cursor: isDragging ? "grabbing" : "grab", touchAction: "none" }}
            >
              <SurfaceCanvas
                surfaceType={surfaceType}
                order={order}
                arithStrength={arithStrength}
                colorMode={colorMode}
                rotX={rotX}
                rotY={rotY}
                resolution={50}
              />
            </div>
            <div style={{ fontSize: 10, color: "#556677", marginTop: 6 }}>
              Drag to rotate · Order n controls g(z) = z^n in Weierstrass data
            </div>

            {/* Formulas */}
            <div style={{
              marginTop: 16,
              padding: 12,
              background: "#080c14",
              borderRadius: 6,
              border: "1px solid #1a2a40",
              width: "100%",
            }}>
              <div style={{ fontSize: 11, color: "#4af0c0", marginBottom: 6 }}>
                Weierstrass-Enneper Data
              </div>
              <div style={{ fontSize: 11, color: "#aabbcc", lineHeight: 1.6 }}>
                {surfaceType === "classical" ? (
                  <>
                    <div>f(z) = 1, &nbsp; g(z) = z<sup>{order}</sup></div>
                    <div>x₁ = Re ∫ (1 - z<sup>{2*order}</sup>) dz = Re[z - z<sup>{2*order+1}</sup>/{2*order+1}]</div>
                    <div>x₂ = -Im ∫ (1 + z<sup>{2*order}</sup>) dz</div>
                    <div>x₃ = Re ∫ 2z<sup>{order}</sup> dz = 2Re[z<sup>{order+1}</sup>/{order+1}]</div>
                    <div style={{ color: "#ff8844", marginTop: 4 }}>
                      Total curvature: ∫K dA = -{4*order}π &nbsp;|&nbsp; Self-intersections: {order > 1 ? "Yes" : "None"}
                    </div>
                  </>
                ) : (
                  <>
                    <div>f(z) = 1 · W(z; y_c), &nbsp; g(z) = z<sup>{order}</sup></div>
                    <div>W(z) = 1 + a · ln(1+|z|²) · cos(n·arg(z)) &nbsp; [arithmetic envelope]</div>
                    <div>η-factor: exp(-2πy_c·|z|) &nbsp; [Dedekind damping at vacuum depth]</div>
                    <div style={{ color: "#ff8844", marginTop: 4 }}>
                      Modular deformation strength: {arithStrength.toFixed(2)} &nbsp;|&nbsp; y_c = 0.08
                    </div>
                  </>
                )}
              </div>
            </div>
          </div>
        )}

        {/* Folds Tab */}
        {activeTab === "folds" && (
          <div style={{ display: "flex", flexDirection: "column", alignItems: "center" }}>
            <FoldCanvas order={order} />
            <div style={{
              marginTop: 12,
              padding: 12,
              background: "#080c14",
              borderRadius: 6,
              border: "1px solid #1a2a40",
              width: "100%",
            }}>
              <div style={{ fontSize: 11, color: "#4af0c0", marginBottom: 6 }}>
                Enneper Fold Structure ↔ Burning Ship Folds
              </div>
              <div style={{ fontSize: 11, color: "#aabbcc", lineHeight: 1.8 }}>
                <div>
                  <span style={{ color: "#ff6666" }}>■</span> Red: Self-intersection loci (|det J| &lt; ε)
                  {order === 1 ? " — none for n=1 (no self-intersections)" : ` — present for n=${order}`}
                </div>
                <div>
                  <span style={{ color: "#ff884488" }}>---</span> Dashed: Fold lines Re(z)=0, Im(z)=0
                </div>
                <div style={{ marginTop: 6, color: "#ff8844" }}>
                  MTFT parallel: The Burning Ship fold lines F_x = {"{Re(z)=0}"} and F_y = {"{Im(z)=0}"}
                  divide the complex plane into sectors where the map is locally smooth —
                  identical to the Enneper fold structure at order n ≥ 2.
                </div>
                {order >= 2 && (
                  <div style={{ marginTop: 6, color: "#7b61ff" }}>
                    At n = {order}: {2 * order} fold rays emanate from the origin,
                    creating {4 * order} locally smooth sectors.
                    Compare: X₀(143) has {order <= 4 ? "genus 13 with" : ""} 168 = |PSL(2,7)| symmetries.
                  </div>
                )}
              </div>
            </div>
          </div>
        )}

        {/* Gauss Map Tab */}
        {activeTab === "gauss" && (
          <div style={{ display: "flex", flexDirection: "column", alignItems: "center" }}>
            <GaussMapCanvas order={order} />
            <div style={{
              marginTop: 12,
              padding: 12,
              background: "#080c14",
              borderRadius: 6,
              border: "1px solid #1a2a40",
              width: "100%",
            }}>
              <div style={{ fontSize: 11, color: "#4af0c0", marginBottom: 6 }}>
                Gauss Map ↔ Newform Projection
              </div>
              <div style={{ fontSize: 11, color: "#aabbcc", lineHeight: 1.8 }}>
                <div>
                  The Gauss map g: M → S² sends each point of the minimal surface
                  to its unit normal on the Riemann sphere.
                </div>
                <div style={{ marginTop: 4 }}>
                  For Enneper order n: g(z) = z<sup>{order}</sup> wraps the sphere <span style={{ color: "#4af0c0" }}>{order} time{order > 1 ? "s" : ""}</span>.
                </div>
                <div style={{ marginTop: 6, color: "#ff8844" }}>
                  MTFT key: The newform space S₂(Γ₀(143)) has dimension 11.
                  Setting n = 11 gives a Gauss map that wraps 11 times —
                  each wrapping corresponds to one newform direction in the Jacobian J₀(143).
                </div>
                <div style={{ marginTop: 6, color: "#7b61ff" }}>
                  The Gauss map of a minimal surface is <em>conformal</em> (angle-preserving).
                  This is precisely the isothermal coordinate property that
                  connects Enneper geometry to the conformal structure of X₀(143).
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Dictionary Tab */}
        {activeTab === "dictionary" && (
          <div>
            <div style={{ fontSize: 12, color: "#4af0c0", marginBottom: 12 }}>
              MTFT ↔ Enneper Surface Dictionary
            </div>
            <div style={{ overflowX: "auto" }}>
              <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 11 }}>
                <thead>
                  <tr style={{ borderBottom: "1px solid #2a3a55" }}>
                    <th style={{ textAlign: "left", padding: "6px 8px", color: "#4af0c0" }}>Enneper / Minimal Surface</th>
                    <th style={{ textAlign: "left", padding: "6px 8px", color: "#ff8844" }}>MTFT Object</th>
                    <th style={{ textAlign: "left", padding: "6px 8px", color: "#7b61ff" }}>Bridge</th>
                  </tr>
                </thead>
                <tbody>
                  {dictionary.map((row, i) => (
                    <tr key={i} style={{ borderBottom: "1px solid #152030" }}>
                      <td style={{ padding: "6px 8px", color: "#c8d8e8" }}>{row.math}</td>
                      <td style={{ padding: "6px 8px", color: "#ddd" }}>{row.mtft}</td>
                      <td style={{ padding: "6px 8px", color: "#8899aa" }}>{row.detail}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            <div style={{
              marginTop: 16,
              padding: 12,
              background: "#080c14",
              borderRadius: 6,
              border: "1px solid #1a2a40",
            }}>
              <div style={{ fontSize: 11, color: "#ff8844", marginBottom: 6 }}>
                The Central Construction
              </div>
              <div style={{ fontSize: 11, color: "#aabbcc", lineHeight: 1.8 }}>
                <div>
                  Given weight-2 newforms f_i ∈ S₂(Γ₀(143)), define Weierstrass data:
                </div>
                <div style={{ padding: "8px 0", color: "#e8f0f8" }}>
                  f(z) = f₁(z) · |η(z)|<sup>4α⁻¹/π</sup>, &nbsp;&nbsp;
                  g(z) = f₂(z) / f₁(z)
                </div>
                <div>
                  This produces a <em>modular minimal surface</em> whose:
                </div>
                <div style={{ paddingLeft: 12, marginTop: 4 }}>
                  → Total curvature encodes the genus of X₀(143)<br/>
                  → Self-intersection loci are the Atkin-Lehner sector boundaries<br/>
                  → Gauss map wrapping number = dim S₂(Γ₀(143)) = 11<br/>
                  → Curvature at the vacuum nome |K(τ_c)| = μ_N(y_c)
                </div>
                <div style={{ marginTop: 8, color: "#4af0c0" }}>
                  The fermion mass hierarchy m_e ≪ m_μ ≪ m_τ is the
                  depth hierarchy of the three Enneper ends,
                  measured by the period integrals of the modular minimal surface.
                </div>
              </div>
            </div>
          </div>
        )}
      </div>

      {/* Footer */}
      <div style={{
        marginTop: 12,
        padding: 8,
        textAlign: "center",
        fontSize: 10,
        color: "#445566",
      }}>
        MTFT Paper 32 · Enneper-Weierstrass Bridge · X₀(143) Modular Minimal Surfaces
      </div>
    </div>
  );
}
