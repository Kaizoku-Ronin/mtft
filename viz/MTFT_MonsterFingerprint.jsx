import { useState, useEffect, useRef, useCallback } from "react";

// ═══════════════════════════════════════════════════════════════
// MONSTERHASH VISUAL FINGERPRINT
// Every seed produces a unique, deterministic geometric identity.
// Change one character → ~50% of the visual changes (avalanche).
// 13-round SL(2,ℤ)-inspired mixing — MTFT structural constants.
// ═══════════════════════════════════════════════════════════════

const SZ = 280;
const CX = SZ / 2;
const R = SZ / 2 - 12;

const CYAN  = [0, 245, 212];
const BG    = [8, 8, 14];

function monsterHash(input) {
  const bytes = new TextEncoder().encode(input || "MTFT");
  const s = new Uint32Array([
    0x4D544654, 0x8F1BBCDC, 0x00000143, 0x0000000D,
    0xA953FD4E, 0x5C4DD124, 0xE49B69C1, 0x2DE92C6F,
  ]);
  for (const b of bytes) {
    for (let r = 0; r < 13; r++) {
      const i = r & 7;
      s[i] = (s[i] + b + r + 1) | 0;
      s[i] = Math.imul(s[i], 0x5BD1E995) | 0;
      s[i] = s[i] ^ (s[i] >>> 15);
      s[(i + 1) & 7] = (s[(i + 1) & 7] ^ s[i]) | 0;
    }
  }
  for (let r = 0; r < 13; r++) {
    for (let i = 0; i < 8; i++) {
      s[i] = (s[i] ^ s[(i + 3) & 7]) | 0;
      s[i] = Math.imul(s[i], 0x1B873593) | 0;
      s[i] = (s[i] ^ (s[i] >>> 13)) | 0;
      s[(i + 5) & 7] = (s[(i + 5) & 7] + s[i]) | 0;
    }
  }
  return s;
}

function hashHex(s) {
  return Array.from(s).map(v => (v >>> 0).toString(16).padStart(8, "0")).join("");
}

function byte(s, idx) { return ((s[idx >> 2] >>> ((3 - (idx & 3)) * 8)) & 0xFF); }
function norm(s, idx) { return byte(s, idx) / 255; }
function rangeInt(s, idx, min, max) { return min + (byte(s, idx) % (max - min + 1)); }

function hsl(h, s, l) {
  h = ((h % 360) + 360) % 360;
  const a = s * Math.min(l, 1 - l);
  const f = (n) => {
    const k = (n + h / 30) % 12;
    return l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1));
  };
  return `rgb(${Math.round(f(0)*255)},${Math.round(f(8)*255)},${Math.round(f(4)*255)})`;
}

function extractParams(s) {
  const symmetry = rangeInt(s, 0, 3, 13);
  const hue1 = byte(s, 1) * 1.41;
  const hue2 = (hue1 + 60 + byte(s, 2) * 0.94) % 360;
  const hue3 = (hue1 + 180 + byte(s, 3) * 0.47) % 360;
  const layers = rangeInt(s, 4, 3, 6);
  const rotation = norm(s, 5) * Math.PI * 2;
  const outerShape = rangeInt(s, 6, 0, 4); // 0=circle 1=sq 2=tri 3=diamond 4=hex
  const outerCount = symmetry * rangeInt(s, 7, 1, 3);
  const petalWidth = 0.15 + norm(s, 8) * 0.35;
  const petalLength = 0.3 + norm(s, 9) * 0.35;
  const innerRings = rangeInt(s, 10, 1, 3);
  const coreShape = rangeInt(s, 11, 3, 8);
  const coreFill = norm(s, 12) > 0.5;
  const spokeStyle = rangeInt(s, 13, 0, 3);
  const dotSize = 1.5 + norm(s, 14) * 3;
  const webSkip = rangeInt(s, 15, 1, Math.max(2, symmetry - 1));
  const glowIntensity = 0.1 + norm(s, 16) * 0.3;
  const sat1 = 0.5 + norm(s, 17) * 0.5;
  const sat2 = 0.4 + norm(s, 18) * 0.5;
  const lineWeight = 1 + norm(s, 19) * 1.5;

  return {
    symmetry, hue1, hue2, hue3, layers, rotation, outerShape, outerCount,
    petalWidth, petalLength, innerRings, coreShape, coreFill, spokeStyle,
    dotSize, webSkip, glowIntensity, sat1, sat2, lineWeight,
  };
}

function drawFingerprint(canvas, params) {
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const { symmetry: N, hue1, hue2, hue3, rotation, outerCount,
          petalWidth, petalLength, innerRings, coreShape, coreFill,
          spokeStyle, dotSize, webSkip, glowIntensity, sat1, sat2, lineWeight } = params;

  ctx.clearRect(0, 0, SZ, SZ);

  // Background circle
  const bgGrad = ctx.createRadialGradient(CX, CX, 0, CX, CX, R + 10);
  bgGrad.addColorStop(0, hsl(hue1, 0.15, 0.08));
  bgGrad.addColorStop(0.7, hsl(hue1, 0.1, 0.04));
  bgGrad.addColorStop(1, "transparent");
  ctx.fillStyle = bgGrad;
  ctx.beginPath();
  ctx.arc(CX, CX, R + 8, 0, Math.PI * 2);
  ctx.fill();

  // Glow ring
  ctx.strokeStyle = hsl(hue1, 0.8, 0.3 + glowIntensity);
  ctx.lineWidth = 1;
  ctx.globalAlpha = 0.3;
  ctx.beginPath();
  ctx.arc(CX, CX, R, 0, Math.PI * 2);
  ctx.stroke();
  ctx.globalAlpha = 1;

  // Layer 1: Connection web
  ctx.save();
  ctx.translate(CX, CX);
  ctx.rotate(rotation);
  const webR = R * 0.85;
  ctx.strokeStyle = hsl(hue2, sat2, 0.35);
  ctx.lineWidth = lineWeight * 0.5;
  ctx.globalAlpha = 0.25;
  for (let i = 0; i < outerCount; i++) {
    const a1 = (i / outerCount) * Math.PI * 2;
    const j = (i + webSkip) % outerCount;
    const a2 = (j / outerCount) * Math.PI * 2;
    ctx.beginPath();
    ctx.moveTo(Math.cos(a1) * webR, Math.sin(a1) * webR);
    ctx.lineTo(Math.cos(a2) * webR, Math.sin(a2) * webR);
    ctx.stroke();
  }
  ctx.globalAlpha = 1;

  // Layer 2: Outer elements
  const outerR = R * 0.82;
  for (let i = 0; i < outerCount; i++) {
    const angle = (i / outerCount) * Math.PI * 2;
    const ox = Math.cos(angle) * outerR;
    const oy = Math.sin(angle) * outerR;
    ctx.fillStyle = i % N === 0
      ? hsl(hue1, sat1, 0.55)
      : hsl(hue2, sat2 * 0.7, 0.4);
    const sz = i % N === 0 ? dotSize * 1.4 : dotSize;
    drawShape(ctx, ox, oy, sz, params.outerShape, angle);
  }

  // Layer 3: Petals
  const petalR = R * 0.62;
  ctx.lineWidth = lineWeight;
  for (let i = 0; i < N; i++) {
    const angle = (i / N) * Math.PI * 2;
    const px = Math.cos(angle) * petalR;
    const py = Math.sin(angle) * petalR;
    const pLen = petalR * petalLength;
    const pWid = petalR * petalWidth;

    ctx.save();
    ctx.translate(px * 0.1, py * 0.1);
    ctx.rotate(angle);
    ctx.beginPath();
    ctx.moveTo(petalR * 0.3, 0);
    ctx.quadraticCurveTo(petalR * 0.5 + pLen * 0.3, -pWid * 30, petalR * 0.3 + pLen * 60, 0);
    ctx.quadraticCurveTo(petalR * 0.5 + pLen * 0.3, pWid * 30, petalR * 0.3, 0);
    ctx.closePath();
    ctx.fillStyle = hsl(hue3, sat1 * 0.8, 0.3);
    ctx.globalAlpha = 0.4;
    ctx.fill();
    ctx.globalAlpha = 1;
    ctx.strokeStyle = hsl(hue3, sat1, 0.5);
    ctx.lineWidth = lineWeight * 0.7;
    ctx.stroke();
    ctx.restore();
  }

  // Layer 4: Inner rings
  for (let ring = 0; ring < innerRings; ring++) {
    const iR = R * (0.35 - ring * 0.08);
    const iCount = N * (ring + 1);
    ctx.strokeStyle = hsl(hue1, sat1, 0.4 + ring * 0.1);
    ctx.lineWidth = lineWeight * 0.6;
    ctx.globalAlpha = 0.5;
    ctx.beginPath();
    ctx.arc(0, 0, iR, 0, Math.PI * 2);
    ctx.stroke();
    ctx.globalAlpha = 1;

    for (let i = 0; i < iCount; i++) {
      const angle = (i / iCount) * Math.PI * 2;
      const ix = Math.cos(angle) * iR;
      const iy = Math.sin(angle) * iR;
      ctx.fillStyle = hsl(hue1, sat1, 0.5);
      ctx.beginPath();
      ctx.arc(ix, iy, dotSize * 0.5, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  // Layer 5: Spokes
  if (spokeStyle > 0) {
    ctx.strokeStyle = hsl(hue2, sat2, 0.4);
    ctx.lineWidth = lineWeight * 0.4;
    ctx.globalAlpha = 0.3;
    const spokeCount = spokeStyle === 1 ? N : spokeStyle === 2 ? N * 2 : N * 3;
    for (let i = 0; i < spokeCount; i++) {
      const angle = (i / spokeCount) * Math.PI * 2;
      ctx.beginPath();
      ctx.moveTo(Math.cos(angle) * R * 0.15, Math.sin(angle) * R * 0.15);
      ctx.lineTo(Math.cos(angle) * R * 0.55, Math.sin(angle) * R * 0.55);
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
  }

  // Layer 6: Core
  const coreR = R * 0.18;
  ctx.fillStyle = hsl(hue1, sat1, coreFill ? 0.45 : 0.1);
  ctx.strokeStyle = hsl(hue1, sat1, 0.6);
  ctx.lineWidth = lineWeight;
  ctx.beginPath();
  for (let i = 0; i <= coreShape; i++) {
    const angle = (i / coreShape) * Math.PI * 2 - Math.PI / 2;
    const r = coreR * (coreFill ? 1 : (i % 2 === 0 ? 1 : 0.5));
    const px = Math.cos(angle) * r;
    const py = Math.sin(angle) * r;
    i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
  }
  ctx.closePath();
  ctx.fill();
  ctx.stroke();

  // Core dot
  ctx.fillStyle = hsl(hue1, 1, 0.7);
  ctx.beginPath();
  ctx.arc(0, 0, dotSize * 0.7, 0, Math.PI * 2);
  ctx.fill();

  ctx.restore();
}

function drawShape(ctx, x, y, size, shape, angle) {
  ctx.save();
  ctx.translate(x, y);
  ctx.rotate(angle);
  ctx.beginPath();
  if (shape === 0) {
    ctx.arc(0, 0, size, 0, Math.PI * 2);
  } else {
    const sides = shape === 1 ? 4 : shape === 2 ? 3 : shape === 3 ? 4 : 6;
    const startAngle = shape === 3 ? Math.PI / 4 : -Math.PI / 2;
    for (let i = 0; i <= sides; i++) {
      const a = startAngle + (i / sides) * Math.PI * 2;
      const px = Math.cos(a) * size, py = Math.sin(a) * size;
      i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
    }
    ctx.closePath();
  }
  ctx.fill();
  ctx.restore();
}

function avalancheDiff(h1, h2) {
  let bits = 0, total = 0;
  for (let i = 0; i < 8; i++) {
    const xor = (h1[i] ^ h2[i]) >>> 0;
    for (let b = 0; b < 32; b++) {
      if ((xor >>> b) & 1) bits++;
      total++;
    }
  }
  return { flipped: bits, total, pct: ((bits / total) * 100).toFixed(1) };
}

export default function MonsterFingerprint() {
  const [seed, setSeed] = useState("MonsterChain");
  const canvasA = useRef(null);
  const canvasB = useRef(null);
  const [avalanche, setAvalanche] = useState(null);
  const [paramsA, setParamsA] = useState(null);
  const [paramsB, setParamsB] = useState(null);
  const [hexA, setHexA] = useState("");
  const [hexB, setHexB] = useState("");

  const render = useCallback(() => {
    const hashA = monsterHash(seed);
    const modified = seed.length > 0
      ? seed.slice(0, -1) + String.fromCharCode(seed.charCodeAt(seed.length - 1) + 1)
      : "a";
    const hashB = monsterHash(modified);

    const pA = extractParams(hashA);
    const pB = extractParams(hashB);
    setParamsA(pA);
    setParamsB(pB);
    setHexA(hashHex(hashA));
    setHexB(hashHex(hashB));
    setAvalanche(avalancheDiff(hashA, hashB));

    drawFingerprint(canvasA.current, pA);
    drawFingerprint(canvasB.current, pB);
  }, [seed]);

  useEffect(() => { render(); }, [render]);

  const paramLabels = paramsA ? [
    ["Symmetry", paramsA.symmetry, paramsB?.symmetry],
    ["Primary Hue", Math.round(paramsA.hue1) + "°", Math.round(paramsB?.hue1) + "°"],
    ["Secondary Hue", Math.round(paramsA.hue2) + "°", Math.round(paramsB?.hue2) + "°"],
    ["Core Shape", paramsA.coreShape + "-gon", paramsB?.coreShape + "-gon"],
    ["Outer Elements", paramsA.outerCount, paramsB?.outerCount],
    ["Spoke Style", paramsA.spokeStyle, paramsB?.spokeStyle],
  ] : [];

  const modified = seed.length > 0
    ? seed.slice(0, -1) + String.fromCharCode(seed.charCodeAt(seed.length - 1) + 1)
    : "a";

  return (
    <div style={{
      background: `rgb(${BG.join(",")})`, minHeight: "100vh",
      fontFamily: "'JetBrains Mono', 'Courier New', monospace",
      color: "#d0d0d0", padding: "16px 8px",
      display: "flex", flexDirection: "column", alignItems: "center",
    }}>
      <h1 style={{
        fontSize: 17, fontWeight: 700, letterSpacing: 3, margin: 0,
        color: `rgb(${CYAN.join(",")})`,
        textShadow: `0 0 20px rgba(${CYAN.join(",")},0.4)`,
      }}>
        MONSTERHASH FINGERPRINT
      </h1>
      <div style={{ fontSize: 9, color: "#555", letterSpacing: 1, marginTop: 4, marginBottom: 14 }}>
        DETERMINISTIC GEOMETRIC IDENTITY — 13-ROUND SL(2,ℤ) MIXING — 256-BIT AVALANCHE
      </div>

      {/* Input */}
      <div style={{ width: "100%", maxWidth: 620, marginBottom: 16 }}>
        <label style={{ fontSize: 9, color: `rgb(${CYAN.join(",")})`, letterSpacing: 2, display: "block", marginBottom: 4 }}>
          SEED / ADDRESS / WALLET ID
        </label>
        <input
          type="text"
          value={seed}
          onChange={(e) => setSeed(e.target.value)}
          placeholder="Enter any text..."
          style={{
            width: "100%", boxSizing: "border-box",
            background: "rgba(255,255,255,0.04)",
            border: `1px solid rgba(${CYAN.join(",")},0.3)`,
            borderRadius: 6, padding: "10px 14px",
            color: "#eee", fontSize: 14,
            fontFamily: "'JetBrains Mono', monospace",
            outline: "none",
          }}
        />
      </div>

      {/* Fingerprints side by side */}
      <div style={{
        display: "flex", gap: 20, alignItems: "flex-start", flexWrap: "wrap",
        justifyContent: "center", marginBottom: 16,
      }}>
        {/* Identity */}
        <div style={{ textAlign: "center" }}>
          <div style={{ fontSize: 9, letterSpacing: 2, color: `rgb(${CYAN.join(",")})`, marginBottom: 6 }}>
            IDENTITY
          </div>
          <canvas ref={canvasA} width={SZ} height={SZ} style={{
            borderRadius: "50%",
            boxShadow: `0 0 30px rgba(${CYAN.join(",")},0.12)`,
            maxWidth: "min(280px, 42vw)", height: "auto",
          }} />
          <div style={{ fontSize: 8, color: "#444", marginTop: 6, wordBreak: "break-all", maxWidth: 260 }}>
            "{seed}"
          </div>
        </div>

        {/* Avalanche indicator */}
        <div style={{
          display: "flex", flexDirection: "column", alignItems: "center",
          justifyContent: "center", minHeight: 200, gap: 8,
        }}>
          <div style={{ fontSize: 9, letterSpacing: 1, color: "#666" }}>AVALANCHE</div>
          {avalanche && (
            <>
              <div style={{
                fontSize: 28, fontWeight: 700,
                color: Math.abs(parseFloat(avalanche.pct) - 50) < 8
                  ? `rgb(${CYAN.join(",")})`
                  : "#ff6666",
              }}>
                {avalanche.pct}%
              </div>
              <div style={{ fontSize: 9, color: "#666" }}>
                {avalanche.flipped}/{avalanche.total} bits
              </div>
              <div style={{
                width: 40, height: 4, borderRadius: 2,
                background: "rgba(255,255,255,0.1)", overflow: "hidden",
              }}>
                <div style={{
                  width: `${avalanche.pct}%`, height: "100%",
                  background: `rgb(${CYAN.join(",")})`, borderRadius: 2,
                }} />
              </div>
              <div style={{ fontSize: 8, color: "#555", textAlign: "center", maxWidth: 60 }}>
                ideal: 50%
              </div>
            </>
          )}
          <div style={{ fontSize: 18, color: "#333" }}>⟷</div>
        </div>

        {/* Avalanche twin */}
        <div style={{ textAlign: "center" }}>
          <div style={{ fontSize: 9, letterSpacing: 2, color: "#886644", marginBottom: 6 }}>
            AVALANCHE TWIN
          </div>
          <canvas ref={canvasB} width={SZ} height={SZ} style={{
            borderRadius: "50%",
            boxShadow: "0 0 30px rgba(255,200,50,0.08)",
            maxWidth: "min(280px, 42vw)", height: "auto",
          }} />
          <div style={{ fontSize: 8, color: "#444", marginTop: 6, wordBreak: "break-all", maxWidth: 260 }}>
            "{modified}"
          </div>
        </div>
      </div>

      {/* Parameter comparison */}
      {paramLabels.length > 0 && (
        <div style={{
          background: "rgba(255,255,255,0.02)", border: "1px solid rgba(255,255,255,0.06)",
          borderRadius: 8, padding: 14, maxWidth: 620, width: "100%",
        }}>
          <div style={{ fontSize: 9, letterSpacing: 2, color: `rgb(${CYAN.join(",")})`, marginBottom: 8 }}>
            PARAMETER COMPARISON
          </div>
          <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: "4px 12px", fontSize: 10 }}>
            <span style={{ color: "#666" }}>Parameter</span>
            <span style={{ color: `rgb(${CYAN.join(",")})` }}>Identity</span>
            <span style={{ color: "#886644" }}>Twin</span>
            {paramLabels.map(([label, a, b], i) => {
              const changed = String(a) !== String(b);
              return [
                <span key={`l${i}`} style={{ color: "#888" }}>{label}</span>,
                <span key={`a${i}`} style={{ color: "#ccc" }}>{a}</span>,
                <span key={`b${i}`} style={{
                  color: changed ? "#ff8866" : "#666",
                  fontWeight: changed ? 700 : 400,
                }}>{b} {changed ? "△" : ""}</span>,
              ];
            })}
          </div>
        </div>
      )}

      {/* Hash hex */}
      <div style={{
        marginTop: 12, maxWidth: 620, width: "100%",
        background: "rgba(255,255,255,0.02)", border: "1px solid rgba(255,255,255,0.06)",
        borderRadius: 8, padding: 14,
      }}>
        <div style={{ fontSize: 9, letterSpacing: 2, color: `rgb(${CYAN.join(",")})`, marginBottom: 6 }}>
          256-BIT HASH DIGEST
        </div>
        <div style={{ fontSize: 9, color: "#888", wordBreak: "break-all", lineHeight: 1.8 }}>
          <span style={{ color: "#666" }}>ID: </span>
          {hexA.match(/.{1,8}/g)?.map((chunk, i) => {
            const chunkB = hexB.substring(i * 8, i * 8 + 8);
            const same = chunk === chunkB;
            return (
              <span key={i} style={{ color: same ? "#555" : `rgb(${CYAN.join(",")})`, marginRight: 3 }}>
                {chunk}
              </span>
            );
          })}
        </div>
        <div style={{ fontSize: 9, color: "#888", wordBreak: "break-all", lineHeight: 1.8, marginTop: 2 }}>
          <span style={{ color: "#666" }}>AV: </span>
          {hexB.match(/.{1,8}/g)?.map((chunk, i) => {
            const chunkA = hexA.substring(i * 8, i * 8 + 8);
            const same = chunk === chunkA;
            return (
              <span key={i} style={{ color: same ? "#555" : "#886644", marginRight: 3 }}>
                {chunk}
              </span>
            );
          })}
        </div>
      </div>

      <div style={{ marginTop: 14, fontSize: 8, color: "#333", letterSpacing: 1, textAlign: "center", maxWidth: 500 }}>
        MTFT — 13-round absorption + 13-round squeeze — structural constants from X₀(143)
        — symmetry range 3–13 (genus) — every wallet gets a unique, reproducible sigil
      </div>
    </div>
  );
}
