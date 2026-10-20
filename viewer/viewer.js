"use strict";
/* CVT mechanism viewer. Data contract: docs/viewer_contract.md. SI units, world frame:
   primary axis at the origin, +X toward the secondary. Axial cuts: z along the shaft, r radial. */

const INK = "#1c1e21", MUTE = "#5a5f66", RED = "#d71921", LINE = "#a9adb3";
const MONO = 'ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace';
// Drawing conventions (not physics): back-plate thickness, shaft radius, spring pocket, weight size.
const DRAW = { plate: 0.012, shaft: 0.0127, rho: 0.0075, Ls0p: 0.060, spLen0s: 0.030, cup: 0.014 };

const VIRIDIS = ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#fde725"];
const LUT = (() => {
  const c = document.createElement("canvas"); c.width = 256; c.height = 1;
  const x = c.getContext("2d"), g = x.createLinearGradient(0, 0, 256, 0);
  VIRIDIS.forEach((col, i) => g.addColorStop(i / (VIRIDIS.length - 1), col));
  x.fillStyle = g; x.fillRect(0, 0, 256, 1);
  const d = x.getImageData(0, 0, 256, 1).data, out = [];
  for (let i = 0; i < 256; i++) out.push(`rgb(${d[4 * i]},${d[4 * i + 1]},${d[4 * i + 2]})`);
  return out;
})();
const cmap = f => LUT[Math.max(0, Math.min(255, Math.round(f * 255)))];

const $ = id => document.getElementById(id);
const S = { D: null, t: 0, playing: false, speed: 1, last: 0 };
let ser = {}, fr = {}, G = {}, plots = [], tRange = [0, 1], tLim = [0, 1000], derived = {};

/* ---------- data ---------- */
function bisect(a, x) {           // largest i with a[i] <= x
  let lo = 0, hi = a.length - 1;
  while (lo < hi) { const m = (lo + hi + 1) >> 1; if (a[m] <= x) lo = m; else hi = m - 1; }
  return lo;
}
function lerpAt(ts, ys, t) {
  if (t <= ts[0]) return ys[0];
  if (t >= ts[ts.length - 1]) return ys[ys.length - 1];
  const i = bisect(ts, t), a = (t - ts[i]) / (ts[i + 1] - ts[i]);
  return ys[i] + a * (ys[i + 1] - ys[i]);
}
const sv = (k, t) => (ser[k] ? lerpAt(ser.t, ser[k], t) : NaN);

function prepare(D) {
  S.D = D; ser = D.series; fr = D.frames; G = D.geometry;
  G.tb = Math.tan(G.beta_deg * Math.PI / 180);
  G.rim_p = G.rp_mech_max + 0.5 * G.belt_thickness + 0.001;
  G.rim_s = G.rs_mech_max + 0.5 * G.belt_thickness + 0.001;
  const r_eff = ser.r_eff_p && ser.r_eff_s;
  const rp = r_eff ? ser.r_eff_p : ser.rstar_p, rs = r_eff ? ser.r_eff_s : ser.rstar_s;
  derived = {
    ratio: rp.map((v, i) => rs[i] / v),
    speed: ser.v ? ser.v.map(v => v * 3.6) : null,
    rpm: ser.w_e.map(w => w * 30 / Math.PI),
    xp: ser.x_p.map(v => v * 1e3), xs: ser.x_s.map(v => v * 1e3),
    sp: ser.slip_p ? ser.slip_p.map(v => v * 100) : null,
    ss: ser.slip_s ? ser.slip_s.map(v => v * 100) : null,
    mp: ser.cap_p ? ser.tq_p.map((v, i) => 100 * (1 - Math.abs(v) / ser.cap_p[i])) : null,
    ms: ser.cap_s ? ser.tq_s_fix.map((v, i) => 100 * (1 - Math.abs(v + ser.tq_s_mov[i]) / ser.cap_s[i])) : null,
  };
  let mx = 0; for (const row of fr.tension) for (const v of row) if (v > mx) mx = v;
  const step = mx > 400 ? 100 : 50;
  G.Tmax = Math.ceil(mx / step) * step;
  tRange = [ser.t[0], ser.t[ser.t.length - 1]];
  const m = D.meta || {};
  $("meta").textContent = `${m.run || "run"}  ${(m.commit || "").slice(0, 7)}  tune:${m.tune || "-"}  dt ${m.dt || "-"}`;
  buildReadouts(); buildPlots(); setTime(tRange[0]);
  $("msg").classList.add("hide");
}

/* ---------- readouts ---------- */
const RO = [
  ["ENGINE", t => sv("w_e", t) * 30 / Math.PI, "rpm", 0],
  ["SPEED", t => sv("v", t) * 3.6, "km/h", 1],
  ["RATIO", t => sv("r_eff_s", t) / sv("r_eff_p", t) || sv("rstar_s", t) / sv("rstar_p", t), "", 3],
  ["F CLAMP P", t => sv("Fa_p", t), "N", 0],
  ["F CLAMP S", t => sv("Fa_s", t), "N", 0],
  ["T TIGHT", t => sv("T_tight", t), "N", 0],
  ["T SLACK", t => sv("T_slack", t), "N", 0],
  ["BELT TEMP", t => sv("T_belt", t), "°C", 1],
];
function buildReadouts() {
  const box = $("readouts"); box.innerHTML = "";
  const items = RO.map(r => ({ l: r[0], f: r[1], u: r[2], d: r[3] }));
  const sh = S.D.shafts || {};
  Object.keys(sh).filter(k => k !== "t").slice(0, 4).forEach(k =>
    items.push({ l: k.replace(/_/g, " ").toUpperCase(), f: t => lerpAt(sh.t, sh[k], t), u: "N", d: 0 }));
  S.ro = items.map(it => {
    const e = document.createElement("div"); e.className = "ro";
    e.innerHTML = `<div class="l">${it.l}</div><div class="v"><span></span><small>${it.u}</small></div>`;
    box.appendChild(e); return { ...it, el: e.querySelector("span") };
  });
}
function updateReadouts(t) {
  for (const r of S.ro) { const v = r.f(t); r.el.textContent = Number.isFinite(v) ? v.toFixed(r.d) : "—"; }
}

/* ---------- canvas helpers ---------- */
function fit(cv) {
  const dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  const c = cv.getContext("2d"); c.setTransform(dpr, 0, 0, dpr, 0, 0); c.clearRect(0, 0, w, h);
  return { c, w, h };
}
function hatch(color, dir) {
  const p = document.createElement("canvas"); p.width = p.height = 8;
  const x = p.getContext("2d"); x.strokeStyle = color; x.lineWidth = 1; x.beginPath();
  if (dir > 0) { x.moveTo(0, 8); x.lineTo(8, 0); x.moveTo(-1, 1); x.lineTo(1, -1); x.moveTo(7, 9); x.lineTo(9, 7); }
  else { x.moveTo(0, 0); x.lineTo(8, 8); x.moveTo(-1, 7); x.lineTo(1, 9); x.moveTo(7, -1); x.lineTo(9, 1); }
  x.stroke(); return x.createPattern(p, "repeat");
}
let HATCH = null;
function text(c, s, x, y, o = {}) {
  c.font = `${o.w || 500} ${o.s || 11}px ${MONO}`; c.fillStyle = o.col || INK;
  c.textAlign = o.al || "left"; c.textBaseline = o.bl || "alphabetic"; c.fillText(s, x, y);
}
function poly(c, pts, V, fill, stroke, lw = 1.2) {
  c.beginPath(); pts.forEach(([a, b], i) => (i ? c.lineTo(V.x(a), V.y(b)) : c.moveTo(V.x(a), V.y(b)))); c.closePath();
  if (fill) { c.fillStyle = fill; c.fill(); }
  if (stroke) { c.strokeStyle = stroke; c.lineWidth = lw; c.lineJoin = "round"; c.stroke(); }
}
function steel(c, V, pts, kind) {
  // base gradient + hatch + outline; kind 0 fixed, 1 movable, 2 shaft/cover
  const ys = pts.map(p => V.y(p[1])), g = c.createLinearGradient(0, Math.min(...ys), 0, Math.max(...ys));
  const cols = [["#aeb4bb", "#8c929a"], ["#c3c8ce", "#9aa1a9"], ["#9ea4ab", "#7d838b"]][kind];
  g.addColorStop(0, cols[0]); g.addColorStop(1, cols[1]);
  poly(c, pts, V, g, null);
  c.save(); c.globalAlpha = 0.5; poly(c, pts, V, HATCH[kind === 1 ? 1 : 0], null); c.restore();
  poly(c, pts, V, null, "#2c3035", 1.4);
}
function coil(c, V, z0, z1, r0, amp, turns, col, lw = 2.4) {
  const n = 160; c.lineJoin = "round"; c.lineCap = "round";
  for (const [w, cc] of [[lw + 2, "rgba(20,22,25,0.55)"], [lw, col]]) {
    c.beginPath();
    for (let i = 0; i <= n; i++) {
      const f = i / n, z = z0 + (z1 - z0) * f, r = r0 + amp * Math.sin(2 * Math.PI * turns * f);
      i ? c.lineTo(V.x(z), V.y(r)) : c.moveTo(V.x(z), V.y(r));
    }
    c.strokeStyle = cc; c.lineWidth = w; c.stroke();
  }
}
function arrow(c, x0, y0, x1, y1, col = RED, lw = 1.5) {
  const a = Math.atan2(y1 - y0, x1 - x0), h = 6;
  c.strokeStyle = col; c.fillStyle = col; c.lineWidth = lw; c.beginPath(); c.moveTo(x0, y0); c.lineTo(x1, y1); c.stroke();
  c.beginPath(); c.moveTo(x1, y1); c.lineTo(x1 - h * Math.cos(a - 0.4), y1 - h * Math.sin(a - 0.4));
  c.lineTo(x1 - h * Math.cos(a + 0.4), y1 - h * Math.sin(a + 0.4)); c.closePath(); c.fill();
}

/* ---------- state at time t ---------- */
function stateAt(t) {
  const ft = fr.t, i = Math.min(ft.length - 1, bisect(ft, t));
  const i2 = Math.min(ft.length - 1, i + 1), a = i2 === i ? 0 : (t - ft[i]) / (ft[i2] - ft[i]);
  const near = a < 0.5 ? i : i2, L = k => fr[k][i] + a * (fr[k][i2] - fr[k][i]);
  const rp = sv("r_eff_p", t), rs = sv("r_eff_s", t);
  return {
    fi: near, xp: L("x_p"), xs: L("x_s"), th_e: L("theta_e"), th_s: L("theta_s"), rw: L("r_w"),
    rp: Math.max(G.rp_min, Number.isFinite(rp) ? rp : sv("rstar_p", t)),
    rs: Number.isFinite(rs) ? rs : sv("rstar_s", t),
    belt: fr.belt[near], ten: fr.tension[near], Tm: sv("T_max", t),
  };
}

/* ---------- side view ---------- */
function drawSide(st) {
  const { c, w, h } = fit($("cv-side"));
  const x0 = -G.rim_p, x1 = G.C + G.rim_s, yh = Math.max(G.rim_p, G.rim_s);
  const padT = 34, padB = 56, pad = 18;
  const sc = Math.min((w - 2 * pad) / (x1 - x0), (h - padT - padB) / (2 * yh));
  const ox = w / 2 - sc * (x0 + x1) / 2, oy = padT + (h - padT - padB) / 2;
  const V = { x: x => ox + sc * x, y: y => oy - sc * y, s: sc };

  c.strokeStyle = "rgba(28,30,33,0.35)"; c.setLineDash([8, 5, 2, 5]); c.lineWidth = 1;
  c.beginPath(); c.moveTo(V.x(x0 - 0.01), V.y(0)); c.lineTo(V.x(x1 + 0.01), V.y(0)); c.stroke(); c.setLineDash([]);

  const sheave = (cx, rim, hub, rp, th) => {
    const X = V.x(cx), Y = V.y(0), g = c.createRadialGradient(X, Y, hub * sc, X, Y, rim * sc);
    g.addColorStop(0, "#e2e4e7"); g.addColorStop(1, "#a3a8af");
    c.beginPath(); c.arc(X, Y, rim * sc, 0, 7); c.fillStyle = g; c.fill();
    c.strokeStyle = "#464b52"; c.lineWidth = 3.5; c.stroke();
    c.beginPath(); c.arc(X, Y, (rim - 0.004) * sc, 0, 7); c.strokeStyle = "rgba(255,255,255,0.45)"; c.lineWidth = 1; c.stroke();
    for (let k = 0; k < 6; k++) {   // rotating spokes (rotation sense: counter-clockwise for +theta)
      const a = th + k * Math.PI / 3, ca = Math.cos(a), sa = -Math.sin(a);
      c.lineCap = "butt";
      for (const [lw, col] of [[7, "#6b7078"], [4.5, "#c8ccd1"]]) {
        c.strokeStyle = col; c.lineWidth = lw; c.beginPath();
        c.moveTo(X + ca * hub * sc, Y + sa * hub * sc); c.lineTo(X + ca * (rim - 0.006) * sc, Y + sa * (rim - 0.006) * sc); c.stroke();
      }
    }
    c.beginPath(); c.arc(X, Y, hub * sc, 0, 7);
    const gh = c.createRadialGradient(X - hub * sc * 0.3, Y - hub * sc * 0.3, 2, X, Y, hub * sc);
    gh.addColorStop(0, "#8a9098"); gh.addColorStop(1, "#3b3f45"); c.fillStyle = gh; c.fill();
    c.strokeStyle = "#1c1e21"; c.lineWidth = 1.5; c.stroke();
    c.beginPath(); c.arc(X, Y, rp * sc, 0, 7); c.setLineDash([3, 3]);
    c.strokeStyle = "rgba(28,30,33,0.7)"; c.lineWidth = 1; c.stroke(); c.setLineDash([]);
    const ka = th;  // keyway tick
    c.strokeStyle = "#e8e9eb"; c.lineWidth = 2; c.beginPath();
    c.moveTo(X + Math.cos(ka) * hub * 0.35 * sc, Y - Math.sin(ka) * hub * 0.35 * sc);
    c.lineTo(X + Math.cos(ka) * hub * 0.9 * sc, Y - Math.sin(ka) * hub * 0.9 * sc); c.stroke();
  };
  sheave(0, G.rim_p, G.rp_min, st.rp, st.th_e);
  sheave(G.C, G.rim_s, G.rs_mech_min, st.rs, st.th_s);

  // belt: dark band, then tension-coloured elements
  const B = st.belt, n = B.length, bw = G.belt_thickness * sc;
  c.lineCap = "butt"; c.lineJoin = "round";
  c.beginPath(); B.forEach((p, i) => (i ? c.lineTo(V.x(p[0]), V.y(p[1])) : c.moveTo(V.x(p[0]), V.y(p[1])))); c.closePath();
  c.strokeStyle = "#16181b"; c.lineWidth = bw + 3; c.stroke();
  for (let k = 0; k < n; k++) {
    const p = B[k], q = B[(k + 1) % n];
    c.beginPath(); c.moveTo(V.x(p[0]), V.y(p[1])); c.lineTo(V.x(q[0]), V.y(q[1]));
    c.strokeStyle = cmap(st.ten[k] / G.Tmax); c.lineWidth = bw; c.stroke();
  }
  c.beginPath(); B.forEach((p, i) => (i ? c.lineTo(V.x(p[0]), V.y(p[1])) : c.moveTo(V.x(p[0]), V.y(p[1])))); c.closePath();
  c.strokeStyle = "rgba(255,255,255,0.28)"; c.lineWidth = 1; c.stroke();

  // flyweights at r_w, 120 deg spacing, rotating with the primary
  for (let k = 0; k < G.n_weights; k++) {
    const a = st.th_e + 2 * Math.PI * k / G.n_weights, X = V.x(st.rw * Math.cos(a)), Y = V.y(st.rw * Math.sin(a));
    c.save(); c.translate(X, Y); c.rotate(-a);
    const wl = 0.016 * sc, ww = 0.012 * sc;
    c.beginPath(); c.roundRect(-wl / 2, -ww / 2, wl, ww, 2); c.fillStyle = RED; c.fill();
    c.strokeStyle = "#fff"; c.lineWidth = 1.5; c.stroke(); c.restore();
  }
  c.beginPath(); c.arc(V.x(0), V.y(0), G.r_w0 * sc, 0, 7); c.setLineDash([2, 4]);
  c.strokeStyle = "rgba(215,25,33,0.55)"; c.lineWidth = 1; c.stroke(); c.setLineDash([]);

  text(c, `PRIMARY  rₚ ${(st.rp * 1e3).toFixed(1)} mm`, V.x(0), V.y(0) + (G.rim_p + 0.012) * sc, { al: "center", w: 700, s: 11 });
  text(c, `SECONDARY  rₛ ${(st.rs * 1e3).toFixed(1)} mm`, V.x(G.C), V.y(0) + (G.rim_s + 0.012) * sc, { al: "center", w: 700, s: 11 });
  text(c, `${(sv("w_e", S.t) * 30 / Math.PI).toFixed(0)} rpm`, V.x(0), V.y(0) - (G.rim_p + 0.007) * sc, { al: "center", col: MUTE });
  text(c, `C = ${(G.C * 1e3).toFixed(0)} mm`, V.x(G.C / 2), V.y(0) - 8, { al: "center", col: MUTE });

  // scale bar and colour bar
  const bx = 16, by = h - 18;
  c.fillStyle = INK; c.fillRect(bx, by, 0.05 * sc, 2); c.fillRect(bx, by - 4, 1, 10); c.fillRect(bx + 0.05 * sc, by - 4, 1, 10);
  text(c, "50 mm", bx, by - 8, { s: 10, col: MUTE });
  const cw = Math.min(260, w * 0.34), cx = w - cw - 18, cy = h - 30;
  const gr = c.createLinearGradient(cx, 0, cx + cw, 0); VIRIDIS.forEach((col, i) => gr.addColorStop(i / 8, col));
  c.fillStyle = gr; c.fillRect(cx, cy, cw, 10); c.strokeStyle = INK; c.lineWidth = 1; c.strokeRect(cx + .5, cy + .5, cw - 1, 9);
  for (let i = 0; i <= 4; i++) {
    const x = cx + cw * i / 4; c.fillStyle = INK; c.fillRect(x, cy + 10, 1, 4);
    text(c, `${Math.round(G.Tmax * i / 4)}`, x, cy + 25, { al: i === 0 ? "left" : i === 4 ? "right" : "center", s: 10, col: MUTE });
  }
  text(c, "BELT TENSION [N]", cx, cy - 5, { s: 10, w: 700 });
}

/* ---------- axial cuts ---------- */
function cutView(cv, zmin, zmax, rmax) {
  const { c, w, h } = fit(cv), pad = 14, top = 30, res = Math.min(150, w * 0.32), aw = w - res;
  const sc = Math.min((aw - 2 * pad) / (zmax - zmin), (h - top - pad) / (2 * rmax));
  const ox = aw / 2 - sc * (zmin + zmax) / 2, oy = top + (h - top - pad) / 2;
  return { c, w, h, rx: aw, V: { x: z => ox + sc * z, y: r => oy - sc * r, s: sc, oy } };
}
function info(c, x, y, rows) {
  for (const [l, v, col] of rows) {
    text(c, l, x, y, { s: 9, col: MUTE }); text(c, v, x, y + 15, { s: 13, w: 700, col: col || INK }); y += 38;
  }
  return y;
}
const mirror = pts => pts.map(([z, r]) => [z, -r]);
const both = (c, V, pts, kind) => { steel(c, V, pts, kind); steel(c, V, mirror(pts), kind); };

function drawBeltSection(c, V, zF, zM, rp, st) {
  const t = G.belt_thickness, r1 = rp - t / 2, r2 = rp + t / 2;
  const pts = [[zF(r1), r1], [zF(r2), r2], [zM(r2), r2], [zM(r1), r1]];
  const col = cmap(Math.min(1, (Number.isFinite(st.Tm) ? st.Tm : G.Tmax) / G.Tmax));
  for (const P of [pts, mirror(pts)]) { poly(c, P, V, col, "#16181b", 1.6); }
}
function drawAxis(c, V, z0, z1) {
  c.strokeStyle = "rgba(28,30,33,0.5)"; c.setLineDash([10, 4, 2, 4]); c.lineWidth = 1; c.beginPath();
  c.moveTo(V.x(z0), V.y(0)); c.lineTo(V.x(z1), V.y(0)); c.stroke(); c.setLineDash([]);
}
function dimV(c, V, z, r0, r1, label, col = RED) {
  const x = V.x(z); arrow(c, x, V.y(r0), x, V.y(r1), col); arrow(c, x, V.y(r1), x, V.y(r0), col);
  text(c, label, x + 5, (V.y(r0) + V.y(r1)) / 2 + 4, { col, w: 700, s: 10 });
}

function drawPrimary(st) {
  const b = G.belt_width, tb = G.tb, ro = G.rim_p, rpm = G.rp_min, P = DRAW.plate, rs = DRAW.shaft;
  const zF = r => -(b / 2 + tb * (r - rpm)), zMx = (r, xp) => b / 2 + tb * (r - rpm) - xp;
  const zbOpen = zMx(ro, -G.gap) + P, zcup = zbOpen + DRAW.Ls0p, zend = zcup + DRAW.cup;
  const zmin = zF(ro) - P - 0.006, zmax = zend + 0.006;
  const { c, w, h, V, rx } = cutView($("cv-pri"), zmin, zmax, ro + 0.012);
  const u = st.xp + G.gap, zM = r => zMx(r, st.xp), zb = zM(ro) + P;
  drawAxis(c, V, zmin, zmax);

  // shaft
  steel(c, V, [[zmin + 0.003, -rs], [zend, -rs], [zend, rs], [zmin + 0.003, rs]], 2);
  // fixed half
  const fix = [[zF(ro) - P, rs], [zF(ro) - P, ro], [zF(ro), ro], [zF(rpm), rpm], [zF(rpm), rs]];
  both(c, V, fix, 0);
  // movable half (+ sleeve on the shaft)
  const mov = [[zM(rpm), rs + 0.003], [zM(rpm), rpm], [zM(ro), ro], [zb, ro], [zb, rs + 0.003]];
  both(c, V, mov, 1);
  both(c, V, [[zM(rpm), rs], [zb + 0.012, rs], [zb + 0.012, rs + 0.003], [zM(rpm), rs + 0.003]], 1);
  // open-position outline of the movable half
  c.save(); c.setLineDash([4, 3]); c.strokeStyle = RED; c.lineWidth = 1;
  for (const s of [1, -1]) {
    c.beginPath(); c.moveTo(V.x(zMx(rpm, -G.gap)), V.y(s * rpm)); c.lineTo(V.x(zMx(ro, -G.gap)), V.y(s * ro));
    c.lineTo(V.x(zbOpen), V.y(s * ro)); c.stroke();
  }
  c.restore();

  drawBeltSection(c, V, zF, zM, st.rp, st);

  // cover with ramp profile; flyweight (circle) on the ramp; spring
  const rho = DRAW.rho, ramp = G.ramp, B = ramp.map(([ui, ri], i) => {
    const j0 = Math.max(0, i - 1), j1 = Math.min(ramp.length - 1, i + 1);
    const kap = (ramp[j1][1] - ramp[j0][1]) / (ramp[j1][0] - ramp[j0][0]), n = Math.hypot(1, kap);
    return [zbOpen + rho - ui + rho * kap / n, ri + rho / n];
  });
  const rcup = 0.032, rtop = B[B.length - 1][1] + 0.01;
  const cover = [[zcup, rs], [zcup, rcup], [B[0][0], rcup], ...B, [B[B.length - 1][0], rtop], [zend, rtop], [zend, rs]];
  both(c, V, cover, 2);
  const zc = zb + rho;
  c.beginPath(); c.arc(V.x(zc), V.y(st.rw), rho * V.s, 0, 7);
  const gw = c.createRadialGradient(V.x(zc) - 3, V.y(st.rw) - 3, 1, V.x(zc), V.y(st.rw), rho * V.s);
  gw.addColorStop(0, "#ff6b70"); gw.addColorStop(1, RED); c.fillStyle = gw; c.fill();
  c.strokeStyle = "#16181b"; c.lineWidth = 1.6; c.stroke();
  const rsp = 0.019;
  for (const s of [1, -1]) coil(c, V, zb + 0.012, zcup, s * rsp, 0.0045, 9, "#d9dde2");

  // annotations
  const Fsp = G.spring_p.F0 + G.spring_p.k * u;
  arrow(c, V.x(zbOpen), V.y(ro + 0.006), V.x(zb), V.y(ro + 0.006));
  c.beginPath(); B.forEach(([z, r], i) => (i ? c.lineTo(V.x(z), V.y(r)) : c.moveTo(V.x(z), V.y(r))));
  c.strokeStyle = RED; c.lineWidth = 2; c.stroke();
  info(c, rx + 8, 52, [["TRAVEL u = x_p + g", `${(u * 1e3).toFixed(2)} mm`, RED],
    ["SPRING FORCE", `${Fsp.toFixed(0)} N`], ["FLYWEIGHT RADIUS", `${(st.rw * 1e3).toFixed(1)} mm`, RED],
    ["BELT RADIUS", `${(st.rp * 1e3).toFixed(1)} mm`], ["SPRING RATE", `${(G.spring_p.k / 1e3).toFixed(2)} N/mm`]]);
  text(c, "fixed", V.x(zF(ro) - P / 2), V.y(ro) - 5, { al: "center", s: 10, col: MUTE });
  text(c, "movable", V.x(zb - P / 2), V.y(-ro) + 12, { al: "center", s: 10, col: MUTE });
  text(c, "cover + ramp", V.x(zend - 0.005), V.y(rtop) - 5, { al: "center", s: 10, col: MUTE });
  text(c, "cut through flyweight 1", V.x(zend), V.y(-rtop) + 14, { al: "right", s: 9, col: MUTE });
}

function drawSecondary(st) {
  const b = G.belt_width, tb = G.tb, ro = G.rim_s, rhub = G.rs_mech_min, P = DRAW.plate, rs = DRAW.shaft, rc = G.r_cam;
  const zF = r => -(b / 2 + tb * (r - G.rs_mech_max)), zMx = (r, xs) => b / 2 + tb * (r - G.rs_mech_max) + xs;
  const zb0 = zMx(ro, 0) + P, zcup = zb0 + G.s_s + DRAW.spLen0s, zend = zcup + DRAW.cup;
  const zmin = zF(ro) - P - 0.006, zmax = zend + 0.006;
  const { c, w, h, V, rx } = cutView($("cv-sec"), zmin, zmax, ro + 0.012);
  const zM = r => zMx(r, st.xs), zb = zM(ro) + P;
  drawAxis(c, V, zmin, zmax);

  steel(c, V, [[zmin + 0.003, -rs], [zend, -rs], [zend, rs], [zmin + 0.003, rs]], 2);
  both(c, V, [[zF(ro) - P, rs], [zF(ro) - P, ro], [zF(ro), ro], [zF(rhub), rhub], [zF(rhub), rs]], 0);
  both(c, V, [[zM(rhub), rs + 0.003], [zM(rhub), rhub], [zM(ro), ro], [zb, ro], [zb, rs + 0.003]], 1);
  both(c, V, [[zM(rhub), rs], [zb + 0.012, rs], [zb + 0.012, rs + 0.003], [zM(rhub), rs + 0.003]], 1);
  // closed-position outline
  c.save(); c.setLineDash([4, 3]); c.strokeStyle = RED; c.lineWidth = 1;
  for (const s of [1, -1]) {
    c.beginPath(); c.moveTo(V.x(zMx(rhub, 0)), V.y(s * rhub)); c.lineTo(V.x(zMx(ro, 0)), V.y(s * ro)); c.lineTo(V.x(zb0), V.y(s * ro)); c.stroke();
  }
  c.restore();
  drawBeltSection(c, V, zF, zM, st.rs, st);

  // cam sleeve (on the movable half) with slot at the roller, spring cup (fixed to shaft)
  const tw = 0.003, zroll = zcup - 0.020;
  for (const s of [1, -1]) {
    const sl = [[zb, rc - tw], [zroll - 0.007, rc - tw], [zroll - 0.007, rc + tw], [zb, rc + tw]];
    steel(c, V, s > 0 ? sl : mirror(sl), 1);
  }
  both(c, V, [[zcup, rs], [zcup, rc + 0.012], [zend, rc + 0.012], [zend, rs]], 2);
  // roller on its pin at r_cam, fixed to the cup
  for (const s of [1, -1]) {
    c.beginPath(); c.arc(V.x(zroll), V.y(s * rc), 0.0055 * V.s, 0, 7); c.fillStyle = RED; c.fill(); c.strokeStyle = "#16181b"; c.lineWidth = 1.6; c.stroke();
    steel(c, V, [[zroll, s * rc - 0.0015], [zcup, s * rc - 0.0015], [zcup, s * rc + 0.0015], [zroll, s * rc + 0.0015]], 2);
  }
  // linear spring and torsion spring (outside the sleeve)
  for (const s of [1, -1]) {
    coil(c, V, zb + 0.012, zcup, s * 0.0235, 0.0045, 8, "#d9dde2");
    coil(c, V, zb + 0.002, zcup - 0.002, s * (rc + 0.0085), 0.0022, 14, "#7a9bc0", 2);
  }

  // annotations
  const hx = interp(G.helix, st.xs), th = G.spring_s.theta_t0_deg + hx * 180 / Math.PI;
  const Fl = G.spring_s.F0 + G.spring_s.k * st.xs;
  arrow(c, V.x(zb0), V.y(ro + 0.006), V.x(zb), V.y(ro + 0.006));
  info(c, rx + 8, 52, [["TRAVEL x_s", `${(st.xs * 1e3).toFixed(2)} mm`, RED],
    ["LINEAR SPRING", `${Fl.toFixed(0)} N`], ["TORSION ANGLE", `${th.toFixed(1)}°`, "#3f6592"],
    ["BELT RADIUS", `${(st.rs * 1e3).toFixed(1)} mm`]]);
  text(c, "fixed", V.x(zF(ro) - P / 2), V.y(ro) - 5, { al: "center", s: 10, col: MUTE });
  text(c, "movable + helix sleeve", V.x(zb), V.y(-ro) + 12, { al: "center", s: 10, col: MUTE });

  // inset: unrolled helix, circumferential s = r_cam h against axial x (equal aspect)
  const sz = Math.min(110, w - rx - 24), ix = rx + 10, iy = h - sz - 30, hmax = G.helix[G.helix.length - 1][1] * rc;
  const sc = sz / Math.max(G.s_s, hmax), X = x => ix + 4 + sc * x, Y = s => iy + sz - 4 - sc * s;
  c.fillStyle = "rgba(232,233,235,0.7)"; c.fillRect(ix - 6, iy - 18, sz + 12, sz + 36);
  c.strokeStyle = LINE; c.strokeRect(ix - 5.5, iy - 17.5, sz + 11, sz + 35);
  text(c, "HELIX  h(x)", ix, iy - 5, { s: 9, w: 700 });
  c.beginPath(); G.helix.forEach(([x, hh], i) => (i ? c.lineTo(X(x), Y(hh * rc)) : c.moveTo(X(x), Y(hh * rc))));
  c.strokeStyle = INK; c.lineWidth = 2; c.stroke();
  c.setLineDash([2, 3]); c.strokeStyle = MUTE; c.lineWidth = 1; c.beginPath();
  c.moveTo(X(st.xs), Y(0)); c.lineTo(X(st.xs), Y(hx * rc)); c.lineTo(X(0), Y(hx * rc)); c.stroke(); c.setLineDash([]);
  c.beginPath(); c.arc(X(st.xs), Y(hx * rc), 4.5, 0, 7); c.fillStyle = RED; c.fill(); c.strokeStyle = "#fff"; c.lineWidth = 1.5; c.stroke();
  text(c, "x (axial)", ix + sz, iy + sz + 12, { al: "right", s: 9, col: MUTE });
  text(c, "r·h", ix, iy + sz + 12, { s: 9, col: MUTE });
}
function interp(tab, x) {
  if (x <= tab[0][0]) return tab[0][1];
  for (let i = 1; i < tab.length; i++) if (x <= tab[i][0]) {
    const a = (x - tab[i - 1][0]) / (tab[i][0] - tab[i - 1][0]); return tab[i - 1][1] + a * (tab[i][1] - tab[i - 1][1]);
  }
  return tab[tab.length - 1][1];
}

/* ---------- plots ---------- */
const PCOL = { a: "#23262b", b: "#2f6aa8", c: RED };
function buildPlots() {
  const host = $("plots"); host.innerHTML = ""; plots.forEach(p => p.u.destroy()); plots = [];
  const T = ser.t, d = derived;
  const defs = [
    ["ENGINE RPM vs SPEED", "km/h", d.speed || T, [["rpm", d.rpm, PCOL.a]], "rpm"],
    ["RATIO", "s", T, [["r_s / r_p", d.ratio, PCOL.a]], ""],
    ["SHEAVE TRAVEL", "s", T, [["x_p [mm]", d.xp, PCOL.a], ["x_s [mm]", d.xs, PCOL.b]], "mm"],
    ["BELT TENSION", "s", T, [["T tight", ser.T_tight, PCOL.a], ["T slack", ser.T_slack, PCOL.b], ["T max", ser.T_max, PCOL.c]], "N"],
    ["SLIP", "s", T, [["primary", d.sp, PCOL.a], ["secondary", d.ss, PCOL.b]], "%"],
    ["TRACTION MARGIN", "s", T, [["primary", d.mp, PCOL.a], ["secondary", d.ms, PCOL.b]], "% of capacity"],
  ];
  for (const [title, xl, xs, lines, yl] of defs) {
    const L = lines.filter(l => l[1]); if (!L.length) continue;
    const el = document.createElement("div"); el.className = "plot";
    el.innerHTML = `<h3>${title}<i>${yl}</i></h3><div></div>`; host.appendChild(el);
    const ax = { stroke: INK, grid: { stroke: "rgba(28,30,33,0.12)", width: 1 }, ticks: { stroke: "rgba(28,30,33,0.4)", width: 1 }, font: `10px ${MONO}`, size: 40 };
    const opts = {
      width: el.clientWidth - 16, height: 200,
      scales: { x: { time: false } }, legend: { live: false },
      cursor: { drag: { x: false, y: false, setScale: false }, points: { show: false } },
      axes: [{ ...ax, size: 28, label: xl, labelFont: `10px ${MONO}`, labelSize: 14 }, { ...ax, size: 52 }],
      series: [{}, ...L.map(l => ({ label: l[0], stroke: l[2], width: 1.6 }))],
      hooks: { draw: [u => cursorMark(u, xs)], ready: [u => seekHooks(u)] },
    };
    const u = new uPlot(opts, [xs, ...L.map(l => l[1])], el.children[1]);
    plots.push({ u, xs, el });
  }
}
function cursorMark(u, xs) {
  if (!S.D) return;
  const c = u.ctx, t = S.t, xv = xs === ser.t ? t : lerpAt(ser.t, xs, t), x = u.valToPos(xv, "x", true);
  c.save(); c.strokeStyle = RED; c.lineWidth = 1.5;
  if (xs === ser.t) { c.beginPath(); c.moveTo(x, u.bbox.top); c.lineTo(x, u.bbox.top + u.bbox.height); c.stroke(); }
  for (let i = 1; i < u.series.length; i++) {
    const y = u.valToPos(lerpAt(ser.t, u.data[i], t), "y", true);
    c.beginPath(); c.arc(x, y, 4 * devicePixelRatio, 0, 7); c.fillStyle = RED; c.fill(); c.strokeStyle = "#fff"; c.lineWidth = 1.2 * devicePixelRatio; c.stroke();
  }
  c.restore();
}
function seekHooks(u) {
  const go = e => { if (u.cursor.idx != null) setTime(ser.t[u.cursor.idx]); };
  u.over.addEventListener("click", go);
  u.over.addEventListener("pointermove", e => { if (e.buttons & 1) go(e); });
}

/* ---------- transport ---------- */
let raf = 0;
function setTime(t) {
  S.t = Math.max(tRange[0], Math.min(tRange[1], t));
  const st = stateAt(S.t);
  $("tread").textContent = `t ${S.t.toFixed(4)} s  f ${st.fi + 1}/${fr.t.length}`;
  $("scrub").value = Math.round(10000 * (S.t - tRange[0]) / (tRange[1] - tRange[0]));
  updateReadouts(S.t);
  if (!raf) raf = requestAnimationFrame(() => { raf = 0; render(); });
}
function render() {
  if (!S.D) return;
  const st = stateAt(S.t);
  drawSide(st); drawPrimary(st); drawSecondary(st);
  plots.forEach(p => p.u.redraw(false, false));
}
function tick(now) {
  if (!S.playing) return;
  const dt = (now - S.last) / 1000; S.last = now;
  let t = S.t + dt * S.speed;
  if (t >= tRange[1]) { t = tRange[1]; setPlaying(false); }
  setTime(t); requestAnimationFrame(tick);
}
function setPlaying(on) {
  S.playing = on; $("btn-play").textContent = on ? "PAUSE" : "PLAY"; $("btn-play").classList.toggle("on", on);
  if (on) { if (S.t >= tRange[1]) S.t = tRange[0]; S.last = performance.now(); requestAnimationFrame(tick); }
}
function step(n) { const ft = fr.t, i = Math.max(0, Math.min(ft.length - 1, stateAt(S.t).fi + n)); setPlaying(false); setTime(ft[i]); }

/* ---------- loading ---------- */
function showMsg(s, err) { const m = $("msg"); m.textContent = s; m.className = "msg" + (err ? " err" : ""); }
async function loadUrl(u) {
  try { const r = await fetch(u); if (!r.ok) throw new Error(`${u}: HTTP ${r.status}`); prepare(await r.json()); }
  catch (e) { showMsg("Load failed: " + e.message, true); }
}
function init() {
  HATCH = [hatch("rgba(20,22,25,0.75)", 1), hatch("rgba(20,22,25,0.75)", -1)];
  $("btn-play").onclick = () => setPlaying(!S.playing);
  $("btn-prev").onclick = () => step(-1); $("btn-next").onclick = () => step(1);
  $("speed").onchange = e => { S.speed = +e.target.value; };
  $("scrub").oninput = e => { setTime(tRange[0] + (tRange[1] - tRange[0]) * e.target.value / 10000); };
  addEventListener("keydown", e => {
    if (e.target.tagName === "SELECT") return;
    if (e.code === "Space") { e.preventDefault(); setPlaying(!S.playing); }
    else if (e.code === "ArrowRight") step(e.shiftKey ? 10 : 1);
    else if (e.code === "ArrowLeft") step(e.shiftKey ? -10 : -1);
  });
  addEventListener("resize", () => { plots.forEach(p => p.u.setSize({ width: p.el.clientWidth - 16, height: 200 })); render(); });
  addEventListener("dragover", e => { e.preventDefault(); document.body.classList.add("drag"); });
  addEventListener("dragleave", () => document.body.classList.remove("drag"));
  addEventListener("drop", e => {
    e.preventDefault(); document.body.classList.remove("drag");
    const f = e.dataTransfer.files[0]; if (!f) return;
    f.text().then(s => prepare(JSON.parse(s))).catch(err => showMsg("Load failed: " + err.message, true));
  });
  const run = new URLSearchParams(location.search).get("run");
  if (run) loadUrl(run);
}
document.fonts && document.fonts.load('700 20px "Doto"').catch(() => {});
init();
