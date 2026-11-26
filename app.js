/* ============================================================
   app.js — 交互、Canvas 渲染、图例、反应面板
   依赖 data.js / chem.js
   ============================================================ */
(function () {
  "use strict";
  const D = PBData, C = PBChem, byId = D.byId;
  const $ = id => document.getElementById(id);
  const $$ = sel => Array.prototype.slice.call(document.querySelectorAll(sel));
  window.addEventListener("error", function (e) {
    const w = document.getElementById("warnbox");
    if (w) w.innerHTML = '<div style="color:#b84a4a">运行错误：' + (e.message || e.error) + ' @ ' + (e.filename || "") + ':' + (e.lineno || "") + '</div>';
  });
  const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
  const PI = Math.PI;

  /* ================= 通用格式 ================= */
  const SUBD = "₀₁₂₃₄₅₆₇₈₉";
  const SUPD = "⁰¹²³⁴⁵⁶⁷⁸⁹";
  function uniSub(n) { return String(n).split("").map(c => SUBD[+c] || c).join(""); }
  function uniSup(s) {
    return String(s).split("").map(c => {
      if (c === "+") return "⁺";
      if (c === "-") return "⁻";
      if (c === "−") return "⁻";
      if (SUPD[+c] !== undefined && /[0-9]/.test(c)) return SUPD[+c];
      return c;
    }).join("");
  }
  function fmtConc(c) {
    if (!(c > 0)) return "—";
    if (c >= 0.01 && c < 100) {
      let s = c >= 1 ? c.toPrecision(3) : c.toPrecision(3);
      s = s.replace(/(\.\d*?)0+$/, "$1").replace(/\.$/, "");
      return s + " M";
    }
    const e = Math.floor(Math.log10(c));
    const m = c / Math.pow(10, e);
    return (Math.abs(m - 1) < 0.005 ? "1.0" : m.toFixed(2)) + "×10" + uniSup(String(e)) + " M";
  }
  function fmtP(p) {
    if (!(p > 0)) return "—";
    if (p >= 0.01 && p < 100) return p.toPrecision(3).replace(/\.?0+$/, "") + " bar";
    const e = Math.floor(Math.log10(p));
    const m = p / Math.pow(10, e);
    return (Math.abs(m - 1) < 0.005 ? "1.0" : m.toFixed(2)) + "×10" + uniSup(String(e)) + " bar";
  }
  function fmtNum(v, d) { return (Math.abs(v) < 5e-7 ? "0" : v.toFixed(d == null ? 3 : d)).replace("-", "−"); }
  function fmtE(v) { return fmtNum(v, 3) + " V"; }
  function fmtG(v) { return fmtNum(v, 1) + " kJ·mol⁻¹"; }
  function fmtLogK(v) {
    if (!isFinite(v)) return "—";
    if (Math.abs(v) >= 100) {
      const e = Math.floor(Math.log10(Math.abs(v)));
      return (v < 0 ? "−" : "") + (Math.abs(v) / Math.pow(10, e)).toFixed(2) + "×10" + uniSup(String(e));
    }
    return fmtNum(v, 2);
  }
  function oxLabel(ox) {
    if (ox == null) return "—";
    if (Math.abs(ox - Math.round(ox)) < 1e-9) {
      const n = Math.round(ox);
      return n === 0 ? "0" : n > 0 ? "+" + n : "−" + Math.abs(n);
    }
    if (Math.abs(ox - 8 / 3) < 1e-9) return "8/3";
    return String(Math.round(ox * 100) / 100).replace("-", "−");
  }
  function stateName(st) { return st === "aq" ? "aq" : st === "s" ? "s" : st === "g" ? "g" : st === "l" ? "l" : ""; }

  /* ================= 化学式解析 / 渲染 ================= */
  const pfCache = Object.create(null);
  function parseFormula(s) {
    if (pfCache[s]) return pfCache[s];
    let i = 0;
    function num() {
      let n = "";
      while (i < s.length && /[0-9]/.test(s[i])) n += s[i++];
      return n ? parseInt(n, 10) : 1;
    }
    function seq(stopAtParen) {
      const out = [];
      while (i < s.length) {
        const c = s[i];
        if (c === ")" && stopAtParen) { i++; break; }
        if (c === "(") { i++; const inner = seq(true); const n = num(); out.push({ t: "g", inner, n }); continue; }
        if (c === "^") {
          i++; let q = "";
          while (i < s.length && /[0-9+\-−]/.test(s[i])) q += s[i++];
          out.push({ t: "q", q: q.replace("−", "-") });
          continue;
        }
        if (c === "·") {
          i++; let m = "";
          if (s[i] === "x") { m = "x"; i++; }
          else { while (i < s.length && /[0-9.]/.test(s[i])) m += s[i++]; }
          out.push({ t: "hyd", m: m || "1" });
          continue;
        }
        if (/[A-Z]/.test(c)) {
          let sym = c; i++;
          if (i < s.length && /[a-z]/.test(s[i])) sym += s[i++];
          const n = num();
          out.push({ t: "e", sym, n });
          continue;
        }
        i++;
      }
      return out;
    }
    const tree = seq(false);
    const q = tree.find(t => t.t === "q");
    const body = tree.filter(t => t.t !== "q");
    const out = { body, charge: q ? q.q : null };
    pfCache[s] = out;
    return out;
  }
  function bodyHasPoly(body) {
    return body.length > 1 || body.some(t => t.t === "g");
  }
  function renderTree(items, mode, plainStyle) {
    let out = "";
    const sub = n => mode === "html" ? "<sub>" + n + "</sub>" : mode === "unicode" ? uniSub(n) : String(n);
    for (const it of items) {
      if (it.t === "e") out += it.sym + (it.n > 1 ? sub(it.n) : "");
      else if (it.t === "g") out += "(" + renderTree(it.inner, mode, plainStyle) + ")" + (it.n > 1 ? sub(it.n) : "");
      else if (it.t === "hyd") out += "·" + it.m + "H" + sub(2) + "O";
    }
    return out;
  }
  function renderFormula(pf, mode, plainStyle) {
    mode = mode || state.fmode;
    plainStyle = plainStyle || state.fplain;
    if (pf === "e^-" || pf === "e-") {
      if (mode === "html") return "e<sup>−</sup>";
      if (mode === "unicode") return "e⁻";
      return "e-";
    }
    const tree = parseFormula(pf);
    const body = renderTree(tree.body, mode, plainStyle);
    if (!tree.charge) return body;
    const q = tree.charge.replace("-", "−");
    if (mode === "html") return body + "<sup>" + q + "</sup>";
    if (mode === "unicode") return body + uniSup(tree.charge);
    if (plainStyle === "spaced") return body + " " + q;
    if (plainStyle === "paren" && bodyHasPoly(tree.body)) return "(" + body + ")" + q;
    return body + q;
  }
  const chemHTML = pf => renderFormula(pf, "html");
  const chemUni = pf => renderFormula(pf, "unicode");
  const chemPlain = pf => renderFormula(pf, "plain", state.fplain);
  function formulaFor(sp, mode) { return renderFormula(sp.formula, mode || state.fmode); }
  function formulaById(id, mode) {
    if (id === "e-") return mode === "html" ? "e<sup>−</sup>" : mode === "unicode" ? "e⁻" : "e-";
    if (id.startsWith("@")) {
      const sym = id.replace(/^@el:/, "").replace(/^@/, "");
      const f = renderFormula(sym, mode || "html");
      return mode === "html" ? f + "<sub>(s)</sub>" : f + "(s)";
    }
    const sp = byId[id];
    return sp ? renderFormula(sp.formula, mode || "html") : id;
  }
  function formulaByIdUni(id) { return formulaById(id, "unicode"); }

  /* ================= 状态 ================= */
  const state = {
    main: "Fe", selected: new Set(), customC: {}, customP: {},
    env: {}, filters: { ox: new Set(), cat: new Set() },
    Ctot: 1e-3, Penv: 1, V0: 0, couple: true,
    mode: "C", Cdisp: 1e-5, hoverPx: 8,
    xAxis: "pH", pHfix: 7,
    xmin: 0, xmax: 14, emin: -1.2, emax: 1.6,
    palMode: "intrinsic", palSeed: 0, customColors: {}, randomMap: {},
    fmode: "html", fplain: "compact", flegend: "unicode",
    showBoundary: true, showWater: true, showMedium: true, showAllRedox: false,
    showLabels: true, showGrid: true,
    sys: null, sysKey: null, polygons: [], boundaries: [], boundaryChains: [], labels: [], groups: null,
    field: null, fieldKey: null,
    hover: null, activePreset: null,
    renderMs: 0,
  };

  const THEMES = {
    quiet: {
      plotTop: "#fdfcff", plotBottom: "#faf7fd", plot: "#fdfcff",
      grid: "#e9e1f3", gridMinor: "#f4eef9", axis: "#6d5a88", text: "#4a4158",
      muted: "#8b7f9c", line: "#d9cde9", water: "#b3a0c8", medium: "#a58fc4",
      boundary: "#6f5d88", labelBg: "rgba(255,255,255,0.86)", paper: "#ffffff",
      ink: "#342c40", accent: "#6f4fa0", dot: "#ffffff",
      bgRgb: [253, 252, 255], light: true,
    },
  };
  const theme = () => THEMES.quiet;

  /* ================= 颜色 ================= */
  function hash(s) {
    let h = 2166136261;
    for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); }
    return h >>> 0;
  }
  function hslRgb(h, s, l) {
    h = ((h % 360) + 360) % 360 / 360;
    const f = n => {
      const k = (n + h * 12) % 12;
      const a = s * Math.min(l, 1 - l);
      return Math.round(255 * (l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1))));
    };
    return [f(0), f(8), f(4)];
  }
  function hexRgb(c) {
    if (!c) return [150, 140, 170];
    if (Array.isArray(c)) return c;
    const m = String(c).match(/rgba?\((\d+),\s*(\d+),\s*(\d+)/);
    if (m) return [+m[1], +m[2], +m[3]];
    let s = String(c).replace("#", "");
    if (s.length === 3) s = s.split("").map(x => x + x).join("");
    return [parseInt(s.slice(0, 2), 16), parseInt(s.slice(2, 4), 16), parseInt(s.slice(4, 6), 16)];
  }
  function rgbStr(a) { return "rgb(" + a[0] + "," + a[1] + "," + a[2] + ")"; }
  function toHex(c) {
    if (Array.isArray(c)) return "#" + c.map(v => clamp(Math.round(v), 0, 255).toString(16).padStart(2, "0")).join("");
    if (String(c)[0] === "#") return String(c).slice(0, 7);
    const m = String(c).match(/rgb\((\d+),(\d+),(\d+)\)/);
    return m ? toHex([+m[1], +m[2], +m[3]]) : "#8b7f9c";
  }
  function mixRgb(a, b, t) { return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]; }
  function lightenRgb(a, t) { return mixRgb(a, [255, 255, 255], t); }
  function darkenRgb(a, t) { return mixRgb(a, [10, 8, 16], t); }

  // 无色/近白物种没有“真实颜色”，用代表性淡彩区分；有真实颜色的离子不在此表。
  const COLOR_OVERRIDES = {
    "Fe": "#9a9ca0", "Mg": "#b8bcc0", "Ag": "#c8c8d0", "Pb": "#8a8c94",
    "Zn": "#b0b4b8", "Zn2+": "#c8ccd0", "ZnO": "#d8d8d0", "Zn(OH)2": "#f2f2ee",
    "Zn(OH)42-": "#c8d8e0", "ZnCO3": "#e8e8d8", "ZnS": "#d8d0a0",
    "Zn(NH3)42+": "#c8d0e8", "Zn(CN)42-": "#d8d0e0", "ZnCl42-": "#d0d8d8",
    "V": "#b8bcc0", "Cr": "#b0b4b8", "Mn": "#a8acb0",
    "Al": "#b8bcc0", "Al3+": "#7fa8d8", "Al(OH)3": "#d8b878", "AlO(OH)": "#c8a868",
    "Al2O3": "#b8a0d8", "Al(OH)4-": "#78b898", "AlF63-": "#98b8b0",
    "Cl-": "#7fa8c8", "Cl2": "#b8cc60", "Cl2g": "#a8c050", "HClO": "#e0b878",
    "ClO-": "#e8a888", "HCl": "#a8a8c0", "HClg": "#b8b8c8", "AgCl": "#e8e8e0",
    "NH4+": "#7fa8d8", "NH3": "#78c8a8", "NH3g": "#98d0b8", "N2": "#b8b8a8",
    "HNO2": "#d88898", "NO2-": "#e8a868", "NO": "#a8a8b8", "NO2": "#c88868", "NO3-": "#a8b8c8",
    "H2S": "#88a8c8", "H2Sg": "#98b8d8", "HS-": "#78b8b8", "S2-": "#88c888",
    "SO2": "#e8b878", "SO2g": "#e8c888", "HSO3-": "#e8a8a8", "SO32-": "#e8c888",
    "S2O32-": "#98c8b8", "S4O62-": "#b8a8d8", "S2O42-": "#88b8d8",
    "S2O52-": "#b8c8a8", "S2O72-": "#a8b8c8", "S2O82-": "#b8b8d8",
    "HSO4-": "#a8b8c8", "SO42-": "#98a8b8",
    "H3PO4": "#e8c898", "H2PO4-": "#e0d0a0", "HPO42-": "#c8d8a8", "PO43-": "#a8d8b8",
    "H3PO3": "#e8b8a8", "H2PO3-": "#e0c8a8", "HPO32-": "#c8d8b8",
    "H3PO2": "#e8d0a8", "H2PO2-": "#d8d8a8", "PH3": "#b8b8c8",
    "CO2": "#b8c8d8", "CO2g": "#c8d8e8", "H2CO3": "#c8d8e8", "HCO3-": "#a8c8d8",
    "CO32-": "#98b8d0", "CO": "#b8b8b8", "HCOOH": "#d8b8a8", "HCOO-": "#c8a898",
    "CH4": "#d0d0d0", "CH3OH": "#d8d0c0", "CH3OHg": "#d8d8c8",
    "Fe2+": "#78c898", "Fe(OH)2": "#b8d8c0", "Fe(OH)3": "#b87858",
    "FeO(OH)": "#a87858", "Fe3O4": "#5a5a62",
    "Mn2+": "#f0b8d0", "Mn(OH)2": "#d8d8c8", "MnO": "#8a9a78",
    "MnO(OH)": "#9a7a5a", "Mn2O3": "#6a4a3a", "MnO2": "#3a3a3a",
    "MnO42-": "#58b878", "MnO4-": "#a04aa0",
    "Mg2+": "#9fb8d0", "Mg(OH)2": "#d8e0e8", "MgO": "#c8ccd8",
    "Ca2+": "#b8c8d8", "Ba2+": "#b8c8d8", "Ag+": "#c8c8d8",
    "Pb2+": "#9fb0c0", "Pb(OH)2": "#d8dce4", "PbO": "#d8b858", "PbO2": "#8a6a50",
    "I-": "#b8a0c8", "I3-": "#a07040", "I2": "#6a4a7a", "IO3-": "#e0d0a0", "AgI": "#e8d870",
    "H2O2": "#a8cce8", "OH-": "#9fc0e8", "O2": "#a8c8f0",
  };
  const BOUNDARY_COLORS = {
    redox: "#7189bd", acidbase: "#9b7fc0", precipitation: "#ab8b6d",
    complex: "#74a08c", phase: "#9c92aa", stable: "#6f5d88",
  };
  function baseColorFor(id) {
    const sp = byId[id];
    if (state.customColors[id]) return state.customColors[id];
    if (state.palMode === "intrinsic") {
      if (COLOR_OVERRIDES[id]) return COLOR_OVERRIDES[id];
      if (sp && sp.color) {
        const rgb = hexRgb(sp.color);
        const mx = Math.max(...rgb), mn = Math.min(...rgb);
        const sat = mx ? (mx - mn) / mx : 0;
        const lum = (rgb[0] * 0.299 + rgb[1] * 0.587 + rgb[2] * 0.114) / 255;
        if (sat >= 0.18 || lum < 0.55) return sp.color;
      }
    }
    if (state.palMode === "random") {
      const h = state.randomMap[id] != null ? state.randomMap[id] : hash(id + "#" + state.palSeed) % 360;
      return rgbStr(hslRgb(h, 0.38, 0.72));
    }
    const h = (hash(id) + (sp && sp.main ? sp.main.charCodeAt(0) * 7 : 0)) % 360;
    return rgbStr(hslRgb(h, 0.38, 0.72));
  }
  function colorFor(id) {
    const th = theme();
    const base = hexRgb(baseColorFor(id));
    let fill, line, text;
    const lum = (base[0] * 0.299 + base[1] * 0.587 + base[2] * 0.114) / 255;
    if (th.light) {
      fill = lightenRgb(base, 0.70);
      line = darkenRgb(base, 0.34);
      text = darkenRgb(base, 0.52);
    } else {
      fill = mixRgb(base, th.bgRgb, 0.52);
      line = lightenRgb(base, 0.38);
      text = lightenRgb(base, 0.60);
    }
    // 深色/黑色相区稍微提亮，避免糊成一片
    if (lum < 0.22) fill = mixRgb(fill, th.light ? [210, 206, 216] : [72, 66, 84], th.light ? 0.45 : 0.35);
    // 白色沉淀（AgCl、Fe(OH)₂ 等）加一点灰绿底，避免与背景融为一体
    if (lum > 0.86) {
      fill = mixRgb(fill, th.light ? [214, 220, 216] : [60, 58, 66], th.light ? 0.34 : 0.30);
      line = mixRgb(line, th.light ? [132, 138, 134] : [170, 178, 174], 0.45);
    }
    return { fill: rgbStr(fill), rgbFill: fill, line: rgbStr(line), text: rgbStr(text), trueColor: toHex(base) };
  }
  function recolorRandom() {
    const ids = state.sys ? state.sys.selectedMain.slice() : [];
    if (!ids.length) { state.randomMap = {}; return; }
    const adj = new Map(ids.map(id => [id, new Set()]));
    for (const s of state.boundaries) {
      if (adj.has(s.a) && adj.has(s.b) && s.a !== s.b) { adj.get(s.a).add(s.b); adj.get(s.b).add(s.a); }
    }
    const order = ids.slice().sort((a, b) => adj.get(b).size - adj.get(a).size);
    const seed = state.palSeed || 0;
    const hues = [];
    for (let i = 0; i < 12; i++) hues.push((i * 30 + seed * 37) % 360);
    hues.sort((a, b) => hash("h" + a + "#" + seed) - hash("h" + b + "#" + seed));
    const map = {};
    for (const id of order) {
      const used = [...adj.get(id)].map(n => map[n]).filter(v => v != null);
      let best = hues[0], bestScore = -1;
      for (const h of hues) {
        let mind = 360;
        for (const u of used) { let d = Math.abs(h - u); d = Math.min(d, 360 - d); if (d < mind) mind = d; }
        const score = used.length ? mind : 360;
        if (score > bestScore) { bestScore = score; best = h; }
      }
      map[id] = best;
    }
    state.randomMap = map;
  }

  /* ================= 主题 / 画布 ================= */
  function applyTheme() { document.body.dataset.theme = "quiet"; }
  const plot = $("plot"), overlay = $("overlay"), ctx = plot.getContext("2d"), octx = overlay.getContext("2d");
  let geom = null;
  const PAD = { left: 64, right: 24, top: 18, bottom: 48 };

  function setupCanvas() {
    const wrap = $("canvas-wrap");
    const cssW = Math.max(320, Math.round(wrap.clientWidth || 900));
    const cssH = Math.round(clamp(cssW * 0.615, 420, 720));
    const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
    plot.style.height = cssH + "px";
    plot.width = Math.round(cssW * dpr);
    plot.height = Math.round(cssH * dpr);
    overlay.style.height = cssH + "px";
    overlay.width = Math.round(cssW * dpr);
    overlay.height = Math.round(cssH * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    octx.setTransform(dpr, 0, 0, dpr, 0, 0);
    geom = {
      W: cssW, H: cssH, dpr,
      plotX: PAD.left, plotY: PAD.top,
      plotW: cssW - PAD.left - PAD.right,
      plotH: cssH - PAD.top - PAD.bottom,
    };
  }
  const X = x => geom.plotX + (x - state.xmin) / (state.xmax - state.xmin) * geom.plotW;
  const Y = e => geom.plotY + (state.emax - e) / (state.emax - state.emin) * geom.plotH;
  const invX = px => state.xmin + (px - geom.plotX) / geom.plotW * (state.xmax - state.xmin);
  const invY = py => state.emax - (py - geom.plotY) / geom.plotH * (state.emax - state.emin);

  function roundRect(c, x, y, w, h, r) {
    r = Math.min(r, w / 2, h / 2);
    c.beginPath();
    c.moveTo(x + r, y);
    c.arcTo(x + w, y, x + w, y + h, r);
    c.arcTo(x + w, y + h, x, y + h, r);
    c.arcTo(x, y + h, x, y, r);
    c.arcTo(x, y, x + w, y, r);
    c.closePath();
  }

  /* ================= 体系 / 相区计算 ================= */
  function envKey() {
    return Object.keys(state.env).sort().map(k => {
      const e = state.env[k];
      return k + ":" + (e.on ? 1 : 0) + ":" + (e.c == null ? "-" : e.c.toPrecision(6)) + ":" + (e.ratio == null ? 1 : e.ratio);
    }).join(";");
  }
  function makeCfg() {
    return {
      main: state.main, selected: [...state.selected], env: state.env,
      Ctot: state.Ctot, Penv: state.Penv, xAxis: state.xAxis, pHfix: state.pHfix,
      couple: state.couple, customC: state.customC, customP: state.customP,
    };
  }
  function systemKey() {
    return JSON.stringify([state.main, [...state.selected].sort(), envKey(), state.Ctot.toPrecision(6),
      state.Penv, state.xAxis, state.pHfix, state.couple, JSON.stringify(state.customC)]);
  }

  function onBorder(A, B, eps) {
    eps = eps || 1e-7;
    const onX = (v) => Math.abs(v - state.xmin) < eps || Math.abs(v - state.xmax) < eps;
    const onE = (v) => Math.abs(v - state.emin) < eps || Math.abs(v - state.emax) < eps;
    return (onX(A[0]) && onX(B[0])) || (onE(A[1]) && onE(B[1]));
  }
  function pointInPoly(x, y, poly) {
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const xi = poly[i][0], yi = poly[i][1], xj = poly[j][0], yj = poly[j][1];
      const inter = ((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
      if (inter) inside = !inside;
    }
    return inside;
  }
  function distToSeg(px, py, x1, y1, x2, y2) {
    const dx = x2 - x1, dy = y2 - y1;
    const l2 = dx * dx + dy * dy;
    if (l2 < 1e-12) return Math.hypot(px - x1, py - y1);
    let t = ((px - x1) * dx + (py - y1) * dy) / l2;
    t = clamp(t, 0, 1);
    return Math.hypot(px - (x1 + t * dx), py - (y1 + t * dy));
  }

  function buildBoundaries(sys, polys) {
    const segs = [], seen = new Set();
    const eps = 1e-4;
    for (const p of polys) {
      const poly = p.poly;
      if (!poly || poly.length < 3) continue;
      const cx = p.centroid[0], cy = p.centroid[1];
      for (let i = 0; i < poly.length; i++) {
        const A = poly[i], B = poly[(i + 1) % poly.length];
        if (onBorder(A, B)) continue;
        const mx = (A[0] + B[0]) / 2, my = (A[1] + B[1]) / 2;
        let nx = -(B[1] - A[1]), ny = B[0] - A[0];
        const len = Math.hypot(nx, ny);
        if (len < 1e-9) continue;
        nx /= len; ny /= len;
        if ((cx - mx) * nx + (cy - my) * ny < 0) { nx = -nx; ny = -ny; }
        const pin = [mx + nx * eps, my + ny * eps];
        const pout = [mx - nx * eps, my - ny * eps];
        const winIn = sys.winner(pin[0], pin[1]);
        if (winIn !== p.id) continue;
        const winOut = sys.winner(pout[0], pout[1]);
        if (!winOut || winOut === p.id) continue;
        const k1 = A[0].toFixed(4) + "," + A[1].toFixed(4);
        const k2 = B[0].toFixed(4) + "," + B[1].toFixed(4);
        const pair = [p.id, winOut].sort().join("|");
        const key = pair + "|" + [k1, k2].sort().join("|");
        if (seen.has(key)) continue;
        seen.add(key);
        const spA = byId[p.id], spB = byId[winOut];
        segs.push({
          a: p.id, b: winOut, x1: A[0], e1: A[1], x2: B[0], e2: B[1],
          type: sys.classify(p.id, winOut),
          solid: (spA && spA.state === "s") || (spB && spB.state === "s"),
        });
      }
    }
    return segs;
  }

  function groupBoundaryChains(segs) {
    const chains = [];
    const tol = 2e-3;
    const near = (p, q) => Math.abs(p[0] - q[0]) < tol && Math.abs(p[1] - q[1]) < tol;
    for (const s of segs) {
      const pair = [s.a, s.b].sort().join("|");
      let ch = null;
      for (const c of chains) {
        if (c.pair !== pair) continue;
        const first = c.pts[0], last = c.pts[c.pts.length - 1];
        if (near(last, [s.x1, s.e1])) { c.pts.push([s.x2, s.e2]); ch = c; break; }
        if (near(last, [s.x2, s.e2])) { c.pts.push([s.x1, s.e1]); ch = c; break; }
        if (near(first, [s.x2, s.e2])) { c.pts.unshift([s.x1, s.e1]); ch = c; break; }
        if (near(first, [s.x1, s.e1])) { c.pts.unshift([s.x2, s.e2]); ch = c; break; }
      }
      if (!ch) {
        ch = { pair, a: s.a, b: s.b, type: s.type, solid: s.solid, pts: [[s.x1, s.e1], [s.x2, s.e2]], segs: [] };
        chains.push(ch);
      }
      ch.segs.push(s);
      s.chain = ch;
    }
    return chains;
  }

  function buildLabels(sys, polys) {
    const bySpecies = new Map();
    for (const p of polys) {
      if (!bySpecies.has(p.id)) bySpecies.set(p.id, []);
      bySpecies.get(p.id).push(p);
    }
    const cand = [];
    for (const [id, list] of bySpecies) {
      const sp = byId[id];
      let total = 0, cx = 0, cy = 0;
      for (const p of list) {
        total += p.area;
        cx += p.centroid[0] * p.area;
        cy += p.centroid[1] * p.area;
      }
      cx /= total || 1; cy /= total || 1;
      let rx0 = Infinity, rx1 = -Infinity, ry0 = Infinity, ry1 = -Infinity;
      for (const p of list) for (const q of p.poly) {
        rx0 = Math.min(rx0, q[0]); rx1 = Math.max(rx1, q[0]);
        ry0 = Math.min(ry0, q[1]); ry1 = Math.max(ry1, q[1]);
      }
      const bcx = (rx0 + rx1) / 2, bcy = (ry0 + ry1) / 2;
      list.sort((a, b) => b.area - a.area);
      const hostB = list.find(p => pointInPoly(bcx, bcy, p.poly));
      const hostC = list.find(p => pointInPoly(cx, cy, p.poly));
      const host = hostB || hostC ||
        list.slice().sort((a, b) => Math.hypot(a.centroid[0] - bcx, a.centroid[1] - bcy) - Math.hypot(b.centroid[0] - bcx, b.centroid[1] - bcy))[0];
      const targetX = hostB ? bcx : (hostC ? cx : host.centroid[0]);
      const targetY = hostB ? bcy : (hostC ? cy : host.centroid[1]);
      const hostPoly = host.poly;
      let bx0 = Infinity, bx1 = -Infinity, by0 = Infinity, by1 = -Infinity;
      for (const q of hostPoly) {
        bx0 = Math.min(bx0, q[0]); bx1 = Math.max(bx1, q[0]);
        by0 = Math.min(by0, q[1]); by1 = Math.max(by1, q[1]);
      }
      const pts = [];
      for (let ix = 1; ix <= 13; ix++) {
        for (let iy = 1; iy <= 13; iy++) {
          const x = bx0 + (bx1 - bx0) * ix / 14;
          const y = by0 + (by1 - by0) * iy / 14;
          if (!pointInPoly(x, y, hostPoly)) continue;
          const sx = X(x), sy = Y(y);
          let d = Infinity;
          for (let i = 0; i < hostPoly.length; i++) {
            const a = hostPoly[i], b = hostPoly[(i + 1) % hostPoly.length];
            d = Math.min(d, distToSeg(sx, sy, X(a[0]), Y(a[1]), X(b[0]), Y(b[1])));
          }
          const dc = Math.hypot(sx - X(targetX), sy - Y(targetY));
          pts.push({ x, y, sx, sy, d, score: d - dc * 0.08 });
        }
      }
      // 区域面积加权质心优先，保证标签落在相区中部而不是最大条带的边缘
      if (pointInPoly(targetX, targetY, hostPoly)) {
        const sx = X(targetX), sy = Y(targetY);
        let d = Infinity;
        for (let i = 0; i < hostPoly.length; i++) {
          const a = hostPoly[i], b = hostPoly[(i + 1) % hostPoly.length];
          d = Math.min(d, distToSeg(sx, sy, X(a[0]), Y(a[1]), X(b[0]), Y(b[1])));
        }
        pts.push({ x: targetX, y: targetY, sx, sy, d, score: d + 1000 });
      }
      pts.sort((a, b) => b.score - a.score);
      const suffix = sp.state === "aq" ? "" : "(" + sp.state + ")";
      cand.push({ id, pts, area: total, text: chemUni(sp.formula) + suffix, color: colorFor(id).text, sp });
    }
    cand.sort((a, b) => b.area - a.area);
    const placed = [];
    const labels = [];
    for (const c of cand) {
      const w = Math.max(26, c.text.length * 7.3 + 12), h = 18;
      let chosen = null;
      for (const p of c.pts) {
        if (p.sx - w / 2 < geom.plotX + 3 || p.sx + w / 2 > geom.plotX + geom.plotW - 3) continue;
        if (p.sy - h / 2 < geom.plotY + 3 || p.sy + h / 2 > geom.plotY + geom.plotH - 3) continue;
        let hit = false;
        for (const r of placed) {
          if (Math.abs(p.sx - r.x) < (w + r.w) / 2 + 5 && Math.abs(p.sy - r.y) < (h + r.h) / 2 + 4) { hit = true; break; }
        }
        if (!hit) { chosen = p; break; }
      }
      if (chosen) {
        placed.push({ x: chosen.sx, y: chosen.sy, w, h });
        labels.push({ id: c.id, x: chosen.sx, y: chosen.sy, text: c.text, color: c.color, sp: c.sp });
      }
    }
    return labels;
  }

  /* ================= pCl 横轴的 A 模式回退 =================
     解析多边形在“同斜率、多物种”的 pCl 情形下会丢相区；
     改用逐列下包络的矩形条带，边界/标签仍由同一批多边形生成。 */
  function buildPclPolygons(sys) {
    const N = 260;
    const out = [];
    const dx = (state.xmax - state.xmin) / N;
    for (let k = 0; k < N; k++) {
      const x0 = state.xmin + k * dx, x1 = x0 + dx, xm = (x0 + x1) / 2;
      const env = sys.envelope(xm, state.emin, state.emax);
      for (const seg of env) {
        const h = seg.e1 - seg.e0;
        if (h <= 1e-9) continue;
        out.push({
          id: seg.id, x0, x1,
          poly: [[x0, seg.e0], [x1, seg.e0], [x1, seg.e1], [x0, seg.e1]],
          area: dx * h,
          centroid: [xm, (seg.e0 + seg.e1) / 2],
        });
      }
    }
    return out;
  }
  function computePolygons(sys) {
    if (state.xAxis !== "pH") return buildPclPolygons(sys);
    return sys.dominance(state.xmin, state.xmax, state.emin, state.emax, 180);
  }

  /* ================= 场渲染（B/C） ================= */
  function buildField(sys, mode) {
    const nx = mode === "C" ? 220 : 200;
    const ny = Math.round(nx * geom.plotH / geom.plotW);
    const off = document.createElement("canvas");
    off.width = nx; off.height = ny;
    const octx2 = off.getContext("2d");
    const img = octx2.createImageData(nx, ny);
    const data = img.data;
    const th = theme();
    const bg = th.bgRgb;
    const Cdisp = state.Cdisp;
    for (let i = 0; i < nx; i++) {
      const x = state.xmin + (i + 0.5) / nx * (state.xmax - state.xmin);
      const mu = sys.muBasis(x);
      const G0 = {};
      for (const id of sys.selectedMain) G0[id] = sys.baseEnergy(id, x, mu);
      for (let j = 0; j < ny; j++) {
        const E = state.emax - (j + 0.5) / ny * (state.emax - state.emin);
        const sp = sys.speciate(x, E, { Ctot: state.Ctot, mu, G0, keepGas: true, maxIter: 12 });
        let r = bg[0], g = bg[1], b = bg[2];
        if (mode === "B") {
          let cur = [r, g, b];
          for (const p of sp.present) {
            const col = colorFor(p.id).rgbFill;
            let a = 0;
            if (p.gas) {
              if (p.p > 0) a = 0.55 * clamp(p.p / state.Penv, 0, 1);
            } else if (p.solid) {
              if (p.frac > 1e-4) a = 0.40;
            } else {
              if (p.conc > Cdisp) a = 0.36;
            }
            if (a > 0) cur = mixRgb(cur, col, a);
          }
          if (sp.solid && sp.solid.frac > 0) cur = mixRgb(cur, [255, 255, 255], Math.min(0.25, sp.solid.frac * 0.22));
          r = cur[0]; g = cur[1]; b = cur[2];
        } else {
          const items = [];
          let sum = 0;
          for (const p of sp.present) {
            if (p.gas) continue;
            let w = 0;
            if (p.solid) {
              if (p.frac > 1e-4) w = p.frac;
            } else {
              if (p.conc >= Cdisp) w = p.frac;
            }
            if (w > 0) { items.push([p.id, w]); sum += w; }
          }
          let acc = bg.slice();
          if (sum > 0) {
            let mix = [0, 0, 0];
            for (const [id, w] of items) {
              const col = colorFor(id).rgbFill;
              mix[0] += col[0] * w / sum; mix[1] += col[1] * w / sum; mix[2] += col[2] * w / sum;
            }
            acc = mixRgb(bg, mix, Math.min(1, sum));
          }
          // 气相按分压着色：p → P_env 时完全覆盖液相/固相颜色。
          const gases = sp.present.filter(p => p.gas && p.p > 0).sort((a, b) => a.p - b.p);
          for (const gp of gases) {
            const pg = clamp(gp.p / state.Penv, 0, 1);
            if (pg <= 1e-6) continue;
            acc = mixRgb(acc, colorFor(gp.id).rgbFill, pg);
          }
          if (sp.solid && sp.solid.frac > 0) acc = mixRgb(acc, [255, 255, 255], Math.min(0.25, sp.solid.frac * 0.20));
          r = acc[0]; g = acc[1]; b = acc[2];
        }
        const idx = (j * nx + i) * 4;
        data[idx] = clamp(Math.round(r), 0, 255);
        data[idx + 1] = clamp(Math.round(g), 0, 255);
        data[idx + 2] = clamp(Math.round(b), 0, 255);
        data[idx + 3] = 255;
      }
    }
    octx2.putImageData(img, 0, 0);
    return { canvas: off, nx, ny, mode };
  }

  /* ================= 绘制 ================= */
  function drawEmpty() {
    const th = theme();
    ctx.fillStyle = th.muted;
    ctx.font = '13px "Segoe UI","Microsoft YaHei",sans-serif';
    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    ctx.fillText("勾选左侧物种，开始绘制相图", geom.plotX + geom.plotW / 2, geom.plotY + geom.plotH / 2);
  }

  function drawGridAndAxes() {
    const th = theme();
    const x0 = geom.plotX, y0 = geom.plotY, w = geom.plotW, h = geom.plotH;
    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0, w, h); ctx.clip();
    // 网格
    if (state.showGrid) {
      const xStepMinor = state.xAxis === "pH" ? 1 : 0.5;
      const xStepMajor = state.xAxis === "pH" ? 2 : 1;
      const eStepMinor = 0.1, eStepMajor = 0.5;
      ctx.lineWidth = 1;
      for (let x = Math.ceil(state.xmin / xStepMinor) * xStepMinor; x <= state.xmax + 1e-9; x += xStepMinor) {
        const major = Math.abs(x / xStepMajor - Math.round(x / xStepMajor)) < 1e-6;
        ctx.strokeStyle = major ? th.grid : th.gridMinor;
        ctx.globalAlpha = major ? 0.75 : 0.55;
        ctx.beginPath(); ctx.moveTo(X(x), y0); ctx.lineTo(X(x), y0 + h); ctx.stroke();
      }
      for (let e = Math.ceil(state.emin / eStepMinor) * eStepMinor; e <= state.emax + 1e-9; e += eStepMinor) {
        const major = Math.abs(e / eStepMajor - Math.round(e / eStepMajor)) < 1e-6;
        ctx.strokeStyle = major ? th.grid : th.gridMinor;
        ctx.globalAlpha = major ? 0.75 : 0.55;
        ctx.beginPath(); ctx.moveTo(x0, Y(e)); ctx.lineTo(x0 + w, Y(e)); ctx.stroke();
      }
      ctx.globalAlpha = 1;
    }
    ctx.restore();

    // 轴刻度
    ctx.save();
    ctx.fillStyle = th.muted;
    ctx.font = '10px "JetBrains Mono",Consolas,monospace';
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    const xStep = state.xAxis === "pH" ? 2 : 1;
    for (let x = Math.ceil(state.xmin / xStep) * xStep; x <= state.xmax + 1e-9; x += xStep) {
      const sx = X(x);
      if (sx < x0 - 1 || sx > x0 + w + 1) continue;
      ctx.strokeStyle = th.line; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(sx, y0 + h); ctx.lineTo(sx, y0 + h + 4); ctx.stroke();
      ctx.fillText(state.xAxis === "pH" ? String(Math.round(x * 10) / 10) : String(Math.round(x * 10) / 10), sx, y0 + h + 7);
    }
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    const eStep = 0.5;
    for (let e = Math.ceil(state.emin / eStep) * eStep; e <= state.emax + 1e-9; e += eStep) {
      const sy = Y(e);
      if (sy < y0 - 1 || sy > y0 + h + 1) continue;
      ctx.strokeStyle = th.line; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x0 - 4, sy); ctx.lineTo(x0, sy); ctx.stroke();
      ctx.fillText(fmtNum(e, 1), x0 - 7, sy);
    }
    ctx.restore();

    // 轴标题
    ctx.save();
    ctx.fillStyle = th.text;
    ctx.font = '600 11.5px "Segoe UI","Microsoft YaHei",sans-serif';
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.fillText(state.xAxis === "pH" ? "pH" : state.xAxis === "pCl" ? "lg[Cl⁻]" : "lg[NH₃]", x0 + w / 2, y0 + h + 26);
    ctx.translate(15, y0 + h / 2);
    ctx.rotate(-PI / 2);
    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    ctx.fillText("E / V (vs SHE)", 0, 0);
    ctx.restore();

    // 边框
    ctx.strokeStyle = th.line; ctx.lineWidth = 1;
    ctx.strokeRect(x0 + 0.5, y0 + 0.5, w - 1, h - 1);
    // 25°C 角标
    ctx.fillStyle = th.muted; ctx.font = '9px "JetBrains Mono",Consolas,monospace';
    ctx.textAlign = "right"; ctx.textBaseline = "top";
    ctx.fillText("25 °C", x0 + w - 7, y0 + 6);
  }

  function drawMediumLines() {
    const th = theme();
    const x0 = geom.plotX, y0 = geom.plotY, w = geom.plotW, h = geom.plotH;
    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0, w, h); ctx.clip();
    const nernst = C.NERNST;
    const pH = state.xAxis === "pH" ? null : state.pHfix;

    if (state.showWater) {
      ctx.strokeStyle = th.water; ctx.lineWidth = 1; ctx.setLineDash([5, 4]); ctx.globalAlpha = 0.62;
      for (const [name, fn] of [
        ["H⁺/H₂", x => -nernst * x],
        ["O₂/H₂O", x => 1.229 - nernst * x],
      ]) {
        ctx.beginPath();
        for (let k = 0; k <= 40; k++) {
          const x = state.xmin + (state.xmax - state.xmin) * k / 40;
          const e = state.xAxis === "pH" ? fn(x) : fn(pH);
          if (k === 0) ctx.moveTo(X(x), Y(e)); else ctx.lineTo(X(x), Y(e));
        }
        ctx.stroke();
      }
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
      // 标签
      ctx.font = '9px "Cambria Math","Segoe UI",sans-serif';
      ctx.textAlign = "right"; ctx.textBaseline = "middle";
      for (const [name, fn] of [["H⁺/H₂", x => -nernst * x], ["O₂/H₂O", x => 1.229 - nernst * x]]) {
        const e = state.xAxis === "pH" ? fn(state.xmax - 0.2) : fn(pH);
        if (e < state.emin || e > state.emax) continue;
        const sx = X(state.xmax - 0.25), sy = Y(e);
        const tw = ctx.measureText(name).width;
        ctx.fillStyle = th.plot; ctx.globalAlpha = 0.75;
        ctx.fillRect(sx - tw - 3, sy - 7, tw + 5, 14);
        ctx.globalAlpha = 0.8; ctx.fillStyle = th.water;
        ctx.fillText(name, sx, sy);
      }
      ctx.globalAlpha = 1;
    }

    if (state.showMedium) {
      // pKa 线
      if (state.xAxis === "pH" && state.sys) {
        const pks = state.sys.pKaList();
        ctx.font = '8px "Cambria Math","Segoe UI",sans-serif';
        ctx.textAlign = "center"; ctx.textBaseline = "top";
        let row = 0;
        const placed = [];
        for (const pk of pks) {
          if (pk.pH < state.xmin || pk.pH > state.xmax) continue;
          const sx = X(pk.pH);
          ctx.strokeStyle = th.medium; ctx.globalAlpha = 0.33; ctx.setLineDash([1, 3]); ctx.lineWidth = 1;
          ctx.beginPath(); ctx.moveTo(sx, y0); ctx.lineTo(sx, y0 + h); ctx.stroke();
          ctx.setLineDash([]);
          const label = chemUni(byId[pk.from].formula) + "/" + chemUni(byId[pk.to].formula);
          const tw = ctx.measureText(label).width;
          if (placed.some(r => Math.abs(r.x - sx) < (r.w + tw) / 2 + 6 && r.row === row)) row = 1 - row;
          const ty = y0 + 6 + row * 13 + 4;
          ctx.globalAlpha = 0.75;
          ctx.fillStyle = th.plot;
          ctx.fillRect(sx - tw / 2 - 2, ty - 2, tw + 4, 11);
          ctx.fillStyle = th.medium;
          ctx.fillText(label, sx, ty);
          placed.push({ x: sx, w: tw, row });
        }
        ctx.globalAlpha = 1;
      }
      // 电对线
      const redox = D.ENV.filter(c => c.redox && state.env[c.id] && state.env[c.id].on);
      ctx.font = '9px "Cambria Math","Segoe UI",sans-serif';
      ctx.textAlign = "right"; ctx.textBaseline = "middle";
      for (const comp of redox) {
        if (!state.sys) continue;
        const midX = state.xAxis === "pH" ? (state.xmin + state.xmax) / 2 : state.pHfix;
        const info = state.sys.redoxLine({ ...comp.redox, id: comp.id }, midX);
        if (!info || !isFinite(info.E)) continue;
        const slope = info.r.slope || 0;
        const eAt = x => info.E + slope * (x - midX);
        const ex1 = state.xAxis === "pH" ? eAt(state.xmin) : info.E;
        const ex2 = state.xAxis === "pH" ? eAt(state.xmax) : info.E;
        if ((ex1 < state.emin && ex2 < state.emin) || (ex1 > state.emax && ex2 > state.emax)) continue;
        ctx.strokeStyle = th.medium; ctx.globalAlpha = 0.42; ctx.setLineDash([2, 4]); ctx.lineWidth = 1;
        ctx.beginPath();
        if (state.xAxis === "pH") { ctx.moveTo(X(state.xmin), Y(ex1)); ctx.lineTo(X(state.xmax), Y(ex2)); }
        else { ctx.moveTo(X(state.xmin), Y(info.E)); ctx.lineTo(X(state.xmax), Y(info.E)); }
        ctx.stroke(); ctx.setLineDash([]);
        const ey = state.xAxis === "pH" ? ex2 : info.E;
        if (ey >= state.emin && ey <= state.emax) {
          const sx = X(state.xmax - 0.25), sy = Y(ey);
          const label = comp.label;
          const tw = ctx.measureText(label).width;
          ctx.globalAlpha = 0.75; ctx.fillStyle = th.plot;
          ctx.fillRect(sx - tw - 3, sy - 7, tw + 5, 14);
          ctx.globalAlpha = 0.8; ctx.fillStyle = th.medium;
          ctx.fillText(label, sx, sy);
        }
      }
      ctx.globalAlpha = 1;
    }

    // 所有电对线（含亚稳）
    if (state.showAllRedox && state.sys) {
      const ids = state.sys.selectedMain;
      let n = 0;
      ctx.strokeStyle = th.medium; ctx.globalAlpha = 0.16; ctx.setLineDash([1, 4]); ctx.lineWidth = 1;
      for (let i = 0; i < ids.length && n < 36; i++) {
        for (let j = i + 1; j < ids.length && n < 36; j++) {
          const r = state.sys.balance(ids[i], ids[j]);
          if (!r || !(r.nE > 0) || r.nE > 8 || r.E0 == null) continue;
          const ex1 = state.xAxis === "pH" ? r.E0 + (r.slope || 0) * state.xmin : r.E0 + (r.slope || 0) * state.pHfix;
          const ex2 = state.xAxis === "pH" ? r.E0 + (r.slope || 0) * state.xmax : r.E0 + (r.slope || 0) * state.pHfix;
          if ((ex1 < state.emin && ex2 < state.emin) || (ex1 > state.emax && ex2 > state.emax)) continue;
          ctx.beginPath();
          ctx.moveTo(X(state.xmin), Y(ex1)); ctx.lineTo(X(state.xmax), Y(ex2));
          ctx.stroke(); n++;
        }
      }
      ctx.setLineDash([]); ctx.globalAlpha = 1;
    }
    ctx.restore();
  }

  function drawFills() {
    const th = theme();
    ctx.save();
    ctx.beginPath(); ctx.rect(geom.plotX, geom.plotY, geom.plotW, geom.plotH); ctx.clip();
    if (state.mode === "A" || !state.field) {
      if (!state.sys || !state.sys.selectedMain.length) { drawEmpty(); ctx.restore(); return; }
      // A 模式：精确多边形
      for (const p of state.polygons) {
        const col = colorFor(p.id);
        ctx.beginPath();
        for (let i = 0; i < p.poly.length; i++) {
          const q = p.poly[i];
          const sx = X(q[0]), sy = Y(q[1]);
          if (i === 0) ctx.moveTo(sx, sy); else ctx.lineTo(sx, sy);
        }
        ctx.closePath();
        ctx.fillStyle = col.fill;
        ctx.fill();
        // 用同色细描边覆盖相邻条带接缝
        ctx.strokeStyle = col.fill;
        ctx.lineWidth = 0.7;
        ctx.stroke();
      }
      } else {
      // B/C 场
      try {
        ctx.imageSmoothingEnabled = true;
        ctx.globalAlpha = state.mode === "B" ? 0.94 : 1;
        ctx.drawImage(state.field.canvas, geom.plotX, geom.plotY, geom.plotW, geom.plotH);
        ctx.globalAlpha = 1;
      } catch (e) {
        ctx.restore(); return;
      }
    }
    ctx.restore();
  }

  function drawBoundaries() {
    if (!state.showBoundary) return;
    ctx.save();
    ctx.beginPath(); ctx.rect(geom.plotX, geom.plotY, geom.plotW, geom.plotH); ctx.clip();
    ctx.lineCap = "round"; ctx.lineJoin = "round";
    for (const ch of state.boundaryChains) {
      if (!ch.pts.length) continue;
      let visible = false;
      ctx.beginPath();
      for (let i = 0; i < ch.pts.length; i++) {
        const p = ch.pts[i];
        const sx = X(p[0]), sy = Y(p[1]);
        if (i === 0) ctx.moveTo(sx, sy); else ctx.lineTo(sx, sy);
        if (sx >= geom.plotX - 2 && sx <= geom.plotX + geom.plotW + 2 && sy >= geom.plotY - 2 && sy <= geom.plotY + geom.plotH + 2) visible = true;
      }
      if (!visible) continue;
      const col = BOUNDARY_COLORS[ch.type] || BOUNDARY_COLORS.stable;
      ctx.strokeStyle = col;
      ctx.globalAlpha = ch.solid ? 0.78 : 0.62;
      ctx.lineWidth = ch.solid ? 1.25 : 1.0;
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
    ctx.restore();
  }

  function drawLabels() {
    if (!state.showLabels) return;
    ctx.save();
    ctx.beginPath(); ctx.rect(geom.plotX, geom.plotY, geom.plotW, geom.plotH); ctx.clip();
    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    for (const L of state.labels) {
      if (L.x < geom.plotX || L.x > geom.plotX + geom.plotW) continue;
      if (L.y < geom.plotY || L.y > geom.plotY + geom.plotH) continue;
      ctx.font = '600 12px "Cambria Math","Segoe UI","Microsoft YaHei",serif';
      const w = ctx.measureText(L.text).width;
      ctx.fillStyle = theme().labelBg;
      roundRect(ctx, L.x - w / 2 - 4, L.y - 9, w + 8, 18, 5);
      ctx.fill();
      ctx.fillStyle = L.color;
      ctx.fillText(L.text, L.x, L.y + 0.5);
    }
    ctx.restore();
  }

  function draw() {
    if (!geom) setupCanvas();
    const th = theme();
    // 背景
    ctx.clearRect(0, 0, geom.W, geom.H);
    ctx.fillStyle = th.plot;
    ctx.fillRect(0, 0, geom.W, geom.H);
    const g = ctx.createLinearGradient(0, geom.plotY, 0, geom.plotY + geom.plotH);
    g.addColorStop(0, th.plotTop); g.addColorStop(1, th.plotBottom);
    ctx.fillStyle = g;
    ctx.fillRect(geom.plotX, geom.plotY, geom.plotW, geom.plotH);

    // 相区
    drawFills();
    if (state.sys && state.sys.selectedMain.length) {
      // 网格在填色之上、相界线之下
      drawGridAndAxes();
      drawMediumLines();
      drawBoundaries();
      drawLabels();
    } else {
      drawGridAndAxes();
    }
    drawHighlight();
  }

  /* ================= 悬停高亮 / overlay ================= */
  function drawHighlight() {
    if (!geom) return;
    octx.clearRect(0, 0, geom.W, geom.H);
    const h = state.hover;
    if (!h) return;
    octx.save();
    octx.beginPath(); octx.rect(geom.plotX, geom.plotY, geom.plotW, geom.plotH); octx.clip();
    if (h.kind === "boundary") {
      const ch = h.seg.chain;
      const pts = ch ? ch.pts : [[h.seg.x1, h.seg.e1], [h.seg.x2, h.seg.e2]];
      octx.strokeStyle = theme().accent; octx.globalAlpha = 0.95; octx.lineWidth = 2.6; octx.lineCap = "round"; octx.lineJoin = "round";
      octx.beginPath();
      for (let i = 0; i < pts.length; i++) {
        const p = pts[i], sx = X(p[0]), sy = Y(p[1]);
        if (i === 0) octx.moveTo(sx, sy); else octx.lineTo(sx, sy);
      }
      octx.stroke();
    } else if (h.kind === "point") {
      const x = X(h.x), y = Y(h.e);
      octx.strokeStyle = theme().accent; octx.globalAlpha = 0.9; octx.lineWidth = 1.5;
      octx.beginPath(); octx.arc(x, y, 6, 0, PI * 2); octx.stroke();
      octx.globalAlpha = 0.5;
      octx.beginPath(); octx.moveTo(x - 10, y); octx.lineTo(x - 4, y); octx.moveTo(x + 4, y); octx.lineTo(x + 10, y);
      octx.moveTo(x, y - 10); octx.lineTo(x, y - 4); octx.moveTo(x, y + 4); octx.lineTo(x, y + 10); octx.stroke();
    }
    octx.restore();
  }

  /* ================= 悬停信息 ================= */
  const tip = $("tooltip");
  function nearestBoundary(px, py) {
    let best = null, bd = Infinity;
    for (const s of state.boundaries) {
      const x1 = X(s.x1), y1 = Y(s.e1), x2 = X(s.x2), y2 = Y(s.e2);
      if (Math.min(x1, x2) - 20 > px || Math.max(x1, x2) + 20 < px) continue;
      if (Math.min(y1, y2) - 20 > py || Math.max(y1, y2) + 20 < py) continue;
      const d = distToSeg(px, py, x1, y1, x2, y2);
      if (d < bd) { bd = d; best = s; }
    }
    return bd <= state.hoverPx ? best : null;
  }
  const TYPE_NAME = { redox: "氧化还原", acidbase: "酸碱", precipitation: "沉淀", complex: "配位", phase: "相变", stable: "稳定相界" };
  function reactionHTML(r) {
    const side = arr => arr.map(x => {
      const coef = x.coef > 1.0001 ? '<span style="font-family:var(--mono);font-size:.86em;color:var(--muted)">' + Math.round(x.coef) + '</span>&thinsp;' : "";
      return coef + formulaById(x.id, state.fmode);
    }).join(" + ");
    return side(r.left) + ' <span style="color:var(--accent);font-weight:600">⇌</span> ' + side(r.right);
  }
  function reactionForTip(seg) {
    const r = state.sys.balance(seg.a, seg.b);
    if (!r) return null;
    const pH = state.xAxis === "pH" ? (seg.x1 + seg.x2) / 2 : state.pHfix;
    const useAlk = pH > 7 && r.hasH;
    return { r, shown: useAlk ? state.sys.alkalineForm(r) : r, pH, useAlk };
  }
  function showTip(html, px, py) {
    tip.innerHTML = html;
    tip.style.display = "block";
    const wrap = $("canvas-wrap");
    const tw = tip.offsetWidth, th = tip.offsetHeight;
    let left = px + 15, top = py + 15;
    if (left + tw > wrap.clientWidth - 6) left = px - tw - 15;
    if (top + th > wrap.clientHeight - 6) top = py - th - 15;
    tip.style.left = Math.max(6, left) + "px";
    tip.style.top = Math.max(6, top) + "px";
  }
  function hideTip() { tip.style.display = "none"; }

  function onMove(e) {
    if (!geom || !state.sys) return;
    const rect = plot.getBoundingClientRect();
    const px = e.clientX - rect.left, py = e.clientY - rect.top;
    if (px < geom.plotX || px > geom.plotX + geom.plotW || py < geom.plotY || py > geom.plotY + geom.plotH) {
      state.hover = null; hideTip(); drawHighlight(); return;
    }
    const x = invX(px), ee = invY(py);
    const seg = nearestBoundary(px, py);
    if (seg) {
      const info = reactionForTip(seg);
      state.hover = { kind: "boundary", seg };
      drawHighlight();
      if (!info) return hideTip();
      const spA = byId[seg.a], spB = byId[seg.b];
      const cA = colorFor(seg.a), cB = colorFor(seg.b);
      const meta = [];
      if (info.r.nE > 0) meta.push("E° = <b>" + fmtE(info.r.E0) + "</b>");
      if (info.r.slope != null && Math.abs(info.r.slope) > 1e-5) meta.push("斜率 = <b>" + fmtNum(info.r.slope, 3) + " V/pH</b>");
      if (info.r.pHstd != null) meta.push("pH 边界 = <b>" + fmtNum(info.r.pHstd, 2) + "</b>");
      meta.push("log K = <b>" + fmtLogK(info.r.logK) + "</b>");
      const html =
        '<div class="thead"><span class="dot" style="background:' + cA.trueColor + '"></span>' + chemUni(spA.formula) +
        ' <span style="color:var(--muted);font-weight:400">⇌</span> <span class="dot" style="background:' + cB.trueColor + '"></span>' + chemUni(spB.formula) +
        '<span class="tbadge">' + (TYPE_NAME[seg.type] || "相界") + (seg.solid ? " · 固相界" : "") + '</span></div>' +
        '<div class="tmeta">pH ' + fmtNum(info.pH, 2) + ' · E ' + fmtNum(ee, 3) + ' V</div>' +
        '<div class="teq">' + reactionHTML(info.shown) + '</div>' +
        '<div class="tmeta">' + meta.join(" · ") + '</div>' +
        '<div class="tnote">由 ΔGf° 现场配平' + (info.useAlk ? "；碱性写法：H⁺ + OH⁻ ⇌ H₂O" : "") + '。' +
        (spA.derived || spB.derived ? " 含「推得」数据，见物种行。" : " 数据：NBS 1982。") + '</div>';
      return showTip(html, px, py);
    }

    const win = state.sys.winner(x, ee);
    if (!win) { state.hover = null; hideTip(); drawHighlight(); return; }
    state.hover = { kind: "point", x, e: ee };
    drawHighlight();
    const sp = byId[win];
    const col = colorFor(win);
    const g = state.sys.gA(win, x, ee);
    const meta = ["Ĝ = <b>" + fmtG(g) + "</b>"];
    const parts = [];
    let body = "";
    if (state.mode === "A") {
      const lab = sp.state === "aq" ? "相区浓度" : sp.state === "s" ? "固相 · 承载 C_tot" : sp.state === "g" ? "气相分压" : "相";
      const val = sp.state === "g" ? fmtP(state.Penv)
        : fmtConc(state.customC[win] != null ? Math.pow(10, state.customC[win]) : state.Ctot);
      body += '<div class="trow"><span>' + lab + '</span><b>' + val + '</b></div>';
    } else {
      const spec = state.sys.speciate(x, ee, { Ctot: state.Ctot, keepGas: true, maxIter: 18 });
      const list = spec.present.slice().sort((a, b) => (b.frac || 0) - (a.frac || 0));
      for (const p of list.slice(0, 5)) {
        const s2 = byId[p.id];
        if (!s2) continue;
        const c2 = colorFor(p.id);
        let val = "";
        if (p.gas) val = "p = " + fmtP(p.p);
        else if (p.solid) val = "固相分数 " + (p.frac * 100).toFixed(1) + "%";
        else val = fmtConc(p.conc);
        parts.push('<span style="display:inline-flex;align-items:center;gap:3px"><i style="width:7px;height:7px;border-radius:50%;background:' + c2.trueColor + '"></i>' + formulaFor(s2, state.fmode) + ' <span style="color:var(--muted)">' + val + '</span></span>');
      }
      if (spec.solid) body += '<div class="trow"><span>固相</span><b>' + formulaById(spec.solid.id, state.fmode) + ' · ' + (spec.solid.frac * 100).toFixed(1) + '%</b></div>';
      if (spec.gas) {
        body += '<div class="trow"><span>气相分压</span><b>' + fmtP(spec.gas.p) + '</b></div>';
        if (state.V0 > 0 || spec.gas.p >= state.Penv * 0.999) {
          try {
            const gs = state.sys.gasState(x, ee, { Ctot: state.Ctot, V0: state.V0, Vsol: 1 });
            if (gs && gs.Vg > 0) body += '<div class="trow"><span>气相体积</span><b>' + gs.Vg.toFixed(3) + ' L</b></div>';
          } catch (e2) { /* ignore */ }
        }
      }
      body += '<div class="mini-legend">' + parts.join("") + '</div>';
    }
    const html =
      '<div class="thead"><span class="dot" style="background:' + col.trueColor + '"></span>' + formulaFor(sp, state.fmode) +
      '<span class="tbadge">' + (sp.state === "s" ? "固相" : sp.state === "g" ? "气相" : sp.state === "l" ? "液相" : "溶液") + '</span></div>' +
      '<div class="tmeta">pH ' + fmtNum(x, 2) + ' · E ' + fmtNum(ee, 3) + ' V · ' + meta.join(" · ") + '</div>' +
      '<div class="trow"><span>' + sp.name + '</span><b>' + (sp.derived ? "推得" : "NBS") + '</b></div>' +
      body;
    showTip(html, px, py);
  }
  function onLeave() { state.hover = null; hideTip(); drawHighlight(); }
  function hoverAt(x, e) {
    if (!geom) return;
    const rect = plot.getBoundingClientRect();
    onMove({ clientX: rect.left + X(x), clientY: rect.top + Y(e) });
  }

  /* ================= 图例 / 反应面板 ================= */
  function renderLegend() {
    const el = $("legend");
    el.innerHTML = "";
    const ids = [...state.selected].filter(id => byId[id] && byId[id].el[state.main] > 0 && !byId[id].solvent);
    if (!ids.length) {
      el.innerHTML = '<span class="muted" style="font-size:11px">尚未勾选物种。</span>';
      return;
    }
    const dominant = new Set(state.polygons.map(p => p.id));
    for (const id of ids) {
      const sp = byId[id];
      const col = colorFor(id);
      const chip = document.createElement("span");
      chip.className = "lg" + (dominant.has(id) ? "" : " dim");
      chip.title = sp.name + " · " + (sp.derived ? sp.note || "推得" : "NBS 1982");
      const sw = document.createElement("input");
      sw.type = "color"; sw.className = "sw"; sw.value = col.trueColor;
      sw.title = "自定义颜色";
      sw.addEventListener("input", () => {
        state.customColors[id] = sw.value;
        scheduleRender(0);
      });
      const fm = document.createElement("span");
      fm.className = "fm chem";
      fm.innerHTML = renderFormula(sp.formula, state.flegend);
      const nm = document.createElement("span");
      nm.className = "nm"; nm.textContent = sp.name;
      const tag = document.createElement("span");
      tag.className = "tag";
      tag.textContent = sp.state === "aq" ? "" : "(" + sp.state + ")";
      chip.appendChild(sw); chip.appendChild(fm); chip.appendChild(tag); chip.appendChild(nm);
      el.appendChild(chip);
    }
  }

  function renderReactions() {
    const el = $("reactions");
    const prevScroll = el.scrollTop;
    el.innerHTML = "";
    const groups = state.groups;
    if (!groups) { el.innerHTML = '<div class="rxn-empty">勾选物种后，这里列出现场配平的反应式。</div>'; return; }
    const defs = [
      { key: "stable", name: "稳定相界", desc: "相邻优势相之间", limit: 48 },
      { key: "redox", name: "氧化还原", desc: "含 e⁻ 的半反应", limit: 40 },
      { key: "acidbase", name: "酸碱", desc: "质子转移", limit: 32 },
      { key: "precipitation", name: "沉淀 / 溶解", desc: "固相与溶液", limit: 32 },
      { key: "complex", name: "配位", desc: "配合物生成 / 解离", limit: 32 },
      { key: "phase", name: "相变 / 水合", desc: "同组成固相之间", limit: 16 },
    ];
    let total = 0;
    for (const d of defs) total += (groups[d.key] || []).length;
    let shown = 0;
    for (const d of defs) {
      const arr = groups[d.key] || [];
      if (!arr.length) continue;
      const sec = document.createElement("section");
      sec.className = "rxn-group";
      const h = document.createElement("h3");
      h.innerHTML = d.name + ' <span class="n">' + arr.length + '</span><span class="desc">' + d.desc + '</span>';
      sec.appendChild(h);
      for (const r of arr.slice(0, d.limit)) {
        shown++;
        const spA = byId[r.a], spB = byId[r.b];
        const pH = state.xAxis === "pCl" ? state.pHfix : (r.at != null ? r.at : (r.pHstd != null ? r.pHstd : 7));
        const useAlk = pH > 7 && r.hasH;
        const rr = useAlk ? state.sys.alkalineForm(r) : r;
        const card = document.createElement("article");
        card.className = "rxn " + r.type;
        const meta = [];
        if (r.nE > 0) meta.push("E° <b>" + fmtE(r.E0) + "</b>");
        if (r.slope != null && Math.abs(r.slope) > 1e-5) meta.push("斜率 <b>" + fmtNum(r.slope, 3) + " V/pH</b>");
        if (r.at != null) {
          meta.push((state.xAxis === "pCl" ? "相界 lg[Cl⁻]" : "相界 pH") + " <b>" + fmtNum(r.at, 2) + "</b>");
          if (r.atE != null) meta.push("相界 E <b>" + fmtNum(r.atE, 3) + " V</b>");
        } else if (r.pHstd != null) {
          meta.push((Math.abs(Math.abs(r.cH) - 1) < 0.01 ? "pKa" : "pH 边界") + " <b>" + fmtNum(r.pHstd, 2) + "</b>");
        }
        if (r.at != null && r.type === "acidbase" && Math.abs(Math.abs(r.cH) - 1) < 0.01) {
          meta.push("pKa（标准） <b>" + fmtNum(-r.logK, 2) + "</b>");
        }
        meta.push("log K <b>" + fmtLogK(r.logK) + "</b>");
        meta.push('<span class="tag">' + (TYPE_NAME[r.type] || r.type) + '</span>');
        if (useAlk) meta.push('<span class="tag">碱性写法</span>');
        const eq = document.createElement("div");
        eq.className = "eq";
        eq.innerHTML = reactionHTML(rr);
        const md = document.createElement("div");
        md.className = "meta";
        md.innerHTML = meta.join(" · ");
        card.appendChild(eq); card.appendChild(md);
        sec.appendChild(card);
      }
      if (arr.length > d.limit) {
        const more = document.createElement("div");
        more.className = "rxn-empty";
        more.textContent = "仅显示前 " + d.limit + " 条；缩小物种选择范围可减少反应数。";
        sec.appendChild(more);
      }
      el.appendChild(sec);
    }
    if (!shown) el.innerHTML = '<div class="rxn-empty">当前选择下没有可列出的反应。尝试勾选更多物种或启用配体。</div>';
    $("reaction-count").textContent = total + " 条 · 显示 " + shown;
    el.scrollTop = prevScroll;
  }

  /* ================= 图表元信息 / 警告 ================= */
  function renderMeta() {
    const ligandNames = D.ENV.filter(c => state.env[c.id] && state.env[c.id].on && !c.redox
      && (c.groups.includes("ligand") || c.groups.includes("precipitant"))
      && !(c.id === "Cl-" && state.xAxis === "pCl")
      && !(c.id === "NH3" && state.xAxis === "pNH3")).slice(0, 2).map(c => c.label);
    const suffix = ligandNames.length ? "–" + ligandNames.join("–") : "";
    $("chart-title").innerHTML = state.main + suffix + "–H₂O 体系";
    const nPha = new Set(state.polygons.map(p => p.id)).size;
    const lines = [
      "25 °C · 理想溶液",
      "C<sub>tot</sub> = " + fmtConc(state.Ctot) + " · P<sub>env</sub> = " + fmtP(state.Penv) + (state.V0 > 0 ? " · V₀ = " + state.V0.toFixed(1) + " L" : ""),
      "模式 " + state.mode + " · " + nPha + " 个优势相" + (state.xAxis !== "pH" ? " · pH 固定 " + fmtNum(state.pHfix, 2) : ""),
    ];
    $("chart-meta").innerHTML = lines.join("<br>");
    const w = [];
    if (!state.sys) w.push("正在初始化计算…");
    else {
      if (state.sys.unbalanced.length) w.push("以下物种无法在当前基组中配平：" + state.sys.unbalanced.join("、"));
      const solvent = state.selected && [...state.selected].some(id => byId[id] && byId[id].solvent);
      if (solvent) w.push("H₂O 作为溶剂不参与相区比较。");
      if (!state.selected.size) w.push("未勾选物种；从左侧「典型物种」开始。");
      const filtered = [...state.selected].filter(id => byId[id] && byId[id].el[state.main] > 0 && !byId[id].solvent && !state.sys.selectedMain.includes(id));
      if (filtered.length) w.push("未参与相区比较：" + filtered.map(id => chemUni(byId[id].formula)).join("、"));
      if (Object.keys(state.customC).length) w.push("已启用 " + Object.keys(state.customC).length + " 个物种的自定义浓度；总浓度约束部分解除。");
    }
    $("warnbox").innerHTML = w.map(s => '<div style="margin:2px 0">' + s + '</div>').join("");
  }

  /* ================= 主渲染 ================= */
  let renderTimer = null, rafId = null;
  function scheduleRender(delay) {
    if (renderTimer) clearTimeout(renderTimer);
    renderTimer = setTimeout(() => {
      if (rafId) cancelAnimationFrame(rafId);
      rafId = requestAnimationFrame(() => { rafId = null; renderAll(); });
    }, delay == null ? 35 : delay);
  }

  function renderAll() {
    if (!geom) setupCanvas();
    const t0 = performance.now();
    const sk = systemKey();
    if (sk !== state.sysKey) {
      state.sys = C.createSystem(makeCfg());
      state.sysKey = sk;
      state.polygons = [];
      state.boundaries = [];
      state.labels = [];
      state.field = null;
      state.fieldKey = null;
    }
    const sys = state.sys;
    if (!sys) return;
    const polyKey = JSON.stringify([sk, state.xmin, state.xmax, state.emin, state.emax, 180]);
    if (!state.polygons.length && sys.selectedMain.length) {
      state.polygons = computePolygons(sys);
      state.boundaries = buildBoundaries(sys, state.polygons);
      state.boundaryChains = groupBoundaryChains(state.boundaries);
      state.labels = buildLabels(sys, state.polygons);
      state.groups = sys.reactionList({ boundaries: state.boundaries });
      if (state.palMode === "random") recolorRandom();
      state._polyKey = polyKey;
      state.fieldKey = null;
      state.field = null;
    } else if (!sys.selectedMain.length) {
      state.polygons = []; state.boundaries = []; state.boundaryChains = []; state.labels = [];
      state.groups = sys.reactionList({ boundaries: [] });
      state.field = null; state.fieldKey = null;
    }
    // 坐标轴范围变化但体系未变：重算多边形
    if (state._polyKey !== polyKey) {
      state._polyKey = polyKey;
      if (sys.selectedMain.length) {
        state.polygons = computePolygons(sys);
        state.boundaries = buildBoundaries(sys, state.polygons);
        state.labels = buildLabels(sys, state.polygons);
        state.groups = sys.reactionList({ boundaries: state.boundaries });
        state.fieldKey = null; state.field = null;
      }
    }
    if (state.mode !== "A" && sys.selectedMain.length) {
      const fk = JSON.stringify([state.sysKey, state.mode, state._polyKey, state.Cdisp, state.Penv, state.Ctot,
        state.palMode, state.palSeed, state.customColors]);
      if (fk !== state.fieldKey) {
        $("render-status").textContent = "计算场…";
        state.field = buildField(sys, state.mode);
        state.fieldKey = fk;
      }
    } else {
      state.field = null; state.fieldKey = null;
    }
    draw();
    renderLegend();
    renderReactions();
    renderMeta();
    state.renderMs = performance.now() - t0;
    $("render-status").textContent = state.renderMs > 20 ? Math.round(state.renderMs) + " ms" : "";
  }

  /* ================= UI 构建 ================= */
  const DEFAULT_SELECTED = {
    Mg: ["Mg", "Mg2+", "Mg(OH)2", "MgO"],
    Al: ["Al", "Al3+", "Al(OH)3", "Al(OH)4-"],
    P: ["H3PO4", "H2PO4-", "HPO42-", "PO43-", "P4", "PH3"],
    S: ["S8", "H2S", "HS-", "S2-", "SO2", "HSO3-", "SO32-", "S2O32-"],
    Cl: ["Cl-", "Cl2", "Cl2g", "HClO", "ClO-"],
    I: ["I-", "I3-", "I2", "IO3-"],
    N: ["N2", "NH4+", "NH3", "HNO2", "NO2-", "NO"],
    C: ["CO2", "HCO3-", "CO32-", "CO", "HCOOH", "HCOO-"],
    V: ["V", "V2+", "V3+", "VO2+", "VO2plus", "V2O5"],
    Cr: ["Cr", "Cr3+", "Cr(OH)3", "Cr2O3", "CrO42-", "Cr2O72-"],
    Mn: ["Mn", "Mn2+", "Mn(OH)2", "MnO2", "MnO42-", "MnO4-"],
    Fe: ["Fe", "Fe2+", "Fe3+", "Fe(OH)2", "Fe(OH)3", "Fe3O4"],
    Cu: ["Cu", "Cu+", "Cu2+", "Cu2O", "CuO", "Cu(OH)2"],
    Ag: ["Ag", "Ag+", "Ag2O"],
    Zn: ["Zn", "Zn2+", "ZnO", "Zn(OH)2"],
    Pb: ["Pb", "Pb2+", "PbO", "PbO2"],
  };
  const E_RANGES = {
    Mg: [-2.8, 0.8], Al: [-2.4, 0.8], P: [-2.0, 1.8], S: [-1.6, 2.0],
    Cl: [-0.4, 2.2], I: [-0.4, 1.8], N: [-1.2, 2.0],
    C: [-1.4, 1.8], V: [-2.0, 2.2], Cr: [-1.6, 2.2], Mn: [-2.0, 2.4],
    Fe: [-1.2, 1.8], Cu: [-0.8, 1.6], Ag: [-0.4, 1.6], Zn: [-1.8, 1.2], Pb: [-1.0, 2.0],
  };
  const ENV_GROUPS = [
    { key: "acidbase", name: "酸碱介质", hint: "按物料选择；浓度随 pH 质子化",
      match: c => c.groups.includes("acidbase") },
    { key: "ligand", name: "配体 / 沉淀剂", hint: "与主元素形成配合物或难溶物",
      match: c => c.groups.includes("ligand") || c.groups.includes("precipitant") },
    { key: "redox", name: "氧化还原介质", hint: "用于书写完整反应式与电对线",
      match: c => c.groups.includes("redox") },
  ];
  const MODE_DESC = {
    A: "每点取 Ĝ 最小者；相区内该形态浓度 = 总浓度。相界是直线/曲线的精确分割。",
    B: "浓度高于展示阈值的形态一律叠色（1:1），不按摩尔分数；相界仍保留。",
    C: "按物料平衡求真实分布，按浓度权重混色；固相计浊度、气相按分压着色。",
  };

  function oxOf(sp, main) {
    if (sp.oxByEl && sp.oxByEl[main] != null) return sp.oxByEl[main];
    if (sp.main === main) return sp.ox;
    return sp.ox;
  }
  function findEnv(id) { return D.ENV.find(c => c.id === id); }

  function buildElements() {
    const el = $("elems");
    el.innerHTML = "";
    for (const e of D.MAIN_ELEMENTS) {
      const b = document.createElement("button");
      b.className = "elbtn" + (e === state.main ? " active" : "");
      b.dataset.el = e;
      b.innerHTML = '<span class="sym">' + e + '</span><span class="cn">' + (D.ELEMENT_NAMES[e] || e) + '</span>';
      b.addEventListener("click", () => switchMain(e));
      el.appendChild(b);
    }
  }
  function buildPresets() {
    const el = $("presets");
    el.innerHTML = "";
    for (const pr of D.PRESETS) {
      const b = document.createElement("button");
      b.className = "preset" + (state.activePreset === pr.id ? " active" : "");
      b.textContent = pr.short + " · " + pr.name.split("·")[1].trim();
      b.title = pr.note;
      b.addEventListener("click", () => loadPreset(pr));
      el.appendChild(b);
    }
  }
  function renderFilters() {
    const sps = D.speciesFor(state.main);
    const oxSet = new Set();
    for (const sp of sps) {
      const ox = oxOf(sp, state.main);
      if (ox != null) oxSet.add(ox);
    }
    const oxList = [...oxSet].sort((a, b) => a - b);
    const oxEl = $("ox-filter");
    oxEl.innerHTML = "";
    for (const ox of oxList) {
      const key = String(ox);
      const b = document.createElement("button");
      b.className = "chip" + (state.filters.ox.has(key) ? " active" : "");
      b.textContent = oxLabel(ox);
      b.addEventListener("click", () => {
        if (state.filters.ox.has(key)) state.filters.ox.delete(key); else state.filters.ox.add(key);
        renderFilters(); renderSpeciesList(); updateCounts();
      });
      oxEl.appendChild(b);
    }
    const catEl = $("cat-filter");
    catEl.innerHTML = "";
    for (const cat of D.CATEGORIES) {
      const b = document.createElement("button");
      b.className = "chip cat" + (state.filters.cat.has(cat) ? " active" : "");
      b.textContent = cat;
      b.addEventListener("click", () => {
        if (state.filters.cat.has(cat)) state.filters.cat.delete(cat); else state.filters.cat.add(cat);
        renderFilters(); renderSpeciesList(); updateCounts();
      });
      catEl.appendChild(b);
    }
  }
  function visibleSpecies() {
    let list = D.speciesFor(state.main).slice();
    if (state.filters.ox.size) list = list.filter(sp => state.filters.ox.has(String(oxOf(sp, state.main))));
    if (state.filters.cat.size) list = list.filter(sp => (sp.cats || []).some(c => state.filters.cat.has(c)));
    list.sort((a, b) => {
      const oa = oxOf(a, state.main), ob = oxOf(b, state.main);
      if (oa == null && ob != null) return 1;
      if (oa != null && ob == null) return -1;
      if (oa != ob) return oa - ob;
      const st = s => s === "aq" ? 0 : s === "s" ? 1 : s === "g" ? 2 : 3;
      if (st(a.state) !== st(b.state)) return st(a.state) - st(b.state);
      return a.name.localeCompare(b.name, "zh");
    });
    return list;
  }
  function renderSpeciesList() {
    const el = $("species-list");
    const prevScroll = el.scrollTop;
    el.innerHTML = "";
    const list = visibleSpecies();
    for (const sp of list) {
      const id = sp.id;
      const selected = state.selected.has(id);
      const row = document.createElement("label");
      row.className = "sp-row" + (selected ? " selected" : "");
      if (sp.solvent) row.title = "溶剂，不参与相区比较";
      const cb = document.createElement("input");
      cb.type = "checkbox"; cb.checked = selected;
      if (sp.solvent) cb.disabled = true;
      cb.addEventListener("change", () => toggleSpecies(id, cb.checked));
      const fw = document.createElement("div");
      fw.className = "fwrap";
      const fm = document.createElement("span");
      fm.className = "fm chem";
      fm.innerHTML = formulaFor(sp, state.fmode);
      const nm = document.createElement("span");
      nm.className = "nm"; nm.textContent = sp.name;
      fw.appendChild(fm); fw.appendChild(nm);
      const dg = document.createElement("div");
      dg.className = "dg"; dg.textContent = fmtNum(sp.g, 1); dg.title = "ΔGf° / kJ·mol⁻¹";
      const badges = document.createElement("div");
      badges.className = "badges";
      for (const c of sp.cats || []) {
        const b = document.createElement("span"); b.className = "badge"; b.textContent = c; badges.appendChild(b);
      }
      const st = document.createElement("span"); st.className = "badge state"; st.textContent = stateName(sp.state); badges.appendChild(st);
      if (sp.derived) { const b = document.createElement("span"); b.className = "badge derived"; b.textContent = "推得"; b.title = sp.note || "由 logβ/Ksp/Henry 推得"; badges.appendChild(b); }
      if (sp.solvent) { const b = document.createElement("span"); b.className = "badge"; b.textContent = "溶剂"; badges.appendChild(b); }
      row.appendChild(cb); row.appendChild(fw); row.appendChild(dg); row.appendChild(badges);
      // 自定义浓度
      if (selected && sp.state === "aq") {
        const box = document.createElement("div");
        box.className = "sp-custom";
        const has = state.customC[id] != null;
        const row1 = document.createElement("div");
        row1.className = "row";
        const btn = document.createElement("button");
        btn.textContent = has ? "自定义浓度" : "＋ 自定义浓度";
        btn.addEventListener("click", (ev) => {
          ev.preventDefault(); ev.stopPropagation();
          if (state.customC[id] != null) delete state.customC[id];
          else state.customC[id] = Math.log10(state.Ctot);
          renderSpeciesList(); updateCounts(); scheduleRender();
        });
        const val = document.createElement("span");
        val.className = "v";
        val.textContent = has ? fmtConc(Math.pow(10, state.customC[id])) : "随 C_tot";
        row1.appendChild(btn); row1.appendChild(val);
        if (has) {
          const rm = document.createElement("button");
          rm.textContent = "×"; rm.title = "恢复随总浓度";
          rm.addEventListener("click", (ev) => { ev.preventDefault(); ev.stopPropagation(); delete state.customC[id]; renderSpeciesList(); updateCounts(); scheduleRender(); });
          row1.appendChild(rm);
        }
        box.appendChild(row1);
        if (has) {
          const row2 = document.createElement("div");
          row2.className = "row";
          const r = document.createElement("input");
          r.type = "range"; r.min = "-12"; r.max = "1.7482"; r.step = "0.0001"; r.value = state.customC[id];
          r.addEventListener("input", () => {
            const v = snapLog(+r.value);
            r.value = v;
            state.customC[id] = v;
            val.textContent = fmtConc(Math.pow(10, v));
            scheduleRender(45);
          });
          row2.appendChild(r);
          box.appendChild(row2);
        }
        row.appendChild(box);
      }
      el.appendChild(row);
    }
    el.scrollTop = prevScroll;
  }
  function snapLog(v, snaps) {
    snaps = snaps || [0, -2, -5];
    for (const s of snaps) if (Math.abs(v - s) < 0.06) return s;
    return v;
  }
  function updateCounts() {
    const total = D.speciesFor(state.main).filter(s => !s.solvent).length;
    const sel = [...state.selected].filter(id => byId[id] && !byId[id].solvent).length;
    $("main-count").textContent = total + " 物种";
    $("sp-count").textContent = sel + " / " + total;
    $("env-count").textContent = (Object.keys(state.env).filter(k => state.env[k].on
      && !(k === "Cl-" && state.xAxis === "pCl") && !(k === "NH3" && state.xAxis === "pNH3")).length
      + (state.xAxis !== "pH" ? 1 : 0)) + " 项";
    $("species-note").innerHTML = "勾选带配体/沉淀剂的物种时，会自动启用对应物料。" +
      (Object.keys(state.customC).length ? " 已启用自定义浓度 " + Object.keys(state.customC).length + " 项。" : "");
  }

  function envHidden(c) {
    return c.id === "OH-"
      || (c.id === "Cl-" && state.xAxis === "pCl")
      || (c.id === "NH3" && state.xAxis === "pNH3");
  }
  function envFormulaHTML(comp) {
    if (comp.redox) return formulaById(comp.redox.ox) + ' <span class="muted">/</span> ' + formulaById(comp.redox.red);
    return formulaFor({ formula: comp.formula }, state.fmode);
  }
  function pHControlRow() {
    const row = document.createElement("div");
    row.className = "env-row";
    const top = document.createElement("div"); top.className = "top";
    const lbl = document.createElement("span"); lbl.className = "lbl"; lbl.innerHTML = 'OH⁻ <span class="muted">·</span> pH 控制';
    const tag = document.createElement("span"); tag.className = "tag"; tag.textContent = "固定 pH";
    const v = document.createElement("span"); v.className = "v mono"; v.style.marginLeft = "auto"; v.style.color = "var(--accent2)"; v.style.fontSize = "10px";
    v.textContent = fmtNum(state.pHfix, 1);
    top.appendChild(lbl); top.appendChild(tag); top.appendChild(v); row.appendChild(top);
    const line = document.createElement("div"); line.className = "sliderline";
    const r = document.createElement("input"); r.type = "range"; r.min = "0"; r.max = "14"; r.step = "0.1"; r.value = state.pHfix;
    r.addEventListener("input", () => {
      state.pHfix = +r.value; v.textContent = fmtNum(state.pHfix, 1);
      $("xaxis-note").textContent = "pH 固定 " + fmtNum(state.pHfix, 1);
      scheduleRender(60);
    });
    line.appendChild(r); row.appendChild(line);
    const note = document.createElement("div"); note.className = "note";
    note.textContent = "非 pH 横轴时由 OH⁻ 活度固定 pH：a(OH⁻) = 10^(pH − 14)。";
    row.appendChild(note);
    return row;
  }

  function renderEnvGroups() {
    const el = $("env-groups");
    el.innerHTML = "";
    const seen = new Set();
    let any = false;
    for (const g of ENV_GROUPS) {
      const comps = D.ENV.filter(c => g.match(c) && !envHidden(c));
      const visible = comps.filter(c => state.env[c.id] && state.env[c.id].on && !seen.has(c.id));
      const showPH = g.key === "acidbase" && state.xAxis !== "pH";
      if (!visible.length && !showPH) continue;
      any = true;
      for (const c of visible) seen.add(c.id);
      const det = document.createElement("details");
      det.className = "env-group"; det.open = true;
      const sum = document.createElement("summary");
      sum.innerHTML = '<span class="chev">▶</span><span>' + g.name + '</span><span class="gcount">' + (visible.length + (showPH ? 1 : 0)) + ' 项</span>';
      det.appendChild(sum);
      const body = document.createElement("div");
      body.className = "env-body";
      if (showPH) body.appendChild(pHControlRow());
      for (const comp of visible) body.appendChild(envRow(comp));
      det.appendChild(body);
      el.appendChild(det);
    }
    if (!any) el.innerHTML = '<div class="env-empty">未启用任何物料；点击上方「筛选物料」添加。</div>';
  }

  function envRow(comp) {
    const st = state.env[comp.id];
    const row = document.createElement("div");
    row.className = "env-row";
    const top = document.createElement("div");
    top.className = "top";
    const cb = document.createElement("input"); cb.type = "checkbox"; cb.checked = !!st.on;
    cb.addEventListener("change", () => setEnvOn(comp.id, cb.checked));
    const lbl = document.createElement("span"); lbl.className = "lbl chem"; lbl.innerHTML = envFormulaHTML(comp);
    const rm = document.createElement("button"); rm.className = "rm"; rm.textContent = "×"; rm.title = "关闭";
    rm.addEventListener("click", () => setEnvOn(comp.id, false));
    top.appendChild(cb); top.appendChild(lbl);
    if (!comp.redox) {
      const tag = document.createElement("span"); tag.className = "tag";
      tag.textContent = comp.special === "pH" ? "pH 控制" : "总浓度";
      top.appendChild(tag);
    }
    top.appendChild(rm);
    row.appendChild(top);
    if (comp.note) { const n = document.createElement("div"); n.className = "note"; n.textContent = comp.note; row.appendChild(n); }
    if (comp.special === "pH") {
      const n = document.createElement("div"); n.className = "note"; n.textContent = "a(OH⁻) = 10^(pH − 14)，由 pH 轴决定。"; row.appendChild(n);
    } else {
      const line = document.createElement("div"); line.className = "sliderline";
      const r = document.createElement("input"); r.type = "range"; r.min = "-8"; r.max = "1.7482"; r.step = "0.0001";
      const logC = Math.log10(st.c != null ? st.c : (comp.defaultC || 0.01));
      r.value = logC;
      const v = document.createElement("span"); v.className = "v"; v.textContent = fmtConc(Math.pow(10, logC));
      r.addEventListener("input", () => {
        const lv = snapLog(+r.value);
        r.value = lv;
        state.env[comp.id].c = Math.pow(10, lv);
        v.textContent = fmtConc(state.env[comp.id].c);
        scheduleRender(45);
      });
      line.appendChild(r); line.appendChild(v); row.appendChild(line);
      if (comp.redox) {
        const line2 = document.createElement("div"); line2.className = "sliderline";
        const r2 = document.createElement("input"); r2.type = "range"; r2.min = "-3"; r2.max = "3"; r2.step = "0.1"; r2.value = st.ratio != null ? Math.log10(st.ratio) : 0;
        const v2 = document.createElement("span"); v2.className = "v";
        const ratioLabel = () => { const rr = Math.pow(10, +r2.value); v2.textContent = rr >= 1 ? rr.toFixed(1) + " : 1" : "1 : " + (1 / rr).toFixed(1); };
        ratioLabel();
        r2.addEventListener("input", () => { state.env[comp.id].ratio = Math.pow(10, +r2.value); ratioLabel(); scheduleRender(45); });
        line2.appendChild(r2); line2.appendChild(v2); row.appendChild(line2);
        const ratio = document.createElement("div"); ratio.className = "ratio"; ratio.textContent = "上：氧化态 : 还原态（浓度比）";
        row.appendChild(ratio);
      }
    }
    return row;
  }

  /* ================= 环境弹窗 ================= */
  function openEnvDialog() {
    const dlg = $("env-dialog");
    renderEnvDialog();
    dlg.classList.add("open");
    dlg.setAttribute("aria-hidden", "false");
  }
  function closeEnvDialog() {
    const dlg = $("env-dialog");
    dlg.classList.remove("open");
    dlg.setAttribute("aria-hidden", "true");
    renderEnvGroups();
  }
  function renderEnvDialog() {
    const el = $("env-dialog-body");
    const prev = el.scrollTop;
    el.innerHTML = "";
    for (const g of ENV_GROUPS) {
      const sec = document.createElement("section");
      sec.className = "dlg-group";
      const h = document.createElement("h3"); h.textContent = g.name + " · " + g.hint;
      sec.appendChild(h);
      const comps = D.ENV.filter(c => g.match(c) && !envHidden(c));
      if (g.key === "acidbase" && state.xAxis !== "pH") {
        const item = document.createElement("div");
        item.className = "dlg-item";
        const main = document.createElement("div"); main.className = "main";
        const fm = document.createElement("span"); fm.className = "fm chem"; fm.innerHTML = 'OH⁻ <span class="muted">· pH 控制</span>';
        const note = document.createElement("span"); note.className = "note";
        note.textContent = "非 pH 横轴时用它固定 pH；a(OH⁻) = 10^(pH − 14)。";
        main.appendChild(fm); main.appendChild(document.createElement("br")); main.appendChild(note);
        const line = document.createElement("div"); line.className = "sl";
        const r = document.createElement("input"); r.type = "range"; r.min = "0"; r.max = "14"; r.step = "0.1"; r.value = state.pHfix;
        const v = document.createElement("span"); v.className = "v"; v.textContent = "pH " + fmtNum(state.pHfix, 1);
        r.addEventListener("input", () => {
          state.pHfix = +r.value; v.textContent = "pH " + fmtNum(state.pHfix, 1);
          if ($("xaxis-note")) $("xaxis-note").textContent = "pH 固定 " + fmtNum(state.pHfix, 1);
          scheduleRender(60);
        });
        line.appendChild(r); line.appendChild(v);
        item.appendChild(main); item.appendChild(line);
        sec.appendChild(item);
      }
      for (const comp of comps) {
        const st = state.env[comp.id] || { on: false, c: comp.defaultC, ratio: comp.ratio || 1 };
        const item = document.createElement("div");
        item.className = "dlg-item"; item.dataset.dlgId = comp.id;
        const cb = document.createElement("input"); cb.type = "checkbox"; cb.checked = !!st.on;
        cb.addEventListener("change", () => { setEnvOn(comp.id, cb.checked); renderEnvDialog(); });
        const main = document.createElement("div");
        main.className = "main";
        const fm = document.createElement("span"); fm.className = "fm chem"; fm.innerHTML = envFormulaHTML(comp);
        const note = document.createElement("span"); note.className = "note"; note.textContent = comp.note || "";
        main.appendChild(fm); if (comp.note) { main.appendChild(document.createElement("br")); main.appendChild(note); }
        item.appendChild(cb); item.appendChild(main);
        if (!comp.redox) {
          const tag = document.createElement("span"); tag.className = "tag"; tag.textContent = comp.special === "pH" ? "pH" : "";
          item.appendChild(tag);
        }
        if (st.on && !comp.special) {
          const line = document.createElement("div"); line.className = "sl";
          const r = document.createElement("input"); r.type = "range"; r.min = "-8"; r.max = "1.7482"; r.step = "0.0001";
          r.value = Math.log10(st.c != null ? st.c : (comp.defaultC || 0.01));
          const v = document.createElement("span"); v.className = "v"; v.textContent = fmtConc(Math.pow(10, +r.value));
          r.addEventListener("input", () => {
            const lv = snapLog(+r.value); r.value = lv;
            state.env[comp.id].c = Math.pow(10, lv);
            v.textContent = fmtConc(state.env[comp.id].c);
            scheduleRender(45);
          });
          line.appendChild(r); line.appendChild(v); item.appendChild(line);
          if (comp.redox) {
            const line2 = document.createElement("div"); line2.className = "sl";
            const r2 = document.createElement("input"); r2.type = "range"; r2.min = "-3"; r2.max = "3"; r2.step = "0.1";
            r2.value = Math.log10(st.ratio || 1);
            const v2 = document.createElement("span"); v2.className = "ratio-v";
            const lab = () => { const rr = Math.pow(10, +r2.value); v2.textContent = "氧化:还原 " + (rr >= 1 ? rr.toFixed(1) + ":1" : "1:" + (1 / rr).toFixed(1)); };
            lab();
            r2.addEventListener("input", () => { state.env[comp.id].ratio = Math.pow(10, +r2.value); lab(); scheduleRender(45); });
            line2.appendChild(r2); line2.appendChild(v2); item.appendChild(line2);
          }
        }
        sec.appendChild(item);
      }
      el.appendChild(sec);
    }
    el.scrollTop = prev;
  }

  /* ================= 选择 / 环境联动 ================= */
  let syncing = false;
  function enableLigandDirect(lid) {
    if (lid === "OH-") return;
    const comp = findEnv(lid);
    if (!comp) return;
    if (state.env[lid]) { state.env[lid].on = true; return; }
    state.env[lid] = { on: true, c: comp.defaultC != null ? comp.defaultC : 0.01, ratio: comp.ratio != null ? comp.ratio : 1 };
  }
  function syncLigandsForSelection() {
    if (syncing) return;
    syncing = true;
    for (const id of state.selected) {
      const sp = byId[id];
      if (!sp) continue;
      for (const lid of sp.ligs || []) enableLigandDirect(lid);
    }
    syncing = false;
  }
  function maybeDisableLigand(lid) {
    if (lid === "OH-" || !state.env[lid] || !state.env[lid].on) return;
    for (const id of state.selected) {
      const sp = byId[id];
      if (sp && sp.ligs && sp.ligs.includes(lid)) return;
    }
    const related = D.speciesFor(state.main).filter(sp => sp.ligs && sp.ligs.includes(lid) && sp.cats.includes("难溶物"));
    if (related.length && !related.some(sp => state.selected.has(sp.id))) {
      state.env[lid].on = false;
    }
  }
  function setEnvOn(id, on, opts) {
    opts = opts || {};
    const comp = findEnv(id);
    if (!comp) return;
    if (on) {
      if (!state.env[id]) state.env[id] = { on: true, c: comp.defaultC != null ? comp.defaultC : 0.01, ratio: comp.ratio != null ? comp.ratio : 1 };
      else state.env[id].on = true;
      if (!opts.noAutoSelect && (comp.groups.includes("ligand") || comp.groups.includes("precipitant"))) {
        const related = D.speciesFor(state.main).filter(sp => sp.ligs && sp.ligs.includes(id) && sp.cats.includes("难溶物"));
        for (const sp of related) state.selected.add(sp.id);
      }
    } else {
      if (state.env[id]) state.env[id].on = false;
      for (const sid of [...state.selected]) {
        const sp = byId[sid];
        if (sp && sp.ligs && sp.ligs.includes(id)) { state.selected.delete(sid); delete state.customC[sid]; }
      }
    }
    if (!opts.silent) {
      syncLigandsForSelection();
      refreshStructural();
      if ($("env-dialog").classList.contains("open")) renderEnvDialog();
    }
  }
  function isHydroxide(sp) {
    return sp.state === "s" && sp.cats.includes("难溶物") && sp.name.indexOf("氢氧化") === 0;
  }
  function isChlorideSpecies(sp) {
    return (sp.el.Cl || 0) > 0 && (sp.state === "s" || sp.cats.includes("络合物"));
  }
  function autoSelectForAxis(axis) {
    let n = 0;
    for (const sp of D.speciesFor(state.main)) {
      if (sp.solvent) continue;
      const hit = axis === "pH" ? isHydroxide(sp)
        : axis === "pCl" ? isChlorideSpecies(sp)
        : axis === "pNH3" ? sp.ligs.includes("NH3") : false;
      if (hit && !state.selected.has(sp.id)) { state.selected.add(sp.id); n++; }
    }
    if (n) syncLigandsForSelection();
    return n;
  }
  function switchAxis(next) {
    if (next === state.xAxis) return;
    if (next === "pCl" && state.main === "Cl") return;
    if (next === "pNH3" && state.main === "N") return;
    if (state.xAxis === "pH") state._pHrange = [state.xmin, state.xmax];
    if (next === "pCl" && state._pClrange == null) state._pClrange = [-4, 1.5];
    if (next === "pNH3" && state._pNH3range == null) state._pNH3range = [-6, 1];
    state.xAxis = next;
    if (next === "pH") {
      state.xmin = state._pHrange ? state._pHrange[0] : 0;
      state.xmax = state._pHrange ? state._pHrange[1] : 14;
    } else if (next === "pCl") {
      state.xmin = state._pClrange[0]; state.xmax = state._pClrange[1];
    } else {
      state.xmin = state._pNH3range[0]; state.xmax = state._pNH3range[1];
    }
    autoSelectForAxis(next);
    state._polyKey = null; state.fieldKey = null;
    refreshStructural();
  }
  function toggleSpecies(id, on) {
    const sp = byId[id];
    if (!sp) return;
    if (on) {
      state.selected.add(id);
      for (const lid of sp.ligs || []) enableLigandDirect(lid);
    } else {
      state.selected.delete(id);
      delete state.customC[id];
      for (const lid of sp.ligs || []) maybeDisableLigand(lid);
    }
    syncLigandsForSelection();
    refreshStructural();
  }

  /* ================= 预设 / 主元素 ================= */
  function loadPreset(pr) {
    state.main = pr.main;
    state.activePreset = pr.id;
    state.selected = new Set(pr.selected.filter(id => byId[id] && byId[id].el[pr.main] > 0 && !byId[id].solvent));
    state.customC = {};
    state.filters = { ox: new Set(), cat: new Set() };
    state.env = {};
    for (const id of pr.env || []) {
      const comp = findEnv(id);
      if (!comp) continue;
      state.env[id] = { on: true, c: comp.defaultC != null ? comp.defaultC : 0.01, ratio: comp.ratio != null ? comp.ratio : 1 };
    }
    state.Ctot = pr.Ctot;
    state.mode = "C";
    state.xAxis = "pH"; state.pHfix = 7;
    state.xmin = pr.xmin != null ? pr.xmin : 0;
    state.xmax = pr.xmax != null ? pr.xmax : 14;
    state.emin = pr.emin != null ? pr.emin : (E_RANGES[pr.main] || [-1.2, 1.8])[0];
    state.emax = pr.emax != null ? pr.emax : (E_RANGES[pr.main] || [-1.2, 1.8])[1];
    state.customColors = {};
    state.randomMap = {};
    state.sysKey = null; state._polyKey = null; state.fieldKey = null;
    state.polygons = []; state.boundaries = []; state.boundaryChains = []; state.labels = []; state.groups = null;
    autoSelectForAxis("pH");
    syncLigandsForSelection();
    buildElements(); buildPresets(); refreshStructural();
    scheduleRender(0);
  }
  function switchMain(el) {
    if (el === state.main) return;
    state.main = el;
    state.activePreset = null;
    state.selected = new Set((DEFAULT_SELECTED[el] || []).filter(id => byId[id] && !byId[id].solvent));
    state.customC = {};
    state.filters = { ox: new Set(), cat: new Set() };
    const range = E_RANGES[el] || [-1.2, 1.8];
    state.emin = range[0]; state.emax = range[1];
    state.xmin = 0; state.xmax = 14; state.xAxis = "pH"; state.pHfix = 7;
    state.sysKey = null; state._polyKey = null; state.fieldKey = null;
    state.polygons = []; state.boundaries = []; state.boundaryChains = []; state.labels = []; state.groups = null;
    autoSelectForAxis("pH");
    syncLigandsForSelection();
    buildElements(); buildPresets(); refreshStructural();
    scheduleRender(0);
  }
  function autoScale() {
    const r = E_RANGES[state.main] || [-1.2, 1.8];
    state.emin = r[0]; state.emax = r[1];
    $("emin").value = state.emin.toFixed(1); $("emax").value = state.emax.toFixed(1);
    scheduleRender(0);
  }

  /* ================= 控件同步 ================= */
  function syncControls() {
    $("ctot").value = Math.log10(state.Ctot);
    $("ctotval").textContent = fmtConc(state.Ctot);
    $("penv").value = Math.log10(state.Penv);
    $("penvval").textContent = fmtP(state.Penv);
    $("v0").value = state.V0;
    $("v0val").textContent = state.V0.toFixed(1) + " L";
    $("cdisp").value = Math.log10(state.Cdisp);
    $("cdispval").textContent = fmtConc(state.Cdisp);
    for (const b of $$("#xaxis-btns button")) {
      const axis = b.dataset.axis;
      b.classList.toggle("active", axis === state.xAxis);
      b.disabled = (axis === "pCl" && state.main === "Cl") || (axis === "pNH3" && state.main === "N");
      b.title = axis === "pCl" && state.main === "Cl" ? "主元素为 Cl 时，Cl⁻ 浓度不是独立横轴"
        : axis === "pNH3" && state.main === "N" ? "主元素为 N 时，NH₃ 浓度不是独立横轴" : "";
    }
    $("xaxis-note").textContent = state.xAxis === "pH" ? "" : "pH 固定 " + fmtNum(state.pHfix, 1);
    $("fmode").value = state.fmode;
    $("fplain").value = state.fplain;
    $("flegend").value = state.flegend;
    $("couple").checked = state.couple;
    $("show-boundary").checked = state.showBoundary;
    $("show-water").checked = state.showWater;
    $("show-medium").checked = state.showMedium;
    $("show-all-redox").checked = state.showAllRedox;
    $("show-labels").checked = state.showLabels;
    $("show-grid").checked = state.showGrid;
    $("emin").value = state.emin.toFixed(1);
    $("emax").value = state.emax.toFixed(1);
    for (const b of $$("#mode-btns button")) b.classList.toggle("active", b.dataset.mode === state.mode);
    for (const b of $$("#palette-btns button")) b.classList.toggle("active", b.dataset.pal === state.palMode);
    $("mode-desc").textContent = MODE_DESC[state.mode] || "";
    $("ctot-note").innerHTML = "滑条在 1 M / 0.01 M / 10⁻⁵ M 处磁吸。" +
      (Object.keys(state.customC).length
        ? " 已启用 " + Object.keys(state.customC).length + " 个物种自定义浓度。"
        : " 所有物种共用该分析浓度。");
  }
  function refreshStructural() {
    renderFilters();
    renderSpeciesList();
    renderEnvGroups();
    updateCounts();
    syncControls();
    scheduleRender();
  }

  /* ================= 事件绑定 ================= */
  function bindEvents() {
    $("ctot").addEventListener("input", () => {
      const lv = snapLog(+$("ctot").value);
      $("ctot").value = lv;
      state.Ctot = Math.pow(10, lv);
      $("ctotval").textContent = fmtConc(state.Ctot);
      $("ctot-note").innerHTML = "滑条在 1 M / 0.01 M / 10⁻⁵ M 处磁吸。" +
        (Object.keys(state.customC).length ? " 已启用 " + Object.keys(state.customC).length + " 个物种自定义浓度。" : " 所有物种共用该分析浓度。");
      scheduleRender(60);
    });
    $("sp-select-all").addEventListener("click", (e) => {
      e.preventDefault();
      for (const sp of visibleSpecies()) {
        if (sp.solvent) continue;
        state.selected.add(sp.id);
        for (const lid of sp.ligs || []) enableLigandDirect(lid);
      }
      syncLigandsForSelection(); refreshStructural();
    });
    $("sp-clear").addEventListener("click", (e) => {
      e.preventDefault();
      state.selected.clear(); state.customC = {};
      refreshStructural();
    });
    $("env-manage").addEventListener("click", () => openEnvDialog());
    $("env-dialog-close").addEventListener("click", closeEnvDialog);
    $("env-dialog-done").addEventListener("click", closeEnvDialog);
    const backdrop = document.querySelector(".modal-backdrop");
    if (backdrop) backdrop.addEventListener("click", closeEnvDialog);
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape" && $("env-dialog").classList.contains("open")) closeEnvDialog();
    });
    $("env-dialog-clear").addEventListener("click", () => {
      for (const k in state.env) state.env[k].on = false;
      state.selected.clear(); state.customC = {};
      refreshStructural(); renderEnvDialog();
    });
    $("env-dialog-common").addEventListener("click", () => {
      for (const id of ["Cl-", "CO32-", "NH3"]) {
        const comp = findEnv(id);
        if (comp) setEnvOn(id, true, { noAutoSelect: true, silent: true });
      }
      syncLigandsForSelection();
      refreshStructural(); renderEnvDialog();
    });
    $("penv").addEventListener("input", () => {
      state.Penv = Math.pow(10, +$("penv").value);
      $("penvval").textContent = fmtP(state.Penv);
      scheduleRender(60);
    });
    $("v0").addEventListener("input", () => {
      state.V0 = +$("v0").value;
      $("v0val").textContent = state.V0.toFixed(1) + " L";
      scheduleRender(60);
    });
    $("couple").addEventListener("change", () => { state.couple = $("couple").checked; scheduleRender(0); });
    for (const b of $$("#mode-btns button")) {
      b.addEventListener("click", () => {
        state.mode = b.dataset.mode;
        for (const x of $$("#mode-btns button")) x.classList.toggle("active", x === b);
        $("mode-desc").textContent = MODE_DESC[state.mode] || "";
        scheduleRender(0);
      });
    }
    $("cdisp").addEventListener("input", () => {
      state.Cdisp = Math.pow(10, +$("cdisp").value);
      $("cdispval").textContent = fmtConc(state.Cdisp);
      scheduleRender(80);
    });
    for (const b of $$("#xaxis-btns button")) {
      b.addEventListener("click", () => switchAxis(b.dataset.axis));
    }
    for (const b of $$("#palette-btns button")) {
      b.addEventListener("click", () => {
        state.palMode = b.dataset.pal;
        for (const x of $$("#palette-btns button")) x.classList.toggle("active", x === b);
        if (state.palMode === "random") {
          state.palSeed = (state.palSeed + 1) % 997;
          recolorRandom();
        }
        scheduleRender(0);
      });
    }
    $("recolor").addEventListener("click", () => {
      state.palSeed = (state.palSeed + 1) % 997;
      state.palMode = "random";
      for (const x of $$("#palette-btns button")) x.classList.toggle("active", x.dataset.pal === "random");
      recolorRandom();
      scheduleRender(0);
    });
    $("clear-colors").addEventListener("click", () => { state.customColors = {}; scheduleRender(0); });
    $("fmode").addEventListener("change", () => { state.fmode = $("fmode").value; renderSpeciesList(); renderEnvGroups(); renderLegend(); renderReactions(); });
    $("fplain").addEventListener("change", () => { state.fplain = $("fplain").value; renderSpeciesList(); renderLegend(); renderReactions(); });
    $("flegend").addEventListener("change", () => { state.flegend = $("flegend").value; renderLegend(); });
    const toggleMap = {
      "show-boundary": "showBoundary", "show-water": "showWater", "show-medium": "showMedium",
      "show-all-redox": "showAllRedox", "show-labels": "showLabels", "show-grid": "showGrid",
    };
    for (const id in toggleMap) {
      $(id).addEventListener("change", () => { state[toggleMap[id]] = $(id).checked; scheduleRender(0); });
    }
    $("emin").addEventListener("change", () => {
      const v = +$("emin").value; if (isFinite(v) && v < state.emax) { state.emin = v; state._polyKey = null; scheduleRender(0); }
    });
    $("emax").addEventListener("change", () => {
      const v = +$("emax").value; if (isFinite(v) && v > state.emin) { state.emax = v; state._polyKey = null; scheduleRender(0); }
    });
    $("autoscale").addEventListener("click", autoScale);
    $("reset-all").addEventListener("click", () => {
      const pr = D.PRESETS.find(p => p.main === state.main);
      if (pr) loadPreset(pr);
      else { switchMain(state.main); }
    });
    $("reaction-toggle").addEventListener("click", () => {
      const el = $("reactions");
      const exp = el.classList.toggle("expanded");
      $("reaction-toggle").textContent = exp ? "收起" : "展开";
    });
    $("save-png").addEventListener("click", savePNG);
    plot.addEventListener("mousemove", onMove);
    plot.addEventListener("mouseleave", onLeave);
    window.addEventListener("resize", debounce(() => {
      const oldW = geom ? geom.W : 0;
      setupCanvas();
      if (Math.abs(geom.W - oldW) > 1) { state._polyKey = null; scheduleRender(0); }
    }, 120));
  }
  function debounce(fn, ms) {
    let t = null;
    return function () { if (t) clearTimeout(t); t = setTimeout(fn, ms); };
  }

  /* ================= 保存 PNG ================= */
  function savePNG() {
    if (!state.sys) return;
    const dpr = geom.dpr;
    const cssW = geom.W, cssH = geom.H;
    const margin = 34;
    const ids = [...state.selected].filter(id => byId[id] && byId[id].el[state.main] > 0 && !byId[id].solvent).slice(0, 18);
    const legendRows = Math.max(1, Math.ceil(ids.length / 3));
    const outW = Math.round((cssW + margin * 2) * dpr);
    const outH = Math.round((66 + cssH + 22 + legendRows * 21 + 34) * dpr);
    const out = document.createElement("canvas");
    out.width = outW; out.height = outH;
    const c = out.getContext("2d");
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.fillStyle = "#ffffff";
    c.fillRect(0, 0, cssW + margin * 2, outH / dpr);
    c.fillStyle = "#342c40";
    c.font = '650 22px "Segoe UI","Microsoft YaHei",sans-serif';
    c.textAlign = "left"; c.textBaseline = "alphabetic";
    c.fillText(state.main + " Pourbaix 图", margin, 31);
    c.fillStyle = "#8b7f9c";
    c.font = '11px "Segoe UI","Microsoft YaHei",sans-serif';
    c.fillText("25 °C · 理想溶液 · C_tot = " + fmtConc(state.Ctot) + " · 模式 " + state.mode + " · E / V vs SHE", margin, 50);
    c.drawImage(plot, margin, 66, cssW, cssH);
    // 图例
    let lx = margin, ly = 66 + cssH + 22;
    c.font = '11px "Cambria Math","Segoe UI",sans-serif';
    for (const id of ids) {
      const sp = byId[id];
      const col = colorFor(id);
      c.fillStyle = col.fill;
      c.fillRect(lx, ly - 9, 11, 11);
      c.fillStyle = "#4a4158";
      const txt = chemUni(sp.formula) + (sp.state === "aq" ? "" : "(" + sp.state + ")");
      c.fillText(txt, lx + 15, ly);
      lx += 15 + c.measureText(txt).width + 18;
      if (lx > cssW - 60) { lx = margin; ly += 21; }
    }
    c.fillStyle = "#8b7f9c";
    c.font = '10px "JetBrains Mono",Consolas,monospace';
    c.fillText("数据：NBS 1982 / 推得 · 反应式由 ΔGf° 现场配平", margin, outH / dpr - 12);
    const a = document.createElement("a");
    a.download = "pourbaix-" + state.main + "-" + state.mode + ".png";
    a.href = out.toDataURL("image/png");
    a.click();
  }

  /* ================= 初始化 ================= */
  function init() {
    applyTheme();
    buildElements();
    buildPresets();
    bindEvents();
    const params = new URLSearchParams(location.search);
    const presetId = params.get("preset");
    const main = params.get("main");
    const mode = params.get("mode");
    if (presetId) {
      const pr = D.PRESETS.find(p => p.id === presetId);
      if (pr) loadPreset(pr);
      else { state.main = main || "Fe"; state.selected = new Set(DEFAULT_SELECTED[state.main] || []); syncLigandsForSelection(); refreshStructural(); }
    } else if (main && D.MAIN_ELEMENTS.includes(main)) {
      state.main = main;
      state.selected = new Set(DEFAULT_SELECTED[main] || []);
      state.emin = (E_RANGES[main] || [-1.2, 1.8])[0];
      state.emax = (E_RANGES[main] || [-1.2, 1.8])[1];
      state.xAxis = "pH";
      autoSelectForAxis("pH");
      buildElements(); buildPresets();
      syncLigandsForSelection(); refreshStructural();
    } else {
      loadPreset(D.PRESETS[0]);
    }
    if (mode === "A" || mode === "B" || mode === "C") state.mode = mode;
    const envParam = params.get("env");
    if (envParam) {
      for (const token of envParam.split(",")) {
        const parts = token.split(":");
        const id = parts[0], logC = parts[1];
        const comp = findEnv(id);
        if (!comp) continue;
        state.env[id] = {
          on: true,
          c: logC != null ? Math.pow(10, +logC) : (comp.defaultC != null ? comp.defaultC : 0.01),
          ratio: comp.ratio != null ? comp.ratio : 1,
        };
      }
    }
    const selParam = params.get("sel");
    if (selParam) {
      const ids = selParam.replace(/ /g, "+").split(",").map(s => s.trim()).filter(id => byId[id] && byId[id].el[state.main] > 0);
      if (ids.length) state.selected = new Set(ids);
      state.customC = {};
    }
    const ctotParam = params.get("ctot");
    if (ctotParam != null) {
      const v = +ctotParam;
      if (isFinite(v)) state.Ctot = Math.pow(10, Math.max(-12, Math.min(1.7482, v)));
    }
    const xaxisParam = params.get("xaxis");
    if (xaxisParam === "pCl" && state.main !== "Cl") {
      state.xAxis = "pCl"; state.xmin = -4; state.xmax = 1.5;
      autoSelectForAxis("pCl");
    } else if (xaxisParam === "pNH3" && state.main !== "N") {
      state.xAxis = "pNH3"; state.xmin = -6; state.xmax = 1;
      autoSelectForAxis("pNH3");
    } else if (xaxisParam === "pH") {
      state.xAxis = "pH"; state.xmin = 0; state.xmax = 14;
    }
    const phfixParam = params.get("phfix");
    if (phfixParam != null && isFinite(+phfixParam)) state.pHfix = Math.max(0, Math.min(14, +phfixParam));
    const palParam = params.get("pal");
    if (["intrinsic", "random"].includes(palParam)) state.palMode = palParam;
    const fmodeParam = params.get("fmode");
    if (["html", "unicode", "plain"].includes(fmodeParam)) state.fmode = fmodeParam;
    const fplainParam = params.get("fplain");
    if (["compact", "spaced", "paren"].includes(fplainParam)) state.fplain = fplainParam;
    const flegendParam = params.get("flegend");
    if (["html", "unicode", "plain"].includes(flegendParam)) state.flegend = flegendParam;
    syncLigandsForSelection();
    refreshStructural();
    applyTheme();
    syncControls();
    setupCanvas();
    renderAll();
    if (params.get("scroll") === "display") {
      const sec = $("sec-display");
      if (sec) sec.scrollIntoView({ block: "start" });
    }
    // 首次绘制后再校正一次尺寸（字体/滚动条可能改变宽度）
    setTimeout(() => { setupCanvas(); state._polyKey = null; renderAll(); }, 180);
    if (params.get("dialog") === "1") openEnvDialog();
    const setenvParam = params.get("setenv");
    if (setenvParam) {
      for (const id of setenvParam.split(",")) setEnvOn(id, true);
      renderAll();
    }
    const hoverParam = params.get("hover");
    if (hoverParam) {
      const hv = hoverParam.split(",").map(Number);
      if (hv.length === 2 && isFinite(hv[0]) && isFinite(hv[1])) hoverAt(hv[0], hv[1]);
    }
  }

  // 供调试/自动化使用
  window.PBApp = { state, loadPreset, switchMain, renderAll, setupCanvas, savePNG, setEnvOn, openEnvDialog, D, C };

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();