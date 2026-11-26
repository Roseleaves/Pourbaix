/* ============================================================
   chem.js — 纯函数计算层（无 DOM，可在 Node 中单测）
   ------------------------------------------------------------
   参考态分解：物种 → 主元素标准态 + H⁺ + H₂O + e⁻ + 配体参考态。
   分解只负责“配平”，物种自身的 ΔGf° 负责“能量”，二者分离，
   因此不会把不同氧化态的配体混进酸碱反应。
   Ĝ_i(pH,E) = [ΔGf°_i − Σ ν_b μ_b(pH)] / n_i + Q_i E
   Q_i = ν_e F / n_i
   ============================================================ */
var PBChem = (function (D) {
  "use strict";
  const R = 8.314462618;                 // J·mol⁻¹·K⁻¹
  const T = 298.15;
  const F = 96485.33212;                 // C·mol⁻¹
  const F_KJ = F / 1000;                 // 96.485 kJ·V⁻¹·mol⁻¹
  const RTln10 = D.RTln10;               // 5.708 kJ·mol⁻¹
  const NERNST = RTln10 / F_KJ;          // 0.05916 V
  const RGAS_LBAR = 0.0831446;           // L·bar·mol⁻¹·K⁻¹
  const EPS = 1e-9;
  const TOL = 1e-5;

  const byId = D.byId;
  const GAS_PAIR = Object.create(null);
  for (const sp of D.DB) if (sp.aqPair) GAS_PAIR[sp.aqPair] = sp.id;

  /* ---------------- 线性代数 ---------------- */
  function solveLin(A, b) {
    const n = b.length;
    const M = A.map((r, i) => r.slice().concat([b[i]]));
    for (let c = 0; c < n; c++) {
      let p = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      if (Math.abs(M[p][c]) < 1e-10) return null;
      if (p !== c) { const t = M[c]; M[c] = M[p]; M[p] = t; }
      const piv = M[c][c];
      for (let r = 0; r < n; r++) {
        if (r === c) continue;
        const f = M[r][c] / piv;
        if (Math.abs(f) < 1e-15) continue;
        for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
      }
    }
    return M.map((r, i) => r[n] / M[i][i]);
  }

  function solveLS(A, b) {
    const m = A.length, n = A[0].length;
    const AtA = Array.from({ length: n }, () => new Array(n).fill(0));
    const Atb = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        let s = 0;
        for (let k = 0; k < m; k++) s += A[k][i] * A[k][j];
        AtA[i][j] = AtA[j][i] = s;
      }
      let s = 0;
      for (let k = 0; k < m; k++) s += A[k][i] * b[k];
      Atb[i] = s;
    }
    const x = solveLin(AtA, Atb);
    if (!x) return null;
    let maxRes = 0;
    for (let k = 0; k < m; k++) {
      let r = -b[k];
      for (let j = 0; j < n; j++) r += A[k][j] * x[j];
      maxRes = Math.max(maxRes, Math.abs(r));
    }
    return { x, maxRes };
  }

  function rankCols(cols, m) {
    const n = cols.length;
    if (!n) return 0;
    const M = Array.from({ length: m }, (_, i) => cols.map(c => c[i]));
    let rank = 0;
    for (let c = 0; c < n && rank < m; c++) {
      let p = rank;
      for (let i = rank + 1; i < m; i++) if (Math.abs(M[i][c]) > Math.abs(M[p][c])) p = i;
      if (Math.abs(M[p][c]) < 1e-9) continue;
      if (p !== rank) { const t = M[rank]; M[rank] = M[p]; M[p] = t; }
      const piv = M[rank][c];
      for (let i = 0; i < m; i++) {
        if (i === rank) continue;
        const f = M[i][c] / piv;
        if (Math.abs(f) < 1e-15) continue;
        for (let k = c; k < n; k++) M[i][k] -= f * M[rank][k];
      }
      rank++;
    }
    return rank;
  }

  function gcd(a, b) { a = Math.abs(Math.round(a)); b = Math.abs(Math.round(b)); while (b) { const t = a % b; a = b; b = t; } return a || 1; }

  function integerize(coeffs) {
    for (let L = 1; L <= 60; L++) {
      if (coeffs.every(c => Math.abs(c * L - Math.round(c * L)) < 1e-5)) {
        return { L, vals: coeffs.map(c => Math.round(c * L)) };
      }
    }
    return { L: 1, vals: coeffs.map(c => Math.round(c)) };
  }

  /* ---------------- 酸碱家族 ---------------- */
  function familyInfo(key) {
    const fam = D.FAMILIES[key];
    if (!fam) return null;
    const list = fam.list;
    const pKa = [];
    for (let j = 1; j < list.length; j++) {
      let dG = byId[list[j]].g - byId[list[j - 1]].g;
      if (j === 1 && fam.hydrate) dG -= byId["H2O"].g;
      pKa.push(dG / RTln10);
    }
    return { key, list, pKa, hydrate: !!fam.hydrate, ref: fam.ref };
  }

  // α_j(pH)：家族中第 j 个形态占总分析浓度的分数
  function familyFractions(key, pH) {
    const info = familyInfo(key);
    if (!info) return null;
    const { list, pKa } = info;
    const n = list.length - 1;
    const logRatio = new Array(list.length).fill(0); // [A_j]/[A_n]
    for (let j = n - 1; j >= 0; j--) {
      logRatio[j] = logRatio[j + 1] + (pKa[j] - pH);
    }
    let sum = 0;
    const raw = logRatio.map(lr => Math.pow(10, Math.max(-300, Math.min(300, lr))));
    // 用最大值缩放，避免溢出
    const mx = Math.max(...raw);
    for (let j = 0; j < raw.length; j++) { raw[j] /= mx; sum += raw[j]; }
    return raw.map(v => v / sum);
  }

  function ligandAlpha(refId, pH) {
    const key = D.LIG_FAMILY[refId];
    if (!key) return 1;
    if (refId === "OH-") return Math.pow(10, pH - 14);
    const info = familyInfo(key);
    if (!info) return 1;
    const frac = familyFractions(key, pH);
    const idx = info.list.indexOf(refId);
    if (idx < 0) return 0;
    return frac[idx];
  }

  /* ---------------- 体系构建 ---------------- */
  function normalizeEnv(env) {
    const out = {};
    if (!env) return out;
    for (const id in env) {
      const e = env[id];
      if (e === true) out[id] = { on: true, c: null, ratio: 1 };
      else out[id] = { on: !!e.on, c: e.c != null ? e.c : null, ratio: e.ratio != null ? e.ratio : 1 };
    }
    return out;
  }

  function createSystem(cfg) {
    const main = cfg.main;
    const selected = (cfg.selected || []).filter(id => byId[id]);
    const env = normalizeEnv(cfg.env);
    const Ctot = Math.max(cfg.Ctot != null ? cfg.Ctot : 1e-3, 1e-30);
    const Penv = Math.max(cfg.Penv != null ? cfg.Penv : 1, 1e-6);
    const xAxis = cfg.xAxis || "pH";
    const pHfix = cfg.pHfix != null ? cfg.pHfix : 7;
    const couple = cfg.couple !== false;
    const customC = cfg.customC || {};
    const customP = cfg.customP || {};

    // 当前主元素下的物种（排除溶剂 H2O）
    const allSel = selected
      .map(id => byId[id])
      .filter(sp => sp && sp.el[main] > 0 && !sp.solvent);
    const solventIds = selected.filter(id => byId[id] && byId[id].solvent);

    // 分析浓度与元素集合
    const elements = new Set([main, "H", "O"]);
    for (const sp of allSel) for (const e in sp.el) elements.add(e);
    const elementList = [...elements];
    const rows = elementList.concat(["z"]);
    const mRows = rows.length;

    function vecOf(sp) {
      return rows.map(r => r === "z" ? sp.z : (sp.el[r] || 0));
    }

    // 配体参考态 → 环境物料
    const refToEnv = {};
    for (const comp of D.ENV) {
      let refId = comp.id;
      if (comp.family && D.FAMILIES[comp.family]) refId = D.FAMILIES[comp.family].ref;
      if (!refToEnv[refId]) refToEnv[refId] = comp;
      if (!refToEnv[comp.id]) refToEnv[comp.id] = comp;
    }

    // 从环境读取配体总浓度
    function envFor(refId) { return refToEnv[refId] || null; }

    function ligandActivity(refId, pH) {
      if (refId === "OH-") return Math.pow(10, pH - 14);
      const comp = envFor(refId);
      let c = null, on = false;
      if (comp) {
        const st = env[comp.id];
        if (st && st.on) { on = true; c = st.c != null ? st.c : comp.defaultC; }
      }
      if (c == null) c = comp && comp.defaultC != null ? comp.defaultC : 1e-3;
      c = Math.max(c, 1e-30);
      if (comp && comp.family) {
        const key = comp.family === "water" ? "water" : comp.family;
        const frac = familyFractions(key, pH);
        const info = familyInfo(key);
        if (info && frac) {
          const idx = info.list.indexOf(refId);
          if (idx >= 0) return Math.max(c * frac[idx], 1e-300);
        }
      }
      return c;
    }

    // 配体活度：如果环境里没有该物料，仍按 1e-3 M 处理，避免“关掉配体就无解”
    function ligandMu(refId, pH) {
      const sp = byId[refId];
      if (!sp) return 0;
      const a = refId === "OH-" ? Math.pow(10, pH - 14) : ligandActivity(refId, pH);
      return sp.g + RTln10 * Math.log10(Math.max(a, 1e-300));
    }

    // 基组
    const basis = [];
    const basisIds = new Set();
    const basisCols = [];

    function addBasis(comp, kind, meta) {
      if (!comp || basisIds.has(comp.id)) return false;
      const col = vecOf(comp);
      const before = rankCols(basisCols, mRows);
      const after = rankCols(basisCols.concat([col]), mRows);
      if (after <= before) return false;
      basis.push(Object.assign({ kind }, meta || {}, {
        id: comp.id, el: comp.el, z: comp.z, g: comp.g, state: comp.state, comp,
      }));
      basisIds.add(comp.id);
      basisCols.push(col);
      return true;
    }

    // 固定基准
    const mainRef = { id: "@" + main, el: { [main]: 1 }, z: 0, g: 0, state: "ref", formula: main };
    addBasis(mainRef, "main");
    addBasis(byId["H+"], "proton");
    addBasis(byId["H2O"], "water");
    addBasis(byId["e-"], "electron");

    // 候选配体/参考态
    const candidates = [];
    const candIds = new Set();
    function pushCand(sp, kind, meta) {
      if (!sp || candIds.has(sp.id)) return;
      if (sp.el[main] > 0) return;               // 主元素物种不能当配体参考
      candIds.add(sp.id);
      candidates.push(Object.assign({ kind, comp: sp }, meta || {}));
    }

    for (const comp of D.ENV) {
      const st = env[comp.id];
      if (!st || !st.on) continue;
      let refId = comp.id;
      if (comp.family && D.FAMILIES[comp.family]) refId = D.FAMILIES[comp.family].ref;
      if (refId === "OH-" || refId === "H2O" || refId === "H+" || refId === "e-") continue;
      const sp = byId[refId];
      if (sp) pushCand(sp, "ligand", { envId: comp.id, family: comp.family || null, refId });
    }
    for (const sp of allSel) {
      for (const lid of (sp.ligs || [])) {
        if (lid === "OH-" || lid === "H2O" || lid === "H+" || lid === "e-") continue;
        const L = byId[lid];
        if (!L) continue;
        const comp = envFor(lid);
        pushCand(L, "ligand", { envId: comp ? comp.id : null, family: comp ? comp.family : null, refId: lid });
      }
    }
    // 额外元素默认参考态
    for (const e of elementList) {
      if (e === main || e === "H" || e === "O") continue;
      const refId = D.DEFAULT_REF[e];
      const sp = refId ? byId[refId] : null;
      if (sp && sp.el[main] === 0) pushCand(sp, "default", { element: e });
    }
    // 元素标准态兜底（最后才用）
    for (const e of elementList) {
      if (e === main || e === "H" || e === "O") continue;
      const id = "@el:" + e;
      if (candIds.has(id)) continue;
      candIds.add(id);
      candidates.push({
        kind: "element", comp: { id, el: { [e]: 1 }, z: 0, g: 0, state: "ref", formula: e },
        element: e,
      });
    }

    for (const cand of candidates) {
      if (basisIds.has(cand.comp.id)) continue;
      if (cand.comp.el[main] > 0) continue;
      addBasis(cand.comp, cand.kind, cand);
    }

    // 分解
    const A = [];
    for (let k = 0; k < mRows; k++) A.push(basis.map(b => rows[k] === "z" ? b.z : (b.el[rows[k]] || 0)));
    const decomp = {};
    const unbalanced = [];
    const selectedMain = [];
    const order = {};
    allSel.forEach((sp, idx) => { order[sp.id] = idx; });

    function decompose(id) {
      const sp = byId[id];
      const b = vecOf(sp);
      const sol = solveLS(A, b);
      if (!sol || sol.maxRes > 1e-4) return null;
      const nMain = sp.el[main];
      if (!(nMain > 0)) return null;
      const nu = {};
      let nuE = 0;
      basis.forEach((bb, i) => {
        const c = Math.abs(sol.x[i]) < 1e-7 ? 0 : sol.x[i];
        if (c) nu[bb.id] = c;
        if (bb.kind === "electron") nuE = c;
      });
      return { id, sp, nu, nMain, nuE, residual: sol.maxRes, Q: nuE * F_KJ / nMain };
    }

    for (const sp of allSel) {
      const d = decompose(sp.id);
      if (!d) { unbalanced.push(sp.id); continue; }
      decomp[sp.id] = d;
      selectedMain.push(sp.id);
    }
    // 若仍有残差，尝试加入元素标准态
    if (unbalanced.length) {
      let changed = false;
      for (const id of unbalanced.slice()) {
        const sp = byId[id];
        for (const e in sp.el) {
          if (e === main || e === "H" || e === "O") continue;
          const cand = candidates.find(c => c.kind === "element" && c.element === e);
          if (cand && !basisIds.has(cand.comp.id)) {
            if (addBasis(cand.comp, cand.kind, cand)) changed = true;
          }
        }
      }
      if (changed) {
        // rebuild A
        A.length = 0;
        for (let k = 0; k < mRows; k++) A.push(basis.map(b => rows[k] === "z" ? b.z : (b.el[rows[k]] || 0)));
        for (const id of unbalanced.slice()) {
          const d = decompose(id);
          if (d) {
            decomp[id] = d;
            selectedMain.push(id);
            unbalanced.splice(unbalanced.indexOf(id), 1);
          }
        }
      }
    }

    /* ---------------- 化学势 ---------------- */
    function pHat(x) { return xAxis === "pH" ? x : pHfix; }

    function muBasis(x, pH) {
      const p = pH != null ? pH : pHat(x);
      const m = {};
      for (const b of basis) {
        if (b.kind === "main" || b.kind === "element") { m[b.id] = 0; continue; }
        if (b.kind === "proton") { m[b.id] = -RTln10 * p; continue; }
        if (b.kind === "water") { m[b.id] = byId["H2O"].g; continue; }
        if (b.kind === "electron") { m[b.id] = 0; continue; }
        // 若横轴为 lg[Cl-]，且该配体是 Cl-，用横轴值覆盖
        if (xAxis === "pCl" && b.id === "Cl-") {
          // 横轴 x = lg[Cl⁻]，故 a(Cl⁻) = 10^x
          m[b.id] = b.g + RTln10 * x;
          continue;
        }
        if (xAxis === "pNH3" && b.id === "NH3") {
          // 横轴 x = lg[NH₃]（游离氨活度），pH 另由 OH⁻ 控制行固定
          m[b.id] = b.g + RTln10 * x;
          continue;
        }
        if (b.kind === "default") {
          const comp = envFor(b.id);
          const st = comp ? env[comp.id] : null;
          m[b.id] = st && st.on ? ligandMu(b.id, p) : b.g;
          continue;
        }
        m[b.id] = ligandMu(b.id, p);
      }
      return m;
    }

    function baseEnergy(id, x, mu) {
      const d = decomp[id];
      if (!d) return Infinity;
      const M = mu || muBasis(x);
      let e = d.sp.g;
      for (const bid in d.nu) {
        if (bid === "e-" || bid.startsWith("@el:") || bid === "@" + main) continue;
        e -= d.nu[bid] * M[bid];
      }
      return e / d.nMain;
    }

    function activityTerm(id, x, kind) {
      const d = decomp[id];
      const sp = d.sp;
      if (kind === "std") return 0;
      if (sp.state === "aq") {
        const c = Math.max(customC[id] != null ? Math.pow(10, customC[id]) : Ctot / d.nMain, 1e-30);
        return RTln10 * Math.log10(c) / d.nMain;
      }
      if (sp.state === "g") {
        const p = Math.max(customP[id] != null ? customP[id] : Penv, 1e-30);
        return RTln10 * Math.log10(p) / d.nMain;
      }
      return 0;
    }

    // A 模式：每相浓度 = 总浓度（或自定义浓度）
    function gA(id, x, E, mu) {
      const d = decomp[id];
      if (!d) return Infinity;
      return baseEnergy(id, x, mu) + activityTerm(id, x, "A") + d.Q * E;
    }
    // B/C 模式：标准态形成能（浓度由平衡求解）
    function gStd(id, x, E, mu) {
      const d = decomp[id];
      if (!d) return Infinity;
      return baseEnergy(id, x, mu) + d.Q * E;
    }

    function winner(x, E, mu) {
      let best = null, bv = Infinity;
      const M = mu || muBasis(x);
      for (const id of selectedMain) {
        const g = gA(id, x, E, M);
        if (g < bv - 1e-9 || (Math.abs(g - bv) < 1e-9 && best && order[id] < order[best])) {
          bv = g; best = id;
        }
      }
      return best;
    }

    function winnerStd(x, E, mu) {
      let best = null, bv = Infinity;
      const M = mu || muBasis(x);
      for (const id of selectedMain) {
        const g = gStd(id, x, E, M);
        if (g < bv - 1e-9 || (Math.abs(g - bv) < 1e-9 && best && order[id] < order[best])) {
          bv = g; best = id;
        }
      }
      return best;
    }

    /* ---------------- 下包络（A 模式精确区间） ---------------- */
    function envelope(x, emin, emax, mu) {
      const M = mu || muBasis(x);
      const lines = [];
      for (const id of selectedMain) {
        const d = decomp[id];
        lines.push({
          id, a: baseEnergy(id, x, M) + activityTerm(id, x, "A"), q: d.Q,
        });
      }
      if (!lines.length) return [];
      const cuts = [emin, emax];
      for (let i = 0; i < lines.length; i++) {
        for (let j = i + 1; j < lines.length; j++) {
          const dq = lines[i].q - lines[j].q;
          if (Math.abs(dq) < 1e-12) continue;
          const e = (lines[j].a - lines[i].a) / dq;
          if (e > emin + 1e-9 && e < emax - 1e-9) cuts.push(e);
        }
      }
      cuts.sort((a, b) => a - b);
      const uniq = [];
      for (const c of cuts) if (!uniq.length || Math.abs(c - uniq[uniq.length - 1]) > 1e-7) uniq.push(c);
      const out = [];
      for (let k = 0; k < uniq.length - 1; k++) {
        const e0 = uniq[k], e1 = uniq[k + 1];
        if (e1 - e0 < 1e-9) continue;
        const em = (e0 + e1) / 2;
        let best = null, bv = Infinity;
        for (const L of lines) {
          const g = L.a + L.q * em;
          if (g < bv - 1e-9 || (Math.abs(g - bv) < 1e-9 && best && order[L.id] < order[best])) { bv = g; best = L.id; }
        }
        if (best) {
          if (out.length && out[out.length - 1].id === best) out[out.length - 1].e1 = e1;
          else out.push({ id: best, e0, e1 });
        }
      }
      return out;
    }

    /* ---------------- 相区多边形 ---------------- */
    function clipPoly(poly, f, keepLE) {
      if (poly.length < 3) return [];
      const out = [];
      for (let i = 0; i < poly.length; i++) {
        const A0 = poly[i], B0 = poly[(i + 1) % poly.length];
        const fa = f(A0[0], A0[1]), fb = f(B0[0], B0[1]);
        const ina = keepLE ? fa <= 1e-12 : fa >= -1e-12;
        const inb = keepLE ? fb <= 1e-12 : fb >= -1e-12;
        if (ina) out.push(A0);
        if (ina !== inb) {
          const t = fa / (fa - fb);
          out.push([A0[0] + t * (B0[0] - A0[0]), A0[1] + t * (B0[1] - A0[1])]);
        }
      }
      return out;
    }

    function polyArea(poly) {
      let s = 0;
      for (let i = 0; i < poly.length; i++) {
        const a = poly[i], b = poly[(i + 1) % poly.length];
        s += a[0] * b[1] - b[0] * a[1];
      }
      return Math.abs(s) / 2;
    }

    function polygonCentroid(poly) {
      let a = 0, cx = 0, cy = 0;
      for (let i = 0; i < poly.length; i++) {
        const p = poly[i], q = poly[(i + 1) % poly.length];
        const cr = p[0] * q[1] - q[0] * p[1];
        a += cr; cx += (p[0] + q[0]) * cr; cy += (p[1] + q[1]) * cr;
      }
      if (Math.abs(a) < 1e-12) {
        let x = 0, y = 0;
        for (const p of poly) { x += p[0]; y += p[1]; }
        return [x / poly.length, y / poly.length];
      }
      a *= 0.5;
      return [cx / (6 * a), cy / (6 * a)];
    }

    function dominance(xmin, xmax, emin, emax, N) {
      const n = N || 180;
      const xs = [];
      for (let k = 0; k <= n; k++) xs.push(xmin + (xmax - xmin) * k / n);
      const mus = xs.map(x => muBasis(x));
      const out = [];
      for (let k = 0; k < n; k++) {
        const x0 = xs[k], x1 = xs[k + 1];
        const M0 = mus[k], M1 = mus[k + 1];
        const A0 = {}, A1 = {}, Q = {};
        for (const id of selectedMain) {
          A0[id] = baseEnergy(id, x0, M0) + activityTerm(id, x0, "A");
          A1[id] = baseEnergy(id, x1, M1) + activityTerm(id, x1, "A");
          Q[id] = decomp[id].Q;
        }
        for (const id of selectedMain) {
          let poly = [[x0, emin], [x1, emin], [x1, emax], [x0, emax]];
          for (const jd of selectedMain) {
            if (jd === id) continue;
            const dq = Q[id] - Q[jd];
            if (Math.abs(dq) > 1e-10) {
              const e0 = (A0[jd] - A0[id]) / dq;
              const e1 = (A1[jd] - A1[id]) / dq;
              const f = (x, E) => E - (e0 + (e1 - e0) * (x - x0) / (x1 - x0));
              poly = clipPoly(poly, f, dq > 0);
            } else {
              const d0 = A0[id] - A0[jd], d1 = A1[id] - A1[jd];
              if (d0 > TOL && d1 > TOL) { poly = []; break; }
              if (d0 < -TOL && d1 < -TOL) continue;
              if (Math.abs(d0) < TOL && Math.abs(d1) < TOL) {
                if (order[id] > order[jd]) { poly = []; break; }
                continue;
              }
              const xc = x0 + (x1 - x0) * (-d0) / (d1 - d0);
              if (d1 > d0) poly = clipPoly(poly, (x) => x - xc, true);       // 保持 x <= xc
              else poly = clipPoly(poly, (x) => xc - x, true);              // 保持 x >= xc
            }
            if (poly.length < 3) break;
          }
          if (poly.length >= 3) {
            const area = polyArea(poly);
            if (area > 1e-8) out.push({ id, x0, x1, poly, area, centroid: polygonCentroid(poly) });
          }
        }
      }
      return out;
    }

    /* ---------------- 反应配平 ---------------- */
    const reactionCache = new Map();
    function balance(idA, idB, opts) {
      opts = opts || {};
      const key = idA + ">" + idB + "|" + (opts.orient || "");
      if (reactionCache.has(key)) return reactionCache.get(key);
      const dA = decomp[idA], dB = decomp[idB];
      if (!dA || !dB) return null;
      const r = dA.nMain / dB.nMain;
      const raw = [1, r];
      const rawBasis = {};
      const ids = new Set(Object.keys(dA.nu).concat(Object.keys(dB.nu)));
      for (const bid of ids) {
        const c = (dA.nu[bid] || 0) - r * (dB.nu[bid] || 0);
        if (Math.abs(c) > 1e-8) rawBasis[bid] = c;
      }
      const coeffs = raw.concat(Object.values(rawBasis));
      const ints = integerize(coeffs);
      const L = ints.L;
      const coefA = ints.vals[0], coefB = ints.vals[1];
      let k = 2;
      const basisCoef = {};
      for (const bid in rawBasis) { basisCoef[bid] = ints.vals[k++]; }

      let left = [{ id: idA, coef: coefA }];
      let right = [{ id: idB, coef: coefB }];
      for (const bid in basisCoef) {
        const c = basisCoef[bid];
        if (!c) continue;
        if (c > 0) right.push({ id: bid, coef: c });
        else left.push({ id: bid, coef: -c });
      }

      // 负系数归一到另一侧；并消去两侧相同物种（防御性归一化）
      for (const x of left.filter(x => x.coef < -1e-9)) {
        left.splice(left.indexOf(x), 1);
        const t = right.find(y => y.id === x.id);
        if (t) t.coef += -x.coef; else right.push({ id: x.id, coef: -x.coef });
      }
      for (const x of right.filter(x => x.coef < -1e-9)) {
        right.splice(right.indexOf(x), 1);
        const t = left.find(y => y.id === x.id);
        if (t) t.coef += -x.coef; else left.push({ id: x.id, coef: -x.coef });
      }
      const allIds = new Set(left.concat(right).map(x => x.id));
      for (const id of allIds) {
        const l = left.find(x => x.id === id), rr = right.find(x => x.id === id);
        if (l && rr) {
          const m = Math.min(l.coef, rr.coef);
          if (m > 0) { l.coef -= m; rr.coef -= m; }
        }
      }
      left = left.filter(x => x.coef > 1e-9);
      right = right.filter(x => x.coef > 1e-9);

      // 约分
      const allCoefs = left.map(x => x.coef).concat(right.map(x => x.coef));
      const g0 = allCoefs.reduce((a, b) => gcd(a, b), 0);
      if (g0 > 1) { left.forEach(x => x.coef /= g0); right.forEach(x => x.coef /= g0); }

      function eCoef(side) {
        return side.filter(x => x.id === "e-").reduce((a, b) => a + b.coef, 0);
      }
      let eL = eCoef(left), eR = eCoef(right);
      // 统一为“还原方向”：e⁻ 在左
      if (eR > 0 && eL === 0) {
        const t = left; left = right; right = t;
        eL = eCoef(left); eR = eCoef(right);
      }
      // 移除系数为 0 的项
      left = left.filter(x => Math.abs(x.coef) > 1e-9);
      right = right.filter(x => Math.abs(x.coef) > 1e-9);

      function gOf(entry) {
        if (entry.id.startsWith("@")) return 0;
        const sp = byId[entry.id];
        return sp ? sp.g : 0;
      }
      let dG0 = 0;
      for (const x of right) dG0 += x.coef * gOf(x);
      for (const x of left) dG0 -= x.coef * gOf(x);

      const nE = Math.abs(eCoef(left) - eCoef(right));
      const cH = right.filter(x => x.id === "H+").reduce((a, b) => a + b.coef, 0)
               - left.filter(x => x.id === "H+").reduce((a, b) => a + b.coef, 0);
      const ligandIds = new Set(basis.filter(b => b.kind === "ligand").map(b => b.id));
      let ligChange = 0;
      for (const bid of ligandIds) {
        const c = (right.filter(x => x.id === bid).reduce((a, b) => a + b.coef, 0))
                - (left.filter(x => x.id === bid).reduce((a, b) => a + b.coef, 0));
        if (Math.abs(c) > 1e-9) ligChange += Math.abs(c);
      }
      const spA = byId[idA], spB = byId[idB];
      const solid = spA.state === "s" || spA.state === "l" || spB.state === "s" || spB.state === "l";
      const gas = spA.state === "g" || spB.state === "g";
      let type = "phase";
      if (nE > 1e-6) type = "redox";
      else if (ligChange > 1e-9) type = solid ? "precipitation" : "complex";
      else if (Math.abs(cH) > 1e-9) type = "acidbase";

      const logK = -dG0 / RTln10;
      const E0 = nE > 1e-9 ? -dG0 / (nE * F_KJ) : null;
      const slope = nE > 1e-9 ? cH * NERNST / nE : null;
      const pHstd = (nE < 1e-9 && Math.abs(cH) > 1e-9) ? dG0 / (cH * RTln10) : null;

      const out = {
        a: idA, b: idB, left, right, dG0, logK, nE, cH, E0, slope, pHstd,
        type, solid, gas, coefA, coefB, L,
        hasH: left.some(x => x.id === "H+") || right.some(x => x.id === "H+"),
      };
      reactionCache.set(key, out);
      return out;
    }

    function classify(idA, idB) {
      const r = balance(idA, idB);
      return r ? r.type : "phase";
    }

    // 碱性写法：H⁺ + OH⁻ → H₂O
    function alkalineForm(reaction) {
      const left = reaction.left.map(x => ({ ...x }));
      const right = reaction.right.map(x => ({ ...x }));
      function apply(side, c, acidOnRight) {
        // c H+ 在 side 上；H+ = H2O - OH-
        const h2o = side === "left" ? left : right;
        const oh = side === "left" ? right : left;
        function add(arr, id, d) {
          const f = arr.find(x => x.id === id);
          if (f) f.coef += d; else arr.push({ id, coef: d });
        }
        // 移除 H+
        const h = (side === "left" ? left : right).filter(x => x.id === "H+");
        for (const x of h) x.coef = 0;
        // 加 H2O 到同一侧，加 OH- 到另一侧
        add(h2o, "H2O", c);
        add(oh, "OH-", c);
      }
      const hL = left.filter(x => x.id === "H+").reduce((a, b) => a + b.coef, 0);
      const hR = right.filter(x => x.id === "H+").reduce((a, b) => a + b.coef, 0);
      if (hL) apply("left", hL);
      if (hR) apply("right", hR);
      // 消去两侧相同物种
      for (const side of [left, right]) side.forEach(x => { if (Math.abs(x.coef) < 1e-9) x.coef = 0; });
      const ids = new Set(left.concat(right).map(x => x.id));
      for (const id of ids) {
        if (id === "H+" || id === "e-") continue;
        const l = left.find(x => x.id === id), r = right.find(x => x.id === id);
        if (l && r) {
          const m = Math.min(l.coef, r.coef);
          if (m > 0) { l.coef -= m; r.coef -= m; }
        }
      }
      return {
        ...reaction, left: left.filter(x => x.coef > 1e-9), right: right.filter(x => x.coef > 1e-9),
      };
    }

    /* ---------------- 平衡求解（B/C 模式） ---------------- */
    function aqList() { return selectedMain.filter(id => decomp[id].sp.state === "aq"); }
    function solidList() { return selectedMain.filter(id => decomp[id].sp.state === "s" || decomp[id].sp.state === "l"); }
    function gasList() { return selectedMain.filter(id => decomp[id].sp.state === "g"); }

    function solveMu(x, Ct, M, fixedMu) {
      const aq = aqList();
      if (!aq.length) return null;
      const G = {}, n = {};
      const G0 = fixedMu && fixedMu.G0 ? fixedMu.G0 : null;
      let gmin = Infinity, best = null;
      for (const id of aq) {
        G[id] = (G0 ? G0[id] : baseEnergy(id, x, M)) + decomp[id].Q * (fixedMu != null ? fixedMu.E : 0);
        n[id] = decomp[id].nMain;
        if (G[id] < gmin) { gmin = G[id]; best = id; }
      }
      // 若给了 E，G 是标准形成能 + QE；这里 fixedMu 传 {E}
      let mu;
      if (fixedMu && fixedMu.mu != null) mu = fixedMu.mu;
      else mu = gmin + RTln10 * Math.log10(Ct) / Math.max(n[best], 1);
      const target = Math.log10(Ct);
      const maxIter = fixedMu && fixedMu.maxIter ? fixedMu.maxIter : 20;
      for (let it = 0; it < maxIter; it++) {
        let mx = -Infinity;
        const logs = {};
        for (const id of aq) {
          const l = n[id] * (mu - G[id]) / RTln10;
          logs[id] = l;
          if (l > mx) mx = l;
        }
        let sum = 0, sumN = 0, sumN2 = 0;
        for (const id of aq) {
          const c = Math.pow(10, Math.max(-300, logs[id] - mx));
          sum += c;
          sumN += n[id] * c;
          sumN2 += n[id] * n[id] * c;
        }
        if (sum <= 0) { mu += 10; continue; }
        const logSum = mx + Math.log10(sum);
        const f = logSum - target;
        if (Math.abs(f) < 1e-10) break;
        const deriv = sumN2 / (RTln10 * sumN); // d logSum / d mu
        let step = f / Math.max(deriv, 1e-12);
        if (!isFinite(step)) step = 1;
        step = Math.max(-20, Math.min(20, step));
        mu -= step;
      }
      const conc = {};
      let mx = -Infinity;
      for (const id of aq) {
        const l = n[id] * (mu - G[id]) / RTln10;
        if (l > mx) mx = l;
      }
      let sumN = 0;
      for (const id of aq) {
        const c = Math.pow(10, Math.max(-300, n[id] * (mu - G[id]) / RTln10 - mx));
        conc[id] = c;
        sumN += n[id] * c;
      }
      // 上面的 conc 是缩放后的；重新归一化到 Ct
      if (sumN > 0) for (const id of aq) conc[id] *= Ct / sumN;
      return { mu, conc };
    }

    function speciate(x, E, opts) {
      opts = opts || {};
      const M = opts.mu || muBasis(x);
      const G0in = opts.G0 || null;
      const maxIter = opts.maxIter || 20;
      const Ct = Math.max(opts.Ctot != null ? opts.Ctot : Ctot, 1e-30);
      const aq = aqList(), solids = solidList(), gases = gasList();
      const G = {};
      for (const id of selectedMain) G[id] = (G0in ? G0in[id] : baseEnergy(id, x, M)) + decomp[id].Q * E;

      // 没有溶质：取最低固/气
      if (!aq.length) {
        let best = null, bv = Infinity;
        for (const id of solids.concat(gases)) if (G[id] < bv) { bv = G[id]; best = id; }
        const present = best ? [{
          id: best, frac: 1, solid: byId[best].state !== "g", gas: byId[best].state === "g",
          p: byId[best].state === "g" ? Penv : null,
        }] : [];
        return { conc: {}, present, solid: null, gas: null, Csol: 0, dissolved: 0, mu: null };
      }

      // 初解（无固相）
      let sol = solveMu(x, Ct, M, { E, G0: G0in, maxIter });
      if (!sol) return { conc: {}, present: [], solid: null, gas: null, Csol: 0, dissolved: 0, mu: null };
      let mu = sol.mu, conc = sol.conc;
      let dissolved = 0;
      for (const id of aq) dissolved += decomp[id].nMain * conc[id];

      // 固相缓冲：取最低 g 的固相
      let solid = null;
      if (solids.length) {
        let best = null, bv = Infinity;
        for (const id of solids) if (G[id] < bv) { bv = G[id]; best = id; }
        if (best && mu > bv + 1e-9) {
          mu = bv;
          // 按 μ = G_solid 直接计算饱和浓度（绝对浓度，不做归一化）
          const aqIds = aq;
          conc = {};
          let sN = 0;
          for (const id of aqIds) {
            const l = Math.max(-300, Math.min(300, decomp[id].nMain * (mu - G[id]) / RTln10));
            const c = Math.pow(10, l);
            conc[id] = c; sN += decomp[id].nMain * c;
          }
          // 数值兜底：若饱和浓度仍超过总量，按质量平衡缩放
          if (sN > Ct && sN > 0) {
            for (const id of aqIds) conc[id] *= Ct / sN;
            sN = Ct;
          }
          dissolved = sN;
          const frac = Math.max(0, Math.min(1, (Ct - dissolved) / Ct));
          solid = { id: best, frac, amount: (Ct - dissolved) };
        }
      }

      // 气体：按 Henry 限制溶解量。场渲染用 simpleGas 只钳制一次，避免逐点重解。
      let gas = null, gasLoss = 0;
      function scanGas() {
        let loss = 0, gasInfo = null;
        for (const id of aq) {
          const gasId = GAS_PAIR[id] || null;
          if (!gasId) continue;
          const gsp = byId[gasId];
          if (!gsp) continue;
          const H = Math.exp((byId[id].g - gsp.g) * 1000 / (R * T)); // bar/M
          const c = conc[id] || 0;
          const p = H * c;
          if (p > Penv * 1.0000001) {
            const cmax = Penv / H;
            const excess = c - cmax;
            loss += decomp[id].nMain * excess;
            if (!gasInfo || excess > gasInfo.excess) {
              gasInfo = { id: gasId, aqId: id, p: Penv, H, excess, cmax };
            }
            conc[id] = cmax;
          } else if (!gasInfo && p > 0 && opts.keepGas) {
            gasInfo = { id: gasId, aqId: id, p, H, excess: 0, cmax: c };
          }
        }
        return { loss, gasInfo };
      }
      if (opts.simpleGas) {
        const r = scanGas();
        gasLoss = r.loss;
        if (r.gasInfo) gas = { ...r.gasInfo, frac: Math.min(1, r.loss / Ct), amount: r.loss };
      } else {
        for (let iter = 0; iter < 2; iter++) {
          const r = scanGas();
          const loss = r.loss, gasInfo = r.gasInfo;
          if (loss > 0) gasLoss = Math.max(gasLoss, loss);
          if (loss < 1e-12) { if (gasInfo) gas = gasInfo; break; }
          const newCt = Math.max(Ct - loss, 1e-30);
          const sol2 = solveMu(x, newCt, M, { E, G0: G0in, maxIter });
          if (!sol2) break;
          mu = sol2.mu; conc = sol2.conc;
          dissolved = 0;
          for (const id of aq) dissolved += decomp[id].nMain * conc[id];
          if (gasInfo) gas = gasInfo;
        }
      }

      // 重新检查固相（气体损失后可能不再饱和）
      if (solid) {
        const frac = Math.max(0, Math.min(1, (Ct - gasLoss - dissolved) / Ct));
        solid = { ...solid, frac, amount: Math.max(0, Ct - gasLoss - dissolved) };
      } else if (gasLoss > 0) {
        // 无固相但气体逸出：气体本身占主元素的一部分
        if (gas) gas = { ...gas, frac: Math.min(1, gasLoss / Ct), amount: gasLoss };
      }

      const present = [];
      for (const id of aq) {
        const c = conc[id];
        if (c > 0 && decomp[id].sp.state === "aq") {
          present.push({ id, conc: c, frac: decomp[id].nMain * c / Ct, solid: false, gas: false });
        }
      }
      // 气相（若有）
      if (gas) {
        present.push({
          id: gas.id, frac: gas.frac != null ? gas.frac : Math.min(1, (gas.excess || 0) / Ct),
          solid: false, gas: true,
          p: gas.p, aqId: gas.aqId, H: gas.H, conc: gas.cmax,
        });
      }
      // 没有水相前体的气体（N₂、NO、NO₂、CO 等）：
      // 用主元素化学势反算平衡分压，C 模式按 p/P_env 着色，而不是按溶解浓度。
      if (mu != null) {
        const pairedGases = new Set(Object.values(GAS_PAIR));
        const already = new Set(present.filter(p => p.gas).map(p => p.id));
        for (const gid of gases) {
          if (pairedGases.has(gid) || already.has(gid)) continue;
          const d = decomp[gid];
          if (!d || G[gid] == null) continue;
          const lp = d.nMain * (mu - G[gid]) / RTln10;
          if (!isFinite(lp)) continue;
          const pRaw = Math.pow(10, Math.max(-300, Math.min(300, lp)));
          if (pRaw <= 1e-12) continue;
          present.push({
            id: gid, frac: Math.min(1, pRaw / Penv), solid: false, gas: true,
            p: Math.min(Penv, pRaw), rawP: pRaw,
          });
        }
      }
      if (solid) present.push({ id: solid.id, frac: solid.frac, solid: true, gas: false, amount: solid.amount });
      return { conc, present, solid, gas, gasLoss, Csol: dissolved, dissolved, mu };
    }

    /* ---------------- 气体体积（死体积 + Ar 稀释） ---------------- */
    function gasState(x, E, opts) {
      opts = opts || {};
      const sp = speciate(x, E, opts);
      const Vsol = opts.Vsol != null ? opts.Vsol : 1; // L
      const V0 = opts.V0 != null ? opts.V0 : 0;
      const volatiles = [];
      for (const id of aqList()) {
        const gasId = GAS_PAIR[id] || null;
        if (!gasId) continue;
        const gsp = byId[gasId];
        if (!gsp) continue;
        const H = Math.exp((byId[id].g - gsp.g) * 1000 / (R * T));
        const c = (sp.conc[id] || 0);
        if (c > 0) volatiles.push({ aqId: id, gasId, H, C: decomp[id].nMain * c });
      }
      if (!volatiles.length) return { p: {}, Vg: V0, Ar: Penv * V0 / (RGAS_LBAR * T), volatiles };
      const nAr = Penv * V0 / (RGAS_LBAR * T);
      const p = {}, Vg = { v: V0 };
      // 单气体解析 + 多气体迭代
      for (const v of volatiles) p[v.gasId] = Math.min(v.H * v.C, Penv * 0.9);
      for (let it = 0; it < 30; it++) {
        let sumP = 0;
        for (const k in p) sumP += p[k];
        if (sumP >= Penv * 0.999999) {
          const f = Penv * 0.999999 / sumP;
          for (const k in p) p[k] *= f;
          sumP = Penv * 0.999999;
        }
        const vg = nAr > 0 ? nAr * RGAS_LBAR * T / Math.max(Penv - sumP, 1e-12) : 0;
        let changed = false;
        for (const v of volatiles) {
          const denom = Vsol / v.H + vg / (RGAS_LBAR * T);
          const pn = denom > 0 ? v.C * Vsol / denom : 0;
          if (Math.abs(pn - p[v.gasId]) > 1e-9) changed = true;
          p[v.gasId] = Math.max(0, Math.min(pn, Penv * 0.999999));
        }
        Vg.v = nAr > 0 ? nAr * RGAS_LBAR * T / Math.max(Penv - sumP, 1e-12) : 0;
        if (!changed && it > 2) break;
      }
      return { p, Vg: Vg.v, Ar: nAr, volatiles, totals: volatiles.map(v => v.C) };
    }

    /* ---------------- 反应列表 ---------------- */
    function reactionList(opts) {
      opts = opts || {};
      const ids = selectedMain.slice();
      const stablePairs = [];
      if (opts.boundaries) {
        const seen = new Set();
        for (const s of opts.boundaries) {
          const key = [s.a, s.b].sort().join("|");
          if (seen.has(key)) continue;
          seen.add(key);
          stablePairs.push([s.a, s.b, (s.x1 + s.x2) / 2, (s.e1 + s.e2) / 2]);
        }
      }
      const groups = { stable: [], redox: [], acidbase: [], precipitation: [], complex: [], phase: [] };
      const seen = new Set();
      function push(group, r, at, atE) {
        if (!r) return;
        const key = [r.a, r.b].sort().join("|") + "|" + r.type;
        if (seen.has(key)) return;
        seen.add(key);
        if (at != null) r.at = at;
        if (atE != null) r.atE = atE;
        groups[group].push(r);
      }
      for (const [a, b, at, atE] of stablePairs) {
        const r = balance(a, b);
        if (r) push("stable", r, at, atE);
      }
      for (let i = 0; i < ids.length; i++) {
        for (let j = i + 1; j < ids.length; j++) {
          const a = ids[i], b = ids[j];
          const r = balance(a, b);
          if (!r) continue;
          if (r.type === "redox") {
            if (r.nE <= 8) push("redox", r);
          } else if (r.type === "acidbase") {
            push("acidbase", r);
          } else if (r.type === "precipitation") {
            push("precipitation", r);
          } else if (r.type === "complex") {
            push("complex", r);
          } else if (r.type === "phase") {
            push("phase", r);
          }
        }
      }
      // 去重：稳定相界优先
      const stableKeys = new Set(groups.stable.map(r => [r.a, r.b].sort().join("|") + "|" + r.type));
      for (const k of ["redox", "acidbase", "precipitation", "complex", "phase"]) {
        groups[k] = groups[k].filter(r => !stableKeys.has([r.a, r.b].sort().join("|") + "|" + r.type));
      }
      return groups;
    }

    /* ---------------- 介质线 ---------------- */
    function pKaList() {
      const out = [];
      const seen = new Set();
      for (const comp of D.ENV) {
        const st = env[comp.id];
        if (!st || !st.on) continue;
        if (!comp.family || comp.family === "water") continue;
        const info = familyInfo(comp.family);
        if (!info) continue;
        for (let j = 1; j < info.list.length; j++) {
          const key = info.list[j - 1] + "|" + info.list[j];
          if (seen.has(key)) continue;
          seen.add(key);
          out.push({ pH: info.pKa[j - 1], from: info.list[j - 1], to: info.list[j], family: comp.family });
        }
      }
      return out;
    }

    function redoxLine(couple, x) {
      const ox = byId[couple.ox], red = byId[couple.red];
      if (!ox || !red) return null;
      // 用一个小体系配平该电对
      const sub = createSystem({
        main: ox.main || red.main, selected: [ox.id, red.id],
        env: {}, Ctot: 1, Penv: Penv, xAxis: "pH", pHfix: pHat(x),
      });
      const r = sub.balance(ox.id, red.id);
      if (!r || r.nE < 1e-9) return null;
      const st = env[couple.id];
      const c = st && st.c != null ? st.c : 0.01;
      const ratio = st && st.ratio != null ? st.ratio : 1;
      const cOx = c * ratio / (1 + ratio), cRed = c / (1 + ratio);
      let q = 0;
      for (const e of r.right) {
        if (e.id === "H+" || e.id === "H2O" || e.id === "e-") continue;
        const sp = byId[e.id];
        if (!sp) continue;
        if (sp.state === "aq") q += e.coef * Math.log10((sp.id === ox.id ? cOx : sp.id === red.id ? cRed : 0.01));
      }
      for (const e of r.left) {
        if (e.id === "H+" || e.id === "H2O" || e.id === "e-") continue;
        const sp = byId[e.id];
        if (!sp) continue;
        if (sp.state === "aq") q -= e.coef * Math.log10((sp.id === ox.id ? cOx : sp.id === red.id ? cRed : 0.01));
      }
      const E = r.E0 - (r.nE ? NERNST / r.nE * q : 0) + (r.slope || 0) * pHat(x);
      return { E, r, label: couple.id };
    }

    return {
      cfg, main, env, Ctot, Penv, xAxis, pHfix, couple,
      basis, elementList, decomp, selectedMain, allSel, unbalanced,
      muBasis, baseEnergy, gA, gStd, winner, winnerStd,
      envelope, dominance, clipPoly, polyArea, polygonCentroid,
      balance, alkalineForm, classify, reactionList,
      speciate, gasState, aqList, solidList, gasList,
      familyInfo, familyFractions, ligandAlpha, pKaList, redoxLine,
      activityTerm,
    };
  }

  return {
    R, T, F, F_KJ, RTln10, NERNST, RGAS_LBAR,
    solveLin, solveLS, rankCols, gcd, integerize,
    familyInfo, familyFractions, ligandAlpha,
    createSystem,
  };
})(typeof PBData !== "undefined" ? PBData : require("./data.js"));
if (typeof module !== "undefined") module.exports = PBChem;