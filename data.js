/* ============================================================
   data.js — 纯数据层
   25 °C，理想溶液（活度 ≈ 浓度；纯固体/纯水 a = 1）
   ΔGf° 单位 kJ·mol⁻¹。NBS = Wagman et al. (1982)。
   配合物/沉淀的 ΔGf° 由 logβ / Ksp 现场推得，保证
   “数据 → 反应式 → 相图” 三者自洽；note 中记录所用常数。
   ============================================================ */
var PBData = (function () {
  "use strict";
  const RTln10 = 5.708;          // kJ·mol⁻¹ at 298.15 K
  const DG_H2O = -237.14;

  const DB = [];
  const byId = Object.create(null);

  function add(o) {
    if (byId[o.id]) throw new Error("duplicate species id: " + o.id);
    o.cats = o.cats || [];
    o.ligs = o.ligs || [];
    o.color = o.color || null;
    o.derived = !!o.derived;
    o.src = o.src || (o.derived ? "推得" : "NBS 1982");
    byId[o.id] = o;
    DB.push(o);
    return o;
  }

  // 简单物种。ox 为“主元素”的氧化数；多元素物种用 oxByEl 覆盖。
  function S(id, formula, el, z, state, g, ox, main, name, extra) {
    return add(Object.assign({
      id, formula, el, z, state, g, ox, main, name,
      oxByEl: Object.assign({}, (extra && extra.oxByEl) || {}),
    }, extra || {}));
  }

  // 配合物：自由离子 + n 配体 → 配合物，logβ
  function Cx(id, formula, el, z, ox, main, name, baseId, ligands, logB, cats, extra) {
    const base = byId[baseId];
    if (!base) throw new Error("missing complex base " + baseId);
    let g = base.g;
    for (const [lid, n] of ligands) {
      const L = byId[lid];
      if (!L) throw new Error("missing ligand " + lid);
      g += n * L.g;
    }
    g -= RTln10 * logB;
    return add(Object.assign({
      id, formula, el, z, state: "aq", g, ox, main, name,
      cats: cats || ["络合物"], ligs: ligands.map(x => x[0]),
      derived: true,
      note: `由 logβ = ${logB} 推得（${baseId} + ${ligands.map(x => x[1] + " " + x[0]).join(" + ")}）`,
    }, extra || {}));
  }

  // 沉淀：M_aX_b(s) → a M + b X，logKsp
  function Ppt(id, formula, el, z, ox, main, name, ions, logKsp, cats, extra) {
    let g = 0;
    for (const [iid, n] of ions) {
      const ion = byId[iid];
      if (!ion) throw new Error("missing precipitate ion " + iid);
      g += n * ion.g;
    }
    g += RTln10 * logKsp;
    return add(Object.assign({
      id, formula, el, z, state: "s", g, ox, main, name,
      cats: cats || ["难溶物"], derived: true,
      note: `由 logKsp = ${logKsp} 推得（溶解：${formula} ⇌ ${ions.map(x => x[1] + " " + x[0]).join(" + ")}）`,
    }, extra || {}));
  }

  /* ================= 介质与基准 ================= */
  S("H2O", "H2O", { H: 2, O: 1 }, 0, "l", DG_H2O, null, null, "水（溶剂）", {
    cats: ["酸碱"], oxByEl: { H: 1, O: -2 }, color: "#eef6ff", solvent: true,
  });
  S("H+", "H^+", { H: 1 }, 1, "aq", 0, null, null, "溶剂化质子", {
    cats: ["酸碱"], oxByEl: { H: 1 },
  });
  S("OH-", "OH^-", { O: 1, H: 1 }, -1, "aq", -157.24, null, null, "氢氧根", {
    cats: ["酸碱"], oxByEl: { O: -2, H: 1 }, color: "#cfe0ff",
  });
  S("e-", "e^-", {}, -1, "e", 0, null, null, "电子", { cats: [] });
  S("H2", "H2", { H: 2 }, 0, "g", 0, 0, "H", "氢气", {
    cats: ["单质", "易挥发物"], color: "#cfe8ff",
  });
  S("O2", "O2", { O: 2 }, 0, "g", 0, 0, "O", "氧气", {
    cats: ["单质", "易挥发物"], color: "#bcd4ff",
  });
  S("H2O2", "H2O2", { H: 2, O: 2 }, 0, "aq", -134.0, -1, "O", "过氧化氢", {
    cats: ["酸碱"], color: "#e8f4ff",
  });

  /* ================= Mg ================= */
  S("Mg", "Mg", { Mg: 1 }, 0, "s", 0, 0, "Mg", "镁", { cats: ["单质"], color: "#c8ccd0" });
  S("Mg2+", "Mg^2+", { Mg: 1 }, 2, "aq", -454.8, 2, "Mg", "镁离子", {
    cats: ["离子", "酸碱"], color: "#d8dce0",
  });
  S("Mg(OH)2", "Mg(OH)2", { Mg: 1, O: 2, H: 2 }, 0, "s", -833.5, 2, "Mg", "氢氧化镁↓", {
    cats: ["难溶物", "酸碱"], color: "#f4f4f0",
  });
  S("MgO", "MgO", { Mg: 1, O: 1 }, 0, "s", -569.3, 2, "Mg", "氧化镁", {
    cats: ["氧化物", "酸碱"], color: "#f2f2f2",
  });

  /* ================= Al ================= */
  S("Al", "Al", { Al: 1 }, 0, "s", 0, 0, "Al", "铝", { cats: ["单质"], color: "#c4c8cc" });
  S("Al3+", "Al^3+", { Al: 1 }, 3, "aq", -485.0, 3, "Al", "铝离子", {
    cats: ["离子", "酸碱"], color: "#d8dce0",
  });
  S("Al(OH)3", "Al(OH)3", { Al: 1, O: 3, H: 3 }, 0, "s", -1154.9, 3, "Al", "氢氧化铝↓（三水铝石）", {
    cats: ["难溶物", "酸碱"], color: "#f6f6f2",
    note: "NBS 三水铝石（gibbsite）值；无定形相约 -1137.6",
  });
  S("AlO(OH)", "AlO(OH)", { Al: 1, O: 2, H: 1 }, 0, "s", -918.4, 3, "Al", "羟基氧化铝（勃姆石）", {
    cats: ["难溶物", "酸碱"], color: "#ececec",
  });
  S("Al2O3", "Al2O3", { Al: 2, O: 3 }, 0, "s", -1582.3, 3, "Al", "氧化铝", {
    cats: ["氧化物", "酸碱"], color: "#f0f0f0",
  });
  S("Al(OH)4-", "Al(OH)4^-", { Al: 1, O: 4, H: 4 }, -1, "aq", -1305.8, 3, "Al", "四羟基合铝酸根", {
    cats: ["离子", "络合物", "酸碱"], color: "#dce8e0",
  });
  /* ================= P ================= */
  S("P4", "P4", { P: 4 }, 0, "s", 0, 0, "P", "白磷", { cats: ["单质"], color: "#f5f0c0" });
  S("PH3", "PH3", { P: 1, H: 3 }, 0, "g", 13.4, -3, "P", "磷化氢", {
    cats: ["易挥发物", "酸碱"], color: "#e8e8e8",
  });
  S("H3PO4", "H3PO4", { H: 3, P: 1, O: 4 }, 0, "aq", -1142.6, 5, "P", "磷酸", {
    cats: ["酸碱"], color: "#f0e8d8",
  });
  S("H2PO4-", "H2PO4^-", { H: 2, P: 1, O: 4 }, -1, "aq", -1130.2, 5, "P", "磷酸二氢根", {
    cats: ["离子", "酸碱"], color: "#eee6d6",
  });
  S("HPO42-", "HPO4^2-", { H: 1, P: 1, O: 4 }, -2, "aq", -1089.1, 5, "P", "磷酸氢根", {
    cats: ["离子", "酸碱"], color: "#ece4d4",
  });
  S("PO43-", "PO4^3-", { P: 1, O: 4 }, -3, "aq", -1018.7, 5, "P", "磷酸根", {
    cats: ["离子", "酸碱"], color: "#eae2d2",
  });
  S("H3PO3", "H3PO3", { H: 3, P: 1, O: 3 }, 0, "aq", -852.2, 3, "P", "亚磷酸", {
    cats: ["酸碱"], color: "#f0e8d8",
  });
  S("H2PO3-", "H2PO3^-", { H: 2, P: 1, O: 3 }, -1, "aq", -844.8, 3, "P", "亚磷酸二氢根", {
    cats: ["离子", "酸碱"], color: "#eee6d6",
  });
  S("HPO32-", "HPO3^2-", { H: 1, P: 1, O: 3 }, -2, "aq", -806.6, 3, "P", "亚磷酸氢根", {
    cats: ["离子", "酸碱"], color: "#ece4d4",
  });
  S("H3PO2", "H3PO2", { H: 3, P: 1, O: 2 }, 0, "aq", -518.6, 1, "P", "次磷酸", {
    cats: ["酸碱"], color: "#f0e8d8",
  });
  S("H2PO2-", "H2PO2^-", { H: 2, P: 1, O: 2 }, -1, "aq", -511.7, 1, "P", "次磷酸根", {
    cats: ["离子", "酸碱"], color: "#eee6d6",
  });

  /* ================= S ================= */
  S("S8", "S8", { S: 8 }, 0, "s", 0, 0, "S", "环八硫（斜方硫）", {
    cats: ["单质"], color: "#e8d44d",
  });
  S("S", "S", { S: 1 }, 0, "s", 0.2, 0, "S", "硫（亚稳态）", {
    cats: ["单质"], color: "#e8d44d", derived: true,
    note: "与 S₈ 近简并，取 +0.2 kJ·mol⁻¹ 避免退化；亚稳态",
  });
  S("H2S", "H2S", { H: 2, S: 1 }, 0, "aq", -27.8, -2, "S", "硫化氢（溶解）", {
    cats: ["酸碱", "易挥发物"], color: "#e8e8c0",
  });
  S("H2Sg", "H2S", { H: 2, S: 1 }, 0, "g", -33.6, -2, "S", "硫化氢（气）", {
    cats: ["易挥发物", "酸碱"], color: "#e0e0b8", aqPair: "H2S",
  });
  S("HS-", "HS^-", { H: 1, S: 1 }, -1, "aq", 12.0, -2, "S", "硫氢根", {
    cats: ["离子", "酸碱"], color: "#e8e0a8",
  });
  S("S2-", "S^2-", { S: 1 }, -2, "aq", 85.8, -2, "S", "硫离子", {
    cats: ["离子", "酸碱"], color: "#e0d890",
  });
  S("SO2", "SO2", { S: 1, O: 2 }, 0, "aq", -300.2, 4, "S", "二氧化硫（溶解）", {
    cats: ["氧化物", "酸碱", "易挥发物"], color: "#f0e0a0",
  });
  S("SO2g", "SO2", { S: 1, O: 2 }, 0, "g", -300.1, 4, "S", "二氧化硫（气）", {
    cats: ["氧化物", "易挥发物"], color: "#e8d890", aqPair: "SO2",
  });
  S("HSO3-", "HSO3^-", { H: 1, S: 1, O: 3 }, -1, "aq", -527.8, 4, "S", "亚硫酸氢根", {
    cats: ["离子", "酸碱"], color: "#eee0a8",
  });
  S("SO32-", "SO3^2-", { S: 1, O: 3 }, -2, "aq", -486.6, 4, "S", "亚硫酸根", {
    cats: ["离子", "酸碱"], color: "#ece0a0",
  });
  S("HSO4-", "HSO4^-", { H: 1, S: 1, O: 4 }, -1, "aq", -755.9, 6, "S", "硫酸氢根", {
    cats: ["离子", "酸碱"], color: "#eee8d8",
  });
  S("SO42-", "SO4^2-", { S: 1, O: 4 }, -2, "aq", -744.5, 6, "S", "硫酸根", {
    cats: ["离子", "酸碱"], color: "#ece6d6",
  });
  S("S2O32-", "S2O3^2-", { S: 2, O: 3 }, -2, "aq", -522.5, 2, "S", "硫代硫酸根", {
    cats: ["离子", "酸碱"], color: "#e8e0b0",
  });
  S("S4O62-", "S4O6^2-", { S: 4, O: 6 }, -2, "aq", -1060.4, 2.5, "S", "连四硫酸根", {
    cats: ["离子"], color: "#e8e0b0", derived: true, note: "NBS 近似值",
  });
  S("S2O42-", "S2O4^2-", { S: 2, O: 4 }, -2, "aq", -473.0, 3, "S", "连二亚硫酸根（保险粉）", {
    cats: ["离子"], color: "#e8e0b0", derived: true, note: "NBS 近似值",
  });
  S("S2O52-", "S2O5^2-", { S: 2, O: 5 }, -2, "aq", -818.5, 4, "S", "焦亚硫酸根", {
    cats: ["离子"], color: "#e8e0b0", derived: true, note: "近似值",
  });
  S("S2O72-", "S2O7^2-", { S: 2, O: 7 }, -2, "aq", -1251.9, 6, "S", "焦硫酸根", {
    cats: ["离子"], color: "#e8e0b0", derived: true, note: "近似值",
  });
  S("S2O82-", "S2O8^2-", { S: 2, O: 8 }, -2, "aq", -1115.0, 7, "S", "过二硫酸根", {
    cats: ["离子"], color: "#e8e0b0",
  });

  /* ================= Cl ================= */
  S("Cl-", "Cl^-", { Cl: 1 }, -1, "aq", -131.2, -1, "Cl", "氯离子", {
    cats: ["离子"], color: "#d8e8a0",
  });
  S("Cl2", "Cl2", { Cl: 2 }, 0, "aq", 6.94, 0, "Cl", "氯气（溶解）", {
    cats: ["单质", "易挥发物"], color: "#d8e8a0",
  });
  S("Cl2g", "Cl2", { Cl: 2 }, 0, "g", 0, 0, "Cl", "氯气（气）", {
    cats: ["单质", "易挥发物"], color: "#c8d860", aqPair: "Cl2",
  });
  S("HClO", "HClO", { H: 1, Cl: 1, O: 1 }, 0, "aq", -79.9, 1, "Cl", "次氯酸", {
    cats: ["酸碱"], color: "#e8e8b0",
  });
  S("ClO-", "ClO^-", { Cl: 1, O: 1 }, -1, "aq", -36.8, 1, "Cl", "次氯酸根", {
    cats: ["离子", "酸碱"], color: "#e8e0a0",
  });
  S("ClO2g", "ClO2", { Cl: 1, O: 2 }, 0, "g", 120.0, 4, "Cl", "二氧化氯（气）", {
    cats: ["氧化物", "易挥发物"], color: "#e8e060",
  });
  S("ClO3-", "ClO3^-", { Cl: 1, O: 3 }, -1, "aq", -3.4, 5, "Cl", "氯酸根", {
    cats: ["离子"], color: "#e8e8c0",
  });
  S("ClO4-", "ClO4^-", { Cl: 1, O: 4 }, -1, "aq", -8.6, 7, "Cl", "高氯酸根", {
    cats: ["离子"], color: "#eee8d0",
  });
  S("HCl", "HCl", { H: 1, Cl: 1 }, 0, "aq", -120.0, -1, "Cl", "盐酸（形式物种）", {
    cats: ["酸碱"], color: "#e8e8d0", derived: true,
    note: "强酸完全解离；取 -120 而非 -131.2，避免与 Cl⁻ + H⁺ 退化",
  });
  S("HClg", "HCl", { H: 1, Cl: 1 }, 0, "g", -95.3, -1, "Cl", "氯化氢（气）", {
    cats: ["易挥发物", "酸碱"], color: "#e8e8d0", aqPair: "HCl",
  });

  /* ================= I ================= */
  S("I-", "I^-", { I: 1 }, -1, "aq", -51.6, -1, "I", "碘离子", {
    cats: ["离子"], color: "#c8b0d8",
  });
  S("I3-", "I3^-", { I: 3 }, -1, "aq", -51.5, -1 / 3, "I", "三碘离子", {
    cats: ["离子"], color: "#8a5a2a", derived: true,
    note: "NBS 近似；I₂(s) + I⁻ ⇌ I₃⁻ 的 logK 约 0（数据源差异）",
  });
  S("I2", "I2", { I: 2 }, 0, "s", 0, 0, "I", "碘", { cats: ["单质"], color: "#5a3a6a" });
  S("ICl", "ICl", { I: 1, Cl: 1 }, 0, "s", -13.0, 1, "I", "一氯化碘", {
    cats: ["单质", "难溶物"], color: "#7a4a3a", derived: true, note: "近似值",
  });
  S("IO3-", "IO3^-", { I: 1, O: 3 }, -1, "aq", -128.0, 5, "I", "碘酸根", {
    cats: ["离子"], color: "#e8e0c0",
  });

  /* ================= N ================= */
  S("N2", "N2", { N: 2 }, 0, "g", 0, 0, "N", "氮气", {
    cats: ["单质", "易挥发物"], color: "#d0d8e8",
  });
  S("NH4+", "NH4^+", { N: 1, H: 4 }, 1, "aq", -79.3, -3, "N", "铵根", {
    cats: ["离子", "酸碱"], color: "#d8e8e0",
  });
  S("NH3", "NH3", { N: 1, H: 3 }, 0, "aq", -26.5, -3, "N", "氨（溶解）", {
    cats: ["酸碱", "易挥发物"], color: "#d0e8d8",
  });
  S("NH3g", "NH3", { N: 1, H: 3 }, 0, "g", -16.4, -3, "N", "氨（气）", {
    cats: ["易挥发物", "酸碱"], color: "#c8e0d0", aqPair: "NH3",
  });
  S("HNO2", "HNO2", { H: 1, N: 1, O: 2 }, 0, "aq", -46.0, 3, "N", "亚硝酸", {
    cats: ["酸碱"], color: "#e8d8d0",
  });
  S("NO2-", "NO2^-", { N: 1, O: 2 }, -1, "aq", -32.2, 3, "N", "亚硝酸根", {
    cats: ["离子", "酸碱"], color: "#e8d0c0",
  });
  S("NO3-", "NO3^-", { N: 1, O: 3 }, -1, "aq", -111.3, 5, "N", "硝酸根", {
    cats: ["离子"], color: "#e0e0e8",
  });
  S("NO", "NO", { N: 1, O: 1 }, 0, "g", 86.6, 2, "N", "一氧化氮", {
    cats: ["氧化物", "易挥发物"], color: "#c8d0d8",
  });
  S("NO2", "NO2", { N: 1, O: 2 }, 0, "g", 51.3, 4, "N", "二氧化氮", {
    cats: ["氧化物", "易挥发物"], color: "#c08040",
  });

  /* ================= C ================= */
  S("CO2", "CO2", { C: 1, O: 2 }, 0, "aq", -386.0, 4, "C", "二氧化碳（溶解）", {
    cats: ["氧化物", "酸碱", "易挥发物"], color: "#e0e0e0",
  });
  S("CO2g", "CO2", { C: 1, O: 2 }, 0, "g", -394.4, 4, "C", "二氧化碳（气）", {
    cats: ["氧化物", "易挥发物"], color: "#d8d8d8", aqPair: "CO2",
  });
  S("H2CO3", "H2CO3", { H: 2, C: 1, O: 3 }, 0, "aq", -614.84, 4, "C", "碳酸（痕量中间体）", {
    cats: ["酸碱"], color: "#e8e8e8", derived: true,
    note: "由 CO₂(aq) + H₂O ⇌ H₂CO₃ 的 K≈0.035 推得；不作为碳酸主形态",
  });
  S("HCO3-", "HCO3^-", { H: 1, C: 1, O: 3 }, -1, "aq", -586.8, 4, "C", "碳酸氢根", {
    cats: ["离子", "酸碱"], color: "#e0e0d8",
  });
  S("CO32-", "CO3^2-", { C: 1, O: 3 }, -2, "aq", -527.8, 4, "C", "碳酸根", {
    cats: ["离子", "酸碱"], color: "#d8d8d0",
  });
  S("CO", "CO", { C: 1, O: 1 }, 0, "g", -137.2, 2, "C", "一氧化碳", {
    cats: ["氧化物", "易挥发物"], color: "#d0d0d0",
  });
  S("HCOOH", "HCOOH", { H: 2, C: 1, O: 2 }, 0, "aq", -372.3, 2, "C", "甲酸", {
    cats: ["酸碱"], color: "#e8e0d0",
  });
  S("HCOO-", "HCOO^-", { H: 1, C: 1, O: 2 }, -1, "aq", -351.0, 2, "C", "甲酸根", {
    cats: ["离子", "酸碱"], color: "#e0d8c8",
  });
  S("CH4", "CH4", { H: 4, C: 1 }, 0, "g", -50.5, -4, "C", "甲烷", {
    cats: ["易挥发物"], color: "#e0e0e0",
  });
  S("CH3OH", "CH3OH", { H: 4, C: 1, O: 1 }, 0, "aq", -175.3, -2, "C", "甲醇（溶解）", {
    cats: ["易挥发物"], color: "#e8e8e8",
  });
  S("CH3OHg", "CH3OH", { H: 4, C: 1, O: 1 }, 0, "g", -162.3, -2, "C", "甲醇（气）", {
    cats: ["易挥发物"], color: "#e0e0e0", aqPair: "CH3OH",
  });

  /* ================= V ================= */
  S("V", "V", { V: 1 }, 0, "s", 0, 0, "V", "钒", { cats: ["单质"], color: "#b8bcc0" });
  S("V2+", "V^2+", { V: 1 }, 2, "aq", -218.0, 2, "V", "二价钒离子（紫）", {
    cats: ["离子"], color: "#8a5ac0",
  });
  S("V(OH)2", "V(OH)2", { V: 1, O: 2, H: 2 }, 0, "s", -612.3, 2, "V", "氢氧化钒(II)↓", {
    cats: ["难溶物", "酸碱"], color: "#9a7ac0", derived: true, note: "近似值",
  });
  S("V3+", "V^3+", { V: 1 }, 3, "aq", -243.1, 3, "V", "三价钒离子（绿）", {
    cats: ["离子"], color: "#4a9a5a",
  });
  S("V(OH)3", "V(OH)3", { V: 1, O: 3, H: 3 }, 0, "s", -878.0, 3, "V", "氢氧化钒(III)↓", {
    cats: ["难溶物", "酸碱"], color: "#5a9a6a", derived: true, note: "近似值",
  });
  S("V2O3", "V2O3", { V: 2, O: 3 }, 0, "s", -1139.0, 3, "V", "三氧化二钒", {
    cats: ["氧化物", "酸碱"], color: "#4a6a5a",
  });
  S("VO2", "VO2", { V: 1, O: 2 }, 0, "s", -736.0, 4, "V", "二氧化钒", {
    cats: ["氧化物", "酸碱"], color: "#4a5a8a", derived: true, note: "近似值",
  });
  S("VO2+", "VO^2+", { V: 1, O: 1 }, 2, "aq", -447.4, 4, "V", "钒酰离子（蓝）", {
    cats: ["离子"], color: "#3a6ac0",
  });
  S("VO2plus", "VO2^+", { V: 1, O: 2 }, 1, "aq", -588.1, 5, "V", "二氧钒离子（黄）", {
    cats: ["离子"], color: "#d0c040",
  });
  S("V2O5", "V2O5", { V: 2, O: 5 }, 0, "s", -1419.0, 5, "V", "五氧化二钒", {
    cats: ["氧化物", "酸碱"], color: "#d08030",
  });
  S("VO3-", "VO3^-", { V: 1, O: 3 }, -1, "aq", -899.0, 5, "V", "偏钒酸根", {
    cats: ["离子"], color: "#d0a040", derived: true, note: "近似；多聚钒酸根未细分",
  });
  S("VO43-", "VO4^3-", { V: 1, O: 4 }, -3, "aq", -1018.0, 5, "V", "正钒酸根", {
    cats: ["离子"], color: "#d0b050", derived: true, note: "近似",
  });

  /* ================= Cr ================= */
  S("Cr", "Cr", { Cr: 1 }, 0, "s", 0, 0, "Cr", "铬", { cats: ["单质"], color: "#b0b4b8" });
  S("Cr3+", "Cr^3+", { Cr: 1 }, 3, "aq", -215.5, 3, "Cr", "三价铬离子（绿/紫）", {
    cats: ["离子"], color: "#5a8a5a",
  });
  S("Cr2O3", "Cr2O3", { Cr: 2, O: 3 }, 0, "s", -1058.1, 3, "Cr", "氧化铬（绿）", {
    cats: ["氧化物", "酸碱"], color: "#3a7a3a",
  });
  S("Cr(OH)3", "Cr(OH)3", { Cr: 1, O: 3, H: 3 }, 0, "s", -859.6, 3, "Cr", "氢氧化铬↓（新鲜沉淀）", {
    cats: ["难溶物", "酸碱"], color: "#6a9a7a", derived: true,
    note: "由 Ksp(Cr(OH)₃)≈6.3×10⁻³¹ 推得；相对 Cr₂O₃ 为亚稳相",
  });
  S("Cr(OH)4-", "Cr(OH)4^-", { Cr: 1, O: 4, H: 4 }, -1, "aq", -1019.1, 3, "Cr", "四羟基合铬酸根", {
    cats: ["离子", "络合物", "酸碱"], color: "#5a9a6a", derived: true,
    note: "由 Cr(OH)₃(s) + OH⁻ ⇌ Cr(OH)₄⁻ 的 logK≈-0.4 推得",
  });
  S("CrO42-", "CrO4^2-", { Cr: 1, O: 4 }, -2, "aq", -727.8, 6, "Cr", "铬酸根（黄）", {
    cats: ["离子"], color: "#e8d020",
  });
  S("Cr2O72-", "Cr2O7^2-", { Cr: 2, O: 7 }, -2, "aq", -1301.1, 6, "Cr", "重铬酸根（橙）", {
    cats: ["离子"], color: "#e87020",
  });

  /* ================= Mn ================= */
  S("Mn", "Mn", { Mn: 1 }, 0, "s", 0, 0, "Mn", "锰", { cats: ["单质"], color: "#a8acb0" });
  S("Mn2+", "Mn^2+", { Mn: 1 }, 2, "aq", -228.1, 2, "Mn", "二价锰离子（淡粉）", {
    cats: ["离子"], color: "#f0d8e0",
  });
  S("Mn(OH)2", "Mn(OH)2", { Mn: 1, O: 2, H: 2 }, 0, "s", -615.0, 2, "Mn", "氢氧化锰↓（白）", {
    cats: ["难溶物", "酸碱"], color: "#f0ece8",
  });
  S("MnO", "MnO", { Mn: 1, O: 1 }, 0, "s", -362.9, 2, "Mn", "氧化锰(II)", {
    cats: ["氧化物", "酸碱"], color: "#6a8a5a",
  });
  S("MnO(OH)", "MnO(OH)", { Mn: 1, O: 2, H: 1 }, 0, "s", -557.0, 3, "Mn", "羟基氧化锰（水锰矿）", {
    cats: ["氧化物", "难溶物", "酸碱"], color: "#6a4a3a", derived: true, note: "近似值",
  });
  S("Mn2O3", "Mn2O3", { Mn: 2, O: 3 }, 0, "s", -881.2, 3, "Mn", "三氧化二锰", {
    cats: ["氧化物", "酸碱"], color: "#5a3a2a",
  });
  S("MnO2", "MnO2", { Mn: 1, O: 2 }, 0, "s", -465.1, 4, "Mn", "二氧化锰（黑）", {
    cats: ["氧化物", "酸碱"], color: "#2a2a2a",
  });
  S("MnO42-", "MnO4^2-", { Mn: 1, O: 4 }, -2, "aq", -500.8, 6, "Mn", "锰酸根（绿）", {
    cats: ["离子"], color: "#3a9a5a",
  });
  S("MnO4-", "MnO4^-", { Mn: 1, O: 4 }, -1, "aq", -447.2, 7, "Mn", "高锰酸根（紫）", {
    cats: ["离子"], color: "#8a2a8a",
  });

  /* ================= Fe ================= */
  S("Fe", "Fe", { Fe: 1 }, 0, "s", 0, 0, "Fe", "铁", { cats: ["单质"], color: "#9a9ca0" });
  S("Fe2+", "Fe^2+", { Fe: 1 }, 2, "aq", -78.9, 2, "Fe", "亚铁离子（浅绿）", {
    cats: ["离子"], color: "#a8d0a8",
  });
  S("Fe(OH)2", "Fe(OH)2", { Fe: 1, O: 2, H: 2 }, 0, "s", -494.89, 2, "Fe", "氢氧化亚铁↓（白）", {
    cats: ["难溶物", "酸碱"], color: "#e8f0e8",
  });
  S("FeO", "FeO", { Fe: 1, O: 1 }, 0, "s", -244.59, 2, "Fe", "氧化亚铁（黑）", {
    cats: ["氧化物", "酸碱"], color: "#2a2a2e",
  });
  S("Fe3O4", "Fe3O4", { Fe: 3, O: 4 }, 0, "s", -1015.4, 8 / 3, "Fe", "四氧化三铁（磁铁矿，黑）", {
    cats: ["氧化物", "酸碱"], color: "#1a1a1e",
  });
  S("Fe3+", "Fe^3+", { Fe: 1 }, 3, "aq", -4.7, 3, "Fe", "铁离子（黄）", {
    cats: ["离子"], color: "#d0a040",
  });
  S("Fe(OH)3", "Fe(OH)3", { Fe: 1, O: 3, H: 3 }, 0, "s", -696.5, 3, "Fe", "氢氧化铁↓（红褐）", {
    cats: ["难溶物", "酸碱"], color: "#8a4a2a",
  });
  S("FeO(OH)", "FeO(OH)", { Fe: 1, O: 2, H: 1 }, 0, "s", -488.6, 3, "Fe", "羟基氧化铁（针铁矿）", {
    cats: ["难溶物", "酸碱"], color: "#9a5a2a", derived: true, note: "近似值",
  });
  S("Fe2O3", "Fe2O3", { Fe: 2, O: 3 }, 0, "s", -742.2, 3, "Fe", "氧化铁（赤铁矿，红棕）", {
    cats: ["氧化物", "酸碱"], color: "#9a3a2a",
  });
  S("FeO42-", "FeO4^2-", { Fe: 1, O: 4 }, -2, "aq", -316.5, 6, "Fe", "高铁酸根（紫红）", {
    cats: ["离子"], color: "#a03060", derived: true,
    note: "由 FeO₄²⁻/Fe³⁺ 的 E°≈2.20 V 推得",
  });

  /* ================= Cu ================= */
  S("Cu", "Cu", { Cu: 1 }, 0, "s", 0, 0, "Cu", "铜（紫红）", { cats: ["单质"], color: "#c07040" });
  S("Cu+", "Cu^+", { Cu: 1 }, 1, "aq", 50.0, 1, "Cu", "亚铜离子", {
    cats: ["离子"], color: "#d0a080",
  });
  S("Cu2+", "Cu^2+", { Cu: 1 }, 2, "aq", 65.5, 2, "Cu", "铜离子（蓝）", {
    cats: ["离子"], color: "#3a8ad0",
  });
  S("Cu2O", "Cu2O", { Cu: 2, O: 1 }, 0, "s", -146.0, 1, "Cu", "氧化亚铜（红）", {
    cats: ["氧化物", "酸碱"], color: "#b03020",
  });
  S("CuO", "CuO", { Cu: 1, O: 1 }, 0, "s", -129.7, 2, "Cu", "氧化铜（黑）", {
    cats: ["氧化物", "酸碱"], color: "#1a1a1e",
  });
  S("Cu(OH)2", "Cu(OH)2", { Cu: 1, O: 2, H: 2 }, 0, "s", -359.5, 2, "Cu", "氢氧化铜↓（蓝）", {
    cats: ["难溶物", "酸碱"], color: "#4a9ad0",
  });
  S("Cu2S", "Cu2S", { Cu: 2, S: 1 }, 0, "s", -86.2, 1, "Cu", "硫化亚铜（辉铜矿，黑）", {
    cats: ["难溶物"], color: "#2a2a2e", ligs: ["S2-"],
  });
  S("Cu2S2", "Cu2S2", { Cu: 2, S: 2 }, 0, "s", -107.2, 1, "Cu", "二硫化二铜", {
    cats: ["难溶物"], color: "#2a2a2e", ligs: ["S2-"], derived: true, note: "近似值",
  });
  S("CuFeS2", "CuFeS2", { Cu: 1, Fe: 1, S: 2 }, 0, "s", -190.6, 1, "Cu", "黄铜矿", {
    cats: ["难溶物"], color: "#c8a040", derived: true, note: "近似值",
  });
  S("Cu2(OH)2CO3", "Cu2(OH)2CO3", { Cu: 2, O: 5, H: 2, C: 1 }, 0, "s", -1051.0, 2, "Cu", "碱式碳酸铜（铜绿）", {
    cats: ["难溶物", "酸碱"], color: "#3a9a7a", ligs: ["CO32-"], derived: true, note: "近似值",
  });

  /* ================= Ag ================= */
  S("Ag", "Ag", { Ag: 1 }, 0, "s", 0, 0, "Ag", "银（白）", { cats: ["单质"], color: "#d0d0d8" });
  S("Ag+", "Ag^+", { Ag: 1 }, 1, "aq", 77.1, 1, "Ag", "银离子", {
    cats: ["离子"], color: "#d8d8e0",
  });
  S("Ag2O", "Ag2O", { Ag: 2, O: 1 }, 0, "s", -11.2, 1, "Ag", "氧化银（褐黑）", {
    cats: ["氧化物", "酸碱"], color: "#3a2a2a",
  });

  /* ================= Pb ================= */
  S("Pb", "Pb", { Pb: 1 }, 0, "s", 0, 0, "Pb", "铅（灰蓝）", { cats: ["单质"], color: "#8a8c94" });
  S("Pb2+", "Pb^2+", { Pb: 1 }, 2, "aq", -24.4, 2, "Pb", "铅离子", {
    cats: ["离子"], color: "#b0b8c0",
  });
  S("PbO", "PbO", { Pb: 1, O: 1 }, 0, "s", -187.9, 2, "Pb", "氧化铅（黄）", {
    cats: ["氧化物", "酸碱"], color: "#d0a030",
  });
  S("PbO2", "PbO2", { Pb: 1, O: 2 }, 0, "s", -217.3, 4, "Pb", "二氧化铅（褐）", {
    cats: ["氧化物", "酸碱"], color: "#4a3020",
  });

  /* ================= 环境离子 / 配体 ================= */
  S("Ca2+", "Ca^2+", { Ca: 1 }, 2, "aq", -553.6, 2, "Ca", "钙离子", {
    cats: ["离子"], color: "#d8d8e0",
  });
  S("Ba2+", "Ba^2+", { Ba: 1 }, 2, "aq", -560.8, 2, "Ba", "钡离子", {
    cats: ["离子"], color: "#d8d8e0",
  });
  S("Br-", "Br^-", { Br: 1 }, -1, "aq", -104.0, -1, "Br", "溴离子", {
    cats: ["离子"], color: "#e0d0b0",
  });
  S("F-", "F^-", { F: 1 }, -1, "aq", -278.8, -1, "F", "氟离子", {
    cats: ["离子", "酸碱"], color: "#d8e8d8",
  });
  S("HF", "HF", { H: 1, F: 1 }, 0, "aq", -296.8, -1, "F", "氢氟酸", {
    cats: ["酸碱"], color: "#e0e8e0",
  });
  S("CN-", "CN^-", { C: 1, N: 1 }, -1, "aq", 172.4, 2, "C", "氰根", {
    cats: ["离子", "络合物", "酸碱"], color: "#e0d8d8",
  });
  S("HCN", "HCN", { H: 1, C: 1, N: 1 }, 0, "aq", 119.8, 2, "C", "氢氰酸", {
    cats: ["酸碱", "易挥发物"], color: "#e8d8d8",
  });
  S("SCN-", "SCN^-", { S: 1, C: 1, N: 1 }, -1, "aq", 92.7, -2, "S", "硫氰根", {
    cats: ["离子", "络合物"], color: "#e8d8c0",
  });
  S("CH3COO-", "CH3COO^-", { C: 2, H: 3, O: 2 }, -1, "aq", -369.4, 0, "C", "醋酸根", {
    cats: ["离子", "酸碱"], color: "#e8e0d8",
  });
  S("CH3COOH", "CH3COOH", { C: 2, H: 4, O: 2 }, 0, "aq", -396.5, 0, "C", "醋酸", {
    cats: ["酸碱"], color: "#eee6de",
  });
  S("C2O42-", "C2O4^2-", { C: 2, O: 4 }, -2, "aq", -674.0, 3, "C", "草酸根", {
    cats: ["离子", "酸碱", "络合物"], color: "#e8e0d0",
  });
  S("HC2O4-", "HC2O4^-", { H: 1, C: 2, O: 4 }, -1, "aq", -697.63, 3, "C", "草酸氢根", {
    cats: ["离子", "酸碱"], color: "#ece4d4", derived: true,
    note: "由 pKa₂(草酸)=4.14 与 C₂O₄²⁻ 的 ΔGf° 推得",
  });
  S("H2C2O4", "H2C2O4", { H: 2, C: 2, O: 4 }, 0, "aq", -704.77, 3, "C", "草酸", {
    cats: ["酸碱"], color: "#f0e8d8", derived: true,
    note: "由 pKa₁(草酸)=1.25 推得",
  });
  S("Zn2+", "Zn^2+", { Zn: 1 }, 2, "aq", -147.1, 2, "Zn", "锌离子", {
    cats: ["离子"], color: "#c8ccd0",
  });
  S("Zn", "Zn", { Zn: 1 }, 0, "s", 0, 0, "Zn", "锌", { cats: ["单质"], color: "#b0b4b8" });
  S("ZnO", "ZnO", { Zn: 1, O: 1 }, 0, "s", -318.3, 2, "Zn", "氧化锌（锌华）", {
    cats: ["氧化物", "酸碱"], color: "#e0e0d8",
  });
  S("Zn(OH)2", "Zn(OH)2", { Zn: 1, O: 2, H: 2 }, 0, "s", -553.5, 2, "Zn", "氢氧化锌↓（白）", {
    cats: ["难溶物", "酸碱"], color: "#f2f2ee",
    note: "取无定形/活性值；晶质相略低，属数据源差异",
  });

  /* ================= 配合物（按 logβ 推得） ================= */
  Cx("AlF63-", "AlF6^3-", { Al: 1, F: 6 }, -3, 3, "Al", "六氟合铝酸根（冰晶石）",
    "Al3+", [["F-", 6]], 19.7, ["络合物"], { color: "#d8e8e0" });
  Cx("FeSCN2+", "FeSCN^2+", { Fe: 1, S: 1, C: 1, N: 1 }, 2, 3, "Fe",
    "硫氰合铁(III)（血红）", "Fe3+", [["SCN-", 1]], 2.3, ["络合物"], { color: "#c02020" });
  Cx("Fe(CN)64-", "Fe(CN)6^4-", { Fe: 1, C: 6, N: 6 }, -4, 2, "Fe",
    "亚铁氰根（黄血盐）", "Fe2+", [["CN-", 6]], 35.4, ["络合物"], { color: "#e8d060" });
  Cx("Fe(CN)63-", "Fe(CN)6^3-", { Fe: 1, C: 6, N: 6 }, -3, 3, "Fe",
    "铁氰根（赤血盐）", "Fe3+", [["CN-", 6]], 43.0, ["络合物"], { color: "#c04030" });
  Cx("CuCl2-", "CuCl2^-", { Cu: 1, Cl: 2 }, -1, 1, "Cu",
    "二氯合铜(I)", "Cu+", [["Cl-", 2]], 5.7, ["络合物"], { color: "#b0b070" });
  Cx("CuCl42-", "CuCl4^2-", { Cu: 1, Cl: 4 }, -2, 2, "Cu",
    "四氯合铜(II)（黄绿）", "Cu2+", [["Cl-", 4]], 5.6, ["络合物"], { color: "#7aa030" });
  Cx("Cu(NH3)2+", "Cu(NH3)2^+", { Cu: 1, N: 2, H: 6 }, 1, 1, "Cu",
    "二氨合铜(I)", "Cu+", [["NH3", 2]], 10.8, ["络合物"], { color: "#6a8ac0" });
  Cx("Cu(NH3)42+", "Cu(NH3)4^2+", { Cu: 1, N: 4, H: 12 }, 2, 2, "Cu",
    "四氨合铜(II)（深蓝）", "Cu2+", [["NH3", 4]], 13.32, ["络合物"], { color: "#1a3ac0" });
  Cx("Cu(OH)42-", "Cu(OH)4^2-", { Cu: 1, O: 4, H: 4 }, -2, 2, "Cu",
    "四羟基合铜酸根（蓝）", "Cu2+", [["OH-", 4]], 15.5, ["络合物", "酸碱"], { color: "#3a7ac0" });
  Cx("AgCl2-", "AgCl2^-", { Ag: 1, Cl: 2 }, -1, 1, "Ag",
    "二氯合银(I)", "Ag+", [["Cl-", 2]], 5.25, ["络合物"], { color: "#c8c8c0" });
  Cx("Ag(NH3)2+", "Ag(NH3)2^+", { Ag: 1, N: 2, H: 6 }, 1, 1, "Ag",
    "二氨合银（银氨溶液）", "Ag+", [["NH3", 2]], 7.22, ["络合物"], { color: "#c8d8c0" });
  Cx("Ag(CN)2-", "Ag(CN)2^-", { Ag: 1, C: 2, N: 2 }, -1, 1, "Ag",
    "二氰合银(I)", "Ag+", [["CN-", 2]], 21.0, ["络合物"], { color: "#c8c0d8" });
  Cx("Ag(S2O3)23-", "Ag(S2O3)2^3-", { Ag: 1, S: 4, O: 6 }, -3, 1, "Ag",
    "二硫代硫酸合银(I)", "Ag+", [["S2O32-", 2]], 13.5, ["络合物"], { color: "#d0d0b0" });
  Cx("Zn(OH)42-", "Zn(OH)4^2-", { Zn: 1, O: 4, H: 4 }, -2, 2, "Zn",
    "四羟基合锌酸根", "Zn2+", [["OH-", 4]], 15.5, ["络合物", "酸碱"], { color: "#c8d8e0" });
  Cx("Zn(NH3)42+", "Zn(NH3)4^2+", { Zn: 1, N: 4, H: 12 }, 2, 2, "Zn",
    "四氨合锌(II)", "Zn2+", [["NH3", 4]], 9.1, ["络合物"], { color: "#c8d0e8" });
  Cx("Zn(CN)42-", "Zn(CN)4^2-", { Zn: 1, C: 4, N: 4 }, -2, 2, "Zn",
    "四氰合锌(II)", "Zn2+", [["CN-", 4]], 16.7, ["络合物"], { color: "#d8d0e0" });
  Cx("ZnCl42-", "ZnCl4^2-", { Zn: 1, Cl: 4 }, -2, 2, "Zn",
    "四氯合锌(II)", "Zn2+", [["Cl-", 4]], 0.2, ["络合物"], { color: "#d0d8d8" });

  /* ================= 沉淀（按 Ksp 推得） ================= */
  Ppt("AgCl", "AgCl", { Ag: 1, Cl: 1 }, 0, 1, "Ag", "氯化银↓（白）",
    [["Ag+", 1], ["Cl-", 1]], -9.75, ["难溶物"], {
      color: "#f4f4f0", ligs: ["Cl-"], oxByEl: { Ag: 1, Cl: -1 },
    });
  Ppt("AgBr", "AgBr", { Ag: 1, Br: 1 }, 0, 1, "Ag", "溴化银↓（淡黄）",
    [["Ag+", 1], ["Br-", 1]], -12.30, ["难溶物"], {
      color: "#e8e0a0", ligs: ["Br-"], oxByEl: { Ag: 1, Br: -1 },
    });
  Ppt("AgI", "AgI", { Ag: 1, I: 1 }, 0, 1, "Ag", "碘化银↓（黄）",
    [["Ag+", 1], ["I-", 1]], -16.08, ["难溶物"], {
      color: "#e8e070", ligs: ["I-"], oxByEl: { Ag: 1, I: -1 },
    });
  Ppt("Ag2S", "Ag2S", { Ag: 2, S: 1 }, 0, 1, "Ag", "硫化银↓（黑）",
    [["Ag+", 2], ["S2-", 1]], -49.2, ["难溶物"], {
      color: "#1a1a1e", ligs: ["S2-"], oxByEl: { Ag: 1, S: -2 },
    });
  Ppt("CuCl", "CuCl", { Cu: 1, Cl: 1 }, 0, 1, "Cu", "氯化亚铜↓（白）",
    [["Cu+", 1], ["Cl-", 1]], -6.77, ["难溶物"], {
      color: "#ecece8", ligs: ["Cl-"],
    });
  Ppt("CuI", "CuI", { Cu: 1, I: 1 }, 0, 1, "Cu", "碘化亚铜↓（白）",
    [["Cu+", 1], ["I-", 1]], -11.96, ["难溶物"], {
      color: "#ece8e0", ligs: ["I-"],
    });
  Ppt("PbI2", "PbI2", { Pb: 1, I: 2 }, 0, 2, "Pb", "碘化铅↓（金黄）",
    [["Pb2+", 1], ["I-", 2]], -8.0, ["难溶物"], {
      color: "#e8c840", ligs: ["I-"], oxByEl: { Pb: 2, I: -1 },
    });
  Ppt("PbSO4", "PbSO4", { Pb: 1, S: 1, O: 4 }, 0, 2, "Pb", "硫酸铅↓（白）",
    [["Pb2+", 1], ["SO42-", 1]], -7.8, ["难溶物"], {
      color: "#f0f0f0", ligs: ["SO42-"], oxByEl: { Pb: 2, S: 6 },
    });
  Ppt("PbCl2", "PbCl2", { Pb: 1, Cl: 2 }, 0, 2, "Pb", "氯化铅↓（白）",
    [["Pb2+", 1], ["Cl-", 2]], -4.8, ["难溶物"], {
      color: "#f0f0f0", ligs: ["Cl-"], oxByEl: { Pb: 2, Cl: -1 },
    });
  Ppt("Pb(OH)2", "Pb(OH)2", { Pb: 1, O: 2, H: 2 }, 0, 2, "Pb", "氢氧化铅↓（白）",
    [["Pb2+", 1], ["OH-", 2]], -19.85, ["难溶物", "酸碱"], {
      color: "#ececec", ligs: ["OH-"],
    });
  Ppt("PbCO3", "PbCO3", { Pb: 1, C: 1, O: 3 }, 0, 2, "Pb", "碳酸铅（白铅矿）↓",
    [["Pb2+", 1], ["CO32-", 1]], -13.1, ["难溶物"], {
      color: "#f0f0e8", ligs: ["CO32-"],
    });
  Ppt("PbS", "PbS", { Pb: 1, S: 1 }, 0, 2, "Pb", "硫化铅（方铅矿）↓",
    [["Pb2+", 1], ["S2-", 1]], -27.5, ["难溶物"], {
      color: "#2a2a30", ligs: ["S2-"], oxByEl: { Pb: 2, S: -2 },
    });
  Ppt("MgCO3", "MgCO3", { Mg: 1, C: 1, O: 3 }, 0, 2, "Mg", "碳酸镁（菱镁矿）↓",
    [["Mg2+", 1], ["CO32-", 1]], -5.17, ["难溶物"], {
      color: "#f4f4f0", ligs: ["CO32-"],
    });
  Ppt("FeCO3", "FeCO3", { Fe: 1, C: 1, O: 3 }, 0, 2, "Fe", "碳酸亚铁（菱铁矿）↓",
    [["Fe2+", 1], ["CO32-", 1]], -10.24, ["难溶物"], {
      color: "#d8c8a8", ligs: ["CO32-"],
    });
  Ppt("ZnCO3", "ZnCO3", { Zn: 1, C: 1, O: 3 }, 0, 2, "Zn", "碳酸锌（菱锌矿）↓",
    [["Zn2+", 1], ["CO32-", 1]], -10.0, ["难溶物"], {
      color: "#e8e8d8", ligs: ["CO32-"],
    });
  Ppt("ZnS", "ZnS", { Zn: 1, S: 1 }, 0, 2, "Zn", "硫化锌（闪锌矿）↓",
    [["Zn2+", 1], ["S2-", 1]], -24.7, ["难溶物"], {
      color: "#d8d0a0", ligs: ["S2-"], oxByEl: { Zn: 2, S: -2 },
    });
  Ppt("MnCO3", "MnCO3", { Mn: 1, C: 1, O: 3 }, 0, 2, "Mn", "碳酸锰↓",
    [["Mn2+", 1], ["CO32-", 1]], -10.66, ["难溶物"], {
      color: "#e8d0d0", ligs: ["CO32-"],
    });
  Ppt("CaMnO4", "CaMnO4", { Ca: 1, Mn: 1, O: 4 }, 0, 6, "Mn", "锰酸钙↓",
    [["Ca2+", 1], ["MnO42-", 1]], -4.0, ["难溶物"], {
      color: "#5a8a6a", ligs: ["Ca2+"], oxByEl: { Ca: 2, Mn: 6 },
      note: "近似：logKsp≈-4；多聚/水合形态未细分",
    });
  S("FeS", "FeS", { Fe: 1, S: 1 }, 0, "s", -100.4, 2, "Fe", "硫化亚铁↓（黑）", {
    cats: ["难溶物"], color: "#2a2a2a", ligs: ["S2-"],
  });
  S("FeS2", "FeS2", { Fe: 1, S: 2 }, 0, "s", -166.9, 2, "Fe", "二硫化亚铁（黄铁矿）", {
    cats: ["难溶物"], color: "#c8a840", ligs: ["S2-"], derived: true,
    note: "NBS 值；黄铁矿中硫以 S₂²⁻ 形式存在，未进一步细分",
  });
  Ppt("CaCO3", "CaCO3", { Ca: 1, C: 1, O: 3 }, 0, 2, "Ca", "碳酸钙（石灰石）↓",
    [["Ca2+", 1], ["CO32-", 1]], -8.48, ["难溶物"], {
      color: "#f4f4f0", ligs: ["CO32-"],
    });
  Ppt("BaCO3", "BaCO3", { Ba: 1, C: 1, O: 3 }, 0, 2, "Ba", "碳酸钡↓",
    [["Ba2+", 1], ["CO32-", 1]], -8.29, ["难溶物"], {
      color: "#f4f4f0", ligs: ["CO32-"],
    });
  Ppt("BaSO4", "BaSO4", { Ba: 1, S: 1, O: 4 }, 0, 6, "S", "硫酸钡（钡餐）↓",
    [["Ba2+", 1], ["SO42-", 1]], -9.96, ["难溶物"], {
      color: "#f6f6f6", ligs: ["Ba2+"], oxByEl: { S: 6, Ba: 2 },
    });

  /* ================= 酸碱家族（最质子化 → 参考态） ================= */
  // hydrate: 第一步为 A0 + H2O ⇌ A1 + H+（CO2、SO2 等）
  const FAMILIES = {
    water:     { list: ["H2O", "OH-"], hydrate: false, ref: "OH-" },
    carbonate: { list: ["CO2", "HCO3-", "CO32-"], hydrate: true, ref: "CO32-" },
    ammonia:   { list: ["NH4+", "NH3"], hydrate: false, ref: "NH3" },
    acetate:   { list: ["CH3COOH", "CH3COO-"], hydrate: false, ref: "CH3COO-" },
    sulfide:   { list: ["H2S", "HS-", "S2-"], hydrate: false, ref: "S2-" },
    sulfite:   { list: ["SO2", "HSO3-", "SO32-"], hydrate: true, ref: "SO32-" },
    sulfate:   { list: ["HSO4-", "SO42-"], hydrate: false, ref: "SO42-" },
    phosphate: { list: ["H3PO4", "H2PO4-", "HPO42-", "PO43-"], hydrate: false, ref: "PO43-" },
    phosphite: { list: ["H3PO3", "H2PO3-", "HPO32-"], hydrate: false, ref: "HPO32-" },
    hypophosphite: { list: ["H3PO2", "H2PO2-"], hydrate: false, ref: "H2PO2-" },
    fluoride:  { list: ["HF", "F-"], hydrate: false, ref: "F-" },
    cyanide:   { list: ["HCN", "CN-"], hydrate: false, ref: "CN-" },
    oxalate:   { list: ["H2C2O4", "HC2O4-", "C2O42-"], hydrate: false, ref: "C2O42-" },
    hypochlorite: { list: ["HClO", "ClO-"], hydrate: false, ref: "ClO-" },
    nitrite:   { list: ["HNO2", "NO2-"], hydrate: false, ref: "NO2-" },
    formate:   { list: ["HCOOH", "HCOO-"], hydrate: false, ref: "HCOO-" },
  };

  // 参考配体 → 家族 key
  const LIG_FAMILY = {
    "OH-": "water",
    "CO32-": "carbonate",
    "HCO3-": "carbonate",
    "NH3": "ammonia",
    "CH3COO-": "acetate",
    "S2-": "sulfide",
    "HS-": "sulfide",
    "SO32-": "sulfite",
    "HSO3-": "sulfite",
    "SO42-": "sulfate",
    "HPO42-": "phosphate",
    "H2PO4-": "phosphate",
    "PO43-": "phosphate",
    "HPO32-": "phosphite",
    "H2PO3-": "phosphite",
    "H2PO2-": "hypophosphite",
    "F-": "fluoride",
    "CN-": "cyanide",
    "C2O42-": "oxalate",
    "HC2O4-": "oxalate",
    "ClO-": "hypochlorite",
    "NO2-": "nitrite",
    "HCOO-": "formate",
  };

  /* ================= 化学环境物料 ================= */
  // groups: acidbase / precipitant / ligand / redox
  const ENV = [
    { id: "OH-", label: "OH⁻", formula: "OH^-", groups: ["acidbase", "ligand"], family: "water",
      defaultC: null, special: "pH", note: "活度由 pH 决定" },
    { id: "CO32-", label: "CO₃²⁻", formula: "CO3^2-", groups: ["acidbase", "precipitant"], family: "carbonate",
      defaultC: 0.01, note: "总碳酸浓度；随 pH 质子化" },
    { id: "NH3", label: "NH₃", formula: "NH3", groups: ["acidbase", "ligand"], family: "ammonia",
      defaultC: 0.01, note: "总氨浓度；随 pH 质子化" },
    { id: "CH3COOH", label: "CH₃COOH", formula: "CH3COOH", groups: ["acidbase"], family: "acetate",
      defaultC: 0.01, note: "总醋酸浓度" },
    { id: "S2-", label: "S²⁻", formula: "S^2-", groups: ["precipitant"], family: "sulfide",
      defaultC: 0.001, note: "总硫化物浓度；随 pH 质子化" },
    { id: "Ca2+", label: "Ca²⁺", formula: "Ca^2+", groups: ["precipitant"], defaultC: 0.01 },
    { id: "Ba2+", label: "Ba²⁺", formula: "Ba^2+", groups: ["precipitant"], defaultC: 0.01 },
    { id: "Ag+", label: "Ag⁺", formula: "Ag^+", groups: ["precipitant"], defaultC: 0.001 },
    { id: "Cl-", label: "Cl⁻", formula: "Cl^-", groups: ["ligand"], defaultC: 0.1 },
    { id: "CN-", label: "CN⁻", formula: "CN^-", groups: ["ligand"], family: "cyanide", defaultC: 0.001 },
    { id: "SCN-", label: "SCN⁻", formula: "SCN^-", groups: ["ligand"], defaultC: 0.001 },
    { id: "C2O42-", label: "C₂O₄²⁻", formula: "C2O4^2-", groups: ["ligand"], family: "oxalate", defaultC: 0.01 },
    { id: "S2O32-", label: "S₂O₃²⁻", formula: "S2O3^2-", groups: ["ligand"], defaultC: 0.01 },
    { id: "F-", label: "F⁻", formula: "F^-", groups: ["ligand"], family: "fluoride", defaultC: 0.01 },
    { id: "I-", label: "I⁻", formula: "I^-", groups: ["ligand"], defaultC: 0.001, note: "AgI/PbI₂ 等沉淀剂" },
    { id: "Br-", label: "Br⁻", formula: "Br^-", groups: ["ligand"], defaultC: 0.001, note: "AgBr 沉淀剂" },
    { id: "SO42-", label: "SO₄²⁻", formula: "SO4^2-", groups: ["ligand"], family: "sulfate", defaultC: 0.01, note: "PbSO₄/BaSO₄ 沉淀剂" },
    { id: "I2/I-", label: "I₂/I⁻", formula: "I2", groups: ["redox"], redox: { ox: "I2", red: "I-", n: 2 },
      defaultC: 0.01, ratio: 1 },
    { id: "Fe3+/Fe2+", label: "Fe³⁺/Fe²⁺", formula: "Fe^3+", groups: ["redox"],
      redox: { ox: "Fe3+", red: "Fe2+", n: 1 }, defaultC: 0.01, ratio: 1 },
    { id: "Cr2O72-/Cr3+", label: "Cr₂O₇²⁻/Cr³⁺", formula: "Cr2O7^2-", groups: ["redox"],
      redox: { ox: "Cr2O72-", red: "Cr3+", n: 6 }, defaultC: 0.001, ratio: 1 },
    { id: "MnO2/Mn2+", label: "MnO₂/Mn²⁺", formula: "MnO2", groups: ["redox"],
      redox: { ox: "MnO2", red: "Mn2+", n: 2, solidOx: true }, defaultC: 0.01, ratio: 1 },
    { id: "Zn2+/Zn", label: "Zn²⁺/Zn", formula: "Zn^2+", groups: ["redox"],
      redox: { ox: "Zn2+", red: "Zn", n: 2, solidRed: true }, defaultC: 0.01, ratio: 1 },
    { id: "C2O42-/CO2", label: "C₂O₄²⁻/CO₂", formula: "C2O4^2-", groups: ["redox"],
      redox: { ox: "C2O42-", red: "CO2", n: 2 }, defaultC: 0.01, ratio: 1, note: "用于书写完整反应式" },
  ];

  /* ================= 元素物种表（严格按需求文档） ================= */
  const ELEMENT_SPECIES = {
    Mg: ["Mg", "Mg2+", "Mg(OH)2", "MgO", "MgCO3"],
    Al: ["Al", "Al3+", "Al(OH)3", "AlO(OH)", "Al2O3", "Al(OH)4-", "AlF63-"],
    P: ["H3PO4", "H2PO4-", "HPO42-", "PO43-", "H3PO3", "H2PO3-", "HPO32-", "H3PO2", "H2PO2-", "P4", "PH3"],
    S: ["S", "S8", "H2S", "H2Sg", "HS-", "S2-", "SO2", "SO2g", "HSO3-", "SO32-", "HSO4-", "SO42-",
        "S2O32-", "S4O62-", "S2O42-", "S2O52-", "S2O72-", "S2O82-", "BaSO4"],
    Cl: ["Cl-", "Cl2", "Cl2g", "HClO", "ClO-", "ClO2g", "ClO3-", "ClO4-", "AgCl", "HCl", "HClg"],
    I: ["I-", "I3-", "I2", "ICl", "IO3-", "AgI", "PbI2"],
    N: ["N2", "NH4+", "NH3", "NH3g", "HNO2", "NO2-", "NO3-", "NO", "NO2"],
    C: ["CO2", "CO2g", "H2CO3", "HCO3-", "CO32-", "CO", "HCOOH", "HCOO-", "CH4", "CH3OH", "CH3OHg"],
    V: ["V", "V2+", "V(OH)2", "V3+", "V(OH)3", "V2O3", "VO2", "VO2+", "VO2plus", "V2O5", "VO3-", "VO43-"],
    Cr: ["Cr", "Cr3+", "Cr2O3", "Cr(OH)3", "Cr(OH)4-", "CrO42-", "Cr2O72-"],
    Mn: ["Mn", "Mn2+", "MnCO3", "Mn(OH)2", "MnO", "MnO(OH)", "Mn2O3", "MnO2", "MnO42-", "CaMnO4", "MnO4-"],
    Fe: ["Fe", "Fe2+", "Fe(OH)2", "FeO", "FeS", "FeS2", "FeCO3", "Fe(CN)64-", "Fe3O4", "Fe3+", "Fe(OH)3",
         "FeO(OH)", "Fe2O3", "FeSCN2+", "Fe(CN)63-", "FeO42-"],
    Cu: ["Cu", "Cu+", "Cu2O", "Cu2S", "Cu2S2", "CuFeS2", "CuCl", "CuI", "CuCl2-", "Cu(NH3)2+",
         "Cu2+", "CuO", "Cu(OH)2", "Cu(OH)42-", "Cu(NH3)42+", "CuCl42-", "Cu2(OH)2CO3"],
    Ag: ["Ag", "Ag+", "Ag2O", "AgCl", "AgBr", "AgI", "Ag2S", "AgCl2-", "Ag(NH3)2+", "Ag(CN)2-", "Ag(S2O3)23-"],
    Zn: ["Zn", "Zn2+", "ZnO", "Zn(OH)2", "Zn(OH)42-", "ZnCO3", "ZnS",
         "Zn(NH3)42+", "Zn(CN)42-", "ZnCl42-"],
    Pb: ["Pb", "Pb2+", "Pb(OH)2", "PbO", "PbCl2", "PbSO4", "PbCO3", "PbS", "PbI2", "PbO2"],
  };

  const ELEMENT_NAMES = {
    Mg: "镁", Al: "铝", P: "磷", S: "硫", Cl: "氯", I: "碘", N: "氮",
    C: "碳", V: "钒", Cr: "铬", Mn: "锰", Fe: "铁", Cu: "铜", Ag: "银", Pb: "铅",
    Ca: "钙", Ba: "钡", Zn: "锌", H: "氢", F: "氟", Br: "溴",
  };

  const MAIN_ELEMENTS = ["Mg", "Al", "P", "S", "Cl", "I", "N", "C", "V", "Cr", "Mn", "Fe", "Cu", "Ag", "Zn", "Pb"];

  /* ================= 默认参考态 ================= */
  // 当物种含有非主元素、又没有被启用的配体覆盖时，用它作参考态。
  const DEFAULT_REF = {
    Ca: "Ca2+", Ba: "Ba2+", Ag: "Ag+", Zn: "Zn2+", Fe: "Fe2+", Cu: "Cu2+", Pb: "Pb2+",
    F: "F-", Br: "Br-", I: "I-", Cl: "Cl-",
    C: "CO32-", S: "SO42-", N: "NH3", P: "PO43-",
  };

  /* ================= 预设 ================= */
  const PRESETS = [
    {
      id: "fe-rust", name: "铁 · 锈蚀", short: "Fe",
      main: "Fe", Ctot: 1e-3, mode: "A",
      selected: ["Fe", "Fe2+", "Fe3+", "Fe(OH)2", "Fe(OH)3", "Fe3O4"],
      env: [], xmin: 0, xmax: 14, emin: -1.2, emax: 1.6,
      note: "锈蚀初期以 Fe(OH)₃ 为主；FeO(OH)/Fe₂O₃ 可在物种列表中勾选",
    },
    {
      id: "al-amphoteric", name: "铝 · 酸碱两性", short: "Al",
      main: "Al", Ctot: 1e-3, mode: "A",
      selected: ["Al", "Al3+", "Al(OH)3", "Al(OH)4-"],
      env: [], xmin: 0, xmax: 14, emin: -2.4, emax: 1.0,
      note: "酸性区 Al³⁺，近中性 Al(OH)₃，强碱区 Al(OH)₄⁻",
    },
    {
      id: "mn-competition", name: "锰 · 竞争反应", short: "Mn",
      main: "Mn", Ctot: 1e-3, mode: "A",
      selected: ["Mn", "Mn2+", "Mn(OH)2", "MnO(OH)", "MnO2", "MnO42-", "MnO4-"],
      env: [], xmin: 0, xmax: 14, emin: -1.5, emax: 2.2,
      note: "Mn²⁺ → Mn(OH)₂/MnO₂ 的 pH 竞争；高电位出现 MnO₄⁻/MnO₄²⁻",
    },
    {
      id: "cl-room", name: "氯 · 常温歧化", short: "Cl",
      main: "Cl", Ctot: 1e-2, mode: "A",
      selected: ["Cl-", "Cl2", "Cl2g", "HClO", "ClO-"],
      env: [], xmin: 0, xmax: 14, emin: -0.2, emax: 1.8,
      note: "只列 Cl(−1) 与 Cl(+1)；ClO₃⁻/ClO₄⁻ 可在物种列表中主动勾选",
    },
    {
      id: "n-nitrite", name: "氮 · 亚硝酸–氨", short: "N",
      main: "N", Ctot: 1e-2, mode: "A",
      selected: ["N2", "NH4+", "NH3", "NH3g", "HNO2", "NO2-", "NO"],
      env: [], xmin: 0, xmax: 14, emin: -0.8, emax: 1.6,
      note: "亚硝酸–氨体系制 N₂；不列 NO₃⁻（可在物种列表中主动勾选）",
    },
    {
      id: "s-rich", name: "硫 · 亚硫酸–硫代硫酸", short: "S",
      main: "S", Ctot: 1e-2, mode: "A",
      selected: ["S8", "H2S", "H2Sg", "HS-", "S2-", "SO2", "SO2g", "HSO3-", "SO32-",
                "S2O32-", "S4O62-", "S2O42-"],
      env: [], xmin: 0, xmax: 14, emin: -0.8, emax: 1.2,
      note: "硫–亚硫酸–硫代硫酸体系；不列硫酸根（可在物种列表中主动勾选）",
    },
  ];

  /* ================= 归类 ================= */
  const CATEGORIES = ["单质", "氧化物", "酸碱", "离子", "难溶物", "络合物", "易挥发物"];

  // 工具：取某元素在主元素列表中的物种
  function speciesFor(el) {
    const ids = ELEMENT_SPECIES[el] || [];
    return ids.map(id => byId[id]).filter(Boolean);
  }

  return {
    DB, byId, MAIN_ELEMENTS, ELEMENT_SPECIES, ELEMENT_NAMES, CATEGORIES,
    ENV, FAMILIES, LIG_FAMILY, DEFAULT_REF, PRESETS, RTln10, DG_H2O,
    speciesFor,
  };
})();
if (typeof module !== "undefined") module.exports = PBData;