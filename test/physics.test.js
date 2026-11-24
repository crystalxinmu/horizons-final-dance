import { test } from 'node:test';
import assert from 'node:assert/strict';
import { traceRay } from '../src/physics/geodesic.js';
import {
  buildTrajectory, circularOrbit, barrierRadius, remnant, pnFlux, pnFluxTaylor, peakFlux,
  OMEGA_QNM, TAU_QNM, T_UNIT,
} from '../src/physics/eob.js';

// ————————————————— 光线：Schwarzschild 零测地线 —————————————————

const BH = [{ pos: [0, 0, 0], mass: 1 }];           // M = 1
const SHADOW = 3 * Math.sqrt(3);                    // 3√3 M ≈ 5.196，临界瞄准参数
const R_EMIT = 4000;

function ray(b, opts) { return traceRay([-R_EMIT, b, 0], [1, 0, 0], BH, { escapeRadius: 1.5 * R_EMIT, ...opts }); }
const deflection = (q) => Math.acos(q[3] / Math.hypot(q[3], q[4], q[5]));

/** 精确偏折角的独立参考：δ = 2∫₀^{u₀} du/√(b⁻² − u² + 2u³) − π（换元消去转折点的平方根奇异性） */
function exactDeflection(b) {
  const f = (u) => 1 / (b * b) - u * u + 2 * u ** 3;
  let lo = 0, hi = 1.5 / b;
  for (let i = 0; i < 200; i++) { const m = 0.5 * (lo + hi); if (f(m) > 0) lo = m; else hi = m; }
  const u0 = 0.5 * (lo + hi);
  const g = (th) => {
    const s = Math.sin(th), c = Math.cos(th);
    if (Math.abs(c) < 1e-7) return 2 * Math.sqrt(u0) / Math.sqrt(2 * u0 - 6 * u0 * u0);
    return (2 * u0 * s * c) / Math.sqrt(f(u0 * s * s));
  };
  const N = 20000, h = Math.PI / 2 / N;
  let sum = g(0) + g(Math.PI / 2);
  for (let i = 1; i < N; i++) sum += (i % 2 ? 4 : 2) * g(i * h);
  return (2 * sum * h) / 3 - Math.PI;
}

test('光线方程复现 Schwarzschild 阴影半径 b_c = 3√3 M', () => {
  let lo = 4, hi = 6;                               // lo 被捕获，hi 逃逸
  for (let i = 0; i < 40; i++) { const mid = 0.5 * (lo + hi); if (ray(mid).captured) lo = mid; else hi = mid; }
  const bc = 0.5 * (lo + hi);
  assert.ok(Math.abs(bc / SHADOW - 1) < 2e-3, `b_c = ${bc}, 期望 ${SHADOW}`);
});

test('偏折角与精确椭圆积分求积一致（弱场到强场）', () => {
  for (const b of [6, 10, 20, 40, 80]) {
    const r = ray(b);
    assert.ok(!r.captured);
    const d = deflection(r.q), ref = exactDeflection(b);
    assert.ok(Math.abs(d / ref - 1) < 3e-4, `b=${b}: 数值 ${d}, 精确 ${ref}`);
  }
});

test('弱场极限回到 4GM/(c²b)', () => {
  const b = 320, d = deflection(ray(b).q);
  assert.ok(Math.abs(d / (4 / b) - 1) < 0.01, `δ=${d}, 4/b=${4 / b}`);
});

test('瞄准参数小于视界半径的光线必被捕获', () => {
  for (const b of [0.5, 1, 1.9, 2.9, 5.1]) assert.ok(ray(b).captured, `b=${b} 应被捕获`);
});

test('积分步长自适应：步长缩小 16 倍结果保持不变', () => {
  const a = deflection(ray(6.5, { stepScale: 1 }).q), b = deflection(ray(6.5, { stepScale: 0.0625 }).q);
  assert.ok(Math.abs(a / b - 1) < 5e-5, `${a} vs ${b}`);
});

// ————————————————— 轨道：EOB 哈密顿量 —————————————————

test('test-mass 极限：圆轨道角动量 L² = r²/(r−3)，圆轨道精确还原 Schwarzschild', () => {
  const nu = 1e-6;                                 // 足够小以逼近 test-mass，又不至因浮点抵消损失精度
  for (const r of [6.5, 8, 15, 40, 200]) {
    const c = circularOrbit(r, nu), u = 1 / r;
    assert.ok(Math.abs((c.pp * c.pp) / (r * r / (r - 3)) - 1) < 1e-4, `L² at r=${r}`);    // 偏差为 O(ν)
    assert.ok(Math.abs(c.om / Math.pow(r, -1.5) - 1) < 1e-4, `Ω at r=${r}`);
    assert.ok(Math.abs((c.H - 1) / nu / ((1 - 2 * u) / Math.sqrt(1 - 3 * u) - 1) - 1) < 5e-7, `E at r=${r}`);
  }
});

test('等质量圆轨道能量与 2PN 级数一致到 3PN 量级', () => {
  const nu = 0.25;
  for (const [r, tol] of [[100, 6e-4], [40, 4e-3]]) {
    const c = circularOrbit(r, nu), x = Math.pow(c.om, 2 / 3);
    const pn = -0.5 * nu * x * (1 - (9 + nu) / 12 * x + (-81 - 57 * nu + nu * nu) / 24 * x * x);
    assert.ok(Math.abs((c.H - 1) / pn - 1) < tol, `r=${r}: E=${c.H - 1}, 2PN=${pn}`);
  }
});

test('特征半径：ν→0 时 ISCO = 6M、光子环 = 3M；等质量没有势垒', () => {
  const tm = barrierRadius(1e-9);
  assert.ok(Math.abs(tm.isco - 6) < 0.01 && Math.abs(tm.lr - 3) < 0.01, JSON.stringify(tm));
  const eq = barrierRadius(0.25);
  assert.ok(Number.isNaN(eq.isco) && eq.lr > 2.5 && eq.lr < 3.2, JSON.stringify(eq));
  const q4 = barrierRadius(0.16);
  assert.ok(q4.isco > 4.5 && q4.isco < 5.5, JSON.stringify(q4));
});

test('能量通量：牛顿极限、与 Taylor 级数的一致性、峰值以 NR 光度封顶', () => {
  const nu = 0.25;
  assert.ok(Math.abs(pnFlux(1e-10, nu) / ((32 / 5) * nu * nu * 1e-50) - 1) < 1e-6);
  assert.ok(Math.abs(pnFlux(0.01, nu) / pnFluxTaylor(0.01, nu) - 1) < 2e-3);   // 与 2.5PN 级数展开到 t⁵ 一致
  assert.ok(Math.abs(pnFlux(0.02, nu) / pnFluxTaylor(0.02, nu) - 1) < 2e-2);
  const peak = peakFlux(nu);
  assert.ok(Math.abs(peak / 0.0010 - 1) < 1e-9);            // 等质量 ≈ 3.6×10⁴⁹ W
  for (const x of [0.1, 0.15, 0.2, 0.25, 0.3]) assert.ok(pnFlux(x, nu) > 0 && pnFlux(x, nu) <= peak && Number.isFinite(pnFlux(x, nu)));
  assert.ok(pnFlux(0.3, nu) === peak);
});

test('能量自洽：引力波光度的全程积分接近残骸拟合给出的辐射总能量（模型的已知低估约 30%）', () => {
  for (const q of [1, 2, 4, 10]) {
    const t = buildTrajectory({ q }), T = t.tbl;
    let E = 0, fPeak = 0;
    for (let k = 1; k < t.count; k++) {
      E += 0.5 * (T.flux[k] + T.flux[k - 1]) * (T.tau[k] - T.tau[k - 1]);
      fPeak = Math.max(fPeak, T.flux[k]);
    }
    assert.ok(E / t.eRad > 0.6 && E / t.eRad < 1.0, `q=${q}: ∫F dτ = ${E}, ε = ${t.eRad}`);
    assert.ok(fPeak > 0 && fPeak <= peakFlux(t.nu) * 1.0001, `q=${q}: 峰值 ${fPeak} 超出封顶`);
  }
});

test('并合残骸：辐射能量与最终自旋符合数值相对论', () => {
  const eq = remnant(0.25);
  assert.ok(Math.abs(eq.eRad - 0.0483) < 1e-9 && Math.abs(eq.spin - 0.687) < 0.01, `χ_f=${eq.spin}`);
  const q4 = remnant(0.16);
  assert.ok(q4.eRad > 0.015 && q4.eRad < 0.035 && q4.spin > 0.4 && q4.spin < 0.6);      // NR（q=4）：ΔE ≈ 2.5%，χ_f ≈ 0.48
});

// ————————————————— 整条轨道 —————————————————

test('等质量旋近时间与四极辐射 leading-order 一致', () => {
  const t = buildTrajectory({ q: 1 }), nu = 0.25;
  const lead = (5 / 256) * (10 ** 4 - 6 ** 4) / nu;     // t = 5a⁴/(256ν)，从 a = 10M 到 6M，以 M 计
  // 到达 r = 6M 的物理时间：在表里找第一次 r < 6
  let tau6 = NaN;
  for (let k = 0; k < t.count; k++) if (t.tbl.r[k] <= 6) { tau6 = t.tbl.tau[k]; break; }
  const ratio = tau6 / lead;
  assert.ok(ratio > 0.75 && ratio < 1.1, `τ(r=6M)=${tau6}, leading-order=${lead}, 比值=${ratio}`);
});

for (const q of [1, 1.5, 2, 3, 4, 10]) {
  test(`轨迹表 q=${q}：单调、有限，并合后残骸权重为 1，铃宕按 QNM 衰减`, () => {
    const t = buildTrajectory({ q }), { tbl, count } = t;
    for (const col of ['tau', 'r', 'phi', 'om', 'w', 'h']) {
      for (let k = 0; k < count; k++) assert.ok(Number.isFinite(tbl[col][k]), `${col}[${k}]`);
    }
    for (let k = 1; k < count; k++) assert.ok(tbl.tau[k] >= tbl.tau[k - 1] - 1e-6);
    for (let k = 0; k < count; k++) assert.ok(tbl.w[k] >= 0 && tbl.w[k] <= 1);
    const iRing = Math.floor(t.sRing * t.rate);
    assert.ok(t.sTrig > 15 && t.sTrig < 45 && t.sRing > t.sTrig && t.sEnd > t.sRing + 5, `s: ${t.sTrig}, ${t.sRing}, ${t.sEnd}`);
    assert.equal(tbl.w[Math.min(iRing + 1, count - 1)], 1);
    assert.equal(tbl.r[count - 1], 0);
    // 取半个 QNM 周期内的最大振幅，避免恰逢过零点
    const amp = (ta, tb) => { let m = 0; for (let k = 0; k < count; k++) if (tbl.tau[k] >= ta && tbl.tau[k] <= tb) m = Math.max(m, Math.abs(tbl.h[k])); return m; };
    const a0 = amp(t.tauRing, t.tauRing + 9), aEnd = amp(tbl.tau[count - 1] - 9, tbl.tau[count - 1]);
    assert.ok(a0 > 0.05 && aEnd < 0.03 * a0, `铃宕振幅 ${a0} → ${aEnd}`);
    // 过零点半周期对应的物理时间应为 π/(2·Mω_QNM)：4νx·cos(2φ + ωτ)
    const zeros = [];
    for (let k = iRing + 2; k < count - 1 && zeros.length < 4; k++) if (tbl.h[k] * tbl.h[k - 1] < 0) zeros.push(tbl.tau[k]);
    const half = (zeros[3] - zeros[1]) / 2;
    assert.ok(Math.abs(half / (Math.PI / OMEGA_QNM) - 1) < 0.02, `QNM 半周期 ${half}，期望 ${Math.PI / OMEGA_QNM}`);
    assert.ok(Math.abs(TAU_QNM - 1 / 0.0889621576) < 1e-6);
  });
}

test('真实单位：60 M☉ 系统的并合时间在百毫秒量级', () => {
  assert.ok(Math.abs(T_UNIT - 4.9255e-6) < 2e-9, `t_unit=${T_UNIT}`);
  const t = buildTrajectory({ q: 1 });
  const seconds = t.tauTrig * T_UNIT * 60;
  assert.ok(seconds > 0.05 && seconds < 0.5, `t_merge = ${seconds} s`);
});

test('轨迹对象导出渲染所需的量：并合时刻、铃宕相位、特征半径都是有限数', () => {
  for (const q of [1, 4, 10]) {
    const t = buildTrajectory({ q });
    for (const k of ['tauTrig', 'tauRing', 'xTrig', 'phiRing', 'sTrig', 'sRing', 'sEnd', 'eRad', 'mFin', 'spin']) assert.ok(Number.isFinite(t[k]), `q=${q}: ${k} 不是有限数`);
    assert.ok(Number.isFinite(t.bar.stop) && t.bar.stop > 1);
  }
});
