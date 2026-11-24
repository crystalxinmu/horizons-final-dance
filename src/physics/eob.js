/**
 * 双黑洞轨道动力学（几何单位 G = c = 1，总质量 M 归一为 1，时间以 GM/c³ 计）。
 *
 * 旋近与坠落：有效单体（EOB）哈密顿量（Buonanno & Damour 2000；3PN 势取自 Damour, Jaranowski & Schäfer 2000）
 *   H  = √(1 + 2ν(H_eff − 1))
 *   H_eff² = A(r)·(1 + p_φ²/r² + p_r²/B(r) + Q(r))
 *   A(u) = 1 − 2u + 2νu³ + a₄νu⁴，B(u) = (1 − 2u)⁻¹，Q(u) = 2ν(4 − 3ν)u³
 *   p_φ ≡ P_φ/(μM)，p_r ≡ P_r/(μc)；运动方程中 ν 被 μ 的归一化消去。
 * 辐射反作用：圆轨道 2.5PN 能量通量
 *   F = (32/5) ν² x⁵ [1 + f₁x + f_{3/2}x^{3/2} + f₂x² + f_{5/2}x^{5/2}]，x = (MHΩ)^{2/3}
 *   ṗ_φ = −F/(νΩ)：只施加切向反作用，于是 Ḣ = −F 严格成立。
 * 铃宕：Schwarzschild l = m = 2, n = 0 准正规模式，Mω = 0.373672 − 0.088963 i。
 */

export const G = 6.67430e-11;                   // m³ kg⁻¹ s⁻²
export const C = 299792458;                     // m s⁻¹
export const MSUN = 1.98892e30;                 // kg
export const T_UNIT = (G * MSUN) / (C ** 3);    // 1 M☉ 对应的几何时间 ≈ 4.9255e-6 s
export const L_UNIT = C ** 5 / G;               // c⁵/G ≈ 3.628e52 W

const A4 = 94 / 3 - (41 * Math.PI * Math.PI) / 32;          // 4PN 势系数 ≈ 18.688（微扰论值）
const PEAK_FLUX = 0.0010;      // 等质量并合的峰值引力波光度 / (c⁵/G)：数值相对论结果，约 3.6×10⁴⁹ W
export const OMEGA_QNM = 0.3736717085;          // Re(Mω)
export const DAMP_QNM = 0.0889621576;           // −Im(Mω)
export const TAU_QNM = 1 / DAMP_QNM;            // ≈ 11.24 M
export const ISCO = 6;                          // 最内稳定圆轨道 /M

const smooth = (x) => { x = Math.min(Math.max(x, 0), 1); return x * x * (3 - 2 * x); };

/** 并合残骸：辐射能量分数与最终自旋（对非自旋、准圆轨道数值相对论结果的拟合） */
export function remnant(nu) {
  const eRad = 0.0483 * Math.pow(4 * nu, 1.5);                         // q = 1 时 4.8%
  const spin = Math.sqrt(12) * nu - 3.871 * nu * nu + 4.028 * nu ** 3; // Rezzolla et al. 2008
  return { eRad, mFin: 1 - eRad, spin };
}

/** 圆轨道 2.5PN 能量通量的 Taylor 级数（M = 1）；x ≳ 0.2 时发散，仅供对照 */
export function pnFluxTaylor(x, nu) {
  const f = fluxCoeffs(nu), s = Math.sqrt(x);
  return (32 / 5) * nu * nu * x ** 5 * (1 + x * (f.f1 + s * (f.f15 + s * (f.f2 + s * f.f25))));
}

/**
 * 能量通量（M = 1），Padé[2/3] 重求和（Damour, Iyer & Sathyaprakash 2000 的 P-近似思想）：
 *   f(t) = (1 + n₁t + n₂t²)/(1 + d₁t + d₂t² + d₃t³)，t = √x，与 2.5PN 级数展开到 t⁵ 一致。
 * Taylor 级数在 x ≳ 0.2 发散；Padé 把高阶项重求和为一个极点 x_pole，
 * ν→0 时 x_pole ≈ 0.216，落在 test-mass 光子环 1/6 之外——通量在轨道坠落前急剧上升。
 * 自极点起通量封顶于数值相对论的峰值光度 PEAK_FLUX·(4ν)²。对等质量系统，全程积分 ∫F dτ ≈ 0.035 M，
 * 约为数值相对论辐射总能量（0.048 M）的七成：本模型的高光度阶段偏短，这是已知的系统偏差。
 */
export function pnFlux(x, nu) {
  const m = fluxModel(nu), t = Math.sqrt(x);
  if (t >= m.pole) return m.peak;                   // 自极点起：真实通量有限，取 NR 峰值引力波光度
  const { f, d1, d2, d3 } = m;
  const pade = (1 + d1 * t + (d2 + f.f1) * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t);
  return Math.min((32 / 5) * nu * nu * x ** 5 * pade, m.peak);
}

function fluxCoeffs(nu) {                         // 2.5PN 通量因子 f(x) = 1 + f₁x + f_{3/2}x^{3/2} + f₂x² + f_{5/2}x^{5/2}
  return {
    f1: -1247 / 336 - (35 / 12) * nu,
    f15: 4 * Math.PI,
    f2: -44711 / 9072 + (9271 / 504) * nu + (65 / 576) * nu * nu,
    f25: -Math.PI * (8191 / 672 + (535 / 24) * nu),
  };
}

/** 峰值引力波光度（单位 c⁵/G），随对称质量比按 (4ν)² 标度 */
export const peakFlux = (nu) => PEAK_FLUX * Math.pow(4 * nu, 2);

const fluxModels = new Map();
function fluxModel(nu) {                            // Padé 系数与极点位置，按 ν 缓存
  const hit = fluxModels.get(nu);
  if (hit) return hit;
  const f = fluxCoeffs(nu);
  const det = (f.f2 - f.f1 * f.f1) * f.f1 - f.f15 * f.f15;
  const d1 = (-f.f25 * f.f1 + f.f1 * f.f1 * f.f15 + f.f15 * f.f2) / det;
  const d2 = (f.f1 * f.f1 * f.f2 - f.f2 * f.f2 + f.f15 * f.f25 - f.f1 * f.f15 * f.f15) / det;
  const d3 = -f.f15 - f.f1 * d1;
  let pole = Infinity;                              // 分母的最小正实根，t = √x
  for (let t = 0.0005; t < 4; t += 0.0005) {
    if (1 + d1 * t + d2 * t * t + d3 * t * t * t <= 0) { pole = t; break; }
  }
  const model = { f, d1, d2, d3, pole, peak: peakFlux(nu) };
  fluxModels.set(nu, model);
  return model;
}

/** EOB 势函数：A（度规势）、A' = dA/dr、B、Q；u = M/r */
function potential(u, nu) {
  const u2 = u * u, u3 = u2 * u;
  return {
    A: 1 - 2 * u + 2 * nu * u3 + A4 * nu * u3 * u,      // 3PN 结果；ν→0 时精确还原 1 − 2u
    Ap: u2 * (2 - 6 * nu * u2 - 4 * A4 * nu * u3),      // dA/dr = −u²·dA/du
    B: 1 / (1 - 2 * u),
    Q: 2 * nu * (4 - 3 * nu) * u3,
  };
}

/** EOB 哈密顿量及其偏导：H（以 Mc² 计）、ṙ、Ω = dφ/dt、∂H/∂r */
export function hamiltonian(r, pr, pp, nu) {
  const { A, Ap, B, Q } = potential(1 / r, nu);
  const r2 = r * r;
  const S = 1 + (pp * pp) / r2 + (pr * pr) / B + Q;
  const Heff = Math.sqrt(A * S);
  const H = Math.sqrt(Math.max(1 + 2 * nu * (Heff - 1), 1e-12));
  const k = 1 / (2 * H * Heff);
  const dHdr = k * (Ap * S + A * ((2 * (pr * pr) - (2 * (pp * pp)) / r) / r2 - (3 * Q) / r));
  return { H, rdot: 2 * k * A * pr * B, om: (2 * k * A * pp) / r2, dHdr };
}

/**
 * 分离 r 处圆轨道的 { p_φ, H, Ω }：p_r = 0 且 ∂H/∂r = 0，即 ∂[A(1 + p_φ²/r² + Q)]/∂r = 0。
 * 取 A = 1 − 2u、Q = 0 时还原 test-mass 的精确角动量 L² = r²/(r − 3)。
 * r 已越过势垒（不存在圆轨道）时返回 NaN。
 */
export function circularOrbit(r, nu) {
  const { A, Ap, Q } = potential(1 / r, nu);
  const Gp = Ap * (1 + Q) - (3 * A * Q) / r;          // G = A(1+Q) 对 r 的导数
  const den = 2 * A - Ap * r;
  if (den <= 0) return { pp: NaN, H: NaN, om: NaN };
  const pp = Math.sqrt((Gp * r ** 3) / den);
  const { H, om } = hamiltonian(r, 0, pp, nu);
  return { pp, H, om };
}

/**
 * 有效度规的两个特征半径：
 *  - 势垒顶端（最内稳定圆轨道 ISCO）：圆轨道角动量的极小点，ν ≥ ν_bif ≈ 0.2285 时不存在；
 *  - 光子环：圆轨道角速度的极大点，ν→0 时给出 r = 3M。
 */
export function barrierRadius(nu) {
  const step = 0.002, r0 = 2.2;
  let a = circularOrbit(r0, nu).pp, b = circularOrbit(r0 + step, nu).pp;
  let isco = NaN, lr = NaN, omMax = -1;
  for (let k = 2; r0 + k * step < 9; k++) {
    const r = r0 + k * step, c = circularOrbit(r, nu);
    if (c.om > omMax) { omMax = c.om; lr = r; }
    if (isFinite(a) && isFinite(c.pp) && b < a && b <= c.pp) isco = r - step;   // 局部极小
    a = b; b = c.pp;
  }
  return { isco, lr, stop: Number.isFinite(isco) ? isco : lr };
}

/**
 * 预积分整条轨道，并按“屏幕时间”均匀重采样。
 * 屏幕时钟是非线性的：旋近阶段令屏幕角速度 V(Ω) = vMin·(Ω/Ω₀)^β，并合之后恒为 wPost（单位 M/s）。
 * 这样毫秒量级的并合与铃宕才看得见；物理量本身完全不受回放时钟影响。
 */
export function buildTrajectory({ q = 1, r0 = 10, vMin = 0.9, beta = 0.5, wPost = 6, tail = 66, onStep = null } = {}) {
  const m1 = q / (1 + q), m2 = 1 / (1 + q), nu = m1 * m2;
  const R_TRIG = 2.4;   // 两视界即将接触：开始并合过渡
  const R_RING = 0.9;   // 共同视界已经形成：进入铃宕
  const R_FALL = 2.9;   // 自此不再信任 PN 势，见 fB 的终盘分支
  const RATE = 240;     // 屏幕时间采样率 /s
  const bar = barrierRadius(nu);
  // 屏幕角速度 V(Ω) = V₀·(Ω/Ω₀)^β。圈数 ∝ 1/ν，所以 V₀ 随 ν 减小而提高，使各质量比的观看时长相近
  const V0 = vMin * Math.sqrt(0.25 / nu);
  const W_MERGE = 3;    // 并合瞬间的慢镜头回放速率 dτ/ds（M/s）
  let om0b = 1;         // 稍后设为 Ω₀^β，使初始屏幕角速度恰为 V₀
  const screenRate = (om) => (Math.pow(om, 1 - beta) * om0b) / V0;     // 旋近阶段的 ds/dτ = Ω / V(Ω)

  const orbit = (r) => circularOrbit(r, nu);
  const dEdr = (r) => { const d = 1e-5 * r; return (orbit(r + d).H - orbit(r - d).H) / (2 * d); };
  const adiabatic = (r) => {                       // 准圆轨道 + 能量平衡
    const { H, om } = orbit(r);
    const flux = pnFlux(Math.pow(om * H, 2 / 3), nu);
    return { rdot: -flux / dEdr(r), om, flux };
  };
  const fA = (t, y, d) => {                        // 阶段 A 的导数，y = [r, φ, s]
    const a = adiabatic(Math.max(y[0], bar.stop));
    d[0] = a.rdot; d[1] = a.om; d[2] = screenRate(a.om);
  };

  // ——— 阶段 A：绝热旋近 ———
  const rows = [];                                  // τ, r, φ, Ω, F, 屏幕时间 s
  let tau = 0, y = [r0, 0, 0], a = adiabatic(r0);
  om0b = Math.pow(a.om, beta);
  for (;;) {
    rows.push(tau, y[0], y[1], a.om, a.flux, y[2]);
    // 更早地交给哈密顿量 + 辐射反作用；没有势垒（ν ≥ ν_bif）时以光子环为界
    if (Math.abs(a.rdot) >= 0.05 || y[0] <= bar.stop * 1.7) break;
    if (!(tau < 1e5)) throw new Error('旋近积分未能在预期时间内到达势垒');
    const h = Math.min(0.02 / a.om, (0.01 * y[0]) / Math.max(Math.abs(a.rdot), 1e-12), 1);
    y = rk4(fA, tau, y, h);
    tau += h;
    a = adiabatic(y[0]);
  }

  // ——— 阶段 B：纯哈密顿量坠落，径向动量与绝热解的径向速度衔接（ṙ = A·p_r·B/(H·H_eff)）———
  const rs = y[0], phase0 = y[1], s0 = y[2], tauStop = tau;   // 阶段 A 结束时的状态
  const circ = orbit(rs), pot = potential(1 / rs, nu);
  const Heff = Math.sqrt(pot.A * (1 + (circ.pp * circ.pp) / (rs * rs) + pot.Q));
  const pr0 = (a.rdot * circ.H * Heff) / (pot.A * pot.B);
  // stage 0：辐射反作用驱动（到势垒顶端/光子环）；1：纯哈密顿量坠落；2：并合过渡；3：铃宕
  const ctx = { stage: 0, tauTrig: NaN, tauRing: NaN, wTrig: wPost, xTrig: 0.3, phiRing: 0, omFall: NaN };
  const cur = { om: 0, flux: 0 };
  const peak = peakFlux(nu);
  // 并合瞬间（两个视界融为一体，物理上只有约 1M）放慢到 W_MERGE，随后在 4M 内恢复到 wPost
  const warp = (t) => ctx.stage < 3
    ? ctx.wTrig + (W_MERGE - ctx.wTrig) * smooth((t - ctx.tauTrig) / 0.4)
    : W_MERGE + (wPost - W_MERGE) * smooth((t - ctx.tauRing) / 4);
  const fB = (t, yy, d) => {                        // y = [r, φ, p_r, p_φ, s]
    if (ctx.stage === 3) {                          // 铃宕：轨道已终结，只推进回放时钟
      d[0] = d[1] = d[2] = d[3] = 0; d[4] = 1 / warp(t); cur.om = 0;
      cur.flux = peak * Math.exp(-2 * (t - ctx.tauRing) * DAMP_QNM);
      return;
    }
    if (ctx.stage >= 1 && yy[0] < R_FALL) {         // 终盘：PN 展开的势在 r ≲ 3M 已不可信，
      const omFall = ctx.omFall > 0 ? ctx.omFall : hamiltonian(yy[0], yy[2], yy[3], nu).om;
      d[0] = -0.35; d[1] = omFall; d[2] = d[3] = 0;       // 这里以恒定径向速率并入，角速度锁定为进入时的值
      d[4] = ctx.stage <= 1 ? screenRate(omFall) : 1 / warp(t);
      cur.om = omFall; return;
    }
    const { H, rdot, om, dHdr } = hamiltonian(yy[0], yy[2], yy[3], nu);
    d[0] = rdot; d[1] = om; d[2] = -dHdr;
    // 辐射反作用驱动阶段取 Padé 值（已封顶），之后的坠落与并合过渡处于光度峰值
    cur.flux = ctx.stage === 0 ? pnFlux(Math.pow(om * H, 2 / 3), nu) : peak;
    d[3] = ctx.stage === 0 ? -cur.flux / (nu * om) : 0;
    d[4] = ctx.stage <= 1 ? screenRate(om) : 1 / warp(t);
    cur.om = om;
  };
  const d = [0, 0, 0, 0, 0];
  y = [rs, phase0, pr0, circ.pp, s0];
  fB(tau, y, d);
  for (;;) {
    rows.push(tau, y[0], y[1], cur.om, cur.flux, y[4]);
    if (onStep) onStep({ stage: ctx.stage, tau, r: y[0], pp: y[3], pr: y[2], rdot: d[0], om: cur.om });
    if (ctx.stage === 3 && tau >= ctx.tauRing + tail) break;
    if (!(tau < 1e5)) throw new Error('EOB 坠落积分发散');
    let h = ctx.stage === 3 ? 2
      : Math.min(0.02 / Math.max(cur.om, 1e-6), (0.01 * Math.max(y[0], 0.3)) / Math.max(Math.abs(d[0]), 1e-9));
    let next = rk4(fB, tau, y, h);
    for (let tries = 0; !next.every(Number.isFinite) && tries < 40; tries++) { h *= 0.5; next = rk4(fB, tau, y, h); }
    if (!next.every(Number.isFinite)) throw new Error(`EOB 坠落积分发散：stage=${ctx.stage} τ=${tau.toFixed(3)} r=${y[0].toFixed(4)} p_r=${y[2].toExponential(3)} p_φ=${y[3].toFixed(4)}`);
    y = next; tau += h;
    fB(tau, y, d);
    if (ctx.stage === 0 && y[0] <= bar.stop) { ctx.stage = 1; ctx.omFall = cur.om; fB(tau, y, d); }
    else if (ctx.stage === 1 && y[0] < R_FALL && !(ctx.omFall > 0)) { ctx.omFall = cur.om; fB(tau, y, d); }
    else if (ctx.stage === 1 && y[0] < R_TRIG) {
      ctx.stage = 2; ctx.tauTrig = tau; ctx.xTrig = Math.pow(cur.om, 2 / 3);
      ctx.wTrig = 1 / screenRate(cur.om);           // 此刻的 dτ/ds，保证屏幕时钟连续
      fB(tau, y, d);
    } else if (ctx.stage === 2 && y[0] < R_RING) {
      ctx.stage = 3; ctx.tauRing = tau; ctx.phiRing = y[1];
      y[0] = 0; y[2] = 0;                            // 共同视界之后轨道已无意义
      fB(tau, y, d);
    }
  }

  // ——— 按屏幕时间均匀重采样 ———
  const n = rows.length / 6;
  const sEnd = rows[n * 6 - 1];
  const count = Math.ceil(sEnd * RATE) + 2;
  const T = {
    tau: new Float32Array(count), r: new Float32Array(count), phi: new Float32Array(count),
    om: new Float32Array(count), flux: new Float32Array(count), w: new Float32Array(count), h: new Float32Array(count),
  };
  const screenAt = (t) => {                          // 物理 τ → 屏幕时间
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (rows[mid * 6] < t) lo = mid; else hi = mid; }
    return rows[lo * 6 + 5];
  };
  const sRing = screenAt(ctx.tauRing), sTrig = screenAt(ctx.tauTrig);
  const hRing = 4 * nu * ctx.xTrig;                  // 铃宕起始振幅 = 并合瞬间的四极振幅
  let j = 0;
  const cols = [T.tau, T.r, T.phi, T.om, T.flux];
  for (let k = 0; k < count; k++) {
    const s = Math.min(k / RATE, sEnd);
    while (j < n - 2 && rows[(j + 1) * 6 + 5] < s) j++;
    const i0 = j * 6, i1 = i0 + 6;
    const span = rows[i1 + 5] - rows[i0 + 5];
    const f = span > 0 ? Math.min(Math.max((s - rows[i0 + 5]) / span, 0), 1) : 0;
    for (let c = 0; c < 5; c++) cols[c][k] = rows[i0 + c] + f * (rows[i1 + c] - rows[i0 + c]);
    if (T.tau[k] >= ctx.tauRing) T.flux[k] = peak * Math.exp(-2 * (T.tau[k] - ctx.tauRing) * DAMP_QNM);
    T.w[k] = smooth((R_TRIG - T.r[k]) / (R_TRIG - R_RING));
    const dt = T.tau[k] - ctx.tauRing;
    T.h[k] = dt < 0
      ? 4 * nu * Math.pow(T.om[k], 2 / 3) * Math.cos(2 * T.phi[k])
      : hRing * Math.exp(-dt / TAU_QNM) * Math.cos(2 * ctx.phiRing + OMEGA_QNM * dt);
  }

  const { eRad, mFin, spin } = remnant(nu);
  return {
    q, m1, m2, nu, rate: RATE, count, tbl: T, eRad, mFin, spin, bar,
    tauStop, tauTrig: ctx.tauTrig, tauRing: ctx.tauRing, xTrig: ctx.xTrig, phiRing: ctx.phiRing,
    sTrig, sRing, sEnd,
    /** 屏幕时间 → 状态（线性插值） */
    sample(s) {
      const t = Math.min(Math.max(s, 0), sEnd) * RATE;
      const i = Math.min(Math.floor(t), count - 2), f = t - i;
      const L = (a) => a[i] + f * (a[i + 1] - a[i]);
      const om = L(T.om);
      return { tau: L(T.tau), r: L(T.r), phi: L(T.phi), om, flux: L(T.flux), w: L(T.w), h: L(T.h), x: Math.pow(om, 2 / 3) };
    },
  };
}

function rk4(f, tau, y, h) {
  const n = y.length, k1 = y.slice(), k2 = y.slice(), k3 = y.slice(), k4 = y.slice(), t = y.slice();
  for (let s = 0; s < 4; s++) {
    const c = s === 3 ? 1 : 0.5;
    for (let i = 0; i < n; i++) t[i] = s === 0 ? y[i] : y[i] + c * h * (s === 1 ? k1[i] : s === 2 ? k2[i] : k3[i]);
    f(tau + c * h, t, s === 0 ? k1 : s === 1 ? k2 : s === 2 ? k3 : k4);
  }
  const out = y.slice();
  for (let i = 0; i < n; i++) out[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
  return out;
}
