/**
 * 光线追踪方程的 CPU 参考实现，与着色器里的方程逐项一致，供测试校验与文档引用。
 *
 * Schwarzschild 时空中的零测地线精确满足 Binet 方程
 *     u'' + u = 3M u²      (u = 1/r，撇号表示对轨道平面内的极角求导)
 * 写成三维形式即
 *     x'' = −3 M h² x / r⁵ ,   h² = |x × x'|
 * 其中仿射参数以远处 |x'| = 1 归一。它给出的轨迹与真测地线完全相同：
 * 弱场极限回到偏折角 4GM/(c²b)，临界瞄准参数为 3√3 M（即阴影半径）。
 * 多个黑洞按 test-mass 方式叠加加速度，这是双黑洞度规的一阶近似（忽略黑洞间的背反应）。
 */

/** 状态 q = [x, y, z, vx, vy, vz] 对仿射参数的导数 */
export function derivative(q, bodies, out) {
  out[0] = q[3]; out[1] = q[4]; out[2] = q[5];
  let ax = 0, ay = 0, az = 0;
  for (const b of bodies) {
    const dx = q[0] - b.pos[0], dy = q[1] - b.pos[1], dz = q[2] - b.pos[2];
    const r2 = dx * dx + dy * dy + dz * dz;
    const cx = dy * q[5] - dz * q[4], cy = dz * q[3] - dx * q[5], cz = dx * q[4] - dy * q[3];
    const k = (-3 * b.mass * (cx * cx + cy * cy + cz * cz)) / (r2 * r2 * Math.sqrt(r2));
    ax += k * dx; ay += k * dy; az += k * dz;
  }
  out[3] = ax; out[4] = ay; out[5] = az;
}

function nearestHorizon(q, bodies) {
  let d = Infinity;
  for (const b of bodies) d = Math.min(d, Math.hypot(q[0] - b.pos[0], q[1] - b.pos[1], q[2] - b.pos[2]) - 2 * b.mass);
  return d;
}

/**
 * 从 p 沿单位方向 v 积分光线。
 * @returns {{captured:boolean, escaped:boolean, q:number[], steps:number}}
 */
export function traceRay(p, v, bodies, { maxSteps = 200000, escapeRadius = 3e4, stepScale = 1 } = {}) {
  const q = [...p, ...v];
  const k1 = [0, 0, 0, 0, 0, 0], k2 = [0, 0, 0, 0, 0, 0], k3 = [0, 0, 0, 0, 0, 0], k4 = [0, 0, 0, 0, 0, 0];
  const t = [0, 0, 0, 0, 0, 0];
  let captured = false, escaped = false, steps = 0;
  for (let i = 0; i < maxSteps; i++) {
    steps = i;
    if (nearestHorizon(q, bodies) <= 0) { captured = true; break; }
    if (Math.hypot(q[0], q[1], q[2]) > escapeRadius && q[0] * q[3] + q[1] * q[4] + q[2] * q[5] > 0) { escaped = true; break; }
    const dmin = Math.max(nearestHorizon(q, bodies), 0);
    const h = Math.min(stepScale * (0.14 * dmin + 0.01), 0.25 * Math.hypot(q[0], q[1], q[2]));
    const K = [k1, k2, k3, k4];
    for (let s = 0; s < 4; s++) {
      const f = s === 1 || s === 2 ? 0.5 : 1;
      for (let j = 0; j < 6; j++) t[j] = s === 0 ? q[j] : q[j] + f * h * K[s - 1][j];
      derivative(t, bodies, K[s]);
    }
    for (let j = 0; j < 6; j++) q[j] += (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
  }
  return { captured, escaped, q, steps };
}
