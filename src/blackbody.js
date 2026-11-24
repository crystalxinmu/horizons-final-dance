import * as THREE from 'three';

/**
 * CIE 1931 2° 标准观察者配色函数的多瓣高斯拟合
 * （Wyman, Sloan & Shirley 2013, "Simple Analytic Approximations to the CIE XYZ Color Matching Functions"）。
 * λ 以 nm 计。
 */
function cmf(lam) {
  const g = (mu, s1, s2) => { const t = (lam - mu) / (lam < mu ? s1 : s2); return Math.exp(-0.5 * t * t); };
  return [
    1.056 * g(599.8, 37.9, 31.0) + 0.362 * g(442.0, 16.0, 26.7) - 0.065 * g(501.1, 20.4, 26.2),  // x̄(λ)
    0.821 * g(568.8, 46.9, 40.5) + 0.286 * g(530.9, 16.3, 31.1),                            // ȳ(λ)
    1.217 * g(437.0, 11.8, 36.0) + 0.681 * g(459.0, 26.0, 13.8),                            // z̄(λ)
  ];
}

/**
 * 生成“色温 → 线性 sRGB 色度”的一维查找表。
 * 每个温度下对普朗克黑体谱在可见波段做数值积分：
 *     B(λ, T) ∝ 1 / (λ⁵ (e^{hc/(λkT)} − 1))，hc/k = 1.4388×10⁷ nm·K
 * 再经 CIE 配色函数得到 XYZ，最后用 sRGB (D65) 矩阵变到线性 RGB。色度以最大通道归一化。
 */
export function buildBlackbodyTexture(size = 1024, tMin = 800, tMax = 40000) {
  const data = new Uint16Array(size * 4);
  const toHalf = THREE.DataUtils.toHalfFloat;
  for (let i = 0; i < size; i++) {
    const T = tMin * Math.pow(tMax / tMin, i / (size - 1));
    let X = 0, Y = 0, Z = 0;
    for (let lam = 380; lam <= 780; lam += 2) {
      const B = 1 / (lam ** 5 * Math.expm1(1.4388e7 / (lam * T)));
      const [cx, cy, cz] = cmf(lam);
      X += B * cx; Y += B * cy; Z += B * cz;
    }
    const R = Math.max(3.2406 * X - 1.5372 * Y - 0.4986 * Z, 0);
    const G = Math.max(-0.9689 * X + 1.8758 * Y + 0.0415 * Z, 0);
    const Bc = Math.max(0.0557 * X - 0.2040 * Y + 1.0570 * Z, 0);
    const peak = Math.max(R, G, Bc, 1e-30);
    data.set([toHalf(R / peak), toHalf(G / peak), toHalf(Bc / peak), toHalf(0)], i * 4);
  }
  const tex = new THREE.DataTexture(data, size, 1, THREE.RGBAFormat, THREE.HalfFloatType);
  tex.minFilter = tex.magFilter = THREE.LinearFilter;
  tex.wrapS = tex.wrapT = THREE.ClampToEdgeWrapping;
  tex.colorSpace = THREE.NoColorSpace;
  tex.needsUpdate = true;
  return tex;
}
