/**
 * 《视界终舞》的屏幕空间着色器：对每个像素在 Schwarzschild 时空中积分零测地线（广义相对论光线追踪）。
 *
 * 光线方程（几何单位）：Schwarzschild 零测地线精确满足 Binet 方程
 *     u'' + u = 3Mu²      (u = 1/r)
 * 等价的三维形式
 *     x'' = −3 M h² (x − x_c) / r⁵ ,   h² = |(x − x_c) × x'|
 * 仿射参数以无穷远处 |x'| = 1 归一，于是 x' 正是归一化到 E = 1 的光子坐标动量，可直接用于多普勒因子。
 * 多个黑洞按 test-mass 方式叠加加速度——双黑洞度规的一阶近似（见 README“已知近似”）。
 * 同一方程的 CPU 版本在 test/ 中与精确椭圆积分求积逐点核对到 3×10⁻⁴。
 */

export const vertexShader = /* glsl */ `
varying vec2 vUv;
void main() {
  vUv = uv;
  gl_Position = vec4(position.xy, 0.0, 1.0);
}
`;

export const fragmentShader = /* glsl */ `
varying vec2 vUv;

uniform vec2  uResolution;
uniform vec3  uCamPos;
uniform mat3  uCamBasis;       // 列向量：right, up, forward
uniform float uTanHalfFov;
uniform float uTau;            // 仿真时间 τ（单位 M）
uniform float uStepScale;
uniform int   uMaxSteps;
uniform float uEscRadius;
uniform int   uSupersample;       // 每个像素的光线数（N×N）

uniform vec4  uBh[3];          // xyz 位置, w 质量；[0] [1] 两个黑洞，[2] 并合残骸
uniform vec2  uQnm;            // 残骸 l = m = 2 准正规模式的四极形变：振幅、相位

uniform sampler2D uBlackbody;  // 色温 → 线性 RGB 色度（普朗克谱 × CIE 配色函数）

uniform vec3  uDiskC[4];       // 盘中心
uniform vec4  uDiskP[4];       // r_in, r_out, 中心质量, 权重
uniform vec3  uDiskS[4];       // 温度尺度(K)、纹理种子、类型（0 微盘及并合后内盘 / 1 环双星盘）
uniform vec3  uDiskV[4];       // 盘中心（黑洞）的坐标速度
uniform vec2  uSpiral;         // 环双星盘 m = 2 螺旋密度波：图案相位、缠绕数
uniform vec3  uPulse;          // 并合后向外传播的声波脉冲：半径、幅度、宽度

const float TAU = 6.28318530718;

bool  gInside = false;         // 本步的采样点是否有落入视界的
float gDmin = 1e9;             // 本步采样点到最近视界距离的最小值
float gRmin = 1e9;             // 本步采样点到最近黑洞中心距离的最小值

// ————————————————— 噪声 —————————————————

uint hashU(uint x) {                                   // 低差异的整数散列（Thomas Wang / PCG 风格）
  x ^= x >> 16; x *= 0x7feb352du;
  x ^= x >> 15; x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}
float hash13(ivec3 c) {                                // 对整数格点均匀散列到 [0, 1)
  uint h = hashU(uint(c.x)) ^ hashU(uint(c.y) + 0x68bc21ebu) ^ hashU(uint(c.z) + 0x02e5be93u);
  return float(hashU(h)) * (1.0 / 4294967296.0);
}
vec3 hash33(ivec3 c) {
  return vec3(hash13(c), hash13(c + ivec3(1, 7, 13)), hash13(c + ivec3(29, 5, 11)));
}
float vnoise(vec3 p) {
  ivec3 i = ivec3(floor(p));
  vec3 f = fract(p);
  f = f * f * (3.0 - 2.0 * f);
  return mix(
    mix(mix(hash13(i),                   hash13(i + ivec3(1, 0, 0)), f.x),
        mix(hash13(i + ivec3(0, 1, 0)),   hash13(i + ivec3(1, 1, 0)), f.x), f.y),
    mix(mix(hash13(i + ivec3(0, 0, 1)),   hash13(i + ivec3(1, 0, 1)), f.x),
        mix(hash13(i + ivec3(0, 1, 1)),   hash13(i + ivec3(1, 1, 1)), f.x), f.y),
    f.z);
}
float fbm(vec3 p) {
  float s = 0.0, a = 0.5;
  for (int i = 0; i < 3; i++) { s += a * vnoise(p); p = p * 2.07 + vec3(1.7, 9.2, 3.1); a *= 0.5; }
  return s;
}

// ————————————————— 引力：零测地线 —————————————————

/** 残骸的 l = m = 2, n = 0 准正规模式：赤道面内的四极形变 sin²θ·cos 2(φ − ψ)，同时作用于视界与光子球 */
float qnmDeform(vec3 d) {
  float q = ((d.x * d.x - d.z * d.z) * cos(uQnm.y) - 2.0 * d.x * d.z * sin(uQnm.y))
          * (d.x * d.x + d.z * d.z) / (dot(d, d) * dot(d, d));
  return 1.0 + uQnm.x * q;
}

float effRadius(vec3 d, int k) {
  return k == 2 ? length(d) / max(qnmDeform(d), 0.25) : length(d);
}

/** x'' = −Σ 3 M_i h_i² (x − x_i) / r_i⁵，h_i = |(x − x_i) × x'| */
vec3 accel(vec3 p, vec3 v) {
  vec3 a = vec3(0.0);
  for (int k = 0; k < 3; k++) {
    float m = uBh[k].w;
    if (m <= 0.0) continue;
    vec3 d = p - uBh[k].xyz;
    float r = effRadius(d, k);
    gInside = gInside || (r <= 2.0 * m);
    gDmin = min(gDmin, max(r - 2.0 * m, 0.0));         // 步长由到视界的距离决定
    gRmin = min(gRmin, r);
    vec3 c = cross(d, v);
    float r2 = r * r;
    a += (-3.0 * m * dot(c, c) / (r2 * r2 * r)) * d;
  }
  return a;
}

bool insideAny(vec3 p) {
  for (int k = 0; k < 3; k++) {
    if (uBh[k].w > 0.0 && effRadius(p - uBh[k].xyz, k) <= 2.0 * uBh[k].w) return true;
  }
  return false;
}

// ————————————————— 黑体色与背景星空 —————————————————

vec3 blackbody(float T) {
  float u = (log2(max(T, 800.0)) - 9.643856) / 5.643856;          // 800 K … 40000 K
  return texture2D(uBlackbody, vec2(clamp(u, 0.0, 1.0), 0.5)).rgb;
}

/** 程序化星表：立方体 six-face 网格，每格至多一颗星；小于一个像素时按面亮度守恒变暗而非闪烁 */
vec3 starfield(vec3 d) {
  vec3 ad = abs(d);
  vec2 uv; float face;
  if (ad.x >= ad.y && ad.x >= ad.z) { uv = d.yz / ad.x; face = d.x > 0.0 ? 1.0 : 2.0; }
  else if (ad.y >= ad.z)            { uv = d.zx / ad.y; face = d.y > 0.0 ? 3.0 : 4.0; }
  else                              { uv = d.xy / ad.z; face = d.z > 0.0 ? 5.0 : 6.0; }
  vec2 g = (uv * 0.5 + 0.5) * 18.0;
  vec2 id0 = floor(g);
  float px = max(fwidth(g.x), fwidth(g.y));
  vec3 col = vec3(0.0);
  for (int j = -1; j <= 1; j++) for (int i = -1; i <= 1; i++) {
    vec2 cid = id0 + vec2(float(i), float(j));
    vec3 h = hash33(ivec3(ivec2(cid), int(face)));
    if (h.z < 0.40) {                                              // 每格约 40% 的概率有一颗星
      vec2 dp = cid + 0.2 + 0.6 * h.xy - g;
      float b = hash13(ivec3(ivec2(cid), int(face) + 11));
      float bright = min(0.10 + 6.5 * pow(b, 16.0), 2.4);              // 亮度呈重尾分布，并给最亮的星封顶
      float src = 0.034, rad = max(src, 0.85 * px);
      col += blackbody(2600.0 + 14000.0 * h.x * h.x) * bright * exp(-dot(dp, dp) / (rad * rad)) * (src * src) / (rad * rad);
    }
  }
  return col;
}

vec3 background(vec3 d) {
  vec3 gn = normalize(vec3(0.36, 0.80, 0.48));                    // 银道面法向
  float b = abs(dot(d, gn));
  float band = exp(-b * b * 22.0);
  float neb = fbm(d * 4.0 + vec3(3.1, 1.7, 5.2));
  float lanes = fbm(d * 11.0 + vec3(7.0, 2.0, 1.0));              // 尘埃暗带
  vec3 col = vec3(0.026, 0.028, 0.042) * band * (0.30 + 2.2 * neb * neb) * (0.55 + 0.45 * lanes);
  col += vec3(0.0022, 0.0026, 0.0038);                             // 极弱的底光
  return col + starfield(d) * (0.75 + 0.9 * band);
}

// ————————————————— 吸积盘 —————————————————

/** 开普勒差速旋转的湍流纹理。窗口长度取 0.4 个轨道周期、时间上交错的两套采样交叉淡化，避免剪切把纹理越拉越细 */
float diskNoise(float r, float phi, float om, float rin, float seed) {
  float T = 0.4 * TAU / om;
  float t1 = mod(uTau, T), t2 = mod(uTau + 0.5 * T, T);
  float w1 = 1.0 - smoothstep(0.0, 1.0, t1 / T);
  float w2 = 1.0 - smoothstep(0.0, 1.0, t2 / T);
  vec3 base = vec3(0.0, 0.0, log(r / rin) * 6.5 + seed);
  float a1 = phi - om * t1 + 3.0 * seed, a2 = phi - om * t2 + 3.0 * seed;
  float n1 = fbm(base + vec3(cos(a1), sin(a1), 0.0) * 0.9);
  float n2 = fbm(base + vec3(cos(a2), sin(a2), 0.0) * 0.9);
  return mix(n1, n2, w2 / max(w1 + w2, 1e-4));
}

/** 对盘平面 y = 0 上的一个交点着色。返回是否命中 */
bool diskShade(int i, vec3 hit, vec3 v, out vec3 emis, out float alpha) {
  emis = vec3(0.0); alpha = 0.0;
  float weight = uDiskP[i].w;
  if (weight <= 0.002) return false;
  vec3 rel = hit - uDiskC[i];
  float r = length(rel.xz);
  float rin = uDiskP[i].x, rout = uDiskP[i].y, m = uDiskP[i].z;
  if (r < rin || r > rout) return false;
  bool circ = uDiskS[i].z > 0.5;

  float om = sqrt(m / (r * r * r));                    // Ω = √(M/r³)：Schwarzschild 中圆轨道角速度的严格结果
  float phi = atan(rel.z, rel.x);
  float n = diskNoise(r, phi, om, rin, uDiskS[i].y);

  // 密度包络：内缘在 ISCO 处陡起，外缘平滑消失
  float env = circ
    ? smoothstep(rin, rin * 1.55, r) * (1.0 - smoothstep(rout * 0.45, rout, r)) * (0.8 + 0.6 * n)
    : smoothstep(rin, rin * 1.05, r) * (1.0 - smoothstep(rout * 0.6, rout, r)) * sqrt(rin / r) * (0.30 + 2.2 * n * n);
  if (circ) {                                           // 双星 m = 2 势激发的螺旋密度波，在空腔边缘最强
    float spiral = cos(2.0 * (phi - uSpiral.x) + uSpiral.y * log(r / rin));
    env *= 1.0 + 0.45 * spiral * exp(-(r - rin) / (0.5 * rin));
  }
  if (uPulse.y > 0.0) env *= 1.0 + uPulse.y * exp(-pow((r - uPulse.x) / uPulse.z, 2.0));   // 并合后的声波脉冲

  // Shakura–Sunyaev / Novikov–Thorne 薄盘的温度分布：T ∝ [(M/r³)(1 − √(r_in/r))]^{1/4}
  float Tn = pow(m / (r * r * r) * (1.0 - sqrt(rin / r)), 0.25);

  // 引力红移与多普勒：g = ν_obs/ν_em = √(1 − 3M/r) / (1 − Ω·L)，L 为光子绕盘轴的角动量（E = 1 归一）
  vec3 k = -v;                                         // 光子实际的传播方向（朝观测者）
  float L = rel.x * k.z - rel.z * k.x;
  float g = sqrt(max(1.0 - 3.0 * m / r, 0.0)) / max(1.0 - om * L, 0.06);
  vec3 V = uDiskV[i];                                  // 黑洞轨道运动的狭义相对论多普勒
  float v2 = dot(V, V);
  if (v2 > 1e-9) g /= sqrt(1.0 - v2) * (1.0 - dot(V, k));
  g = min(g, 4.0);

  float T = uDiskS[i].x * Tn * g;                      // 观测色温 = 发射温度 × g
  float lum = pow(T / 6400.0, 1.7);                    // 相对论增亮：T_obs = g·T（显示用压缩指数，真实黑体总辐射 ∝ g⁴）
  emis = weight * env * lum * blackbody(T) * (circ ? 3.5 : 1.4);   // 环双星盘温度低，另乘增益以保持可见
  alpha = clamp(weight * env * (circ ? 1.1 : 1.6), 0.0, 1.0);
  return true;
}

// ————————————————— 主循环 —————————————————

vec3 cameraRay(vec2 ndc) {                             // ndc ∈ [−1, 1]²，y 向上
  return normalize(uCamBasis * vec3(ndc.x * (uResolution.x / uResolution.y) * uTanHalfFov, ndc.y * uTanHalfFov, 1.0));
}

/** 从相机沿 dir 反向追踪一条光线，返回线性 HDR 亮度 */
vec3 trace(vec3 dir) {
  vec3 p = uCamPos, v = dir;                           // |v| = 1 即以 E = 1 归一
  vec3 col = vec3(0.0);
  float T = 1.0;                                       // 透射率
  float esc2 = uEscRadius * uEscRadius;

  for (int step = 0; step < 2000; step++) {
    if (step >= uMaxSteps) break;                      // 步数耗尽：只在光子环附近发生，保持黑色
    gInside = false; gDmin = 1e9; gRmin = 1e9;
    // 四阶 Runge–Kutta；步长与到最近视界的距离成正比，与 test/ 中的 CPU 版本一致
    vec3 a1 = accel(p, v);
    float h = uStepScale * min(0.14 * gDmin + 0.01, 0.25 * gRmin);
    vec3 v2 = v + 0.5 * h * a1;  vec3 a2 = accel(p + 0.5 * h * v,  v2);
    vec3 v3 = v + 0.5 * h * a2;  vec3 a3 = accel(p + 0.5 * h * v2, v3);
    vec3 v4 = v + h * a3;        vec3 a4 = accel(p + h * v3,       v4);
    vec3 pn = p + (h / 6.0) * (v + 2.0 * v2 + 2.0 * v3 + v4);
    vec3 vn = v + (h / 6.0) * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
    if (gInside || insideAny(pn)) break;               // 落入视界
    if (p.y * pn.y < 0.0) {                            // 穿过盘平面 y = 0
      vec3 hit = mix(p, pn, p.y / (p.y - pn.y));
      for (int i = 0; i < 4; i++) {
        vec3 e; float a;
        if (diskShade(i, hit, v, e, a)) { col += T * e; T *= 1.0 - a; }
      }
      if (T < 0.01) break;
    }
    p = pn; v = vn;
    if (dot(p, p) > esc2 && dot(p, v) > 0.0) {         // 逃逸：采样无穷远处的星空
      col += T * background(normalize(v));
      break;
    }
  }
  return col;
}

void main() {
  vec2 px = 1.0 / uResolution;
  vec2 base = gl_FragCoord.xy * px * 2.0 - 1.0;
  vec3 sum = vec3(0.0);
  for (int j = 0; j < 2; j++) for (int i = 0; i < 2; i++) {
    if (i >= uSupersample || j >= uSupersample) continue;
    vec2 off = (vec2(float(i), float(j)) + 0.5) / float(uSupersample) - 0.5;   // 子像素偏移（像素单位）
    sum += trace(cameraRay(base + off * 2.0 * px));
  }
  gl_FragColor = vec4(sum / float(uSupersample * uSupersample), 1.0);
}
`;
