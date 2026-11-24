/**
 * 《视界终舞》主程序：Three.js 渲染管线、OrbitControls、回放时钟，以及把轨道状态写进着色器。
 *
 * 长度单位：总质量对应的几何长度 GM/c²（下文记作“M”）；时间单位：GM/c³。
 * 画面的几何形状与总质量无关；只有把时间换算成秒、把频率换算成赫兹时才用到 M☉。
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import { buildTrajectory, ISCO, T_UNIT, L_UNIT, OMEGA_QNM, DAMP_QNM } from './physics/eob.js';
import { vertexShader, fragmentShader } from './shaders/raytracer.js';
import { buildBlackbodyTexture } from './blackbody.js';
import { createUI } from './ui.js';

const FOV = 46;
const R_CAVITY = 2.0 * 10;       // 环双星盘空腔半径：约为初始间距的 2 倍
const C_SOUND = 0.035;           // 环双星盘的声速（单位 c），对应几何薄盘 H/r ≈ 0.05–0.1
const QNM_AMPLITUDE = 0.16;      // 铃宕开始时视界的四极形变幅度
export const HOLD = 2.5;         // 终幕停留秒数，之后循环
const DIST_INSPIRAL = 62;        // 自动运镜：旋近期与并合期的相机距离
const DIST_MERGE = 23;
const PRESETS = {                // 画质：内部分辨率倍率、最大积分步数、积分步长倍率、超采样
  low: { scale: 0.45, steps: 420, step: 0.6, ss: 1 },
  medium: { scale: 0.7, steps: 700, step: 0.8, ss: 1 },
  high: { scale: 0.5, steps: 1100, step: 1.0, ss: 2 },       // 精细档：每个像素 2×2 条光线，用超采样换抗锯齿
};

const clamp = (x, a, b) => Math.min(Math.max(x, a), b);
const smoothstep = (a, b, x) => { const t = clamp((x - a) / (b - a), 0, 1); return t * t * (3 - 2 * t); };
const params = new URLSearchParams(location.search);
const num = (key, fallback) => { const v = parseFloat(params.get(key)); return Number.isFinite(v) ? v : fallback; };

/** 薄盘温度分布 T ∝ [(M/r³)(1 − √(r_in/r))]^{1/4} 的峰值，位于 r = (49/36) r_in */
const tNormPeak = (m, rin) => Math.pow((0.05666 * m) / rin ** 3, 0.25);
const T_SCALE = 6500 / tNormPeak(0.5, ISCO * 0.5);   // 显示色温的绝对标度：等质量微盘的峰值 ≈ 6500 K

const state = {
  q: clamp(num('q', 1), 1, 12), massMsun: num('M', 60) || 60, rate: 1,
  quality: PRESETS[params.get('quality')] ? params.get('quality') : 'auto',
  showDisk: params.get('disk') !== '0', autoCam: params.get('cam') !== '0',
  playing: params.get('paused') !== '1', s: Math.max(num('t', 0), 0),
};
let traj = buildTrajectory({ q: state.q });

// ————————————————— 渲染管线 —————————————————

const stage = document.getElementById('stage');
let renderer;
try {
  renderer = new THREE.WebGLRenderer({ antialias: false, powerPreference: 'high-performance' });
} catch (err) {
  document.getElementById('fatalMsg').textContent = String(err?.message ?? err);
  document.getElementById('fatal').hidden = false;
  throw err;
}
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.15;
// 着色器编译失败时把驱动的报错显示出来，而不是留下一块黑屏
renderer.debug.onShaderError = (gl, program, vs, fs) => {
  const log = [gl.getProgramInfoLog(program), gl.getShaderInfoLog(vs), gl.getShaderInfoLog(fs)].filter(Boolean).join('\n');
  document.getElementById('fatalMsg').textContent = `着色器编译失败：\n${log}`;
  document.getElementById('fatal').hidden = false;
};

const uniformArray = (Ctor, n) => ({ value: Array.from({ length: n }, () => new Ctor()) });
const uniforms = {
  uResolution: { value: new THREE.Vector2() },
  uCamPos: { value: new THREE.Vector3() },
  uCamBasis: { value: new THREE.Matrix3() },
  uTanHalfFov: { value: Math.tan(THREE.MathUtils.degToRad(FOV / 2)) },
  uTau: { value: 0 },
  uStepScale: { value: 0.8 },
  uMaxSteps: { value: 700 },
  uSupersample: { value: 1 },
  uEscRadius: { value: 300 },
  uBh: uniformArray(THREE.Vector4, 3),
  uQnm: { value: new THREE.Vector2() },
  uBlackbody: { value: buildBlackbodyTexture() },
  uDiskC: uniformArray(THREE.Vector3, 4),
  uDiskP: uniformArray(THREE.Vector4, 4),
  uDiskS: uniformArray(THREE.Vector3, 4),
  uDiskV: uniformArray(THREE.Vector3, 4),
  uSpiral: { value: new THREE.Vector2() },
  uPulse: { value: new THREE.Vector3(0, 0, 1) },
};

// 整个画面来自一个全屏片元着色器：正交相机 + 单位四边形
const scene = new THREE.Scene();
const quad = new THREE.Mesh(
  new THREE.PlaneGeometry(2, 2),
  new THREE.ShaderMaterial({ uniforms, vertexShader, fragmentShader, depthTest: false, depthWrite: false }),
);
quad.frustumCulled = false;
scene.add(quad);
const ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);

const composer = new EffectComposer(renderer);
composer.addPass(new RenderPass(scene, ortho));
const bloom = new UnrealBloomPass(new THREE.Vector2(1, 1), 0.7, 0.45, 1.1);
// UnrealBloomPass 对孤立的高亮像素会在最粗一级模糊上留下方形亮块：限制送进泛光的亮度上限
{
  const mat = bloom.materialHighPassFilter, src = 'vec4 texel = texture2D( tDiffuse, vUv );';
  if (mat.fragmentShader.includes(src)) {
    mat.fragmentShader = mat.fragmentShader.replace(src, `${src}\n\ttexel.rgb = min( texel.rgb, vec3( 8.0 ) );`);
    mat.needsUpdate = true;
  }
}
composer.addPass(bloom);
composer.addPass(new OutputPass());

// 观众操控的透视相机只提供位置与朝向，由 OrbitControls 驱动
const cam = new THREE.PerspectiveCamera(FOV, 1, 0.1, 5000);
const controls = new OrbitControls(cam, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.06;
controls.rotateSpeed = 0.7;
controls.minDistance = 5;
controls.maxDistance = 600;
controls.autoRotate = true;
controls.autoRotateSpeed = 0.22;
{
  const el = THREE.MathUtils.degToRad(15), az = 0.55;               // 略高于盘面，带一点侧视以显出多普勒不对称
  cam.position.set(Math.sin(az) * Math.cos(el), Math.sin(el), Math.cos(az) * Math.cos(el)).multiplyScalar(DIST_INSPIRAL);
  controls.update();
}
stage.appendChild(renderer.domElement);

let renderScale = 0.7;
function resize() {
  const w = window.innerWidth, h = window.innerHeight;
  const dpr = Math.min(window.devicePixelRatio || 1, 1.5);
  renderer.setPixelRatio(dpr);
  renderer.setSize(w, h);
  cam.aspect = w / h;
  cam.updateProjectionMatrix();
  composer.setPixelRatio(dpr * renderScale);
  composer.setSize(w, h);
  uniforms.uResolution.value.set(composer.renderTarget1.width, composer.renderTarget1.height);
}
window.addEventListener('resize', resize);
if (params.has('scale')) renderScale = clamp(num('scale', 0.5), 0.05, 2);   // 主要供截图与测试用
resize();

// ————————————————— 把轨道状态写进着色器 —————————————————

function updateWorld(st) {
  const { m1, m2, mFin } = traj;
  const w = st.w, r = st.r, phi = st.phi, om = st.om;
  const c = Math.cos(phi), s = Math.sin(phi);
  const r1 = r * m2, r2 = r * m1;                                    // 各自到质心的距离

  const [b0, b1, b2] = uniforms.uBh.value;
  b0.set(-r1 * c, 0, -r1 * s, m1 * (1 - w));                         // 并合过渡中，质量平滑地转给残骸
  b1.set(r2 * c, 0, r2 * s, m2 * (1 - w));
  b2.set(0, 0, 0, mFin * w);

  const dt = st.tau - traj.tauRing;
  const shape = dt >= 0 ? Math.exp(-dt * DAMP_QNM) : w;              // 四极形变：过渡中长起，之后按 QNM 衰减
  uniforms.uQnm.value.set(QNM_AMPLITUDE * shape, 2 * traj.phiRing + OMEGA_QNM * Math.max(dt, 0));

  const P = uniforms.uDiskP.value, C = uniforms.uDiskC.value, S = uniforms.uDiskS.value, V = uniforms.uDiskV.value;
  for (const p of P) p.w = 0;
  if (state.showDisk) {
    // 微盘：外缘被伴星的潮汐力限制在约 0.45a；收缩到 ISCO 附近时，间隙里的气体已被耗尽
    const mini = (i, mass, cx, cz, vx, vz) => {
      const rin = ISCO * mass, rout = 0.45 * r;
      C[i].set(cx, 0, cz); V[i].set(vx, 0, vz);
      P[i].set(rin, rout, mass, smoothstep(0.10, 0.45, rout / rin - 1) * (1 - w));
      S[i].set(T_SCALE, 0.3 + 1.7 * i, 0);
    };
    mini(0, m1, -r1 * c, -r1 * s, r1 * om * s, -r1 * om * c);       // 盘随黑洞平动；v = d(pos)/dτ
    mini(1, m2, r2 * c, r2 * s, -r2 * om * s, r2 * om * c);

    // 并合后的内盘：两个微盘的残余气体落向残骸，在几十 M 内建立
    C[2].set(0, 0, 0); V[2].set(0, 0, 0);
    P[2].set(ISCO * mFin, R_CAVITY * 0.62, mFin, w * smoothstep(0, 14, dt));
    S[2].set(T_SCALE, 4.1, 0);

    // 环双星盘：内缘被双星的潮汐转矩截断。粘滞时标远长于并合时标，所以内缘不随并合移动
    C[3].set(0, 0, 0); V[3].set(0, 0, 0);
    P[3].set(R_CAVITY, R_CAVITY * 2.0, 1, 1);
    S[3].set(T_SCALE, 8.3, 1);
    uniforms.uSpiral.value.set(phi, 5.0);

    // 并合瞬间中心质量骤减约 5%：环双星盘内缘被抛起，一道声波向外传播
    const after = Math.max(dt, 0);
    uniforms.uPulse.value.set(R_CAVITY * 1.2 + C_SOUND * after, dt > 0 ? 0.55 * Math.exp(-after / 70) : 0, 1.4 + 0.012 * after);
  }
  uniforms.uTau.value = st.tau;
}

// ————————————————— 观测量的真实单位换算 —————————————————
// 总质量只在这里出现：1 M☉ 的几何时间 GM/c³ ≈ 4.93 μs

function measure(st) {
  const tScale = T_UNIT * state.massMsun;                            // 秒 / M
  const ringing = st.tau >= traj.tauRing;
  // 旋近阶段四极辐射的频率是轨道频率的两倍；铃宕时应变以准正规模式的角频率 ω_QNM 本身振荡
  const fGw = ringing ? OMEGA_QNM / (2 * Math.PI * tScale) : st.om / (Math.PI * tScale);
  return {
    st, ringing, fGw,
    fOrb: ringing ? NaN : fGw / 2,
    lum: st.flux * L_UNIT,                             // 轨迹表给出的通量：旋近 → 峰值 → 铃宕衰减
    tPhys: st.tau * tScale,
  };
}

// ————————————————— 回放时钟与运镜 —————————————————

const clock = new THREE.Clock(false);
let lastUserInput = -1e9, autoDist = DIST_INSPIRAL, slowdown = 0, prevTau = NaN;

function advance(ds, trajEnd) {
  state.s += ds;
  if (state.s > trajEnd + HOLD) state.s = 0;                         // 循环
}
const fadeAt = (s, end) => Math.max(1 - smoothstep(0, 0.7, s), smoothstep(end + HOLD - 0.6, end + HOLD, s));

/** 并合前约 6 秒把镜头缓缓推近；用户操作后 8 秒内不接管 */
function applyAutoCamera(dt) {
  if (!state.autoCam || clock.elapsedTime - lastUserInput < 8) {
    autoDist = cam.position.distanceTo(controls.target);
    return;
  }
  const target = DIST_INSPIRAL + (DIST_MERGE - DIST_INSPIRAL) * smoothstep(traj.sTrig - 6, traj.sTrig, state.s);
  autoDist += (target - autoDist) * (1 - Math.exp(-dt / 1.1));
  const dir = cam.position.clone().sub(controls.target).normalize();
  cam.position.copy(controls.target).addScaledVector(dir, autoDist);
}

let ema = 16, frameCount = 0;
function adaptQuality(dt) {
  const preset = PRESETS[state.quality] ?? PRESETS.medium;
  uniforms.uStepScale.value = preset.step;
  uniforms.uMaxSteps.value = preset.steps;
  uniforms.uSupersample.value = preset.ss;
  if (params.has('scale')) return;                                   // 固定分辨率时不做自适应
  if (state.quality !== 'auto') {
    if (renderScale !== preset.scale) { renderScale = preset.scale; resize(); }
    return;
  }
  ema = ema * 0.92 + Math.min(dt * 1000, 60) * 0.08;                 // 帧时间的指数平均（毫秒）
  if (frameCount % 30 !== 0) return;
  const next = ema > 25 ? Math.max(renderScale * 0.9, 0.3) : ema < 13.5 ? Math.min(renderScale * 1.05, 1.0) : renderScale;
  if (next !== renderScale) { renderScale = next; resize(); }
}
controls.addEventListener('start', () => { lastUserInput = clock.elapsedTime; });

// ————————————————— 界面 —————————————————

const app = {
  state,
  getTraj: () => traj,
  setQ(q) { state.q = q; traj = buildTrajectory({ q }); state.s = 0; },
  setMass(m) { state.massMsun = m; },
  setRate(r) { state.rate = r; },
  setQuality(name) { state.quality = name; },
  setShowDisk(on) { state.showDisk = on; },
  setAutoCam(on) { state.autoCam = on; },
  togglePlay() { state.playing = !state.playing; },
  reset() { state.s = 0; },
  seek(s) { state.s = clamp(s, 0, traj.sEnd + HOLD); },
};
const ui = createUI(app);
if (params.get('ui') === '0') document.body.classList.add('clean');

// ————————————————— 主循环 —————————————————

function frame() {
  requestAnimationFrame(frame);
  const dt = Math.min(clock.getDelta(), 0.05);
  frameCount++;
  if (state.playing) advance(dt * state.rate, traj.sEnd);

  const st = traj.sample(Math.min(state.s, traj.sEnd));
  const scale = T_UNIT * state.massMsun;
  if (st.tau > prevTau) slowdown = (dt * state.rate) / ((st.tau - prevTau) * scale);   // 回放比真实时间慢多少倍
  prevTau = state.playing ? st.tau : NaN;

  controls.update();
  applyAutoCamera(dt);
  cam.updateMatrixWorld();
  const e = cam.matrixWorld.elements;
  uniforms.uCamBasis.value.set(e[0], e[4], -e[8], e[1], e[5], -e[9], e[2], e[6], -e[10]);   // 列：right, up, forward
  uniforms.uCamPos.value.copy(cam.position);
  uniforms.uEscRadius.value = Math.max(260, cam.position.length() * 2.5);
  updateWorld(st);

  composer.render();
  adaptQuality(dt);
  ui.update({ ...measure(st), slow: slowdown, fade: fadeAt(state.s, traj.sEnd), screen: state.s, distance: cam.position.length() });
}

document.documentElement.dataset.ready = '1';
window.__horizons = { state, getTraj: () => traj, ui, renderer, uniforms, bloom };   // 自动化测试与调试入口
clock.start();
frame();
