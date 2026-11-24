/**
 * 界面：数据面板、时间轴、引力波波形图和控制项。只负责显示与输入，不含物理计算。
 */

const $ = (id) => document.getElementById(id);
const clamp01 = (x) => Math.min(Math.max(x, 0), 1);
const fmtHz = (f) => (Number.isFinite(f) ? (f >= 1000 ? `${(f / 1000).toFixed(2)} kHz` : `${f.toFixed(1)} Hz`) : '—');
const fmtTime = (t) => {
  if (t < 1e-6) return `${(t * 1e9).toFixed(1)} ns`;
  if (t < 1e-3) return `${(t * 1e6).toFixed(t < 1e-4 ? 1 : 0)} μs`;
  if (t < 1) return `${(t * 1e3).toFixed(1)} ms`;
  return `${t.toFixed(2)} s`;
};
/** 1.3×10⁵⁰ 这样的功率/能量表示 */
const pow10 = (x) => {
  if (!(x > 0) || !Number.isFinite(x)) return '0';
  const e = Math.floor(Math.log10(x));
  return `${(x / 10 ** e).toFixed(2)}×10<sup>${e}</sup>`;
};
const expTimes = (n) => (n >= 100 ? `${Math.round(n).toLocaleString('en-US')}×` : `${n.toFixed(1)}×`);

const IDS = [
  'vPlay', 'vReset', 'vTimeline', 'vFill', 'vMarkTrig', 'vMarkRing', 'vPhase',
  'vTime', 'vSlow', 'vSep', 'vFgw', 'vForb', 'vSpeed', 'vLum', 'vDist',
  'vQ', 'vQOut', 'vMass', 'vRate', 'vRateOut', 'vQuality', 'vDisk', 'vAutoCam', 'vStatic', 'vHide', 'vHelp',
  'wave', 'about', 'vAboutClose', 'fade',
];

export function createUI(app) {
  const el = Object.fromEntries(IDS.map((id) => [id, $(id)]));
  const traj = app.getTraj;
  const ctx = el.wave.getContext('2d');
  const WAVE_WINDOW = 12;                                 // 波形图显示最近多少秒
  let dragging = false, lastText = 0, hidden = false;

  // ——— 控制项 ———
  const syncPlay = () => { el.vPlay.textContent = app.state.playing ? '❚❚' : '▶'; el.vPlay.setAttribute('aria-label', app.state.playing ? '暂停' : '播放'); };
  el.vPlay.addEventListener('click', () => { app.togglePlay(); syncPlay(); });
  el.vReset.addEventListener('click', () => app.reset());
  el.vQ.addEventListener('input', () => {
    app.setQ(Math.exp((Math.log(12) * parseFloat(el.vQ.value)) / 100));           // 对数刻度 1 … 12
    el.vQOut.textContent = `1 : ${app.state.q < 9.95 ? app.state.q.toFixed(1) : Math.round(app.state.q)}`;
    syncTimeline(); syncStatic();
  });
  el.vMass.addEventListener('change', () => { app.setMass(parseFloat(el.vMass.value)); syncStatic(); });
  el.vRate.addEventListener('input', () => {
    app.setRate(2 ** parseFloat(el.vRate.value));
    el.vRateOut.textContent = `×${app.state.rate < 10 ? app.state.rate.toFixed(2).replace(/\.?0+$/, '') : app.state.rate}`;
  });
  el.vQuality.addEventListener('change', () => app.setQuality(el.vQuality.value));
  el.vDisk.addEventListener('change', () => app.setShowDisk(el.vDisk.checked));
  el.vAutoCam.addEventListener('change', () => app.setAutoCam(el.vAutoCam.checked));
  const toggleBare = () => { hidden = !hidden; document.body.classList.toggle('bare', hidden); };
  el.vHide.addEventListener('click', toggleBare);
  el.vHelp.addEventListener('click', () => el.about.showModal());
  el.vAboutClose.addEventListener('click', () => el.about.close());
  el.about.addEventListener('click', (ev) => { if (ev.target === el.about) el.about.close(); });   // 点击遮罩关闭
  window.addEventListener('keydown', (ev) => {
    if (ev.target instanceof HTMLInputElement || ev.target instanceof HTMLSelectElement) return;
    if (ev.code === 'Space') { ev.preventDefault(); app.togglePlay(); syncPlay(); }
    else if (ev.key === 'r' || ev.key === 'R') app.reset();
    else if (ev.key === 'h' || ev.key === 'H') toggleBare();
  });
  // 时间轴：按下或拖动跳转
  const seekFrom = (ev) => {
    const rect = el.vTimeline.getBoundingClientRect();
    app.seek(clamp01((ev.clientX - rect.left) / rect.width) * traj().sEnd);
  };
  el.vTimeline.addEventListener('pointerdown', (ev) => { dragging = true; el.vTimeline.setPointerCapture(ev.pointerId); seekFrom(ev); });
  el.vTimeline.addEventListener('pointermove', (ev) => { if (dragging) seekFrom(ev); });
  el.vTimeline.addEventListener('pointerup', () => { dragging = false; });
  el.vTimeline.addEventListener('pointercancel', () => { dragging = false; });
  el.vQ.value = (Math.log(app.state.q) / Math.log(12)) * 100;
  el.vQOut.textContent = `1 : ${app.state.q < 9.95 ? app.state.q.toFixed(1) : Math.round(app.state.q)}`;
  syncPlay();

  // ——— 不随播放变化、只随质量比与总质量变化的信息 ———
  function syncStatic() {
    const t = traj(), M = app.state.massMsun, c2 = 1.7871e47;                     // 1 M☉c² ≈ 1.79×10⁴⁷ J
    const eJoule = t.eRad * M * c2;
    const nu = t.nu;
    el.vStatic.innerHTML = [
      ['对称质量比 ν = m₁m₂/M²', nu.toFixed(3)],
      ['初始间距', `10 M = ${(10 * 1.4766 * M).toFixed(0)} km`],
      ['并合时刻（真实时间）', fmtTime(t.tauTrig * 4.9255e-6 * M)],
      ['辐射的能量', `${(t.eRad * M).toFixed(2)} M☉c² ≈ ${pow10(eJoule)} J`],
      ['残骸', `M<sub>f</sub> = ${t.mFin.toFixed(3)} M，χ<sub>f</sub> = ${t.spin.toFixed(3)}`],
    ].map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
  }
  function syncTimeline() {
    const t = traj();
    el.vMarkTrig.style.left = `${(t.sTrig / t.sEnd) * 100}%`;
    el.vMarkRing.style.left = `${(t.sRing / t.sEnd) * 100}%`;
    const [labInspiral, labTrig, labRing] = el.vTimeline.parentElement.querySelectorAll('.ticks span');
    labInspiral.style.left = '0';
    labTrig.style.left = `${(t.sTrig / t.sEnd) * 100}%`;
  }
  syncStatic();
  syncTimeline();

  function fitCanvas() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    el.wave.width = Math.max(1, Math.round(el.wave.clientWidth * dpr));
    el.wave.height = Math.max(1, Math.round(el.wave.clientHeight * dpr));
  }
  window.addEventListener('resize', fitCanvas);
  fitCanvas();

  /** 滚动的引力波应变 h(t)，l = m = 2 模式 */
  function drawWave(screen) {
    const t = traj(), W = el.wave.width, H = el.wave.height;
    if (W < 4) return;
    ctx.clearRect(0, 0, W, H);
    ctx.strokeStyle = 'rgba(255,255,255,0.12)';
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(0, H / 2); ctx.lineTo(W, H / 2); ctx.stroke();
    const i1 = Math.min(Math.floor(screen * t.rate), t.count - 1), span = Math.round(WAVE_WINDOW * t.rate);
    const i0 = Math.max(i1 - span, 0), x = (i) => W - ((i1 - i) / span) * W, h = t.tbl.h;
    if (t.sTrig * t.rate >= i0 && t.sTrig * t.rate <= i1) {                     // 并合时刻
      const xm = x(t.sTrig * t.rate);
      ctx.strokeStyle = 'rgba(255,177,92,0.55)';
      ctx.beginPath(); ctx.moveTo(xm, 0); ctx.lineTo(xm, H); ctx.stroke();
    }
    ctx.strokeStyle = '#6fe3ff';
    ctx.lineWidth = Math.max(1.2, el.wave.width / 520);
    ctx.beginPath();
    for (let i = i0; i <= i1; i++) { const y = H / 2 - h[i] * (H * 0.44 / 0.34); if (i === i0) ctx.moveTo(x(i), y); else ctx.lineTo(x(i), y); }
    ctx.stroke();
  }

  function phaseName(st) {
    const t = traj();
    if (st.tau >= t.tauRing) return '铃宕';
    if (st.tau >= t.tauTrig) return '并合';
    return st.r <= t.bar.stop ? '坠落' : '旋近';
  }

  /** 每帧调用 */
  function update(info) {
    const now = performance.now();
    el.fade.style.opacity = info.fade.toFixed(3);
    el.vFill.style.width = `${clamp01(info.screen / traj().sEnd) * 100}%`;
    drawWave(info.screen);
    const t = traj(), scale = 4.9255e-6 * app.state.massMsun;
    el.vPhase.textContent = phaseName(info.st);
    el.vPhase.dataset.stage = info.st.tau >= t.tauRing ? 'ring' : info.st.tau >= t.tauTrig ? 'trig' : info.st.r <= t.bar.stop ? 'fall' : 'inspiral';
    if (now - lastText < 90) return;                       // 数值文字以约 11 Hz 刷新，避免闪烁
    lastText = now;
    el.vTime.innerHTML = `${fmtTime(info.tPhys)} <small>= ${info.st.tau.toFixed(0)} GM/c³</small>`;
    el.vSlow.textContent = info.slow > 0 ? `${expTimes(info.slow)} 慢放` : '—';
    el.vSep.textContent = info.ringing ? '—' : `${info.st.r.toFixed(2)} M`;
    el.vForb.textContent = info.ringing ? '—' : fmtHz(info.fOrb);
    el.vFgw.textContent = fmtHz(info.fGw);
    el.vSpeed.textContent = info.ringing ? '—' : `${Math.cbrt(info.st.om).toFixed(3)} c`;
    el.vLum.innerHTML = `${pow10(info.lum)} W<small>${(info.lum / 3.628e52).toFixed(4)} c⁵/G</small>`;
    el.vDist.textContent = `${info.distance.toFixed(1)} M`;
    syncPlay();
  }
  return { update };
}
