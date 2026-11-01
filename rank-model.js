/* ============================================================
   素数あつめ ランキングの「偏差値」（概算）と、公開データの形（v3）
   ・RATING（index.html）と ランキング分析（ranking-analysis.html）の両方で使います
   ・画面を持たない、計算だけの部分です（Node.js でも そのまま動きます）

   偏差値の出し方
     1. 上位の割合     p ＝（順位 − 0.5）÷ 参加人数 N
     2. 標準正規分布で「上から p の所」にある値     z ＝ Φ⁻¹(1 − p)
     3. 偏差値 ＝ 50 ＋ 10 × z
   記録 x は、まず y ＝ g(x) に直します（点数 g(x) ＝ x^(3/4)、タイム g(t) ＝ −√t、累積 g(x) ＝ x^0.3）。
   読めた順位のあいだは、y と z を まっすぐ結んで見つもります。
   読めた順位より下は「予測分布」で出します：y は つりがね形（正規分布）に ちらばる、と考えて
     y ＝ μ ＋ σz、z ＝（y − μ）÷ σ、順位 ≈ N × Q(z) ＋ 0.5（Q は上側確率）
     σ：読めた順位の いちばん下の点 L から、z の幅 0.4 の中の点で「y を z の直線で表した傾き」（最小二乗）。
        点数と累積は「最下位（N 位）の人は 0 くらい」とするので σ ≥ y_L ÷（z_L − z_N）
     μ：直線が L を通るように μ ＝ y_L − σ z_L
   予測平均点・予測標準偏差：参加者 N 人の記録を 1 位から N 位まで並べたもの x̂₁ … x̂_N の 平均と標準偏差です。
     x̂_r は、読めた順位は その記録、読めた順位のあいだは 上の まっすぐ結んだ値、
     読めた順位より下は x̂_r ＝ g⁻¹(μ ＋ σ z_r)（μ ＋ σ z_r が 0 より下なら 0）。
   ============================================================ */
(function (G) {
'use strict';

/* ---------- 標準正規分布 ---------- */
/** 相補誤差関数 erfc(x)。相対誤差は どこでも 1.2×10⁻⁷ より小さい（Numerical Recipes の erfcc） */
function erfc(x) {
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
    t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
/** 上側確率 Q(z) ＝ 1 − Φ(z)（z より上にいる人の割合） */
const upper = z => z === Infinity ? 0 : z === -Infinity ? 1 : 0.5 * erfc(z / Math.SQRT2);
/** 下側確率 Φ(z) */
const normCdf = z => 1 - upper(z);
/** Φ⁻¹(p)：下側確率が p になる z（Acklam の近似。相対誤差 1.15×10⁻⁹） */
function normInv(p) {
  if (!(p > 0)) return -Infinity;
  if (!(p < 1)) return Infinity;
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  const lo = 0.02425;
  if (p < lo) {
    const q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  if (p > 1 - lo) return -normInv(1 - p);
  const q = p - 0.5, r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}
/** 上側確率が q になる z（q が とても小さくても、1 − q の丸めで ずれないようにしています） */
const zOfUpper = q => -normInv(q);
/** 参加 N 人の中の r 位 → z */
const zOfRank = (r, N) => zOfUpper((r - 0.5) / N);
/** z → 偏差値、偏差値 → z */
const tOfZ = z => 50 + 10 * z;
const zOfT = t => (t - 50) / 10;
/** 順位 → 偏差値 */
const tOfRank = (r, N) => tOfZ(zOfRank(r, N));
/** 偏差値 → 順位（小数のまま。ちょうど その偏差値になる順位） */
const rankOfT = (t, N) => N * upper(zOfT(t)) + 0.5;

/* ---------- 記録の見方（y ＝ g(x)） ----------
   pow：点数の 3/4 乗（素数ゲームの点・スコアの合計）。0 点 → 0。予測平均点・予測標準偏差の単位は「点」
   nsqrt：タイムの √（1/2 乗）に マイナスを付けたもの（小さいほど良いので）。単位は「秒」
   pcum：累積の 0.3 乗（生成した乱数の累積・スコアの累積）。0 → 0
   lin・neg・log10：v2 までの形（そのまま・タイムに マイナス・常用対数）。くらべる時だけ使います
   べき（3/4・1/2・0.3）は、2026-11-01 の 1〜1000 位で「上位 200・300・500 位だけで当てはめて、その下を当てる」こと、
   一覧より下の記録（非公開）を当てること、累積は 2026-10-31 に読めた 1480〜10000 位を当てることを ためして、
   いちばん よく当たったものです */
const POW = { score: 0.75, time: 0.5, cum: 0.3 };
const powTR = (k, unit) => ({ f: v => Math.pow(v, k), inv: y => y > 0 ? Math.pow(y, 1 / k) : 0, ok: v => v >= 0, floor: true, unit, pow: k });
const TR = {
  pow: powTR(POW.score, '点'),
  nsqrt: { f: v => -Math.pow(v, POW.time), inv: y => y < 0 ? Math.pow(-y, 1 / POW.time) : 0, ok: v => v > 0, unit: '秒', pow: POW.time },
  pcum: powTR(POW.cum, ''),
  lin: { f: v => v, inv: y => y > 0 ? y : 0, ok: v => v >= 0, floor: true, unit: '点', pow: 1 },
  neg: { f: v => -v, inv: y => -y, ok: v => v > 0, unit: '秒', pow: 1 },
  log10: { f: v => v > 0 ? Math.log10(v) : -Infinity, inv: y => Math.pow(10, y), ok: v => v >= 0, unit: '桁', log: true }
};
/** いちばん下の所の傾きを出すのに使う、z の幅 */
const TAIL_SPAN = 0.4;
/** 確かさの しるし */
const MARK = { exact: '◎', top: '◎', interp: '○', above: '○', extrap: '△', below: '—', over: '—', under: '—' };
const CONF_TEXT = {
  exact: '読めた順位で、ぴったり決まります',
  top: 'いまの1位より上（または同じ）です',
  interp: '読めた順位のあいだを見つもっています',
  above: '読めた順位より上です',
  extrap: '読めた順位より下なので、予測分布から出しています',
  below: '見つもれません'
};

/** 最小二乗の傾き（y を z の直線で表した時の 傾き） */
function slope(list) {
  const n = list.length;
  if (n < 2) return NaN;
  let sz = 0, sy = 0;
  for (const p of list) { sz += p.z; sy += p.y; }
  const mz = sz / n, my = sy / n;
  let num = 0, den = 0;
  for (const p of list) { num += (p.z - mz) * (p.y - my); den += (p.z - mz) * (p.z - mz); }
  return den > 0 ? num / den : NaN;
}

/** 種目のモデルを作ります。
    cat = { name, asc（小さいほど良い）, kind（'cum' は累積）, players（参加人数）, pts: [[順位, 値], …], extra: [[順位, 値], …] }
    opt.N で参加人数を変えられます。opt.extra === false で 一覧の外の点を使いません。opt.tr で 記録の見方を変えられます */
function build(cat, opt) {
  opt = opt || {};
  if (!cat || !Array.isArray(cat.pts) || !cat.pts.length) return null;
  const N = opt.N > 0 ? Math.round(opt.N) : (cat.players > 0 ? Math.round(cat.players) : 0);
  if (!N) return null;
  const asc = !!cat.asc;
  const trName = TR[opt.tr] ? opt.tr : (asc ? 'nsqrt' : (cat.kind === 'cum' || /累積/.test(String(cat.name || '')) ? 'pcum' : 'pow'));
  const tr = TR[trName];
  const good = (a, b) => asc ? a <= b : a >= b;       // a は b と同じか、それより良い
  /* 一覧の点（並びが さかさまの所があれば、上の値にそろえます） */
  const seen = new Set(), list = [];
  for (const p of cat.pts) {
    const r = Number(p[0]), v = Number(p[1]);
    if (!Number.isInteger(r) || r < 1 || r > N || !isFinite(v) || !tr.ok(v) || seen.has(r)) continue;
    seen.add(r); list.push({ r, v, src: 'list' });
  }
  if (!list.length) return null;
  list.sort((a, b) => a.r - b.r);
  for (let i = 1; i < list.length; i++) if (!good(list[i - 1].v, list[i].v)) list[i].v = list[i - 1].v;
  /* 一覧の外の点（画面の「あなたの順位」の欄など）。前後の点と並びが合う時だけ使います */
  const extra = [];
  if (opt.extra !== false) {
    for (const e of (Array.isArray(cat.extra) ? cat.extra : [])) {
      const r = Math.round(Number(e && e[0])), v = Number(e && e[1]);
      if (!(r >= 1 && r <= N) || !isFinite(v) || !tr.ok(v) || seen.has(r)) continue;
      let above = null, below = null;
      for (const p of list) { if (p.r < r) above = p; else { below = p; break; } }
      if ((above && !good(above.v, v)) || (below && !good(v, below.v))) continue;
      seen.add(r); extra.push({ r, v, src: 'extra' });
    }
  }
  const all = list.concat(extra).sort((a, b) => a.r - b.r);
  const pts = [];
  for (const p of all) {
    /* 一覧の外の点どうしで並びが合わない時は、あとの点を使いません */
    if (pts.length && !good(pts[pts.length - 1].v, p.v)) {
      if (p.src === 'extra') continue;
      while (pts.length && pts[pts.length - 1].src === 'extra' && !good(pts[pts.length - 1].v, p.v)) pts.pop();
    }
    pts.push(p);
  }
  for (const p of pts) { p.z = zOfRank(p.r, N); p.y = tr.f(p.v); }
  let contig = 0;
  for (const p of list) { if (p.r === contig + 1) contig = p.r; else if (p.r > contig + 1) break; }
  /* いちばん下の所の傾き：下から点を足していき、z の幅が TAIL_SPAN に とどくまでの点で出します */
  const last = pts[pts.length - 1], tail = [];
  for (let i = pts.length - 1; i >= 0; i--) {
    tail.push(pts[i]);
    if (tail.length >= 2 && pts[i].z - last.z >= TAIL_SPAN) break;
  }
  /* 予測分布：y ＝ μ ＋ σz。σ は いちばん下の所の傾き、μ は 直線が いちばん下の点 L を通るように決めます。
     点数と累積（0 より下が無いもの）は「最下位（N 位）の人は 0 くらい」とするので、
     直線が z_N より上で 0 に とどくくらいの傾きは要ります：σ ≥ y_L ÷（z_L − z_N） */
  const zBottom = zOfRank(N, N);
  let sigma = slope(tail), clamp = '';
  if (!(sigma > 0) || !isFinite(sigma)) sigma = slope(pts);
  if (tr.floor && last.y > 0 && last.z > zBottom) {
    const lo = last.y / (last.z - zBottom);
    if (!(sigma >= lo) || !isFinite(sigma)) { sigma = lo; clamp = 'low'; }
  }
  if (!(sigma > 0) || !isFinite(sigma)) sigma = null;
  const mu = sigma ? last.y - last.z * sigma : null;
  const listMax = list[list.length - 1].r, extraUsed = pts.filter(p => p.src === 'extra').length;
  const m = {
    name: cat.name || '', N, asc, tr: trName, f: tr.f, inv: tr.inv, ok: tr.ok, floor: !!tr.floor, log: !!tr.log, unit: tr.unit,
    pts, last, contig, mu, sigma, b: sigma, clamp, listMax, extraUsed, tailFrom: tail[tail.length - 1].r, tailN: tail.length,
    zTop: zOfRank(1, N), zBottom, playersFrom: cat.playersFrom || null
  };
  /* 予測平均点・予測標準偏差は、はじめて使う時に 計算します（参加人数ぶん 足すので） */
  let dist;
  Object.defineProperty(m, 'dist', { enumerable: true, get() { if (dist === undefined) dist = sigma ? moments(m) : null; return dist; } });
  return m;
}

/** 読めた順位の中の z の所の y（まっすぐ結んだ値）。z が いちばん上の点より上なら いちばん上の点の y */
function yInside(m, z) {
  const P = m.pts;
  if (z >= P[0].z) return P[0].y;
  let lo = 0, hi = P.length - 1;            // P[lo].z > z >= P[hi].z
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (P[mid].z > z) lo = mid; else hi = mid; }
  const a = P[lo], b = P[hi];
  return b.z === a.z ? b.y : a.y + (z - a.z) / (b.z - a.z) * (b.y - a.y);
}
/** 参加者 N 人の「予測した記録の並び」x̂₁ … x̂_N の 平均・標準偏差（記録の見方が log10 の時だけ、y ＝ log₁₀ のまま）。
    読めた順位 r ≤ r_L：読めた記録（抜けている順位は まっすぐ結んだ値）
    読めた順位より下 r > r_L：x̂_r ＝ g⁻¹(μ ＋ σ z_r)、z_r ＝ Φ⁻¹(1 − (r − 0.5) ÷ N)（点数と累積は 0 より下を 0） */
function moments(m) {
  const N = m.N, P = m.pts, L = m.last;
  let n = 0, mean = 0, M2 = 0, zero = 0, worst = null, i = 0;
  for (let r = 1; r <= N; r++) {
    let y;
    if (r <= L.r) {
      while (P[i].r < r) i++;
      y = P[i].r === r ? P[i].y : (i === 0 ? P[0].y : yInside(m, zOfRank(r, N)));
    } else y = m.mu + m.sigma * zOfRank(r, N);
    let x = m.log ? y : m.inv(y);
    if (m.floor && !(x > 0)) { x = 0; if (r > L.r) zero++; }
    if (!isFinite(x)) continue;
    n++; const d = x - mean; mean += d / n; M2 += d * (x - mean);
    if (r === N) worst = x;
  }
  if (!n) return null;
  const sd = Math.sqrt(M2 / n);
  const med = valueAtZ(m, 0);
  if (m.log) return { mean: Math.pow(10, mean), sd, logMean: mean, median: med, unit: '桁', log: true, N, read: L.r, worst: worst == null ? null : Math.pow(10, worst) };
  return { mean, sd, median: med, unit: m.unit, N, read: L.r, worst, zero: zero / N, zeroN: zero };
}
/** z の所の 記録（読めた順位の中は まっすぐ結んだ値、下は 予測分布。点数は 0 点より下を 0 点） */
function valueAtZ(m, z) {
  if (!m) return null;
  const L = m.last;
  let y;
  if (z >= L.z) y = yInside(m, z);
  else if (m.sigma) y = m.mu + m.sigma * Math.max(z, m.zBottom);
  else return null;
  const x = m.inv(y);
  return m.floor && !(x > 0) ? 0 : x;
}
/** 予測分布の中での偏差値 ＝ 50 ＋ 10 ×（y − μ）÷ σ（y ＝ g(記録)。読めた順位より下では、偏差値と同じになります） */
function classicT(m, value) {
  if (!m || !m.sigma) return null;
  const v = Number(value);
  if (!isFinite(v) || !m.ok(v)) return null;
  return tOfZ((m.f(v) - m.mu) / m.sigma);
}
/** ふつうの偏差値 ＝ 50 ＋ 10 ×（記録 − 予測平均点）÷ 予測標準偏差（タイムは 速いほど上、累積は 対数で）。
    くらべるための参考です（点数の ちらばりが つりがね形でないと、上位% と合いません） */
function plainT(m, value) {
  const D = m && m.dist;
  if (!D || !(D.sd > 0)) return null;
  const v = Number(value);
  if (!isFinite(v) || !m.ok(v)) return null;
  if (D.log) return v > 0 ? tOfZ((Math.log10(v) - D.logMean) / D.sd) : null;
  return tOfZ((m.asc ? D.mean - v : v - D.mean) / D.sd);
}

/** その記録の 順位・上位%・偏差値。
    conf：exact（◎ ぴったり）／top（◎ 1位以上）／interp（○ 読めた順位のあいだ）／extrap（△ 読めた順位より下）／below（見つもれない） */
function locate(m, value) {
  if (!m) return null;
  const v = Number(value);
  if (!isFinite(v) || !m.ok(v)) return null;
  const y = m.f(v), P = m.pts, N = m.N;
  const out = (rank, z, conf, more) => Object.assign({ rank, z, T: tOfZ(z), top: rank / N * 100, conf, mark: MARK[conf], N }, more || {});
  /* 1) いちばん上の点と同じか、それより良い */
  if (y >= P[0].y) {
    if (P[0].r === 1) return out(1, zOfRank(1, N), y > P[0].y ? 'top' : 'exact');
    return out(P[0].r, P[0].z, 'above', { to: P[0].r });
  }
  /* 2) 読めた点のあいだ：P[j-1] より悪く、P[j] と同じか良い所 */
  if (y >= P[P.length - 1].y) {
    let lo = 0, hi = P.length - 1;          // P[lo].y > y >= P[hi].y
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (P[mid].y > y) lo = mid; else hi = mid; }
    const a = P[lo], b = P[hi];
    /* すぐ上の順位まで読めている時と、読めた点と同じ記録の時は、順位がぴったり決まります */
    if (b.r === a.r + 1 || y === b.y) return out(b.r, b.z, 'exact');
    const t = (a.y - y) / (a.y - b.y);
    const z = a.z + t * (b.z - a.z);
    const rank = Math.min(b.r, Math.max(a.r + 1, Math.round(N * upper(z) + 0.5)));
    return out(rank, z, 'interp', { from: a.r, to: b.r });
  }
  /* 3) 読めた点より下：予測分布で出します。z ＝（y − μ）÷ σ（最下位で止めます） */
  const L = m.last;
  if (!m.sigma) return { rank: null, z: null, T: null, top: null, conf: 'below', mark: MARK.below, N, below: L.r, rankMin: L.r + 1, tMax: tOfZ(L.z) };
  let z = y === -Infinity ? -Infinity : (y - m.mu) / m.sigma;
  if (!(z >= m.zBottom)) z = m.zBottom;
  if (z > L.z) z = L.z;
  const rank = Math.min(N, Math.max(L.r + 1, Math.round(N * upper(z) + 0.5)));
  return out(rank, z, 'extrap', { from: L.r });
}

/** その偏差値に とどくのに要る記録（タイムは「これより速い」、ほかは「これより上」）。
    value が null の時：over（1位でも とどかない）／under（いちばん下でも とどく） */
function valueForT(m, T) {
  if (!m || !isFinite(T)) return null;
  const z = zOfT(T), P = m.pts, N = m.N;
  const rk = N * upper(z) + 0.5;
  const res = (value, conf, more) => Object.assign({ T, rank: rk, value, conf, mark: MARK[conf] }, more || {});
  if (z > P[0].z) return res(null, 'over', { best: P[0].v, tBest: tOfZ(P[0].z) });
  if (z >= P[P.length - 1].z) {
    let lo = 0, hi = P.length - 1;          // P[lo].z >= z > P[hi].z（ちょうど同じ時は lo 側）
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (P[mid].z >= z) lo = mid; else hi = mid; }
    const a = P[lo], b = P[hi];
    if (z === a.z) return res(a.v, 'exact', { rank: a.r });
    if (z === b.z) return res(b.v, 'exact', { rank: b.r });
    if (b.r === a.r + 1) {
      /* 抜けのない所：a.r 位に並べば、その偏差値に とどきます */
      return res(a.v, 'exact', { rank: a.r });
    }
    const y = a.y + (z - a.z) / (b.z - a.z) * (b.y - a.y);
    return res(m.inv(y), 'interp', { from: a.r, to: b.r });
  }
  if (!m.sigma) return res(null, 'below', { below: m.last.r });
  if (z < m.zBottom) return res(null, 'under');
  let value = m.inv(m.mu + z * m.sigma), zero = false;
  if (m.floor && !(value > 0)) { value = 0; zero = true; }
  return res(value, 'extrap', { from: m.last.r, zero });
}

/* ---------- くらべ：「参加者 全員の記録が 正規分布」と考えた時 ----------
   上位 k 人の記録 x₁ ≥ … ≥ x_k と、残りの N − k 人は「x_k より下」ということから、
   正規分布 N(μ, σ²) の μ・σ を 最尤推定で決めます（打ち切りのある正規分布）：
     対数尤度 ℓ(μ, σ) ＝ −k log σ − Σ (x_i − μ)² ÷ (2σ²) ＋ (N − k) log Φ((x_k − μ) ÷ σ)（＋ 定数）
   タイムは −t で（速いほど大きく）考えます。累積は（何万倍も ちがうので）くらべません */
function logPhi(z) {
  if (z > -30) return Math.log(Math.max(1e-300, normCdf(z)));
  return -z * z / 2 - Math.log(-z) - 0.5 * Math.log(2 * Math.PI);   // とても小さい時（漸近形）
}
/** ネルダー–ミード法（f を小さくする x を さがします） */
function nelderMead(f, x0, step, iters) {
  const n = x0.length;
  let S = [x0.slice()].concat(x0.map((_, i) => { const p = x0.slice(); p[i] += step[i]; return p; })).map(p => ({ p, f: f(p) }));
  for (let it = 0; it < (iters || 2000); it++) {
    S.sort((a, b) => a.f - b.f);
    if (it > 50 && Math.abs(S[n].f - S[0].f) <= 1e-11 * (1 + Math.abs(S[0].f))) break;
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let d = 0; d < n; d++) c[d] += S[i].p[d] / n;
    const at = t => c.map((v, d) => v + t * (v - S[n].p[d]));
    const r = at(1), fr = f(r);
    if (fr < S[0].f) { const e = at(2), fe = f(e); S[n] = fe < fr ? { p: e, f: fe } : { p: r, f: fr }; }
    else if (fr < S[n - 1].f) S[n] = { p: r, f: fr };
    else {
      const k = at(-0.5), fk = f(k);
      if (fk < S[n].f) S[n] = { p: k, f: fk };
      else S = S.map((q, i) => { if (!i) return q; const p = q.p.map((v, d) => (v + S[0].p[d]) / 2); return { p, f: f(p) }; });
    }
  }
  S.sort((a, b) => a.f - b.f);
  return S[0];
}
/** 「全員が正規分布」の μ・σ（上位から抜けなく読めた所 1〜k 位を使います）。m は build の結果 */
function fitNormalAll(m) {
  if (!m || m.log || m.tr === 'pcum') return null;
  const xs = [];
  for (const p of m.pts) { if (p.r !== xs.length + 1) break; xs.push(m.asc ? -p.v : p.v); }
  const k = xs.length, N = m.N;
  if (k < 3 || k >= N) return null;
  const xk = xs[k - 1];
  const nll = ([mu, ls]) => {
    const sd = Math.exp(ls); let s = k * ls;
    for (const x of xs) { const e = (x - mu) / sd; s += e * e / 2; }
    return s - (N - k) * logPhi((xk - mu) / sd);
  };
  /* はじめの値：x を z の直線で表した時の 切片と傾き */
  const P = xs.map((x, i) => ({ z: zOfRank(i + 1, N), y: x }));
  const s0 = slope(P) > 0 ? slope(P) : 1, my = P.reduce((a, p) => a + p.y, 0) / k, mz = P.reduce((a, p) => a + p.z, 0) / k;
  const best = nelderMead(nll, [my - s0 * mz, Math.log(s0)], [s0 * 0.3, 0.2], 3000);
  const mu = best.p[0], sd = Math.exp(best.p[1]);
  /* 点数は 0 点より下が無いので、0 点より下の人は 0 点として 平均を出します（タイムは そのまま） */
  const a = mu / sd, floorMean = m.asc ? null : mu * normCdf(a) + sd * Math.exp(-a * a / 2) / Math.sqrt(2 * Math.PI);
  return { mu: m.asc ? -mu : mu, sd, k, floorMean, below0: m.asc ? null : normCdf(-a) };
}
/** 上位 K 位まで だけで当てはめて、その下の読めた順位を当てた時の 偏差値のずれ。
    how：'v3'（この計算）／'normal'（全員が正規分布）。返す値 { n, mae（平均の大きさ）, max, end（いちばん下の点の ずれ） } */
function holdout(cat, K, how) {
  const full = build(cat, { extra: false });
  if (!full) return null;
  const fitPts = cat.pts.filter(p => p[0] <= K);
  const m = build(Object.assign({}, cat, { pts: fitPts, extra: [] }), { N: full.N });
  if (!m) return null;
  let nf = null;
  if (how === 'normal') { nf = fitNormalAll(m); if (!nf) return null; }
  let n = 0, sa = 0, mx = 0, end = null;
  for (const p of full.pts) {
    if (p.r <= K) continue;
    const T = tOfRank(p.r, full.N);
    const Tp = how === 'normal' ? tOfZ((m.asc ? nf.mu - p.v : p.v - nf.mu) / nf.sd) : locate(m, p.v).T;
    if (Tp == null || !isFinite(Tp)) continue;
    const e = Tp - T; n++; sa += Math.abs(e); mx = Math.max(mx, Math.abs(e)); end = e;
  }
  return n ? { n, mae: sa / n, max: mx, end } : null;
}

/** 参加人数が読めていない種目に、前の週（日付が前で いちばん新しい週）の人数を入れます。
    入れた種目には playersFrom（人数を読んだ週の日付）を付けます。weeks は readWeek の結果（{ date, cats }）の配列です */
function fillPlayers(weeks) {
  const ws = (Array.isArray(weeks) ? weeks : []).filter(w => w && w.cats && w.date).slice().sort((a, b) => a.date < b.date ? -1 : a.date > b.date ? 1 : 0);
  let filled = 0;
  for (let i = 1; i < ws.length; i++) {
    for (const name of Object.keys(ws[i].cats)) {
      const c = ws[i].cats[name];
      if (c.players > 0) continue;
      for (let k = i - 1; k >= 0; k--) {
        const o = ws[k].cats[name];
        if (o && o.players > 0) { c.players = o.players; c.playersFrom = o.playersFrom || ws[k].date; filled++; break; }
      }
    }
  }
  return filled;
}

/* ---------- 公開データ（rankings/）の形 ---------- */
/** ランキング集計の JSON（ranking-ocr.html の書き出し・公開データ どちらでも）を読みます。
    返す種目名は JSON に書かれている名前のままです。opt.self === false で「あなたの順位」の欄を使いません */
function readWeek(j, opt) {
  opt = opt || {};
  if (!j || typeof j !== 'object' || !j.categories || typeof j.categories !== 'object') return null;
  const bad = new Set();
  [].concat(Array.isArray(j.checking) ? j.checking : [], Array.isArray(j.excluded) ? j.excluded : [])
    .forEach(c => { if (c && c.category != null) bad.add(c.category + ':' + c.rank); });
  const cats = {};
  for (const name of Object.keys(j.categories)) {
    const c = j.categories[name];
    if (!c || typeof c !== 'object' || !c.ranks || typeof c.ranks !== 'object') continue;
    const pts = [];
    for (const k of Object.keys(c.ranks)) {
      const r = Number(k), v = Number(c.ranks[k]);
      if (!Number.isInteger(r) || r < 1 || !isFinite(v) || v < 0 || bad.has(name + ':' + r)) continue;
      pts.push([r, v]);
    }
    if (!pts.length) continue;
    pts.sort((a, b) => a[0] - b[0]);
    const extra = [];
    if (Array.isArray(c.extra)) c.extra.forEach(e => { if (Array.isArray(e) && isFinite(+e[0]) && isFinite(+e[1])) extra.push([+e[0], +e[1]]); });
    if (opt.self !== false && c.self && +c.self.rank > 0 && isFinite(+c.self.score)) extra.push([+c.self.rank, +c.self.score]);
    const asc = c.order === 'asc' || c.unit === '秒';
    cats[name] = {
      name, asc, kind: asc ? 'time' : (/累積/.test(name) ? 'cum' : 'int'),
      players: +c.players > 0 ? Math.round(+c.players) : null, pts, extra
    };
  }
  const d = String(j.date || '');
  return { date: /^\d{4}-\d{2}-\d{2}$/.test(d) ? d : null, cats };
}
/** 公開する形にします（checking! の残った順位は入れません。名前は もともと入っていません）。
    opt.self：画面の「あなたの順位」の欄を、一覧の外の点（extra）として入れるか */
function toPublic(j, opt) {
  opt = opt || {};
  const w = readWeek(j, { self: false });
  if (!w || !w.date) throw new Error('日付の入った ランキング集計の JSON ではありません');
  const names = Object.keys(w.cats);
  if (!names.length) throw new Error('読み取れた順位が1つもありません');
  const out = {
    app: '素数あつめ ランキング集計', format: 2, public: true, date: w.date,
    createdAt: j.createdAt || null, images: j.images || null,
    note: '読み取れた順位と点数だけが入っています（名前は入っていません）。extra は、一覧の外で読めた［順位, 点数］です。',
    categories: {}
  };
  for (const name of names) {
    const c = w.cats[name], src = j.categories[name] || {};
    const ranks = {};
    c.pts.forEach(([r, v]) => { ranks[r] = v; });
    const o = {
      unit: src.unit || (c.asc ? '秒' : '点'), order: c.asc ? 'asc' : 'desc', players: c.players,
      count: c.pts.length, top: c.pts[c.pts.length - 1][0], ranks
    };
    const ex = c.extra.filter(e => !(e[0] in ranks));
    if (opt.self && src.self && +src.self.rank > 0 && isFinite(+src.self.score) && !(+src.self.rank in ranks)) ex.push([+src.self.rank, +src.self.score]);
    if (ex.length) o.extra = ex;
    out.categories[name] = o;
  }
  return out;
}
/** 公開データを文字にします（1種目を1行にするので、週ごとのちがいが見やすくなります） */
function stringify(pub) {
  const lines = ['{'];
  Object.keys(pub).filter(k => k !== 'categories').forEach(k => lines.push(`  ${JSON.stringify(k)}: ${JSON.stringify(pub[k])},`));
  lines.push('  "categories": {');
  const names = Object.keys(pub.categories || {});
  names.forEach((n, i) => {
    const c = pub.categories[n];
    const head = Object.keys(c).filter(k => k !== 'ranks').map(k => `${JSON.stringify(k)}: ${JSON.stringify(c[k])}`);
    head.push(`"ranks": ${JSON.stringify(c.ranks)}`);
    lines.push(`    ${JSON.stringify(n)}: {${head.join(', ')}}${i < names.length - 1 ? ',' : ''}`);
  });
  lines.push('  }', '}');
  return lines.join('\n') + '\n';
}
/** rankings/index.json に、その週を足します（同じ日付は入れかえます） */
function addToIndex(idx, pub, file) {
  const weeks = (idx && Array.isArray(idx.weeks) ? idx.weeks : []).filter(w => w && w.date && w.date !== pub.date);
  const names = Object.keys(pub.categories || {});
  weeks.push({
    date: pub.date, file: file || (pub.date + '.json'), categories: names.length,
    ranks: names.reduce((a, n) => a + (pub.categories[n].count || 0), 0)
  });
  weeks.sort((a, b) => a.date < b.date ? -1 : a.date > b.date ? 1 : 0);
  return {
    app: '素数あつめ ランキング集計', note: 'ランキングの全部の順位（ゲームの画面から読み取ったもの）の一覧です。週ごとのファイルは、この index.json と同じ所にあります。',
    latest: weeks[weeks.length - 1].date, weeks
  };
}
/** その日付に使う週（同じ日付、なければ その日より前で いちばん新しい週。それも無ければ いちばん新しい週） */
function pickWeek(idx, date) {
  const ws = idx && Array.isArray(idx.weeks) ? idx.weeks.filter(w => w && w.date) : [];
  if (!ws.length) return null;
  const sorted = ws.slice().sort((a, b) => a.date < b.date ? -1 : 1);
  if (!date) return sorted[sorted.length - 1];
  const before = sorted.filter(w => w.date <= date);
  return before.length ? before[before.length - 1] : sorted[sorted.length - 1];
}

G.RankModel = {
  VERSION: 3, TAIL_SPAN, POW, TR, MARK, CONF_TEXT,
  erfc, upper, normCdf, normInv, zOfRank, tOfZ, zOfT, tOfRank, rankOfT,
  build, moments, valueAtZ, locate, valueForT, classicT, plainT, fillPlayers, fitNormalAll, holdout,
  readWeek, toPublic, stringify, addToIndex, pickWeek
};
if (typeof module !== 'undefined' && module.exports) module.exports = G.RankModel;
})(typeof globalThis !== 'undefined' ? globalThis : this);
