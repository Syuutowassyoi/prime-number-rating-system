/* ============================================================
   素数あつめ ランキングの「偏差値」（概算）と、公開データの形（v4）
   ・RATING（index.html）と ランキング分析（ranking-analysis.html）の両方で使います
   ・画面を持たない、計算だけの部分です（Node.js でも そのまま動きます）

   偏差値の出し方
     1. 上位の割合     p ＝（順位 − 0.5）÷ 参加人数 N
     2. 標準正規分布で「上から p の所」にある値     z ＝ Φ⁻¹(1 − p)
     3. 偏差値 ＝ 50 ＋ 10 × z
   読めた順位のあいだは、記録を y ＝ g(x) に直して（点数 x^(3/4)、タイム −√t、累積 x^0.3）y と z を まっすぐ結びます。

   読めた順位より下（予測分布。v4）
     点数・累積（0 より下が無いもの）：
       記録 x(z) ＝ x_L ×（(z − z_N) ÷ (z_L − z_N)）^κ
         L は 読めた いちばん下の順位、z_N は 最下位（N 位）の z。
         z ＝ z_L で ちょうど x_L になり、最下位の人で はじめて 0 になります（とちゅうで 0 には なりません）。
       κ（下がり方の強さ）：L から z の幅 0.4 の中の点で、ln(x ÷ x_L) を ln((z − z_N) ÷ (z_L − z_N)) の
         原点を通る直線で表した 傾き κ̂（最小二乗）。ただし 点数は 4/3 以上、累積は 1/0.3 以上（κ ＝ max(κ_min, κ̂)）。
         κ ＝ 4/3 は「x^(3/4) が z の直線で、最下位で 0」と同じ形（v3 の決まりと同じ）です。
     タイム：y ＝ −√t を z の直線 y ＝ μ ＋ σz で表します（σ は L から z の幅 0.4 の中の点の最小二乗、直線は L を通ります）。
   予測平均点・予測標準偏差：参加者 N 人の記録を 1 位から N 位まで並べたもの x̂₁ … x̂_N の 平均と標準偏差です
     （読めた順位は その記録、あいだは まっすぐ結んだ値、下は 上の予測分布。足し算は 誤差を補正して 2 回に分けて行います）。

   v3 → v4 で変えた所
     v3 は「x^(3/4)（累積は x^0.3）が z の直線」を そのまま下へ のばしたので、直線が 最下位より ずっと上で 0 に とどき、
     「偏差値 50〜60 に要る記録が 0」になる種目がありました（スコアの累積は 78% の人が 0 点、など）。
     v4 は 最下位の人だけが 0 になるようにして、L での 下がり方は 読めた点に合わせました。
     いろいろな形（v3・x^λ の いろいろな λ・対数正規・対数の2次式・ロジット・…）を、
     「上位 200・300・500 位だけで当てはめて その下を当てる」「累積の 1〜10000 位で 深い所を当てる」
     「一覧より下の記録（非公開）を当てる」で くらべて、いちばん よく当たり、0 にも ならない形を えらびました。
   標準正規分布：Φ⁻¹ は Wichura の AS241（相対誤差 1×10⁻¹⁶ ほど）、erfc は 級数と連分数（相対誤差 1×10⁻¹⁵ ほど）です。
   ============================================================ */
(function (G) {
'use strict';

/* ---------- 標準正規分布（高い精度） ---------- */
const SQRT_PI = Math.sqrt(Math.PI);
/** erf(x)（|x| が小さい所）：erf(x) ＝ 2/√π · e^(−x²) · Σ 2ⁿ x^(2n+1) ÷ (1·3·5·…·(2n+1))（項は ぜんぶ正なので 打ち消しが起きません） */
function erfSeries(x) {
  const x2 = x * x;
  let term = x, sum = x;
  for (let n = 0; n < 500; n++) {
    term *= 2 * x2 / (2 * n + 3);
    sum += term;
    if (term < sum * 1e-17) break;
  }
  return 2 / SQRT_PI * Math.exp(-x2) * sum;
}
/** erfc(x)（x ≥ 1.5）：連分数 erfc(x) ＝ e^(−x²)/√π · 1/(x ＋ (1/2)/(x ＋ 1/(x ＋ (3/2)/(x ＋ …))))（Lentz の方法） */
function erfcCF(x) {
  const tiny = 1e-300;
  let f = x, C = x, Dd = 0;
  for (let n = 1; n < 5000; n++) {
    const a = n / 2;
    Dd = x + a * Dd; if (Math.abs(Dd) < tiny) Dd = tiny;
    C = x + a / C; if (Math.abs(C) < tiny) C = tiny;
    Dd = 1 / Dd;
    const delta = C * Dd;
    f *= delta;
    if (Math.abs(delta - 1) < 1e-16) break;
  }
  return Math.exp(-x * x) / SQRT_PI / f;
}
/** 相補誤差関数 erfc(x) */
function erfc(x) {
  if (x !== x) return NaN;
  if (x < 0) return 2 - erfc(-x);
  if (x < 1.5) return 1 - erfSeries(x);
  if (x > 27) return 0;
  return erfcCF(x);
}
/** 上側確率 Q(z) ＝ 1 − Φ(z)（z より上にいる人の割合） */
const upper = z => z === Infinity ? 0 : z === -Infinity ? 1 : 0.5 * erfc(z / Math.SQRT2);
/** 下側確率 Φ(z) */
const normCdf = z => z === -Infinity ? 0 : z === Infinity ? 1 : 0.5 * erfc(-z / Math.SQRT2);
/** Φ⁻¹(p)：下側確率が p になる z（Wichura, Algorithm AS241 PPND16。相対誤差 約 1×10⁻¹⁶） */
const AS = {
  a: [3.3871328727963666080e0, 1.3314166789178437745e+2, 1.9715909503065514427e+3, 1.3731693765509461125e+4,
      4.5921953931549871457e+4, 6.7265770927008700853e+4, 3.3430575583588128105e+4, 2.5090809287301226727e+3],
  b: [1, 4.2313330701600911252e+1, 6.8718700749205790830e+2, 5.3941960214247511077e+3,
      2.1213794301586595867e+4, 3.9307895800092710610e+4, 2.8729085735721942674e+4, 5.2264952788528545610e+3],
  c: [1.42343711074968357734e0, 4.63033784615654529590e0, 5.76949722146069140550e0, 3.64784832476320460504e0,
      1.27045825245236838258e0, 2.41780725177450611770e-1, 2.27238449892691845833e-2, 7.74545014278341407640e-4],
  d: [1, 2.05319162663775882187e0, 1.67638483018380384940e0, 6.89767334985100004550e-1,
      1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4, 1.05075007164441684324e-9],
  e: [6.65790464350110377720e0, 5.46378491116411436990e0, 1.78482653991729133580e0, 2.96560571828504891230e-1,
      2.65321895265761230930e-2, 1.24266094738807843860e-3, 2.71155556874348757815e-5, 2.01033439929228813265e-7],
  f: [1, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2,
      7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7, 2.04426310338993978564e-15]
};
const poly = (k, r) => { let s = k[7]; for (let i = 6; i >= 0; i--) s = s * r + k[i]; return s; };
function normInv(p) {
  if (!(p > 0)) return -Infinity;
  if (!(p < 1)) return Infinity;
  const q = p - 0.5;
  if (Math.abs(q) <= 0.425) { const r = 0.180625 - q * q; return q * poly(AS.a, r) / poly(AS.b, r); }
  let r = q < 0 ? p : 1 - p;
  r = Math.sqrt(-Math.log(r));
  let x;
  if (r <= 5) { r -= 1.6; x = poly(AS.c, r) / poly(AS.d, r); }
  else { r -= 5; x = poly(AS.e, r) / poly(AS.f, r); }
  return q < 0 ? -x : x;
}
/** 上側確率が q になる z（q が とても小さくても、1 − q の丸めで ずれないようにしています） */
const zOfUpper = q => -normInv(q);
/** 参加 N 人の中の r 位 → z（上半分は上側確率、下半分は下側確率から出すので、どちらの端でも 丸めで ずれません） */
const zOfRank = (r, N) => (r - 0.5 <= N / 2) ? -normInv((r - 0.5) / N) : normInv((N - r + 0.5) / N);
/** z → 偏差値、偏差値 → z */
const tOfZ = z => 50 + 10 * z;
const zOfT = t => (t - 50) / 10;
/** 順位 → 偏差値 */
const tOfRank = (r, N) => tOfZ(zOfRank(r, N));
/** 偏差値 → 順位（小数のまま。ちょうど その偏差値になる順位） */
const rankOfT = (t, N) => N * upper(zOfT(t)) + 0.5;

/* ---------- 記録の見方（y ＝ g(x)。読めた順位のあいだを結ぶ時と、タイムの予測分布に使います） ----------
   pow：点数の 3/4 乗（素数ゲームの点・スコアの合計）。0 点 → 0
   nsqrt：タイムの √（1/2 乗）に マイナスを付けたもの（小さいほど良いので）。単位は「秒」
   pcum：累積の 0.3 乗（生成した乱数の累積・スコアの累積）。0 → 0
   lin・neg・log10：v2 までの形（くらべる時だけ使います） */
const POW = { score: 0.75, time: 0.5, cum: 0.3 };
const powTR = (k, unit, mid) => ({ f: v => Math.pow(v, k), inv: y => y > 0 ? Math.pow(y, 1 / k) : 0, ok: v => v >= 0, floor: true, mid: !!mid, unit, pow: k });
const TR = {
  pow: powTR(POW.score, '点', true),
  nsqrt: { f: v => -Math.pow(v, POW.time), inv: y => y < 0 ? Math.pow(-y, 1 / POW.time) : 0, ok: v => v > 0, unit: '秒', pow: POW.time },
  pcum: powTR(POW.cum, ''),
  lin: { f: v => v, inv: y => y > 0 ? y : 0, ok: v => v >= 0, floor: true, mid: true, unit: '点', pow: 1 },
  neg: { f: v => -v, inv: y => -y, ok: v => v > 0, unit: '秒', pow: 1 },
  log10: { f: v => v > 0 ? Math.log10(v) : -Infinity, inv: y => Math.pow(10, y), ok: v => v >= 0, unit: '桁', log: true }
};
/** いちばん下の所の形を出すのに使う、z の幅 */
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
/** 予測分布の形の名前 */
const TAIL_TEXT = {
  low: '4/3 乗の形（x^(3/4) が z の直線で、最下位で 0。v3 と同じ）',
  heavy: '読めた いちばん下の所の 下がり方に合わせた形',
  cumlow: '1/0.3 乗の形（x^0.3 が z の直線で、最下位で 0）',
  line: 'z の直線（タイム）',
  none: '読めた順位が最下位まで あります'
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
    opt.N で参加人数を変えられます。opt.extra === false で 一覧の外の点を使いません。opt.tr で 記録の見方を変えられます。
    opt.tail で 読めた順位より下の形を変えられます：'v4'（ふだん）／'v3'（前の版）／'log'（対数正規） */
function build(cat, opt) {
  opt = opt || {};
  if (!cat || !Array.isArray(cat.pts) || !cat.pts.length) return null;
  const N = opt.N > 0 ? Math.round(opt.N) : (cat.players > 0 ? Math.round(cat.players) : 0);
  if (!N) return null;
  const asc = !!cat.asc;
  const trName = TR[opt.tr] ? opt.tr : (asc ? 'nsqrt' : (cat.kind === 'cum' || /累積/.test(String(cat.name || '')) ? 'pcum' : 'pow'));
  const tr = TR[trName];
  const tailKind = opt.tail === 'v3' || opt.tail === 'log' ? opt.tail : 'v4';
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
  /* 一覧の外の点（画面の「あなたの順位」の欄、近い日の 深い順位 など）。前後の点と並びが合う時だけ使います */
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
  /* いちばん下の所：下から点を足していき、z の幅が TAIL_SPAN に とどくまでの点 */
  const last = pts[pts.length - 1], tail = [];
  for (let i = pts.length - 1; i >= 0; i--) {
    tail.push(pts[i]);
    if (tail.length >= 2 && pts[i].z - last.z >= TAIL_SPAN) break;
  }
  const zBottom = zOfRank(N, N), D = last.z - zBottom;
  const listMax = list[list.length - 1].r, extraUsed = pts.filter(p => p.src === 'extra').length;
  const m = {
    name: cat.name || '', N, asc, tr: trName, f: tr.f, inv: tr.inv, ok: tr.ok, floor: !!tr.floor, log: !!tr.log, unit: tr.unit, pow: tr.pow,
    pts, last, contig, listMax, extraUsed, tailFrom: tail[tail.length - 1].r, tailN: tail.length,
    zTop: zOfRank(1, N), zBottom, D, tailKind, playersFrom: cat.playersFrom || null,
    borrowFrom: cat.borrowFrom || null, borrowN: cat.borrowN || 0,
    mu: null, sigma: null, b: null, clamp: '', kappa: null, kappaHat: null, kappaMin: null, regime: '', hasTail: false
  };
  if (!(last.r < N) || !(D > 0)) {
    /* 読めた順位が最下位まで ある */
    m.regime = 'none'; m.hasTail = true;
    m.tailX = () => last.v; m.tailZ = () => zBottom;
  } else if (tr.floor && tailKind === 'v4') {
    /* v4：x(z) ＝ x_L ×（(z − z_N) ÷ D）^κ、κ ＝ max(κ_min, κ̂) */
    const xL = last.v, kMin = 1 / tr.pow;
    let sab = 0, saa = 0;
    for (const p of tail) {
      if (p === last || !(p.v > 0) || !(p.z > last.z) || !(xL > 0)) continue;
      const a = Math.log((p.z - zBottom) / D), b = Math.log(p.v / xL);
      sab += a * b; saa += a * a;
    }
    const kHat = saa > 0 ? sab / saa : NaN;
    const kappa = isFinite(kHat) && kHat > kMin ? kHat : kMin;
    m.kappa = kappa; m.kappaHat = isFinite(kHat) ? kHat : null; m.kappaMin = kMin;
    m.regime = kappa > kMin ? 'heavy' : (trName === 'pcum' ? 'cumlow' : 'low');
    m.hasTail = true;
    if (xL > 0) {
      m.tailX = z => z <= zBottom ? 0 : xL * Math.pow((Math.min(z, last.z) - zBottom) / D, kappa);
      m.tailZ = x => x > 0 ? zBottom + D * Math.pow(x / xL, 1 / kappa) : zBottom;
    } else { m.tailX = () => 0; m.tailZ = () => zBottom; }
    /* 参考：L での 接線（y ＝ g(x) の空間で）。σ ＝ dy/dz、μ ＝ y_L − σ z_L */
    if (xL > 0) { m.sigma = tr.pow * kappa * Math.pow(xL, tr.pow) / D; m.mu = last.y - m.sigma * last.z; m.b = m.sigma; }
  } else if (tr.floor && tailKind === 'log') {
    /* くらべる用：ln x を z の直線で表して のばします（0 には なりません） */
    const lt = tail.filter(p => p.v > 0).map(p => ({ z: p.z, y: Math.log(p.v) }));
    let s = slope(lt);
    if (!(s > 0) || !isFinite(s)) s = slope(pts.filter(p => p.v > 0).map(p => ({ z: p.z, y: Math.log(p.v) })));
    if (s > 0 && isFinite(s) && last.v > 0) {
      const lx = Math.log(last.v);
      m.sigma = s; m.mu = lx - s * last.z; m.b = s; m.regime = 'log'; m.hasTail = true;
      m.tailX = z => Math.exp(lx - s * (last.z - Math.max(z, zBottom)));
      m.tailZ = x => x > 0 ? last.z - (lx - Math.log(x)) / s : zBottom;
    }
  } else {
    /* タイム（と v3）：y ＝ μ ＋ σz。σ は いちばん下の所の傾き、μ は 直線が L を通るように決めます */
    let sigma = slope(tail), clamp = '';
    if (!(sigma > 0) || !isFinite(sigma)) sigma = slope(pts);
    if (tailKind === 'v3') {
      /* v3 の決まり：最下位は 0 くらい（σ ≥ y_L ÷（z_L − z_N））、点数は まん中が 0 点より上（σ ≤ y_L ÷ z_L） */
      if (tr.floor && last.y > 0 && last.z > zBottom) {
        const lo = last.y / (last.z - zBottom);
        if (!(sigma >= lo) || !isFinite(sigma)) { sigma = lo; clamp = 'low'; }
      }
      if (tr.mid && last.y > 0 && last.z > 0 && sigma > last.y / last.z) { sigma = last.y / last.z; clamp = 'high'; }
    }
    if (sigma > 0 && isFinite(sigma)) {
      const mu = last.y - last.z * sigma;
      m.sigma = sigma; m.mu = mu; m.b = sigma; m.clamp = clamp; m.regime = tailKind === 'v3' ? 'v3' : 'line'; m.hasTail = true;
      m.tailX = z => { const x = tr.inv(mu + sigma * Math.max(z, zBottom)); return tr.floor && !(x > 0) ? 0 : x; };
      m.tailZ = x => (tr.f(x) - mu) / sigma;
    }
  }
  /* 予測平均点・予測標準偏差は、はじめて使う時に 計算します（参加人数ぶん 足すので） */
  let dist;
  Object.defineProperty(m, 'dist', { enumerable: true, get() { if (dist === undefined) dist = m.hasTail ? moments(m) : null; return dist; } });
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
/** 誤差を補正して足します（Neumaier の方法） */
function nsum(arr) {
  let s = 0, c = 0;
  for (let k = 0; k < arr.length; k++) {
    const x = arr[k], t = s + x;
    c += Math.abs(s) >= Math.abs(x) ? (s - t) + x : (x - t) + s;
    s = t;
  }
  return s + c;
}
/** 参加者 N 人の「予測した記録の並び」x̂₁ … x̂_N の 平均・標準偏差（記録の見方が log10 の時だけ、y ＝ log₁₀ のまま）。
    読めた順位 r ≤ r_L：読めた記録（抜けている順位は まっすぐ結んだ値）
    読めた順位より下 r > r_L：予測分布 x̂_r ＝ x(z_r)、z_r ＝ Φ⁻¹(1 − (r − 0.5) ÷ N)
    平均は 補正つきの和 ÷ N、標準偏差は √(Σ(x − 平均)² ÷ N)（2 回に分けて足すので、けたの大きい累積でも ずれません） */
function moments(m) {
  const N = m.N, P = m.pts, L = m.last;
  const xs = new Float64Array(N);
  let zero = 0, i = 0;
  for (let r = 1; r <= N; r++) {
    let x;
    if (r <= L.r) {
      while (P[i].r < r) i++;
      if (P[i].r === r) x = m.log ? P[i].y : P[i].v;
      else { const y = i === 0 ? P[0].y : yInside(m, zOfRank(r, N)); x = m.log ? y : m.inv(y); }
    } else {
      const v = m.tailX(zOfRank(r, N));
      x = m.log ? m.f(v) : v;
    }
    if (m.floor && !m.log && !(x > 0)) { x = 0; if (r > L.r) zero++; }
    if (!isFinite(x)) x = 0;
    xs[r - 1] = x;
  }
  const mean = nsum(xs) / N;
  const dev = new Float64Array(N);
  for (let k = 0; k < N; k++) { const d = xs[k] - mean; dev[k] = d * d; }
  const sd = Math.sqrt(nsum(dev) / N);
  const worst = xs[N - 1];
  const at = z => valueAtZ(m, z);
  const q = { top25: at(0.6744897501960817), median: at(0), top75: at(-0.6744897501960817) };
  if (m.log) return { mean: Math.pow(10, mean), sd, logMean: mean, median: q.median, unit: '桁', log: true, N, read: L.r, worst: Math.pow(10, worst), q };
  return { mean, sd, median: q.median, q, unit: m.unit, N, read: L.r, worst, zero: zero / N, zeroN: zero };
}
/** z の所の 記録（読めた順位の中は まっすぐ結んだ値、下は 予測分布） */
function valueAtZ(m, z) {
  if (!m) return null;
  const L = m.last;
  if (z >= L.z) { const x = m.inv(yInside(m, z)); return m.floor && !(x > 0) ? 0 : x; }
  if (!m.hasTail) return null;
  return m.tailX(Math.max(z, m.zBottom));
}
/** 予測分布の中での偏差値（読めた順位より下では、偏差値と同じになります） */
function classicT(m, value) {
  if (!m || !m.hasTail) return null;
  const v = Number(value);
  if (!isFinite(v) || !m.ok(v)) return null;
  return tOfZ(m.tailZ(v));
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
  /* 3) 読めた点より下：予測分布で出します（最下位で止めます） */
  const L = m.last;
  if (!m.hasTail) return { rank: null, z: null, T: null, top: null, conf: 'below', mark: MARK.below, N, below: L.r, rankMin: L.r + 1, tMax: tOfZ(L.z) };
  let z = m.tailZ(v);
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
  if (!m.hasTail) return res(null, 'below', { below: m.last.r });
  if (z < m.zBottom) return res(null, 'under');
  let value = m.tailX(z), zero = false;
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
    how：'v4'（この計算）／'v3'（前の版）／'log'（対数正規）／'normal'（全員が正規分布）。
    返す値 { n, mae（平均の大きさ）, bias（平均。＋は 予測が上すぎ）, max, end（いちばん下の点の ずれ） } */
function holdout(cat, K, how) {
  const full = build(cat, { extra: false });
  if (!full) return null;
  const fitPts = cat.pts.filter(p => p[0] <= K);
  const tail = how === 'v3' || how === 'log' ? how : 'v4';
  const m = build(Object.assign({}, cat, { pts: fitPts, extra: [] }), { N: full.N, tail });
  if (!m) return null;
  let nf = null;
  if (how === 'normal') { nf = fitNormalAll(m); if (!nf) return null; }
  let n = 0, sa = 0, sb = 0, mx = 0, end = null;
  for (const p of full.pts) {
    if (p.r <= K) continue;
    const T = tOfRank(p.r, full.N);
    const Tp = how === 'normal' ? tOfZ((m.asc ? nf.mu - p.v : p.v - nf.mu) / nf.sd) : (locate(m, p.v) || {}).T;
    if (Tp == null || !isFinite(Tp)) continue;
    const e = Tp - T; n++; sa += Math.abs(e); sb += e; mx = Math.max(mx, Math.abs(e)); end = e;
  }
  return n ? { n, mae: sa / n, bias: sb / n, max: mx, end } : null;
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

/* ---------- data.json（ボーダー）と 全部の順位（rankings/）を あわせる ---------- */
/** タイムの種目（小さいほど良い） */
const TIME_CATS = new Set(['1から順に', '2から順に', '1桁素数', '2の累乗数']);
const TOTAL_EX = 'スコアの合計 全種（素数判定を除外）', RAND = '生成した乱数の累積', SCORE = 'スコアの累積';
/** 'YYYY-MM-DD' どうしの 日数の差（b − a） */
function dayDiff(a, b) {
  const p = s => { const m = /^(\d{4})-(\d{2})-(\d{2})$/.exec(String(s || '')); return m ? Date.UTC(+m[1], +m[2] - 1, +m[3]) : NaN; };
  return Math.round((p(b) - p(a)) / 86400000);
}
/** data.json の1週（{ date, borders, totalBorders, lsBorders, stats }）を、readWeek と同じ形にします。
    素数ゲームは 1・5・10・50・100・200 位、十種競技総合は 1〜1000 位、累積は 1〜10000 位 などの ボーダーです。
    参加人数は、累積だけ 統計情報（stats）にあります */
function fromDataWeek(dw) {
  if (!dw || typeof dw !== 'object' || !/^\d{4}-\d{2}-\d{2}$/.test(String(dw.date || ''))) return null;
  const cats = {};
  const add = (name, obj, players) => {
    if (!obj || typeof obj !== 'object') return;
    const asc = TIME_CATS.has(name), pts = [];
    for (const k of Object.keys(obj)) {
      const r = Number(k), v = Number(obj[k]);
      if (Number.isInteger(r) && r >= 1 && isFinite(v) && v > 0) pts.push([r, v]);
    }
    if (!pts.length) return;
    pts.sort((a, b) => a[0] - b[0]);
    cats[name] = { name, asc, kind: asc ? 'time' : (/累積/.test(name) ? 'cum' : 'int'), players: players > 0 ? Math.round(players) : null, pts, extra: [] };
  };
  for (const n of Object.keys(dw.borders || {})) add(n, dw.borders[n], null);
  add(TOTAL_EX, dw.totalBorders, null);
  const st = dw.stats || {}, ls = dw.lsBorders || {};
  add(RAND, ls.rand, st.rand && st.rand.players && +st.rand.players.n);
  add(SCORE, ls.score, st.score && st.score.players && +st.score.players.n);
  return { date: dw.date, cats, src: 'data' };
}
/** 全部の順位の週（readWeek の結果）と data.json の週（fromDataWeek の結果）を、日付ごとに あわせます。
    同じ日付・同じ順位は 全部の順位のほうを使います。参加人数は、その日に読めた人数 → 前の週の人数（多いほう）の順です。
    返す週：{ date, cats, hasRank, hasData }。種目には own（その日の点の数）・ownMax（その日の いちばん下の順位）が付きます */
function unify(rankWeeks, dataWeeks) {
  const by = new Map();
  const get = d => { if (!by.has(d)) by.set(d, { date: d, cats: {}, hasRank: false, hasData: false }); return by.get(d); };
  const clone = c => Object.assign({}, c, { pts: c.pts.map(p => [p[0], p[1]]), extra: (c.extra || []).map(p => [p[0], p[1]]) });
  for (const w of rankWeeks || []) {
    if (!w || !w.date || !w.cats) continue;
    const u = get(w.date); u.hasRank = true;
    for (const n of Object.keys(w.cats)) { u.cats[n] = clone(w.cats[n]); u.cats[n].playersSrc = w.cats[n].players > 0 ? 'rank' : null; }
  }
  for (const w of dataWeeks || []) {
    if (!w || !w.date || !w.cats) continue;
    const u = get(w.date); u.hasData = true;
    for (const n of Object.keys(w.cats)) {
      const c = w.cats[n], t = u.cats[n];
      if (!t) { u.cats[n] = clone(c); u.cats[n].playersSrc = c.players > 0 ? 'data' : null; continue; }
      const have = new Set(t.pts.map(p => p[0]));
      for (const p of c.pts) if (!have.has(p[0])) t.pts.push([p[0], p[1]]);
      t.pts.sort((a, b) => a[0] - b[0]);
      if (!(t.players > 0) && c.players > 0) { t.players = c.players; t.playersSrc = 'data'; }
    }
  }
  const weeks = Array.from(by.values()).sort((a, b) => a.date < b.date ? -1 : a.date > b.date ? 1 : 0);
  /* 参加人数：無ければ 前の週から。累積の 統計情報の人数は 古いままのことがあるので、前の週（全部の順位で読めた人数）のほうが多ければ そちらを使います */
  const lastKnown = {};
  for (const u of weeks) {
    for (const n of Object.keys(u.cats)) {
      const c = u.cats[n], prev = lastKnown[n];
      c.own = c.pts.length; c.ownMax = c.pts.length ? c.pts[c.pts.length - 1][0] : 0;
      if (!(c.players > 0)) { if (prev) { c.players = prev.n; c.playersFrom = prev.d; } }
      else if (c.playersSrc === 'data' && prev && prev.n > c.players) { c.players = prev.n; c.playersFrom = prev.d; }
      if (c.players > 0 && !c.playersFrom && (!prev || c.players >= prev.n || c.playersSrc === 'rank')) lastKnown[n] = { n: c.players, d: u.date };
    }
  }
  return weeks;
}
/** 近い日（maxDays 日以内）の週に もっと深い順位まで あれば、その順位を「一覧の外の点」として借ります（偏差値の計算用）。
    借りた種目には borrowFrom（借りた週の日付）・borrowN（借りた点の数）が付きます。週の点そのものは変えません */
function borrowDeep(weeks, maxDays) {
  const lim = maxDays > 0 ? maxDays : 28;
  for (const u of weeks || []) {
    for (const n of Object.keys(u.cats)) {
      const c = u.cats[n], own = c.ownMax || (c.pts.length ? c.pts[c.pts.length - 1][0] : 0);
      let best = null;
      for (const v of weeks) {
        if (v === u) continue;
        const o = v.cats[n]; if (!o || !(o.ownMax > own)) continue;
        const dd = Math.abs(dayDiff(u.date, v.date));
        if (!(dd <= lim)) continue;
        if (!best || dd < best.dd || (dd === best.dd && v.date < best.v.date)) best = { v, o, dd };
      }
      if (!best) continue;
      const add = best.o.pts.filter(p => p[0] > own).slice(0, best.o.own);
      if (!add.length) continue;
      c.extra = (c.extra || []).concat(add.map(p => [p[0], p[1]]));
      c.borrowFrom = best.v.date; c.borrowN = add.length; c.borrowMax = add[add.length - 1][0];
    }
  }
  return weeks;
}

/* ---------- 週ごとの変化（マイルストーン・不動日数・ヒートマップ） ---------- */
/** 細かく見る順位（マイルストーン） */
const MILESTONES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90,
  100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000];
/** その種目・その順位の 週ごとの記録 [{ date, v }]（その日に読めた点だけ。借りた点は使いません） */
function seriesAt(weeks, name, r) {
  const out = [];
  for (const u of weeks || []) {
    const c = u.cats && u.cats[name]; if (!c) continue;
    const own = c.own != null ? c.pts.slice(0, c.own) : c.pts;
    let lo = 0, hi = own.length - 1, hit = null;
    while (lo <= hi) { const mid = (lo + hi) >> 1; const pr = own[mid][0]; if (pr === r) { hit = own[mid]; break; } if (pr < r) lo = mid + 1; else hi = mid - 1; }
    if (hit) out.push({ date: u.date, v: hit[1] });
  }
  return out;
}
/** 不動日数：いまの記録が はじめて見えた日から、最後に見た日までの日数（その間 ずっと同じ記録なら「少なくとも この日数は動いていない」）。
    返す値 { days（不動日数）, since（その記録が はじめて見えた日）, until（最後に見た日）, prevChange（その前に ちがう記録を見た日）,
             maxDays（いちばん長い不動日数）, changes（変わった回数）, n（見た回数） } */
function fudoOf(series) {
  const s = (series || []).filter(p => p && p.date && isFinite(p.v));
  if (!s.length) return null;
  let start = 0, maxDays = 0, changes = 0;
  for (let i = 1; i < s.length; i++) {
    if (s[i].v !== s[i - 1].v) { maxDays = Math.max(maxDays, dayDiff(s[start].date, s[i - 1].date)); start = i; changes++; }
  }
  const last = s[s.length - 1], days = dayDiff(s[start].date, last.date);
  maxDays = Math.max(maxDays, days);
  return { days, since: s[start].date, until: last.date, prevChange: start > 0 ? s[start - 1].date : null, maxDays, changes, n: s.length, value: last.v };
}
/** マイルストーンの表（1 種目）。行：{ r, last（最後の記録）, date, prev, prevDate, d（前回からの差）, rel（前回からの割合）,
    first, firstDate, perWeek（1週あたりの差）, fudo（fudoOf）, n, gapNext（次のマイルストーンとの差） } */
function milestoneRows(weeks, name, ranks) {
  const rs = Array.isArray(ranks) ? ranks : MILESTONES;
  const rows = rs.map(r => {
    const s = seriesAt(weeks, name, r);
    if (!s.length) return { r, n: 0 };
    const L = s[s.length - 1], P = s.length > 1 ? s[s.length - 2] : null, F = s[0];
    const span = dayDiff(F.date, L.date);
    return {
      r, n: s.length, series: s, last: L.v, date: L.date, prev: P ? P.v : null, prevDate: P ? P.date : null,
      d: P ? L.v - P.v : null, rel: P && P.v ? (L.v - P.v) / Math.abs(P.v) : null,
      first: F.v, firstDate: F.date, perWeek: span > 0 ? (L.v - F.v) / span * 7 : null, fudo: fudoOf(s)
    };
  });
  for (let i = 0; i < rows.length; i++) {
    const a = rows[i], b = rows.slice(i + 1).find(x => x.n && x.date === a.date);
    if (a.n && b) a.gapNext = { r: b.r, d: a.last - b.last };
  }
  return rows;
}
/** 不動日数のランキング（種目 × マイルストーン）。opt.names（種目）・opt.ranks（順位）・opt.min（何日以上か） */
function fudoRanking(weeks, opt) {
  opt = opt || {};
  const names = opt.names || Array.from(new Set([].concat(...(weeks || []).map(u => Object.keys(u.cats || {})))));
  const rs = opt.ranks || MILESTONES, out = [];
  for (const n of names) for (const r of rs) {
    const f = fudoOf(seriesAt(weeks, n, r));
    if (!f || f.n < 2 || !(f.days >= (opt.min || 0))) continue;
    out.push(Object.assign({ name: n, r }, f));
  }
  out.sort((a, b) => b.days - a.days || b.n - a.n || a.r - b.r || (a.name < b.name ? -1 : 1));
  return out;
}
/** ヒートマップ用：となりあう2つの日（読めた点がある日）のあいだの、順位ごとの変化。
    返す値 [{ from, to, days, cells: Map(順位 → { a（前）, b（後）, d（差）, rel（割合）, good（良くなったか） }) }] */
function heatPairs(weeks, name, opt) {
  opt = opt || {};
  const ws = (weeks || []).filter(u => u.cats && u.cats[name] && (u.cats[name].own || u.cats[name].pts.length) >= (opt.minPts || 1));
  const out = [];
  for (let i = 1; i < ws.length; i++) {
    const A = ws[i - 1].cats[name], B = ws[i].cats[name];
    const ma = new Map((A.own != null ? A.pts.slice(0, A.own) : A.pts).map(p => [p[0], p[1]]));
    const cells = new Map(), asc = !!B.asc;
    for (const p of (B.own != null ? B.pts.slice(0, B.own) : B.pts)) {
      if (!ma.has(p[0])) continue;
      const a = ma.get(p[0]), b = p[1], d = b - a;
      cells.set(p[0], { a, b, d, rel: a ? d / Math.abs(a) : null, good: asc ? d < 0 : d > 0 });
    }
    if (cells.size) out.push({ from: ws[i - 1].date, to: ws[i].date, days: dayDiff(ws[i - 1].date, ws[i].date), cells });
  }
  return out;
}

G.RankModel = {
  VERSION: 4, TAIL_SPAN, POW, TR, MARK, CONF_TEXT, TAIL_TEXT, MILESTONES, TIME_CATS,
  erfc, upper, normCdf, normInv, zOfUpper, zOfRank, tOfZ, zOfT, tOfRank, rankOfT,
  build, moments, valueAtZ, locate, valueForT, classicT, plainT, fillPlayers, fitNormalAll, holdout,
  readWeek, toPublic, stringify, addToIndex, pickWeek,
  dayDiff, fromDataWeek, unify, borrowDeep, seriesAt, fudoOf, milestoneRows, fudoRanking, heatPairs
};
if (typeof module !== 'undefined' && module.exports) module.exports = G.RankModel;
})(typeof globalThis !== 'undefined' ? globalThis : this);
