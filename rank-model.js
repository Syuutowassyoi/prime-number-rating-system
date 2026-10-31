/* ============================================================
   素数あつめ ランキングの「偏差値」（概算）と、公開データの形
   ・RATING（index.html）と ランキング分析（ranking-analysis.html）の両方で使います
   ・画面を持たない、計算だけの部分です（Node.js でも そのまま動きます）

   偏差値の出し方
     1. 上位の割合     p ＝（順位 − 0.5）÷ 参加人数
     2. 標準正規分布で「上から p の所」にある値     z ＝ Φ⁻¹(1 − p)
     3. 偏差値 ＝ 50 ＋ 10 × z
   読めた順位のあいだは、点数と z をまっすぐ結んで見つもります。
   読めた順位より下は「予測分布」（つりがね形。予測平均点 μ・予測標準偏差 σ）で出します：
     z ＝（点数 − μ）÷ σ、偏差値 ＝ 50 ＋ 10 ×（点数 − μ）÷ σ、順位 ≈ 参加人数 × Q(z)
   μ と σ は、読めた順位のいちばん下の所に 直線「点数 ＝ μ ＋ σz」を当てはめて出します。
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

/* ---------- 点数の見方 ----------
   lin：そのまま（素数ゲームの点・スコアの合計）。予測平均点・予測標準偏差の単位は「点」
   neg：マイナスを付ける（タイム。小さいほど良いので）。単位は「秒」
   log10：常用対数（生成した乱数の累積・スコアの累積。上と下で 何万倍もちがうので）。
          予測平均は「中央値」、予測標準偏差の単位は「桁」（1 桁 ＝ 10 倍） */
const TR = {
  lin: { f: v => v, inv: y => y, ok: v => v >= 0 },
  neg: { f: v => -v, inv: y => -y, ok: v => v > 0 },
  log10: { f: v => v > 0 ? Math.log10(v) : -Infinity, inv: y => Math.pow(10, y), ok: v => v >= 0 }
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
  extrap: '読めた順位より下なので、予測分布（予測平均点と予測標準偏差）から出しています',
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
    opt.N で参加人数を変えられます。opt.extra === false で 一覧の外の点を使いません */
function build(cat, opt) {
  opt = opt || {};
  if (!cat || !Array.isArray(cat.pts) || !cat.pts.length) return null;
  const N = opt.N > 0 ? Math.round(opt.N) : (cat.players > 0 ? Math.round(cat.players) : 0);
  if (!N) return null;
  const asc = !!cat.asc;
  const trName = opt.tr || (asc ? 'neg' : (cat.kind === 'cum' || /累積/.test(String(cat.name || '')) ? 'log10' : 'lin'));
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
  /* 予測分布（つりがね形）：いちばん下の所の点に 直線「y ＝ μ ＋ σz」を当てはめます。
     傾き σ が予測標準偏差、z ＝ 0（まん中の人）の所の μ が予測平均点です（直線は いちばん下の点を通します）。
     点数（lin）は 0点より下が無いので、2つの決まりで σ の はばを決めます：
       ・最下位（N 位）の人は 0点か、それより下になる  → σ ≥ y_L ÷（z_L − z_N）
       ・まん中の人は 0点より下にならない（z_L > 0 の時）→ σ ≤ y_L ÷ z_L */
  let sd = slope(tail), clamp = '';
  if (!(sd > 0) || !isFinite(sd)) sd = slope(pts);
  const zBottom = zOfRank(N, N);
  if (trName === 'lin' && last.y > 0) {
    const lo = last.y / (last.z - zBottom), hi = last.z > 0 ? last.y / last.z : Infinity;
    if (!(sd > 0) || !isFinite(sd)) { sd = lo; clamp = 'low'; }
    else if (sd < lo) { sd = lo; clamp = 'low'; }
    else if (sd > hi) { sd = hi; clamp = 'high'; }
  }
  if (!(sd > 0) || !isFinite(sd)) sd = null;
  /* 「まん中 ≥ 0点」の決まりで止めた時は、予測平均点は ちょうど 0点です（−0 にならないように） */
  const mean = !sd ? null : clamp === 'high' ? 0 : last.y - last.z * sd;
  /* 画面に出す形：点（lin）・秒（neg）・中央値と桁（log10） */
  const dist = !sd ? null : trName === 'lin' ? { mean, sd, unit: '点' }
    : trName === 'neg' ? { mean: -mean, sd, unit: '秒' }
    : { mean: Math.pow(10, mean), sd, unit: '桁', log: true };
  const listMax = list[list.length - 1].r, extraUsed = pts.filter(p => p.src === 'extra').length;
  return {
    name: cat.name || '', N, asc, tr: trName, f: tr.f, inv: tr.inv, ok: tr.ok, pts, last, contig,
    sd, mean, b: sd, dist, clamp, listMax, extraUsed, tailFrom: tail[tail.length - 1].r,
    zTop: zOfRank(1, N), zBottom, zFloor: zBottom
  };
}
/** 予測分布だけで出した偏差値 ＝ 50 ＋ 10 ×（点数 − 予測平均点）÷ 予測標準偏差（タイムは 速いほど上） */
function classicT(m, value) {
  if (!m || !m.sd) return null;
  const v = Number(value);
  if (!isFinite(v) || !m.ok(v)) return null;
  return tOfZ((m.f(v) - m.mean) / m.sd);
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
  /* 3) 読めた点より下：予測分布で出します。z ＝（点数 − 予測平均点）÷ 予測標準偏差（最下位で止めます） */
  const L = m.last;
  if (!m.sd) return { rank: null, z: null, T: null, top: null, conf: 'below', mark: MARK.below, N, below: L.r, rankMin: L.r + 1, tMax: tOfZ(L.z) };
  let z = y === -Infinity ? -Infinity : (y - m.mean) / m.sd;
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
  if (!m.sd) return res(null, 'below', { below: m.last.r });
  if (z < m.zBottom) return res(null, 'under');
  let value = m.inv(m.mean + z * m.sd), zero = false;
  if (m.tr === 'lin' && value <= 0) { value = 0; zero = true; }
  return res(value, 'extrap', { from: m.last.r, zero });
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
  VERSION: 2, TAIL_SPAN, MARK, CONF_TEXT,
  erfc, upper, normCdf, normInv, zOfRank, tOfZ, zOfT, tOfRank, rankOfT,
  build, locate, valueForT, classicT, readWeek, toPublic, stringify, addToIndex, pickWeek
};
if (typeof module !== 'undefined' && module.exports) module.exports = G.RankModel;
})(typeof globalThis !== 'undefined' ? globalThis : this);
