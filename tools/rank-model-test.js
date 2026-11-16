#!/usr/bin/env node
/* ============================================================
   偏差値の計算（rank-model.js）を確かめる道具です（Node.js だけで動きます）

     node tools/rank-model-test.js                 … rankings/ の いちばん新しい週で確かめます
     node tools/rank-model-test.js <週の JSON>     … ランキング集計の JSON（公開データでも）で確かめます

   確かめること
     1. 標準正規分布の計算（Φ と Φ⁻¹）が合っているか
     2. 読めた順位の点数を入れると、その順位になるか
     3. 偏差値 → 記録 → 偏差値 と もどしても、もとの偏差値に とどくか
     4. 記録が良いほど、偏差値が下がらないか
     5. 予測平均点・予測標準偏差が出るか（0 より大きいか）、順位ごとに足した値が 積分で出した値と合うか
        点数・累積は「0 と予測される人」が 最下位の 1 人だけか（v4）、偏差値 20〜90 に要る記録が 0 にならないか
     6. （見るだけ）上位 200・500・1000 位だけで当てはめて その下を当てた時の ずれ（v4 ／ v3 ／ 全員が正規分布）
     7. （見るだけ）一覧の外の点を かくして、一覧だけで当てると どれくらい ずれるか
     8. data.json と あわせた週（unify）・近い日の深い順位を借りる（borrowDeep）・不動日数（fudoOf）が 合っているか
   ============================================================ */
'use strict';
const fs = require('fs'), path = require('path');
const M = require(path.join(__dirname, '..', 'rank-model.js'));

function loadWeek(arg) {
  if (arg) return JSON.parse(fs.readFileSync(arg, 'utf8'));
  const dir = path.join(__dirname, '..', 'rankings');
  const idx = JSON.parse(fs.readFileSync(path.join(dir, 'index.json'), 'utf8'));
  const ent = M.pickWeek(idx);
  if (!ent) throw new Error('rankings/index.json に週がありません');
  return JSON.parse(fs.readFileSync(path.join(dir, ent.file), 'utf8'));
}

let bad = 0;
const ng = (...a) => { bad++; if (bad <= 20) console.log('  NG', ...a); };
const near = (a, b, e) => Math.abs(a - b) <= e;

/* 1. 標準正規分布（高い精度：相対誤差 1e-14 より小さいか） */
[[M.normCdf(1.959964), 0.975, 1e-6], [M.normInv(0.975), 1.959963984540054, 1e-14], [M.normInv(1e-6), -4.753424308822899, 1e-13],
 [M.normInv(1e-10), -6.361340902404056, 1e-13], [M.normInv(0.3), -0.5244005127080407, 1e-15],
 [M.upper(3), 0.0013498980316300946, 1e-17], [M.tOfRank(5000.5, 10000), 50, 1e-12]].forEach(([a, b, e], i) => { if (!near(a, b, e)) ng('normal', i, a, b); });
[[0.5, 0.4795001221869535], [1, 0.15729920705028513], [2, 0.004677734981047266], [3, 2.209049699858544e-05],
 [5, 1.5374597944280349e-12], [10, 2.088487583762545e-45]].forEach(([x, v]) => { if (!(Math.abs(M.erfc(x) - v) <= v * 1e-14)) ng('erfc', x, M.erfc(x), v); });
for (const p of [1e-300, 1e-20, 1e-9, 1e-6, 0.001, 0.02, 0.3, 0.5, 0.7, 0.98, 0.999999]) {
  const x = M.normInv(p), back = M.normCdf(x), phi = Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI);
  if (!(Math.abs(back - p) / phi <= 1e-13 * Math.max(1, Math.abs(x)))) ng('roundtrip', p, back);
}
if (!near(M.zOfRank(1, 80145), -M.zOfRank(80145, 80145), 1e-12)) ng('zOfRank symmetric');

const j = loadWeek(process.argv[2]);
const w = M.readWeek(j);
console.log(`週：${w.date}　種目：${Object.keys(w.cats).length}　（rank-model.js v${M.VERSION}）`);
/* 参加人数が読めていない種目は、前の週の人数を使います（rankings/ に前の週がある時） */
if (!process.argv[2]) {
  const dir = path.join(__dirname, '..', 'rankings');
  const idx = JSON.parse(fs.readFileSync(path.join(dir, 'index.json'), 'utf8'));
  const older = (idx.weeks || []).filter(x => x.date < w.date).map(x => M.readWeek(JSON.parse(fs.readFileSync(path.join(dir, x.file), 'utf8'))));
  const n = M.fillPlayers(older.concat([w]));
  if (n) console.log(`参加人数が読めていない ${n} 種目は、前の週の人数を使います：` + Object.keys(w.cats).filter(k => w.cats[k].playersFrom).map(k => `${k}（${w.cats[k].playersFrom}）`).join('・'));
}
/** 予測平均点を 順位ごとに足すのとは別のやり方（z の積分）で出して、くらべます */
function integralMean(m) {
  const N = m.N, L = m.last;
  let s = 0;
  for (let r = 1; r <= L.r; r++) { const x = M.valueAtZ(m, M.zOfRank(r, N)); s += m.log ? Math.log10(x) : x; }
  /* 読めた順位より下：上側確率 q を (r_L ÷ N, 1) で 細かく区切って 中点で足します */
  const q0 = L.r / N, K = 20000;
  for (let i = 0; i < K; i++) {
    const q = q0 + (1 - q0) * (i + 0.5) / K;
    const z = Math.max(m.zBottom, -M.normInv(q));
    let x = M.valueAtZ(m, z); if (m.floor && !(x > 0)) x = 0;
    s += (m.log ? Math.log10(x) : x) * (1 - q0) * N / K;
  }
  return m.log ? Math.pow(10, s / N) : s / N;
}
const fmt = (c, v) => c.asc ? v.toFixed(3) + '秒' : Math.round(v).toLocaleString('en-US');
const fmtD = d => d.log ? `幾何平均 ${Math.round(d.mean).toLocaleString('en-US')} ／ 標準偏差 ${d.sd.toFixed(2)} 桁`
  : `平均 ${d.mean.toLocaleString('en-US', { minimumFractionDigits: 2, maximumFractionDigits: 2 })} ／ 標準偏差 ${d.sd.toLocaleString('en-US', { minimumFractionDigits: 2, maximumFractionDigits: 2 })} ${d.unit}`;
const val = [], dist = [], hold = [];
for (const name of Object.keys(w.cats)) {
  const c = w.cats[name], m = M.build(c);
  if (!m) { console.log(`  ${name}：参加人数が無いので とばします`); continue; }
  /* 5. 予測分布があるか（予測標準偏差 > 0、予測平均点 > 0）、順位ごとに足した値と 積分の値が合うか */
  const D = m.dist;
  if (!D || !m.hasTail || !(D.sd > 0) || !(D.mean > 0)) ng('no predicted distribution', name, D && D.mean, D && D.sd);
  else {
    /* v4：0 と予測されるのは 最下位の 1 人だけ。偏差値 20〜90 に要る記録は 0 にならない */
    if (m.floor && D.zeroN > 1) ng('zero mass', name, D.zeroN);
    if (m.floor) for (let T = 20; T <= 90; T += 5) { const r = M.valueForT(m, T); if (r && r.conf !== 'over' && !(r.value > 0)) ng('zero value for T', name, T, r && r.value); }
    const im = integralMean(m);
    if (!(Math.abs(im - D.mean) <= Math.max(1e-6, Math.abs(D.mean) * 2e-3))) ng('mean (sum vs integral)', name, D.mean, im);
    dist.push(`  ${name.padEnd(20)} 参加 ${String(m.N).padStart(6)}${m.playersFrom ? '*' : ' '} ${fmtD(D)}` +
      (D.log ? `（中央値 ${Math.round(D.median).toLocaleString('en-US')}）` : `（中央値 ${fmt(c, D.median)}）`) +
      (m.kappa ? `　κ ＝ ${m.kappa.toFixed(3)}（${m.regime === 'heavy' ? '読めた所に合わせた形' : '4/3 乗などの ゆるい形'}）` : ''));
  }
  /* 6. 上位 K 位だけで当てはめて、その下を当てる */
  if (!m.log) {
    const cell = K => { const a = M.holdout(c, K, 'v4'), b = M.holdout(c, K, 'v3'), n = M.holdout(c, K, 'normal'); return a && b ? [a.mae, b.mae, n ? n.mae : null] : null; };
    const h2 = cell(200), h5 = cell(500), h10 = cell(1000);
    if (h2 || h5) hold.push({ name, h2, h5, h10 });
  }
  /* 2. 読めた点 → その順位 */
  for (const p of m.pts) {
    const L = M.locate(m, p.v);
    const tie = m.pts.some(q => q.r < p.r && q.v === p.v);
    if (!L || (L.rank !== p.r && !(tie && L.rank < p.r))) ng('rank', name, p.r, L && L.rank);
  }
  /* 3. 偏差値 → 記録 → 偏差値 */
  for (let T = 30; T <= 95; T += 0.5) {
    const r = M.valueForT(m, T);
    if (!r || r.value == null) continue;
    const L = M.locate(m, r.value);
    if (!L || !(L.T >= T - 1e-6)) ng('T roundtrip', name, T, r.value, L && L.T);
  }
  /* 4. 記録が良いほど 偏差値は下がらない */
  let prev = -Infinity;
  const lo = c.asc ? m.pts[0].v * 0.8 : 0, hi = c.asc ? m.pts[m.pts.length - 1].v * 3 : m.pts[0].v * 1.05;
  for (let k = 0; k <= 400; k++) {
    const v = c.asc ? hi - (hi - lo) * k / 400 : lo + (hi - lo) * k / 400;
    const L = M.locate(m, v);
    if (v > 0 && (!L || L.T == null || L.rank == null)) { ng('no rank', name, v); continue; }
    if (!L || L.T == null) continue;
    if (L.T < prev - 1e-9) ng('monotone', name, v, L.T, prev);
    prev = L.T;
    /* 読めた順位より下は「50 ＋ 10 ×（点 − 予測平均点）÷ 予測標準偏差」と同じになるか（最下位で止めた所は のぞく） */
    if (L.conf === 'extrap' && L.rank < m.N && Math.abs(L.T - M.classicT(m, v)) > 1e-9) ng('classic T', name, v, L.T, M.classicT(m, v));
  }
  /* 5. 一覧の外の点を かくして当てる */
  const m0 = M.build(c, { extra: false });
  for (const p of m.pts.filter(q => q.src === 'extra')) {
    const L = m0 ? M.locate(m0, p.v) : null;
    const T = M.tOfRank(p.r, m.N);
    val.push(`  ${name.padEnd(20)} ${String(p.r).padStart(5)}位 ${fmt(c, p.v).padStart(14)}  本当 ${T.toFixed(1)}  一覧だけ ` +
      (!L ? '—' : L.T == null ? `${L.tMax.toFixed(1)} より下（${T < L.tMax ? '合う' : 'はずれ'}）` : `${L.T.toFixed(1)}${L.mark}（ずれ ${Math.abs(L.T - T) < 0.05 ? '±0.0' : (L.T - T > 0 ? '+' : '') + (L.T - T).toFixed(1)}）`));
  }
}
if (dist.length) console.log('予測平均点・予測標準偏差（参加者 全員の 予測した記録の 平均と標準偏差。* は前の週の人数）：\n' + dist.join('\n'));
if (hold.length) {
  const f = x => x == null ? '   —' : x.toFixed(2).padStart(5);
  const avg = k => { const a = hold.filter(h => h[k]); const mean = j => { const b = a.filter(h => h[k][j] != null); return b.length ? b.reduce((s, h) => s + h[k][j], 0) / b.length : null; }; return a.length ? [mean(0), mean(1), mean(2)] : [null, null, null]; };
  console.log('上位 K 位だけで当てはめて、その下の読めた順位を当てた時の 偏差値のずれ（平均）：v4 ／ v3 ／ 全員が正規分布');
  const k10 = hold.some(h => h.h10);
  const c3 = h => `${f(h && h[0])} ／${f(h && h[1])} ／${f(h && h[2])}`;
  hold.forEach(h => console.log(`  ${h.name.padEnd(20)} K=200：${c3(h.h2)}    K=500：${c3(h.h5)}${k10 ? `    K=1000：${c3(h.h10)}` : ''}`));
  const a2 = avg('h2'), a5 = avg('h5'), a10 = avg('h10');
  console.log(`  ${'（種目の平均）'.padEnd(18)} K=200：${c3(a2)}    K=500：${c3(a5)}${k10 && a10[0] != null ? `    K=1000：${c3(a10)}` : ''}`);
}
if (val.length) console.log('一覧の外の点を かくして当てた結果：\n' + val.join('\n'));
/* 8. data.json と あわせる・深い順位を借りる・不動日数 */
{
  const root = path.join(__dirname, '..');
  const dj = JSON.parse(fs.readFileSync(path.join(root, 'data.json'), 'utf8'));
  const idx = JSON.parse(fs.readFileSync(path.join(root, 'rankings', 'index.json'), 'utf8'));
  const rws = (idx.weeks || []).map(x => M.readWeek(JSON.parse(fs.readFileSync(path.join(root, 'rankings', x.file), 'utf8')), { self: false }));
  const dws = (dj.weeks || []).map(M.fromDataWeek).filter(Boolean);
  const U = M.unify(rws, dws);
  const dates = U.map(u => u.date);
  if (new Set(dates).size !== dates.length) ng('unify: duplicate dates');
  for (const u of U) for (const n of Object.keys(u.cats)) {
    const c = u.cats[n];
    for (let i = 1; i < c.pts.length; i++) if (!(c.pts[i][0] > c.pts[i - 1][0])) { ng('unify: ranks not sorted', u.date, n); break; }
  }
  /* 同じ日付の data.json のボーダーと 全部の順位は 同じ値のはず */
  for (const dw of dj.weeks) {
    const u = U.find(x => x.date === dw.date); if (!u || !u.hasRank) continue;
    for (const n of Object.keys(dw.borders || {})) for (const r of Object.keys(dw.borders[n])) {
      const c = u.cats[n]; const p = c && c.pts.find(q => q[0] === +r);
      if (+dw.borders[n][r] > 0 && (!p || p[1] !== +dw.borders[n][r])) ng('unify: border mismatch', dw.date, n, r, p && p[1], dw.borders[n][r]);
    }
  }
  M.borrowDeep(U, 28);
  const last = U[U.length - 1];
  const lastEv = last.cats['素因数分解・初級'];
  console.log(`あわせた週：${U.length} 週（${U[0].date} 〜 ${last.date}）。いちばん新しい週の 素因数分解・初級：自分の点 ${lastEv.own}・借りた点 ${lastEv.borrowN || 0}（${lastEv.borrowFrom || '—'}）・参加 ${lastEv.players}${lastEv.playersFrom ? `（${lastEv.playersFrom} の人数）` : ''}`);
  for (const n of Object.keys(last.cats)) { const m = M.build(last.cats[n]); if (!m) continue; const D = m.dist; if (m.floor && D.zeroN > 1) ng('latest week zero mass', n, D.zeroN); }
  /* 不動日数 */
  const f1 = M.fudoOf([{ date: '2026-10-01', v: 5 }, { date: '2026-10-08', v: 5 }, { date: '2026-10-15', v: 6 }, { date: '2026-10-22', v: 6 }, { date: '2026-11-05', v: 6 }]);
  if (!f1 || f1.days !== 21 || f1.since !== '2026-10-15' || f1.maxDays !== 21 || f1.changes !== 1 || f1.prevChange !== '2026-10-08') ng('fudoOf', JSON.stringify(f1));
  const rk = M.fudoRanking(U, {});
  if (!rk.length || rk.some((x, i) => i && x.days > rk[i - 1].days)) ng('fudoRanking order');
  console.log(`不動日数ランキング（上位 3）：${rk.slice(0, 3).map(x => `${x.name} ${x.r}位 ${x.days}日（${x.since}〜）`).join('　')}`);
  const hp = M.heatPairs(U, '素因数分解・初級');
  if (!hp.length) ng('heatPairs empty');
}
console.log(bad ? `\n確かめ：${bad} か所 おかしい所があります` : '\n確かめ：ぜんぶ合っています');
process.exit(bad ? 1 : 0);
