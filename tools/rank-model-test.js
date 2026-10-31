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
     5. （見るだけ）一覧の外の点を かくして、一覧だけで当てると どれくらい ずれるか
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

/* 1. 標準正規分布 */
[[M.normCdf(1.959964), 0.975, 1e-6], [M.normInv(0.975), 1.959964, 1e-6], [M.normInv(1e-6), -4.753424, 1e-5],
 [M.upper(3), 0.0013499, 1e-7], [M.tOfRank(5000.5, 10000), 50, 1e-6]].forEach(([a, b, e], i) => { if (!near(a, b, e)) ng('normal', i, a, b); });
for (const p of [1e-9, 1e-6, 0.001, 0.02, 0.3, 0.5, 0.7, 0.98, 0.999999]) {
  if (!near(M.normCdf(M.normInv(p)), p, Math.max(1e-12, p * 2e-7))) ng('roundtrip', p);
}

const j = loadWeek(process.argv[2]);
const w = M.readWeek(j);
console.log(`週：${w.date}　種目：${Object.keys(w.cats).length}`);
const fmt = (c, v) => c.asc ? v.toFixed(3) + '秒' : Math.round(v).toLocaleString('en-US');
const fmtD = d => d.log ? `中央値 ${Math.round(d.mean).toLocaleString('en-US')} ／ σ ${d.sd.toFixed(2)} 桁`
  : `平均 ${d.mean.toLocaleString('en-US', { maximumFractionDigits: 2 })} ／ σ ${d.sd.toLocaleString('en-US', { maximumFractionDigits: 2 })} ${d.unit}`;
const val = [], dist = [];
for (const name of Object.keys(w.cats)) {
  const c = w.cats[name], m = M.build(c);
  if (!m) { console.log(`  ${name}：参加人数が無いので とばします`); continue; }
  /* 予測分布があるか（予測標準偏差 > 0、予測平均点は 点数なら 0 以上） */
  if (!m.dist || !(m.sd > 0)) ng('no predicted distribution', name);
  else {
    if (m.tr === 'lin' && m.mean < -1e-6) ng('negative predicted mean', name, m.mean);
    dist.push(`  ${name.padEnd(20)} 参加 ${String(m.N).padStart(6)}  ${fmtD(m.dist)}${m.clamp ? `（${m.clamp === 'low' ? '最下位≈0点' : 'まん中≥0点'} の決まりで調整）` : ''}`);
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
if (dist.length) console.log('予測分布（予測平均点・予測標準偏差）：\n' + dist.join('\n'));
if (val.length) console.log('一覧の外の点を かくして当てた結果：\n' + val.join('\n'));
console.log(bad ? `\n確かめ：${bad} か所 おかしい所があります` : '\n確かめ：ぜんぶ合っています');
process.exit(bad ? 1 : 0);
