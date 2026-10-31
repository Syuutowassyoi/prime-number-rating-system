#!/usr/bin/env node
/* ============================================================
   ranking-ocr.html の読み取りを、パソコンで確かめるための道具です。
   Node.js だけで動きます（追加のインストールはいりません）。

   使い方（リポジトリのいちばん上で）:
     node tools/ranking-ocr-dev.js test  <画像フォルダ>
       … フォルダの PNG を全部読み取って、結果と checking! を表示します。
         data.json の最新の週とも くらべます。
     node tools/ranking-ocr-dev.js train <画像フォルダ> [種目表.json]
       … 読み取りに使う「お手本」を、そのフォルダの画像から作り直して
         ranking-ocr.html に書き込みます（アプリの見た目が変わった時用）。
         種目表.json は { "IMG_0001.png": "g01", ... } の形で、
         種目名がまだ読めない画像の種目を教える時だけ使います。
   ============================================================ */
'use strict';
const fs = require('fs');
const path = require('path');
const zlib = require('zlib');

const ROOT = path.join(__dirname, '..');
const HTML = path.join(ROOT, 'ranking-ocr.html');

/* ---------- ranking-ocr.html から、読み取りのしくみとお手本を取り出します ---------- */
function loadOCR() {
  const html = fs.readFileSync(HTML, 'utf8');
  const core = /<script id="ocr-core">([\s\S]*?)<\/script>/.exec(html);
  const tpl = /<script type="application\/json" id="ocr-templates">([\s\S]*?)<\/script>/.exec(html);
  if (!core || !tpl) throw new Error('ranking-ocr.html の中に ocr-core / ocr-templates が見つかりません');
  new Function(core[1])();
  const O = globalThis.PrimeOCR;
  const templates = JSON.parse(tpl[1]);
  O.setTemplates(templates);
  return { O, html, templates };
}

/* ---------- PNG を読む（8bit・インターレースなし） ---------- */
function decodePNG(buf) {
  if (buf.readUInt32BE(0) !== 0x89504e47) throw new Error('PNG ではありません');
  let pos = 8, w = 0, h = 0, depth = 0, ctype = 0, inter = 0;
  const idat = [];
  while (pos < buf.length) {
    const len = buf.readUInt32BE(pos), type = buf.toString('latin1', pos + 4, pos + 8);
    const data = buf.subarray(pos + 8, pos + 8 + len);
    if (type === 'IHDR') { w = data.readUInt32BE(0); h = data.readUInt32BE(4); depth = data[8]; ctype = data[9]; inter = data[12]; }
    else if (type === 'IDAT') idat.push(data);
    else if (type === 'IEND') break;
    pos += 12 + len;
  }
  const ch = { 0: 1, 2: 3, 4: 2, 6: 4 }[ctype];
  if (depth !== 8 || inter !== 0 || !ch) throw new Error('この形の PNG には対応していません');
  const raw = zlib.inflateSync(Buffer.concat(idat));
  const stride = w * ch, out = new Uint8ClampedArray(w * h * 4);
  let prev = new Uint8Array(stride), cur = new Uint8Array(stride);
  for (let y = 0; y < h; y++) {
    const base = y * (stride + 1), f = raw[base];
    for (let i = 0; i < stride; i++) {
      const a = i >= ch ? cur[i - ch] : 0, b = prev[i], c = i >= ch ? prev[i - ch] : 0;
      let v = raw[base + 1 + i];
      if (f === 1) v += a;
      else if (f === 2) v += b;
      else if (f === 3) v += (a + b) >> 1;
      else if (f === 4) {
        const p = a + b - c, pa = Math.abs(p - a), pb = Math.abs(p - b), pc = Math.abs(p - c);
        v += (pa <= pb && pa <= pc) ? a : (pb <= pc ? b : c);
      }
      cur[i] = v & 255;
    }
    for (let x = 0; x < w; x++) {
      const o = (y * w + x) * 4, i = x * ch;
      if (ch <= 2) { out[o] = out[o + 1] = out[o + 2] = cur[i]; out[o + 3] = ch === 2 ? cur[i + 1] : 255; }
      else { out[o] = cur[i]; out[o + 1] = cur[i + 1]; out[o + 2] = cur[i + 2]; out[o + 3] = ch === 4 ? cur[i + 3] : 255; }
    }
    const t = prev; prev = cur; cur = t;
  }
  return { width: w, height: h, data: out };
}

function listPNG(dir) {
  return fs.readdirSync(dir).filter(f => /\.png$/i.test(f))
    .sort((a, b) => a.localeCompare(b, 'ja', { numeric: true }))
    .map(f => path.join(dir, f));
}

/* ---------- test ---------- */
function cmdTest(dir) {
  const { O } = loadOCR();
  const files = listPNG(dir);
  const items = [];
  let ms = 0;
  for (const [i, f] of files.entries()) {
    const r = O.analyze(decodePNG(fs.readFileSync(f)));
    ms += r.ms;
    items.push({ fi: i, result: r });
    const cat = r.categoryId ? O.CAT_BY_ID[r.categoryId].name : '（種目がわかりません）';
    const iss = r.entries.filter(e => e.issues.length).length;
    console.log(`${path.basename(f)}  ${cat}${r.categorySure ? '' : '（自信なし）'}  ${r.layout}  行 ${r.entries.length}  切れ ${r.rejected.length}  気になる行 ${iss}  人数 ${r.players || '—'}  ${r.ms}ms`);
  }
  const M = O.mergeResults(items);
  let chk = 0;
  console.log('\n--- 種目ごと ---');
  for (const c of O.CATS) {
    const m = M[c.id];
    if (!m) continue;
    chk += m.checking.length;
    console.log(`${c.name}：${m.rows.length} 行（一覧は 1〜${m.listMax}位）抜け ${m.missing.length}  checking! ${m.checking.length}`);
    for (const s of m.checking) console.log(`   checking! #${s.rank} ${s.text}  ${s.flags.join(',')}`);
  }
  /* data.json の最新の週と、同じ順位どうしをくらべます */
  const dj = path.join(ROOT, 'data.json');
  if (fs.existsSync(dj)) {
    const data = JSON.parse(fs.readFileSync(dj, 'utf8'));
    const wk = data.weeks[data.weeks.length - 1];
    const byName = {};
    for (const id in M) byName[M[id].name] = M[id];
    const val = (name, r) => { const m = byName[name]; const s = m && m.map[r]; return s && !isNaN(s.score) ? s.score : undefined; };
    let same = 0, diff = 0;
    const cmp = (label, name, obj) => {
      for (const [r, v] of Object.entries(obj || {})) {
        const got = val(name, +r);
        if (got === undefined) continue;
        if (Math.abs(got - v) < 1e-9) same++;
        else { diff++; console.log(`   ちがい：${label} ${r}位  data.json=${v}  画像=${got}`); }
      }
    };
    console.log(`\n--- data.json（${wk.date} の週）とくらべます ---`);
    for (const [name, b] of Object.entries(wk.borders || {})) cmp(name, name, b);
    cmp('totalBorders', 'スコアの合計 全種（素数判定を除外）', wk.totalBorders);
    cmp('lsBorders.rand', '生成した乱数の累積', (wk.lsBorders || {}).rand);
    cmp('lsBorders.score', 'スコアの累積', (wk.lsBorders || {}).score);
    console.log(`同じ ${same} ／ ちがう ${diff}（画像に写っている順位だけ くらべています）`);
  }
  console.log(`\n${files.length} 枚、読み取り時間 合計 ${ms}ms、checking! ${chk} 件`);
}

/* ---------- train ---------- */
function cmdTrain(dir, mapFile) {
  const { O, html, templates } = loadOCR();
  const map = mapFile ? JSON.parse(fs.readFileSync(mapFile, 'utf8')) : {};
  const D = O._dev;
  /* 1つの文字につき 1200 個まで（多すぎる時は、まんべんなく選びます） */
  const bags = { score: {}, rank: {}, info: {} }, seen = { score: {}, rank: {}, info: {} };
  const CAP = 1200;
  let seed = 12345;
  const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  const add = (kind, ch) => {
    if (!ch.vec || ch.c.length !== 1 || ch.unsure) return;
    const bag = bags[kind][ch.c] = bags[kind][ch.c] || [];
    const n = seen[kind][ch.c] = (seen[kind][ch.c] || 0) + 1;
    if (bag.length < CAP) bag.push(Array.from(ch.vec));
    else { const k = Math.floor(rnd() * n); if (k < CAP) bag[k] = Array.from(ch.vec); }
  };
  const prints = {};
  for (const f of listPNG(dir)) {
    const img = decodePNG(fs.readFileSync(f));
    const forced = map[path.basename(f)] || null;
    const r = O.analyze(img, { dev: true, categoryId: forced });
    if (r.categoryId) {
      const hdr = D.findHeader(img);
      const p = hdr && D.titlePrint(img, hdr);
      if (p && !p.prefix && (forced || r.categorySure)) (prints[r.categoryId] = prints[r.categoryId] || []).push(p);
    }
    /* checking! の出ない、形の正しい行だけを お手本に使います */
    for (const row of r._rows || []) {
      const e = r.entries.find(x => x.y === Math.round(row.yc) && x.col === row.col);
      if (!e || e.issues.length || !e.scoreOk) continue;
      row.sc.chars.forEach(ch => add('score', ch));
      if (row.rk && row.rk.chars) row.rk.chars.forEach(ch => add('rank', ch));
    }
    if (r._info) r._info.chars.forEach(ch => add('info', ch));
  }
  const dist = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) { const d = a[i] - b[i]; s += d * d; } return s / a.length; };
  const b64 = v => Buffer.from(Uint8Array.from(v, x => Math.max(0, Math.min(255, Math.round(x * 255))))).toString('base64');
  const out = { score: [], rank: [], info: [], titles: [] };
  for (const kind of ['score', 'rank', 'info']) {
    const fresh = new Set(Object.keys(bags[kind]));
    /* 新しい画像に出てこなかった文字は、前のお手本をそのまま残します */
    out[kind] = (templates[kind] || []).filter(t => !fresh.has(t.c));
    for (const c of fresh) {
      const cl = [];
      for (const v of bags[kind][c]) {
        let best = null, bd = 1e9;
        for (const k of cl) { const d = dist(v, k.rep); if (d < bd) { bd = d; best = k; } }
        if (best && bd < 0.009) best.m.push(v); else cl.push({ rep: v, m: [v] });
      }
      for (const k of cl) {
        const mean = new Array(k.rep.length).fill(0);
        k.m.forEach(v => v.forEach((x, j) => { mean[j] += x; }));
        out[kind].push({ c, v: b64(mean.map(x => x / k.m.length)) });
      }
    }
    console.log(`${kind}: ${[...fresh].map(c => `${c}(${seen[kind][c]})`).join(' ')} → お手本 ${out[kind].length} 個`);
  }
  /* 種目名の指紋：種目ごとに いちばん代表的なもの（ほかとの ちがいの合計が いちばん小さいもの）と、
     それから いちばん はなれたもの（見え方のちがいを おぼえるため）を残します */
  const titles = [];
  for (const id of Object.keys(prints)) {
    const ps = prints[id];
    let best = ps[0], bs = Infinity;
    for (const a of ps) { let s2 = 0; for (const b of ps) if (a !== b) s2 += D.titleDistance(a, b); if (s2 < bs) { bs = s2; best = a; } }
    titles.push(best);
    let far = null, fd = 0;
    for (const b of ps) { const dd = D.titleDistance(b, best); if (dd > fd) { fd = dd; far = b; } }
    if (far && fd > 0.02) titles.push(far);
    best.id = id; if (far) far.id = id;
  }
  const fresh = new Set(Object.keys(prints));
  out.titles = (templates.titles || []).filter(t => !fresh.has(t.id))
    .concat(titles.map(p => ({ id: p.id, w: p.w, h: p.h, a: Buffer.from(p.a).toString('base64') })));
  console.log(`種目名の指紋：${fresh.size} 種目を作り直し、合計 ${out.titles.length} 個`);
  const next = html.replace(/(<script type="application\/json" id="ocr-templates">)[\s\S]*?(<\/script>)/,
    (m0, a, b) => a + JSON.stringify(out) + b);
  fs.writeFileSync(HTML, next);
  console.log('ranking-ocr.html に書き込みました。test でもう一度確かめてください。');
}

const [cmd, dir, extra] = process.argv.slice(2);
if (cmd === 'test' && dir) cmdTest(dir);
else if (cmd === 'train' && dir) cmdTrain(dir, extra);
else {
  console.log('使い方：node tools/ranking-ocr-dev.js test <画像フォルダ>');
  console.log('　　　　node tools/ranking-ocr-dev.js train <画像フォルダ> [種目表.json]');
  process.exit(1);
}
