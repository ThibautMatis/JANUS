// Interface de JANUS : chargement de la référence, analyse des .ab1, résultats et chromatogramme.

import { parseABIF } from './abif.js';
import { analyzeRead, mergeCalls, modelSupport } from './decompose.js';
import { fetchGenBank, loadGene, referenceFromRecord, referenceFromFasta, referenceFromGenBank, preferredTranscript, listTranscripts } from './reference.js';
import { complement } from './sequence.js';
import { TraceView, snapshot } from './trace-view.js';
import { CLASSES, classOf, conclusion, geneOf } from './report.js';
import { NAME, releaseLabel } from './version.js';

document.getElementById('release').textContent = releaseLabel();
import { proteinName, translate, THREE } from './protein.js';
import { variantLocation, parseVariantC, variantKey, hgvsNames } from './hgvs.js';
import { directCheck, explainDiscrepancy, explainShift } from './check.js';
import { genomicMapping, gnomadId, gnomadFrequencies, gnomadUrl } from './gnomad.js';

const $ = id => document.getElementById(id);
const state = { ref: null, ngRecord: null, localGenBank: null, results: [], failed: [], merged: [] };

const TYPE_LABEL = { snv: 'SNV', del: 'délétion', dup: 'duplication', ins: 'insertion', delins: 'délétion-insertion' };

function setStatus(id, text, error = false) {
  const el = $(id);
  el.textContent = text;
  el.hidden = false;
  el.classList.toggle('error', error);
}

function option(value, text) {
  const o = document.createElement('option');
  o.value = value;
  o.textContent = text;
  return o;
}

function cell(tr, content, className) {
  const td = document.createElement('td');
  if (content instanceof Node) td.append(content); else td.textContent = content ?? '';
  if (className) td.className = className;
  tr.append(td);
  return td;
}

function tag(text, kind) {
  const s = document.createElement('span');
  s.className = 'tag ' + kind;
  s.textContent = text;
  return s;
}

// ---------------------------------------------------------------------------
// 1. Référence

function describeReference(ref) {
  const lines = [`Référence : ${ref.accession} (${ref.seq.length.toLocaleString('fr-FR')} nt) · ${ref.source}`];
  if (ref.definition) lines.push(ref.definition);
  if (ref.transcript) {
    lines.push(`Transcrit : ${ref.transcript.id} · ${ref.transcript.exons.length} exons` + (ref.cds ? ` · CDS ${ref.cds.cdsStart}–${ref.cds.cdsEnd} du transcrit` : ''));
  } else {
    lines.push('Aucun transcrit sélectionné : nomenclature g. uniquement.');
  }
  if (ref.exonNumbering) lines.push(`Exons : ${ref.exonNumbering}.`);
  for (const w of ref.warnings) lines.push('⚠ ' + w);
  return lines.join('\n');
}

function setReference(ref) {
  state.ref = ref;
  markQuickGene(String(geneOf(ref) || '').toUpperCase());
  setStatus('refStatus', describeReference(ref));
  if (state.results.length) setStatus('readStatus', 'Référence modifiée : relancez l’analyse.');
}

function fillTranscripts(transcripts, selected = '') {
  const select = $('nmSelect');
  select.replaceChildren();
  select.append(option('', '— aucun (nomenclature g.) —'));
  for (const t of transcripts) {
    const label = `${t.id}${t.product ? ' — ' + t.product.slice(0, 80) : ''}${/MANE Select/i.test(t.note) ? ' · MANE Select' : ''}`;
    select.append(option(t.id, label));
  }
  select.value = selected;
  select.disabled = !transcripts.length;
}

/** Affiche une référence chargée depuis NCBI et sa liste de transcrits. */
function showLoaded(ref) {
  state.ngRecord = ref.record;
  state.localGenBank = null;
  fillTranscripts(ref.transcripts, ref.transcript?.id || '');
  setReference(ref);
}

async function loadGeneClick() {
  const gene = $('gene').value.trim();
  if (!gene) return setStatus('refStatus', 'Saisissez un symbole de gène (ex. BRCA2).', true);
  try {
    showLoaded(await loadGene(gene, text => setStatus('refStatus', text)));
    try { localStorage.setItem('janus.gene', gene); } catch { /* facultatif */ }
  } catch (e) {
    setStatus('refStatus', `Chargement impossible : ${e.message}. Vous pouvez saisir une accession NG_ précise ou importer un fichier local.`, true);
  }
}
$('loadGene').addEventListener('click', loadGeneClick);
$('gene').addEventListener('keydown', e => { if (e.key === 'Enter') $('loadGene').click(); });

// Accès rapide : un clic sur un gène fréquent le charge directement.
const QUICK_GENES = ['BRCA1', 'BRCA2', 'PALB2', 'RAD51C', 'RAD51D', 'PTEN', 'CDH1', 'MSH2', 'MSH6', 'MLH1'];
function markQuickGene(gene) {
  for (const b of $('geneChips').querySelectorAll('button')) b.setAttribute('aria-pressed', String(b.textContent === gene));
}
for (const gene of QUICK_GENES) {
  const b = document.createElement('button');
  b.type = 'button';
  b.textContent = gene;
  b.title = `Charger ${gene}`;
  b.setAttribute('aria-pressed', 'false');
  b.addEventListener('click', async () => {
    $('gene').value = gene;
    const chips = [...$('geneChips').querySelectorAll('button')];
    chips.forEach(c => { c.disabled = true; });
    try { await loadGeneClick(); } finally { chips.forEach(c => { c.disabled = false; }); }
  });
  $('geneChips').append(b);
}
try { const last = localStorage.getItem('janus.gene'); if (last) $('gene').value = last; } catch { /* facultatif */ }

$('loadNg').addEventListener('click', async () => {
  const acc = $('ngAccession').value.trim().toUpperCase();
  setStatus('refStatus', `Chargement de ${acc}…`);
  try {
    const record = await fetchGenBank(acc);
    const t = preferredTranscript(record, listTranscripts(record));
    showLoaded(await referenceFromRecord(record, t?.id || null));
  } catch (e) {
    setStatus('refStatus', `Chargement impossible : ${e.message}`, true);
  }
});

// Changer de transcrit reconstruit la référence (le NM_ est téléchargé une seule fois puis mémorisé).
$('nmSelect').addEventListener('change', async () => {
  const nm = $('nmSelect').value || null;
  try {
    if (state.localGenBank) setReference(referenceFromGenBank(state.localGenBank, nm));
    else if (state.ngRecord) {
      setStatus('refStatus', nm ? `Chargement du transcrit ${nm}…` : 'Construction de la référence…');
      setReference(await referenceFromRecord(state.ngRecord, nm));
    }
  } catch (e) {
    setStatus('refStatus', `Erreur : ${e.message}`, true);
  }
});

$('loadRefFile').addEventListener('click', async () => {
  const file = $('refFile').files[0];
  if (!file) return setStatus('refStatus', 'Choisissez un fichier GenBank ou FASTA.', true);
  try {
    const text = await file.text();
    if (/^LOCUS/m.test(text)) {
      const ref = referenceFromGenBank(text);
      state.localGenBank = text;
      state.ngRecord = null;
      fillTranscripts(ref.transcripts);
      setReference(ref);
      if (ref.transcripts.length) setStatus('refStatus', describeReference(ref) + '\nChoisissez un transcrit pour les coordonnées c.');
    } else {
      setReference(referenceFromFasta(text));
    }
  } catch (e) {
    setStatus('refStatus', `Fichier refusé : ${e.message}`, true);
  }
});

// ---------------------------------------------------------------------------
// 2. Analyse

const options = () => ({
  ratio: Number($('optRatio').value) || 0.25,
  minSignal: Number($('optSignal').value) || 0.08,
  snrMin: Number($('optSnr').value) || 4,
});

$('analyze').addEventListener('click', async () => {
  const files = [...$('abFiles').files];
  if (!state.ref) return setStatus('readStatus', 'Chargez d’abord une référence (étape 1).', true);
  if (!files.length) return setStatus('readStatus', 'Choisissez au moins un fichier .ab1.', true);
  setStatus('readStatus', 'Analyse en cours…');
  await new Promise(r => setTimeout(r, 20));

  const results = [], failed = [];
  for (const file of files) {
    try {
      const read = parseABIF(await file.arrayBuffer());
      read.name = file.name;
      results.push(analyzeRead(read, state.ref, options()));
    } catch (e) {
      const reason = e.message.startsWith(file.name) ? e.message.slice(file.name.length).replace(/^\s*:\s*/, '') : e.message;
      failed.push({ name: file.name, reason });
    }
  }
  state.results = results;
  state.failed = failed;
  $('gnomadStatus').hidden = true;
  // Une lecture illisible n'apporte que du bruit : ni appels, ni confirmation.
  state.merged = mergeCalls(results.filter(r => r.quality.verdict !== 'illisible'), state.ref);
  // Contrôle indépendant : bases du séquenceur alignées directement sur le génome.
  for (const r of results) r.check = r.quality.verdict === 'illisible' ? null : directCheck(r, state.ref);
  for (const v of state.merged) { v.protein = proteinName(v, state.ref); v.location = variantLocation(v, state.ref.mapper); }
  renderReads();
  renderVariants();
  renderTraceSelect();
  evaluateSought();
  resetReport();
  setStatus('readStatus', readSummary(), !results.length);
  // Variant recherché absent : chromatogramme positionné sur sa position, sur la première lecture informative.
  const s = state.sought;
  if (s && !s.found) {
    const first = s.reads.find(x => x.center !== null && x.gain !== null) || s.reads.find(x => x.center !== null);
    if (first) { $('traceRead').value = String(state.results.indexOf(first.result)); showTrace(state.results.indexOf(first.result), first.center); }
  }
});

// ---------------------------------------------------------------------------
// Variant recherché : trouvé ou non, et ce que dit chaque lecture à sa position.

/** Base de la lecture la plus proche d'une position génomique (null si la lecture ne la couvre pas). */
function centerFor(result, g) {
  let center = null, best = Infinity;
  for (const [k, a] of annotationsOf(result).map) if (a.g && Math.abs(a.g - g) < best) { best = Math.abs(a.g - g); center = k; }
  return best <= 2 ? center : null;
}

function evaluateSought() {
  const text = $('sought').value.trim();
  state.sought = null;
  if (!text) { $('soughtStatus').hidden = true; renderVariants(); return; }
  if (!state.ref) return renderSoughtError('chargez d’abord une référence (étape 1)');
  let parsed;
  try { parsed = parseVariantC(text, state.ref); } catch (e) { return renderSoughtError(e.message); }
  const v = parsed.variant;
  Object.assign(v, { key: variantKey(v), names: hgvsNames(v, state.ref), protein: proteinName(v, state.ref), location: variantLocation(v, state.ref.mapper) });
  const s = { variant: v, transcript: parsed.transcript, found: null, near: null, reads: [] };
  if (state.results.length) {
    s.found = state.merged.find(m => m.key === v.key) || null;
    // Indel de même bilan à quelques bases : probablement le même événement, nommé autrement.
    const net = x => x.ins.length - x.del.length;
    if (!s.found && v.type !== 'snv') s.near = state.merged.find(m => m.type !== 'snv' && net(m) === net(v) && Math.abs(m.pos - v.pos) <= 10) || null;
    for (const r of state.results.filter(x => x.quality.verdict !== 'illisible')) {
      const called = s.found?.reads.some(x => x.result === r) ? 'direct' : s.found?.supportedBy.some(x => x.name === r.name) ? 'modèle' : null;
      s.reads.push({ result: r, called, gain: modelSupport(r, v, state.ref.seq), center: centerFor(r, v.pos + 1) });
    }
  }
  state.sought = s;
  renderSought();
  renderVariants();
}

function renderSoughtError(message) {
  const box = $('soughtStatus');
  box.hidden = false;
  box.className = 'sought-box error';
  box.replaceChildren(el('b', {}, 'Variant recherché non valide : '), message);
}

/** Verdict d'une lecture sur le variant recherché. */
function soughtReadVerdict(s, x) {
  if (x.called === 'direct') return 'variant appelé';
  if (x.called === 'modèle') return 'variant confirmé par modélisation (lecture débutant dans la zone décalée)';
  if (x.gain === null) return x.center === null ? 'position non couverte' : 'non informative à cette position (trop près d’une extrémité)';
  const strong = s.variant.type === 'snv' ? 2 : 20;
  return x.gain >= strong ? 'signal compatible avec le variant, mais non appelé : à vérifier sur le chromatogramme' : 'pas de signal du variant à cette position';
}

function renderSought() {
  const s = state.sought, box = $('soughtStatus');
  box.hidden = false;
  const v = s.variant;
  const name = `${v.names.c}${v.protein?.short ? ` ${v.protein.short}` : ''}`;
  const status = !state.results.length ? tag('sera recherché à l’analyse', 'warn')
    : s.found ? tag(`trouvé — ${s.found.status}`, s.found.status === 'confirmé F/R' ? 'ok' : 'warn')
    : tag('non retrouvé', 'bad');
  box.className = 'sought-box' + (!state.results.length ? '' : s.found ? ' found' : ' absent');
  const lines = [];
  if (s.transcript && state.ref.transcript && s.transcript !== state.ref.transcript.id) {
    lines.push(el('li', {}, `⚠ Saisi sur ${s.transcript}, interprété sur ${state.ref.transcript.id} : vérifier que la numérotation c. est la même.`));
  }
  if (s.found) {
    const fr = s.found.reads.map(x => x.alleleFraction).filter(Number.isFinite).map(f => `${Math.round(100 * f)} %`).join(' / ');
    lines.push(el('li', {}, `${s.found.genotype}${fr ? `, fraction ${fr}` : ''}`));
  }
  if (s.near) lines.push(el('li', {}, `≈ indel de même taille appelé à proximité : ${s.near.names.c || s.near.names.g} — vérifier s’il s’agit du même événement décrit autrement.`));
  for (const x of s.reads) {
    const go = el('button', { type: 'button', class: 'ghost' }, 'Voir');
    go.addEventListener('click', () => {
      const i = state.results.indexOf(x.result);
      $('traceRead').value = String(i);
      showTrace(i, x.center);
      $('step-trace').scrollIntoView({ behavior: 'smooth', block: 'start' });
    });
    lines.push(el('li', {}, `${x.result.name} (${x.result.strand === '+' ? 'Forward' : 'Reverse'}) : ${soughtReadVerdict(s, x)}`, x.center !== null ? go : null));
  }
  box.replaceChildren(
    el('h3', {}, 'Variant recherché : ', el('span', { class: 'mono' }, name), v.location ? el('span', { class: 'note' }, `(${v.location})`) : null, status),
    lines.length ? el('ul', {}, lines) : '');
}

$('sought').addEventListener('change', () => {
  evaluateSought();
  if (!state.results.length) return;
  // Le variant recherché, s'il est trouvé, est coché pour le rapport ; les autres choix sont conservés.
  const f = state.sought?.found;
  if (f && report.choices.has(f.key)) report.choices.get(f.key).checked = true;
  renderReportTable();
  updateConclusion();
  renderReportPreview();
});
$('sought').addEventListener('keydown', e => { if (e.key === 'Enter') $('sought').dispatchEvent(new Event('change')); });

/** Orientation devinée d'après le nom de fichier (lectures non localisées). */
function strandFromName(name) {
  if (/forward|fwd/i.test(name)) return '+';
  if (/reverse|rev/i.test(name)) return '-';
  const m = name.match(/(?<![A-Za-z])([FR])(?![A-Za-z])/i);
  return m ? (m[1].toUpperCase() === 'F' ? '+' : '-') : null;
}

/** Bilan des lectures : combien sont exploitables, par brin, et ce qui manque pour confirmer. */
function readSummary() {
  const usable = state.results.filter(r => r.quality.verdict !== 'illisible');
  const lines = [`${usable.length} lecture(s) exploitable(s) sur ${state.results.length + state.failed.length}, ${state.merged.length} variant(s).`];
  for (const [strand, label] of [['+', 'Forward'], ['-', 'Reverse']]) {
    const ok = usable.some(r => r.strand === strand);
    const bad = [...state.results.filter(r => r.strand === strand && r.quality.verdict === 'illisible').map(r => r.name),
      ...state.failed.filter(f => strandFromName(f.name) === strand).map(f => f.name)];
    if (!ok) lines.push(`⚠ Aucune lecture ${label} exploitable${bad.length ? ` (${bad.join(', ')} illisible)` : ''} : les variants ne pourront pas être confirmés sur les deux brins.`);
  }
  for (const f of state.failed) lines.push(`⚠ ${f.name} : illisible — ${f.reason}`);
  return lines.join('\n');
}

function regionLabel(r) {
  const m = state.ref.mapper;
  const g = `g.${r.refStart}_${r.refEnd}`;
  return m ? `${g} (c.${m.c(r.refStart)}_${m.c(r.refEnd)})` : g;
}

const QUALITY_TAG = { bonne: ['exploitable', 'ok'], médiocre: ['exploitable avec prudence', 'warn'], illisible: ['illisible', 'bad'] };

function renderReads() {
  const table = $('readTable'), body = table.querySelector('tbody');
  body.replaceChildren();
  table.hidden = !state.results.length && !state.failed.length;
  for (const r of state.results) {
    const tr = document.createElement('tr');
    cell(tr, r.name, 'mono');
    cell(tr, r.strand === '+' ? 'Forward (+)' : 'Reverse (−)');
    const [label, kind] = QUALITY_TAG[r.quality.verdict];
    const q = cell(tr, tag(label, kind));
    const details = [...r.quality.reasons, `signal/bruit ${r.quality.snr.toFixed(0)}`, r.quality.meanQ !== null ? `Q moyen ${Math.round(r.quality.meanQ)}` : null].filter(Boolean);
    const note = document.createElement('div');
    note.className = 'note';
    note.textContent = details.join(' · ');
    q.append(note);
    const ck = checkSummary(r);
    const cc = cell(tr, ck ? tag(ck.ok ? 'concordant' : 'à vérifier', ck.ok ? 'ok' : 'warn') : '—');
    if (ck) cc.append(el('div', { class: 'note' }, ck.lines.join(' · ')));
    cell(tr, regionLabel(r), 'mono');
    cell(tr, `${r.alignedBases} / ${r.usableBases}`);
    cell(tr, `${Math.round(100 * r.mixedFraction)} %`);
    const shift = r.partial ? Math.abs(r.partial.ins.length - r.partial.del.length) : 0;
    cell(tr, r.quality.verdict === 'illisible' ? 'ignorée' : String(r.variants.length) + (shift ? ` · décalage d’environ ${shift} nt dès le début de la lecture (taille exacte donnée par l’autre lecture)` : ''));
    body.append(tr);
  }
  for (const f of state.failed) {
    const tr = document.createElement('tr');
    cell(tr, f.name, 'mono');
    const strand = strandFromName(f.name);
    cell(tr, strand ? `${strand === '+' ? 'Forward' : 'Reverse'} ? (d’après le nom)` : '—');
    const q = cell(tr, tag('illisible', 'bad'));
    const note = document.createElement('div');
    note.className = 'note';
    note.textContent = f.reason;
    q.append(note);
    cell(tr, '—'); cell(tr, '—'); cell(tr, '—'); cell(tr, '—'); cell(tr, 'ignorée');
    body.append(tr);
  }
}

function statusTag(v) {
  if (v.status === 'confirmé F/R') return tag('confirmé F/R', 'ok');
  if (v.status === 'une seule lecture') return tag('une seule lecture', 'warn');
  return tag('non retrouvé sur l’autre lecture', 'bad');
}

// Tri du tableau des variants : clic sur un en-tête (une seconde fois pour inverser l'ordre).
const sort = { key: 'position', dir: 1 };
const STATUS_ORDER = { 'confirmé F/R': 0, 'une seule lecture': 1, 'non retrouvé sur l’autre lecture': 2 };
const meanFraction = v => {
  const f = v.reads.map(x => x.alleleFraction).filter(Number.isFinite);
  return f.length ? f.reduce((a, b) => a + b, 0) / f.length : null;
};
// Valeur de tri ; null = inconnue, toujours rangée en fin de tableau.
const SORT_VALUE = {
  position: v => v.pos,
  protein: v => { const m = v.protein?.short?.match(/(\d+)/); return m ? Number(m[1]) : null; },
  type: v => TYPE_LABEL[v.type] || v.type,
  genotype: v => v.genotype,
  fraction: meanFraction,
  status: v => STATUS_ORDER[v.status] ?? 3,
  gnomad: v => (v.gnomad?.freq ? v.gnomad.freq.af : v.gnomad?.freq === null ? 0 : null),
  reads: v => v.reads.length + v.supportedBy.length,
};

function sortedVariants() {
  const value = SORT_VALUE[sort.key];
  return state.merged.slice().sort((a, b) => {
    const x = value(a), y = value(b);
    if (x === null || y === null) return x === y ? a.pos - b.pos : x === null ? 1 : -1;
    const d = typeof x === 'string' ? x.localeCompare(y, 'fr') : x - y;
    return d ? sort.dir * d : a.pos - b.pos;
  });
}

$('variantTable').querySelectorAll('th[data-sort]').forEach(th => th.querySelector('button').addEventListener('click', () => {
  const key = th.dataset.sort;
  sort.dir = sort.key === key ? -sort.dir : 1;
  sort.key = key;
  $('variantTable').querySelectorAll('th[data-sort]').forEach(h => h.removeAttribute('aria-sort'));
  th.setAttribute('aria-sort', sort.dir === 1 ? 'ascending' : 'descending');
  renderVariants();
}));

function renderVariants() {
  const body = $('variantTable').querySelector('tbody');
  body.replaceChildren();
  $('exportCsv').disabled = !state.merged.length;
  $('gnomadButton').disabled = !state.merged.length;
  $('variantDetail').hidden = true;
  if (!state.merged.length) {
    const tr = document.createElement('tr');
    cell(tr, state.results.length ? 'Aucun variant détecté. Un résultat vide n’exclut pas un variant.' : 'Aucune analyse.', 'empty').colSpan = 10;
    body.append(tr);
    return;
  }
  for (const v of sortedVariants()) {
    const tr = document.createElement('tr');
    tr.className = 'clickable' + (state.sought && state.sought.variant.key === v.key ? ' sought' : '');
    cell(tr, v.names.c ? v.names.c.split(':')[1] : '—', 'mono');
    cell(tr, v.location || '—');
    const p = cell(tr, v.protein?.short || '—', 'mono');
    p.title = v.protein?.effect || '';
    cell(tr, v.names.g.split(':')[1], 'mono');
    cell(tr, TYPE_LABEL[v.type] || v.type);
    cell(tr, v.genotype);
    const fractions = v.reads.map(x => x.alleleFraction).filter(Number.isFinite);
    cell(tr, fractions.length ? fractions.map(f => `${Math.round(100 * f)} %`).join(' / ') : '—');
    cell(tr, statusTag(v));
    cell(tr, gnomadCell(v));
    cell(tr, [...v.reads.map(x => `${x.name} (${x.strand})`), ...v.supportedBy.map(x => `${x.name} (${x.strand}, par modélisation)`)].join(', '), 'mono');
    tr.addEventListener('click', () => {
      body.querySelectorAll('tr').forEach(x => x.classList.remove('selected'));
      tr.classList.add('selected');
      showVariant(v);
    });
    body.append(tr);
  }
}

function showVariant(v) {
  const seq = state.ref.seq, flank = 25;
  const left = seq.slice(Math.max(0, v.pos - flank), v.pos);
  const right = seq.slice(v.pos + v.del.length, v.pos + v.del.length + flank);
  const width = Math.max(v.del.length, v.ins.length);
  const pad = s => s + '-'.repeat(width - s.length);
  const lines = [
    v.names.c || v.names.g,
    v.protein?.name ? `${v.protein.name} · ${v.protein.effect}` : v.protein?.effect ? `Conséquence : ${v.protein.effect}` : '',
    v.location ? `Localisation : ${v.location}` : '',
    v.names.c ? v.names.g : '',
    v.gnomad?.id ? `GRCh38 ${v.gnomad.id} · gnomAD v4 : ${formatFrequency(v.gnomad)}` : '',
    `${TYPE_LABEL[v.type] || v.type} · ${v.genotype} · ${v.status}`,
    '',
    `Référence  ${left.toLowerCase()}${pad(v.del)}${right.toLowerCase()}`,
    `Variant    ${left.toLowerCase()}${pad(v.ins)}${right.toLowerCase()}`,
    '',
    ...v.reads.map(x => `${x.name} (${x.strand}) : base ${x.readBase} · ${x.evidence}${Number.isFinite(x.alleleFraction) ? ` · fraction ${Math.round(100 * x.alleleFraction)} %` : ''}${x.flags.length ? ' · ' + x.flags.join(', ') : ''}`),
    ...v.supportedBy.map(x => `${x.name} (${x.strand}) : non appelé directement, mais signal nettement mieux expliqué avec ce variant (gain d’ajustement ${Math.round(x.gain)})`),
    ...(v.missing.length ? [`Non retrouvé sur : ${v.missing.join(', ')} (lecture informative à cette position)`] : []),
  ].filter((l, i) => l !== '' || i > 5);
  const detail = $('variantDetail');
  detail.textContent = lines.join('\n');
  detail.hidden = false;

  const first = v.reads[0];
  const idx = state.results.indexOf(first.result);
  $('traceRead').value = String(idx);
  showTrace(idx, first.readBase - 1);
  $('step-trace').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

$('exportCsv').addEventListener('click', () => {
  const header = ['HGVS_c', 'localisation', 'HGVS_p', 'consequence', 'HGVS_g', 'GRCh38', 'gnomAD_AF', 'gnomAD_AC_AN', 'type', 'genotype', 'statut', 'lectures', 'fractions', 'arguments'];
  const rows = sortedVariants().map(v => [ // dans l’ordre affiché
    v.names.c || '', v.location || '', v.protein?.name || '', v.protein?.effect || '', v.names.g, v.gnomad?.id || '',
    v.gnomad?.freq ? v.gnomad.freq.af : v.gnomad?.freq === null ? 0 : '',
    v.gnomad?.freq ? `${v.gnomad.freq.ac}/${v.gnomad.freq.an}` : '',
    v.type, v.genotype, v.status,
    v.reads.map(x => `${x.name} (${x.strand})`).join('; '),
    v.reads.map(x => Number.isFinite(x.alleleFraction) ? x.alleleFraction.toFixed(2) : '').join('; '),
    v.reads.map(x => x.evidence).join('; '),
  ]);
  const csv = [header, ...rows].map(r => r.map(x => `"${String(x).replace(/"/g, '""')}"`).join(',')).join('\r\n');
  const url = URL.createObjectURL(new Blob(['﻿' + csv], { type: 'text/csv;charset=utf-8' }));
  const a = document.createElement('a');
  a.href = url;
  a.download = 'janus_variants.csv';
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});

// ---------------------------------------------------------------------------
// 4. Chromatogramme

const view = new TraceView($('trace'), { onSelect: k => describeBase(k), label: col => hoverLabel(col) });

/** Étiquette de survol : position c. (sinon g.) en titre, puis g., base de la lecture, génome et allèles. */
function hoverLabel(col) {
  const r = view.result;
  const cOf = g => (state.ref.mapper ? `c.${state.ref.mapper.c(g)}` : null);
  if (col.k === null) {
    const { g, genome } = col.gap;
    return [cOf(g) || `g.${g}`, [cOf(g) ? `g.${g}` : null, `génome ${genome}`, 'aucun pic en face (base absente de la lecture)'].filter(Boolean).join(' · ')];
  }
  const k = col.k, a = view.annotationMap?.get(k);
  if (!r || !a?.g) return [`base ${k + 1}`, k < r?.trim[0] || k >= r?.trim[1] ? 'hors de la zone analysée' : a?.genome === '-' ? 'base lue absente du génome (insertion)' : 'non alignée sur la référence'];
  const alleles = a.a1 === a.a2 ? `${a.a1}` : `${a.a1}/${a.a2}`;
  // Bases dans le sens de la lecture, comme sous les pics : complémentaires du transcrit sur une Reverse.
  const strand = r.strand === '-' ? ' (brin −)' : '';
  const cd = view.context?.(a.g)?.codon;
  const codon = cd ? `codon ${cd.n} (${cd.aa3}, base ${cd.phase + 1}/3)` : null;
  return [cOf(a.g) || `g.${a.g}`, [cOf(a.g) ? `g.${a.g}` : null, `base ${k + 1}`, `génome ${a.genome}`, `allèles ${alleles}${strand}`, codon].filter(Boolean).join(' · ')];
}

/**
 * Contexte d'une position génomique pour le bandeau du chromatogramme : exon / intron (et, dans un
 * intron, l'exon le plus proche et son côté à l'écran), codon et acide aminé de référence.
 */
function contextFor(result) {
  const ref = state.ref, m = ref.mapper;
  if (!m) return null;
  if (!ref._protein && ref.cds) {
    const spliced = m.exons.map(([s, e]) => ref.seq.slice(s - 1, e)).join('');
    ref._protein = translate(spliced.slice(ref.cds.cdsStart - 1));
  }
  const reverse = result.strand === '-';
  return g => {
    const region = m.region(g);
    let nearest = null;
    if (region.kind === 'intron') {
      // c.N+d : exon en amont (positions génomiques plus petites) ; c.N-d : exon en aval.
      const upstream = /\+\d+$/.test(m.c(g));
      nearest = { n: upstream ? region.n : region.next, left: upstream !== reverse };
    }
    let codon = null;
    const t = m.exonic(g);
    if (t !== null && ref.cds && ref._protein && t >= ref.cds.cdsStart && t <= ref.cds.cdsEnd) {
      const i = t - ref.cds.cdsStart, n = Math.floor(i / 3), aa = ref._protein[n] || 'X';
      codon = { n: n + 1, phase: i % 3, aa, aa3: aa === '*' ? 'Ter' : THREE[aa] || 'Xaa' };
    }
    return { region, nearest, codon };
  };
}

/**
 * Annotations d'une lecture, dans le sens de la lecture, à partir de l'alignement de JANUS sur
 * l'allèle de référence :
 *  - map : base lue → { genome, a1, a2, g } (génome « - » pour une base lue absente du génome) ;
 *  - gaps : base lue → bases du génome sans pic, à afficher juste après elle (ligne Génome continue).
 */
function annotationsOf(result) {
  const map = new Map(), gaps = new Map(), flip = result.strand === '-';
  const show = b => (flip && b !== '-' ? complement(b) : b);
  let pending = [], prev = null;
  for (const c of result.cols) {
    if (c.q === null) { pending.push({ genome: show(c.ref), g: result.winStart + c.t + 1 }); continue; }
    const k = result.profile[c.q].index;
    map.set(k, { genome: show(c.ref), a1: show(c.a1), a2: show(c.a2), g: c.t !== null ? result.winStart + c.t + 1 : null, mixed: c.mixed });
    if (pending.length && prev !== null) {
      // Forward : entre la base lue précédente et celle-ci. Reverse : la lecture parcourt la référence à
      // l'envers, les bases sautées viennent après celle-ci, dans l'ordre inverse.
      if (flip) gaps.set(k, pending.reverse()); else gaps.set(prev, pending);
    }
    pending = [];
    prev = k;
  }
  return { map, gaps };
}

/** Bilan du contrôle direct d'une lecture : discordances et alertes non expliquées par les variants retenus. */
function checkSummary(r) {
  const ck = r.check;
  if (!ck) return null;
  const pos = g => (state.ref.mapper ? `c.${state.ref.mapper.c(g)}` : `g.${g}`);
  const unexplained = ck.discrepancies.filter(d => !explainDiscrepancy(d, r, state.merged));
  const shifts = ck.mapping.ranges.filter(x => !explainShift(x, r, state.merged));
  const lines = [`identité ${(100 * ck.identity).toLocaleString('fr-FR', { maximumFractionDigits: 1 })} % sur ${ck.compared} bases`];
  const alerts = [];
  if (!ck.strandAgrees) alerts.push('orientation différente de celle retenue par JANUS');
  for (const s of ck.skippedInMixed) alerts.push(`JANUS a sauté ${s.bases} base(s) du génome en ${pos(s.gStart)} juste avant une zone de doubles pics : vérifier que l’indel n’est pas nommé à l’envers`);
  for (const s of shifts) alerts.push(`positionnement différent de JANUS sur ${s.n} bases autour de ${pos(s.gStart)} (zone propre)`);
  if (unexplained.length) {
    const list = unexplained.slice(0, 5).map(d => `${pos(d.g)} ${d.type === 'délétion' ? `${d.ref} absent` : d.type === 'insertion' ? `+${d.read}` : `${d.ref}>${d.read}`}`);
    alerts.push(`${unexplained.length} discordance(s) lecture / génome non expliquée(s) par les variants retenus : ${list.join(', ')}${unexplained.length > 5 ? '…' : ''}`);
  }
  if (!alerts.length) lines.push(ck.discrepancies.length ? `${ck.discrepancies.length} discordance(s), toutes expliquées par les variants retenus` : 'aucune discordance');
  return { ok: !alerts.length, lines: [...lines, ...alerts] };
}

function renderTraceSelect() {
  const select = $('traceRead');
  select.replaceChildren();
  state.results.forEach((r, i) => select.append(option(String(i), `${r.name} (${r.strand})`)));
  if (state.results.length) showTrace(0);
}

function showTrace(i, center = null) {
  const r = state.results[i];
  if (!r) return;
  // Repères : variants appelés sur cette lecture et position du variant recherché.
  const markers = r.variants.map(v => v.readBase - 1);
  const soughtAt = state.sought?.reads.find(x => x.result === r)?.center;
  if (soughtAt !== undefined && soughtAt !== null) markers.push(soughtAt);
  if (view.result !== r) {
    view.result = r;
    const ann = annotationsOf(r);
    view.annotationMap = ann.map;
    view.context = contextFor(r);
    view.setRead(r.read, ann.map, markers, ann.gaps);
  } else {
    view.markers = markers;
    view.draw();
  }
  if (center !== null) {
    view.centerOn(center);
    describeBase(center);
  }
}

function describeBase(k) {
  const r = view.result;
  if (!r) return;
  const a = view.annotationMap.get(k);
  let text = `${r.name} · base ${k + 1} : ${r.read.bases[k]}`;
  if (a?.g) {
    text += ` · g.${a.g}`;
    if (state.ref.mapper) text += ` · c.${state.ref.mapper.c(a.g)}`;
    text += ` · génome ${a.genome} · allèle 1 ${a.a1} · allèle 2 ${a.a2}`;
  } else if (k < r.trim[0] || k >= r.trim[1]) {
    text += ' · hors de la zone analysée (extrémité de lecture)';
  }
  setStatus('baseInfo', text);
}

$('traceRead').addEventListener('change', e => showTrace(Number(e.target.value)));
$('zoomIn').addEventListener('click', () => view.zoom(0.6));
$('zoomOut').addEventListener('click', () => view.zoom(1.6));
$('zoomAll').addEventListener('click', () => view.read && view.show(0, view.read.bases.length));

// ---------------------------------------------------------------------------
// Fréquences gnomAD (à la demande : envoie les positions des variants à Ensembl et gnomAD)

function formatFrequency(g) {
  if (g.error) return `indisponible (${g.error})`;
  if (g.freq === undefined) return 'non interrogé';
  if (g.freq === null) return 'absent';
  const { af, ac, an } = g.freq;
  const pct = af === 0 ? '0' : 100 * af >= 1e-4 ? (100 * af).toLocaleString('fr-FR', { maximumSignificantDigits: 2 }) : (100 * af).toExponential(1).replace('.', ',');
  return `${pct} % (${ac.toLocaleString('fr-FR')} / ${an.toLocaleString('fr-FR')})`;
}

function gnomadCell(v) {
  if (!v.gnomad) return '—';
  if (!v.gnomad.id) return v.gnomad.error ? `indisponible` : '—';
  const a = document.createElement('a');
  a.href = gnomadUrl(v.gnomad.id);
  a.target = '_blank';
  a.rel = 'noopener noreferrer';
  a.textContent = formatFrequency(v.gnomad);
  a.title = `${v.gnomad.id} — ouvrir dans gnomAD`;
  a.addEventListener('click', e => e.stopPropagation());
  return a;
}

$('gnomadButton').addEventListener('click', async () => {
  if (!state.merged.length) return;
  const button = $('gnomadButton');
  button.disabled = true;
  setStatus('gnomadStatus', 'Correspondance avec le génome de référence GRCh38 (Ensembl)…');
  try {
    const map = await genomicMapping(state.ref);
    for (const v of state.merged) v.gnomad = { id: gnomadId(v, map) };
    const ids = state.merged.map(v => v.gnomad.id).filter(Boolean);
    setStatus('gnomadStatus', `Interrogation de gnomAD v4 pour ${ids.length} variant(s)…`);
    try {
      const freqs = await gnomadFrequencies(ids);
      for (const v of state.merged) if (v.gnomad.id) v.gnomad.freq = freqs.get(v.gnomad.id);
      setStatus('gnomadStatus', `gnomAD v4 (exomes + génomes), chromosome ${map.chrom} (GRCh38), référence sur le brin ${map.strand === "+" ? "+" : "−"}. « absent » = variant non observé dans gnomAD. Cliquer sur une fréquence ouvre la page gnomAD.`);
    } catch (e) {
      // API gnomAD inaccessible : on garde les liens vers le site, qui reste consultable.
      for (const v of state.merged) if (v.gnomad.id) v.gnomad.error = 'API gnomAD inaccessible';
      setStatus('gnomadStatus', `L’API gnomAD n’a pas répondu (${e.message}). Les liens de la colonne gnomAD ouvrent quand même chaque variant sur le site.`, true);
    }
  } catch (e) {
    for (const v of state.merged) v.gnomad = { id: null, error: e.message };
    setStatus('gnomadStatus', `Correspondance GRCh38 impossible : ${e.message}`, true);
  }
  renderVariants();
  renderReportPreview();
  button.disabled = false;
});

// ---------------------------------------------------------------------------
// 5. Rapport

const report = { choices: new Map(), edited: false, images: new Map() };

function el(tag, props = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (k === 'class') node.className = v; else if (k in node) node[k] = v; else node.setAttribute(k, v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined && c !== '') node.append(c instanceof Node ? c : String(c));
  return node;
}

/** Nouvelle analyse : variants confirmés F/R cochés par défaut, à classer. */
function resetReport() {
  report.choices = new Map(state.merged.map(v => [v.key, { checked: v.status === 'confirmé F/R' || state.sought?.found === v, classId: '' }]));
  report.edited = false;
  report.images = new Map();
  renderReportTable();
  updateConclusion();
  renderReportPreview();
}

const selectedForReport = () => state.merged
  .filter(v => report.choices.get(v.key)?.checked)
  .map(v => ({ variant: v, classId: report.choices.get(v.key).classId }));

function renderReportTable() {
  const body = $('reportTable').querySelector('tbody');
  body.replaceChildren();
  if (!state.merged.length) {
    body.append(el('tr', {}, el('td', { colSpan: 4, class: 'empty' }, state.results.length ? 'Aucun variant détecté : le rapport conclura à l’absence de variant dans la région analysée.' : 'Lancez d’abord une analyse.')));
    return;
  }
  for (const v of state.merged) {
    const choice = report.choices.get(v.key);
    const box = el('input', { type: 'checkbox', checked: choice.checked });
    box.setAttribute('aria-label', `Inclure ${v.names.c || v.names.g} dans le rapport`);
    const select = el('select', {}, el('option', { value: '' }, '— à classer —'), CLASSES.map(c => el('option', { value: c.id }, c.long)));
    select.value = choice.classId;
    const refresh = () => {
      choice.checked = box.checked;
      choice.classId = select.value;
      select.classList.toggle('missing', choice.checked && !choice.classId);
      updateConclusion();
      renderReportPreview();
    };
    box.addEventListener('change', refresh);
    select.addEventListener('change', refresh);
    select.classList.toggle('missing', choice.checked && !choice.classId);
    body.append(el('tr', {},
      el('td', {}, box),
      el('td', { class: 'mono' }, `${v.names.c ? v.names.c.split(':')[1] : v.names.g.split(':')[1]} ${v.protein?.short || ''}`),
      el('td', {}, statusTag(v)),
      el('td', {}, select)));
  }
}

function updateConclusion(force = false) {
  if (report.edited && !force) return;
  const absent = state.sought && !state.sought.found ? [state.sought.variant] : [];
  $('repConclusion').value = state.ref && state.results.length ? conclusion(selectedForReport(), state.ref, state.results, absent) : '';
  report.edited = false;
}

$('repConclusion').addEventListener('input', () => { report.edited = true; renderReportPreview(); });
$('repRegenerate').addEventListener('click', () => { updateConclusion(true); renderReportPreview(); });
for (const id of ['repSample', 'repAuthor']) $(id).addEventListener('input', () => renderReportPreview());

/**
 * Bandeau au-dessus d'une capture du rapport : sens de la lecture par rapport au gène. Une lecture
 * Reverse est affichée dans son propre sens : bases complémentaires, gène orienté de droite à gauche.
 */
function directionBadge(result) {
  const forward = result.strand === '+';
  return el('div', { class: `direction ${forward ? 'forward' : 'reverse'}` },
    el('span', { class: 'badge' }, forward ? 'FORWARD' : 'REVERSE'),
    el('span', { class: 'arrow', 'aria-hidden': 'true' }, forward ? '5′ ━━━━━━━━━▶ 3′' : '3′ ◀━━━━━━━━━ 5′'),
    el('span', { class: 'text' }, forward
      ? 'lu dans le sens du gène (de gauche à droite)'
      : 'brin complémentaire : bases complémentaires du gène, qui va de droite à gauche'));
}

/**
 * Bases de la référence génomique concernées par un variant normalisé (1-based, bornes incluses) :
 * bases supprimées ou substituées, bases dupliquées, ou les deux bases qui encadrent une insertion.
 */
function alteredSpan(v) {
  if (v.del.length) return { from: v.pos + 1, to: v.pos + v.del.length };
  if (v.type === 'dup') return { from: v.pos - v.ins.length + 1, to: v.pos };
  return { from: v.pos, to: v.pos + 1 };
}

/** Capture centrée sur le variant pour chaque lecture qui le soutient (appel direct ou modélisation). */
function variantFigures(v) {
  const sources = [
    ...v.reads.map(x => ({ result: x.result, center: x.readBase - 1, how: 'appel direct' })),
    ...v.supportedBy.map(x => {
      const result = state.results.find(r => r.name === x.name);
      if (!result) return null;
      // Base de la lecture la plus proche de la position du variant sur la référence.
      let center = null, best = Infinity;
      for (const [k, a] of annotationsOf(result).map) if (a.g && Math.abs(a.g - (v.pos + 1)) < best) { best = Math.abs(a.g - (v.pos + 1)); center = k; }
      return center === null ? null : { result, center, how: 'confirmé par modélisation' };
    }),
  ].filter(Boolean);
  return sources.map(({ result, center, how }) => {
    const key = `${v.key}|${result.name}`;
    if (!report.images.has(key)) {
      const ann = annotationsOf(result);
      report.images.set(key, snapshot(result.read, ann.map, [center], center, { gaps: ann.gaps, context: contextFor(result), box: alteredSpan(v) }));
    }
    const pos = state.ref.mapper ? `c.${state.ref.mapper.c(v.pos + 1)}` : `g.${v.pos + 1}`;
    return el('figure', {},
      directionBadge(result),
      el('img', { src: report.images.get(key), alt: `Chromatogramme ${result.name} autour de ${pos}` }),
      el('figcaption', {}, `${result.name} (${result.strand === '+' ? 'Forward' : 'Reverse'}, ${how}) — base ${center + 1} surlignée, bases de référence concernées encadrées en rouge ; lignes : lecture, génome, allèle 1, allèle 2 (sens de lecture).`));
  });
}

/**
 * Rapport : variant recherché. S'il n'est pas retrouvé, capture de chaque lecture informative centrée
 * sur sa position, avec ce que dit la lecture (absence de signal, ou signal à vérifier).
 */
function soughtSection() {
  const s = state.sought;
  if (!s || !state.results.length) return [];
  const v = s.variant;
  const name = `${v.names.c}${v.protein?.short ? ` ${v.protein.short}` : ''}`;
  const out = [el('h2', {}, 'Variant recherché'),
    el('p', {}, el('span', { class: 'mono' }, name), v.location ? ` (${v.location})` : '', ' — ', s.found ? `trouvé (${s.found.status})` : 'non retrouvé')];
  if (s.found) return out; // ses captures figurent avec les variants retenus
  for (const x of s.reads) {
    if (x.center === null) continue;
    const key = `recherche|${v.key}|${x.result.name}`;
    if (!report.images.has(key)) {
      const ann = annotationsOf(x.result);
      report.images.set(key, snapshot(x.result.read, ann.map, [x.center], x.center, { gaps: ann.gaps, context: contextFor(x.result), box: alteredSpan(v) }));
    }
    out.push(el('figure', {},
      directionBadge(x.result),
      el('img', { src: report.images.get(key), alt: `Chromatogramme ${x.result.name} à la position ${v.names.c}` }),
      el('figcaption', {}, `${x.result.name} (${x.result.strand === '+' ? 'Forward' : 'Reverse'}) — position du variant recherché surlignée, bases de référence concernées encadrées en rouge : ${soughtReadVerdict(s, x)}.`)));
  }
  return out;
}

/** Texte où chaque occurrence du symbole du gène est mise en italique (convention de nomenclature). */
function italicGene(text, gene) {
  if (!gene) return [text];
  const re = new RegExp(`(?<![A-Za-z0-9])${gene.replace(/[.*+?^${}()|[\]\\-]/g, '\\$&')}(?![A-Za-z0-9])`, 'g');
  const out = [];
  let last = 0;
  for (const m of text.matchAll(re)) {
    if (m.index > last) out.push(text.slice(last, m.index));
    out.push(el('i', { class: 'gene' }, m[0]));
    last = m.index + m[0].length;
  }
  if (last < text.length) out.push(text.slice(last));
  return out;
}

function renderReportPreview() {
  const view = $('reportView');
  const ready = Boolean(state.ref && state.results.length);
  $('repPrint').disabled = !ready;
  view.hidden = !ready;
  if (!ready) return;
  const ref = state.ref, gene = geneOf(ref);
  const selected = selectedForReport();
  const date = new Date().toLocaleDateString('fr-FR', { day: '2-digit', month: 'long', year: 'numeric' });

  const variants = selected.length
    ? selected.map(({ variant: v, classId }) => {
      const cls = classOf(classId);
      const fractions = v.reads.map(x => x.alleleFraction).filter(Number.isFinite).map(f => `${Math.round(100 * f)} %`).join(' / ');
      const rows = [
        ['Variant', v.names.c || v.names.g],
        ['Localisation', v.location || '—'],
        ['Protéine (prédite)', v.protein?.name ? `${v.protein.name} — ${v.protein.effect}` : v.protein?.effect || '—'],
        ['Génomique', v.names.g + (v.gnomad?.id ? ` · GRCh38 ${v.gnomad.id}` : '')],
        ['Génotype', `${v.genotype}${fractions ? ` (fraction ${fractions})` : ''}`],
        ['Confirmation', `${v.status} — ${[...v.reads.map(x => `${x.name} (${x.strand})`), ...v.supportedBy.map(x => `${x.name} (${x.strand}, modélisation)`)].join(', ')}`],
        ['gnomAD v4', v.gnomad ? formatFrequency(v.gnomad) : 'non interrogé'],
        ['Classification', cls ? cls.long : 'à préciser'],
      ];
      return el('div', { class: 'variant' },
        el('h3', { class: 'mono' }, v.names.c ? v.names.c.split(':')[1] : v.names.g.split(':')[1], v.protein?.short ? ` ${v.protein.short}` : ''),
        el('table', {}, el('tbody', {}, rows.map(([k, val]) => el('tr', {}, el('th', {}, k), el('td', { class: k === 'Classification' && cls ? `class-${cls.id}` : k === 'Variant' ? 'mono' : '' }, val))))),
        variantFigures(v));
    })
    : [el('p', {}, 'Aucun variant retenu.')];

  view.replaceChildren(
    el('img', { class: 'report-logo', src: 'img/janus-logo.png', alt: 'JANUS' }),
    el('h1', {}, 'Analyse par séquençage Sanger — gène ', el('i', { class: 'gene' }, gene)),
    el('p', { class: 'meta' }, [
      $('repSample').value.trim() ? `Échantillon : ${$('repSample').value.trim()}` : null,
      `Date : ${date}`,
      $('repAuthor').value.trim() ? `Analyse : ${$('repAuthor').value.trim()}` : null,
    ].filter(Boolean).join(' · ')),
    el('h2', {}, 'Conclusion'),
    el('p', { class: 'conclusion' }, italicGene($('repConclusion').value.trim() || '—', gene)),
    ...soughtSection(),
    el('h2', {}, selected.length > 1 ? 'Variants retenus' : 'Variant retenu'),
    ...variants,
    el('p', { class: 'disclaimer' }, `Analyse réalisée avec ${NAME} ${releaseLabel()} (déconvolution d’allèles de chromatogrammes Sanger, méthode inspirée de Tracy, Rausch et al., BMC Genomics 2020), outil de recherche non validé pour le diagnostic. Les nomenclatures p. sont prédites. Tout variant doit être confirmé visuellement sur le chromatogramme et par une méthode validée ; l’absence de variant détecté n’exclut pas un variant non détectable par cette méthode (mosaïque, grande délétion, allèle non amplifié).`),
  );
}

$('repPrint').addEventListener('click', () => {
  const missing = selectedForReport().filter(x => !x.classId);
  if (missing.length && !window.confirm(`${missing.length} variant(s) retenu(s) sans classification. Imprimer quand même ?`)) return;
  renderReportPreview();
  window.print();
});
