// Interface de SVIB : chargement de la référence, analyse des .ab1, résultats et chromatogramme.

import { parseABIF } from './abif.js';
import { analyzeRead, mergeCalls } from './decompose.js';
import { fetchGenBank, loadGene, referenceFromRecord, referenceFromFasta, referenceFromGenBank, preferredTranscript, listTranscripts } from './reference.js';
import { complement } from './sequence.js';
import { TraceView } from './trace-view.js';
import { proteinName } from './protein.js';
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
  for (const w of ref.warnings) lines.push('⚠ ' + w);
  return lines.join('\n');
}

function setReference(ref) {
  state.ref = ref;
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

$('loadGene').addEventListener('click', async () => {
  const gene = $('gene').value.trim();
  if (!gene) return setStatus('refStatus', 'Saisissez un symbole de gène (ex. BRCA2).', true);
  try {
    showLoaded(await loadGene(gene, text => setStatus('refStatus', text)));
    try { localStorage.setItem('svib.gene', gene); } catch { /* facultatif */ }
  } catch (e) {
    setStatus('refStatus', `Chargement impossible : ${e.message}. Vous pouvez saisir une accession NG_ précise ou importer un fichier local.`, true);
  }
});
$('gene').addEventListener('keydown', e => { if (e.key === 'Enter') $('loadGene').click(); });
try { const last = localStorage.getItem('svib.gene'); if (last) $('gene').value = last; } catch { /* facultatif */ }

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
  for (const v of state.merged) v.protein = proteinName(v, state.ref);
  renderReads();
  renderVariants();
  renderTraceSelect();
  setStatus('readStatus', readSummary(), !results.length);
});

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
    cell(tr, regionLabel(r), 'mono');
    cell(tr, `${r.alignedBases} / ${r.usableBases}`);
    cell(tr, `${Math.round(100 * r.mixedFraction)} %`);
    const shift = r.partial ? Math.abs(r.partial.ins.length - r.partial.del.length) : 0;
    cell(tr, r.quality.verdict === 'illisible' ? 'ignorée' : String(r.variants.length) + (shift ? ` · décalage de ${shift} nt dès le début de la lecture` : ''));
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
    cell(tr, '—'); cell(tr, '—'); cell(tr, '—'); cell(tr, 'ignorée');
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
    cell(tr, state.results.length ? 'Aucun variant détecté. Un résultat vide n’exclut pas un variant.' : 'Aucune analyse.', 'empty').colSpan = 9;
    body.append(tr);
    return;
  }
  for (const v of sortedVariants()) {
    const tr = document.createElement('tr');
    tr.className = 'clickable';
    cell(tr, v.names.c ? v.names.c.split(':')[1] : '—', 'mono');
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
  ].filter((l, i) => l !== '' || i > 4);
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
  const header = ['HGVS_c', 'HGVS_p', 'consequence', 'HGVS_g', 'GRCh38', 'gnomAD_AF', 'gnomAD_AC_AN', 'type', 'genotype', 'statut', 'lectures', 'fractions', 'arguments'];
  const rows = sortedVariants().map(v => [ // dans l’ordre affiché
    v.names.c || '', v.protein?.name || '', v.protein?.effect || '', v.names.g, v.gnomad?.id || '',
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
  a.download = 'svib_variants.csv';
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});

// ---------------------------------------------------------------------------
// 4. Chromatogramme

const view = new TraceView($('trace'), { onSelect: k => describeBase(k), label: k => hoverLabel(k) });

/** Étiquette de survol : position c. (sinon g.) en titre, puis g., base de la lecture et allèles. */
function hoverLabel(k) {
  const r = view.result, a = view.annotationMap?.get(k);
  if (!r || !a?.g) return [`base ${k + 1}`, k < r?.trim[0] || k >= r?.trim[1] ? 'hors de la zone analysée' : 'non alignée sur la référence'];
  const c = state.ref.mapper ? `c.${state.ref.mapper.c(a.g)}` : null;
  const alleles = a.a1 === a.a2 ? `${a.a1}` : `${a.a1}/${a.a2}`;
  // Bases dans le sens de la lecture, comme sous les pics : complémentaires du transcrit sur une Reverse.
  const strand = r.strand === '-' ? ' (brin −)' : '';
  return [c || `g.${a.g}`, [c ? `g.${a.g}` : null, `base ${k + 1}`, `réf. ${a.ref}`, `allèles ${alleles}${strand}`].filter(Boolean).join(' · ')];
}

/** Annotations par base d'origine : référence et allèles, dans le sens de la lecture. */
function annotationsOf(result) {
  const map = new Map(), flip = result.strand === '-';
  const show = b => (flip && b !== '-' ? complement(b) : b);
  for (const c of result.cols) {
    if (c.q === null) continue;
    const index = result.profile[c.q].index;
    map.set(index, { ref: show(c.ref), a1: show(c.a1), a2: show(c.a2), g: c.t !== null ? result.winStart + c.t + 1 : null });
  }
  return map;
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
  if (view.result !== r) {
    view.result = r;
    view.annotationMap = annotationsOf(r);
    view.setRead(r.read, view.annotationMap, r.variants.map(v => v.readBase - 1));
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
    text += ` · réf. ${a.ref} · allèle 1 ${a.a1} · allèle 2 ${a.a2}`;
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
  button.disabled = false;
});
