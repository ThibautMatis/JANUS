// Interface de SVIB : chargement de la référence, analyse des .ab1, résultats et chromatogramme.

import { parseABIF } from './abif.js';
import { analyzeRead, mergeCalls } from './decompose.js';
import { fetchGenBank, loadGene, referenceFromRecord, referenceFromFasta, referenceFromGenBank, preferredTranscript, listTranscripts } from './reference.js';
import { complement } from './sequence.js';
import { TraceView } from './trace-view.js';

const $ = id => document.getElementById(id);
const state = { ref: null, ngRecord: null, localGenBank: null, results: [], merged: [] };

const TYPE_LABEL = { snv: 'SNV', del: 'délétion', dup: 'duplication', ins: 'insertion', delins: 'délétion-insertion' };

function setStatus(id, text, error = false) {
  const el = $(id);
  el.textContent = text;
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

  const results = [], errors = [];
  for (const file of files) {
    try {
      const read = parseABIF(await file.arrayBuffer());
      read.name = file.name;
      results.push(analyzeRead(read, state.ref, options()));
    } catch (e) {
      errors.push(e.message.startsWith(file.name) ? e.message : `${file.name} : ${e.message}`);
    }
  }
  state.results = results;
  state.merged = mergeCalls(results, state.ref);
  renderReads();
  renderVariants();
  renderTraceSelect();
  const summary = `${results.length} lecture(s) analysée(s), ${state.merged.length} variant(s).`;
  setStatus('readStatus', [summary, ...errors.map(e => '⚠ ' + e)].join('\n'), errors.length > 0 && !results.length);
});

function regionLabel(r) {
  const m = state.ref.mapper;
  const g = `g.${r.refStart}_${r.refEnd}`;
  return m ? `${g} (c.${m.c(r.refStart)}_${m.c(r.refEnd)})` : g;
}

function renderReads() {
  const table = $('readTable'), body = table.querySelector('tbody');
  body.replaceChildren();
  table.hidden = !state.results.length;
  for (const r of state.results) {
    const tr = document.createElement('tr');
    cell(tr, r.name, 'mono');
    cell(tr, r.strand === '+' ? 'Forward (+)' : 'Reverse (−)');
    cell(tr, regionLabel(r), 'mono');
    cell(tr, `${r.alignedBases} / ${r.usableBases}`);
    cell(tr, `${Math.round(100 * r.mixedFraction)} %`);
    const shift = r.partial ? Math.abs(r.partial.ins.length - r.partial.del.length) : 0;
    cell(tr, String(r.variants.length) + (shift ? ` · décalage de ${shift} nt dès le début de la lecture` : ''));
    body.append(tr);
  }
}

function statusTag(v) {
  if (v.status === 'confirmé F/R') return tag('confirmé F/R', 'ok');
  if (v.status === 'une seule lecture') return tag('une seule lecture', 'warn');
  return tag('non retrouvé sur l’autre lecture', 'bad');
}

function renderVariants() {
  const body = $('variantTable').querySelector('tbody');
  body.replaceChildren();
  $('exportCsv').disabled = !state.merged.length;
  $('variantDetail').hidden = true;
  if (!state.merged.length) {
    const tr = document.createElement('tr');
    cell(tr, state.results.length ? 'Aucun variant détecté. Un résultat vide n’exclut pas un variant.' : 'Aucune analyse.', 'empty').colSpan = 7;
    body.append(tr);
    return;
  }
  for (const v of state.merged) {
    const tr = document.createElement('tr');
    tr.className = 'clickable';
    cell(tr, v.names.c ? v.names.c.split(':')[1] : '—', 'mono');
    cell(tr, v.names.g.split(':')[1], 'mono');
    cell(tr, TYPE_LABEL[v.type] || v.type);
    cell(tr, v.genotype);
    const fractions = v.reads.map(x => x.alleleFraction).filter(Number.isFinite);
    cell(tr, fractions.length ? fractions.map(f => `${Math.round(100 * f)} %`).join(' / ') : '—');
    cell(tr, statusTag(v));
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
    v.names.c ? v.names.g : '',
    `${TYPE_LABEL[v.type] || v.type} · ${v.genotype} · ${v.status}`,
    '',
    `Référence  ${left.toLowerCase()}${pad(v.del)}${right.toLowerCase()}`,
    `Variant    ${left.toLowerCase()}${pad(v.ins)}${right.toLowerCase()}`,
    '',
    ...v.reads.map(x => `${x.name} (${x.strand}) : base ${x.readBase} · ${x.evidence}${Number.isFinite(x.alleleFraction) ? ` · fraction ${Math.round(100 * x.alleleFraction)} %` : ''}${x.flags.length ? ' · ' + x.flags.join(', ') : ''}`),
    ...v.supportedBy.map(x => `${x.name} (${x.strand}) : non appelé directement, mais signal nettement mieux expliqué avec ce variant (gain d’ajustement ${Math.round(x.gain)})`),
    ...(v.missing.length ? [`Non retrouvé sur : ${v.missing.join(', ')} (lecture informative à cette position)`] : []),
  ].filter((l, i) => l !== '' || i > 2);
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
  const header = ['HGVS_c', 'HGVS_g', 'type', 'genotype', 'statut', 'lectures', 'fractions', 'arguments'];
  const rows = state.merged.map(v => [
    v.names.c || '', v.names.g, v.type, v.genotype, v.status,
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

const view = new TraceView($('trace'), { onSelect: k => describeBase(k) });

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
