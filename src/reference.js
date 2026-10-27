// Références : lecture GenBank (texte), construction d'une référence génomique RefSeqGene (NG_)
// avec la structure exonique des transcrits NM_, requêtes NCBI E-utilities.

import { transcriptMapper } from './hgvs.js';

/** Analyse une localisation GenBank (join, complement, order, <, >) en intervalles 1-based. */
export function parseLocation(loc) {
  // Ignore les segments situés dans d'autres entrées (ex. J00194.1:100..202).
  const clean = loc.replace(/\s/g, '').replace(/[A-Za-z_]+\d+(?:\.\d+)?:[<>]?\d+(?:\.\.[<>]?\d+)?/g, '');
  const complement = /^complement\(/.test(clean);
  const intervals = [];
  for (const m of clean.matchAll(/<?(\d+)(?:\.\.>?(\d+))?/g)) {
    const s = Number(m[1]), e = Number(m[2] ?? m[1]);
    intervals.push([Math.min(s, e), Math.max(s, e)]);
  }
  return { intervals, complement };
}

/** Analyse un fichier GenBank texte (un ou plusieurs enregistrements). */
export function parseGenBank(text) {
  const records = [];
  for (const chunk of text.replace(/\r/g, '').split(/^\/\/\s*$/m)) {
    if (!/^LOCUS/m.test(chunk)) continue;
    const accession = chunk.match(/^VERSION\s+(\S+)/m)?.[1] || chunk.match(/^ACCESSION\s+(\S+)/m)?.[1] || '';
    const definition = (chunk.match(/^DEFINITION\s+([\s\S]*?)\n(?=\S)/m)?.[1] || '').replace(/\s+/g, ' ').trim();
    const origin = chunk.split(/^ORIGIN.*$/m)[1] || '';
    const seq = origin.replace(/[^A-Za-z]/g, '').toUpperCase();

    const features = [];
    const ft = chunk.match(/^FEATURES.*\n([\s\S]*?)^(?:ORIGIN|CONTIG|BASE COUNT)/m)?.[1] || '';
    let current = null, lastQualifier = null;
    for (const line of ft.split('\n')) {
      const head = line.match(/^ {5}(\S+)\s+(\S.*)$/);
      if (head) {
        current = { key: head[1], location: head[2].trim(), qualifiers: {} };
        features.push(current);
        lastQualifier = null;
        continue;
      }
      if (!current) continue;
      const body = line.replace(/^ {21}/, '');
      const q = body.match(/^\/(\w+)(?:=(.*))?$/);
      if (q) {
        lastQualifier = q[1];
        const value = (q[2] ?? '').replace(/^"/, '');
        (current.qualifiers[lastQualifier] ??= []).push(value);
      } else if (lastQualifier) {
        const list = current.qualifiers[lastQualifier];
        list[list.length - 1] += (/^[A-Z]+$/.test(body.trim()) ? '' : ' ') + body.trim();
      } else {
        current.location += body.trim(); // localisation sur plusieurs lignes
      }
    }
    for (const f of features) {
      for (const k of Object.keys(f.qualifiers)) f.qualifiers[k] = f.qualifiers[k].map(v => v.replace(/"$/, '').trim());
      Object.assign(f, parseLocation(f.location));
    }
    records.push({ accession, definition, seq, features });
  }
  if (!records.length) throw new Error('Aucun enregistrement GenBank reconnu');
  return records;
}

/** Transcrits (mRNA avec /transcript_id) décrits dans un enregistrement génomique. */
export function listTranscripts(record) {
  const out = [];
  for (const f of record.features) {
    if (f.key !== 'mRNA') continue;
    const id = f.qualifiers.transcript_id?.[0];
    if (!id) continue;
    out.push({
      id,
      exons: f.intervals,
      strand: f.complement ? '-' : '+',
      gene: f.qualifiers.gene?.[0] || '',
      product: f.qualifiers.product?.[0] || '',
      note: (f.qualifiers.note || []).join(' '),
    });
  }
  return out;
}

/** CDS (1-based, coordonnées du transcrit) d'un enregistrement NM_. */
export function cdsOfTranscript(record) {
  const cds = record.features.find(f => f.key === 'CDS');
  if (!cds || cds.complement || !cds.intervals.length) return null;
  return { cdsStart: cds.intervals[0][0], cdsEnd: cds.intervals.at(-1)[1] };
}

/**
 * CDS d'un transcrit déduit des annotations de l'enregistrement génomique (secours si le NM_ n'est
 * pas téléchargeable) : CDS dont toutes les jonctions internes coïncident avec celles des exons.
 */
export function cdsFromGenomic(record, transcript) {
  const mapper = transcriptMapper({ exons: transcript.exons, cdsStart: 1, cdsEnd: 1 });
  const ex = mapper.exons;
  for (const f of record.features) {
    if (f.key !== 'CDS' || f.complement !== (transcript.strand === '-') || !f.intervals.length) continue;
    const iv = f.intervals;
    const start = mapper.exonic(iv[0][0]), end = mapper.exonic(iv.at(-1)[1]);
    if (start === null || end === null) continue;
    let ok = true;
    for (let k = 0; k + 1 < iv.length && ok; k++) {
      const e = ex.findIndex(([a, b]) => b === iv[k][1]);
      ok = e >= 0 && ex[e + 1]?.[0] === iv[k + 1][0];
    }
    if (ok && end - start + 1 === iv.reduce((n, [a, b]) => n + b - a + 1, 0)) return { cdsStart: start, cdsEnd: end, proteinId: f.qualifiers.protein_id?.[0] || '' };
  }
  return null;
}

/**
 * Référence génomique prête pour l'analyse.
 * transcript : { id, exons, strand } (facultatif) ; cds : { cdsStart, cdsEnd } (facultatif) ;
 * nmSeq : séquence NM_ pour vérifier que les exons concaténés lui correspondent.
 */
export function buildReference({ accession, seq, definition = '', transcript = null, cds = null, nmSeq = null, source = '' }) {
  if (!seq || !/^[ACGTN]+$/.test(seq)) throw new Error('Séquence de référence invalide (A/C/G/T/N attendus)');
  const ref = { accession, seq, definition, source, transcript: null, mapper: null, warnings: [] };
  if (!transcript) return ref;
  if (transcript.strand !== '+') {
    ref.warnings.push(`Le transcrit ${transcript.id} est sur le brin − de ${accession} : coordonnées c. non calculées.`);
    return ref;
  }
  const spliced = transcript.exons.map(([s, e]) => seq.slice(s - 1, e)).join('');
  if (nmSeq && spliced !== nmSeq) {
    ref.warnings.push(`Les exons de ${transcript.id} dans ${accession} ne reproduisent pas exactement la séquence NM_ (${spliced.length} vs ${nmSeq.length} nt) : vérifier les coordonnées c.`);
  }
  ref.transcript = transcript;
  if (cds) {
    ref.mapper = transcriptMapper({ exons: transcript.exons, ...cds });
    ref.cds = cds;
    const start = ref.mapper.exons.length && spliced.slice(cds.cdsStart - 1, cds.cdsStart + 2);
    if (start !== 'ATG') ref.warnings.push(`Le CDS de ${transcript.id} ne commence pas par ATG dans la référence : vérifier.`);
  } else {
    ref.warnings.push(`CDS de ${transcript.id} inconnu : seule la nomenclature g. est fournie.`);
  }
  return ref;
}

/** Référence à partir d'un FASTA collé ou importé (sans annotation : nomenclature g. seule). */
export function referenceFromFasta(text) {
  const lines = text.replace(/^﻿/, '').trim().split(/\r?\n/);
  if (!lines[0]?.startsWith('>')) throw new Error('En-tête FASTA absent');
  const accession = lines[0].slice(1).trim().split(/\s+/)[0] || 'référence';
  const seq = lines.slice(1).join('').replace(/\s/g, '').toUpperCase();
  return buildReference({ accession, seq, definition: lines[0].slice(1).trim(), source: 'FASTA manuel' });
}

// ---------------------------------------------------------------------------
// NCBI E-utilities (références publiques uniquement : aucune donnée patient dans les URL).

const EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

async function eutils(endpoint, params) {
  const url = new URL(EUTILS + endpoint);
  for (const [k, v] of Object.entries(params)) url.searchParams.set(k, String(v));
  url.searchParams.set('tool', 'svib');
  const response = await fetch(url, { mode: 'cors', cache: 'no-store', credentials: 'omit', referrerPolicy: 'no-referrer' });
  if (!response.ok) throw new Error('NCBI HTTP ' + response.status);
  const body = await response.text();
  if (body.length > 40e6) throw new Error('Réponse NCBI trop volumineuse');
  if (/^\s*<(!DOCTYPE html|html)/i.test(body)) throw new Error('NCBI a renvoyé une page HTML');
  return body;
}

/** Recherche des enregistrements RefSeqGene (NG_) d'un gène humain. */
export async function searchRefSeqGene(gene) {
  if (!/^[A-Za-z0-9][A-Za-z0-9.-]{0,29}$/.test(gene)) throw new Error('Symbole de gène invalide');
  const term = `${gene}[Gene Name] AND refseqgene[Keyword] AND "Homo sapiens"[Organism]`;
  const search = JSON.parse(await eutils('esearch.fcgi', { db: 'nuccore', term, retmode: 'json', retmax: 20 })).esearchresult;
  const ids = search?.idlist || [];
  if (!ids.length) return [];
  const summary = JSON.parse(await eutils('esummary.fcgi', { db: 'nuccore', id: ids.join(','), retmode: 'json' })).result;
  return (summary.uids || [])
    .map(id => ({ accession: summary[id]?.accessionversion || '', title: summary[id]?.title || '' }))
    .filter(x => /^NG_\d+\.\d+$/.test(x.accession));
}

/** Télécharge un enregistrement GenBank texte et vérifie l'accession retournée. */
export async function fetchGenBank(accession) {
  if (!/^(NG|NM|NC|NR)_\d+\.\d+$/.test(accession)) throw new Error('Accession versionnée attendue (ex. NG_012772.3)');
  const [record] = parseGenBank(await eutils('efetch.fcgi', { db: 'nuccore', id: accession, rettype: 'gbwithparts', retmode: 'text' }));
  if (record.accession !== accession) throw new Error(`NCBI a retourné ${record.accession} au lieu de ${accession}`);
  return record;
}

/** Charge une référence NG_ et, si demandé, un transcrit NM_ (CDS et vérification de séquence). */
export async function loadRefSeqGene(ngAccession, nmAccession = null) {
  const ng = await fetchGenBank(ngAccession);
  const transcripts = listTranscripts(ng);
  let transcript = null, cds = null, nmSeq = null;
  const warnings = [];
  if (nmAccession) {
    transcript = transcripts.find(t => t.id === nmAccession);
    if (!transcript) throw new Error(`${nmAccession} n’est pas décrit dans ${ngAccession} (transcrits : ${transcripts.map(t => t.id).join(', ') || 'aucun'})`);
    try {
      const nm = await fetchGenBank(nmAccession);
      cds = cdsOfTranscript(nm);
      nmSeq = nm.seq;
    } catch (e) {
      cds = cdsFromGenomic(ng, transcript);
      warnings.push(`NM_ non téléchargé (${e.message}) : CDS déduit des annotations de ${ngAccession}.`);
    }
  }
  const ref = buildReference({ accession: ng.accession, seq: ng.seq, definition: ng.definition, transcript, cds, nmSeq, source: 'NCBI RefSeqGene' });
  ref.warnings.unshift(...warnings);
  ref.transcripts = transcripts;
  ref.record = ng;
  return ref;
}

/** Référence à partir d'un fichier GenBank local (NG_ téléchargé au préalable). */
export function referenceFromGenBank(text, transcriptId = null) {
  const [record] = parseGenBank(text);
  const transcripts = listTranscripts(record);
  const transcript = transcriptId ? transcripts.find(t => t.id === transcriptId) : null;
  const cds = transcript ? cdsFromGenomic(record, transcript) : null;
  const ref = buildReference({ accession: record.accession, seq: record.seq, definition: record.definition, transcript, cds, source: 'GenBank local' });
  ref.transcripts = transcripts;
  ref.record = record;
  return ref;
}
