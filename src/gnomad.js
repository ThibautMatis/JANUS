// Fréquences gnomAD (v4, GRCh38) des variants appelés.
// 1. Correspondance NG_ → chromosome : coordonnées du gène (Ensembl REST), séquence de la région,
//    puis repérage de la séquence NG_ dans le chromosome (brin et décalage).
// 2. Conversion du variant au format gnomAD (chrom-pos-ref-alt, normalisé à gauche, base d'ancrage).
// 3. Requête GraphQL à l'API gnomAD.
// Seules des positions de variants et un symbole de gène sont envoyés, jamais les chromatogrammes.

import { revcomp } from './sequence.js';

const ENSEMBL = 'https://rest.ensembl.org';
const GNOMAD = 'https://gnomad.broadinstitute.org/api';
export const GNOMAD_DATASET = 'gnomad_r4';

async function getJSON(url, init) {
  const response = await fetch(url, { credentials: 'omit', referrerPolicy: 'no-referrer', ...init });
  if (!response.ok) throw new Error(`${new URL(url).host} HTTP ${response.status}`);
  return response.json();
}

/** Symbole du gène de la référence (qualifiant /gene du transcrit, ou « (SYMBOLE) » dans la définition). */
export function geneSymbol(ref) {
  return ref.transcript?.gene || ref.definition?.match(/\(([A-Za-z0-9-]+)\)/)?.[1] || null;
}

/**
 * Repère la référence NG_ dans une séquence chromosomique `chrSeq` commençant à la position `chrStart`
 * (1-based). Retourne { strand, toChr(ng0) → position chromosomique 1-based } ou null.
 */
export function anchorReference(ngSeq, chrSeq, chrStart, anchors) {
  const k = 40, found = [];
  for (const a of anchors) {
    if (a < 0 || a + k > ngSeq.length) continue;
    const kmer = ngSeq.slice(a, a + k);
    if (/N/.test(kmer)) continue;
    const plus = chrSeq.indexOf(kmer);
    if (plus >= 0 && chrSeq.indexOf(kmer, plus + 1) < 0) { found.push({ strand: '+', offset: chrStart + plus - a }); continue; }
    const minus = chrSeq.indexOf(revcomp(kmer));
    // NG base a (début du k-mer) ↔ dernière base du k-mer complémentaire sur le chromosome.
    if (minus >= 0) found.push({ strand: '-', offset: chrStart + minus + k - 1 + a });
  }
  if (found.length < 2) return null;
  // Correspondance retenue : celle que partagent le plus d'ancres (au moins 2, et la majorité).
  const best = found.map(f => ({ ...f, n: found.filter(g => g.strand === f.strand && g.offset === f.offset).length })).sort((x, y) => y.n - x.n)[0];
  if (best.n < 2 || best.n <= found.length / 2) return null;
  const { strand, offset } = best;
  return { strand, toChr: ng0 => (strand === '+' ? offset + ng0 : offset - ng0) };
}

const mappings = new WeakMap();

/** Correspondance NG_ → GRCh38 pour une référence (mise en cache pour la session). */
export async function genomicMapping(ref) {
  if (mappings.has(ref)) return mappings.get(ref);
  const symbol = geneSymbol(ref);
  if (!symbol) throw new Error('symbole du gène introuvable dans la référence');
  const gene = await getJSON(`${ENSEMBL}/lookup/symbol/homo_sapiens/${encodeURIComponent(symbol)}?content-type=application/json`);
  if (!gene?.seq_region_name) throw new Error(`gène ${symbol} inconnu d'Ensembl`);
  const margin = 30000;
  const start = Math.max(1, gene.start - margin), end = gene.end + margin;
  if (end - start > 5e6) throw new Error('région génomique trop grande');
  const response = await fetch(`${ENSEMBL}/sequence/region/human/${gene.seq_region_name}:${start}..${end}:1?content-type=text/plain`, { credentials: 'omit', referrerPolicy: 'no-referrer' });
  if (!response.ok) throw new Error(`rest.ensembl.org HTTP ${response.status}`);
  const chrSeq = (await response.text()).trim().toUpperCase();
  // Points d'ancrage : débuts d'exons du transcrit, sinon répartis le long de la référence.
  const anchors = ref.transcript
    ? ref.transcript.exons.map(([s]) => s - 1)
    : [0.2, 0.4, 0.6, 0.8].map(x => Math.floor(ref.seq.length * x));
  const map = anchorReference(ref.seq, chrSeq, start, anchors);
  if (!map) throw new Error(`${ref.accession} non retrouvée sur le chromosome ${gene.seq_region_name} (GRCh38)`);
  const result = { chrom: String(gene.seq_region_name), chrSeq, chrStart: start, ...map };
  mappings.set(ref, result);
  return result;
}

/**
 * Identifiant gnomAD (chrom-pos-ref-alt, VCF, GRCh38) d'un variant {pos, del, ins} de la référence NG_.
 * `map` : { chrom, chrSeq, chrStart, strand, toChr }.
 */
export function gnomadId(v, map) {
  const { chrSeq, chrStart, strand, toChr, chrom } = map;
  let pos0, del, ins; // index 0-based dans chrSeq, allèles en orientation chromosomique
  if (strand === '+') {
    pos0 = toChr(v.pos) - chrStart; del = v.del; ins = v.ins;
  } else {
    // Bloc NG [pos, pos+len) → chromosome [toChr(pos+len-1), toChr(pos)] ; insertion avant la base NG pos
    // → après la base chromosomique toChr(pos).
    pos0 = (v.del.length ? toChr(v.pos + v.del.length - 1) : toChr(v.pos) + 1) - chrStart;
    del = revcomp(v.del); ins = revcomp(v.ins);
  }
  if (pos0 < 1 || pos0 + del.length > chrSeq.length) return null;
  if (chrSeq.slice(pos0, pos0 + del.length) !== del) return null; // la référence ne correspond pas
  // Normalisation à gauche (convention VCF / gnomAD).
  while (del.length && ins.length && del.at(-1) === ins.at(-1)) { del = del.slice(0, -1); ins = ins.slice(0, -1); }
  while (del.length && ins.length && del[0] === ins[0]) { del = del.slice(1); ins = ins.slice(1); pos0++; }
  if (!del.length || !ins.length) {
    let allele = del || ins;
    while (pos0 > 1 && chrSeq[pos0 - 1] === allele.at(-1)) { allele = chrSeq[pos0 - 1] + allele.slice(0, -1); pos0--; }
    if (del.length) del = allele; else ins = allele;
    const anchor = chrSeq[pos0 - 1];
    return `${chrom}-${chrStart + pos0 - 1}-${anchor + del}-${anchor + ins}`;
  }
  return `${chrom}-${chrStart + pos0}-${del}-${ins}`;
}

/**
 * Fréquences gnomAD v4 d'une liste d'identifiants. Retourne Map(id → { ac, an, af, exome, genome } | null),
 * null signifiant « absent de gnomAD ».
 */
export async function gnomadFrequencies(ids) {
  const out = new Map();
  for (let i = 0; i < ids.length; i += 20) {
    const batch = ids.slice(i, i + 20);
    const fields = batch.map((id, j) => `v${j}: variant(variantId: ${JSON.stringify(id)}, dataset: ${GNOMAD_DATASET}) { variant_id exome { ac an } genome { ac an } }`).join('\n');
    const body = await getJSON(GNOMAD, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query: `{\n${fields}\n}` }),
    });
    if (!body.data && body.errors?.length && !body.errors.every(e => /not found/i.test(e.message))) {
      throw new Error('gnomAD : ' + body.errors[0].message);
    }
    batch.forEach((id, j) => {
      const v = body.data?.[`v${j}`];
      if (!v) { out.set(id, null); return; }
      const ac = (v.exome?.ac || 0) + (v.genome?.ac || 0), an = (v.exome?.an || 0) + (v.genome?.an || 0);
      out.set(id, { ac, an, af: an ? ac / an : 0, exome: v.exome, genome: v.genome });
    });
  }
  return out;
}

export const gnomadUrl = id => `https://gnomad.broadinstitute.org/variant/${id}?dataset=${GNOMAD_DATASET}`;
