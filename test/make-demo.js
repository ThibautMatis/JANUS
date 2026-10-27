// Génère un jeu de démonstration entièrement synthétique (aucune donnée réelle) :
// demo/NG_DEMO.gb (gène fictif) et deux chromatogrammes portant une duplication hétérozygote.
//   node test/make-demo.js
import { writeFileSync, mkdirSync } from 'node:fs';
import { simulate, writeABIF, rng, randomSeq } from './synth.js';
import { revcomp } from '../src/sequence.js';

const SEQ = randomSeq(6000, rng(2024));
const exons = [[201, 400], [1001, 1300], [2001, 2600], [3001, 3800]];
const cds = [[249, 400], [1001, 1300], [2001, 2600], [3001, 3151]]; // 1203 nt = 401 codons
// Force un ATG au début du CDS et un codon stop à la fin, pour un transcrit plausible.
const chars = SEQ.split('');
chars.splice(248, 3, 'A', 'T', 'G');
chars.splice(3148, 3, 'T', 'A', 'A');
// Cadre de lecture ouvert : les codons stop internes du CDS sont remplacés (TGG, Trp).
{
  const pos = cds.flatMap(([a, b]) => Array.from({ length: b - a + 1 }, (_, i) => a - 1 + i));
  for (let i = 0; i + 3 < pos.length - 3; i += 3) {
    const codon = pos.slice(i, i + 3).map(p => chars[p]).join('');
    if (['TAA', 'TAG', 'TGA'].includes(codon)) pos.slice(i, i + 3).forEach((p, j) => { chars[p] = 'TGG'[j]; });
  }
}
const seq = chars.join('');

const join = iv => `join(${iv.map(([a, b]) => `${a}..${b}`).join(',')})`;
const origin = [];
for (let i = 0; i < seq.length; i += 60) {
  const chunk = seq.slice(i, i + 60).toLowerCase().match(/.{1,10}/g).join(' ');
  origin.push(String(i + 1).padStart(9) + ' ' + chunk);
}
const gb = `LOCUS       NG_DEMO               ${seq.length} bp    DNA     linear   SYN 01-JAN-2026
DEFINITION  Gène fictif DEMO1, séquence aléatoire pour démonstration de SVIB.
ACCESSION   NG_DEMO
VERSION     NG_DEMO.1
FEATURES             Location/Qualifiers
     source          1..${seq.length}
                     /organism="synthetic construct"
     gene            201..3800
                     /gene="DEMO1"
     mRNA            ${join(exons)}
                     /gene="DEMO1"
                     /product="demo protein"
                     /transcript_id="NM_DEMO.1"
     CDS             ${join(cds)}
                     /gene="DEMO1"
                     /protein_id="NP_DEMO.1"
ORIGIN
${origin.join('\n')}
//
`;

// Duplication hétérozygote d'une base dans l'exon 3, amplicon g.1850–2650.
let g = 2300;
while (seq[g - 1] === seq[g]) g++;
const wt = seq.slice(1849, 2650);
const i = g - 1850 + 1;
const mut = wt.slice(0, i) + seq[g - 1] + wt.slice(i);
const F = simulate([wt.slice(20), mut.slice(20)], [0.5, 0.5], { seed: 11, length: 700, sample: 'DEMO_F' });
const R = simulate([revcomp(wt).slice(20), revcomp(mut).slice(20)], [0.5, 0.5], { seed: 12, length: 700 });

mkdirSync('demo', { recursive: true });
writeFileSync('demo/NG_DEMO.gb', gb);
writeFileSync('demo/DEMO_F.ab1', Buffer.from(writeABIF({ ...F, sample: 'DEMO_F' })));
writeFileSync('demo/DEMO_R.ab1', Buffer.from(writeABIF({ ...R, sample: 'DEMO_R' })));
// CDS au nt 49 du transcrit ; exon 3 (g.2001) = nt 501 du transcrit = c.453.
console.log(`Variant simulé : NG_DEMO.1:g.${g}dup = NM_DEMO.1:c.${g - 2001 + 453}dup`);
