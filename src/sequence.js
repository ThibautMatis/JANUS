// Utilitaires de séquence nucléotidique.

export const DNA = 'ACGT';

const COMP = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', R: 'Y', Y: 'R', K: 'M', M: 'K', S: 'S', W: 'W', B: 'V', V: 'B', D: 'H', H: 'D', '-': '-' };

export const complement = b => COMP[b] || 'N';

export function revcomp(s) {
  let out = '';
  for (let i = s.length - 1; i >= 0; i--) out += complement(s[i]);
  return out;
}

const IUPAC = { AG: 'R', CT: 'Y', GT: 'K', AC: 'M', CG: 'S', AT: 'W' };

/** Code IUPAC pour deux bases (identiques → la base elle-même). */
export function iupac(a, b) {
  if (a === b) return a;
  return IUPAC[[a, b].sort().join('')] || 'N';
}
