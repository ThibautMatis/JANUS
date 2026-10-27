// Extraction, à chaque base appelée, de la hauteur de pic des quatre canaux
// puis des bases primaire et secondaire (principe de Tracy, Rausch et al. 2020).

import { DNA, complement } from './sequence.js';

/** Hauteur du pic d'un canal près de `center` : maximum local le plus proche, sinon valeur au centre. */
function peakHeight(trace, center, half) {
  const lo = Math.max(1, center - half), hi = Math.min(trace.length - 2, center + half);
  let best = -1, bestDist = Infinity;
  for (let x = lo; x <= hi; x++) {
    if (trace[x] >= trace[x - 1] && trace[x] >= trace[x + 1] && trace[x] > 0) {
      const d = Math.abs(x - center);
      if (d < bestDist || (d === bestDist && trace[x] > trace[best])) { best = x; bestDist = d; }
    }
  }
  return best >= 0 ? trace[best] : trace[Math.min(trace.length - 1, Math.max(0, center))];
}

/**
 * Profil par base : hauteurs [A,C,G,T], base primaire, secondaire et rapport secondaire/primaire.
 * `index` = indice de la base dans la lecture d'origine (utile après inversion).
 */
export function peakProfile(read) {
  const { positions, traces } = read;
  const out = [];
  for (let k = 0; k < positions.length; k++) {
    const p = positions[k];
    const left = k > 0 ? p - positions[k - 1] : (positions[k + 1] ?? p + 12) - p;
    const right = k + 1 < positions.length ? positions[k + 1] - p : left;
    const half = Math.max(1, Math.min(8, Math.round(Math.min(left, right) * 0.4)));
    const h = DNA.split('').map(b => peakHeight(traces[b], p, half));
    const order = [0, 1, 2, 3].sort((a, b) => h[b] - h[a]);
    const h1 = h[order[0]], h2 = h[order[1]];
    out.push({
      index: k,
      h,
      primary: h1 > 0 ? DNA[order[0]] : 'N',
      secondary: h2 > 0 ? DNA[order[1]] : 'N',
      ratio: h1 > 0 ? h2 / h1 : 0,
      total: h.reduce((a, b) => a + b, 0),
    });
  }
  return out;
}

/** Profil retourné-complémenté (pour travailler dans l'orientation de la référence). */
export function reverseProfile(profile) {
  return profile.slice().reverse().map(x => ({
    ...x,
    h: [x.h[3], x.h[2], x.h[1], x.h[0]],
    primary: complement(x.primary),
    secondary: complement(x.secondary),
  }));
}

/**
 * Bornes [début, fin) de la partie exploitable : on retire les extrémités où l'intensité du signal
 * s'effondre. On ne s'appuie pas sur la qualité du basecaller, qui chute aussi dans les zones
 * mixtes d'un indel hétérozygote — c'est justement le signal que l'on veut garder.
 */
export function trimRange(profile, { minLeft = 15, window = 15, minFraction = 0.12 } = {}) {
  const n = profile.length;
  if (n < 2 * window + minLeft) return [Math.min(minLeft, n), n];
  const totals = profile.map(x => x.total);
  const sorted = totals.slice().sort((a, b) => a - b);
  const ref = sorted[Math.floor(sorted.length * 0.75)] || 1;
  const smooth = [];
  let sum = 0;
  for (let i = 0; i < n; i++) {
    sum += totals[i];
    if (i >= window) sum -= totals[i - window];
    smooth.push(sum / Math.min(i + 1, window));
  }
  // smooth[i] = moyenne de [i-window+1, i]
  let start = minLeft;
  while (start < n - window && smooth[start + window - 1] < minFraction * ref) start++;
  let end = n;
  while (end > start + window && smooth[end - 1] < minFraction * ref) end--;
  return [start, end];
}
