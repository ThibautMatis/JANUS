// Alignement local à pénalités de gap affines (Gotoh) avec fonction de score libre,
// et localisation rapide d'une lecture sur une longue référence par k-mers.

const NEG = -1e9;

/**
 * Aligne une requête de longueur n sur une cible de longueur m.
 * score(i, j) : score de l'appariement requête[i] / cible[j].
 * Un gap de longueur L coûte gapOpen + gapExtend * L (valeurs négatives).
 * Retourne { score, qStart, qEnd, tStart, tEnd, ops } avec ops = [{ op: 'M'|'I'|'D', q, t }]
 * ('I' = base de la requête absente de la cible, 'D' = base de la cible absente de la requête).
 */
export function alignLocal(n, m, score, { gapOpen = -8, gapExtend = -2 } = {}) {
  if (n * m > 4e7) throw new Error('Alignement trop grand');
  const w = m + 1;
  const tb = new Uint8Array((n + 1) * w);
  let Mp = new Float64Array(w).fill(NEG), Xp = new Float64Array(w).fill(NEG), Yp = new Float64Array(w).fill(NEG);
  let Mc = new Float64Array(w), Xc = new Float64Array(w), Yc = new Float64Array(w);
  Mp[0] = NEG;
  let best = 0, bi = 0, bj = 0;
  const go = gapOpen + gapExtend, ge = gapExtend;

  for (let i = 1; i <= n; i++) {
    Mc[0] = NEG; Xc[0] = NEG; Yc[0] = NEG;
    const row = i * w;
    for (let j = 1; j <= m; j++) {
      // M : fin sur un appariement.
      const dm = Mp[j - 1], dx = Xp[j - 1], dy = Yp[j - 1];
      let prev = 0, from = 0;
      if (dm > prev) { prev = dm; from = 1; }
      if (dx > prev) { prev = dx; from = 2; }
      if (dy > prev) { prev = dy; from = 3; }
      const mv = prev + score(i - 1, j - 1);
      // X : base de requête i en face d'un gap (vient de la ligne précédente).
      let xv = Mp[j] + go, xf = 0;
      if (Xp[j] + ge > xv) { xv = Xp[j] + ge; xf = 1; }
      if (Yp[j] + go > xv) { xv = Yp[j] + go; xf = 2; }
      // Y : base de cible j en face d'un gap (vient de la colonne précédente).
      let yv = Mc[j - 1] + go, yf = 0;
      if (Yc[j - 1] + ge > yv) { yv = Yc[j - 1] + ge; yf = 1; }
      if (Xc[j - 1] + go > yv) { yv = Xc[j - 1] + go; yf = 2; }
      Mc[j] = mv; Xc[j] = xv; Yc[j] = yv;
      tb[row + j] = from | (xf << 2) | (yf << 4);
      if (mv > best) { best = mv; bi = i; bj = j; }
    }
    [Mp, Mc] = [Mc, Mp]; [Xp, Xc] = [Xc, Xp]; [Yp, Yc] = [Yc, Yp];
  }

  const ops = [];
  let i = bi, j = bj, state = 1; // 1 = M, 2 = X, 3 = Y
  while (i > 0 && j > 0) {
    const t = tb[i * w + j];
    if (state === 1) {
      ops.push({ op: 'M', q: i - 1, t: j - 1 });
      const from = t & 3;
      i--; j--;
      if (from === 0) break;
      state = from;
    } else if (state === 2) {
      ops.push({ op: 'I', q: i - 1, t: j });
      const xf = (t >> 2) & 3;
      i--;
      state = xf === 0 ? 1 : xf === 1 ? 2 : 3;
    } else {
      ops.push({ op: 'D', q: i, t: j - 1 });
      const yf = (t >> 4) & 3;
      j--;
      state = yf === 0 ? 1 : yf === 1 ? 3 : 2;
    }
  }
  ops.reverse();
  const first = ops[0], last = ops[ops.length - 1];
  return {
    score: best,
    ops,
    qStart: first ? first.q : 0,
    qEnd: bi,
    tStart: first ? first.t : 0,
    tEnd: bj,
  };
}

/** Alignement local de deux chaînes (match/mismatch simples ; N neutre). */
export function alignStrings(query, target, { match = 3, mismatch = -4, gapOpen = -8, gapExtend = -2 } = {}) {
  return alignLocal(query.length, target.length, (i, j) => {
    const a = query[i], b = target[j];
    if (a === 'N' || b === 'N') return 0;
    return a === b ? match : mismatch;
  }, { gapOpen, gapExtend });
}

/** Index de k-mers d'une séquence. */
export function kmerIndex(seq, k = 12) {
  const index = new Map();
  for (let i = 0; i + k <= seq.length; i++) {
    const kmer = seq.slice(i, i + k);
    if (kmer.includes('N')) continue;
    const list = index.get(kmer);
    if (list) { if (list.length < 8) list.push(i); } else index.set(kmer, [i]);
  }
  return index;
}

/**
 * Diagonale (position cible − position requête) la plus soutenue par des k-mers communs.
 * Retourne { diagonal, votes } ou null.
 */
export function locate(query, index, k = 12) {
  const votes = new Map();
  for (let i = 0; i + k <= query.length; i++) {
    const hits = index.get(query.slice(i, i + k));
    if (!hits || hits.length > 4) continue;
    for (const t of hits) {
      // Regroupe les diagonales proches (indels) par paquets de 16.
      const d = Math.round((t - i) / 16);
      votes.set(d, (votes.get(d) || 0) + 1);
    }
  }
  let best = null;
  for (const [d, v] of votes) {
    const total = v + (votes.get(d - 1) || 0) + (votes.get(d + 1) || 0);
    if (!best || total > best.votes) best = { diagonal: d * 16, votes: total };
  }
  return best;
}
