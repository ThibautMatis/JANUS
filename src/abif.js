// Lecteur ABIF (.ab1) sans dépendance. Les fichiers sont non fiables : toutes les lectures sont bornées.

const TYPE_BYTE = 1, TYPE_CHAR = 2, TYPE_WORD = 3, TYPE_SHORT = 4, TYPE_LONG = 5;

export function parseABIF(buffer) {
  const view = new DataView(buffer);
  const bytes = new Uint8Array(buffer);
  const ascii = new TextDecoder('latin1');

  const check = (p, n) => {
    if (!Number.isSafeInteger(p) || !Number.isSafeInteger(n) || p < 0 || n < 0 || p + n > buffer.byteLength) {
      throw new Error('Fichier ABI tronqué ou invalide');
    }
  };
  const u16 = p => (check(p, 2), view.getUint16(p, false));
  const u32 = p => (check(p, 4), view.getUint32(p, false));

  check(0, 34);
  if (ascii.decode(bytes.subarray(0, 4)) !== 'ABIF') throw new Error('Signature ABIF absente : ce n’est pas un fichier .ab1');

  // L'entrée racine du répertoire commence à l'octet 6.
  const count = u32(18), offset = u32(26);
  if (count > 10000) throw new Error('Répertoire ABIF trop grand');
  check(offset, count * 28);

  const entries = new Map();
  for (let i = 0; i < count; i++) {
    const p = offset + i * 28;
    const name = ascii.decode(bytes.subarray(p, p + 4));
    const num = u32(p + 4), type = u16(p + 8), size = u16(p + 10), n = u32(p + 12), dataSize = u32(p + 16);
    const at = dataSize <= 4 ? p + 20 : u32(p + 20);
    check(at, dataSize);
    entries.set(name + num, { type, size, n, dataSize, at });
  }

  const has = (name, num) => entries.has(name + num);
  const entry = (name, num) => {
    const e = entries.get(name + num);
    if (!e) throw new Error('Tag ABI absent : ' + name + num);
    return e;
  };
  const chars = (name, num) => {
    const e = entry(name, num);
    let s = e.at, len = e.dataSize;
    // Les chaînes pString (type 18) commencent par un octet de longueur.
    if (e.type === 18 && len > 0) { s += 1; len = Math.min(len - 1, bytes[e.at]); }
    return ascii.decode(bytes.subarray(s, s + len)).replace(/\0/g, '');
  };
  const numbers = (name, num) => {
    const e = entry(name, num);
    if (e.n > 2e7 || e.size < 1 || e.n * e.size > e.dataSize) throw new Error('Tag numérique ABI invalide : ' + name + num);
    const out = new Array(e.n);
    for (let i = 0; i < e.n; i++) {
      const p = e.at + i * e.size;
      if (e.type === TYPE_SHORT) out[i] = view.getInt16(p, false);
      else if (e.type === TYPE_WORD || e.size === 2) out[i] = view.getUint16(p, false);
      else if (e.type === TYPE_LONG) out[i] = view.getInt32(p, false);
      else if (e.size === 4) out[i] = view.getUint32(p, false);
      else if (e.type === TYPE_CHAR || e.type === TYPE_BYTE || e.size === 1) out[i] = bytes[p];
      else throw new Error('Type ABI non géré : ' + name + num);
    }
    return out;
  };

  // Version 2 = corrigée par l'utilisateur / le basecaller KB ; version 1 = appel d'origine.
  let bases = has('PBAS', 2) ? chars('PBAS', 2) : '';
  let positions = has('PLOC', 2) ? numbers('PLOC', 2) : [];
  if (!bases || positions.length !== bases.length) {
    bases = chars('PBAS', 1);
    positions = numbers('PLOC', 1);
  }
  bases = bases.toUpperCase();
  if (!/^[ACGTNRYKMSWBDHV]*$/.test(bases)) throw new Error('Séquence PBAS invalide');
  if (positions.length !== bases.length) throw new Error('PBAS et PLOC de longueurs différentes');

  const order = chars('FWO_', 1).slice(0, 4).toUpperCase();
  if (!/^[ACGT]{4}$/.test(order) || new Set(order).size !== 4) throw new Error('Ordre des canaux FWO_ invalide');

  // DATA 9–12 = signal analysé (DATA 1–4 = brut).
  const traces = {};
  for (let i = 0; i < 4; i++) traces[order[i]] = Float32Array.from(numbers('DATA', 9 + i), x => Math.max(0, x));
  const length = traces.A.length;
  if (!length || Object.values(traces).some(t => t.length !== length)) throw new Error('Longueurs des canaux incohérentes');
  if (positions.some(p => p < 0 || p >= length)) throw new Error('Positions PLOC hors du signal');

  let quality = [];
  if (has('PCON', 2)) quality = numbers('PCON', 2);
  else if (has('PCON', 1)) quality = numbers('PCON', 1);
  if (quality.length !== bases.length) quality = [];

  let sample = '';
  try { if (has('SMPL', 1)) sample = chars('SMPL', 1); } catch { /* facultatif */ }

  return { bases, positions, traces, quality, order, sample };
}
