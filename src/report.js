// Rapport : classification des variants retenus et conclusion rédigée automatiquement.

/** Classes ACMG/AMP (Richards et al. 2015), numérotées comme en pratique française. */
export const CLASSES = [
  { id: '1', label: 'bénin', long: 'Bénin (classe 1)' },
  { id: '2', label: 'probablement bénin', long: 'Probablement bénin (classe 2)' },
  { id: '3', label: 'de signification incertaine (VUS)', long: 'Variant de signification incertaine – VUS (classe 3)' },
  { id: '4', label: 'probablement pathogène', long: 'Probablement pathogène (classe 4)' },
  { id: '5', label: 'pathogène', long: 'Pathogène (classe 5)' },
];
export const classOf = id => CLASSES.find(c => c.id === id) || null;

/** Symbole du gène de la référence. */
export const geneOf = ref => ref.transcript?.gene || ref.definition?.match(/\(([A-Za-z0-9-]+)\)/)?.[1] || ref.accession;

/** Région analysée (lectures exploitables), en c. si le transcrit est connu. */
export function analysedRegion(ref, results) {
  const usable = results.filter(r => r.quality?.verdict !== 'illisible');
  if (!usable.length) return null;
  const start = Math.min(...usable.map(r => r.refStart)), end = Math.max(...usable.map(r => r.refEnd));
  return ref.mapper ? `c.${ref.mapper.c(start)}_${ref.mapper.c(end)}` : `g.${start}_${end}`;
}

/**
 * Conclusion : une phrase par variant retenu, {variant, classId}, dans l'ordre fourni ;
 * en l'absence de variant retenu, absence de variant dans la région analysée.
 */
export function conclusion(selected, ref, results) {
  const gene = geneOf(ref);
  const transcript = ref.transcript?.id;
  if (!selected.length) {
    const region = analysedRegion(ref, results);
    return `Absence de variant détecté dans la région analysée du gène ${gene}${transcript ? ` (${transcript}` : ''}${region ? `${transcript ? ', ' : ' ('}${region})` : transcript ? ')' : ''}.`;
  }
  return selected.map(({ variant: v, classId }) => {
    const cls = classOf(classId);
    const name = v.names.c || v.names.g;
    const p = v.protein?.short ? ` ${v.protein.short}` : '';
    const state = v.genotype === 'homozygote' ? 'à l’état homozygote' : v.genotype === 'hétérozygote' ? 'à l’état hétérozygote' : '';
    const confirm = v.status === 'confirmé F/R' ? '' : ' Variant non confirmé sur les deux brins : à vérifier.';
    return `Détection${state ? ` ${state}` : ''} du variant ${cls ? cls.label : '[classification à préciser]'} ${name}${p} dans le gène ${gene}${cls ? ` (classe ${cls.id})` : ''}.${confirm}`;
  }).join('\n');
}
