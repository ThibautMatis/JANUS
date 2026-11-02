// Affichage d'un chromatogramme sur <canvas> : déplacement, zoom, sélection d'une base,
// avec sous les pics la séquence du génome et les deux allèles reconstruits.
//
// La ligne Génome est continue : chaque base du génome a sa colonne. Une base du génome sans pic en
// face (délétion par rapport à la lecture) reçoit une colonne insérée entre deux pics : le tracé s'y
// écarte d'autant, comme dans un alignement, et les lignes de la lecture y portent un « – ».

// G est lu dans --trace-g pour rester visible en thème sombre.
const COLORS = { A: '#1a8a4f', C: '#2463c7', G: '#30353b', T: '#d0383a' };
const OTHER = '#8c60a7';
const MONO = 'ui-monospace, SFMono-Regular, Menlo, monospace';

export class TraceView {
  /**
   * onSelect(k) : clic sur une base de la lecture ;
   * label(colonne) → [titre, détail] au survol (colonne = { k } pour une base lue, { gap } sinon) ;
   * context(g) → { region: { kind, n }, side, codon: { n, phase, aa } } pour la position génomique g
   * (bandeau exon / intron, acides aminés de référence, limites des codons), ou null.
   */
  constructor(canvas, { onSelect, label, context } = {}) {
    this.canvas = canvas;
    this.onSelect = onSelect;
    this.label = label;
    this.context = context;
    this.hover = null;
    this.read = null;
    this.start = 0;
    this.count = 45;
    this.annotations = new Map();
    this.gaps = new Map();
    this.markers = [];
    this.bindEvents();
    // Au redimensionnement, on garde la même densité de bases par pixel autour du centre de la vue.
    let width = canvas.clientWidth;
    new ResizeObserver(() => {
      const w = canvas.clientWidth;
      if (this.read && width && w && w !== width) {
        const center = this.start + this.count / 2, count = this.count * w / width;
        this.count = Math.max(8, Math.min(this.read.bases.length, Math.round(count)));
        this.start = Math.max(0, Math.min(this.read.bases.length - this.count, Math.round(center - this.count / 2)));
      }
      width = w;
      this.draw();
    }).observe(canvas);
  }

  /**
   * read : { bases, positions, traces } ;
   * annotations : Map(indexBase → { genome, a1, a2 }) — génome « - » pour une base lue absente du génome ;
   * gaps : Map(indexBase → [{ genome, g }]) — bases du génome sans pic, placées après cette base lue ;
   * markers : [indexBase] surlignés.
   */
  setRead(read, annotations = new Map(), markers = [], gaps = new Map()) {
    Object.assign(this, { read, annotations, markers, gaps });
    this.show(0, Math.min(this.comfortable(), read.bases.length));
  }

  show(start, count) {
    if (!this.read) return;
    const n = this.read.bases.length;
    this.count = Math.max(8, Math.min(n, Math.round(count)));
    this.start = Math.max(0, Math.min(n - this.count, Math.round(start)));
    this.draw();
  }

  /** Nombre de bases lisibles à la largeur actuelle (≈ 16 px par base). */
  comfortable() { return Math.max(15, Math.min(60, Math.floor((this.canvas.clientWidth - 70) / 16))); }

  centerOn(index, count = this.comfortable()) { this.show(index - count / 2, count); }
  zoom(factor, anchor = 0.5) {
    const next = Math.max(8, Math.round(this.count * factor));
    this.show(this.start + anchor * (this.count - next), next);
  }

  /**
   * Coordonnée « virtuelle » d'un point du tracé : position en échantillons, augmentée d'un espacement
   * de base pour chaque colonne de génome insérée avant lui.
   */
  virtualScale() {
    const { positions } = this.read;
    const diffs = [];
    for (let k = 1; k < positions.length; k++) diffs.push(positions[k] - positions[k - 1]);
    diffs.sort((a, b) => a - b);
    const spacing = diffs[Math.floor(diffs.length / 2)] || 12;
    const cuts = [...this.gaps.entries()]
      .filter(([k]) => k + 1 < positions.length)
      .map(([k, cols]) => ({ k, at: (positions[k] + positions[k + 1]) / 2, n: cols.length }))
      .sort((a, b) => a.at - b.at);
    let before = 0;
    for (const c of cuts) { c.before = before; before += c.n; }
    const v = s => {
      let extra = 0;
      for (const c of cuts) { if (c.at < s) extra = c.before + c.n; else break; }
      return s + extra * spacing;
    };
    return { v, cuts, spacing };
  }

  draw() {
    const { canvas, read } = this;
    // width / height / dpr fixes : rendu hors écran (capture pour le rapport).
    const dpr = this.dpr ?? (window.devicePixelRatio || 1);
    const W = this.width ?? canvas.clientWidth, H = this.height ?? canvas.clientHeight;
    if (!W || !H) return;
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    if (!read) return;

    const css = getComputedStyle(canvas);
    const colors = { ...COLORS, G: css.getPropertyValue('--trace-g').trim() || COLORS.G };
    const muted = css.getPropertyValue('--muted').trim() || '#647587';
    const accent = css.getPropertyValue('--accent').trim() || '#147b95';
    const rows = this.annotations.size ? 3 : 0;
    const labelH = 18, rowH = 15, bottom = 26 + labelH + rows * rowH;
    // Bandeau exon / intron et ligne des acides aminés en haut, si le transcrit est connu.
    const top = this.context ? 36 : 0;
    const plotH = H - bottom - 8 - top;

    const end = Math.min(read.bases.length, this.start + this.count);
    const x0 = Math.max(0, read.positions[this.start] - 8);
    const x1 = Math.min(read.traces.A.length - 1, read.positions[end - 1] + 8);
    const { v, cuts, spacing } = this.virtualScale();
    const v0 = v(x0), v1 = v(x1), span = Math.max(1, v1 - v0);
    const X = u => ((u - v0) / span) * (W - 70) + 60;
    let ymax = 1;
    for (const b of 'ACGT') for (let i = x0; i <= x1; i++) if (read.traces[b][i] > ymax) ymax = read.traces[b][i];

    // Colonnes visibles : bases lues et bases du génome insérées entre elles.
    const layout = [];
    for (let k = this.start; k < end; k++) {
      layout.push({ k, x: X(v(read.positions[k])) });
      const cut = cuts.find(c => c.k === k);
      if (cut && k + 1 < end) {
        const vStart = cut.at + cut.before * spacing;
        this.gaps.get(k).forEach((gap, i) => layout.push({ k: null, gap, x: X(vStart + (i + 0.5) * spacing) }));
      }
    }
    this.layout = layout;

    // Variants signalés.
    ctx.fillStyle = accent + '22';
    const half = Math.max(4, (W - 70) / Math.max(1, layout.length) / 2);
    for (const c of layout) if (c.k !== null && this.markers.includes(c.k)) ctx.fillRect(c.x - half, 0, 2 * half, H - 6);

    ctx.lineWidth = 1.4;
    for (const b of 'ACGT') {
      ctx.beginPath();
      ctx.strokeStyle = colors[b];
      const t = read.traces[b];
      let prevV = null;
      for (let i = x0; i <= x1; i++) {
        const vi = v(i), x = X(vi), y = top + 8 + plotH - (t[i] / ymax) * plotH;
        // Colonne de génome insérée : le tracé s'interrompt au lieu d'être étiré en diagonale.
        if (prevV === null || vi - prevV > 1.5) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        prevV = vi;
      }
      ctx.stroke();
    }

    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    const perCol = (W - 70) / Math.max(1, layout.length);
    if (this.context) this.drawContext(ctx, layout, { W, H, top, perCol, muted, accent, rowsBottom: H - 20 });
    const step = perCol < 3 ? 50 : perCol < 6 ? 20 : perCol < 12 ? 10 : 5;
    const showLetters = perCol >= 9;
    const yBase = H - bottom + labelH / 2 + 4;
    const letter = (b, x, y, color, bold = false) => {
      ctx.font = `${bold ? '700 ' : ''}12px ${MONO}`;
      ctx.fillStyle = color;
      ctx.fillText(b, x, y);
    };
    for (const col of layout) {
      const { x } = col;
      if (showLetters) {
        if (col.k === null) {
          // Base du génome sans pic en face : génome en clair, « – » sur la lecture et les allèles.
          letter('–', x, yBase, muted);
          letter(col.gap.genome, x, yBase + labelH, colors[col.gap.genome] || muted);
          letter('–', x, yBase + labelH + rowH, accent, true);
          letter('–', x, yBase + labelH + 2 * rowH, accent, true);
          continue;
        }
        ctx.font = `600 13px ${MONO}`;
        ctx.fillStyle = colors[read.bases[col.k]] || OTHER;
        ctx.fillText(read.bases[col.k], x, yBase);
        const a = this.annotations.get(col.k);
        if (a) {
          letter(a.genome, x, yBase + labelH, a.genome === '-' ? muted : (colors[a.genome] || muted));
          [a.a1, a.a2].forEach((b, r) => {
            const differs = b !== a.genome;
            letter(b, x, yBase + labelH + (r + 1) * rowH, differs ? accent : (colors[b] || muted), differs);
          });
        }
      }
      if (col.k !== null && (col.k + 1) % step === 0) {
        ctx.font = '10px system-ui, sans-serif';
        ctx.fillStyle = muted;
        ctx.fillText(String(col.k + 1), x, H - 10);
      }
    }
    ctx.textAlign = 'left';
    ctx.font = '11px system-ui, sans-serif';
    ctx.fillStyle = muted;
    ctx.fillText('Lecture', 4, yBase);
    if (rows) ['Génome', 'Allèle 1', 'Allèle 2'].forEach((t, r) => ctx.fillText(t, 4, yBase + labelH + r * rowH));

    if (this.hover !== null && layout[this.hover]) this.drawHover(ctx, layout[this.hover], W, H, css, accent);
  }

  /** Position génomique d'une colonne (base lue alignée ou base du génome insérée), ou null. */
  columnG(col) {
    if (col.k === null) return col.gap.g;
    const a = this.annotations.get(col.k);
    return a && a.genome !== '-' ? a.g : null;
  }

  /**
   * Bandeau exon / intron (avec l'exon le plus proche dans un intron), acides aminés de référence au
   * centre de chaque codon et limites des codons en pointillés discrets.
   */
  drawContext(ctx, layout, { W, top, perCol, muted, accent, rowsBottom }) {
    const ctxs = layout.map(c => { const g = this.columnG(c); return g === null ? null : this.context(g); });
    const half = perCol / 2;
    // Bandeau : segments consécutifs de même région.
    let i = 0;
    while (i < layout.length) {
      const r = ctxs[i]?.region;
      let j = i;
      while (j + 1 < layout.length && ctxs[j + 1]?.region?.kind === r?.kind && ctxs[j + 1]?.region?.n === r?.n) j++;
      if (r) {
        const x0 = Math.max(60, layout[i].x - half), x1 = Math.min(W - 10, layout[j].x + half);
        const exon = r.kind === 'exon';
        ctx.fillStyle = exon ? accent + '33' : muted + '18';
        ctx.fillRect(x0, 2, x1 - x0, 14);
        if (exon) { ctx.fillStyle = accent; ctx.fillRect(x0, 14, x1 - x0, 2); }
        let text = r.kind === 'exon' ? `Exon ${r.n}` : r.kind === 'intron' ? `Intron ${r.n}` : r.kind === 'amont' ? 'En amont du gène' : 'En aval du gène';
        // Dans un intron : exon le plus proche, avec une flèche dans sa direction à l'écran.
        const near = ctxs[Math.round((i + j) / 2)]?.nearest;
        if (r.kind === 'intron' && near) text = near.left ? `← Exon ${near.n} · ${text}` : `${text} · Exon ${near.n} →`;
        ctx.font = `${exon ? '700 ' : ''}11px system-ui, sans-serif`;
        ctx.fillStyle = exon ? accent : muted;
        if (ctx.measureText(text).width < x1 - x0 - 4) ctx.fillText(text, (x0 + x1) / 2, 9);
      }
      i = j + 1;
    }
    // Acides aminés de référence (au milieu du codon) et limites des codons.
    ctx.save();
    ctx.strokeStyle = muted + '55';
    ctx.setLineDash([2, 3]);
    ctx.lineWidth = 1;
    for (let k = 0; k < layout.length; k++) {
      const cd = ctxs[k]?.codon;
      if (!cd) continue;
      if (cd.phase === 1) {
        const name = perCol >= 14 ? `${cd.aa3}${cd.n}` : perCol >= 9 ? cd.aa3 : cd.aa;
        ctx.font = `10px ${MONO}`;
        ctx.fillStyle = cd.aa === '*' ? '#b3261e' : muted;
        ctx.fillText(name, layout[k].x, 27);
      }
      const next = ctxs[k + 1]?.codon;
      if (next && next.n !== cd.n) {
        const xm = (layout[k].x + layout[k + 1].x) / 2;
        ctx.beginPath();
        ctx.moveTo(xm, 20);
        ctx.lineTo(xm, rowsBottom);
        ctx.stroke();
      }
    }
    ctx.restore();
  }

  /** Barre verticale sur la colonne survolée et étiquette (position c., détail). */
  drawHover(ctx, col, W, H, css, accent) {
    const x = col.x;
    ctx.save();
    ctx.strokeStyle = accent;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.moveTo(Math.round(x) + 0.5, 0);
    ctx.lineTo(Math.round(x) + 0.5, H);
    ctx.stroke();
    ctx.setLineDash([]);
    const [title, detail] = this.label ? this.label(col) : [col.k !== null ? `base ${col.k + 1}` : 'génome', ''];
    ctx.font = `700 13px ${MONO}`;
    const w1 = ctx.measureText(title).width;
    ctx.font = '11px system-ui, sans-serif';
    const w2 = detail ? ctx.measureText(detail).width : 0;
    const bw = Math.max(w1, w2) + 16, bh = detail ? 38 : 22;
    const bx = Math.min(Math.max(2, x + 8), W - bw - 2), by = this.context ? 40 : 4; // sous le bandeau
    ctx.fillStyle = css.getPropertyValue('--card').trim() || '#fff';
    ctx.strokeStyle = accent;
    ctx.beginPath();
    ctx.roundRect ? ctx.roundRect(bx, by, bw, bh, 5) : ctx.rect(bx, by, bw, bh);
    ctx.fill();
    ctx.stroke();
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillStyle = css.getPropertyValue('--text').trim() || '#17293a';
    ctx.font = `700 13px ${MONO}`;
    ctx.fillText(title, bx + 8, by + 11);
    if (detail) {
      ctx.font = '11px system-ui, sans-serif';
      ctx.fillStyle = css.getPropertyValue('--muted').trim() || '#647587';
      ctx.fillText(detail, bx + 8, by + 28);
    }
    ctx.restore();
  }

  /** Indice (dans this.layout) de la colonne la plus proche d'une abscisse écran. */
  columnAt(clientX) {
    if (!this.layout?.length) return null;
    const rect = this.canvas.getBoundingClientRect();
    const x = clientX - rect.left;
    let best = 0;
    this.layout.forEach((c, i) => { if (Math.abs(c.x - x) < Math.abs(this.layout[best].x - x)) best = i; });
    return best;
  }

  bindEvents() {
    const c = this.canvas;
    let drag = null;
    c.addEventListener('pointerdown', e => {
      if (!this.read) return;
      drag = { x: e.clientX, start: this.start, moved: false };
      c.setPointerCapture(e.pointerId);
    });
    c.addEventListener('pointermove', e => {
      if (!drag) {
        // Survol à la souris : barre verticale et position de la colonne (pas au doigt, réservé au glisser).
        if (e.pointerType === 'touch' || !this.read) return;
        const i = this.columnAt(e.clientX);
        if (i !== this.hover) { this.hover = i; this.draw(); }
        return;
      }
      if (this.hover !== null) this.hover = null;
      const dx = e.clientX - drag.x;
      if (Math.abs(dx) > 3) drag.moved = true;
      if (drag.moved) this.show(drag.start - (dx / c.clientWidth) * this.count, this.count);
    });
    c.addEventListener('pointerup', e => {
      if (drag && !drag.moved && this.onSelect) {
        const i = this.columnAt(e.clientX);
        const col = i !== null ? this.layout[i] : null;
        if (col && col.k !== null) this.onSelect(col.k);
      }
      drag = null;
    });
    c.addEventListener('pointerleave', () => {
      if (this.hover !== null) { this.hover = null; this.draw(); }
    });
    c.addEventListener('wheel', e => {
      if (!this.read) return;
      e.preventDefault();
      const rect = c.getBoundingClientRect();
      this.zoom(e.deltaY < 0 ? 0.8 : 1.25, (e.clientX - rect.left) / rect.width);
    }, { passive: false });
    c.addEventListener('dblclick', e => {
      const rect = c.getBoundingClientRect();
      this.zoom(0.5, (e.clientX - rect.left) / rect.width);
    });
  }
}

/**
 * Image PNG (data URL) d'un chromatogramme centré sur la base `center`, en couleurs claires,
 * pour le rapport imprimable.
 */
export function snapshot(read, annotations, markers, center, { width = 1000, height = 330, count = 41, gaps = new Map(), context = null } = {}) {
  const canvas = document.createElement('canvas');
  const view = Object.create(TraceView.prototype);
  const n = read.bases.length;
  count = Math.min(n, count);
  Object.assign(view, {
    canvas, read, annotations, markers, gaps, context, width, height, dpr: 2, hover: null, label: null,
    count, start: Math.max(0, Math.min(n - count, Math.round(center - count / 2))),
  });
  view.draw();
  return canvas.toDataURL('image/png');
}
