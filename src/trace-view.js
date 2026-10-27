// Affichage d'un chromatogramme sur <canvas> : déplacement, zoom, sélection d'une base,
// avec sous les pics la référence et les deux allèles reconstruits.

// G est lu dans --trace-g pour rester visible en thème sombre.
const COLORS = { A: '#1a8a4f', C: '#2463c7', G: '#30353b', T: '#d0383a' };
const OTHER = '#8c60a7';

export class TraceView {
  /** onSelect(k) : clic sur une base ; label(k) → [titre, détail] affichés au survol (barre verticale). */
  constructor(canvas, { onSelect, label } = {}) {
    this.canvas = canvas;
    this.onSelect = onSelect;
    this.label = label;
    this.hover = null;
    this.read = null;
    this.start = 0;
    this.count = 45;
    this.annotations = new Map();
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

  /** read : { bases, positions, traces } ; annotations : Map(indexBase → {ref, a1, a2}) ; markers : [indexBase]. */
  setRead(read, annotations = new Map(), markers = []) {
    this.read = read;
    this.annotations = annotations;
    this.markers = markers;
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

  draw() {
    const { canvas, read } = this;
    const dpr = window.devicePixelRatio || 1;
    const W = canvas.clientWidth, H = canvas.clientHeight;
    if (!W || !H) return;
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    if (!read) return;

    const css = getComputedStyle(canvas);
    COLORS.G = css.getPropertyValue('--trace-g').trim() || COLORS.G;
    const muted = css.getPropertyValue('--muted').trim() || '#647587';
    const accent = css.getPropertyValue('--accent').trim() || '#147b95';
    const rows = this.annotations.size ? 3 : 0;
    const labelH = 18, rowH = 15, bottom = 26 + labelH + rows * rowH;
    const plotH = H - bottom - 8;

    const end = Math.min(read.bases.length, this.start + this.count);
    const x0 = Math.max(0, read.positions[this.start] - 8);
    const x1 = Math.min(read.traces.A.length - 1, read.positions[end - 1] + 8);
    const span = Math.max(1, x1 - x0);
    const X = s => ((s - x0) / span) * (W - 70) + 60;
    let ymax = 1;
    for (const b of 'ACGT') for (let i = x0; i <= x1; i++) if (read.traces[b][i] > ymax) ymax = read.traces[b][i];
    this.geometry = { x0, x1, span, end, W };

    // Variants signalés.
    ctx.fillStyle = accent + '22';
    for (const m of this.markers) {
      if (m < this.start || m >= end) continue;
      const half = Math.max(4, ((read.positions[Math.min(m + 1, end - 1)] - read.positions[Math.max(m - 1, this.start)]) / span) * (W - 70) / 4);
      ctx.fillRect(X(read.positions[m]) - half, 0, 2 * half, H - 6);
    }

    ctx.lineWidth = 1.4;
    for (const b of 'ACGT') {
      ctx.beginPath();
      ctx.strokeStyle = COLORS[b];
      const t = read.traces[b];
      for (let i = x0; i <= x1; i++) {
        const x = X(i), y = 8 + plotH - (t[i] / ymax) * plotH;
        if (i === x0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    const perBase = (W - 70) / this.count;
    const step = perBase < 3 ? 50 : perBase < 6 ? 20 : perBase < 12 ? 10 : 5;
    const showLetters = this.count <= (W - 70) / 9;
    const yBase = H - bottom + labelH / 2 + 4;
    for (let k = this.start; k < end; k++) {
      const x = X(read.positions[k]);
      if (showLetters) {
        ctx.font = '600 13px ui-monospace, SFMono-Regular, Menlo, monospace';
        ctx.fillStyle = COLORS[read.bases[k]] || OTHER;
        ctx.fillText(read.bases[k], x, yBase);
        const a = this.annotations.get(k);
        if (a) {
          ctx.font = '12px ui-monospace, SFMono-Regular, Menlo, monospace';
          [a.ref, a.a1, a.a2].forEach((b, r) => {
            const differs = r > 0 && b !== a.ref;
            ctx.fillStyle = differs ? accent : (COLORS[b] || muted);
            if (differs) ctx.font = '700 12px ui-monospace, SFMono-Regular, Menlo, monospace';
            ctx.fillText(b, x, yBase + labelH + r * rowH);
            ctx.font = '12px ui-monospace, SFMono-Regular, Menlo, monospace';
          });
        }
      }
      if ((k + 1) % step === 0) {
        ctx.font = '10px system-ui, sans-serif';
        ctx.fillStyle = muted;
        ctx.fillText(String(k + 1), x, H - 10);
      }
    }
    ctx.textAlign = 'left';
    ctx.font = '11px system-ui, sans-serif';
    ctx.fillStyle = muted;
    ctx.fillText('Lecture', 4, yBase);
    if (rows) ['Réf.', 'Allèle 1', 'Allèle 2'].forEach((t, r) => ctx.fillText(t, 4, yBase + labelH + r * rowH));

    if (this.hover !== null && this.hover >= this.start && this.hover < end) this.drawHover(ctx, X(read.positions[this.hover]), W, H, css, accent);
  }

  /** Barre verticale sur la base survolée et étiquette (position c., détail). */
  drawHover(ctx, x, W, H, css, accent) {
    ctx.save();
    ctx.strokeStyle = accent;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.moveTo(Math.round(x) + 0.5, 0);
    ctx.lineTo(Math.round(x) + 0.5, H);
    ctx.stroke();
    ctx.setLineDash([]);
    const [title, detail] = this.label ? this.label(this.hover) : [`base ${this.hover + 1}`, ''];
    ctx.font = '700 13px ui-monospace, SFMono-Regular, Menlo, monospace';
    const w1 = ctx.measureText(title).width;
    ctx.font = '11px system-ui, sans-serif';
    const w2 = detail ? ctx.measureText(detail).width : 0;
    const bw = Math.max(w1, w2) + 16, bh = detail ? 38 : 22;
    const bx = Math.min(Math.max(2, x + 8), W - bw - 2), by = 4;
    ctx.fillStyle = css.getPropertyValue('--card').trim() || '#fff';
    ctx.strokeStyle = accent;
    ctx.beginPath();
    ctx.roundRect ? ctx.roundRect(bx, by, bw, bh, 5) : ctx.rect(bx, by, bw, bh);
    ctx.fill();
    ctx.stroke();
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillStyle = css.getPropertyValue('--text').trim() || '#17293a';
    ctx.font = '700 13px ui-monospace, SFMono-Regular, Menlo, monospace';
    ctx.fillText(title, bx + 8, by + 11);
    if (detail) {
      ctx.font = '11px system-ui, sans-serif';
      ctx.fillStyle = css.getPropertyValue('--muted').trim() || '#647587';
      ctx.fillText(detail, bx + 8, by + 28);
    }
    ctx.restore();
  }

  baseAt(clientX) {
    const g = this.geometry;
    if (!g || !this.read) return null;
    const rect = this.canvas.getBoundingClientRect();
    const s = g.x0 + ((clientX - rect.left - 60) / (g.W - 70)) * g.span;
    let best = this.start;
    for (let k = this.start; k < g.end; k++) {
      if (Math.abs(this.read.positions[k] - s) < Math.abs(this.read.positions[best] - s)) best = k;
    }
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
        // Survol à la souris : barre verticale et position de la base (pas au doigt, réservé au glisser).
        if (e.pointerType === 'touch' || !this.read) return;
        const k = this.baseAt(e.clientX);
        if (k !== this.hover) { this.hover = k; this.draw(); }
        return;
      }
      if (this.hover !== null) this.hover = null;
      const dx = e.clientX - drag.x;
      if (Math.abs(dx) > 3) drag.moved = true;
      if (drag.moved) this.show(drag.start - (dx / c.clientWidth) * this.count, this.count);
    });
    c.addEventListener('pointerup', e => {
      if (drag && !drag.moved && this.onSelect) {
        const k = this.baseAt(e.clientX);
        if (k !== null) this.onSelect(k);
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
