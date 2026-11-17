// Version affichée dans l'interface et le rapport.
export const NAME = 'JANUS';
export const VERSION = '0.8.10';
export const RELEASE_DATE = '2026-11-17';
export const releaseLabel = () => `v${VERSION} · ${new Date(RELEASE_DATE + 'T12:00:00').toLocaleDateString('fr-FR', { day: 'numeric', month: 'long', year: 'numeric' })}`;
