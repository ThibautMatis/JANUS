// Version affichée dans l'interface et le rapport.
export const NAME = 'JANUS';
export const VERSION = '0.8';
export const RELEASE_DATE = '2026-10-27';
export const releaseLabel = () => `v${VERSION} · ${new Date(RELEASE_DATE + 'T12:00:00').toLocaleDateString('fr-FR', { day: 'numeric', month: 'long', year: 'numeric' })}`;
