<img src="img/janus-logo.png" alt="JANUS — Sanger deconvolution & variant calling" height="90">

# JANUS — *Joint Allele aNalysis Using Sanger*

**v0.8.1 · 28 octobre 2026**

Outil web d’appel de variants sur chromatogrammes Sanger (`.ab1`) : il sépare les deux allèles d’un échantillon hétérozygote, nomme les SNV et les indels (y compris les grands indels) en HGVS `c.`, `p.` et `g.`, et produit un rapport imprimable.

Tout se passe dans le navigateur : les chromatogrammes ne sont jamais envoyés à un serveur.

**Accès :** https://thibautmatis.github.io/JANUS/

> ⚠ Outil de recherche, non validé pour le diagnostic. Tout variant doit être confirmé sur le chromatogramme et par une méthode validée.

## Utilisation

1. **Référence** : taper le gène (ex. `BRCA2`) puis cliquer sur « Charger ». Le RefSeqGene (NG_) et le transcrit de référence sont choisis automatiquement, et la référence reste mémorisée dans le navigateur.
2. **Chromatogrammes** : importer les `.ab1` Forward et Reverse, puis cliquer sur « Analyser ». La lisibilité de chaque lecture est indiquée.
3. **Variants** : nomenclatures `c.`, `p.` et `g.`, génotype, fraction allélique, confirmation F/R. Le tableau se trie par colonne, et le bouton « Fréquences gnomAD » interroge gnomAD v4.
4. **Chromatogramme** : cliquer sur un variant pour centrer la vue dessus. Au survol, la position `c.` de chaque base s’affiche.
5. **Rapport** : choisir les variants, les classer (classes 1 à 5), vérifier la conclusion, puis cliquer sur « Imprimer / PDF ».

Un jeu de démonstration fictif est fourni dans `demo/` (référence `NG_DEMO.gb`, à charger via « Sans connexion à NCBI », et deux lectures portant `c.752dup`).

En local : `python3 -m http.server`, puis ouvrir http://localhost:8000.

## Références

- Méthode de déconvolution inspirée de Tracy : Rausch T, Fritz MH-Y, Untergasser A, Benes V. *BMC Genomics*. 2020;21:230. doi:10.1186/s12864-020-6635-8
- Classification : Richards S et al. *Genet Med*. 2015;17:405-424 (ACMG/AMP).
- Nomenclature : den Dunnen JT et al. *Hum Mutat*. 2016;37:564-569 (HGVS).
