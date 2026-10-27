# SVIB — déconvolution d’allèles et appel de variants Sanger

Outil web 100 % local (JavaScript, sans dépendance) qui lit des chromatogrammes `.ab1`, reconstruit les deux allèles d’un échantillon hétérozygote et nomme les variants (SNV, délétions, duplications, insertions) en HGVS `g.` et `c.`.

> ⚠ **Outil de recherche, non validé pour le diagnostic.** Tout appel doit être vérifié sur le chromatogramme et confirmé par une méthode validée. Un résultat vide n’exclut pas un variant.

## Principe (v8)

L’approche s’inspire de Tracy / Indigo (Rausch et al., *BMC Bioinformatics* 2020) :

1. **Profil de pics** (`src/peaks.js`) : à chaque base, la hauteur des quatre canaux donne un pic primaire, un pic secondaire et leur rapport. Seules les extrémités de faible intensité sont rognées, pas les zones mixtes.
2. **Référence génomique** (`src/reference.js`) : les lectures sur ADN génomique couvrent exons et introns. On travaille donc sur l’enregistrement RefSeqGene (NG_), et le transcrit NM_ fournit les coordonnées `c.` (exons, CDS, introns `c.N+d` / `c.N-d`, UTR `c.-N` / `c.*N`).
3. **Alignement du profil** (`src/decompose.js`, `src/align.js`) : la lecture est localisée par k-mers, orientée, puis alignée (Gotoh, gaps affines). Une base de référence compte comme appariée dès qu’elle a un signal mesurable. L’alignement suit ainsi l’allèle de référence, même après un indel hétérozygote.
4. **Séparation des allèles** : allèle 1 = la référence quand elle est présente, allèle 2 = l’autre pic. L’allèle 2 est ensuite corrigé par un alignement qui accepte les pics faibles.
5. **Appel des variants** : chaque allèle est aligné sur la référence. Les indels sont normalisés en 3′ (règle HGVS), les duplications sont reconnues, et les SNV hétérozygotes doivent dépasser le bruit de fond local.
6. **Ré-ajustement des indels** : on teste tous les placements possibles d’un indel autour du point de rupture, et on garde l’allèle prédit qui reproduit le mieux les profils de pics.
7. **Forward + Reverse** : un variant appelé sur les deux brins est « confirmé F/R ». Un variant non retrouvé sur une lecture qui couvre pourtant la position est signalé.

## Utilisation

Les modules ES nécessitent un serveur HTTP : GitHub Pages convient, ou en local `python3 -m http.server` puis http://localhost:8000.

1. **Référence** : rechercher le gène et charger le NG_ proposé (ex. `NG_012772.3` pour BRCA2), puis choisir le transcrit (ex. `NM_000059.4`). Sans accès à NCBI, importer un fichier GenBank NG_ téléchargé au préalable (format « GenBank (full) »).
2. **Chromatogrammes** : importer les `.ab1` Forward et Reverse, puis cliquer sur « Analyser ».
3. **Variants** : cliquer sur une ligne pour voir le contexte et le chromatogramme centré sur le variant. Les lignes Réf. / Allèle 1 / Allèle 2 sous les pics montrent la déconvolution. Export CSV disponible.

Un jeu de démonstration entièrement synthétique est fourni dans `demo/` (gène fictif, `c.752dup` hétérozygote) : référence `demo/NG_DEMO.gb`, lectures `DEMO_F.ab1` et `DEMO_R.ab1`. Pour le régénérer : `node test/make-demo.js`.

## Tests

```sh
npm test                    # tests unitaires et de bout en bout (Node ≥ 20)
node test/benchmark.js 30   # banc d'essai sur chromatogrammes simulés
```

Banc d’essai actuel (30 variants aléatoires par condition : SNV, duplications, délétions 1–10 nt, insertions 1–5 nt ; lectures F et R simulées) :

| Fraction mutée | Bruit | Exact et confirmé F/R | Trouvé sur ≥ 1 lecture | Faux appels / cas (confirmés F/R) |
|---|---|---|---|---|
| 50 % | faible | 100 % | 100 % | 0,00 (0) |
| 50 % | fort | 97 % | 100 % | 0,07 (0) |
| 35 % | faible | 100 % | 100 % | 0,00 (0) |
| 35 % | fort | 93 % | 100 % | 0,10 (0) |
| 25 % | faible | 93 % | 97 % | 0,03 (0) |
| 25 % | fort | 87 % | 93 % | 0,17 (0) |

Ces chiffres portent sur des **données simulées**. La validation sur de vrais chromatogrammes, avec variants connus et témoins négatifs, reste à faire.

## Limites connues

- Allèle minoritaire < 20–25 % (mosaïque) : sensibilité réduite.
- Indels hétérozygotes composites (deux allèles mutés différents) : non résolus.
- Extrémités de lecture (~20 bases) : pas d’appel d’indel, car le signal en aval est insuffisant pour le confirmer.
- Transcrits sur le brin − de leur NG_ : nomenclature `g.` seulement (cas exceptionnel, les NG_ étant orientés dans le sens du gène).
- Pas encore de conséquence protéique (`p.`).

## Structure

```
index.html, style.css     interface
src/app.js                logique de l'interface
src/trace-view.js         affichage du chromatogramme (canvas)
src/abif.js               lecteur ABIF (.ab1)
src/peaks.js              profils de pics, rognage
src/align.js              alignement local affine, localisation par k-mers
src/decompose.js          déconvolution, appel et fusion F/R des variants
src/hgvs.js               normalisation 3′, coordonnées c., noms HGVS
src/reference.js          GenBank, RefSeqGene, NCBI E-utilities
test/                     tests, simulateur de chromatogrammes, banc d'essai
```
