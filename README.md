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
6. **Recherche directe du décalage** : au point où le signal devient mixte (ou dès le début de la lecture), toutes les délétions (1–120 nt) et insertions (1–40 nt) sont testées ; on garde celle dont le mélange prédit reproduit le mieux les profils de pics. C’est ce qui permet d’appeler les grands indels (ex. MSH2 c.212-25_239del, 53 nt) en un seul variant.
7. **Ré-ajustement local** : les autres groupes d’appels sont re-placés ou écartés de la même façon, et les SNV de l’allèle 2 ne sont gardés que s’ils améliorent l’ajustement.
8. **Forward + Reverse** : un variant appelé sur les deux brins est « confirmé F/R ». Une lecture qui ne l’a pas appelé — par exemple parce qu’elle commence déjà dans la zone décalée — le confirme aussi si son signal est nettement mieux expliqué avec lui (« par modélisation »). Quand deux lectures proposent des placements voisins d’un même indel, on retient celui qui explique le mieux l’ensemble des lectures.

## Utilisation

Les modules ES nécessitent un serveur HTTP : GitHub Pages convient, ou en local `python3 -m http.server` puis http://localhost:8000.

1. **Référence** : taper le symbole du gène (ex. `BRCA2`) puis cliquer sur « Charger ». L’outil trouve le RefSeqGene du gène (NG_012772.3 pour BRCA2) et présélectionne le transcrit de référence (MANE Select s’il est signalé, sinon le premier décrit). Les références déjà utilisées sont mémorisées dans le navigateur (IndexedDB) et se rechargent sans NCBI. Une accession NG_ précise ou un fichier GenBank / FASTA local restent possibles.
2. **Chromatogrammes** : importer les `.ab1` Forward et Reverse, puis cliquer sur « Analyser ».
3. **Variants** : cliquer sur une ligne pour voir le contexte et le chromatogramme centré sur le variant. Les lignes Réf. / Allèle 1 / Allèle 2 sous les pics montrent la déconvolution. Export CSV disponible.

Un jeu de démonstration entièrement synthétique est fourni dans `demo/` (gène fictif, `c.752dup` hétérozygote) : référence `demo/NG_DEMO.gb`, lectures `DEMO_F.ab1` et `DEMO_R.ab1`. Pour le régénérer : `node test/make-demo.js`.

## Tests

```sh
npm test                    # tests unitaires et de bout en bout (Node ≥ 20)
node test/benchmark.js 30   # banc d'essai sur chromatogrammes simulés
```

Banc d’essai actuel (30 variants aléatoires par condition ; lectures F et R simulées) :

| Fraction mutée | Bruit | Petits variants¹ : exact et confirmé F/R | Avec grands indels² : exact et confirmé F/R | Faux appels confirmés F/R |
|---|---|---|---|---|
| 50 % | faible | 100 % | 100 % | 0 |
| 50 % | fort | 100 % | 100 % | 0 |
| 35 % | faible | 100 % | 100 % | 0 |
| 35 % | fort | 100 % | 100 % | 0 |
| 25 % | faible | 97 % | 100 % | 0 |
| 25 % | fort | 90 % | 93 % | ≤ 0,03 / cas |

¹ SNV, duplications 1–3 nt, délétions 1–10 nt, insertions 1–5 nt (`node test/benchmark.js 30`).
² mêmes types plus délétions 11–80 nt et insertions 6–30 nt (`node test/benchmark.js 30 --grands`).

Premier cas réel vérifié : MSH2 c.212-25_239del (délétion hétérozygote de 53 nt englobant un poly-T), appelée exactement et confirmée F/R alors que la lecture Forward commence dans la zone décalée.

Ces chiffres portent sur des **données simulées**. La validation sur de vrais chromatogrammes, avec variants connus et témoins négatifs, reste à faire.

## Limites connues

- Allèle minoritaire < 20–25 % (mosaïque) : sensibilité réduite.
- Indels hétérozygotes composites (deux allèles mutés différents) : non résolus.
- Extrémités de lecture (~20 bases) : pas d’appel d’indel, car le signal en aval est insuffisant pour le confirmer.
- Indels > 120 nt (délétion) ou > 40 nt (insertion) : non recherchés directement.
- L’enregistrement NG_ peut décrire une version antérieure du transcrit (ex. NG_007110.2 → NM_000251.2) : la numérotation c. du CDS est en général identique, mais à vérifier pour les UTR.
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
src/cache.js              mémoire locale des références (IndexedDB)
test/                     tests, simulateur de chromatogrammes, banc d'essai
```
