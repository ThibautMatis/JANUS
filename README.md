<img src="img/janus-logo.png" alt="JANUS — Sanger deconvolution & variant calling" height="90">

# JANUS — *Joint Allele aNalysis Using Sanger*

**v0.8.9 · 10 November 2026**

Web tool for variant calling on Sanger chromatograms (`.ab1`), inspired by Tracy (Rausch T et al. *BMC Genomics*. 2020;21:230. doi:10.1186/s12864-020-6635-8).

**Access:** https://thibautmatis.github.io/JANUS/

> ⚠ Research tool, not validated for diagnostic use. Every variant must be confirmed on the chromatogram and by a validated method.

## Usage

1. **Reference**: click a gene in the quick-access bar (BRCA1, BRCA2, PALB2, RAD51C, RAD51D, PTEN, CDH1, MSH2, MSH6, MLH1), or type the gene symbol and click « Charger » (Load). The RefSeqGene (NG_) and the reference transcript are selected automatically, and the reference is cached in the browser.
2. **Chromatograms**: import the Forward and Reverse `.ab1` files, optionally type the **targeted variant** in c. notation (e.g. `c.1813dup`), then click « Analyser » (Analyse). The targeted variant is reported as found or not found, read by read; if it is absent, the chromatogram is centred on its position and the report concludes « Absence du variant … recherché » (absence of the targeted variant), with the corresponding capture. The readability of each read is shown, together with a **reference sequence check** (« Contrôle Séq. Réf. » column): the bases called by the sequencer are aligned to the genome independently of JANUS, and any discrepancy not explained by the retained variants is flagged.
3. **Variants**: `c.`, `p.` and `g.` nomenclature, exon or intron (conventional numbering for BRCA1: no exon 4), genotype, allele fraction, F/R confirmation. The table can be sorted by any column, and the « Fréquences gnomAD » (gnomAD frequencies) button queries gnomAD v4.
4. **Chromatogram**: click a variant to centre the view on it. Hovering over a base shows its `c.` position. Below the peaks, the **Génome** (Genome) row shows the reference sequence placed by the direct check (in red where it differs from the reference projected by JANUS on a single peak).
5. **Report**: select the variants, classify them (classes 1 to 5), check the conclusion, then click « Imprimer / PDF » (Print / PDF).

A fictitious demonstration dataset is provided in `demo/` (reference `NG_DEMO.gb`, to be loaded via « Sans connexion à NCBI » (without NCBI connection), and two reads carrying `c.752dup`).

Running locally: `python3 -m http.server`, then open http://localhost:8000.

## References

- Deconvolution method inspired by Tracy: Rausch T, Fritz MH-Y, Untergasser A, Benes V. *BMC Genomics*. 2020;21:230. doi:10.1186/s12864-020-6635-8
- Classification: Richards S et al. *Genet Med*. 2015;17:405-424 (ACMG/AMP).
- Nomenclature: den Dunnen JT et al. *Hum Mutat*. 2016;37:564-569 (HGVS).
