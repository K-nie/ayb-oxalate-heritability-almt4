# Marker-based heritability and the aluminum-activated malate transporter ALMT4 anchor soluble-oxalate variation in African yam bean (Sphenostylis stenocarpa)

Analysis code and result tables/figures for marker-based heritability and candidate-gene work on 13 seed-quality and anti-nutritional traits in the African yam bean DArTseq panel. Covers phenotype exploration, trait PCA, single- and multi-trait GBLUP, bootstrap and Bayesian heritability, mixed-model and FarmCPU/mrMLM GWAS, candidate-gene callouts against the Funannotate annotation, the ALMT4 protein-topology and orthology analysis, and breeding-decision tools (GEBV, optimal contribution selection, usefulness criterion, selection index).

This repository contains the **analysis code and derived result tables/figures** for
the study. Raw genotype and phenotype data are archived separately (see Data below);
manuscript drafts are not included.

## Repository layout

```
scripts/    numbered Python and R analysis scripts (shared helpers: _plotstyle, _pheno, _figstyle)
results/    one directory per analysis stage, each with figures/, tables/, and a README.md
            documenting method, inputs, outputs, findings, and caveats
refs/       machine-learning best-practice reference material
```

Each `results/<NN>_*/README.md` is the reproducibility and methods record for that
stage: what it does and why, its inputs, the exact commands and thresholds, the outputs,
and the caveats.

## Reproducing

Scripts run from a single conda environment (`yeast-viz`: numpy, pandas, scikit-learn,
matplotlib, plus the R toolchain for the `.R` steps). Stages are numbered by pipeline order;
run the lower-numbered QC/PCA stage first, then the analysis of interest. Per-stage READMEs
give the exact command and parameters.

## Data

- **Raw DArTseq genotypes** (order DAf18-2580, 105 accessions; post-QC 95 accessions x 1,625
  markers at sample and marker call rate >= 0.90, MAF >= 0.05): Zenodo [10.5281/zenodo.20348832](https://doi.org/10.5281/zenodo.20348832).
- **Reference genome**: *S. stenocarpa* chromosome-scale assembly, ENA [PRJEB57813](https://www.ebi.ac.uk/ena/browser/view/PRJEB57813)
  (Shorinola et al. 2024); Funannotate annotation at Zenodo [10.5281/zenodo.13853757](https://doi.org/10.5281/zenodo.13853757).
- Phenotype data are held by the breeding program and available from the authors on
  reasonable request.

## Headline findings

- Soluble-oxalate variation anchors to the aluminum-activated malate transporter ALMT4 near Ss10:15.4 Mb by proximity to a mixed-model GWAS signal.
- Genomic heritability is estimated per trait with bootstrap and Bayesian intervals; the profile-likelihood interval excludes zero for soluble oxalate.
- Multi-trait joint GBLUP improves predictive ability over single-trait models for the correlated oxalate fractions.
- Breeding-decision layers (optimal contribution selection, usefulness criterion, selection index) are computed on the panel for downstream mate allocation.

## Citation

If you use this code or its outputs, please cite this repository (see `CITATION.cff`)
and the accompanying manuscript.

## License

Code is released under the MIT License (see `LICENSE`).
