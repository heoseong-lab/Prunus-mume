# Prunus mume germplasm: reliability-weighted phenotypic diversity and organic acid variation

R code for the analysis reported in:

> Phenotypic Diversity and Variation in Organic Acid Content Among 204 Japanese Apricot (*Prunus mume*) Germplasm Accessions (manuscript in preparation, 2026).

The pipeline evaluates 20 phenotypic traits scored on 204 *P. mume* accessions over three years (2024 to 2026, 484 observations) and three organic acids quantified in 88 accessions. Its central idea is that the year-to-year reproducibility of each trait is estimated first and then used as that trait's weight in a Gower distance, so that no trait or accession has to be discarded by an arbitrary cutoff.

## Pipeline

| Script | What it does | Manuscript section |
|---|---|---|
| `R/00_setup.R` | Packages, paths, trait dictionary, helper functions (Fleiss' kappa, Krippendorff's alpha, Spearman-Brown conversion, bias-corrected Cramér's V, adjusted eta-squared, Shannon index, adjusted Rand index, plotting theme) | 2.4 to 2.11 |
| `R/01_preprocess.R` | Data curation (rare-level merging, unit correction of 2024 fruit weight, outlier listing), linear mixed models with year, accession and harvest batch, repeatability, BLUPs, prediction-error-variance reliability, accession-level matrix | 2.4 to 2.6 |
| `R/02_diversity.R` | Descriptive statistics, Shannon-Weaver diversity, between-year correlations, categorical agreement (complete and majority agreement, Fleiss' kappa, Krippendorff's alpha), year shift and cumulative link mixed models, bias-corrected association matrix | 2.5, 2.6, 3.1 to 3.3 |
| `R/03_multivariate.R` | Composite reliability and weights, reliability-weighted Gower distance on sqrt(1 - s), PCoA and FAMD, Ward.D2 clustering with silhouette scan, eight sensitivity configurations with ARI, nominal-coding permutations, core collection by maximum dissimilarity with MD/VD/CR/VR | 2.7 to 2.9, 3.4 to 3.6 |
| `R/04_organic_acid.R` | Organic acid ANOVA (Bartlett, Welch), five-class grading against LSD, Tukey HSD for extremes, compositional types, correlations among the organic acids on absolute contents, selection lists | 2.10, 3.7, 3.8 |
| `R/run_all.R` | Runs the five scripts in order | |

`mermaid/Fig2_analysis_pipeline_EN.mermaid` is the source of Figure 2 (analytical pipeline). `R/README_pipeline_ko.md` is a detailed Korean description of the pipeline, the analytical decisions and their rationale.

## Requirements

R 4.3 or later. Packages: tidyverse, lme4, car, agricolae, FactoMineR, factoextra, cluster, ggrepel, ordinal (optional; needed for the cumulative link mixed models), clustMixType and ggtern (optional).

```r
install.packages(c("tidyverse", "lme4", "car", "agricolae", "FactoMineR",
                   "factoextra", "cluster", "ggrepel", "ordinal"))
```

## How to run

1. Place the two input files in the project root (see Data below).
2. Set the working directory to the project root.
3. Run:

```r
source("R/run_all.R", encoding = "UTF-8")
```

or from a terminal:

```bash
Rscript R/run_all.R
```

The full run takes about 20 seconds on a laptop and writes tables to `output/tables/` and figures to `output/figures/`. The random seed is fixed at 2026, and the label placement in the PCoA plot is made deterministic (`seed`, `max.time = Inf`), so repeated runs give identical output.

## Data

The input files are not part of this repository:

- `매실 전체 데이터.csv`: phenotypic records (one row per accession and year; 20 traits).
- `매실_유기산_자원.csv`: organic acid means and standard deviations of three replicates per accession, recorded in mg per 100 g fresh weight. The scripts convert them to mg per g.

Availability of the data is described in the Data Availability Statement of the manuscript.

## Notes on a few analytical decisions

- Gower's coefficient is used as sqrt(1 - s), which is the Euclidean-embeddable form (Gower 1971); with 1 - s, 75% of the PCoA eigenvalues were negative. Ward.D2 squares its input, so it receives sqrt(1 - s), equivalent to ward.D on 1 - s.
- Composite reliability is the Spearman-Brown-converted repeatability for continuous traits and Fleiss' kappa for categorical traits; Krippendorff's alpha and the PEV-based BLUP reliability are reported as comparison indices, and alpha weighting is kept as sensitivity configuration 8 because the 0/1 partial distance of nominal traits lets a single nominal trait dominate the clustering when its weight rises.
- Harvest date is a bulk harvest date and enters only as a harvest-batch random effect, not as an accession trait.

## Citation

Please cite the manuscript above once it is published. A citation entry will be added here.
