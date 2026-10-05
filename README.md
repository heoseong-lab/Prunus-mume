# Reproducibility-weighted phenotypic characterization of *Prunus mume* germplasm

R code for the analysis reported in:

> Jang, H.S.; Heo, S. Reproducibility-Weighted Phenotypic Characterization of Japanese Apricot (*Prunus mume*) Germplasm from Multi-Year Field Records. *???* (under review, 2026).

The pipeline characterizes 202 *P. mume* accessions from 480 accession-year records (2024 to 2026) with a two-stage procedure. In the first stage, quantitative traits are adjusted within each year (fruit traits by harvest batch, the other traits by year). In the second stage, every variable recorded each year is weighted in a Gower distance by its between-year reproducibility: the mean between-year Spearman correlation of the adjusted values for quantitative traits and Fleiss' kappa for categorical descriptors. The weighted distance is clustered (Ward's minimum-variance method), validated by leave-one-year-out comparisons and disjoint splits of the collection, and used together with the descriptor classes to select a core collection with ShinyCore.

## Pipeline

| Script | What it does | Manuscript section |
|---|---|---|
| `R/00_config.R` | Packages, constants, trait dictionary, curation settings, shared functions (Fleiss' kappa, Krippendorff's alpha, adjusted Rand index, weighted Gower distance, ShinyCore selection) | 2.2 to 2.7 |
| `R/01_prepare.R` | Curation, exclusions, reference-year harmonization of the two varietal descriptors, stage-1 within-year adjustment, summary tables | 2.2, 2.3 |
| `R/02_reproducibility.R` | Accession profiles, between-year reproducibility coefficients and weights, bootstrap intervals, year effects on ordinal descriptors (cumulative link mixed models) | 2.3, 2.4, 3.1, 3.2 |
| `R/03_structure.R` | Weighted Gower distance, Lingoes correction, PCoA, Ward clustering with silhouette scan, agreement with descriptors, sensitivity configurations | 2.5, 2.6, 3.3, 3.4 |
| `R/04_validation.R` | Disjoint splits (weight transfer) and leave-one-year-out comparisons | 2.6, 3.4 |
| `R/05_priority.R` | ShinyCore core collection, maximum-dissimilarity and random benchmark sets, sensitivity to class number and coverage targets | 2.7, 3.5 |
| `R/06_figures.R` | Figures 1 to 5 | |
| `R/07_export.R` | Collects the numbers quoted in the manuscript in `output/manuscript_numbers.json` | |
| `R/run_all.R` | Runs the scripts in order and writes `output/sessionInfo.txt` | |

## Requirements

R 4.3 or later. Packages: cluster, ggplot2 (3.5 or later for the legend placement of Figure 2; older versions use a fallback), dplyr, tidyr, readr, jsonlite, ordinal.

```r
install.packages(c("cluster", "ggplot2", "dplyr", "tidyr", "readr", "jsonlite", "ordinal"))
```

## How to run

1. Place the input file `data/phenotypes.csv` in the project (see Data below).
2. From the project root:

```r
source("R/run_all.R", encoding = "UTF-8")
```

or from a terminal:

```bash
Rscript R/run_all.R
```

The full run takes about 2 minutes and writes tables to `output/tables/`, figures to `output/figures/`, and the manuscript numbers to `output/manuscript_numbers.json`. All resampling uses fixed seeds (`SEED = 20261001`), so repeated runs give identical output. The default sizes are 1000 accession bootstraps for the reproducibility coefficients, 100 disjoint splits, and 1000 random benchmark sets; the environment variables `PM_BOOT`, `PM_SPLIT`, and `PM_RANDOM` override them for quick checks (for example `PM_BOOT=200 PM_SPLIT=20 Rscript R/run_all.R`).

## Data

The phenotypic records are not distributed with the code (see the Data Availability Statement of the manuscript). The scripts expect `data/phenotypes.csv`, one row per accession and year, with the columns listed in `data/README.md`.

## Trait dictionary

Seventeen variables enter the analysis.

| Variable | Type | Levels or unit | Role |
|---|---|---|---|
| `fruit_weight` | quantitative | g per fruit | adjusted by harvest batch; weighted |
| `soluble_solids_content` | quantitative | °Brix | adjusted by harvest batch; weighted |
| `titratable_acidity` | quantitative | % | adjusted by harvest batch; weighted |
| `petal_number` | quantitative | count per flower | adjusted by year; weighted; mean ≥ 10 defines double flowers |
| `bloom_doy` | quantitative | day of year, from `full_bloom_date` | adjusted by year; weighted |
| `ground_color` | nominal | as recorded | weighted (Fleiss' kappa) |
| `leaf_base_shape` | nominal | acute, obtuse, rounded | weighted (Fleiss' kappa) |
| `fruit_shape` | nominal | as recorded | varietal characteristic; weight 0; describes the groups |
| `flesh_color` | nominal | as recorded | varietal characteristic; weight 0; describes the groups |
| `skin_color_intensity` | ordinal | absent or very weak, weak, medium, strong, very strong | weighted (Fleiss' kappa) |
| `stone_adherence` | ordinal | clingstone, semi-clingstone, semi-freestone, freestone | weighted (Fleiss' kappa) |
| `leaf_green_intensity` | ordinal | light, medium, dark | weighted (Fleiss' kappa) |
| `leaf_pubescence` | ordinal | absent, weak, medium, strong | weighted (Fleiss' kappa) |
| `leaf_angle` | ordinal | acute, right-angled, moderately obtuse, strongly obtuse | weighted (Fleiss' kappa) |
| `leaf_length` | ordinal | short, medium, long | weighted (Fleiss' kappa) |
| `leaf_width` | ordinal | narrow, medium, broad | weighted (Fleiss' kappa) |
| `flower_diameter` | ordinal | small, medium, large | weighted (Fleiss' kappa) |

Ordinal descriptors are scored on an equally spaced 1 to 9 scale for the distance. The weight of a variable is max(coefficient, 0); negative coefficients give weight 0.

## Data-curation rules

These rules are applied in `R/01_prepare.R` and documented in the output tables; the input file itself is never modified.

- The disease descriptors `brown_rot_severity` and `shot_hole_severity` and the accessions KJM261 and KJM264 are excluded.
- `fruit_weight` recorded in 2024 is divided by 3 (three-fruit totals) to give the per-fruit value.
- `skin_color_intensity` spellings `stong` and `very stong` are corrected; `leaf_green_intensity` recorded as `other` is set to missing.
- `leaf_base_shape` recorded as `truncate` or `cordate` is scored as `rounded`, following the three states of the UPOV test guidelines (acute, obtuse, rounded).
- Full bloom dates later than May and harvest dates later than August, or dates that do not fall in the record's year, are set to missing (`date_curation_log.csv`).
- `fruit_shape` and `flesh_color` are varietal characteristics with one value per accession; the recorded value of the reference year (2024 for fruit shape, 2026 for flesh color, otherwise the nearest listed year) is used for all years (`harmonization_log.csv`). They receive weight 0 in the distance and are used to describe the groups and to select the core collection.
- Harvest batches are defined by year and harvest date; fruit traits are centered on the batch median, petal number and bloom date on the annual median, and the overall median is added back so that the original units are kept.

## Core collection

`R/00_config.R` contains a plain-R port of the covering and thickening algorithm of ShinyCore (Kim et al. 2023, *Plant Methods* 19:106; https://github.com/heoseong/ShinyCore), applied to the descriptor classes (categorical levels, and quantitative profiles cut into 5 equal-width classes). Coverage targets are the ShinyCore defaults (maximum 99%, minimum 98%); the environment variables `PM_COV_MAX`, `PM_COV_MIN`, and `PM_CLASSES` change the targets and the class number.

## Citation

Please cite the manuscript above once it is published. A citation entry will be added here.
