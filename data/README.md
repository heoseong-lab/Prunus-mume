# Input data

The phenotypic records are not distributed with the code; see the Data Availability Statement of the manuscript. The scripts read `data/phenotypes.csv` (UTF-8, one row per accession and year) with the following columns:

| Column | Content |
|---|---|
| `no` | accession identifier (e.g., KJM001) |
| `year` | evaluation year (2024, 2025, 2026) |
| `full_bloom_date`, `maturity_date` | dates (YYYY-MM-DD); the maturity date is the harvest date that defines the harvest batch |
| `fruit_shape`, `ground_color`, `flesh_color`, `leaf_base_shape` | nominal descriptors (level names in English) |
| `skin_color_intensity`, `stone_adherence`, `leaf_green_intensity`, `leaf_pubescence`, `leaf_angle`, `leaf_length`, `leaf_width`, `flower_diameter` | ordinal descriptors (level names as listed in the trait dictionary of the main README) |
| `fruit_weight` | g (2024 values are three-fruit totals and are divided by 3 by the scripts) |
| `petal_number` | petals per flower |
| `soluble_solids_content` | °Brix |
| `titratable_acidity` | % |
| `brown_rot_severity`, `shot_hole_severity` | disease descriptors, excluded from the analysis |

Empty cells and `NA` are read as missing. `R/01_prepare.R` writes `data/phenotypes_harmonized.csv` (the records after reference-year harmonization of fruit shape and flesh color); it is ignored by git together with every other `.csv` file.
