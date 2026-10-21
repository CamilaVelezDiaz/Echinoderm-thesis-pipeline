# Taxonomic Assessment and Biodiversity Analysis of Historical Deep-Sea and Cross-Shelf Echinoderm Collections from northeast Australia and the Adjacent Coral Sea

R workflow for the integrated echinoderm occurrence dataset for northeast Australia and the adjacent Coral Sea. The repository has two stages:

1. **Dataset construction (scripts 01–02):** produces a curated, quality-controlled, Darwin Core Archive of 43,470 echinoderm occurrence records from 13 source datasets, described in the Biodiversity Data Journal Data Descriptor.
2. **Biodiversity analyses (seven scripts):** taxonomic completeness, identification and representativeness indices, richness estimation, distribution maps, temporal trends and specimen provenance, reported in the associated MMB (Research) thesis at James Cook University and manuscript in the x.

## Contents

### Stage 1 – Dataset construction

| File | Purpose |
|------|---------|
| `01_echinoderm_dataset_integration.R` | Loads 13 source datasets, filters to Echinodermata within the study bounding box, deduplicates across sources, writes long-format (Layer A) and wide-format (Layer B) tables. |
| `02_echinoderm_post_processing.R` | Resolves cross-source conflicts (coordinates, year, basis of record, depth, taxonomic name against WoRMS), assigns depth zones, writes `echino_wide.csv` and the Darwin Core Archive. |

### Stage 2 – Analyses (all read `echino_wide.csv`)

| File | Purpose | Main outputs |
|------|---------|--------------|
| `Taxonomic Completeness.R` | Taxonomic resolution (rank reached) across depth, class and source: [ordinal models / GAMs / cutpoints]. | [tables and figures] |
| `Taxonomic Identification Index (TIa).R` | Taxonomic Identification Index (TIa) by depth zone and class. | [ ] |
| `Taxonomic Representative Index (TRi).R` | TRi = CSi / CSe × 100, where CSe is the expected species list built from the Australian Faunal Directory (IMCRA regions, WoRMS-resolved). Option `TRI_SCOPE` runs all records or Australia only. | [ ] |
| `Taxonomic Representative Index-Australia RECORDS ONLY.R` | Option `TRI_SCOPE` runs all records or Australia only. | [ ] |
| `Rarefaction and Richness Estimators.R` | Incidence-based rarefaction and Chao2 using iNEXT with 0.5° grid cells as sampling units, for Shelf, Upper Slope and Lower Slope; the Abyssal zone is descriptive only. Option `RICH_SCOPE` runs all records or Australia only. | [ ] |
| `Rarefaction and Richness Estimators-Australia RECORDS ONLY.R` | Option `RICH_SCOPE` runs all records or Australia only. | [ ] |
| `Distribution Maps.R` | Record maps, sampling density and species richness grids. | [ ] |
| `Temporal Trends.R` | Records, new species and TIa by decade; cumulative species discovery. | [ ] |
| `Species Provenance.R` | Where specimen-based records (obs_quality_flag = 1) are held, by country, institution and decade. | [ ] |

Suggested run order: 01 → 02 → any analysis script (they are independent of each other).

## Requirements

- **R version**: 4.5.1 or higher.
- **Stage 1** (`01`, `02`): `dplyr`, `tibble`, `stringr`, `readr`, `purrr`, `tidyr`, `worrms` (WoRMS API), `sf` and `rnaturalearth` (land QC flag in Section 1d of `02`). The OBIS source is a static snapshot generated once with `robis::occurrence()`, so `robis` is not required to run the pipeline.
- **Stage 2**: `tidyverse`, `iNEXT`, `mgcv`, `lme4`, `ordinal`, `sf`, `rnaturalearth`, `ggplot2` [add any others from your scripts: run `renv::snapshot()` to list them].
- Full package versions are in `sessionInfo.txt` (and `renv.lock`).

## How to run

1. Clone the repository and download the source datasets into `data-raw/` (see Data sources).
2. Optional: `Rscript -e 'renv::restore()'`
3. Build the dataset:
```bash
   Rscript "01_echinoderm_dataset_integration.R"
   Rscript "02_echinoderm_post_processing.R"
```
4. Run any analysis script from the repository root. Each expects `echino_wide.csv` in the working directory and writes its outputs to its own folder (e.g. `TRi_all/`, `richness_incidence_all/`).

## Analysis notes

- Random seeds are fixed (`set.seed(42)`) and bootstrap runs use `nboot = 999`.
- Record counts reflect collecting and digitisation effort, not individuals. Richness is therefore estimated with incidence data (0.5° cells); record-based curves are exploratory only.
- The Abyssal zone (3 species-level records from 2 cells) is too sparse for estimation and is reported descriptively.
- A sensitivity analysis excluding Papua New Guinea and New Caledonia records is available through `TRI_SCOPE` / `RICH_SCOPE`.
