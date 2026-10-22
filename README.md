Code accompanying:

> Potgieter LJ, Gildenhuys C, Gaertner M, Richardson DM. Global invasion patterns and niche dynamics in the shrub genus Cotoneaster (Rosaceae). *Biological Invasions* (submitted 2026).

------------------------------------------------------------------------

# Range size analysis

## Methods

Cotoneaster records were downloaded from the Global Biodiversity Information Facility (GBIF) on 15 May 2026 (<https://doi.org/10.15468/dl.gscv7g>). Data were loaded into and all subsequent analyses were conducted in R version 4.6.0.

Known non-native ranges were obtained from the Global Naturalised Alien Flora (GloNAF) database version 3 available at <doi:10.5281/zenodo.17105725> [@davis]. Known native ranges were obtained from the World Checklist of Vascular Plants (WCVP) [@wcvp:wo2024] via the `rWCVPdata` package [@brown2023].

All range information was standardised according to the World Geographic Scheme for Recording Plant Distributions (WGSRPD) level 3, the "Botanical Country" level [@brummitt2001]. The WGSRPD attempts to supply a standardized set of shape files generally following administrative boundaries to designate global regions for for the classification of range size...

The GBIF records were imported into R and cleaned using the `coordinateCleaner` package [@zizka2019]. Records were tested for zeros, seas, institutions, gbif, equal, centroids, and capitals. All flagged records were removed.

All Cotoneaster species names (as recorded in GBIF) were matched with the WCVP to standardise nomenclature between the different databases. Where an exact match could not be found, fuzzy matching was employed and the match edit distance (the number of edits that had to be made to the original name to find a match in the WCVP) and the match similarity (the edit distance divided by the number of characters in the name) were calculated. Results of the matching were examined and all species names with a match similarity of less than 0.9 were discarded as these were deemed unreliable. All cases of multiple matches (i.e. when multiple candidate names were found to a single species name in GBIF) were resolved algorithmically.

GBIF records were spatially joined with the GloNAF and WCVP ranges using the `sf` package. All records were assigned as "native" (when the record falls within the native range boundaries as recorded by the WCVP), "non-native" (when the record falls within the introduced ranges according to WCVP or in the naturalised or invasive ranges according to GloNAF), "disputed" when it is recorded as native by the WCVP and non-native by GloNAF within the same boundary) and finally NA (when the particular species is not recorded as either "native" or "non-native" in the boundary).

Many of these NA records occur because these databases are are incomplete or imperfect (i.e. species introduced into an area but not get recorded as introduced or naturalised or a species occurs just outside the biogeographically arbitrary boundaries). To address these data gaps, distance to the nearest shapefile boundary and identity of said boundary (i.e. "native" or "non-native") were calculated for all NA variables.

If the NA record occurred within 100km of the nearest boundary it was assigned the nativeness of that boundary. If the NA record occurred further than 100km away from the nearest boundary it was assumed to have "non-native" status. All species distributions were mapped to check the the assigned nativeness.

Only two species had disputed nativeness, these were inspected and resolved as either "native" or "non-native". We do not claim that the assigned nativeness following this process is necessarily correct but rather a near enough approximation given gaps in the available databases and given certain assumptions for the purposes of this study. With this final dataset we mapped the total number of species, number of native species, and number of non-native species in each botanical country.

Area of occupancy was derived by overlaying a grid of one degree grid cells and calculating the total area (km<sup>2</sup>) for each native and non-native range of each species. The area was normalised by dividing by the lowest area. Area was logged for the purposes of the modelling. When a species did not have a non-native range, the area of occupancy was set to zero.

We fit a hurdle model modelling non-native area of occupancy against the native area of occupancy. Zeros were explicitly modelled and assumed to represent.

# Residence time analysis

## Cotoneaster residence time and non-native range size

Companion code for the analysis of residence time and spread in the global review of *Cotoneaster* invasions (manuscript in preparation for *Biological Invasions*).

`residence_time_analysis.R` relates the non-native range size of naturalised *Cotoneaster* species to their minimum residence time, controlling for recording effort. It reproduces the Spearman correlation, the log-linear regression, the partial correlation, the combined model, the occurrence-date control, and Figure 3 reported in the "Residence time and spread" section.

## What the script does

1.  Reads the Seebens et al. (2017) first-records database and, for each *Cotoneaster* species, takes the earliest documented first-record year (excluding pre-Linnaean placeholder dates of 1492 and 1500). Minimum residence time is the reference year (2026) minus that year.
2.  Reads the cleaned GBIF occurrence table and, for each species, counts the distinct WGSRPD level-3 units in which it is recorded as non-native (the range measure) and the number of non-native occurrence records.
3.  Joins the two, retains species with at least 20 non-native records, and fits the statistics and figure.

## Inputs

Neither dataset is redistributed in this repository. Download each from its original provider and point the paths in the `Configuration` block at your local copies.

| File | Source |
|----|----|
| `FirstRecords_dataset_public_v4.0.csv` | Seebens, H. et al. (2017) *Nature Communications* 8, 14435. <https://doi.org/10.1038/ncomms14435> |
| Cleaned GBIF table | GBIF download, genus *Cotoneaster*, 20 July 2026, <https://doi.org/10.15468/dl.gscv7g>, cleaned and nativeness-assigned by the main data-processing script in this repository |

The cleaned GBIF table must contain the columns `wcvp_accepted_name`, `nativeness`, `nn_nativeness`, `tdwg3` and `year`.

## Running it

``` r
# from the repository root, with the two inputs placed under data/
source("residence_time_analysis.R")
```

Missing packages are installed automatically. Outputs are written to `output/`:

-   `figure3_residence_time.png` / `.pdf`
-   `residence_time_data.csv` — per-species residence time, range, records and residuals
-   `residence_time_statistics.txt` — the reported test statistics
-   `sessionInfo.txt` — the R and package versions used

## Parameters

Set in the `Configuration` block at the top of the script:

-   `year_now` reference year for minimum residence time (default 2026)
-   `year_floor` lower bound excluding placeholder first-record dates (default 1500)
-   `min_records` minimum non-native occurrence records to retain a species (default 20)

## Dependencies

`readr`, `dplyr`, `stringr`, `tidyr`, `ggplot2`, `ggrepel`. Written and tested under R 4.3 and later.

## Citation

If you use this code, please cite the manuscript (details to be added on acceptance) and the two data sources above.

## Licence

See `LICENSE`. Code is released under the MIT Licence unless stated otherwise; the input datasets remain under the licences of their original providers.

# Cotoneaster niche comparison and model transferability

## Overview

This repository contains the R code used for the climatic niche comparison and species distribution modelling (SDM) analyses described in Aims 3 and 4 of the manuscript. The script runs end-to-end from a cleaned GBIF occurrence dataset and produces all niche overlap metrics, model evaluation statistics, MESS extrapolation diagnostics, figures, and supplementary tables reported in the paper.

The three focal species analysed are *Cotoneaster acutifolius*, *C. divaricatus*, and *C. horizontalis*. These were selected because they have sufficient native-range occurrence records for reliable model calibration; other taxa in the genus were excluded on this basis.

------------------------------------------------------------------------

## Analyses

| Section | Analysis |
|----|----|
| 1 | Spatial thinning (equal-area grid, 5 km) + sensitivity (1/2/5/10 km) |
| 2 | CHELSA V2.1 predictor loading, land masking, VIF selection |
| 3 | Accessible area delineation (200 km buffer) and background sampling |
| 4 | COUE climatic niche comparison (overlap, equivalency, similarity) |
| 5 | Native-range maxnet SDMs, global projection, P10 threshold |
| 6 | Spatially blocked cross-validation (blockCV), geographic transfer evaluation |
| 7 | Geographically constrained null model comparison |
| 8 | Regional breakdown of transfer performance (TDWG3 units) |
| 9 | MESS extrapolation diagnostics (aggregate and per-variable) |
| 10 | Response curves |
| 11 | Reciprocal transfer test (naturalised → native), full and sample-matched |
| 12–13 | Figures and summary tables |

------------------------------------------------------------------------

## Requirements

### R packages

``` r
install.packages(c(
  "tidyverse", "sf", "terra", "ade4", "ecospat", "usdm",
  "maxnet", "blockCV", "precrec", "predicts",
  "rnaturalearth", "patchwork", "tidytext"
))

# flexsdm is used via namespace (flexsdm::calib_area) — install separately:
install.packages("flexsdm", repos = c("https://sjevelazco.r-universe.dev", getOption("repos")))
```

The code was developed under R 4.4. `ecospat` ≥ 4.0 is required for the COUE workflow used here.

### Input data

Three external data sources are required. None are included in this repository.

| File | Source | Notes |
|----|----|----|
| `gbif_final_species.rds` | Contact authors | Cleaned GBIF occurrence data produced by the main *Cotoneaster* review pipeline. Available on request. |
| CHELSA V2.1 bioclimatic layers (`bio*.tif`) | [chelsa-climate.org](https://chelsa-climate.org/bioclim/) | All 19 bio variables; the script applies land masking internally. |
| `ne_50m_land.shp` | [naturalearthdata.com](https://www.naturalearthdata.com/downloads/50m-physical-vectors/) | 50 m land polygons (`ne_50m_land.zip`). |

------------------------------------------------------------------------

## Setup

Open `niche_models.R` and set the three paths near the top of the script (Section 0):

``` r
out_main   <- "path/to/output"                        # <-- SET THIS
chelsa_dir <- "path/to/CHELSA_V2.1"                   # <-- SET THIS
land_shp   <- "path/to/ne_50m_land/ne_50m_land.shp"  # <-- SET THIS
```

`out_main` must already contain `gbif_final_species.rds`. The sub-directories `figures/` and `tables/` are created automatically.

------------------------------------------------------------------------

## Outputs

All outputs are written to `out_main/`:

-   **`figures/`** — PNG figures (COUE niche panels, binary suitability maps, bio19 reciprocal response curves, supplementary MESS maps, null distributions, capture-by-region)
-   **`tables/`** — CSV files for all supplementary tables (thinning sensitivity, PCA loadings, blocked CV, transfer evaluation, null model replicates, MESS, response curves, reciprocal test, main results table)
-   **`suitability_analysis_objects.rds`** — all intermediate R objects saved for downstream use
-   **`sessionInfo_suitability.txt`** — R session info for reproducibility

------------------------------------------------------------------------

## Key analytical decisions

A number of non-default choices were made during the analysis; each is documented in the code comments and described in the manuscript Methods section.

-   **Equal-area thinning** — thinning is performed on a cylindrical equal-area projection (`+proj=cea +lat_ts=30`) so that native (c. 30°N) and naturalised (c. 50°N) records are filtered at the same true ground distance.
-   **CHELSA land masking** — CHELSA V2.1 carries interpolated values over ocean cells; masking to land before any extraction or background sampling is critical.
-   **bio6 substitution** — minimum temperature of the coldest month (bio6) was retained in place of bio1 and bio4 on ecological grounds, not as a VIF result; this is reported as such in the Methods.
-   **Null model evaluation** — null models are evaluated identically to the fitted model (restricted to the naturalised accessible area); global evaluation inflates null performance and was avoided.
-   **Boyce index scope** — for geographic transfer, the Boyce index is computed over the naturalised accessible area only, not the full global raster.

------------------------------------------------------------------------

## Citation

If you use this code, please cite the associated manuscript:

> Potgieter LJ & Gildenhuys C (2026). A global review of *Cotoneaster* (Rosaceae) invasions: distribution, drivers, and projected spread. *Biological Invasions*. [[https://doi.org/[DOI](https://doi.org/%5BDOI](https://doi.org/%5BDOI%5D(https://doi.org/%5BDOI){.uri}){.uri} TBC]

------------------------------------------------------------------------

## License

Code released under the [MIT License](LICENSE).
