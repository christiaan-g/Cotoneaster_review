# Cotoneaster residence time and non-native range size

Companion code for the analysis of residence time and spread in the global
review of *Cotoneaster* invasions (manuscript in preparation for *Biological
Invasions*).

`residence_time_analysis.R` relates the non-native range size of naturalised
*Cotoneaster* species to their minimum residence time, controlling for
recording effort. It reproduces the Spearman correlation, the log-linear
regression, the partial correlation, the combined model, the occurrence-date
control, and Figure 3 reported in the "Residence time and spread" section.

## What the script does

1. Reads the Seebens et al. (2017) first-records database and, for each
   *Cotoneaster* species, takes the earliest documented first-record year
   (excluding pre-Linnaean placeholder dates of 1492 and 1500). Minimum
   residence time is the reference year (2026) minus that year.
2. Reads the cleaned GBIF occurrence table and, for each species, counts the
   distinct WGSRPD level-3 units in which it is recorded as non-native (the
   range measure) and the number of non-native occurrence records.
3. Joins the two, retains species with at least 20 non-native records, and
   fits the statistics and figure.

## Inputs

Neither dataset is redistributed in this repository. Download each from its
original provider and point the paths in the `Configuration` block at your
local copies.

| File | Source |
|------|--------|
| `FirstRecords_dataset_public_v4.0.csv` | Seebens, H. et al. (2017) *Nature Communications* 8, 14435. https://doi.org/10.1038/ncomms14435 |
| Cleaned GBIF table | GBIF download, genus *Cotoneaster*, 20 July 2026, https://doi.org/10.15468/dl.gscv7g, cleaned and nativeness-assigned by the main data-processing script in this repository |

The cleaned GBIF table must contain the columns `wcvp_accepted_name`,
`nativeness`, `nn_nativeness`, `tdwg3` and `year`.

## Running it

```r
# from the repository root, with the two inputs placed under data/
source("residence_time_analysis.R")
```

Missing packages are installed automatically. Outputs are written to
`output/`:

- `figure3_residence_time.png` / `.pdf`
- `residence_time_data.csv` — per-species residence time, range, records and residuals
- `residence_time_statistics.txt` — the reported test statistics
- `sessionInfo.txt` — the R and package versions used

## Parameters

Set in the `Configuration` block at the top of the script:

- `year_now` reference year for minimum residence time (default 2026)
- `year_floor` lower bound excluding placeholder first-record dates (default 1500)
- `min_records` minimum non-native occurrence records to retain a species (default 20)

## Dependencies

`readr`, `dplyr`, `stringr`, `tidyr`, `ggplot2`, `ggrepel`. Written and tested
under R 4.3 and later.

## Citation

If you use this code, please cite the manuscript (details to be added on
acceptance) and the two data sources above.

## Licence

See `LICENSE`. Code is released under the MIT Licence unless stated otherwise;
the input datasets remain under the licences of their original providers.
