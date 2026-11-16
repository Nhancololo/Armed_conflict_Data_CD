# Armed_conflict_Data_CD

Code for the Bayesian spatio-temporal modelling of fatalities in Cabo Delgado
(17 districts, January 2017 to December 2024) with eight count families and
reversible-jump allocation of covariates.

## Structure

```
Armed_conflict_Data_CD/
  scripts/
    data/        data download and preparation (01 to 06)
    analysis/    functions.R and run_models.R
  data/          inputs and prepared data, see data/README.md
  outputs/       tables, figures and model checkpoints
```

## Run order

1. `scripts/data/01_download_acled.R` to `04_build_panel.R` build the
   district-month panel from the ACLED records.
2. `scripts/data/05*` build the covariates, `05_covariate_registry.R` lists
   the candidates and `05l_merge_covariates.R` merges them.
3. `scripts/data/06_build_bundle.R` writes `data/cabo_delgado_data.rds`.
4. `scripts/analysis/run_models.R` fits the families and writes the results to
   `outputs/`.

Open `Armed_conflict_Data_CD.Rproj` first so that `here` finds the project
root. Packages are installed and loaded with `pacman::p_load()`.

## Options of `run_models.R`

Environment variables change one run without editing the script:
`FAMILIES` (comma-separated subset of families), `RBST_SEED_SHIFT` and
`RBST_OUT_SUFFIX` (independent search chain), `RBST_POOL_DIRS` (pool
searches), `RBST_LINEARIZE`, `RBST_TABLES_ONLY`, `RBST_DESCRIPTIVE`,
`RBST_SKIP_POSTERIOR`, `RBST_MAP_SCALE_MAX` and `RBST_OOS`.

## Data

ACLED records require an account and cannot be redistributed. The other
sources are public and are downloaded by the scripts (MozGIS, CHIRPS,
IBTrACS, US Census Bureau IDB, WorldPop, OpenStreetMap, IOM DTM, food
prices). `05h_earth_engine.js` runs in the Google Earth Engine code editor.
