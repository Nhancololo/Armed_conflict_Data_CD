# data

Place the input files here. The scripts read and write with `here::here()`,
so nothing depends on the working directory.

- `raw/acled/`: ACLED records written by `scripts/data/01_download_acled.R`.
- `raw/covariates/`: downloaded rasters, vector layers and tables of the
  covariate scripts (`05*`).
- `processed/`: panel, boundaries, adjacency and covariate tables.
- `cabo_delgado_data.rds`: analysis object read by the modelling script,
  written by `scripts/data/06_build_bundle.R`.

ACLED records cannot be redistributed. Obtain them with an own ACLED account.
