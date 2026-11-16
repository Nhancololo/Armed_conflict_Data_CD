# Export the district boundaries for Google Earth Engine.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(sf, dplyr, readr, here)

IN <- here::here("data", "processed", "cabo_delgado_admin2_canonical.rds")
out_dir <- here::here("data", "raw", "covariates", "gee")
out <- file.path(out_dir, "cabo_delgado_admin2_canonical.geojson")
QA <- here::here("outputs", "covariates", "05G_gee_boundary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(QA, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(IN)) stop("Missing canonical boundary: ", IN)

d <- readRDS(IN) |> sf::st_make_valid() |> sf::st_transform(4326) |> arrange(district_id) |>
  mutate(district_index = row_number()) |>
  select(district_id, canonical_name, district_index, geometry)
if (nrow(d) != 17L || n_distinct(d$district_id) != 17L) stop("Boundary is not exactly 17 districts.")
if (file.exists(out)) file.remove(out)
sf::st_write(d, out, driver = "GeoJSON", quiet = TRUE)
readr::write_csv(sf::st_drop_geometry(d), file.path(QA,
  "00_gee_boundary_keys.csv"), na = "")
readr::write_csv(tibble(file = out, md5 = unname(tools::md5sum(out))),
  file.path(QA, "01_md5.csv"), na = "")
message("GeoJSON for Earth Engine written: ", out)
message(
  "Upload it as a table asset, then set DISTRICT_ASSET in ",
  "05h_earth_engine.js.")
