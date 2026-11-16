# Download the 2017 road network from OpenStreetMap.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(sf, dplyr, tibble, readr, curl, here)

url <- "https://download.geofabrik.de/africa/mozambique-170101-free.shp.zip"
boundary_cfg <- here::here("data", "processed",
  "cabo_delgado_admin2_canonical.rds")
raw_dir <- here::here("data", "raw", "covariates", "osm_geofabrik_2017")
zip <- file.path(raw_dir, "mozambique-170101-free.shp.zip")
unzip_dir <- file.path(raw_dir, "unzipped")
out <- file.path(raw_dir, "roads_osm_2017_cabo_delgado.rds")
QA <- here::here("outputs", "covariates", "05A2_osm2017")
for (d in c(raw_dir, unzip_dir, QA)) dir.create(d, recursive = TRUE,
  showWarnings = FALSE)
if (!file.exists(boundary_cfg)) stop("Missing canonical boundary.")

if (!file.exists(zip) || file.info(zip)$size < 1e6) {
  message("Downloading historical Geofabrik Mozambique snapshot (~250 MB)...")
  curl::curl_download(url, zip, mode = "wb", quiet = FALSE)
}
if (!length(list.files(unzip_dir, pattern = "\\.shp$", recursive = TRUE))) {
  utils::unzip(zip, exdir = unzip_dir)
}
shps <- list.files(unzip_dir, pattern = "\\.shp$", recursive = TRUE,
  full.names = TRUE)
road_shp <- shps[grepl("roads", basename(shps), ignore.case = TRUE)]
if (length(road_shp) != 1L) {
  readr::write_csv(tibble(file = shps), file.path(QA,
    "ERROR_shapefile_inventory.csv"), na = "")
  stop(
    "Could not uniquely identify Geofabrik roads shapefile. Inspect ",
    "ERROR_shapefile_inventory.csv")
}

roads <- sf::st_read(road_shp, quiet = TRUE,
  stringsAsFactors = FALSE) |> sf::st_make_valid()
d <- readRDS(boundary_cfg) |> sf::st_make_valid()
metric_crs <- 32737
extent <- sf::st_buffer(sf::st_union(sf::st_transform(d, metric_crs)),
  50000) |> sf::st_transform(sf::st_crs(roads))
hit <- lengths(sf::st_intersects(roads, extent)) > 0L
roads <- roads[hit, , drop = FALSE] |> sf::st_transform(4326)
saveRDS(roads, out, compress = "xz")

readr::write_csv(
  tibble(
    source_url = url,
    snapshot = "2017-01-01",
    source_zip = zip,
    roads_shapefile = road_shp,
    n_features = nrow(roads),
    md5_zip = unname(tools::md5sum(zip)),
    output = out
  ),
  file.path(QA, "00_osm2017_manifest.csv"), na = ""
)
if ("fclass" %in% names(roads)) {
  readr::write_csv(sf::st_drop_geometry(roads) |> count(fclass, sort = TRUE),
    file.path(QA, "01_road_class_counts.csv"), na = "")
}
message("Historical OSM 2017 road baseline written: ", out)
