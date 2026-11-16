# Build distance and terrain covariates that do not vary in time.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(
  sf, dplyr, tidyr, tibble, readr, stringr, stringi, purrr, lubridate,
  units, terra, geodata, rnaturalearth, rnaturalearthdata, here
)

crs_m <- 32737
start_cfg <- as.Date("2017-01-01")
end_cfg <- as.Date("2024-12-01")
# Radius units are not guessed. Set manually to "m" or "km" only after checking
# the INCM layer metadata/attributes. NA leaves cell_coverage_share unavailable.
cell_tower_radius_unit_cfg <- NA_character_

# The primary road exposure is intentionally anchored at the beginning of the
# study period. Set false only for an explicit sensitivity run using current
# MozGIS/ANE roads.
require_osm2017_roads_cfg <- TRUE

boundary_file <- here::here("data", "processed",
  "cabo_delgado_admin2_canonical.rds")
adj_file <- here::here("data", "processed",
  "cabo_delgado_adjacency_topological.rds")
moz_dir <- here::here("data", "raw", "covariates", "mozgis")
osm2017 <- here::here("data", "raw", "covariates", "osm_geofabrik_2017",
  "roads_osm_2017_cabo_delgado.rds")
terrain_dir <- here::here("data", "raw", "covariates", "terrain_srtm")
out_static <- here::here("data", "processed", "covariates_static_spatial.csv")
out_mining <- here::here("data", "processed", "covariates_mining_2017_2024.csv")
QA <- here::here("outputs", "covariates", "05I_static_spatial")
for (d in c(terrain_dir, dirname(out_static), QA)) dir.create(d,
  recursive = TRUE, showWarnings = FALSE)
for (p in c(boundary_file,
  adj_file)) if (!file.exists(p)) stop("Missing required input: ", p)

norm <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", " ") |>
    stringr::str_squish()
}

load_sf_rds <- function(name, required = FALSE) {
  p <- file.path(moz_dir, paste0(name, ".rds"))
  if (!file.exists(p)) {
    if (required) stop("Missing MozGIS layer: ", p)
    return(NULL)
  }
  x <- readRDS(p)
  if (!inherits(x, "sf")) stop("Not sf: ", p)
  sf::st_make_valid(x)
}

find_col <- function(nms, patterns) {
  for (pat in patterns) {
    hit <- nms[stringr::str_detect(nms, regex(pat, ignore_case = TRUE))]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

deduplicate_tower_sites <- function(x) {
  if (is.null(x) || !nrow(x)) {
    return(list(data = x, qa = tibble(
      input_rows = 0L,
      output_sites = 0L,
      duplicate_rows_removed = 0L,
      id_field = NA_character_,
      method = "no_data"
    )))
  }
  
  # MozGIS ServiTores2020 exposes SiteID. Prefer it because multiple records
  # can otherwise represent the same physical site. Rows without a usable
  # SiteID are retained; we never collapse nearby coordinates heuristically.
  id_col <- find_col(
    names(x),
    c("^SiteID$", "^site_id$", "^siteid$", "site.*id")
  )
  
  n_in <- nrow(x)
  
  if (!is.na(id_col)) {
    id <- stringr::str_squish(as.character(x[[id_col]]))
    usable <- !is.na(id) & nzchar(id)
    
    keep_with_id <- which(usable & !duplicated(id))
    keep_without_id <- which(!usable)
    keep <- sort(c(keep_with_id, keep_without_id))
    
    out <- x[keep, , drop = FALSE]
    method <- paste0("distinct non-missing ", id_col, "; missing IDs retained")
  } else {
    out <- x
    method <- "no site-ID field found; no deduplication applied"
  }
  
  list(
    data = out,
    qa = tibble(
      input_rows = n_in,
      output_sites = nrow(out),
      duplicate_rows_removed = n_in - nrow(out),
      id_field = if (is.na(id_col)) NA_character_ else id_col,
      method = method
    )
  )
}

parse_any_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, c("POSIXct", "POSIXt"))) return(as.Date(x))
  if (is.numeric(x)) {
    med <- suppressWarnings(stats::median(abs(x), na.rm = TRUE))
    if (is.finite(med) && med > 1e11) return(as.Date(as.POSIXct(x / 1000,
      origin = "1970-01-01", tz = "UTC")))
    if (is.finite(med) && med > 1e8) return(as.Date(as.POSIXct(x,
      origin = "1970-01-01", tz = "UTC")))
  }
  z <- suppressWarnings(lubridate::parse_date_time(as.character(x),
    orders = c("Ymd", "Y-m-d", "dmy", "dmY", "mdY", "m/d/Y",
    "d/m/Y"), quiet = TRUE))
  as.Date(z)
}

# Canonical districts ----------------------------------------------------------
d <- readRDS(boundary_file) |> sf::st_make_valid() |> arrange(district_id)
if (nrow(d) != 17L || n_distinct(d$district_id) != 17L) stop("Expected exactly 17 districts.")
dm <- sf::st_transform(d, crs_m)
rep_pts <- suppressWarnings(sf::st_point_on_surface(dm))
area_km2 <- as.numeric(sf::st_area(dm)) / 1e6
base <- tibble(
  district_id = d$district_id,
  canonical_name = d$canonical_name,
  district_area_km2 = area_km2
)

# Generic point/small-feature district metrics --------------------------------
feature_metrics <- function(x, prefix) {
  if (is.null(x) || !nrow(x)) {
    return(tibble(district_id = d$district_id,
                  !!paste0(prefix, "_count") := NA_integer_,
                  !!paste0("distance_", prefix, "_km") := NA_real_))
  }
  xm <- sf::st_transform(sf::st_make_valid(x), crs_m)
  # One representative point per external feature avoids geometry-boundary
  # double-counting. Counts include only features whose representative point
  # intersects a study district. Nearby features outside Cabo Delgado remain
  # eligible for distance calculations but are never assigned to the nearest
  # district merely because they fall inside the 50-km download buffer.
  xp <- suppressWarnings(sf::st_point_on_surface(xm))
  membership <- sf::st_intersects(xp, dm)
  inside_idx <- which(lengths(membership) > 0L)
  assigned_d <- if (length(inside_idx)) vapply(membership[inside_idx],
    function(z) z[1], integer(1)) else integer()
  count_tab <- tibble(idx = assigned_d) |> count(idx, name = "n")
  counts <- rep(0L, nrow(dm))
  if (nrow(count_tab)) counts[count_tab$idx] <- count_tab$n
  dst <- as.matrix(units::drop_units(sf::st_distance(rep_pts, xp)))
  mind <- if (ncol(dst)) apply(dst, 1, min,
    na.rm = TRUE) / 1000 else rep(NA_real_, nrow(dm))
  tibble(
    district_id = d$district_id,
    !!paste0(prefix, "_count") := counts,
    !!paste0("distance_", prefix, "_km") := as.numeric(mind)
  )
}

# Roads ------------------------------------------------------------------------
road_source <- "MozGIS_current_or_unknown_vintage"

if (require_osm2017_roads_cfg && !file.exists(osm2017)) {
  stop(
    "Historical OSM 2017 roads are required for the primary analysis but are missing: ",
    osm2017,
    "\nRun 05a2_download_roads.R first, or set ",
    "REQUIRE_OSM2017_ROADS <- FALSE only for an explicit sensitivity run."
  )
}

if (file.exists(osm2017)) {
  roads_all <- readRDS(osm2017) |> sf::st_make_valid()
  road_source <- "Geofabrik_OSM_2017-01-01"
  # OSM major-road definition: motorway/trunk/primary and link classes.
  class_col <- find_col(names(roads_all), c("^fclass$", "highway", "class"))
  if (!is.na(class_col)) {
    cls <- norm(roads_all[[class_col]])
    major <- cls %in% c("motorway", "motorway link", "trunk", "trunk link",
      "primary", "primary link")
    roads_primary <- roads_all[major, , drop = FALSE]
  } else {
    stop(
      "OSM 2017 road class field was not identified. ",
      "Cannot construct the prespecified primary-road exposure."
    )
  }
  road_class_layers <- list(all = roads_all, primary = roads_primary)
} else {
  road_names <- c("roads_primary", "roads_secondary", "roads_tertiary",
    "roads_vicinal", "roads_other")
  road_class_layers <- lapply(road_names, load_sf_rds)
  names(road_class_layers) <- sub("roads_", "", road_names)
  road_class_layers <- road_class_layers[!vapply(road_class_layers, is.null,
    logical(1))]
  if (!length(road_class_layers)) stop("No road data. Run 05A or 05A2.")
  roads_all <- bind_rows(road_class_layers)
  roads_primary <- road_class_layers[["primary"]]
}

road_length_by_district <- function(x, label) {
  if (is.null(x) || !nrow(x)) return(tibble(district_id = d$district_id,
    !!paste0("road_km_", label) := NA_real_))
  xm <- sf::st_transform(x, crs_m)
  inter <- suppressWarnings(sf::st_intersection(xm |> select(geometry),
    dm |> select(district_id)))
  if (!nrow(inter)) return(tibble(district_id = d$district_id,
    !!paste0("road_km_", label) := 0))
  z <- inter |> mutate(length_km = as.numeric(sf::st_length(geometry)) / 1000) |>
    sf::st_drop_geometry() |> group_by(district_id) |> summarise(km = sum(length_km, na.rm = TRUE), .groups = "drop")
  tibble(district_id = d$district_id) |> left_join(z, by = "district_id") |>
    mutate(km = replace_na(km, 0)) |> rename(!!paste0("road_km_", label) := km)
}

road_metrics <- base |> select(district_id, district_area_km2)
for (nm in names(road_class_layers)) {
  road_metrics <- road_metrics |> left_join(road_length_by_district(road_class_layers[[nm]], nm), by = "district_id")
}
if (!"road_km_all" %in% names(road_metrics)) {
  road_metrics <- road_metrics |> left_join(road_length_by_district(roads_all,
    "all"), by = "district_id")
}
for (nm in grep("^road_km_", names(road_metrics), value = TRUE)) {
  label <- sub("^road_km_", "", nm)
  road_metrics[[paste0("road_density_", label,
    "_km_km2")]] <- road_metrics[[nm]] / road_metrics$district_area_km2
}
# Standard names expected by registry.
if (!"road_density_primary_km_km2" %in% names(road_metrics)) road_metrics$road_density_primary_km_km2 <- NA_real_
if (!"road_density_all_km_km2" %in% names(road_metrics)) road_metrics$road_density_all_km_km2 <- road_metrics$road_km_all / road_metrics$district_area_km2

road_dist <- function(x) {
  if (is.null(x) || !nrow(x)) return(rep(NA_real_, nrow(dm)))
  xm <- sf::st_transform(x, crs_m)
  as.numeric(units::drop_units(sf::st_distance(rep_pts,
    sf::st_union(xm)))) / 1000
}
road_metrics$distance_primary_road_km <- road_dist(roads_primary)
road_metrics$distance_any_road_km <- road_dist(roads_all)
road_metrics$road_source <- road_source
road_metrics <- road_metrics |> select(-district_area_km2)

# Coast / Tanzania / Nampula / Palma ------------------------------------------
coast <- load_sf_rds("coastline", required = TRUE) |> sf::st_transform(crs_m)
distance_coast <- as.numeric(units::drop_units(sf::st_distance(dm,
  sf::st_union(coast)))) / 1000
coastal <- distance_coast <= 1 # 1-km tolerance for cartographic coastline gaps/islands

# Natural Earth stores the sovereign/admin label as
# "United Republic of Tanzania", so ne_countries(country = "Tanzania") is not
# reliable across rnaturalearth/Natural Earth versions. Select by stable ISO3
# code instead.
world_ne <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

iso3_cols <- intersect(
  c("adm0_a3", "iso_a3", "iso_a3_eh", "sov_a3"),
  names(world_ne)
)

if (!length(iso3_cols)) {
  stop("Natural Earth country table has no recognised ISO3 field.")
}

tanzania_hit <- Reduce(
  `|`,
  lapply(
    iso3_cols,
    function(cc) as.character(world_ne[[cc]]) == "TZA"
  )
)

if (sum(tanzania_hit, na.rm = TRUE) != 1L) {
  readr::write_csv(
    sf::st_drop_geometry(world_ne),
    file.path(QA, "ERROR_naturalearth_country_attributes.csv"),
    na = ""
  )
  
  stop(
    "Could not uniquely identify Tanzania by ISO3 code TZA in Natural Earth."
  )
}

tanzania <- world_ne[tanzania_hit, , drop = FALSE] |>
  sf::st_make_valid() |>
  sf::st_transform(crs_m)

distance_tanzania <- as.numeric(
  units::drop_units(
    sf::st_distance(dm, sf::st_union(tanzania))
  )
) / 1000

touches_tanzania <- distance_tanzania <= 1

provs <- load_sf_rds("provinces_2017", required = TRUE)
char_cols <- names(provs)[vapply(sf::st_drop_geometry(provs), is.character,
  logical(1))]
match_namp <- rep(FALSE, nrow(provs))
for (cc in char_cols) match_namp <- match_namp | stringr::str_detect(norm(provs[[cc]]), "nampula")
if (sum(match_namp) != 1L) {
  readr::write_csv(sf::st_drop_geometry(provs), file.path(QA,
    "ERROR_province_attributes.csv"), na = "")
  stop(
    "Could not uniquely identify Nampula province in MozGIS provinces ",
    "layer.")
}
nampula <- provs[match_namp, ] |> sf::st_transform(crs_m)
distance_nampula <- as.numeric(units::drop_units(sf::st_distance(dm,
  nampula))) / 1000
touches_nampula <- distance_nampula <= 1

palma_i <- which(d$district_id == "CD_015")
if (length(palma_i) != 1L) stop("Palma ID not found.")
distance_palma <- as.numeric(units::drop_units(sf::st_distance(dm, dm[palma_i,
  ]))) / 1000

# Use the district IDs stored with the adjacency object rather than assuming
# that
# sf row order and nb index order are identical. This makes the Palma-neighbour
# indicator robust to future reordering of either object.
adj <- readRDS(adj_file)
nb_raw <- adj$nb_queen_raw
if (is.null(nb_raw)) stop("Adjacency RDS lacks nb_queen_raw.")
if (is.null(adj$district_index) || !"district_id" %in% names(adj$district_index)) {
  stop(
    "Adjacency RDS lacks district_index$district_id; cannot map nb indices ",
    "safely.")
}
adj_ids <- as.character(adj$district_index$district_id)
if (length(adj_ids) != length(nb_raw) || !setequal(adj_ids,
  as.character(d$district_id))) {
  stop("Adjacency district IDs do not match the 17 canonical district IDs.")
}
palma_adj_i <- match("CD_015", adj_ids)
palma_neighbor_by_adj <- vapply(
  seq_along(nb_raw),
  function(i) {
    js <- as.integer(nb_raw[[i]])
    js <- js[js > 0L]
    palma_adj_i %in% js
  },
  logical(1)
)
palma_neighbor_by_adj[palma_adj_i] <- FALSE
palma_neighbor_raw <- palma_neighbor_by_adj[match(as.character(d$district_id),
  adj_ids)]

# Afungi / Mozambique LNG project point. This is an external infrastructure
# location, not an ACLED event coordinate. Coordinates are from the OSM/MapCarta
# feature for Mozambique LNG Project and should be retained with source
# metadata.
afungi_lng <- sf::st_as_sf(
  tibble(name = "Mozambique LNG Project / Afungi", lon = 40.56443,
    lat = -10.82188),
  coords = c("lon", "lat"), crs = 4326
) |> sf::st_transform(crs_m)
distance_afungi_lng <- as.numeric(units::drop_units(sf::st_distance(rep_pts,
  afungi_lng))) / 1000

geo_metrics <- tibble(
  district_id = d$district_id,
  coastal = coastal,
  distance_coast_km = distance_coast,
  touches_tanzania = touches_tanzania,
  distance_tanzania_km = distance_tanzania,
  touches_nampula = touches_nampula,
  distance_nampula_km = distance_nampula,
  palma_neighbor_raw = palma_neighbor_raw,
  distance_palma_km = distance_palma,
  distance_afungi_lng_km = distance_afungi_lng
)

# Infrastructure/state-presence snapshots -------------------------------------
airports <- load_sf_rds("airports_2022")
if (is.null(airports) || !nrow(airports)) airports <- load_sf_rds("aerodromes_transport")
ports <- load_sf_rds("ports")
border_posts <- load_sf_rds("border_posts")
hospitals <- load_sf_rds("hospitals")
schools <- load_sf_rds("education_all")
towers_raw <- load_sf_rds("cell_towers_2020")
tower_dedup <- deduplicate_tower_sites(towers_raw)
towers <- tower_dedup$data

readr::write_csv(
  tower_dedup$qa,
  file.path(QA, "00_cell_tower_site_deduplication.csv"),
  na = ""
)

infra <- base |> select(district_id, district_area_km2)
for (spec in list(
  list(x = airports, prefix = "airport"),
  list(x = ports, prefix = "port"),
  list(x = border_posts, prefix = "border_post"),
  list(x = hospitals, prefix = "hospital"),
  list(x = schools, prefix = "school"),
  list(x = towers, prefix = "cell_tower")
)) {
  infra <- infra |> left_join(feature_metrics(spec$x, spec$prefix),
    by = "district_id")
}
infra <- infra |>
  mutate(
    school_density = school_count / district_area_km2 * 1000,
    hospital_density = hospital_count / district_area_km2 * 1000,
    cell_tower_density = cell_tower_count / district_area_km2 * 1000,
    airport_density = airport_count / district_area_km2 * 1000,
    port_density = port_count / district_area_km2 * 1000,
    border_post_density = border_post_count / district_area_km2 * 1000
  ) |>
  select(-district_area_km2)

# Optional tower coverage share: only if radius unit is explicitly supplied.
cell_coverage_share <- rep(NA_real_, nrow(dm))
if (!is.null(towers) && nrow(towers) && !is.na(cell_tower_radius_unit_cfg)) {
  radius_col <- find_col(names(towers), c("raio", "radius"))
  if (is.na(radius_col)) stop("CELL_TOWER_RADIUS_UNIT set, but no tower radius field found.")
  radius <- suppressWarnings(as.numeric(towers[[radius_col]]))
  if (any(!is.finite(radius),
    na.rm = TRUE)) radius[!is.finite(radius)] <- NA_real_
  if (cell_tower_radius_unit_cfg == "km") radius <- radius * 1000
  if (!cell_tower_radius_unit_cfg %in% c("m",
    "km")) stop("CELL_TOWER_RADIUS_UNIT must be 'm', 'km', or NA.")
  tm <- sf::st_transform(towers, crs_m)
  keep <- which(is.finite(radius) & radius > 0)
  if (length(keep)) {
    cover <- suppressWarnings(sf::st_buffer(tm[keep, ], dist = radius[keep]))
    cu <- sf::st_union(cover)
    for (i in seq_len(nrow(dm))) {
      x <- suppressWarnings(sf::st_intersection(dm[i, ], cu))
      cell_coverage_share[i] <- if (nrow(x)) min(1,
        as.numeric(sf::st_area(x)) / as.numeric(sf::st_area(dm[i, ]))) else 0
    }
  }
}
infra$cell_coverage_share <- cell_coverage_share

# Terrain: SRTM 3 arc-second (~90m) via geodata -------------------------------
#
# geodata::elevation_3s() is tile based and takes longitude + latitude, not a
# country argument. CGIAR-SRTM is distributed in 5 x 5 degree tiles. We
# therefore determine every 5-degree tile intersecting the canonical Cabo
# Delgado extent, request one point at the centre of each tile, merge the
# returned rasters, and only then derive terrain metrics.
#
# This avoids:
#   * silently degrading to the 30-arc-second (~1 km) country product;
#   * assuming that a single SRTM tile covers the entire province;
#   * computing slope/TRI without neighbouring cells at the study boundary.

srtm_tile_centres <- function(x) {
  x_ll <- sf::st_transform(x, 4326)
  bb <- sf::st_bbox(x_ll)
  
  # SRTM 3-second CGIAR tiles are aligned on 5-degree boundaries.
  # Subtract a tiny epsilon from xmax/ymax so an extent lying exactly on a
  # tile boundary does not request an unnecessary neighbouring tile.
  eps <- 1e-10
  
  lon0 <- seq(
    floor(as.numeric(bb["xmin"]) / 5) * 5,
    floor((as.numeric(bb["xmax"]) - eps) / 5) * 5,
    by = 5
  )
  
  lat0 <- seq(
    floor(as.numeric(bb["ymin"]) / 5) * 5,
    floor((as.numeric(bb["ymax"]) - eps) / 5) * 5,
    by = 5
  )
  
  tidyr::crossing(
    tile_lon_min = lon0,
    tile_lat_min = lat0
  ) |>
    mutate(
      request_lon = tile_lon_min + 2.5,
      request_lat = tile_lat_min + 2.5
    ) |>
    arrange(tile_lat_min, tile_lon_min)
}

message("Loading/downloading CGIAR-SRTM 3-arc-second tiles for Cabo Delgado...")

srtm_requests <- srtm_tile_centres(d)

srtm_tiles <- purrr::map2(
  srtm_requests$request_lon,
  srtm_requests$request_lat,
  function(lon, lat) {
    geodata::elevation_3s(
      lon = lon,
      lat = lat,
      path = terrain_dir
    )
  }
)

if (!length(srtm_tiles)) {
  stop("No SRTM tiles were identified for the Cabo Delgado extent.")
}

# Record exactly which tile files were used.
srtm_requests$source_file <- vapply(
  srtm_tiles,
  function(r) paste(terra::sources(r), collapse = ";"),
  character(1)
)

readr::write_csv(
  srtm_requests,
  file.path(QA, "00_srtm_tiles_requested.csv"),
  na = ""
)

# Adjacent CGIAR tiles have the same resolution/grid alignment. Merge preserves
# their original elevation values; no resampling or averaging is introduced.
elev <- if (length(srtm_tiles) == 1L) {
  srtm_tiles[[1]]
} else {
  Reduce(function(x, y) terra::merge(x, y), srtm_tiles)
}

# Hard check that the merged raster actually covers the whole canonical extent.
d_ll <- sf::st_transform(d, 4326)
bb <- sf::st_bbox(d_ll)

if (
  terra::xmin(elev) > as.numeric(bb["xmin"]) ||
  terra::xmax(elev) < as.numeric(bb["xmax"]) ||
  terra::ymin(elev) > as.numeric(bb["ymin"]) ||
  terra::ymax(elev) < as.numeric(bb["ymax"])
) {
  stop("Merged SRTM tiles do not cover the full canonical Cabo Delgado extent.")
}

# Keep a 1-km collar around the study polygons before calculating slope/TRI so
# district-edge cells have neighbouring raster cells available.
terrain_buffer_ll <- sf::st_buffer(
  sf::st_union(dm),
  dist = 1000
) |>
  sf::st_transform(4326)

elev_cd <- terra::crop(
  elev,
  terra::vect(terrain_buffer_ll)
)

vv <- terra::vect(d_ll)

slope <- terra::terrain(
  elev_cd,
  v = "slope",
  unit = "degrees",
  neighbors = 8
)

tri <- terra::terrain(
  elev_cd,
  v = "TRI",
  neighbors = 8
)

extract_stat <- function(r, fun) {
  z <- terra::extract(r, vv, fun = fun, na.rm = TRUE)
  value_col <- setdiff(names(z), "ID")
  if (length(value_col) != 1L) stop("Unexpected terrain extraction schema.")
  as.numeric(z[[value_col]])
}
terrain <- tibble(
  district_id = d$district_id,
  elevation_mean_m = extract_stat(elev_cd, mean),
  elevation_sd_m = extract_stat(elev_cd, stats::sd),
  slope_mean_deg = extract_stat(slope, mean),
  terrain_ruggedness = extract_stat(tri, mean)
)

# Mining concessions: monthly historical exposure from concession/validity dates
mining <- load_sf_rds("mining_concessions")
months <- seq.Date(start_cfg, end_cfg, by = "month")
if (is.null(mining) || !nrow(mining)) {
  warning("No mining concession layer available; mining covariates will be NA.")
  mining_monthly <- tidyr::crossing(district_id = d$district_id,
    year_month = months) |>
    mutate(mining_active_area_share = NA_real_,
      mining_active_license_count = NA_integer_, mining_new_license_12m = NA_integer_)
} else {
  start_col <- find_col(names(mining), c("DataConcessao", "concess", "grant",
    "issue"))
  end_col <- find_col(names(mining), c("DataValidade", "valid", "expir",
    "expiry"))
  id_col <- find_col(names(mining), c("NuLicenca", "licen", "license",
    "objectid"))
  if (is.na(start_col)) {
    readr::write_csv(sf::st_drop_geometry(mining) |> slice_head(n = 50),
      file.path(QA, "ERROR_mining_schema_sample.csv"), na = "")
    stop(
      "Mining concession start/concession-date field not found; do not ",
      "fabricate historical activity.")
  }
  mining$license_start <- parse_any_date(mining[[start_col]])
  mining$license_end <- if (!is.na(end_col)) parse_any_date(mining[[end_col]]) else as.Date(NA)
  mining$license_id <- if (!is.na(id_col)) as.character(mining[[id_col]]) else sprintf("LIC_%06d", seq_len(nrow(mining)))
  
  invalid_dates <- !is.na(mining$license_start) & !is.na(mining$license_end) &
    mining$license_end < mining$license_start
  if (any(invalid_dates)) {
    readr::write_csv(
      sf::st_drop_geometry(mining[invalid_dates, ]),
      file.path(QA, "ERROR_mining_end_before_start.csv"),
      na = ""
    )
    stop("Mining concessions contain end dates before start dates.")
  }
  mm <- sf::st_transform(sf::st_make_valid(mining), crs_m)
  inter <- suppressWarnings(sf::st_intersection(
    mm |> select(license_id, license_start, license_end, geometry),
    dm |> select(district_id)
  ))
  inter$piece_area_m2 <- as.numeric(sf::st_area(inter))
  # Records with no concession date cannot support historical active status.
  readr::write_csv(
    sf::st_drop_geometry(inter) |> filter(is.na(license_start)) |> distinct(license_id, .keep_all = TRUE),
    file.path(QA, "00_mining_missing_start_dates.csv"), na = ""
  )
  inter_valid <- inter |> filter(!is.na(license_start))
  
  one_month <- function(m) {
    mend <- lubridate::ceiling_date(m, "month") - lubridate::days(1)
    m12 <- m %m-% months(11)
    purrr::map_dfr(seq_len(nrow(dm)), function(i) {
      did <- dm$district_id[i]
      z <- inter_valid |> filter(district_id == did)
      active <- z |> filter(license_start <= mend,
        is.na(license_end) | license_end >= m)
      if (nrow(active)) {
        union_geom <- sf::st_union(active)
        active_area <- as.numeric(sf::st_area(union_geom))
        active_count <- n_distinct(active$license_id)
      } else {
        active_area <- 0; active_count <- 0L
      }
      new_count <- z |> sf::st_drop_geometry() |> filter(license_start >= m12,
        license_start <= mend) |> summarise(n = n_distinct(license_id)) |> pull(n)
      tibble(
        district_id = did,
        year_month = m,
        mining_active_area_share = min(1,
          active_area / as.numeric(sf::st_area(dm[i, ]))),
        mining_active_license_count = as.integer(active_count),
        mining_new_license_12m = as.integer(new_count)
      )
    })
  }
  mining_monthly <- bind_rows(lapply(months, one_month)) |> arrange(district_id,
    year_month)
}
if (nrow(mining_monthly) != 1632L) stop("Mining output is not 17 x 96.")

# Final static table ------------------------------------------------------------
static <- base |>
  left_join(road_metrics, by = "district_id") |>
  left_join(geo_metrics, by = "district_id") |>
  left_join(infra, by = "district_id") |>
  left_join(terrain, by = "district_id") |>
  arrange(district_id)
if (nrow(static) != 17L || anyDuplicated(static$district_id)) stop("Static covariate key corrupted.")

readr::write_csv(static, out_static, na = "")
readr::write_csv(mining_monthly, out_mining, na = "")
numeric_summary <- static |>
  summarise(
    across(
      where(is.numeric),
      list(
        missing = ~sum(is.na(.x)),
        min = ~safe_min(.x),
        max = ~safe_max(.x)
      )
    )
  )

readr::write_csv(
  numeric_summary,
  file.path(QA, "01_static_numeric_summary.csv"),
  na = ""
)

readr::write_csv(
  tibble(
    road_source = road_source,
    require_osm2017_roads = require_osm2017_roads_cfg,
    cell_tower_radius_unit = cell_tower_radius_unit_cfg
  ),
  file.path(QA, "02_source_choices.csv"),
  na = ""
)

message("Static spatial covariates written: ", out_static)
message("Monthly mining covariates written: ", out_mining)
message("Road source: ", road_source)
message(
  "School/hospital/tower/airport snapshots are sensitivity variables because ",
  "their vintages can post-date conflict onset."
)
