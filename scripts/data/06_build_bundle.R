# Assemble the analysis data object read by the modelling script.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  sf,
  dplyr,
  tibble,
  readr,
  stringr,
  lubridate,
  units,
  geodata,
  here
)

panel_file <- here::here(
  "data", "processed",
  "cabo_delgado_panel_with_covariates_standardized.rds"
)

boundary_file <- here::here(
  "data", "processed",
  "cabo_delgado_admin2_canonical.rds"
)

adjacency_file <- here::here(
  "data", "processed",
  "cabo_delgado_adjacency_topological.rds"
)

event_file <- here::here(
  "data", "processed",
  "cabo_delgado_event_level_2017_2024.rds"
)

registry_file <- here::here(
  "data", "processed",
  "cabo_delgado_covariate_registry.csv"
)

eligibility_file <- here::here(
  "data", "processed",
  "cabo_delgado_covariate_eligibility_matrix.csv"
)

dictionary_file <- here::here(
  "data", "processed",
  "cabo_delgado_district_month_panel_dictionary.csv"
)

standardization_map_file <- here::here(
  "outputs",
  "covariates",
  "05L_merge_QA",
  "07_standardization_name_map.csv"
)

mozgis_dir <- here::here(
  "data", "raw",
  "covariates",
  "mozgis"
)

out_file <- here::here(
  "data",
  "cabo_delgado_data.rds"
)

required_files <- c(
  panel_file,
  boundary_file,
  adjacency_file,
  event_file,
  registry_file,
  eligibility_file,
  dictionary_file,
  standardization_map_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files)) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

read_first_sf <- function(candidates, label) {
  paths <- file.path(
    mozgis_dir,
    paste0(candidates, ".rds")
  )
  paths <- paths[file.exists(paths)]
  
  if (!length(paths)) {
    stop(
      "Missing strategic spatial layer for ",
      label,
      ": ",
      paste(candidates, collapse = " / ")
    )
  }
  
  x <- readRDS(paths[1])
  
  if (!inherits(x, "sf")) {
    stop(label, " layer is not an sf object: ", paths[1])
  }
  
  x <- sf::st_make_valid(x)
  
  list(
    data = x,
    source_file = paths[1]
  )
}

pick_asset_name_field <- function(
    x,
    patterns
) {
  attrs <- sf::st_drop_geometry(x)
  nms <- names(attrs)
  
  if (!length(nms)) {
    return(NA_character_)
  }
  
  for (pat in patterns) {
    hit <- nms[
      stringr::str_detect(
        stringr::str_to_lower(nms),
        stringr::regex(
          pat,
          ignore_case = TRUE
        )
      )
    ]
    
    for (nm in hit) {
      value <- stringr::str_squish(
        as.character(
          attrs[[nm]]
        )
      )
      
      usable <- !is.na(value) &
        nzchar(value)
      
      if (any(usable)) {
        return(nm)
      }
    }
  }
  
  character_fields <- nms[
    vapply(
      attrs,
      function(z) {
        is.character(z) ||
          is.factor(z)
      },
      logical(1)
    )
  ]
  
  for (nm in character_fields) {
    value <- stringr::str_squish(
      as.character(
        attrs[[nm]]
      )
    )
    
    usable <- !is.na(value) &
      nzchar(value)
    
    if (any(usable)) {
      return(nm)
    }
  }
  
  NA_character_
}

title_case_if_upper <- function(x) {
  out <- stringr::str_squish(x)
  
  upper_like <- !is.na(out) &
    nzchar(out) &
    out ==
    stringr::str_to_upper(out)
  
  out[upper_like] <-
    stringr::str_to_title(
      stringr::str_to_lower(
        out[upper_like]
      )
    )
  
  out
}

prepare_asset_labels <- function(
    x,
    asset_type,
    field_patterns
) {
  name_field <-
    pick_asset_name_field(
      x,
      field_patterns
    )
  
  if (is.na(name_field)) {
    raw_name <- rep(
      NA_character_,
      nrow(x)
    )
  } else {
    raw_name <-
      as.character(
        x[[name_field]]
      )
  }
  
  raw_name <-
    title_case_if_upper(
      raw_name
    )
  
  raw_name[
    is.na(raw_name) |
      !nzchar(raw_name)
  ] <- NA_character_
  
  label <-
    switch(
      asset_type,
      
      "Airport / aerodrome" =
        ifelse(
          is.na(raw_name),
          paste0(
            "Airport / Aerodrome ",
            seq_len(nrow(x))
          ),
          ifelse(
            stringr::str_detect(
              stringr::str_to_lower(raw_name),
              "airport|aeroport|aerodrom|airfield"
            ),
            raw_name,
            paste0(
              raw_name,
              " Airport / Aerodrome"
            )
          )
        ),
      
      "Port" =
        ifelse(
          is.na(raw_name),
          paste0(
            "Port ",
            seq_len(nrow(x))
          ),
          ifelse(
            stringr::str_detect(
              stringr::str_to_lower(raw_name),
              "port|porto"
            ),
            raw_name,
            paste0(
              "Port of ",
              raw_name
            )
          )
        ),
      
      "Border post" =
        ifelse(
          is.na(raw_name),
          paste0(
            "Border Post ",
            seq_len(nrow(x))
          ),
          ifelse(
            stringr::str_detect(
              stringr::str_to_lower(raw_name),
              "border|fronteira|post"
            ),
            raw_name,
            paste0(
              raw_name,
              " Border Post"
            )
          )
        ),
      
      raw_name
    )
  
  x$asset_type <- asset_type
  x$asset_label <- label
  x$asset_label_source_field <-
    if (is.na(name_field)) {
      NA_character_
    } else {
      name_field
    }
  
  x
}

nearest_geodesic_distance_km <- function(
    points_sf,
    target_sf
) {
  if (!nrow(target_sf)) {
    return(
      rep(
        NA_real_,
        nrow(points_sf)
      )
    )
  }
  
  p <- sf::st_transform(
    points_sf,
    4326
  )
  
  z <- sf::st_transform(
    target_sf,
    4326
  )
  
  old_s2 <- sf::sf_use_s2()
  on.exit(
    sf::sf_use_s2(old_s2),
    add = TRUE
  )
  sf::sf_use_s2(TRUE)
  
  d <- sf::st_distance(
    p,
    z
  )
  
  apply(
    as.matrix(
      units::drop_units(d)
    ),
    1,
    min,
    na.rm = TRUE
  ) / 1000
}

panel <- readRDS(panel_file) |>
  as_tibble() |>
  arrange(
    district_id,
    time_id
  )

districts <- readRDS(boundary_file)
adjacency <- readRDS(adjacency_file)
events <- readRDS(event_file)

registry <- readr::read_csv(
  registry_file,
  show_col_types = FALSE
)

eligibility <- readr::read_csv(
  eligibility_file,
  show_col_types = FALSE
)

dictionary <- readr::read_csv(
  dictionary_file,
  show_col_types = FALSE
)

standardization_map <- readr::read_csv(
  standardization_map_file,
  show_col_types = FALSE
)

if (!inherits(districts, "sf")) {
  stop("Canonical boundary must be an sf object.")
}

if (
  nrow(panel) != 1632L ||
  dplyr::n_distinct(panel$district_id) != 17L ||
  dplyr::n_distinct(panel$time_id) != 96L ||
  anyDuplicated(panel$cell_id)
) {
  stop("Final panel is not the approved 17 x 96 district-month design.")
}

if (
  nrow(districts) != 17L ||
  anyDuplicated(districts$district_id)
) {
  stop("Canonical geometry does not contain exactly 17 unique districts.")
}

if (is.null(adjacency$nb_queen_connected)) {
  stop("Approved connected adjacency is missing from the adjacency object.")
}

id_order <- adjacency$district_index$district_id

if (
  is.null(id_order) ||
  !setequal(id_order, districts$district_id) ||
  !setequal(id_order, unique(panel$district_id))
) {
  stop("Panel, geometry, and adjacency district IDs are inconsistent.")
}

districts <- districts[
  match(
    id_order,
    districts$district_id
  ),
]

if (inherits(events, "sf")) {
  events_sf <- events
} else {
  required_event_fields <- c(
    "longitude",
    "latitude"
  )
  
  if (!all(required_event_fields %in% names(events))) {
    stop(
      "Event-level data must contain longitude and latitude."
    )
  }
  
  events_sf <- sf::st_as_sf(
    events,
    coords = c(
      "longitude",
      "latitude"
    ),
    crs = 4326,
    remove = FALSE
  )
}

required_event_fields <- c(
  "event_date",
  "district_id",
  "fatalities",
  "event_type",
  "sub_event_type"
)

if (!all(required_event_fields %in% names(events_sf))) {
  stop(
    "Event-level data are missing required field(s): ",
    paste(
      setdiff(
        required_event_fields,
        names(events_sf)
      ),
      collapse = ", "
    )
  )
}

events_sf <- events_sf |>
  mutate(
    event_date = as.Date(event_date),
    year_month = as.Date(
      lubridate::floor_date(
        event_date,
        "month"
      )
    ),
    fatalities = as.numeric(fatalities)
  )

if (
  any(
    !events_sf$district_id %in%
    id_order
  )
) {
  stop(
    "Event-level data contain district IDs outside the canonical ",
    "17-district domain.")
}

airports_obj <- read_first_sf(
  c(
    "airports_2022",
    "aerodromes_transport"
  ),
  "airports/aerodromes"
)

ports_obj <- read_first_sf(
  "ports",
  "ports"
)

border_posts_obj <- read_first_sf(
  "border_posts",
  "border posts"
)

mining_obj <- read_first_sf(
  "mining_concessions",
  "mining concessions"
)

afungi_lng <- sf::st_as_sf(
  tibble(
    name = "Mozambique LNG Project / Afungi",
    longitude = 40.56443,
    latitude = -10.82188
  ),
  coords = c(
    "longitude",
    "latitude"
  ),
  crs = 4326,
  remove = FALSE
)

geo_context <- panel |>
  distinct(
    district_id,
    coastal,
    touches_tanzania,
    touches_nampula,
    distance_coast_km,
    distance_tanzania_km,
    distance_nampula_km,
    distance_palma_km,
    distance_afungi_lng_km,
    distance_airport_km,
    distance_port_km,
    distance_border_post_km,
    distance_primary_road_km
  )

events_sf <- events_sf |>
  left_join(
    geo_context,
    by = "district_id"
  ) |>
  mutate(
    event_geographic_position =
      case_when(
        touches_tanzania & coastal ~ "Tanzania border + coast",
        touches_tanzania ~ "Tanzania border",
        touches_nampula & coastal ~ "Nampula border + coast",
        touches_nampula ~ "Nampula border",
        coastal ~ "Coast",
        TRUE ~ "Interior"
      ),
    event_threat_level =
      case_when(
        fatalities < 5 ~ "Low",
        fatalities <= 10 ~ "Moderate",
        fatalities <= 20 ~ "Substantial",
        fatalities <= 50 ~ "Severe",
        TRUE ~ "Critical"
      )
  )

events_sf$distance_event_airport_km <-
  nearest_geodesic_distance_km(
    events_sf,
    airports_obj$data
  )

events_sf$distance_event_port_km <-
  nearest_geodesic_distance_km(
    events_sf,
    ports_obj$data
  )

events_sf$distance_event_border_post_km <-
  nearest_geodesic_distance_km(
    events_sf,
    border_posts_obj$data
  )

events_sf$distance_event_afungi_lng_km <-
  nearest_geodesic_distance_km(
    events_sf,
    afungi_lng
  )

events_sf$distance_event_mining_concession_km <-
  nearest_geodesic_distance_km(
    events_sf,
    mining_obj$data
  )

airports_display <- prepare_asset_labels(
  airports_obj$data,
  "Airport / aerodrome",
  c(
    "^name$",
    "^nome$",
    "aeroport",
    "aerodrom",
    "airfield",
    "design",
    "local"
  )
)

ports_display <- prepare_asset_labels(
  ports_obj$data,
  "Port",
  c(
    "^name$",
    "^nome$",
    "porto",
    "port",
    "design",
    "local"
  )
)

border_posts_display <- prepare_asset_labels(
  border_posts_obj$data,
  "Border post",
  c(
    "^name$",
    "^nome$",
    "fronteira",
    "border",
    "posto",
    "post",
    "design",
    "local"
  )
)

afungi_lng$asset_type <- "Afungi LNG"
afungi_lng$asset_label <-
  "Mozambique LNG Project / Afungi"
afungi_lng$asset_label_source_field <-
  "name"

mozambique_admin1 <- geodata::gadm(
  country = "MOZ",
  level = 1,
  path = tempdir()
) |>
  sf::st_as_sf() |>
  sf::st_make_valid()

if (!"NAME_1" %in% names(mozambique_admin1)) {
  stop(
    "GADM ADM1 cartographic context is missing NAME_1."
  )
}

mozambique_admin1 <- mozambique_admin1 |>
  mutate(
    is_cabo_delgado =
      stringr::str_to_lower(
        stringr::str_squish(
          NAME_1
        )
      ) ==
      "cabo delgado"
  )

if (
  sum(
    mozambique_admin1$is_cabo_delgado,
    na.rm = TRUE
  ) != 1L
) {
  stop(
    "Mozambique ADM1 inset did not identify exactly one Cabo Delgado province."
  )
}

strategic_geography <- list(
  airports = airports_display,
  ports = ports_display,
  border_posts = border_posts_display,
  mining_concessions = mining_obj$data,
  afungi_lng = afungi_lng
)

bundle <- list(
  panel = panel,
  districts = districts,
  adjacency = adjacency,
  events = events_sf,
  strategic_geography = strategic_geography,
  mozambique_admin1 = mozambique_admin1,
  registry = registry,
  eligibility = eligibility,
  dictionary = dictionary,
  standardization_map = standardization_map,
  metadata = list(
    bundle_version = "3.0.0",
    analysis_unit = "district_month",
    descriptive_event_unit = "ACLED_event",
    study_start = as.Date("2017-01-01"),
    study_end = as.Date("2024-12-01"),
    n_rows = nrow(panel),
    n_events = nrow(events_sf),
    n_districts = dplyr::n_distinct(panel$district_id),
    n_months = dplyr::n_distinct(panel$year_month),
    primary_outcome = "fatalities_all",
    sensitivity_outcome = "fatalities_political_violence",
    geometry_key = "district_id",
    geometry_level = "canonical_admin2",
    inset_geometry = "GADM_ADM1_visual_context_only",
    adjacency_primary = "nb_queen_connected",
    adjacency_sensitivity = "nb_queen_raw",
    event_coordinate_distances_use = "descriptive_only_geodesic_s2",
    strategic_layer_sources = list(
      airports = airports_obj$source_file,
      ports = ports_obj$source_file,
      border_posts = border_posts_obj$source_file,
      mining_concessions = mining_obj$source_file,
      afungi_lng = "OSM/MapCarta point retained from 05I"
    )
  )
)

class(bundle) <- c(
  "cabo_delgado_analysis_bundle",
  "list"
)

dir.create(
  dirname(out_file),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  bundle,
  out_file,
  compress = "xz"
)

message("Complete analysis bundle written: ", out_file)
message(
  "District-month rows: ",
  nrow(panel),
  "; event rows: ",
  nrow(events_sf),
  "; districts: ",
  dplyr::n_distinct(panel$district_id),
  "; months: ",
  dplyr::n_distinct(panel$year_month)
)
message(
  "Single RDS contains panel, canonical geometry, adjacency, event-level data, ",
  "named strategic spatial layers, Mozambique ADM1 inset geometry, registry, ",
  "eligibility, dictionary and metadata."
)
