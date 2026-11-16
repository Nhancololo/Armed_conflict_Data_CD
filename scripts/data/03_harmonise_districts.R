# Harmonise district names and boundaries with the ACLED locations.

# 0. Packages ----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  dplyr,
  tidyr,
  readr,
  stringr,
  stringi,
  tibble,
  purrr,
  sf,
  terra,
  geodata,
  stringdist,
  httr2,
  here,
  fs
)

# 1. Configuration -----------------------------------------------------------
country_cfg  <- "Mozambique"
province <- "Cabo Delgado"

analysis_start <- as.Date("2017-01-01")
analysis_end   <- as.Date("2024-12-31")

raw_acled_root <- here::here("data", "raw", "acled")
processed_dir  <- here::here("data", "processed")

qa_dir <- here::here("outputs", "admin2_harmonization")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

# Local boundary candidates, in priority order.
# Put an official/preferred boundary here if available.
boundary_candidates <- c(
  here::here("data", "raw", "boundaries", "cabo_delgado_admin2_official.rds"),
  here::here("data", "raw", "boundaries", "cabo_delgado_admin2_official.gpkg"),
  here::here("data", "processed", "cabo_delgado_boundary.rds"),
  here::here("data", "processed", "cabo_delgado_boundary2.rds"),
  here::here("cabo_delgado_boundary.rds"),
  here::here("cabo_delgado_boundary2.rds")
)

# Boundary strategy:
#   1) Prefer an explicitly supplied local official/preferred boundary.
#   2) Otherwise download the public CENACARTA/MozGIS district layer.
#   3) GADM is never used as the canonical 17-district geometry here because its
#      Cabo Delgado ADM2 layer currently contains only 15 units. It is
# downloaded
#      only as a diagnostic so that its names/structure are visible and
# auditable.
allow_mozgis_fallback <- TRUE
run_gadm_diagnostic   <- TRUE

mozgis_district_query_url <- paste0(
  "https://www.mozgis.gov.mz/server/rest/services/",
  "Data/Unidades_Administrativas/FeatureServer/2/query"
)

# 2. Helpers -----------------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

is_blankish <- function(x) {
  is.na(x) | stringr::str_trim(as.character(x)) == ""
}

normalize_geo_name <- function(x) {
  x <- as.character(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- stringr::str_to_lower(x)
  x <- stringr::str_replace_all(x, "&", " e ")
  x <- stringr::str_replace_all(x, "[-_/]+", " ")
  x <- stringr::str_replace_all(x, "[^a-z0-9 ]+", " ")
  x <- stringr::str_squish(x)
  
  # Remove generic administrative words only at the edges.
  x <- stringr::str_replace(x, "^distrito de\\s+", "")
  x <- stringr::str_replace(x, "^district of\\s+", "")
  x <- stringr::str_replace(x, "\\s+district$", "")
  x <- stringr::str_squish(x)
  x
}

read_boundary_any <- function(path) {
  ext <- tolower(fs::path_ext(path))
  
  out <- tryCatch({
    if (ext == "rds") {
      readRDS(path)
    } else {
      sf::st_read(path, quiet = TRUE)
    }
  }, error = function(e) NULL)
  
  if (is.null(out)) return(NULL)
  
  if (inherits(out, "SpatVector")) out <- sf::st_as_sf(out)
  if (!inherits(out, "sf")) return(NULL)
  
  out
}

detect_name_field <- function(x) {
  candidates <- c(
    "NAME_2", "NAME_2_EN", "ADM2_EN", "ADM2_NAME", "admin2",
    "Distrito", "DISTRITO", "distrito",
    "district", "DISTRICT", "District", "name_2", "name"
  )
  hit <- intersect(candidates, names(x))
  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

write_qa <- function(x, filename) {
  readr::write_csv(x, file.path(qa_dir, filename), na = "")
  invisible(x)
}

download_mozgis_districts <- function() {
  message(
    "Downloading public CENACARTA/MozGIS district polygons ",
    "(FeatureServer/2)..."
  )
  
  response <- httr2::request(mozgis_district_query_url) |>
    httr2::req_url_query(
      where = "1=1",
      outFields = "Distrito,DistCodigo,Provincia,ProvCodigo",
      returnGeometry = "true",
      outSR = "4326",
      f = "geojson"
    ) |>
    httr2::req_retry(
      max_tries = 5,
      retry_on_failure = TRUE
    ) |>
    httr2::req_timeout(seconds = 120) |>
    httr2::req_perform()
  
  if (httr2::resp_status(response) != 200L) {
    stop(
      "MozGIS/CENACARTA request failed with HTTP status ",
      httr2::resp_status(response), "."
    )
  }
  
  raw_dir <- here::here("data", "raw", "boundaries", "mozgis_cenacarta")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  
  raw_geojson <- file.path(
    raw_dir,
    paste0(
      "cenacarta_mozambique_districts.geojson"
    )
  )
  
  writeBin(httr2::resp_body_raw(response), raw_geojson)
  
  x <- sf::st_read(raw_geojson, quiet = TRUE)
  
  required <- c("Distrito", "Provincia")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      "MozGIS/CENACARTA GeoJSON is missing expected fields: ",
      paste(missing, collapse = ", ")
    )
  }
  
  x_cd <- x |>
    mutate(
      Provincia = as.character(.data$Provincia),
      Distrito = as.character(.data$Distrito)
    ) |>
    filter(
      normalize_geo_name(.data$Provincia) ==
        normalize_geo_name(province)
    )
  
  if (nrow(x_cd) == 0L) {
    stop(
      "MozGIS/CENACARTA returned no district polygons for ",
      province, "."
    )
  }
  
  raw_names <- x_cd |>
    sf::st_drop_geometry() |>
    transmute(
      source = "CENACARTA/MozGIS",
      Provincia = .data$Provincia,
      Distrito = .data$Distrito,
      DistCodigo = if ("DistCodigo" %in% names(x_cd)) {
        as.character(.data$DistCodigo)
      } else {
        NA_character_
      },
      ProvCodigo = if ("ProvCodigo" %in% names(x_cd)) {
        as.character(.data$ProvCodigo)
      } else {
        NA_character_
      },
      normalized = normalize_geo_name(.data$Distrito)
    ) |>
    distinct() |>
    arrange(.data$Distrito)
  
  write_qa(
    raw_names,
    "00B_MOZGIS_CENACARTA_raw_CaboDelgado_district_names.csv"
  )
  
  message(
    "CENACARTA/MozGIS returned ",
    n_distinct(raw_names$Distrito),
    " Cabo Delgado district labels:"
  )
  message("  ", paste(raw_names$Distrito, collapse = " | "))
  
  list(
    boundary = x_cd,
    raw_geojson = raw_geojson,
    raw_names = raw_names
  )
}

# 3. Canonical 17-district nominal domain ------------------------------------
canonical_districts <- tibble::tribble(
  ~district_id, ~canonical_name,
  "CD_001", "Ancuabe",
  "CD_002", "Balama",
  "CD_003", "Chiúre",
  "CD_004", "Ibo",
  "CD_005", "Macomia",
  "CD_006", "Mecúfi",
  "CD_007", "Meluco",
  "CD_008", "Metuge",
  "CD_009", "Mocímboa da Praia",
  "CD_010", "Montepuez",
  "CD_011", "Mueda",
  "CD_012", "Muidumbe",
  "CD_013", "Namuno",
  "CD_014", "Nangade",
  "CD_015", "Palma",
  "CD_016", "Pemba",
  "CD_017", "Quissanga"
) |>
  mutate(
    canonical_normalized = normalize_geo_name(.data$canonical_name),
    valid_from = analysis_start,
    valid_to = analysis_end,
    nominal_reference = "Mozambique INE; Cabo Delgado has 17 districts"
  )

stopifnot(
  nrow(canonical_districts) == 17L,
  n_distinct(canonical_districts$district_id) == 17L,
  n_distinct(canonical_districts$canonical_normalized) == 17L
)

# Explicit aliases. These are validated spelling/diacritic/admin-label variants;
# add new rows here only after manual verification.
explicit_aliases <- tibble::tribble(
  ~alias,                    ~canonical_name,
  "Chiure",                  "Chiúre",
  "Chiúre",                  "Chiúre",
  "Distrito de Chiure",      "Chiúre",
  "Mecufi",                  "Mecúfi",
  "Mecúfi",                  "Mecúfi",
  "Distrito de Mecufi",      "Mecúfi",
  "Mocimboa da Praia",       "Mocímboa da Praia",
  "Mocímboa da Praia",       "Mocímboa da Praia",
  "Mocimboa-da-Praia",       "Mocímboa da Praia",
  "Mocímboa-da-Praia",       "Mocímboa da Praia",
  "Distrito de Mocimboa da Praia", "Mocímboa da Praia",
  "Cidade de Pemba",         "Pemba",
  "Pemba Cidade",             "Pemba",
  "Distrito de Pemba",        "Pemba"
) |>
  mutate(alias_normalized = normalize_geo_name(.data$alias)) |>
  left_join(
    canonical_districts |>
      select("district_id", "canonical_name"),
    by = "canonical_name"
  )

if (any(is.na(explicit_aliases$district_id))) {
  stop("Internal error: an explicit alias points to a non-canonical district.")
}

# 3B. Source-specific historical boundary-name crosswalk ---------------------
# do not merge these rules into the generic aliases above.
#
# Mozambique INE documents that in 2013:
#   * Cidade de Pemba acquired district status; and
#   * the former district Pemba-Metuge was redesignated as Metuge.
#
# The public CENACARTA/MozGIS district layer used below exposes the 17 Cabo
# Delgado polygons with labels including both "CIDADE DE PEMBA" and "PEMBA",
# but not "METUGE". For THIS boundary_cfg SOURCE only we therefore interpret:
#
#   CIDADE DE PEMBA  -> current Pemba
#   PEMBA            -> current Metuge
#
# This mapping is source-specific and must never be applied to ACLED admin2.
boundary_source_crosswalk <- tibble::tribble(
  ~source_pattern,       ~source_name,         ~canonical_name, ~reason,
  "CENACARTA/MozGIS",   "CIDADE DE PEMBA",    "Pemba",
  "INE: Cidade de Pemba became a district in 2013",
  "CENACARTA/MozGIS",   "PEMBA",              "Metuge",
  "INE: former Pemba-Metuge was redesignated Metuge in 2013"
) |>
  mutate(
    source_name_normalized = normalize_geo_name(.data$source_name)
  ) |>
  left_join(
    canonical_districts |>
      select("district_id", "canonical_name"),
    by = "canonical_name"
  )

if (any(is.na(boundary_source_crosswalk$district_id))) {
  stop("Internal error in boundary_source_crosswalk.")
}

# 4. Locate latest ACLED Cabo Delgado typed snapshot -------------------------
if (!dir.exists(raw_acled_root)) {
  stop("ACLED directory not found: ", raw_acled_root)
}

acled_files <- fs::dir_ls(
  raw_acled_root,
  recurse = TRUE,
  type = "file",
  regexp = "acled_cabo_delgado_typed\\.rds$"
)

if (length(acled_files) == 0L) {
  stop(
    "No Cabo Delgado TYPED ACLED RDS found under ", raw_acled_root,
    ". Run script 01 first."
  )
}

file_info <- fs::file_info(acled_files)
acled_input <- as.character(acled_files[which.max(file_info$modification_time)])

acled <- readRDS(acled_input) |>
  tibble::as_tibble()

if (!all(c("event_id_cnty", "event_date", "admin2") %in% names(acled))) {
  stop("ACLED input must contain event_id_cnty, event_date, and admin2.")
}

acled <- acled |>
  mutate(
    event_date = as.Date(.data$event_date),
    admin2_original = as.character(.data$admin2),
    admin2_normalized = normalize_geo_name(.data$admin2_original)
  )

if (anyDuplicated(acled$event_id_cnty) > 0L) {
  stop(
    "Duplicate event_id_cnty detected. Resolve duplicates before ",
    "harmonization.")
}

# 5. Harmonize ACLED admin2 labels -------------------------------------------
canonical_lookup <- canonical_districts |>
  select(
    .data$district_id,
    .data$canonical_name,
    .data$canonical_normalized
  )

# Observed source-name dictionary first; then join it back to events.
observed_names <- acled |>
  filter(!is_blankish(.data$admin2_original)) |>
  count(.data$admin2_original, .data$admin2_normalized, name = "n_events") |>
  arrange(desc(.data$n_events))

# Exact canonical name.
name_match <- observed_names |>
  left_join(
    canonical_lookup |>
      select(
        district_id_exact = .data$district_id,
        canonical_name_exact = .data$canonical_name
      ),
    by = c("admin2_original" = "canonical_name_exact")
  )

# Normalized canonical name.
name_match <- name_match |>
  left_join(
    canonical_lookup |>
      transmute(
        admin2_normalized = .data$canonical_normalized,
        district_id_norm = .data$district_id,
        canonical_name_norm = .data$canonical_name
      ),
    by = "admin2_normalized"
  )

# Explicit aliases.
alias_lookup <- explicit_aliases |>
  select(
    .data$alias,
    .data$alias_normalized,
    district_id_alias = .data$district_id,
    canonical_name_alias = .data$canonical_name
  )

name_match <- name_match |>
  left_join(
    alias_lookup |>
      select(
        admin2_original = .data$alias,
        district_id_alias_exact = .data$district_id_alias,
        canonical_name_alias_exact = .data$canonical_name_alias
      ),
    by = "admin2_original"
  ) |>
  left_join(
    alias_lookup |>
      distinct(
        admin2_normalized = .data$alias_normalized,
        .keep_all = TRUE
      ) |>
      select(
        .data$admin2_normalized,
        district_id_alias_norm = .data$district_id_alias,
        canonical_name_alias_norm = .data$canonical_name_alias
      ),
    by = "admin2_normalized"
  ) |>
  mutate(
    district_id = coalesce(
      .data$district_id_exact,
      .data$district_id_norm,
      .data$district_id_alias_exact,
      .data$district_id_alias_norm
    ),
    canonical_name = coalesce(
      # Exact match canonical name is just the original name.
      if_else(!is.na(.data$district_id_exact), .data$admin2_original,
        NA_character_),
      .data$canonical_name_norm,
      .data$canonical_name_alias_exact,
      .data$canonical_name_alias_norm
    ),
    match_method = case_when(
      !is.na(.data$district_id_exact) ~ "exact_canonical",
      !is.na(.data$district_id_norm) ~ "normalized_exact",
      !is.na(.data$district_id_alias_exact) ~ "explicit_alias_exact",
      !is.na(.data$district_id_alias_norm) ~ "explicit_alias_normalized",
      TRUE ~ "unresolved"
    ),
    confidence = case_when(
      .data$match_method == "exact_canonical" ~ "high",
      .data$match_method == "normalized_exact" ~ "high",
      str_starts(.data$match_method, "explicit_alias") ~ "high",
      TRUE ~ "unresolved"
    ),
    manual_review = .data$match_method == "unresolved"
  ) |>
  select(
    .data$admin2_original,
    .data$admin2_normalized,
    .data$n_events,
    .data$district_id,
    .data$canonical_name,
    .data$match_method,
    .data$confidence,
    .data$manual_review
  )

# Fuzzy suggestions only for unresolved names.
unresolved_names <- name_match |>
  filter(.data$manual_review)

if (nrow(unresolved_names) > 0L) {
  fuzzy_suggestions <- purrr::map_dfr(
    seq_len(nrow(unresolved_names)),
    function(i) {
      src <- unresolved_names$admin2_normalized[i]
      
      d <- stringdist::stringdist(
        src,
        canonical_districts$canonical_normalized,
        method = "jw"
      )
      
      ord <- order(d)
      keep <- ord[seq_len(min(3L, length(ord)))]
      
      tibble(
        admin2_original = unresolved_names$admin2_original[i],
        admin2_normalized = src,
        candidate_rank = seq_along(keep),
        candidate_district_id = canonical_districts$district_id[keep],
        candidate_canonical_name = canonical_districts$canonical_name[keep],
        jaro_winkler_distance = d[keep],
        auto_assignment = FALSE
      )
    }
  )
} else {
  fuzzy_suggestions <- tibble(
    admin2_original = character(),
    admin2_normalized = character(),
    candidate_rank = integer(),
    candidate_district_id = character(),
    candidate_canonical_name = character(),
    jaro_winkler_distance = numeric(),
    auto_assignment = logical()
  )
}

write_qa(name_match, "01_ACLED_name_crosswalk_observed.csv")
write_qa(fuzzy_suggestions, "02_fuzzy_suggestions_MANUAL_REVIEW_ONLY.csv")

# Attach canonical IDs to all events.
acled_h <- acled |>
  left_join(
    name_match |>
      select(
        .data$admin2_original,
        .data$district_id,
        .data$canonical_name,
        .data$match_method,
        .data$confidence,
        .data$manual_review
      ),
    by = "admin2_original"
  )

# 6. Obtain and harmonize a district boundary --------------------------------
boundary <- NULL
boundary_source <- NA_character_
boundary_input <- NA_character_

existing_candidates <- boundary_candidates[file.exists(boundary_candidates)]

if (length(existing_candidates) > 0L) {
  for (p in existing_candidates) {
    candidate <- read_boundary_any(p)
    if (!is.null(candidate)) {
      nm <- detect_name_field(candidate)
      if (!is.na(nm)) {
        boundary <- candidate
        boundary_source <- "local boundary"
        boundary_input <- p
        break
      }
    }
  }
}

# 6A. GADM DIAGNOSTIC only ---------------------------------------------------
# We inspect GADM explicitly because different spelling is not the only
# possible problem: the source may encode a different administrative structure.
# The GADM geometry is never promoted to canonical geometry in this script.
if (isTRUE(run_gadm_diagnostic)) {
  message("\n---------------- GADM ADM2 DIAGNOSTIC ----------------")
  
  gadm_cache <- here::here("data", "raw", "boundaries", "gadm")
  dir.create(gadm_cache, recursive = TRUE, showWarnings = FALSE)
  
  gadm_diag <- tryCatch({
    gadm2 <- geodata::gadm(
      country = "MOZ",
      level = 2,
      path = gadm_cache
    ) |>
      sf::st_as_sf()
    
    if (!all(c("NAME_1", "NAME_2") %in% names(gadm2))) {
      stop("GADM did not return NAME_1/NAME_2.")
    }
    
    gadm_cd <- gadm2 |>
      filter(
        normalize_geo_name(.data$NAME_1) ==
          normalize_geo_name(province)
      )
    
    gadm_names <- gadm_cd |>
      sf::st_drop_geometry() |>
      transmute(
        source = "GADM",
        GID_1 = if ("GID_1" %in% names(gadm_cd)) {
          as.character(.data$GID_1)
        } else {
          NA_character_
        },
        GID_2 = if ("GID_2" %in% names(gadm_cd)) {
          as.character(.data$GID_2)
        } else {
          NA_character_
        },
        NAME_1 = as.character(.data$NAME_1),
        NAME_2 = as.character(.data$NAME_2),
        normalized = normalize_geo_name(.data$NAME_2)
      ) |>
      distinct() |>
      arrange(.data$NAME_2)
    
    gadm_compare <- canonical_districts |>
      select(
        .data$district_id,
        .data$canonical_name,
        .data$canonical_normalized
      ) |>
      full_join(
        gadm_names |>
          select(
            GADM_NAME_2 = .data$NAME_2,
            GADM_GID_2 = .data$GID_2,
            GADM_normalized = .data$normalized
          ),
        by = c("canonical_normalized" = "GADM_normalized")
      ) |>
      mutate(
        GADM_status = case_when(
          is.na(.data$GADM_NAME_2) ~ "MISSING_FROM_GADM_ADM2",
          is.na(.data$canonical_name) ~ "GADM_ONLY",
          .data$canonical_name == .data$GADM_NAME_2 ~ "EXACT_NAME",
          TRUE ~ "NORMALIZED_NAME_MATCH"
        )
      ) |>
      arrange(.data$district_id, .data$GADM_NAME_2)
    
    write_qa(
      gadm_names,
      "00A_GADM_raw_CaboDelgado_ADM2_names.csv"
    )
    
    write_qa(
      gadm_compare,
      "00A_GADM_vs_17_canonical_districts.csv"
    )
    
    missing_names <- gadm_compare |>
      filter(.data$GADM_status == "MISSING_FROM_GADM_ADM2") |>
      pull(.data$canonical_name)
    
    message(
      "GADM Cabo Delgado ADM2 count: ",
      n_distinct(gadm_names$NAME_2)
    )
    message(
      "GADM NAME_2 values: ",
      paste(gadm_names$NAME_2, collapse = " | ")
    )
    
    if (length(missing_names) > 0L) {
      message(
        "Canonical districts absent from GADM ADM2: ",
        paste(missing_names, collapse = " | ")
      )
    }
    
    list(
      names = gadm_names,
      comparison = gadm_compare,
      boundary = gadm_cd
    )
  }, error = function(e) {
    warning("GADM diagnostic failed: ", conditionMessage(e))
    NULL
  })
  
  message("------------------------------------------------------\n")
}

# 6B. CANONICAL GEOMETRY FALLBACK: PUBLIC CENACARTA / MOZGIS -----------------
if (is.null(boundary) && isTRUE(allow_mozgis_fallback)) {
  mozgis <- download_mozgis_districts()
  
  boundary <- mozgis$boundary
  boundary_source <- "CENACARTA/MozGIS public district FeatureServer"
  boundary_input <- mozgis$raw_geojson
}

if (is.null(boundary)) {
  stop(
    "No usable local boundary found and CENACARTA/MozGIS fallback is disabled ",
    "or failed. Provide a current Cabo Delgado Admin-2 boundary before continuing."
  )
}

boundary <- sf::st_make_valid(boundary)
boundary_name_field <- detect_name_field(boundary)

if (is.na(boundary_name_field)) {
  stop("Could not identify the Admin-2 district-name field in the boundary.")
}

boundary_h <- boundary |>
  mutate(
    boundary_name_original = as.character(.data[[boundary_name_field]]),
    boundary_name_normalized = normalize_geo_name(.data$boundary_name_original)
  ) |>
  filter(!is_blankish(.data$boundary_name_original))

# Always expose the actual names used by whichever polygon source was selected.
selected_boundary_names <- boundary_h |>
  sf::st_drop_geometry() |>
  transmute(
    boundary_source = boundary_source,
    source_name = .data$boundary_name_original,
    normalized_name = .data$boundary_name_normalized
  ) |>
  distinct() |>
  arrange(.data$source_name)

write_qa(
  selected_boundary_names,
  "00C_SELECTED_BOUNDARY_raw_district_names.csv"
)

message(
  "Selected boundary contains ",
  n_distinct(selected_boundary_names$source_name),
  " distinct district labels:"
)
message("  ", paste(selected_boundary_names$source_name, collapse = " | "))

# Boundary name reconciliation.
#
# IMPORTANT hierarchy for boundaries:
#   1) source-specific administrative-history rule
#   2) exact canonical name
#   3) normalized canonical name
#   4) generic explicit alias, exact
#   5) generic explicit alias, normalized
#
# Source-specific rules come first because the string "PEMBA" in the selected
# CENACARTA layer refers to the older non-city district geometry that
# corresponds
# to current Metuge, while ACLED's "Pemba" refers to current Pemba.
boundary_names_base <- boundary_h |>
  st_drop_geometry() |>
  tibble::as_tibble() |>
  distinct(.data$boundary_name_original, .data$boundary_name_normalized)

# Activate only rules whose source_pattern matches the actual selected source.
active_source_crosswalk <- boundary_source_crosswalk |>
  filter(stringr::str_detect(boundary_source, fixed(.data$source_pattern))) |>
  select(
    source_name_normalized,
    district_id_source_specific = district_id,
    canonical_name_source_specific = canonical_name,
    source_specific_reason = reason
  ) |>
  distinct(.data$source_name_normalized, .keep_all = TRUE)

# Generic aliases are also matched after normalization to handle uppercase,
# accents and repeated spaces safely.
boundary_alias_norm <- alias_lookup |>
  transmute(
    boundary_alias_normalized = normalize_geo_name(.data$alias),
    district_id_alias_norm = .data$district_id_alias,
    canonical_name_alias_norm = .data$canonical_name_alias
  ) |>
  distinct(.data$boundary_alias_normalized, .keep_all = TRUE)

boundary_name_crosswalk <- boundary_names_base |>
  left_join(
    active_source_crosswalk,
    by = c(
      "boundary_name_normalized" = "source_name_normalized"
    )
  ) |>
  left_join(
    canonical_lookup |>
      transmute(
        boundary_name_original = .data$canonical_name,
        district_id_exact = .data$district_id,
        canonical_name_exact = .data$canonical_name
      ),
    by = "boundary_name_original"
  ) |>
  left_join(
    canonical_lookup |>
      transmute(
        boundary_name_normalized = .data$canonical_normalized,
        district_id_norm = .data$district_id,
        canonical_name_norm = .data$canonical_name
      ),
    by = "boundary_name_normalized"
  ) |>
  left_join(
    alias_lookup |>
      select(
        boundary_name_original = .data$alias,
        district_id_alias_exact = .data$district_id_alias,
        canonical_name_alias_exact = .data$canonical_name_alias
      ),
    by = "boundary_name_original"
  ) |>
  left_join(
    boundary_alias_norm,
    by = c(
      "boundary_name_normalized" = "boundary_alias_normalized"
    )
  ) |>
  mutate(
    district_id = coalesce(
      .data$district_id_source_specific,
      .data$district_id_exact,
      .data$district_id_norm,
      .data$district_id_alias_exact,
      .data$district_id_alias_norm
    ),
    canonical_name = coalesce(
      .data$canonical_name_source_specific,
      .data$canonical_name_exact,
      .data$canonical_name_norm,
      .data$canonical_name_alias_exact,
      .data$canonical_name_alias_norm
    ),
    match_method = case_when(
      !is.na(.data$district_id_source_specific) ~
        "source_specific_admin_history",
      !is.na(.data$district_id_exact) ~
        "exact_canonical",
      !is.na(.data$district_id_norm) ~
        "normalized_exact",
      !is.na(.data$district_id_alias_exact) ~
        "explicit_alias_exact",
      !is.na(.data$district_id_alias_norm) ~
        "explicit_alias_normalized",
      TRUE ~
        "unresolved"
    ),
    confidence = case_when(
      .data$match_method == "source_specific_admin_history" ~ "high",
      .data$match_method %in% c(
        "exact_canonical",
        "normalized_exact",
        "explicit_alias_exact",
        "explicit_alias_normalized"
      ) ~ "high",
      TRUE ~ "unresolved"
    ),
    manual_review = .data$match_method == "unresolved"
  ) |>
  select(
    boundary_name_original,
    boundary_name_normalized,
    district_id,
    canonical_name,
    match_method,
    confidence,
    source_specific_reason,
    manual_review
  ) |>
  arrange(.data$district_id, .data$boundary_name_original)

write_qa(boundary_name_crosswalk, "03_boundary_name_crosswalk.csv")

message("\nBoundary reconciliation:")
print(tibble::as_tibble(boundary_name_crosswalk), n = Inf, width = Inf)

# Strong QA: 17 selected source labels must resolve one-to-one to 17 canonical
# IDs.
boundary_mapping_check <- boundary_name_crosswalk |>
  summarise(
    source_labels = n(),
    resolved_labels = sum(!is.na(.data$district_id)),
    canonical_ids = n_distinct(.data$district_id, na.rm = TRUE),
    duplicated_canonical_ids = sum(
      duplicated(.data$district_id[!is.na(.data$district_id)])
    )
  )

write_qa(
  boundary_mapping_check,
  "03B_boundary_mapping_one_to_one_check.csv"
)

message(
  "Boundary mapping check: source labels=", boundary_mapping_check$source_labels,
  "; resolved=", boundary_mapping_check$resolved_labels,
  "; unique canonical IDs=", boundary_mapping_check$canonical_ids,
  "; duplicated canonical IDs=", boundary_mapping_check$duplicated_canonical_ids
)

if (
  boundary_mapping_check$source_labels == 17L &&
  boundary_mapping_check$resolved_labels == 17L &&
  boundary_mapping_check$canonical_ids == 17L &&
  boundary_mapping_check$duplicated_canonical_ids == 0L
) {
  message(
    "Boundary mapping PASSED: exact one-to-one 17 -> 17 canonical ",
    "districts.")
}

if (any(boundary_name_crosswalk$manual_review)) {
  warning(
    "Some boundary district names remain unresolved after source-specific, ",
    "canonical, normalized, and explicit-alias matching."
  )
}

if (
  boundary_mapping_check$source_labels != 17L ||
  boundary_mapping_check$resolved_labels != 17L ||
  boundary_mapping_check$canonical_ids != 17L ||
  boundary_mapping_check$duplicated_canonical_ids != 0L
) {
  write_qa(
    boundary_name_crosswalk,
    "03C_boundary_mapping_problem_MANUAL_REVIEW.csv"
  )
  
  stop(
    "Selected boundary did not map one-to-one onto the 17 canonical districts. ",
    "Inspect 03_boundary_name_crosswalk.csv and ",
    "03C_boundary_mapping_problem_MANUAL_REVIEW.csv."
  )
}

boundary_h <- boundary_h |>
  left_join(
    boundary_name_crosswalk |>
      select(
        boundary_name_original,
        boundary_district_id = district_id,
        boundary_canonical_name = canonical_name,
        boundary_match_method = match_method,
        boundary_confidence = confidence,
        boundary_source_specific_reason = source_specific_reason,
        boundary_manual_review = manual_review
      ),
    by = "boundary_name_original"
  )

# Do not silently create a partial canonical geography.
if (any(is.na(boundary_h$boundary_district_id))) {
  unresolved_boundary <- boundary_h |>
    st_drop_geometry() |>
    filter(is.na(.data$boundary_district_id))
  
  write_qa(unresolved_boundary,
    "04_unresolved_boundary_names_MANUAL_REVIEW.csv")
  
  stop(
    "Boundary contains unresolved Admin-2 names. Review ",
    file.path(qa_dir, "04_unresolved_boundary_names_MANUAL_REVIEW.csv"),
    " and add validated aliases before continuing."
  )
}

# Dissolve if the source boundary has multiple pieces/features for the same
# district.
canonical_sf <- boundary_h |>
  transmute(
    district_id = .data$boundary_district_id,
    canonical_name = .data$boundary_canonical_name,
    geometry = geometry
  ) |>
  group_by(.data$district_id, .data$canonical_name) |>
  summarise(.groups = "drop")

canonical_sf <- sf::st_make_valid(canonical_sf)

# Require exact 17-district canonical domain.
canonical_missing_from_boundary <- anti_join(
  canonical_districts,
  st_drop_geometry(canonical_sf),
  by = "district_id"
)

boundary_extra <- anti_join(
  st_drop_geometry(canonical_sf),
  canonical_districts,
  by = "district_id"
)

write_qa(canonical_missing_from_boundary,
  "05_canonical_districts_missing_from_boundary.csv")
write_qa(boundary_extra, "06_boundary_districts_not_in_canonical_domain.csv")

if (
  nrow(canonical_sf) != 17L ||
  nrow(canonical_missing_from_boundary) > 0L ||
  nrow(boundary_extra) > 0L
) {
  stop(
    "Canonical boundary does not resolve to exactly the 17 INE districts. ",
    "Inspect QA files 03-06 before building the panel."
  )
}

# Stable ordering independent of source-row order.
canonical_sf <- canonical_sf |>
  mutate(
    canonical_order = match(
      .data$district_id,
      canonical_districts$district_id
    )
  ) |>
  arrange(.data$canonical_order) |>
  select(-.data$canonical_order)

if (!identical(
  as.character(canonical_sf$district_id),
  as.character(canonical_districts$district_id)
)) {
  stop(
    "Canonical boundary ordering could not be aligned to the 17 district ",
    "IDs.")
}

# Geometry diagnostics.
if (is.na(sf::st_crs(canonical_sf))) {
  stop("Canonical boundary has no CRS.")
}

canonical_sf <- st_make_valid(canonical_sf)

canonical_area <- canonical_sf |>
  st_transform(32737) |>
  mutate(area_km2 = as.numeric(sf::st_area(.data$geometry)) / 1e6) |>
  st_drop_geometry() |>
  select(.data$district_id, .data$area_km2)

canonical_sf <- canonical_sf |>
  left_join(canonical_area, by = "district_id") |>
  mutate(
    boundary_source = boundary_source,
    boundary_input = boundary_input
  )

# 7. Spatial validation of ACLED labels against canonical polygons -----------
if (!all(c("longitude", "latitude") %in% names(acled_h))) {
  stop(
    "ACLED data do not contain longitude/latitude required for spatial ",
    "validation.")
}

pts <- acled_h |>
  filter(!is.na(.data$longitude), !is.na(.data$latitude)) |>
  sf::st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

canonical_wgs84 <- st_transform(canonical_sf, 4326)

# st_intersects is deliberate: points on polygon boundaries can map to >1
# polygon.
intersects_list <- st_intersects(pts, canonical_wgs84)

spatial_diag <- tibble(
  event_id_cnty = pts$event_id_cnty,
  n_polygons_intersected = lengths(intersects_list),
  spatial_district_id = vapply(
    intersects_list,
    function(idx) {
      if (length(idx) == 1L) canonical_wgs84$district_id[idx] else NA_character_
    },
    character(1)
  ),
  spatial_canonical_name = vapply(
    intersects_list,
    function(idx) {
      if (length(idx) == 1L) canonical_wgs84$canonical_name[idx] else NA_character_
    },
    character(1)
  )
)

acled_h <- acled_h |>
  left_join(spatial_diag, by = "event_id_cnty") |>
  mutate(
    spatial_status = case_when(
      is.na(.data$spatial_district_id) & .data$n_polygons_intersected == 0L ~ "outside_all_polygons",
      is.na(.data$spatial_district_id) & .data$n_polygons_intersected > 1L ~ "boundary_ambiguous",
      is.na(.data$district_id) ~ "text_unresolved",
      .data$district_id == .data$spatial_district_id ~ "text_spatial_agree",
      TRUE ~ "text_spatial_disagree"
    ),
    spatial_disagreement = .data$spatial_status == "text_spatial_disagree",
    spatial_disagreement_geo1 = .data$spatial_disagreement &
      if ("geo_precision" %in% names(acled_h)) .data$geo_precision == 1L else FALSE,
    needs_manual_review = coalesce(.data$manual_review, FALSE) |
      .data$spatial_status %in% c(
        "outside_all_polygons",
        "boundary_ambiguous",
        "text_unresolved"
      ) |
      .data$spatial_disagreement_geo1
  )

# IMPORTANT: district_id remains based on the harmonized ACLED admin2 label.
# Point-in-polygon is a validation field, not an automatic replacement.

spatial_summary <- if ("geo_precision" %in% names(acled_h)) {
  acled_h |>
    count(.data$spatial_status, .data$geo_precision, name = "n_events",
      sort = TRUE)
} else {
  acled_h |>
    count(.data$spatial_status, name = "n_events", sort = TRUE)
}

write_qa(spatial_summary, "07_spatial_validation_summary.csv")

spatial_disagreements <- acled_h |>
  filter(.data$spatial_status != "text_spatial_agree") |>
  select(
    any_of(c(
      "event_id_cnty", "event_date", "admin2_original", "district_id",
      "canonical_name", "geo_precision", "location", "latitude", "longitude",
      "spatial_district_id", "spatial_canonical_name", "n_polygons_intersected",
      "spatial_status", "notes"
    ))
  )

write_qa(spatial_disagreements, "08_spatial_disagreements_MANUAL_REVIEW.csv")

# 8. Build final crosswalk and manual-review table ---------------------------
acled_crosswalk <- name_match |>
  transmute(
    source = "ACLED",
    source_name = .data$admin2_original,
    normalized_name = .data$admin2_normalized,
    canonical_name = .data$canonical_name,
    district_id = .data$district_id,
    valid_from = analysis_start,
    valid_to = analysis_end,
    match_method = .data$match_method,
    confidence = .data$confidence,
    source_specific_reason = NA_character_,
    manual_review = .data$manual_review,
    n_source_records = .data$n_events
  )

boundary_crosswalk <- boundary_name_crosswalk |>
  transmute(
    source = boundary_source,
    source_name = .data$boundary_name_original,
    normalized_name = .data$boundary_name_normalized,
    canonical_name = .data$canonical_name,
    district_id = .data$district_id,
    valid_from = analysis_start,
    valid_to = analysis_end,
    match_method = .data$match_method,
    confidence = .data$confidence,
    source_specific_reason = .data$source_specific_reason,
    manual_review = .data$manual_review,
    n_source_records = NA_integer_
  )

district_crosswalk <- bind_rows(acled_crosswalk, boundary_crosswalk) |>
  arrange(.data$source, .data$district_id, .data$source_name)

manual_review <- bind_rows(
  acled_h |>
    filter(.data$needs_manual_review) |>
    transmute(
      record_type = "ACLED_event",
      event_id_cnty = .data$event_id_cnty,
      event_date = .data$event_date,
      source_name = .data$admin2_original,
      district_id = .data$district_id,
      canonical_name = .data$canonical_name,
      geo_precision = if ("geo_precision" %in% names(acled_h)) .data$geo_precision else NA_integer_,
      spatial_district_id = .data$spatial_district_id,
      spatial_canonical_name = .data$spatial_canonical_name,
      issue = .data$spatial_status
    ),
  fuzzy_suggestions |>
    filter(.data$candidate_rank == 1L) |>
    transmute(
      record_type = "unresolved_name_fuzzy_suggestion",
      event_id_cnty = NA_character_,
      event_date = as.Date(NA),
      source_name = .data$admin2_original,
      district_id = NA_character_,
      canonical_name = .data$candidate_canonical_name,
      geo_precision = NA_integer_,
      spatial_district_id = NA_character_,
      spatial_canonical_name = NA_character_,
      issue = paste0(
        "fuzzy suggestion only; JW distance=",
        round(.data$jaro_winkler_distance, 4)
      )
    )
) |>
  arrange(.data$record_type, .data$source_name, .data$event_date)

# 9. Hard checks before writing final harmonized objects ---------------------
n_unresolved_events <- sum(is.na(acled_h$district_id))
n_distinct_ids <- n_distinct(acled_h$district_id[!is.na(acled_h$district_id)])

if (n_unresolved_events > 0L) {
  write_qa(
    acled_h |> filter(is.na(.data$district_id)),
    "09_unresolved_ACLED_events_MANUAL_REVIEW.csv"
  )
  stop(
    n_unresolved_events,
    " ACLED events remain without a canonical district_id. ",
    "Resolve aliases before building the district-month panel."
  )
}

if (n_distinct_ids != 17L) {
  warning(
    "ACLED observed only ", n_distinct_ids,
    " canonical districts in this snapshot. This is not necessarily an error; ",
    "the final panel domain still contains all 17 canonical districts."
  )
}

if (nrow(canonical_sf) != 17L) {
  stop("Canonical sf object does not contain exactly 17 districts.")
}

# 10. Save processed outputs -------------------------------------------------
canonical_table <- canonical_districts |>
  select(
    .data$district_id,
    .data$canonical_name,
    .data$valid_from,
    .data$valid_to,
    .data$nominal_reference
  ) |>
  left_join(
    canonical_sf |>
      st_drop_geometry() |>
      select(
        .data$district_id,
        .data$area_km2,
        .data$boundary_source,
        .data$boundary_input
      ),
    by = "district_id"
  )

readr::write_csv(
  district_crosswalk,
  file.path(processed_dir, "cabo_delgado_admin2_crosswalk.csv"),
  na = ""
)

readr::write_csv(
  canonical_table,
  file.path(processed_dir, "cabo_delgado_admin2_canonical.csv"),
  na = ""
)

saveRDS(
  canonical_sf,
  file.path(processed_dir, "cabo_delgado_admin2_canonical.rds"),
  compress = "xz"
)

saveRDS(
  acled_h,
  file.path(processed_dir, "acled_cabo_delgado_harmonized.rds"),
  compress = "xz"
)

readr::write_csv(
  manual_review,
  file.path(processed_dir, "admin2_manual_review.csv"),
  na = ""
)

# QA copies with run timestamp.
write_qa(district_crosswalk, "10_final_admin2_crosswalk.csv")
write_qa(canonical_table, "11_final_canonical_district_table.csv")
write_qa(manual_review, "12_manual_review_final.csv")

harmonization_summary <- tibble(
  metric = c(
    "acled_input",
    "acled_events",
    "canonical_districts",
    "observed_canonical_districts_in_ACLED",
    "unresolved_ACLED_events",
    "spatial_text_agreements",
    "spatial_text_disagreements",
    "spatial_disagreements_geo_precision_1",
    "outside_all_polygons",
    "boundary_source",
    "boundary_input",
    "canonical_boundary_features"
  ),
  value = as.character(c(
    acled_input,
    nrow(acled_h),
    nrow(canonical_districts),
    n_distinct(acled_h$district_id),
    n_unresolved_events,
    sum(acled_h$spatial_status == "text_spatial_agree", na.rm = TRUE),
    sum(acled_h$spatial_status == "text_spatial_disagree", na.rm = TRUE),
    sum(acled_h$spatial_disagreement_geo1, na.rm = TRUE),
    sum(acled_h$spatial_status == "outside_all_polygons", na.rm = TRUE),
    boundary_source,
    boundary_input,
    nrow(canonical_sf)
  ))
)

write_qa(harmonization_summary, "00_HARMONIZATION_SUMMARY.csv")

capture.output(
  sessionInfo(),
  file = file.path(qa_dir, "R_sessionInfo.txt")
)

# 11. Console report ---------------------------------------------------------
message("\nadmin2 harmonization complete")
message("ACLED events                      : ", nrow(acled_h))
message("Canonical district domain         : ", nrow(canonical_districts))
message("Observed district IDs in ACLED    : ", n_distinct(acled_h$district_id))
message("Unresolved ACLED events           : ", n_unresolved_events)
message("Text/spatial agreements           : ",
        sum(acled_h$spatial_status == "text_spatial_agree", na.rm = TRUE))
message("Text/spatial disagreements        : ",
        sum(acled_h$spatial_status == "text_spatial_disagree", na.rm = TRUE))
message("Disagreements with geo_precision=1: ",
        sum(acled_h$spatial_disagreement_geo1, na.rm = TRUE))
message("Events outside all polygons       : ",
        sum(acled_h$spatial_status == "outside_all_polygons", na.rm = TRUE))
message("Boundary source                   : ", boundary_source)
message("Boundary input                    : ", boundary_input)
message("QA directory                      : ", qa_dir)
message("\n")

print(tibble::as_tibble(canonical_table), n = Inf, width = Inf)
print(tibble::as_tibble(spatial_summary), n = Inf, width = Inf)

message(
  "\nIMPORTANT: district_id was assigned from harmonized ACLED admin2 text. ",
  "Spatial containment was used as validation and was NOT allowed to silently ",
  "overwrite the ACLED district assignment.\n"
)
