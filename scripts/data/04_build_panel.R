# Build the district-month panel of events and fatalities.

# 0. Packages ----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  dplyr,
  tidyr,
  readr,
  tibble,
  lubridate,
  stringr,
  purrr,
  sf,
  spdep,
  here
)

# 1. Configuration -----------------------------------------------------------
panel_start <- as.Date("2017-01-01")
panel_end   <- as.Date("2024-12-31")

expected_districts <- 17L
expected_months    <- 96L
expected_cells     <- expected_districts * expected_months

# Spatial graph settings --------------------------------------------------------
# Preserve the pure Queen graph for audit, but create a second connected graph
# for CAR/MRF use by linking only true Queen islands (card == 0) to their
# nearest district using k=1 on projected point-on-surface coordinates.
adjacency_crs_cfg <- 32737L
connect_queen_islands_cfg <- TRUE
island_k_cfg <- 1L

allowed_disorder_types <- c(
  "Political violence",
  "Demonstrations",
  "Political violence; Demonstrations",
  "Strategic developments"
)

dual_disorder_type <- "Political violence; Demonstrations"

allowed_event_types <- c(
  "Battles",
  "Violence against civilians",
  "Explosions/Remote violence",
  "Riots",
  "Protests",
  "Strategic developments"
)

allowed_spatial_status <- c(
  "text_spatial_agree",
  "text_spatial_disagree",
  "outside_all_polygons",
  "boundary_ambiguous",
  "text_unresolved"
)

processed_dir <- here::here("data", "processed")
qa_dir <- here::here("outputs", "district_month_panel")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

event_file <- file.path(
  processed_dir,
  "acled_cabo_delgado_harmonized.rds"
)

canonical_sf_file <- file.path(
  processed_dir,
  "cabo_delgado_admin2_canonical.rds"
)

canonical_table_file <- file.path(
  processed_dir,
  "cabo_delgado_admin2_canonical.csv"
)

# 2. Helpers -----------------------------------------------------------------
is_blankish <- function(x) {
  is.na(x) | stringr::str_trim(as.character(x)) == ""
}

safe_var_to_mean <- function(x) {
  m <- mean(x)
  if (!is.finite(m) || m == 0) {
    return(NA_real_)
  }
  var(x) / m
}

write_qa <- function(x, filename) {
  readr::write_csv(
    x,
    file.path(qa_dir, filename),
    na = ""
  )
  invisible(x)
}

# 3. Input existence and provenance ------------------------------------------
for (p in c(event_file, canonical_sf_file, canonical_table_file)) {
  if (!file.exists(p)) {
    stop(
      "Required input does not exist: ",
      p,
      "\nRun script 03 successfully before script 04."
    )
  }
}

input_md5 <- tools::md5sum(
  c(event_file, canonical_sf_file, canonical_table_file)
)

input_provenance <- tibble(
  input = names(input_md5),
  md5 = unname(input_md5)
)

# 4. Load inputs -------------------------------------------------------------
events <- readRDS(event_file) |>
  tibble::as_tibble()

district_sf <- readRDS(canonical_sf_file)

district_table <- readr::read_csv(
  canonical_table_file,
  show_col_types = FALSE
) |>
  tibble::as_tibble()

if (!inherits(district_sf, "sf")) {
  stop("Canonical district boundary RDS is not an sf object.")
}

if (is.na(sf::st_crs(district_sf))) {
  stop("Canonical district boundary has no CRS.")
}

# 5. Hard preflight: required columns ----------------------------------------
required_event_fields <- c(
  "event_id_cnty",
  "event_date",
  "district_id",
  "canonical_name",
  "fatalities",
  "disorder_type",
  "event_type",
  "sub_event_type",
  "time_precision",
  "geo_precision",
  "spatial_status"
)

missing_event_fields <- setdiff(
  required_event_fields,
  names(events)
)

if (length(missing_event_fields) > 0L) {
  stop(
    "Harmonized ACLED data are missing required fields: ",
    paste(missing_event_fields, collapse = ", ")
  )
}

required_canonical_fields <- c(
  "district_id",
  "canonical_name"
)

missing_table_fields <- setdiff(
  required_canonical_fields,
  names(district_table)
)

if (length(missing_table_fields) > 0L) {
  stop(
    "Canonical district table is missing: ",
    paste(missing_table_fields, collapse = ", ")
  )
}

missing_sf_fields <- setdiff(
  required_canonical_fields,
  names(district_sf)
)

if (length(missing_sf_fields) > 0L) {
  stop(
    "Canonical sf boundary is missing: ",
    paste(missing_sf_fields, collapse = ", ")
  )
}

# 6. Hard preflight: canonical 17-district domain ----------------------------
district_table_core <- district_table |>
  select(district_id, canonical_name) |>
  distinct()

district_sf_core <- district_sf |>
  sf::st_drop_geometry() |>
  tibble::as_tibble() |>
  select(district_id, canonical_name) |>
  distinct()

if (
  nrow(district_table_core) != expected_districts ||
  dplyr::n_distinct(district_table_core$district_id) != expected_districts ||
  dplyr::n_distinct(district_table_core$canonical_name) != expected_districts
) {
  stop(
    "Canonical district table must contain exactly 17 unique district IDs ",
    "and 17 unique canonical names."
  )
}

if (
  nrow(district_sf_core) != expected_districts ||
  dplyr::n_distinct(district_sf_core$district_id) != expected_districts
) {
  stop(
    "Canonical sf boundary must contain exactly 17 unique district IDs."
  )
}

if (!setequal(
  district_table_core$district_id,
  district_sf_core$district_id
)) {
  stop(
    "district_id sets differ between canonical CSV and canonical sf boundary."
  )
}

canonical_name_check <- district_table_core |>
  left_join(
    district_sf_core,
    by = "district_id",
    suffix = c("_table", "_sf")
  ) |>
  mutate(
    names_equal =
      .data$canonical_name_table ==
      .data$canonical_name_sf
  )

if (
  any(is.na(canonical_name_check$names_equal)) ||
  any(!canonical_name_check$names_equal)
) {
  write_qa(
    canonical_name_check,
    "ERROR_canonical_name_mismatch_table_vs_sf.csv"
  )
  
  stop(
    "canonical_name differs between canonical CSV and sf boundary. ",
    "Inspect ERROR_canonical_name_mismatch_table_vs_sf.csv."
  )
}

# Ensure one feature per district.
sf_feature_counts <- district_sf |>
  sf::st_drop_geometry() |>
  count(.data$district_id, name = "n_features")

if (any(sf_feature_counts$n_features != 1L)) {
  write_qa(
    sf_feature_counts,
    "ERROR_canonical_sf_feature_counts.csv"
  )
  
  stop(
    "Canonical sf must contain exactly one feature per district_id."
  )
}

district_sf <- sf::st_make_valid(district_sf)

# 7. Hard preflight: event IDs, dates, fatalities, names, taxonomy -----------
if (nrow(events) == 0L) {
  stop("Harmonized ACLED event file contains zero rows.")
}

if (anyDuplicated(events$event_id_cnty) > 0L) {
  duplicate_ids <- events |>
    count(.data$event_id_cnty, name = "n") |>
    filter(.data$n > 1L)
  
  write_qa(
    duplicate_ids,
    "ERROR_duplicate_event_id_cnty.csv"
  )
  
  stop(
    "Duplicate event_id_cnty found in harmonized ACLED data."
  )
}

# Date parsing must not silently create NA.
events <- events |>
  mutate(
    event_date_raw = as.character(.data$event_date),
    event_date = as.Date(.data$event_date)
  )

if (any(is.na(events$event_date))) {
  bad_dates <- events |>
    filter(is.na(.data$event_date)) |>
    select(event_id_cnty, event_date_raw)
  
  write_qa(
    bad_dates,
    "ERROR_unparseable_event_dates.csv"
  )
  
  stop(
    "At least one ACLED event_date could not be parsed."
  )
}

# Fatalities must be fully observed, finite, non-negative integer counts.
fatalities_raw <- events$fatalities

events <- events |>
  mutate(
    fatalities = suppressWarnings(
      as.numeric(.data$fatalities)
    )
  )

if (any(is.na(events$fatalities))) {
  bad_fatalities <- events |>
    filter(is.na(.data$fatalities)) |>
    select(event_id_cnty, event_date, fatalities)
  
  write_qa(
    bad_fatalities,
    "ERROR_missing_or_non_numeric_fatalities.csv"
  )
  
  stop(
    "Missing/non-numeric fatalities found. ",
    "Do not use na.rm=TRUE to silently turn missing fatalities into zeros."
  )
}

if (any(!is.finite(events$fatalities))) {
  stop("Non-finite fatalities found.")
}

if (any(events$fatalities < 0)) {
  stop("Negative fatalities found.")
}

if (any(abs(events$fatalities - round(events$fatalities)) > 1e-8)) {
  stop("Non-integer fatalities found; investigate before panel construction.")
}

events <- events |>
  mutate(
    fatalities = as.integer(round(.data$fatalities)),
    time_precision = suppressWarnings(
      as.integer(.data$time_precision)
    ),
    geo_precision = suppressWarnings(
      as.integer(.data$geo_precision)
    ),
    district_id = as.character(.data$district_id),
    canonical_name = as.character(.data$canonical_name),
    disorder_type = as.character(.data$disorder_type),
    event_type = as.character(.data$event_type),
    sub_event_type = as.character(.data$sub_event_type),
    spatial_status = as.character(.data$spatial_status)
  )

if (
  any(is_blankish(events$district_id)) ||
  any(is_blankish(events$canonical_name))
) {
  stop("Missing/blank district_id or canonical_name in harmonized ACLED data.")
}

if (
  any(is_blankish(events$disorder_type)) ||
  any(is_blankish(events$event_type)) ||
  any(is_blankish(events$sub_event_type))
) {
  stop(
    "Missing/blank disorder_type, event_type, or sub_event_type in ACLED data."
  )
}

if (
  any(is.na(events$time_precision)) ||
  any(!events$time_precision %in% 1:3)
) {
  stop("time_precision contains missing or values outside {1,2,3}.")
}

if (
  any(is.na(events$geo_precision)) ||
  any(!events$geo_precision %in% 1:3)
) {
  stop("geo_precision contains missing or values outside {1,2,3}.")
}

unknown_disorder_types <- setdiff(
  sort(unique(events$disorder_type)),
  allowed_disorder_types
)

if (length(unknown_disorder_types) > 0L) {
  stop(
    "Unexpected disorder_type value(s): ",
    paste(unknown_disorder_types, collapse = " | ")
  )
}

unknown_event_types <- setdiff(
  sort(unique(events$event_type)),
  allowed_event_types
)

if (length(unknown_event_types) > 0L) {
  stop(
    "Unexpected event_type value(s): ",
    paste(unknown_event_types, collapse = " | ")
  )
}

unknown_spatial_status <- setdiff(
  sort(unique(events$spatial_status)),
  allowed_spatial_status
)

if (length(unknown_spatial_status) > 0L) {
  stop(
    "Unexpected spatial_status value(s): ",
    paste(unknown_spatial_status, collapse = " | ")
  )
}

if (any(events$spatial_status == "text_unresolved")) {
  stop(
    "At least one event still has spatial_status='text_unresolved'. ",
    "Resolve script 03 before building the panel."
  )
}

# 8. Hard preflight: event district IDs/names vs canonical domain ------------
unknown_ids_all <- setdiff(
  unique(events$district_id),
  district_table_core$district_id
)

if (length(unknown_ids_all) > 0L) {
  stop(
    "Event data contain district IDs outside canonical domain: ",
    paste(unknown_ids_all, collapse = ", ")
  )
}

event_name_check <- events |>
  distinct(.data$district_id, .data$canonical_name) |>
  left_join(
    district_table_core,
    by = "district_id",
    suffix = c("_event", "_canonical")
  ) |>
  mutate(
    names_equal =
      .data$canonical_name_event ==
      .data$canonical_name_canonical
  )

if (
  any(is.na(event_name_check$names_equal)) ||
  any(!event_name_check$names_equal)
) {
  write_qa(
    event_name_check,
    "ERROR_event_canonical_name_mismatch.csv"
  )
  
  stop(
    "At least one event has district_id/canonical_name inconsistent ",
    "with the canonical district table."
  )
}

# 9. Restrict to fixed principal analysis window -----------------------------
events_window <- events |>
  filter(
    .data$event_date >= panel_start,
    .data$event_date <= panel_end
  ) |>
  mutate(
    year_month =
      as.Date(lubridate::floor_date(.data$event_date, "month"))
  )

if (nrow(events_window) == 0L) {
  stop("No ACLED events fall inside the 2017-2024 analysis window.")
}

# 9B. ACLED taxonomy consistency for the documented dual disorder category ---
dual_disorder_events <- events_window |>
  filter(.data$disorder_type == dual_disorder_type)

if (nrow(dual_disorder_events) > 0L) {
  bad_dual <- dual_disorder_events |>
    filter(
      .data$event_type != "Protests" |
        .data$sub_event_type != "Excessive force against protesters"
    )
  
  if (nrow(bad_dual) > 0L) {
    write_qa(
      bad_dual,
      "ERROR_dual_disorder_taxonomy_inconsistent.csv"
    )
    
    stop(
      "At least one 'Political violence; Demonstrations' event is not ",
      "Protests / Excessive force against protesters."
    )
  }
}

bad_excessive_force <- events_window |>
  filter(
    .data$sub_event_type == "Excessive force against protesters",
    .data$disorder_type != dual_disorder_type
  )

if (nrow(bad_excessive_force) > 0L) {
  write_qa(
    bad_excessive_force,
    "ERROR_excessive_force_disorder_type_inconsistent.csv"
  )
  
  stop(
    "At least one Excessive force against protesters event does not have ",
    "disorder_type='Political violence; Demonstrations'."
  )
}

# 10. Create complete 17 x 96 domain before joining outcomes -----------------
month_table <- tibble(
  year_month = seq.Date(
    from = lubridate::floor_date(panel_start, "month"),
    to = lubridate::floor_date(panel_end, "month"),
    by = "month"
  )
) |>
  mutate(
    time_id = dplyr::row_number(),
    year = lubridate::year(.data$year_month),
    month = lubridate::month(.data$year_month),
    month_label = format(.data$year_month, "%Y-%m")
  )

if (nrow(month_table) != expected_months) {
  stop(
    "Expected 96 months but generated ",
    nrow(month_table),
    "."
  )
}

district_index <- district_table_core |>
  arrange(.data$district_id) |>
  mutate(
    district_index = dplyr::row_number()
  )

panel_domain <- tidyr::crossing(
  district_index,
  month_table
) |>
  mutate(
    cell_id = paste(
      .data$district_id,
      .data$month_label,
      sep = "__"
    )
  )

if (nrow(panel_domain) != expected_cells) {
  stop(
    "Expected ",
    expected_cells,
    " domain cells but created ",
    nrow(panel_domain),
    "."
  )
}

if (anyDuplicated(panel_domain$cell_id) > 0L) {
  stop("Duplicate cell_id generated in panel domain.")
}

# 11. Event-level helper indicators ------------------------------------------
events_window <- events_window |>
  mutate(
    # ----------------------------------------------------------------------
    # ACLED disorder-type semantics
    # ----------------------------------------------------------------------
    # "Political violence; Demonstrations" is a documented ACLED category
    # for Excessive force against protesters. It belongs to both broader
    # analytical sets. Political violence and Demonstrations therefore overlap
    # when this dual category is included.
    #
    # Mutually exclusive raw disorder categories:
    is_disorder_political_violence_only =
      .data$disorder_type == "Political violence",
    
    is_disorder_demonstrations_only =
      .data$disorder_type == "Demonstrations",
    
    is_disorder_political_violence_demonstrations =
      .data$disorder_type == dual_disorder_type,
    
    is_strategic_development =
      .data$disorder_type == "Strategic developments",
    
    # Broader overlapping analytical sets:
    is_political_violence =
      .data$is_disorder_political_violence_only |
      .data$is_disorder_political_violence_demonstrations,
    
    is_demonstration =
      .data$is_disorder_demonstrations_only |
      .data$is_disorder_political_violence_demonstrations,
    
    # Event types
    is_battle =
      .data$event_type == "Battles",
    
    is_violence_against_civilians =
      .data$event_type == "Violence against civilians",
    
    is_explosions_remote_violence =
      .data$event_type == "Explosions/Remote violence",
    
    is_riot =
      .data$event_type == "Riots",
    
    is_protest =
      .data$event_type == "Protests",
    
    is_strategic_event =
      .data$event_type == "Strategic developments",
    
    # Reported fatalities
    is_zero_reported_fatalities =
      .data$fatalities == 0L,
    
    is_positive_reported_fatalities =
      .data$fatalities > 0L,
    
    # ACLED civilian-targeting field is optional in case the API schema changes.
    is_civilian_targeting =
      if ("civilian_targeting" %in% names(events_window)) {
        !is.na(.data$civilian_targeting) &
          stringr::str_trim(
            as.character(.data$civilian_targeting)
          ) != ""
      } else {
        FALSE
      },
    
    # Precision
    time_precision_1 =
      .data$time_precision == 1L,
    
    time_precision_2 =
      .data$time_precision == 2L,
    
    time_precision_3 =
      .data$time_precision == 3L,
    
    geo_precision_1 =
      .data$geo_precision == 1L,
    
    geo_precision_2 =
      .data$geo_precision == 2L,
    
    geo_precision_3 =
      .data$geo_precision == 3L,
    
    # Spatial QA; does not change district assignment.
    is_spatial_review =
      .data$spatial_status != "text_spatial_agree",
    
    is_text_spatial_disagree =
      .data$spatial_status == "text_spatial_disagree",
    
    is_outside_all_polygons =
      .data$spatial_status == "outside_all_polygons",
    
    is_boundary_ambiguous =
      .data$spatial_status == "boundary_ambiguous",
    
    is_spatial_review_geo1 =
      .data$is_spatial_review &
      .data$geo_precision == 1L
  )

# 12. Aggregate event-derived outcomes by district-month ---------------------
cell_agg <- events_window |>
  group_by(
    .data$district_id,
    .data$year_month
  ) |>
  summarise(
    # --------------------------------------------------------------------------
    # All ACLED events
    # --------------------------------------------------------------------------
    n_events_all = dplyr::n(),
    
    fatalities_all =
      sum(.data$fatalities),
    
    n_zero_fatality_events_all =
      sum(.data$is_zero_reported_fatalities),
    
    n_positive_fatality_events_all =
      sum(.data$is_positive_reported_fatalities),
    
    # --------------------------------------------------------------------------
    # Disorder-type decomposition
    # --------------------------------------------------------------------------
    # Mutually exclusive RAW ACLED disorder categories. These reconstruct all
    # events exactly and are the correct fields for additive decompositions.
    n_disorder_political_violence_only =
      sum(.data$is_disorder_political_violence_only),
    
    fatalities_disorder_political_violence_only =
      sum(
        .data$fatalities[
          .data$is_disorder_political_violence_only
        ]
      ),
    
    n_disorder_demonstrations_only =
      sum(.data$is_disorder_demonstrations_only),
    
    fatalities_disorder_demonstrations_only =
      sum(
        .data$fatalities[
          .data$is_disorder_demonstrations_only
        ]
      ),
    
    n_disorder_political_violence_demonstrations =
      sum(.data$is_disorder_political_violence_demonstrations),
    
    fatalities_disorder_political_violence_demonstrations =
      sum(
        .data$fatalities[
          .data$is_disorder_political_violence_demonstrations
        ]
      ),
    
    n_events_strategic_developments =
      sum(.data$is_strategic_development),
    
    fatalities_strategic_developments =
      sum(
        .data$fatalities[
          .data$is_strategic_development
        ]
      ),
    
    # Broader overlapping analytical sets.
    n_events_political_violence =
      sum(.data$is_political_violence),
    
    fatalities_political_violence =
      sum(
        .data$fatalities[
          .data$is_political_violence
        ]
      ),
    
    n_zero_fatality_events_political_violence =
      sum(
        .data$is_political_violence &
          .data$is_zero_reported_fatalities
      ),
    
    n_positive_fatality_events_political_violence =
      sum(
        .data$is_political_violence &
          .data$is_positive_reported_fatalities
      ),
    
    n_events_demonstrations =
      sum(.data$is_demonstration),
    
    fatalities_demonstrations =
      sum(
        .data$fatalities[
          .data$is_demonstration
        ]
      ),
    
    # --------------------------------------------------------------------------
    # Event-type decomposition
    # --------------------------------------------------------------------------
    n_battles =
      sum(.data$is_battle),
    
    n_violence_against_civilians =
      sum(.data$is_violence_against_civilians),
    
    n_explosions_remote_violence =
      sum(.data$is_explosions_remote_violence),
    
    n_riots =
      sum(.data$is_riot),
    
    n_protests =
      sum(.data$is_protest),
    
    n_strategic_events =
      sum(.data$is_strategic_event),
    
    n_civilian_targeting_events =
      sum(.data$is_civilian_targeting),
    
    # --------------------------------------------------------------------------
    # Temporal/spatial precision composition
    # --------------------------------------------------------------------------
    n_time_precision_1 =
      sum(.data$time_precision_1),
    
    n_time_precision_2 =
      sum(.data$time_precision_2),
    
    n_time_precision_3 =
      sum(.data$time_precision_3),
    
    n_geo_precision_1 =
      sum(.data$geo_precision_1),
    
    n_geo_precision_2 =
      sum(.data$geo_precision_2),
    
    n_geo_precision_3 =
      sum(.data$geo_precision_3),
    
    # --------------------------------------------------------------------------
    # Spatial QA from scripts 03/03b
    # --------------------------------------------------------------------------
    n_spatial_review_events =
      sum(.data$is_spatial_review),
    
    n_text_spatial_disagree =
      sum(.data$is_text_spatial_disagree),
    
    n_outside_all_polygons =
      sum(.data$is_outside_all_polygons),
    
    n_boundary_ambiguous =
      sum(.data$is_boundary_ambiguous),
    
    n_spatial_review_geo1 =
      sum(.data$is_spatial_review_geo1),
    
    .groups = "drop"
  )

# 13. Join to complete domain and create structural panel zeros --------------
panel <- panel_domain |>
  left_join(
    cell_agg,
    by = c(
      "district_id",
      "year_month"
    )
  )

# Only event-derived counts/summed reported fatalities are zero-filled.
# This is not external-covariate imputation.
zero_fill_fields <- c(
  "n_events_all",
  "fatalities_all",
  "n_zero_fatality_events_all",
  "n_positive_fatality_events_all",
  "n_disorder_political_violence_only",
  "fatalities_disorder_political_violence_only",
  "n_disorder_demonstrations_only",
  "fatalities_disorder_demonstrations_only",
  "n_disorder_political_violence_demonstrations",
  "fatalities_disorder_political_violence_demonstrations",
  "n_events_political_violence",
  "fatalities_political_violence",
  "n_zero_fatality_events_political_violence",
  "n_positive_fatality_events_political_violence",
  "n_events_demonstrations",
  "fatalities_demonstrations",
  "n_events_strategic_developments",
  "fatalities_strategic_developments",
  "n_battles",
  "n_violence_against_civilians",
  "n_explosions_remote_violence",
  "n_riots",
  "n_protests",
  "n_strategic_events",
  "n_civilian_targeting_events",
  "n_time_precision_1",
  "n_time_precision_2",
  "n_time_precision_3",
  "n_geo_precision_1",
  "n_geo_precision_2",
  "n_geo_precision_3",
  "n_spatial_review_events",
  "n_text_spatial_disagree",
  "n_outside_all_polygons",
  "n_boundary_ambiguous",
  "n_spatial_review_geo1"
)

panel <- panel |>
  mutate(
    across(
      all_of(zero_fill_fields),
      ~ tidyr::replace_na(.x, 0)
    )
  ) |>
  mutate(
    # --------------------------------------------------------------------------
    # Three-state distinction for all ACLED events
    # --------------------------------------------------------------------------
    no_recorded_event_all =
      .data$n_events_all == 0L,
    
    recorded_event_zero_fatalities_all =
      .data$n_events_all > 0L &
      .data$fatalities_all == 0L,
    
    positive_fatalities_all =
      .data$fatalities_all > 0L,
    
    event_state_all = case_when(
      .data$n_events_all == 0L ~
        "no_recorded_event",
      
      .data$fatalities_all == 0L ~
        "recorded_event_zero_reported_fatalities",
      
      TRUE ~
        "positive_reported_fatalities"
    ),
    
    # --------------------------------------------------------------------------
    # Same distinction restricted to Political violence
    # --------------------------------------------------------------------------
    no_recorded_political_violence =
      .data$n_events_political_violence == 0L,
    
    recorded_political_violence_zero_fatalities =
      .data$n_events_political_violence > 0L &
      .data$fatalities_political_violence == 0L,
    
    positive_political_violence_fatalities =
      .data$fatalities_political_violence > 0L,
    
    event_state_political_violence = case_when(
      .data$n_events_political_violence == 0L ~
        "no_recorded_political_violence",
      
      .data$fatalities_political_violence == 0L ~
        "recorded_political_violence_zero_reported_fatalities",
      
      TRUE ~
        "positive_reported_political_violence_fatalities"
    )
  ) |>
  arrange(
    .data$district_index,
    .data$time_id
  )

# 14. Hard QA: panel domain --------------------------------------------------
if (nrow(panel) != expected_cells) {
  stop(
    "Panel does not have exactly ",
    expected_cells,
    " rows."
  )
}

if (
  dplyr::n_distinct(panel$district_id) !=
  expected_districts
) {
  stop("Panel does not have exactly 17 district IDs.")
}

if (
  dplyr::n_distinct(panel$time_id) !=
  expected_months
) {
  stop("Panel does not have exactly 96 months.")
}

if (anyDuplicated(panel$cell_id) > 0L) {
  stop("Duplicate district-month cell_id found.")
}

month_count_by_district <- panel |>
  count(
    .data$district_id,
    name = "n_months"
  )

if (
  any(
    month_count_by_district$n_months !=
    expected_months
  )
) {
  stop(
    "At least one district does not have exactly 96 months."
  )
}

if (
  any(
    vapply(
      panel[zero_fill_fields],
      function(x) any(is.na(x)),
      logical(1)
    )
  )
) {
  stop(
    "NA remains in at least one event-derived panel field ",
    "after structural zero filling."
  )
}

# 15. Hard QA: decompositions and state identities ---------------------------
# Event-type categories must reconstruct all ACLED events.
event_type_sum <- with(
  panel,
  n_battles +
    n_violence_against_civilians +
    n_explosions_remote_violence +
    n_riots +
    n_protests +
    n_strategic_events
)

if (!all(event_type_sum == panel$n_events_all)) {
  bad_n <- sum(event_type_sum != panel$n_events_all)
  
  stop(
    "Event-type counts fail to reconstruct n_events_all in ",
    bad_n,
    " district-month cells."
  )
}

# Raw disorder_type categories are mutually exclusive and must reconstruct all
# ACLED events exactly.
disorder_type_exclusive_sum <- with(
  panel,
  n_disorder_political_violence_only +
    n_disorder_demonstrations_only +
    n_disorder_political_violence_demonstrations +
    n_events_strategic_developments
)

if (!all(disorder_type_exclusive_sum == panel$n_events_all)) {
  bad_n <- sum(
    disorder_type_exclusive_sum != panel$n_events_all
  )
  
  stop(
    "Mutually exclusive raw disorder_type counts fail to reconstruct ",
    "n_events_all in ",
    bad_n,
    " district-month cells."
  )
}

fatalities_disorder_exclusive_sum <- with(
  panel,
  fatalities_disorder_political_violence_only +
    fatalities_disorder_demonstrations_only +
    fatalities_disorder_political_violence_demonstrations +
    fatalities_strategic_developments
)

if (
  !all(
    fatalities_disorder_exclusive_sum ==
    panel$fatalities_all
  )
) {
  bad_n <- sum(
    fatalities_disorder_exclusive_sum !=
      panel$fatalities_all
  )
  
  stop(
    "Mutually exclusive raw disorder_type fatalities fail to reconstruct ",
    "fatalities_all in ",
    bad_n,
    " district-month cells."
  )
}

# Broad Political violence and Demonstrations overlap through the dual category.
# Inclusion-exclusion must therefore reconstruct all events:
#   PV + Demonstrations + Strategic - Dual = All.
disorder_type_overlap_identity <- with(
  panel,
  n_events_political_violence +
    n_events_demonstrations +
    n_events_strategic_developments -
    n_disorder_political_violence_demonstrations
)

if (
  !all(
    disorder_type_overlap_identity ==
    panel$n_events_all
  )
) {
  bad_n <- sum(
    disorder_type_overlap_identity !=
      panel$n_events_all
  )
  
  stop(
    "Overlap-aware disorder_type identity fails in ",
    bad_n,
    " district-month cells."
  )
}

fatalities_overlap_identity <- with(
  panel,
  fatalities_political_violence +
    fatalities_demonstrations +
    fatalities_strategic_developments -
    fatalities_disorder_political_violence_demonstrations
)

if (
  !all(
    fatalities_overlap_identity ==
    panel$fatalities_all
  )
) {
  bad_n <- sum(
    fatalities_overlap_identity !=
      panel$fatalities_all
  )
  
  stop(
    "Overlap-aware disorder_type fatality identity fails in ",
    bad_n,
    " district-month cells."
  )
}

# Zero + positive event counts must equal total recorded events.
if (
  !all(
    panel$n_zero_fatality_events_all +
    panel$n_positive_fatality_events_all ==
    panel$n_events_all
  )
) {
  stop(
    "Zero/positive-fatality event counts do not reconstruct n_events_all."
  )
}

if (
  !all(
    panel$n_zero_fatality_events_political_violence +
    panel$n_positive_fatality_events_political_violence ==
    panel$n_events_political_violence
  )
) {
  stop(
    "Political-violence zero/positive-fatality event counts ",
    "do not reconstruct n_events_political_violence."
  )
}

# Time and geo precision categories must be exhaustive.
if (
  !all(
    panel$n_time_precision_1 +
    panel$n_time_precision_2 +
    panel$n_time_precision_3 ==
    panel$n_events_all
  )
) {
  stop(
    "time_precision counts do not reconstruct n_events_all."
  )
}

if (
  !all(
    panel$n_geo_precision_1 +
    panel$n_geo_precision_2 +
    panel$n_geo_precision_3 ==
    panel$n_events_all
  )
) {
  stop(
    "geo_precision counts do not reconstruct n_events_all."
  )
}

# Spatial QA categories: review is disagreement + outside + ambiguous.
if (
  !all(
    panel$n_text_spatial_disagree +
    panel$n_outside_all_polygons +
    panel$n_boundary_ambiguous ==
    panel$n_spatial_review_events
  )
) {
  stop(
    "Spatial-review decomposition is inconsistent."
  )
}

# The three outcome states must be mutually exclusive and exhaustive.
all_state_sum <- with(
  panel,
  as.integer(no_recorded_event_all) +
    as.integer(recorded_event_zero_fatalities_all) +
    as.integer(positive_fatalities_all)
)

if (!all(all_state_sum == 1L)) {
  stop(
    "All-event three-state outcome is not mutually exclusive/exhaustive."
  )
}

pv_state_sum <- with(
  panel,
  as.integer(no_recorded_political_violence) +
    as.integer(recorded_political_violence_zero_fatalities) +
    as.integer(positive_political_violence_fatalities)
)

if (!all(pv_state_sum == 1L)) {
  stop(
    "Political-violence three-state outcome is not mutually exclusive/exhaustive."
  )
}

# 16. Hard QA: event-level ↔ panel-level reconciliation ----------------------
event_reconciliation <- tibble(
  metric = c(
    "events_all",
    "reported_fatalities_all",
    "events_political_violence",
    "reported_fatalities_political_violence",
    "events_political_violence_demonstrations_dual",
    "reported_fatalities_political_violence_demonstrations_dual",
    "spatial_review_events",
    "spatial_review_geo_precision_1"
  ),
  
  event_level = c(
    nrow(events_window),
    
    sum(
      events_window$fatalities
    ),
    
    sum(
      events_window$is_political_violence
    ),
    
    sum(
      events_window$fatalities[
        events_window$is_political_violence
      ]
    ),
    
    sum(
      events_window$is_disorder_political_violence_demonstrations
    ),
    
    sum(
      events_window$fatalities[
        events_window$is_disorder_political_violence_demonstrations
      ]
    ),
    
    sum(
      events_window$is_spatial_review
    ),
    
    sum(
      events_window$is_spatial_review_geo1
    )
  ),
  
  panel_level = c(
    sum(panel$n_events_all),
    
    sum(panel$fatalities_all),
    
    sum(panel$n_events_political_violence),
    
    sum(panel$fatalities_political_violence),
    
    sum(panel$n_disorder_political_violence_demonstrations),
    
    sum(panel$fatalities_disorder_political_violence_demonstrations),
    
    sum(panel$n_spatial_review_events),
    
    sum(panel$n_spatial_review_geo1)
  )
) |>
  mutate(
    equal =
      .data$event_level ==
      .data$panel_level
  )

if (!all(event_reconciliation$equal)) {
  write_qa(
    event_reconciliation,
    "ERROR_event_panel_reconciliation.csv"
  )
  
  stop(
    "Event-level and panel-level totals do not reconcile."
  )
}

# 17. Descriptive diagnostics ------------------------------------------------
panel_summary <- tibble(
  metric = c(
    "districts",
    "months",
    "cells",
    "event_rows_2017_2024",
    "cells_no_recorded_event_all",
    "pct_cells_no_recorded_event_all",
    "cells_recorded_event_zero_fatalities_all",
    "pct_cells_recorded_event_zero_fatalities_all",
    "cells_positive_fatalities_all",
    "pct_cells_positive_fatalities_all",
    "cells_no_recorded_political_violence",
    "pct_cells_no_recorded_political_violence",
    "cells_recorded_political_violence_zero_fatalities",
    "pct_cells_recorded_political_violence_zero_fatalities",
    "cells_positive_political_violence_fatalities",
    "pct_cells_positive_political_violence_fatalities",
    "mean_fatalities_all",
    "variance_fatalities_all",
    "variance_to_mean_fatalities_all",
    "max_fatalities_all",
    "mean_fatalities_political_violence",
    "variance_fatalities_political_violence",
    "variance_to_mean_fatalities_political_violence",
    "max_fatalities_political_violence",
    "dual_disorder_events_2017_2024",
    "dual_disorder_reported_fatalities_2017_2024",
    "spatial_review_events_2017_2024",
    "spatial_review_geo_precision_1_2017_2024",
    "district_month_cells_with_spatial_review_event"
  ),
  
  value = c(
    dplyr::n_distinct(panel$district_id),
    dplyr::n_distinct(panel$time_id),
    nrow(panel),
    nrow(events_window),
    
    sum(panel$no_recorded_event_all),
    mean(panel$no_recorded_event_all) * 100,
    
    sum(panel$recorded_event_zero_fatalities_all),
    mean(panel$recorded_event_zero_fatalities_all) * 100,
    
    sum(panel$positive_fatalities_all),
    mean(panel$positive_fatalities_all) * 100,
    
    sum(panel$no_recorded_political_violence),
    mean(panel$no_recorded_political_violence) * 100,
    
    sum(
      panel$recorded_political_violence_zero_fatalities
    ),
    mean(
      panel$recorded_political_violence_zero_fatalities
    ) * 100,
    
    sum(
      panel$positive_political_violence_fatalities
    ),
    mean(
      panel$positive_political_violence_fatalities
    ) * 100,
    
    mean(panel$fatalities_all),
    var(panel$fatalities_all),
    safe_var_to_mean(panel$fatalities_all),
    max(panel$fatalities_all),
    
    mean(
      panel$fatalities_political_violence
    ),
    var(
      panel$fatalities_political_violence
    ),
    safe_var_to_mean(
      panel$fatalities_political_violence
    ),
    max(
      panel$fatalities_political_violence
    ),
    
    sum(
      panel$n_disorder_political_violence_demonstrations
    ),
    
    sum(
      panel$fatalities_disorder_political_violence_demonstrations
    ),
    
    sum(panel$n_spatial_review_events),
    sum(panel$n_spatial_review_geo1),
    sum(panel$n_spatial_review_events > 0L)
  )
)

district_panel_summary <- panel |>
  group_by(
    .data$district_id,
    .data$canonical_name
  ) |>
  summarise(
    months = dplyr::n(),
    
    months_any_event =
      sum(.data$n_events_all > 0L),
    
    months_political_violence =
      sum(
        .data$n_events_political_violence > 0L
      ),
    
    total_events_all =
      sum(.data$n_events_all),
    
    total_fatalities_all =
      sum(.data$fatalities_all),
    
    total_events_political_violence =
      sum(
        .data$n_events_political_violence
      ),
    
    total_fatalities_political_violence =
      sum(
        .data$fatalities_political_violence
      ),
    
    spatial_review_events =
      sum(.data$n_spatial_review_events),
    
    spatial_review_geo1 =
      sum(.data$n_spatial_review_geo1),
    
    .groups = "drop"
  ) |>
  arrange(.data$district_id)

spatial_panel_summary <- panel |>
  filter(.data$n_spatial_review_events > 0L) |>
  select(
    district_id,
    canonical_name,
    year_month,
    time_id,
    n_events_all,
    n_spatial_review_events,
    n_text_spatial_disagree,
    n_outside_all_polygons,
    n_boundary_ambiguous,
    n_spatial_review_geo1
  ) |>
  arrange(
    desc(.data$n_spatial_review_geo1),
    desc(.data$n_spatial_review_events),
    .data$district_id,
    .data$year_month
  )

# 18. Spatial adjacency from canonical polygons ------------------------------
# We preserve TWO graphs:
#
#   nb_queen_raw
#       Pure Queen contiguity. No artificial links.
#
#   nb_queen_connected
#       Queen graph augmented only for true islands (card == 0), using k=1
#       nearest-neighbour links between projected point-on-surface coordinates.
#       This is the graph intended for CAR/MRF models requiring connectivity.
#
# Every artificial bridge is written to QA. No bridge is created silently.

district_sf_model <- district_sf |>
  mutate(
    canonical_order = match(
      .data$district_id,
      district_index$district_id
    )
  ) |>
  arrange(.data$canonical_order) |>
  select(-canonical_order)

if (!identical(
  as.character(district_sf_model$district_id),
  as.character(district_index$district_id)
)) {
  stop(
    "Canonical sf ordering does not align with district_index."
  )
}

# Project before point-on-surface / distance calculations.
district_sf_metric <- sf::st_transform(
  district_sf_model,
  adjacency_crs_cfg
)

# 18.1 Helpers for nb objects ------------------------------------------------
sanitize_nb <- function(nb, region_ids) {
  for (i in seq_along(nb)) {
    nbs <- as.integer(nb[[i]])
    
    # Remove spdep zero sentinel, self-neighbours and duplicates.
    nbs <- nbs[
      is.finite(nbs) &
        nbs > 0L &
        nbs != i
    ]
    
    nb[[i]] <- if (length(nbs) > 0L) {
      as.integer(sort(unique(nbs)))
    } else {
      0L
    }
  }
  
  attr(nb, "region.id") <- as.character(region_ids)
  nb
}

nb_edge_list <- function(nb, district_ids) {
  purrr::map_dfr(
    seq_along(nb),
    function(i) {
      js <- as.integer(nb[[i]])
      js <- js[js > 0L & js != i]
      
      if (length(js) == 0L) {
        return(
          tibble(
            district_id_i = district_ids[i],
            district_id_j = NA_character_
          )
        )
      }
      
      tibble(
        district_id_i = district_ids[i],
        district_id_j = district_ids[js]
      )
    }
  ) |>
    filter(
      is.na(.data$district_id_j) |
        .data$district_id_i < .data$district_id_j
    )
}

nb_components <- function(nb) {
  clean_list <- lapply(
    nb,
    function(js) {
      js <- as.integer(js)
      js[js > 0L]
    }
  )
  
  component_id <- rep(
    NA_integer_,
    length(clean_list)
  )
  
  component_counter <- 0L
  
  for (i in seq_along(clean_list)) {
    if (!is.na(component_id[i])) {
      next
    }
    
    component_counter <- component_counter + 1L
    queue <- i
    component_id[i] <- component_counter
    
    while (length(queue) > 0L) {
      current <- queue[1]
      queue <- queue[-1]
      
      neighbours <- clean_list[[current]]
      
      if (length(neighbours) == 0L) {
        next
      }
      
      new_nodes <- neighbours[
        is.na(component_id[neighbours])
      ]
      
      if (length(new_nodes) > 0L) {
        component_id[new_nodes] <- component_counter
        queue <- c(queue, new_nodes)
      }
    }
  }
  
  component_sizes <- table(component_id)
  
  list(
    component_id = component_id,
    component_size = as.integer(
      component_sizes[
        as.character(component_id)
      ]
    ),
    n_components = length(component_sizes)
  )
}

# 18.2 Pure Queen graph ------------------------------------------------------
nb_queen_raw <- suppressWarnings(
  spdep::poly2nb(
    district_sf_model,
    queen = TRUE,
    row.names = district_sf_model$district_id
  )
)

nb_queen_raw <- sanitize_nb(
  nb_queen_raw,
  district_sf_model$district_id
)

# Make the raw Queen relation explicitly symmetric, then sanitize again.
nb_queen_raw <- spdep::make.sym.nb(
  nb_queen_raw
)

nb_queen_raw <- sanitize_nb(
  nb_queen_raw,
  district_sf_model$district_id
)

raw_card <- spdep::card(
  nb_queen_raw
)

raw_islands <- which(
  raw_card == 0L
)

raw_components <- nb_components(
  nb_queen_raw
)

raw_edge_list <- nb_edge_list(
  nb_queen_raw,
  district_sf_model$district_id
)

# 18.3 Diagnostic geometry for islands ---------------------------------------
representative_points <- suppressWarnings(
  sf::st_point_on_surface(
    district_sf_metric
  )
)

representative_xy <- sf::st_coordinates(
  representative_points
)

if (
  nrow(representative_xy) !=
  nrow(district_sf_model)
) {
  stop(
    "Could not obtain exactly one representative point per district."
  )
}

if (length(raw_islands) > 0L) {
  raw_island_distance_qa <- purrr::map_dfr(
    raw_islands,
    function(i) {
      other_index <- setdiff(
        seq_len(nrow(district_sf_metric)),
        i
      )
      
      polygon_distances_m <- as.numeric(
        sf::st_distance(
          district_sf_metric[i, ],
          district_sf_metric[other_index, ]
        )
      )
      
      point_distances_m <- sqrt(
        (
          representative_xy[other_index, 1] -
            representative_xy[i, 1]
        )^2 +
          (
            representative_xy[other_index, 2] -
              representative_xy[i, 2]
          )^2
      )
      
      nearest_polygon_pos <- which.min(
        polygon_distances_m
      )
      
      nearest_point_pos <- which.min(
        point_distances_m
      )
      
      nearest_polygon_index <-
        other_index[nearest_polygon_pos]
      
      nearest_point_index <-
        other_index[nearest_point_pos]
      
      tibble(
        island_index = i,
        
        district_id =
          district_sf_model$district_id[i],
        
        canonical_name =
          district_sf_model$canonical_name[i],
        
        nearest_by_polygon_district_id =
          district_sf_model$district_id[
            nearest_polygon_index
          ],
        
        nearest_by_polygon_name =
          district_sf_model$canonical_name[
            nearest_polygon_index
          ],
        
        nearest_polygon_distance_m =
          round(
            polygon_distances_m[
              nearest_polygon_pos
            ],
            1
          ),
        
        nearest_by_point_district_id =
          district_sf_model$district_id[
            nearest_point_index
          ],
        
        nearest_by_point_name =
          district_sf_model$canonical_name[
            nearest_point_index
          ],
        
        nearest_point_surface_distance_m =
          round(
            point_distances_m[
              nearest_point_pos
            ],
            1
          )
      )
    }
  )
} else {
  raw_island_distance_qa <- tibble(
    island_index = integer(),
    district_id = character(),
    canonical_name = character(),
    nearest_by_polygon_district_id = character(),
    nearest_by_polygon_name = character(),
    nearest_polygon_distance_m = numeric(),
    nearest_by_point_district_id = character(),
    nearest_by_point_name = character(),
    nearest_point_surface_distance_m = numeric()
  )
}

# 18.4 Augment only true islands using k=1 nearest neighbour -----------------
nb_queen_connected <- nb_queen_raw

bridge_edges <- tibble(
  island_district_id = character(),
  island_name = character(),
  connected_to_district_id = character(),
  connected_to_name = character(),
  method = character(),
  point_surface_distance_m = numeric(),
  polygon_distance_m = numeric(),
  nearest_by_polygon_district_id = character(),
  bridge_matches_nearest_polygon = logical()
)

if (
  isTRUE(connect_queen_islands_cfg) &&
  length(raw_islands) > 0L
) {
  # We only need the nearest-neighbour indices for the Queen islands.
  # Do not convert the global k=1 relation to an nb object here: knn2nb()
  # legitimately warns that a k=1 graph over all 17 districts is disconnected,
  # although that auxiliary graph is not the model graph and is used only to
  # retrieve each island's closest projected representative point.
  nearest_knn <- spdep::knearneigh(
    representative_xy,
    k = island_k_cfg
  )
  
  for (i in raw_islands) {
    candidates <- as.integer(
      nearest_knn$nn[i, ]
    )
    
    candidates <- candidates[
      candidates > 0L &
        candidates != i
    ]
    
    if (length(candidates) == 0L) {
      stop(
        "No valid nearest neighbour found for Queen island ",
        district_sf_model$district_id[i],
        "."
      )
    }
    
    # k=1 by design.
    j <- candidates[1]
    
    # Record diagnostic distances before adding the edge.
    point_distance_m <- sqrt(
      (
        representative_xy[j, 1] -
          representative_xy[i, 1]
      )^2 +
        (
          representative_xy[j, 2] -
            representative_xy[i, 2]
        )^2
    )
    
    polygon_distance_m <- as.numeric(
      sf::st_distance(
        district_sf_metric[i, ],
        district_sf_metric[j, ]
      )
    )
    
    nb_queen_connected[[i]] <- as.integer(
      sort(
        unique(
          c(
            nb_queen_connected[[i]][
              nb_queen_connected[[i]] > 0L
            ],
            j
          )
        )
      )
    )
    
    # Add reverse edge immediately; make.sym.nb below is a second guarantee.
    nb_queen_connected[[j]] <- as.integer(
      sort(
        unique(
          c(
            nb_queen_connected[[j]][
              nb_queen_connected[[j]] > 0L
            ],
            i
          )
        )
      )
    )
    
    bridge_edges <- bind_rows(
      bridge_edges,
      tibble(
        island_district_id =
          district_sf_model$district_id[i],
        
        island_name =
          district_sf_model$canonical_name[i],
        
        connected_to_district_id =
          district_sf_model$district_id[j],
        
        connected_to_name =
          district_sf_model$canonical_name[j],
        
        method =
          paste0(
            "Queen island + k=",
            island_k_cfg,
            " nearest projected point-on-surface"
          ),
        
        point_surface_distance_m =
          round(point_distance_m, 1),
        
        polygon_distance_m =
          round(polygon_distance_m, 1),
        
        nearest_by_polygon_district_id =
          raw_island_distance_qa$nearest_by_polygon_district_id[
            match(
              district_sf_model$district_id[i],
              raw_island_distance_qa$district_id
            )
          ],
        
        bridge_matches_nearest_polygon =
          district_sf_model$district_id[j] ==
          raw_island_distance_qa$nearest_by_polygon_district_id[
            match(
              district_sf_model$district_id[i],
              raw_island_distance_qa$district_id
            )
          ]
      )
    )
  }
}

if (
  nrow(bridge_edges) > 0L &&
  any(!bridge_edges$bridge_matches_nearest_polygon)
) {
  warning(
    "At least one island bridge selected by projected point-on-surface kNN ",
    "does not match the nearest polygon by boundary-to-boundary distance. ",
    "Inspect 07C_island_bridge_edges.csv before spatial modelling."
  )
}

# Force symmetry after augmentation and clean one final time.
nb_queen_connected <- spdep::make.sym.nb(
  nb_queen_connected
)

nb_queen_connected <- sanitize_nb(
  nb_queen_connected,
  district_sf_model$district_id
)

connected_card <- spdep::card(
  nb_queen_connected
)

connected_islands <- which(
  connected_card == 0L
)

connected_components <- nb_components(
  nb_queen_connected
)

connected_edge_list <- nb_edge_list(
  nb_queen_connected,
  district_sf_model$district_id
)

# 18.5 Hard QA for the model graph -------------------------------------------
if (
  isTRUE(connect_queen_islands_cfg) &&
  length(connected_islands) > 0L
) {
  stop(
    "Island augmentation failed: ",
    length(connected_islands),
    " district(s) still have zero neighbours."
  )
}

# Connecting true islands should normally yield one connected component here.
# If it does not, there is a disconnected multi-district component that should
# be investigated explicitly rather than bridged automatically.
if (
  isTRUE(connect_queen_islands_cfg) &&
  connected_components$n_components > 1L
) {
  stop(
    "The augmented adjacency graph still has ",
    connected_components$n_components,
    " connected components. ",
    "Do not add further automatic bridges; inspect the boundary/components."
  )
}

# 18.6 Summaries -------------------------------------------------------------
adjacency_summary <- district_index |>
  mutate(
    n_queen_neighbors_raw =
      raw_card,
    
    is_queen_island_raw =
      .data$n_queen_neighbors_raw == 0L,
    
    queen_component_id_raw =
      raw_components$component_id,
    
    queen_component_size_raw =
      raw_components$component_size,
    
    n_neighbors_model =
      connected_card,
    
    is_isolate_model =
      .data$n_neighbors_model == 0L,
    
    model_component_id =
      connected_components$component_id,
    
    model_component_size =
      connected_components$component_size
  )

raw_component_summary <- adjacency_summary |>
  group_by(
    .data$queen_component_id_raw
  ) |>
  summarise(
    n_districts = dplyr::n(),
    
    district_ids =
      paste(
        .data$district_id,
        collapse = "; "
      ),
    
    district_names =
      paste(
        .data$canonical_name,
        collapse = "; "
      ),
    
    .groups = "drop"
  ) |>
  arrange(
    desc(.data$n_districts),
    .data$queen_component_id_raw
  )

model_component_summary <- adjacency_summary |>
  group_by(
    .data$model_component_id
  ) |>
  summarise(
    n_districts = dplyr::n(),
    
    district_ids =
      paste(
        .data$district_id,
        collapse = "; "
      ),
    
    district_names =
      paste(
        .data$canonical_name,
        collapse = "; "
      ),
    
    .groups = "drop"
  ) |>
  arrange(
    desc(.data$n_districts),
    .data$model_component_id
  )

message(
  "Pure Queen components               : ",
  raw_components$n_components,
  " [sizes: ",
  paste(
    raw_component_summary$n_districts,
    collapse = ", "
  ),
  "]"
)

message(
  "Pure Queen islands                  : ",
  length(raw_islands)
)

if (length(raw_islands) > 0L) {
  message(
    "Queen island district(s)            : ",
    paste(
      adjacency_summary$canonical_name[
        adjacency_summary$is_queen_island_raw
      ],
      collapse = ", "
    )
  )
}

message(
  "Artificial island bridges added     : ",
  nrow(bridge_edges)
)

if (nrow(bridge_edges) > 0L) {
  message(
    "Bridge(s)                           : ",
    paste(
      paste0(
        bridge_edges$island_name,
        " <-> ",
        bridge_edges$connected_to_name,
        " [",
        bridge_edges$point_surface_distance_m,
        " m representative-point distance; ",
        bridge_edges$polygon_distance_m,
        " m polygon distance]"
      ),
      collapse = " | "
    )
  )
}

message(
  "Model graph connected components    : ",
  connected_components$n_components
)

# Backward-compatible alias:
# nb_queen means the connected model graph from this version onward.
nb_queen <- nb_queen_connected
adjacency_edge_list <- connected_edge_list

# 19. Data dictionary --------------------------------------------------------
panel_dictionary <- tibble::tribble(
  ~variable, ~level, ~definition,
  
  "district_id",
  "district",
  "Stable canonical Cabo Delgado district identifier",
  
  "canonical_name",
  "district",
  "Canonical district label",
  
  "district_index",
  "district",
  "Integer spatial index 1..17",
  
  "year_month",
  "month",
  "First day of calendar month; identifies month",
  
  "time_id",
  "month",
  "Integer temporal index 1..96",
  
  "cell_id",
  "district-month",
  "Unique district_id + YYYY-MM key",
  
  "n_events_all",
  "district-month",
  "Count of all ACLED events",
  
  "fatalities_all",
  "district-month",
  "Sum of ACLED reported fatalities across all events",
  
  "n_zero_fatality_events_all",
  "district-month",
  "Count of recorded events whose ACLED fatalities field equals zero",
  
  "n_positive_fatality_events_all",
  "district-month",
  "Count of recorded events with ACLED fatalities > 0",
  
  "n_disorder_political_violence_only",
  "district-month",
  "Count with raw disorder_type exactly Political violence",
  
  "fatalities_disorder_political_violence_only",
  "district-month",
  "Reported fatalities for raw disorder_type exactly Political violence",
  
  "n_disorder_demonstrations_only",
  "district-month",
  "Count with raw disorder_type exactly Demonstrations",
  
  "fatalities_disorder_demonstrations_only",
  "district-month",
  "Reported fatalities for raw disorder_type exactly Demonstrations",
  
  "n_disorder_political_violence_demonstrations",
  "district-month",
  "Count with raw disorder_type exactly Political violence; Demonstrations",
  
  "fatalities_disorder_political_violence_demonstrations",
  "district-month",
  "Reported fatalities for raw dual disorder_type Political violence; Demonstrations",
  
  "n_events_political_violence",
  "district-month",
  "Broad Political violence count: Political violence only + dual Political violence; Demonstrations",
  
  "fatalities_political_violence",
  "district-month",
  "Reported fatalities in broad Political violence, including dual events",
  
  "n_events_demonstrations",
  "district-month",
  "Broad Demonstrations count: Demonstrations only + dual Political violence; Demonstrations",
  
  "fatalities_demonstrations",
  "district-month",
  "Reported fatalities in broad Demonstrations, including dual events",
  
  "n_events_strategic_developments",
  "district-month",
  "Count with disorder_type == Strategic developments",
  
  "fatalities_strategic_developments",
  "district-month",
  "Sum of reported fatalities for Strategic developments",
  
  "n_battles",
  "district-month",
  "Count of Battles",
  
  "n_violence_against_civilians",
  "district-month",
  "Count of Violence against civilians",
  
  "n_explosions_remote_violence",
  "district-month",
  "Count of Explosions/Remote violence",
  
  "n_riots",
  "district-month",
  "Count of Riots",
  
  "n_protests",
  "district-month",
  "Count of Protests",
  
  "n_strategic_events",
  "district-month",
  "Count of Strategic developments event_type",
  
  "n_civilian_targeting_events",
  "district-month",
  "Count of events with nonblank ACLED civilian_targeting field",
  
  "n_time_precision_1",
  "district-month",
  "Count of events with ACLED time_precision == 1",
  
  "n_time_precision_2",
  "district-month",
  "Count of events with ACLED time_precision == 2",
  
  "n_time_precision_3",
  "district-month",
  "Count of events with ACLED time_precision == 3",
  
  "n_geo_precision_1",
  "district-month",
  "Count of events with ACLED geo_precision == 1",
  
  "n_geo_precision_2",
  "district-month",
  "Count of events with ACLED geo_precision == 2",
  
  "n_geo_precision_3",
  "district-month",
  "Count of events with ACLED geo_precision == 3",
  
  "n_spatial_review_events",
  "district-month",
  "Count of events not classified text_spatial_agree in spatial QA",
  
  "n_text_spatial_disagree",
  "district-month",
  "Count whose point falls inside another canonical district polygon",
  
  "n_outside_all_polygons",
  "district-month",
  "Count whose event point falls outside all canonical polygons",
  
  "n_boundary_ambiguous",
  "district-month",
  "Count intersecting more than one polygon/boundary in QA",
  
  "n_spatial_review_geo1",
  "district-month",
  "Spatial-review event count restricted to geo_precision == 1",
  
  "no_recorded_event_all",
  "district-month",
  "TRUE when no ACLED event was recorded",
  
  "recorded_event_zero_fatalities_all",
  "district-month",
  "TRUE when >=1 event recorded but summed reported fatalities are zero",
  
  "positive_fatalities_all",
  "district-month",
  "TRUE when summed reported fatalities > 0",
  
  "event_state_all",
  "district-month",
  "Three-state all-event occurrence/fatality indicator",
  
  "no_recorded_political_violence",
  "district-month",
  "TRUE when no Political violence event was recorded",
  
  "recorded_political_violence_zero_fatalities",
  "district-month",
  "TRUE when Political violence recorded but summed reported fatalities are zero",
  
  "positive_political_violence_fatalities",
  "district-month",
  "TRUE when Political-violence reported fatalities > 0",
  
  "event_state_political_violence",
  "district-month",
  "Three-state Political-violence occurrence/fatality indicator"
)

# 20. Save only after all hard QA has passed ---------------------------------
event_level_rds <- file.path(
  processed_dir,
  "cabo_delgado_event_level_2017_2024.rds"
)

panel_rds <- file.path(
  processed_dir,
  "cabo_delgado_district_month_panel_2017_2024.rds"
)

panel_csv <- file.path(
  processed_dir,
  "cabo_delgado_district_month_panel_2017_2024.csv"
)

adjacency_rds <- file.path(
  processed_dir,
  "cabo_delgado_adjacency_topological.rds"
)

dictionary_csv <- file.path(
  processed_dir,
  "cabo_delgado_district_month_panel_dictionary.csv"
)

saveRDS(
  events_window,
  event_level_rds,
  compress = "xz"
)

saveRDS(
  panel,
  panel_rds,
  compress = "xz"
)

readr::write_csv(
  panel,
  panel_csv,
  na = ""
)

saveRDS(
  list(
    district_index =
      district_index,
    
    # Pure contiguity graph for audit/sensitivity.
    nb_queen_raw =
      nb_queen_raw,
    
    raw_edge_list =
      raw_edge_list,
    
    raw_component_summary =
      raw_component_summary,
    
    raw_island_distance_qa =
      raw_island_distance_qa,
    
    # Connected graph intended for CAR/MRF use.
    nb_queen_connected =
      nb_queen_connected,
    
    connected_edge_list =
      connected_edge_list,
    
    model_component_summary =
      model_component_summary,
    
    bridge_edges =
      bridge_edges,
    
    # Backward-compatible aliases point to the connected model graph.
    nb_queen =
      nb_queen_connected,
    
    adjacency_edge_list =
      connected_edge_list,
    
    adjacency_summary =
      adjacency_summary,
    
    settings = list(
      queen = TRUE,
      adjacency_crs = adjacency_crs_cfg,
      connect_queen_islands = connect_queen_islands_cfg,
      island_k = island_k_cfg,
      bridge_method =
        "knearneigh on projected point-on-surface coordinates, islands only"
    ),
    
    note = paste(
      "Pure Queen graph preserved separately.",
      "Only zero-neighbour Queen islands are bridged in the model graph.",
      "Bridges are symmetric, use projected point-on-surface k=1,",
      "and are recorded explicitly in bridge_edges.",
      "No automatic bridges are added between non-island disconnected components."
    )
  ),
  adjacency_rds,
  compress = "xz"
)

readr::write_csv(
  panel_dictionary,
  dictionary_csv,
  na = ""
)

# 21. Save QA/provenance -----------------------------------------------------
write_qa(
  input_provenance,
  "00_input_provenance_md5.csv"
)

build_metadata <- tibble(
  field = c(
    "panel_start",
    "panel_end",
    "districts",
    "months",
    "cells",
    "event_rows_2017_2024",
    "reported_fatalities_2017_2024",
    "spatial_review_events_2017_2024",
    "spatial_review_geo1_2017_2024"
  ),
  
  value = as.character(
    c(
      panel_start,
      panel_end,
      expected_districts,
      expected_months,
      expected_cells,
      nrow(events_window),
      sum(events_window$fatalities),
      sum(events_window$is_spatial_review),
      sum(events_window$is_spatial_review_geo1)
    )
  )
)

write_qa(
  build_metadata,
  "00_panel_build_metadata.csv"
)

write_qa(
  month_table,
  "01_month_index.csv"
)

write_qa(
  district_index,
  "02_district_index.csv"
)

write_qa(
  event_reconciliation,
  "03_event_panel_reconciliation.csv"
)

write_qa(
  panel_summary,
  "04_panel_outcome_summary.csv"
)

write_qa(
  district_panel_summary,
  "05_panel_summary_by_district.csv"
)

write_qa(
  adjacency_summary,
  "06_adjacency_summary_raw_and_model.csv"
)

write_qa(
  raw_edge_list,
  "07A_queen_raw_edges.csv"
)

write_qa(
  connected_edge_list,
  "07B_model_connected_edges.csv"
)

write_qa(
  bridge_edges,
  "07C_island_bridge_edges.csv"
)

write_qa(
  raw_component_summary,
  "07D_queen_raw_components.csv"
)

write_qa(
  model_component_summary,
  "07E_model_connected_components.csv"
)

write_qa(
  raw_island_distance_qa,
  "07F_queen_island_distance_diagnostics.csv"
)

fatality_frequency_all <- panel |>
  count(
    .data$fatalities_all,
    name = "cells"
  ) |>
  arrange(.data$fatalities_all)

fatality_frequency_pv <- panel |>
  count(
    .data$fatalities_political_violence,
    name = "cells"
  ) |>
  arrange(.data$fatalities_political_violence)

write_qa(
  fatality_frequency_all,
  "08_fatalities_all_exact_frequency.csv"
)

write_qa(
  fatality_frequency_pv,
  "09_fatalities_political_violence_exact_frequency.csv"
)

write_qa(
  spatial_panel_summary,
  "10_panel_cells_with_spatial_review_events.csv"
)

taxonomy_disorder <- events_window |>
  count(
    .data$disorder_type,
    name = "n_events",
    sort = TRUE
  )

taxonomy_event_type <- events_window |>
  count(
    .data$event_type,
    name = "n_events",
    sort = TRUE
  )

write_qa(
  taxonomy_disorder,
  "11_disorder_type_counts_2017_2024.csv"
)

write_qa(
  taxonomy_event_type,
  "12_event_type_counts_2017_2024.csv"
)

dual_disorder_qa <- events_window |>
  filter(.data$disorder_type == dual_disorder_type) |>
  select(
    any_of(
      c(
        "event_id_cnty",
        "event_date",
        "district_id",
        "canonical_name",
        "event_type",
        "sub_event_type",
        "fatalities",
        "location",
        "geo_precision",
        "time_precision"
      )
    )
  ) |>
  arrange(.data$event_date, .data$event_id_cnty)

write_qa(
  dual_disorder_qa,
  "13_dual_disorder_political_violence_demonstrations.csv"
)

capture.output(
  sessionInfo(),
  file = file.path(
    qa_dir,
    "R_sessionInfo.txt"
  )
)

# 22. Console report ---------------------------------------------------------
message(
  "\ndistrict-month panel complete"
)

message(
  "Panel period                         : ",
  panel_start,
  " to ",
  panel_end
)

message(
  "Canonical districts                  : ",
  dplyr::n_distinct(panel$district_id)
)

message(
  "Months                               : ",
  dplyr::n_distinct(panel$time_id)
)

message(
  "Expected cells                       : ",
  expected_cells
)

message(
  "Actual cells                         : ",
  nrow(panel)
)

message(
  "Event rows in 2017-2024              : ",
  nrow(events_window)
)

message(
  "Spatial-review events 2017-2024      : ",
  sum(panel$n_spatial_review_events)
)

message(
  "Spatial-review geo_precision=1       : ",
  sum(panel$n_spatial_review_geo1)
)

message(
  "All-event no-event cells             : ",
  sum(panel$n_events_all == 0L),
  " / ",
  nrow(panel),
  " (",
  round(
    mean(panel$n_events_all == 0L) * 100,
    2
  ),
  "%)"
)

message(
  "Fatalities-all zero cells            : ",
  sum(panel$fatalities_all == 0L),
  " / ",
  nrow(panel),
  " (",
  round(
    mean(panel$fatalities_all == 0L) * 100,
    2
  ),
  "%)"
)

message(
  "Political-violence no-event cells    : ",
  sum(
    panel$n_events_political_violence == 0L
  ),
  " / ",
  nrow(panel),
  " (",
  round(
    mean(
      panel$n_events_political_violence == 0L
    ) * 100,
    2
  ),
  "%)"
)

message(
  "PV-fatalities zero cells             : ",
  sum(
    panel$fatalities_political_violence == 0L
  ),
  " / ",
  nrow(panel),
  " (",
  round(
    mean(
      panel$fatalities_political_violence == 0L
    ) * 100,
    2
  ),
  "%)"
)

message(
  "Pure Queen islands                  : ",
  sum(adjacency_summary$is_queen_island_raw)
)

message(
  "Artificial island bridges           : ",
  nrow(bridge_edges)
)

message(
  "Model graph remaining isolates       : ",
  sum(adjacency_summary$is_isolate_model)
)

message(
  "Model graph connected components     : ",
  connected_components$n_components
)

message(
  "Panel RDS                            : ",
  panel_rds
)

message(
  "Panel CSV                            : ",
  panel_csv
)

message(
  "QA directory                         : ",
  qa_dir
)

message(
  "\n"
)

message("\n--- Event-to-panel reconciliation ---\n")

print(
  tibble::as_tibble(event_reconciliation),
  n = Inf,
  width = Inf
)

message("\n--- Panel outcome summary ---\n")

print(
  tibble::as_tibble(panel_summary),
  n = Inf,
  width = Inf
)

message("\n--- District summary ---\n")

print(
  tibble::as_tibble(district_panel_summary),
  n = Inf,
  width = Inf
)

message("\n--- Topological adjacency summary ---\n")

print(
  tibble::as_tibble(adjacency_summary),
  n = Inf,
  width = Inf
)

message(
  "\nIMPORTANT:\n",
  "1) No external predictor/covariate was merged or imputed.\n",
  "2) Structural zero filling applies only after the complete 17 x 96 domain ",
  "was created.\n",
  "3) Spatial-review events were retained under harmonized ACLED district_id; ",
  "no point-in-polygon reassignment was performed.\n",
  "4) Multiple outcome definitions remain available intentionally.\n",
  "5) ACLED 'Political violence; Demonstrations' is retained as a dual raw ",
  "category. Broad Political violence and Demonstrations overlap and must not ",
  "be added without subtracting the dual category.\n",
  "6) The pure Queen graph is preserved unchanged for audit/sensitivity.\n",
  "7) For the CAR/MRF model graph, only true Queen islands are linked to ",
  "their k=1 nearest projected point-on-surface neighbour; every bridge is ",
  "stored explicitly in bridge_edges and QA outputs.\n",
  "8) No automatic bridge is added between non-island disconnected components.\n"
)
