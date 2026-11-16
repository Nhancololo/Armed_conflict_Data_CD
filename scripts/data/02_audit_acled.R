# Audit the ACLED records for Cabo Delgado.

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
  lubridate,
  tibble,
  purrr,
  ggplot2,
  scales,
  here,
  fs,
  sf,
  digest
)

# 1. Configuration -----------------------------------------------------------
province_cfg <- "Cabo Delgado"
country_cfg  <- "Mozambique"

# Intended principal analysis window. This is only an audit label here.
main_start <- as.Date("2017-01-01")
main_end   <- as.Date("2024-12-31")

data_root <- here::here("data", "raw", "acled")
audit_dir_cfg <- here::here("outputs", "acled_audit")
dir.create(audit_dir_cfg, recursive = TRUE, showWarnings = FALSE)

# Optional known boundary files. The first valid sf object found will be used
# only for DIAGNOSTIC comparison; it will not define the study domain here.
boundary_candidates <- c(
  here::here("data", "processed", "cabo_delgado_boundary.rds"),
  here::here("data", "processed", "cabo_delgado_boundary2.rds"),
  here::here("cabo_delgado_boundary.rds"),
  here::here("cabo_delgado_boundary2.rds")
)

# 2. Helpers -----------------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

is_blankish <- function(x) {
  if (is.character(x) || is.factor(x)) {
    is.na(x) | stringr::str_trim(as.character(x)) == ""
  } else {
    is.na(x)
  }
}

safe_pct <- function(num, den) {
  ifelse(den > 0, 100 * num / den, NA_real_)
}

normalize_geo_name <- function(x) {
  x |>
    as.character() |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", " ") |>
    stringr::str_squish()
}

write_csv_audit <- function(x, filename) {
  readr::write_csv(x, file.path(audit_dir_cfg, filename), na = "")
  invisible(x)
}

count_semicolon_items <- function(x) {
  x <- as.character(x)
  blank <- is.na(x) | stringr::str_trim(x) == ""
  out <- stringr::str_count(x, fixed(";")) + 1L
  out[blank] <- 0L
  out
}

explode_semicolon_counts <- function(data, var, out_name = "item") {
  if (!var %in% names(data)) return(tibble())

  data |>
    transmute(value = as.character(.data[[var]])) |>
    filter(!is.na(.data$value), str_trim(.data$value) != "") |>
    separate_longer_delim(.data$value, delim = ";") |>
    mutate(value = str_squish(.data$value)) |>
    filter(.data$value != "") |>
    count(.data$value, name = "n_events", sort = TRUE) |>
    rename(!!out_name := .data$value)
}

make_frequency <- function(data, var) {
  if (!var %in% names(data)) return(tibble())
  data |>
    mutate(value = as.character(.data[[var]]),
           value = if_else(is.na(.data$value) | str_trim(.data$value) == "",
                           "<MISSING>", .data$value)) |>
    count(.data$value, name = "n", sort = TRUE) |>
    mutate(percent = 100 * .data$n / sum(.data$n)) |>
    rename(!!var := .data$value)
}

save_plot <- function(plot, filename, width = 10, height = 6) {
  ggplot2::ggsave(
    filename = file.path(audit_dir_cfg, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

# 3. Locate the most recent ACLED Cabo Delgado typed file --------------------
if (!dir.exists(data_root)) {
  stop("ACLED data directory does not exist: ", data_root)
}

cd_files <- fs::dir_ls(
  data_root,
  recurse = TRUE,
  type = "file",
  regexp = "acled_cabo_delgado_typed\\.rds$"
)

input_source <- "Cabo Delgado typed derivative"

if (length(cd_files) > 0L) {
  info <- fs::file_info(cd_files)
  INPUT_FILE <- as.character(cd_files[which.max(info$modification_time)])
  acled <- readRDS(INPUT_FILE)
} else {
  # Fallback: use latest typed Mozambique master and derive Cabo Delgado.
  moz_files <- fs::dir_ls(
    data_root,
    recurse = TRUE,
    type = "file",
    regexp = "acled_mozambique_typed\\.rds$"
  )

  if (length(moz_files) == 0L) {
    stop(
      "No typed ACLED RDS was found under ", data_root,
      ". Run 01_download_ACLED_Mozambique_2017_2026_v2.R first."
    )
  }

  info <- fs::file_info(moz_files)
  INPUT_FILE <- as.character(moz_files[which.max(info$modification_time)])
  acled_moz <- readRDS(INPUT_FILE)

  if (!"admin1" %in% names(acled_moz)) {
    stop("Fallback Mozambique file has no admin1 column.")
  }

  acled <- acled_moz |>
    filter(.data$admin1 == province_cfg)

  input_source <- "Mozambique typed master; Cabo Delgado derived during audit"
  rm(acled_moz)
}

acled <- tibble::as_tibble(acled)

message("\nacled cabo delgado audit")
message("Project root : ", here::here())
message("Input file   : ", INPUT_FILE)
message("Input source : ", input_source)
message("Rows loaded  : ", nrow(acled))
message("Audit output : ", audit_dir_cfg)
message("\n")

if (nrow(acled) == 0L) stop("Cabo Delgado dataset contains zero rows.")

# 4. Basic typing and derived audit-only fields ------------------------------
if (!"event_date" %in% names(acled)) {
  stop("event_date is required for the audit but is absent.")
}

acled <- acled |>
  mutate(
    event_date = as.Date(.data$event_date),
    audit_year = lubridate::year(.data$event_date),
    year_month = as.Date(lubridate::floor_date(.data$event_date,
      unit = "month")),
    in_main_window = .data$event_date >= main_start & .data$event_date <= main_end
  )

if ("fatalities" %in% names(acled)) {
  acled <- acled |>
    mutate(fatalities = suppressWarnings(as.numeric(.data$fatalities)))
}

for (v in intersect(c("time_precision", "geo_precision"), names(acled))) {
  acled[[v]] <- suppressWarnings(as.integer(acled[[v]]))
}

for (v in intersect(c("latitude", "longitude", "timestamp",
                      "population_1km", "population_2km",
                      "population_5km", "population_best"), names(acled))) {
  acled[[v]] <- suppressWarnings(as.numeric(acled[[v]]))
}

if ("timestamp" %in% names(acled)) {
  acled <- acled |>
    mutate(
      api_timestamp_utc = as.POSIXct(
        .data$timestamp,
        origin = "1970-01-01",
        tz = "UTC"
      ),
      api_upload_lag_days = as.numeric(as.Date(.data$api_timestamp_utc) - .data$event_date)
    )
}

# 5. Schema and missingness audit --------------------------------------------
schema_audit <- tibble(
  variable = names(acled),
  class = vapply(acled, function(x) paste(class(x), collapse = "/"),
    character(1)),
  n_missing_or_blank = vapply(acled, function(x) sum(is_blankish(x)),
    integer(1)),
  pct_missing_or_blank = vapply(acled, function(x) mean(is_blankish(x)) * 100,
    numeric(1)),
  n_unique_nonmissing = vapply(
    acled,
    function(x) dplyr::n_distinct(x[!is_blankish(x)]),
    integer(1)
  )
) |>
  arrange(desc(.data$pct_missing_or_blank), .data$variable)

write_csv_audit(schema_audit, "01_schema_missingness.csv")

# 6. Core integrity checks ---------------------------------------------------
id_field <- if ("event_id_cnty" %in% names(acled)) "event_id_cnty" else NA_character_

n_duplicate_ids <- if (!is.na(id_field)) {
  sum(duplicated(acled[[id_field]]))
} else {
  NA_integer_
}

duplicate_id_rows <- if (!is.na(id_field)) {
  acled |>
    add_count(.data[[id_field]], name = "id_frequency") |>
    filter(.data$id_frequency > 1) |>
    arrange(.data[[id_field]], .data$event_date)
} else {
  tibble()
}
write_csv_audit(duplicate_id_rows, "02_duplicate_event_ids.csv")

# Logical QA flags. These do not filter anything.
qa_flags <- acled |>
  mutate(
    flag_wrong_country = if ("country" %in% names(acled)) {
      is.na(.data$country) | .data$country != country_cfg
    } else FALSE,
    flag_wrong_admin1 = if ("admin1" %in% names(acled)) {
      is.na(.data$admin1) | .data$admin1 != province_cfg
    } else FALSE,
    flag_date_before_2017 = .data$event_date < as.Date("2017-01-01"),
    flag_date_future = .data$event_date > Sys.Date(),
    flag_year_mismatch = if ("year" %in% names(acled)) {
      suppressWarnings(as.integer(.data$year)) != .data$audit_year
    } else FALSE,
    flag_bad_time_precision = if ("time_precision" %in% names(acled)) {
      !is.na(.data$time_precision) & !.data$time_precision %in% 1:3
    } else FALSE,
    flag_bad_geo_precision = if ("geo_precision" %in% names(acled)) {
      !is.na(.data$geo_precision) & !.data$geo_precision %in% 1:3
    } else FALSE,
    flag_missing_coordinates = if (all(c("latitude",
      "longitude") %in% names(acled))) {
      is.na(.data$latitude) | is.na(.data$longitude)
    } else FALSE,
    flag_invalid_coordinate_range = if (all(c("latitude",
      "longitude") %in% names(acled))) {
      (!is.na(.data$latitude) & (.data$latitude < -90 | .data$latitude > 90)) |
        (!is.na(.data$longitude) & (.data$longitude < -180 | .data$longitude > 180))
    } else FALSE,
    flag_negative_fatalities = if ("fatalities" %in% names(acled)) {
      !is.na(.data$fatalities) & .data$fatalities < 0
    } else FALSE,
    flag_noninteger_fatalities = if ("fatalities" %in% names(acled)) {
      !is.na(.data$fatalities) & abs(.data$fatalities - round(.data$fatalities)) > 1e-8
    } else FALSE,
    flag_negative_upload_lag = if ("api_upload_lag_days" %in% names(acled)) {
      !is.na(.data$api_upload_lag_days) & .data$api_upload_lag_days < 0
    } else FALSE
  )

flag_cols <- names(qa_flags)[stringr::str_starts(names(qa_flags), "flag_")]

flag_summary <- tibble(
  flag = flag_cols,
  n_flagged = vapply(qa_flags[flag_cols], function(x) sum(x, na.rm = TRUE),
    integer(1))
) |>
  mutate(percent = 100 * .data$n_flagged / nrow(qa_flags)) |>
  arrange(desc(.data$n_flagged), .data$flag)

write_csv_audit(flag_summary, "03_integrity_flag_summary.csv")

flagged_events <- qa_flags |>
  filter(if_any(all_of(flag_cols), ~ .x %in% TRUE))
write_csv_audit(flagged_events, "04_integrity_flagged_events.csv")

# 7. Global audit summary ----------------------------------------------------
summary_metrics <- list(
  input_file = INPUT_FILE,
  audit_timestamp = as.character(Sys.time()),
  rows = nrow(acled),
  columns = ncol(acled),
  unique_event_ids = if (!is.na(id_field)) n_distinct(acled[[id_field]]) else NA_integer_,
  duplicate_event_ids = n_duplicate_ids,
  min_event_date = min(acled$event_date, na.rm = TRUE),
  max_event_date = max(acled$event_date, na.rm = TRUE),
  rows_main_2017_2024 = sum(acled$in_main_window, na.rm = TRUE),
  observed_admin2_names = if ("admin2" %in% names(acled)) {
    n_distinct(acled$admin2[!is_blankish(acled$admin2)])
  } else NA_integer_,
  observed_admin3_names = if ("admin3" %in% names(acled)) {
    n_distinct(acled$admin3[!is_blankish(acled$admin3)])
  } else NA_integer_,
  zero_fatality_events = if ("fatalities" %in% names(acled)) {
    sum(acled$fatalities == 0, na.rm = TRUE)
  } else NA_integer_,
  positive_fatality_events = if ("fatalities" %in% names(acled)) {
    sum(acled$fatalities > 0, na.rm = TRUE)
  } else NA_integer_,
  total_reported_fatalities = if ("fatalities" %in% names(acled)) {
    sum(acled$fatalities, na.rm = TRUE)
  } else NA_real_,
  max_event_fatalities = if ("fatalities" %in% names(acled)) {
    max(acled$fatalities, na.rm = TRUE)
  } else NA_real_
)

summary_table <- tibble(
  metric = names(summary_metrics),
  value = vapply(summary_metrics, function(x) paste(as.character(x),
    collapse = " | "), character(1))
)
write_csv_audit(summary_table, "00_AUDIT_SUMMARY.csv")

# 8. Temporal coverage and time precision ------------------------------------
monthly_counts <- acled |>
  count(.data$year_month, name = "events") |>
  arrange(.data$year_month)

full_month_grid <- tibble(
  year_month = seq.Date(
    from = floor_date(min(acled$event_date, na.rm = TRUE), "month"),
    to = floor_date(max(acled$event_date, na.rm = TRUE), "month"),
    by = "month"
  )
) |>
  left_join(monthly_counts, by = "year_month") |>
  mutate(events = replace_na(.data$events, 0L))
write_csv_audit(full_month_grid, "05_monthly_event_coverage_full_snapshot.csv")

main_month_grid <- tibble(
  year_month = seq.Date(
    from = floor_date(main_start, "month"),
    to = floor_date(main_end, "month"),
    by = "month"
  )
) |>
  left_join(monthly_counts, by = "year_month") |>
  mutate(
    events = replace_na(.data$events, 0L),
    province_month_has_no_recorded_event = .data$events == 0L
  )
write_csv_audit(main_month_grid, "06_monthly_event_coverage_main_2017_2024.csv")

if ("time_precision" %in% names(acled)) {
  time_precision_overall <- acled |>
    count(.data$time_precision, name = "n", sort = TRUE) |>
    mutate(percent = 100 * .data$n / sum(.data$n))
  write_csv_audit(time_precision_overall, "07_time_precision_overall.csv")

  time_precision_year <- acled |>
    count(.data$audit_year, .data$time_precision, name = "n") |>
    group_by(.data$audit_year) |>
    mutate(percent_within_year = 100 * .data$n / sum(.data$n)) |>
    ungroup()
  write_csv_audit(time_precision_year, "08_time_precision_by_year.csv")
}

if ("timestamp" %in% names(acled)) {
  upload_lag_summary <- acled |>
    summarise(
      n = n(),
      timestamp_min_utc = min(.data$api_timestamp_utc, na.rm = TRUE),
      timestamp_max_utc = max(.data$api_timestamp_utc, na.rm = TRUE),
      lag_min_days = min(.data$api_upload_lag_days, na.rm = TRUE),
      lag_p25_days = quantile(.data$api_upload_lag_days, 0.25, na.rm = TRUE),
      lag_median_days = median(.data$api_upload_lag_days, na.rm = TRUE),
      lag_p75_days = quantile(.data$api_upload_lag_days, 0.75, na.rm = TRUE),
      lag_p90_days = quantile(.data$api_upload_lag_days, 0.90, na.rm = TRUE),
      lag_p95_days = quantile(.data$api_upload_lag_days, 0.95, na.rm = TRUE),
      lag_p99_days = quantile(.data$api_upload_lag_days, 0.99, na.rm = TRUE),
      lag_max_days = max(.data$api_upload_lag_days, na.rm = TRUE)
    )
  write_csv_audit(upload_lag_summary, "09_api_upload_timestamp_lag_summary.csv")

  upload_lag_year <- acled |>
    group_by(.data$audit_year) |>
    summarise(
      n = n(),
      median_upload_lag_days = median(.data$api_upload_lag_days, na.rm = TRUE),
      p90_upload_lag_days = quantile(.data$api_upload_lag_days, 0.90,
        na.rm = TRUE),
      max_upload_lag_days = max(.data$api_upload_lag_days, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_audit(upload_lag_year,
    "10_api_upload_timestamp_lag_by_event_year.csv")
}

# 9. Geography, district-name stability, geo precision -----------------------
for (v in intersect(c("admin1", "admin2", "admin3", "location"),
  names(acled))) {
  write_csv_audit(make_frequency(acled, v), paste0("11_frequency_", v, ".csv"))
}

if ("admin2" %in% names(acled)) {
  district_summary <- acled |>
    filter(!is_blankish(.data$admin2)) |>
    group_by(.data$admin2) |>
    summarise(
      events = n(),
      first_event = min(.data$event_date, na.rm = TRUE),
      last_event = max(.data$event_date, na.rm = TRUE),
      years_with_events = n_distinct(.data$audit_year),
      months_with_events = n_distinct(.data$year_month),
      zero_fatality_events = if ("fatalities" %in% names(acled)) sum(.data$fatalities == 0, na.rm = TRUE) else NA_integer_,
      total_fatalities = if ("fatalities" %in% names(acled)) sum(.data$fatalities, na.rm = TRUE) else NA_real_,
      pct_geo_precision_1 = if ("geo_precision" %in% names(acled)) mean(.data$geo_precision == 1, na.rm = TRUE) * 100 else NA_real_,
      pct_geo_precision_3 = if ("geo_precision" %in% names(acled)) mean(.data$geo_precision == 3, na.rm = TRUE) * 100 else NA_real_,
      pct_time_precision_1 = if ("time_precision" %in% names(acled)) mean(.data$time_precision == 1, na.rm = TRUE) * 100 else NA_real_,
      .groups = "drop"
    ) |>
    arrange(.data$admin2)
  write_csv_audit(district_summary, "12_admin2_district_summary.csv")

  district_year <- acled |>
    filter(!is_blankish(.data$admin2)) |>
    count(.data$admin2, .data$audit_year, name = "events") |>
    arrange(.data$admin2, .data$audit_year)
  write_csv_audit(district_year, "13_admin2_by_year.csv")

  # IMPORTANT: diagnostic completeness using names OBSERVED in ACLED only.
  # This is not the canonical study-domain panel.
  observed_districts <- sort(unique(as.character(acled$admin2[!is_blankish(acled$admin2)])))

  diagnostic_grid <- tidyr::expand_grid(
    admin2 = observed_districts,
    year_month = seq.Date(
      floor_date(main_start, "month"),
      floor_date(main_end, "month"),
      by = "month"
    )
  )

  observed_district_month <- acled |>
    filter(.data$in_main_window, !is_blankish(.data$admin2)) |>
    count(.data$admin2, .data$year_month, name = "events")

  diagnostic_grid <- diagnostic_grid |>
    left_join(observed_district_month, by = c("admin2", "year_month")) |>
    mutate(
      events = replace_na(.data$events, 0L),
      no_recorded_event = .data$events == 0L,
      district_set_basis = "ACLED-observed admin2 names only; NOT canonical study domain"
    )
  write_csv_audit(diagnostic_grid,
    "14_DIAGNOSTIC_observed_admin2_month_grid_2017_2024.csv")

  district_month_coverage <- diagnostic_grid |>
    group_by(.data$admin2) |>
    summarise(
      expected_months_2017_2024 = n(),
      months_with_recorded_event = sum(.data$events > 0),
      months_without_recorded_event = sum(.data$events == 0),
      pct_months_with_recorded_event = mean(.data$events > 0) * 100,
      .groups = "drop"
    ) |>
    arrange(.data$admin2)
  write_csv_audit(district_month_coverage,
    "15_DIAGNOSTIC_observed_admin2_month_coverage.csv")

  # Normalization / possible spelling variants.
  name_variants <- tibble(admin2 = observed_districts) |>
    mutate(admin2_normalized = normalize_geo_name(.data$admin2)) |>
    group_by(.data$admin2_normalized) |>
    summarise(
      n_original_spellings = n_distinct(.data$admin2),
      original_spellings = paste(sort(unique(.data$admin2)), collapse = " | "),
      .groups = "drop"
    ) |>
    arrange(desc(.data$n_original_spellings), .data$admin2_normalized)
  write_csv_audit(name_variants, "16_admin2_normalized_name_variants.csv")

  # Conservative fuzzy-name candidates for manual review only.
  if (length(observed_districts) >= 2L) {
    pairs <- t(combn(observed_districts, 2))
    fuzzy_pairs <- tibble(
      admin2_a = pairs[, 1],
      admin2_b = pairs[, 2],
      normalized_a = normalize_geo_name(pairs[, 1]),
      normalized_b = normalize_geo_name(pairs[, 2])
    ) |>
      mutate(edit_distance = mapply(
        function(a, b) as.integer(utils::adist(a, b)),
        .data$normalized_a,
        .data$normalized_b
      )) |>
      filter(.data$edit_distance <= 3L) |>
      arrange(.data$edit_distance, .data$admin2_a, .data$admin2_b)
    write_csv_audit(fuzzy_pairs,
      "17_admin2_fuzzy_name_candidates_MANUAL_REVIEW.csv")
  }
}

if ("geo_precision" %in% names(acled)) {
  geo_precision_overall <- acled |>
    count(.data$geo_precision, name = "n", sort = TRUE) |>
    mutate(percent = 100 * .data$n / sum(.data$n))
  write_csv_audit(geo_precision_overall, "18_geo_precision_overall.csv")

  if ("admin2" %in% names(acled)) {
    geo_precision_district <- acled |>
      count(.data$admin2, .data$geo_precision, name = "n") |>
      group_by(.data$admin2) |>
      mutate(percent_within_admin2 = 100 * .data$n / sum(.data$n)) |>
      ungroup()
    write_csv_audit(geo_precision_district, "19_geo_precision_by_admin2.csv")
  }

  geo_precision_year <- acled |>
    count(.data$audit_year, .data$geo_precision, name = "n") |>
    group_by(.data$audit_year) |>
    mutate(percent_within_year = 100 * .data$n / sum(.data$n)) |>
    ungroup()
  write_csv_audit(geo_precision_year, "20_geo_precision_by_year.csv")
}

# Coordinate diagnostics. Repeated coordinates are expected in ACLED because
# coordinates identify named places/centroids; they are not treated as
# duplicates.
if (all(c("latitude", "longitude") %in% names(acled))) {
  coordinate_summary <- acled |>
    summarise(
      n_events = n(),
      missing_coordinates = sum(is.na(.data$latitude) | is.na(.data$longitude)),
      unique_coordinate_pairs = n_distinct(paste(.data$latitude,
        .data$longitude)),
      latitude_min = min(.data$latitude, na.rm = TRUE),
      latitude_max = max(.data$latitude, na.rm = TRUE),
      longitude_min = min(.data$longitude, na.rm = TRUE),
      longitude_max = max(.data$longitude, na.rm = TRUE)
    )
  write_csv_audit(coordinate_summary, "21_coordinate_summary.csv")

  repeated_coordinates <- acled |>
    filter(!is.na(.data$latitude), !is.na(.data$longitude)) |>
    count(.data$latitude, .data$longitude, name = "n_events", sort = TRUE) |>
    filter(.data$n_events > 1)
  write_csv_audit(repeated_coordinates,
    "22_repeated_coordinate_pairs_NOT_ERROR.csv")

  if (all(c("location", "admin2") %in% names(acled))) {
    location_coord <- acled |>
      filter(!is_blankish(.data$location), !is.na(.data$latitude),
        !is.na(.data$longitude)) |>
      group_by(.data$admin2, .data$location) |>
      summarise(
        events = n(),
        unique_coordinate_pairs = n_distinct(paste(.data$latitude,
          .data$longitude)),
        .groups = "drop"
      ) |>
      arrange(desc(.data$unique_coordinate_pairs), desc(.data$events))
    write_csv_audit(location_coord, "23_location_coordinate_consistency.csv")
  }
}

# 10. Optional comparison with an existing district boundary object ----------
boundary_file <- boundary_candidates[file.exists(boundary_candidates)][1] %||% NA_character_
boundary_used <- FALSE
boundary_name_field <- NA_character_

if (!is.na(boundary_file) && nzchar(boundary_file)) {
  boundary_obj <- tryCatch(readRDS(boundary_file), error = function(e) NULL)

  if (!is.null(boundary_obj) && inherits(boundary_obj, "sf")) {
    candidate_name_fields <- c(
      "NAME_2", "admin2", "ADM2_EN", "ADM2_NAME", "district", "DISTRICT"
    )
    boundary_name_field <- intersect(candidate_name_fields,
      names(boundary_obj))[1] %||% NA_character_

    if (!is.na(boundary_name_field)) {
      boundary_used <- TRUE

      boundary_names <- boundary_obj |>
        st_drop_geometry() |>
        transmute(
          boundary_admin2 = as.character(.data[[boundary_name_field]]),
          normalized = normalize_geo_name(.data$boundary_admin2)
        ) |>
        filter(!is_blankish(.data$boundary_admin2)) |>
        distinct()

      acled_names <- if ("admin2" %in% names(acled)) {
        tibble(acled_admin2 = sort(unique(as.character(acled$admin2[!is_blankish(acled$admin2)])))) |>
          mutate(normalized = normalize_geo_name(.data$acled_admin2))
      } else {
        tibble(acled_admin2 = character(), normalized = character())
      }

      boundary_comparison <- full_join(
        boundary_names,
        acled_names,
        by = "normalized"
      ) |>
        mutate(
          exact_original_match = !is.na(.data$boundary_admin2) &
            !is.na(.data$acled_admin2) &
            .data$boundary_admin2 == .data$acled_admin2,
          normalized_match = !is.na(.data$boundary_admin2) & !is.na(.data$acled_admin2),
          status = case_when(
            is.na(.data$acled_admin2) ~ "BOUNDARY_ONLY",
            is.na(.data$boundary_admin2) ~ "ACLED_ONLY",
            .data$exact_original_match ~ "EXACT_MATCH",
            TRUE ~ "NORMALIZED_MATCH"
          )
        ) |>
        arrange(.data$status, .data$boundary_admin2, .data$acled_admin2)

      write_csv_audit(boundary_comparison,
        "24_boundary_vs_ACLED_admin2_names.csv")

      boundary_summary <- tibble(
        boundary_file = boundary_file,
        boundary_name_field = boundary_name_field,
        boundary_admin2_count = n_distinct(boundary_names$boundary_admin2),
        acled_observed_admin2_count = n_distinct(acled_names$acled_admin2),
        exact_or_normalized_matches = sum(boundary_comparison$normalized_match,
          na.rm = TRUE),
        boundary_only = sum(boundary_comparison$status == "BOUNDARY_ONLY",
          na.rm = TRUE),
        acled_only = sum(boundary_comparison$status == "ACLED_ONLY",
          na.rm = TRUE)
      )
      write_csv_audit(boundary_summary, "25_boundary_comparison_summary.csv")

      # Spatial containment check: diagnostic only, using the existing boundary.
      valid_pts <- acled |>
        filter(!is.na(.data$latitude), !is.na(.data$longitude))

      if (nrow(valid_pts) > 0L) {
        pts_sf <- st_as_sf(
          valid_pts,
          coords = c("longitude", "latitude"),
          crs = 4326,
          remove = FALSE
        )
        boundary_wgs84 <- st_transform(boundary_obj, 4326)
        province_union <- st_union(st_geometry(boundary_wgs84))
        inside <- lengths(st_within(pts_sf, province_union)) > 0

        containment_summary <- tibble(
          boundary_file = boundary_file,
          events_with_coordinates = nrow(pts_sf),
          inside_boundary = sum(inside),
          outside_boundary = sum(!inside),
          pct_inside = mean(inside) * 100
        )
        write_csv_audit(containment_summary,
          "26_coordinate_containment_existing_boundary.csv")

        if (any(!inside)) {
          outside_events <- valid_pts[!inside, , drop = FALSE]
          write_csv_audit(outside_events,
            "27_events_outside_existing_boundary_MANUAL_REVIEW.csv")
        }
      }
    }
  }
}

# 11. Event taxonomy: disorder -> event -> sub-event -------------------------
for (v in intersect(c("disorder_type", "event_type", "sub_event_type",
  "civilian_targeting"), names(acled))) {
  write_csv_audit(make_frequency(acled, v), paste0("28_taxonomy_", v, ".csv"))
}

if (all(c("disorder_type", "event_type", "sub_event_type") %in% names(acled))) {
  taxonomy_crosswalk <- acled |>
    count(.data$disorder_type, .data$event_type, .data$sub_event_type,
      name = "n_events") |>
    arrange(.data$disorder_type, .data$event_type, desc(.data$n_events))
  write_csv_audit(taxonomy_crosswalk, "29_taxonomy_crosswalk_observed.csv")
}

if ("event_type" %in% names(acled)) {
  event_type_year <- acled |>
    count(.data$audit_year, .data$event_type, name = "events") |>
    arrange(.data$audit_year, desc(.data$events))
  write_csv_audit(event_type_year, "30_event_type_by_year.csv")
}

# 12. Actors and interactions ------------------------------------------------
for (v in intersect(c("actor1", "actor2", "inter1", "inter2", "interaction"),
  names(acled))) {
  write_csv_audit(make_frequency(acled, v), paste0("31_actor_field_", v,
    ".csv"))
}

if (any(c("actor1", "actor2") %in% names(acled))) {
  actor_long <- bind_rows(
    if ("actor1" %in% names(acled)) acled |> transmute(actor = as.character(.data$actor1), role = "actor1") else tibble(),
    if ("actor2" %in% names(acled)) acled |> transmute(actor = as.character(.data$actor2), role = "actor2") else tibble()
  ) |>
    filter(!is.na(.data$actor), str_trim(.data$actor) != "")

  actor_any_role <- actor_long |>
    count(.data$actor, name = "actor_appearances", sort = TRUE)
  write_csv_audit(actor_any_role, "32_actors_any_main_role.csv")

  actor_role <- actor_long |>
    count(.data$actor, .data$role, name = "appearances") |>
    arrange(.data$actor, .data$role)
  write_csv_audit(actor_role, "33_actors_by_main_role.csv")
}

for (v in intersect(c("assoc_actor_1", "assoc_actor_2"), names(acled))) {
  out <- explode_semicolon_counts(acled, v, out_name = v)
  write_csv_audit(out, paste0("34_", v, "_exploded.csv"))
}

# 13. Sources and reporting environment --------------------------------------
if ("source_scale" %in% names(acled)) {
  write_csv_audit(make_frequency(acled, "source_scale"),
    "35_source_scale_exact_strings.csv")

  source_scale_year <- acled |>
    count(.data$audit_year, .data$source_scale, name = "n") |>
    group_by(.data$audit_year) |>
    mutate(percent_within_year = 100 * .data$n / sum(.data$n)) |>
    ungroup()
  write_csv_audit(source_scale_year, "36_source_scale_by_year.csv")
}

if ("source" %in% names(acled)) {
  source_counts <- explode_semicolon_counts(acled, "source",
    out_name = "source")
  write_csv_audit(source_counts, "37_sources_exploded.csv")

  source_per_event <- acled |>
    mutate(n_source_labels = count_semicolon_items(.data$source)) |>
    count(.data$n_source_labels, name = "events") |>
    mutate(percent = 100 * .data$events / sum(.data$events)) |>
    arrange(.data$n_source_labels)
  write_csv_audit(source_per_event, "38_number_of_source_labels_per_event.csv")

  monthly_reporting <- acled |>
    mutate(n_source_labels = count_semicolon_items(.data$source)) |>
    group_by(.data$year_month) |>
    summarise(
      events = n(),
      mean_source_labels_per_event = mean(.data$n_source_labels, na.rm = TRUE),
      median_source_labels_per_event = median(.data$n_source_labels,
        na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_audit(monthly_reporting,
    "39_monthly_reporting_source_diagnostics.csv")
}

# 14. Fatalities: zeros, tail, extremes, and relationships to taxonomy -------
if ("fatalities" %in% names(acled)) {
  fatality_summary <- acled |>
    summarise(
      n_events = n(),
      missing_fatalities = sum(is.na(.data$fatalities)),
      zero_fatality_events = sum(.data$fatalities == 0, na.rm = TRUE),
      positive_fatality_events = sum(.data$fatalities > 0, na.rm = TRUE),
      zero_percent = mean(.data$fatalities == 0, na.rm = TRUE) * 100,
      total_fatalities = sum(.data$fatalities, na.rm = TRUE),
      mean = mean(.data$fatalities, na.rm = TRUE),
      variance = var(.data$fatalities, na.rm = TRUE),
      variance_to_mean = var(.data$fatalities,
        na.rm = TRUE) / mean(.data$fatalities, na.rm = TRUE),
      min = min(.data$fatalities, na.rm = TRUE),
      p25 = quantile(.data$fatalities, 0.25, na.rm = TRUE),
      median = median(.data$fatalities, na.rm = TRUE),
      p75 = quantile(.data$fatalities, 0.75, na.rm = TRUE),
      p90 = quantile(.data$fatalities, 0.90, na.rm = TRUE),
      p95 = quantile(.data$fatalities, 0.95, na.rm = TRUE),
      p99 = quantile(.data$fatalities, 0.99, na.rm = TRUE),
      max = max(.data$fatalities, na.rm = TRUE)
    )
  write_csv_audit(fatality_summary, "40_fatality_distribution_summary.csv")

  fatality_frequency <- acled |>
    count(.data$fatalities, name = "events", sort = FALSE) |>
    arrange(.data$fatalities)
  write_csv_audit(fatality_frequency, "41_fatality_exact_frequency.csv")

  fatality_bins <- acled |>
    mutate(
      fatality_bin = case_when(
        is.na(.data$fatalities) ~ "Missing",
        .data$fatalities == 0 ~ "0",
        .data$fatalities == 1 ~ "1",
        .data$fatalities <= 4 ~ "2-4",
        .data$fatalities <= 10 ~ "5-10",
        .data$fatalities <= 20 ~ "11-20",
        .data$fatalities <= 50 ~ "21-50",
        TRUE ~ ">50"
      )
    ) |>
    count(.data$fatality_bin, name = "events") |>
    mutate(percent = 100 * .data$events / sum(.data$events))
  write_csv_audit(fatality_bins, "42_fatality_descriptive_bins.csv")

  if ("event_type" %in% names(acled)) {
    fatalities_event_type <- acled |>
      group_by(.data$event_type) |>
      summarise(
        events = n(),
        zero_events = sum(.data$fatalities == 0, na.rm = TRUE),
        zero_percent = mean(.data$fatalities == 0, na.rm = TRUE) * 100,
        positive_events = sum(.data$fatalities > 0, na.rm = TRUE),
        total_fatalities = sum(.data$fatalities, na.rm = TRUE),
        mean_fatalities = mean(.data$fatalities, na.rm = TRUE),
        max_fatalities = max(.data$fatalities, na.rm = TRUE),
        .groups = "drop"
      ) |>
      arrange(desc(.data$events))
    write_csv_audit(fatalities_event_type, "43_fatalities_by_event_type.csv")
  }

  if ("sub_event_type" %in% names(acled)) {
    fatalities_sub_event_type <- acled |>
      group_by(.data$sub_event_type) |>
      summarise(
        events = n(),
        zero_percent = mean(.data$fatalities == 0, na.rm = TRUE) * 100,
        total_fatalities = sum(.data$fatalities, na.rm = TRUE),
        mean_fatalities = mean(.data$fatalities, na.rm = TRUE),
        max_fatalities = max(.data$fatalities, na.rm = TRUE),
        .groups = "drop"
      ) |>
      arrange(desc(.data$events))
    write_csv_audit(fatalities_sub_event_type,
      "44_fatalities_by_sub_event_type.csv")
  }

  extreme_fields <- intersect(
    c(
      "event_id_cnty", "event_date", "disorder_type", "event_type", "sub_event_type",
      "actor1", "actor2", "admin2", "admin3", "location", "latitude", "longitude",
      "geo_precision", "time_precision", "source", "source_scale", "fatalities", "notes"
    ),
    names(acled)
  )

  top_extremes <- acled |>
    arrange(desc(.data$fatalities), .data$event_date) |>
    select(all_of(extreme_fields)) |>
    slice_head(n = 100)
  write_csv_audit(top_extremes, "45_top_100_fatality_events_MANUAL_REVIEW.csv")
}

# 15. Population-exposure fields returned by ACLED (event-level diagnostic only)
pop_fields <- intersect(
  c("population_1km", "population_2km", "population_5km", "population_best"),
  names(acled)
)

if (length(pop_fields) > 0L) {
  population_audit <- map_dfr(pop_fields, function(v) {
    x <- acled[[v]]
    tibble(
      variable = v,
      n = length(x),
      missing = sum(is.na(x)),
      pct_missing = mean(is.na(x)) * 100,
      min = suppressWarnings(min(x, na.rm = TRUE)),
      median = suppressWarnings(median(x, na.rm = TRUE)),
      mean = suppressWarnings(mean(x, na.rm = TRUE)),
      max = suppressWarnings(max(x, na.rm = TRUE))
    )
  })
  write_csv_audit(population_audit,
    "46_ACLED_event_level_population_fields.csv")
}

# 16. Snapshot fingerprint for future revision detection ---------------------
if (!is.na(id_field)) {
  hash_fields <- intersect(
    c(
      "event_date", "time_precision", "disorder_type", "event_type", "sub_event_type",
      "actor1", "assoc_actor_1", "inter1", "actor2", "assoc_actor_2", "inter2",
      "interaction", "civilian_targeting", "admin1", "admin2", "admin3", "location",
      "latitude", "longitude", "geo_precision", "source", "source_scale", "notes",
      "fatalities", "tags"
    ),
    names(acled)
  )

  message(
    "Creating event-level snapshot fingerprints for future revision ",
    "checks...")
  record_hash <- vapply(
    seq_len(nrow(acled)),
    function(i) {
      digest::digest(
        as.list(acled[i, hash_fields, drop = FALSE]),
        algo = "xxhash64",
        serialize = TRUE
      )
    },
    character(1)
  )

  fingerprint <- tibble(
    event_id_cnty = as.character(acled[[id_field]]),
    event_date = acled$event_date,
    api_timestamp_utc = if ("api_timestamp_utc" %in% names(acled)) acled$api_timestamp_utc else as.POSIXct(NA),
    record_hash = record_hash
  ) |>
    arrange(.data$event_id_cnty)

  write_csv_audit(fingerprint,
    "47_snapshot_fingerprint_for_future_revision_comparison.csv")
}

# 17. Figures for visual audit -----------------------------------------------
p_monthly <- full_month_grid |>
  ggplot(aes(x = .data$year_month, y = .data$events)) +
  geom_line(linewidth = 0.5) +
  labs(
    title = "ACLED events in Cabo Delgado by month",
    x = NULL,
    y = "Recorded events"
  ) +
  theme_minimal(base_size = 12)
save_plot(p_monthly, "Figure_01_monthly_events.png", width = 11, height = 6)

if ("event_type" %in% names(acled)) {
  event_type_plot_df <- acled |>
    filter(.data$in_main_window) |>
    count(.data$audit_year, .data$event_type, name = "events")

  p_types <- event_type_plot_df |>
    ggplot(aes(x = factor(.data$audit_year), y = .data$events,
      fill = .data$event_type)) +
    geom_col() +
    labs(
      title = "Event-type composition by year, 2017-2024",
      x = "Year",
      y = "Recorded events",
      fill = "Event type"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_types, "Figure_02_event_type_by_year_2017_2024.png", width = 12,
    height = 7)
}

if ("fatalities" %in% names(acled)) {
  p_fatal <- acled |>
    ggplot(aes(x = log1p(.data$fatalities))) +
    geom_histogram(bins = 50) +
    labs(
      title = "Distribution of reported fatalities per ACLED event",
      subtitle = "x-axis uses log(1 + fatalities)",
      x = "log(1 + fatalities)",
      y = "Events"
    ) +
    theme_minimal(base_size = 12)
  save_plot(p_fatal, "Figure_03_log1p_fatalities.png", width = 9, height = 6)
}

if ("admin2" %in% names(acled)) {
  district_plot_df <- acled |>
    filter(.data$in_main_window, !is_blankish(.data$admin2)) |>
    count(.data$admin2, name = "events") |>
    arrange(.data$events) |>
    mutate(admin2 = factor(.data$admin2, levels = .data$admin2))

  p_district <- district_plot_df |>
    ggplot(aes(x = .data$admin2, y = .data$events)) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Recorded ACLED events by observed Admin 2, 2017-2024",
      x = NULL,
      y = "Recorded events"
    ) +
    theme_minimal(base_size = 11)
  save_plot(p_district, "Figure_04_events_by_observed_admin2_2017_2024.png",
    width = 9, height = 8)
}

# Save a compact R object containing audit metadata, not a modified analytical
# dataset.
audit_metadata <- list(
  created_at = Sys.time(),
  project_root = here::here(),
  input_file = INPUT_FILE,
  input_source = input_source,
  audit_dir = audit_dir_cfg,
  country = country_cfg,
  province = province_cfg,
  main_window = c(as.character(main_start), as.character(main_end)),
  rows = nrow(acled),
  observed_date_range = c(
    as.character(min(acled$event_date, na.rm = TRUE)),
    as.character(max(acled$event_date, na.rm = TRUE))
  ),
  duplicate_event_ids = n_duplicate_ids,
  boundary_used = boundary_used,
  boundary_file = if (boundary_used) boundary_file else NA_character_,
  boundary_name_field = boundary_name_field
)
saveRDS(audit_metadata, file.path(audit_dir_cfg, "audit_metadata.rds"))

capture.output(
  sessionInfo(),
  file = file.path(audit_dir_cfg, "R_sessionInfo.txt")
)

# 20. Console report ---------------------------------------------------------
message("\naudit complete")
message("Events audited            : ", nrow(acled))
message("Observed event IDs        : ",
  if (!is.na(id_field)) n_distinct(acled[[id_field]]) else "event_id_cnty absent")
message("Duplicate event IDs       : ", n_duplicate_ids)
message("Observed date range       : ", min(acled$event_date, na.rm = TRUE),
  " to ", max(acled$event_date, na.rm = TRUE))
message("Rows in main 2017-2024    : ", sum(acled$in_main_window, na.rm = TRUE))
message("Observed Admin 2 names    : ",
  if ("admin2" %in% names(acled)) n_distinct(acled$admin2[!is_blankish(acled$admin2)]) else NA)
message("Boundary comparison       : ", ifelse(boundary_used, "performed",
  "not performed"))
message("Integrity flags with n>0  : ", sum(flag_summary$n_flagged > 0))
message("Audit output directory    : ", audit_dir_cfg)
message("\n")

print(summary_table, n = Inf)
print(flag_summary, n = Inf)

if (exists("district_summary")) print(district_summary, n = Inf)
if (exists("time_precision_overall")) print(time_precision_overall, n = Inf)
if (exists("geo_precision_overall")) print(geo_precision_overall, n = Inf)
if (exists("fatality_summary")) print(fatality_summary, n = Inf)

