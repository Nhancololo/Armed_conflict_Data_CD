# Review ACLED events whose coordinates fall outside their district.

# 0. Packages ----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  dplyr,
  tidyr,
  readr,
  tibble,
  stringr,
  sf,
  ggplot2,
  here,
  fs
)

# 1. Configuration -----------------------------------------------------------
analysis_start <- as.Date("2017-01-01")
analysis_end   <- as.Date("2024-12-31")

# Cabo Delgado is covered by UTM zone 37S.
distance_crs <- 32737

processed_dir <- here::here("data", "processed")

event_file <- file.path(
  processed_dir,
  "acled_cabo_delgado_harmonized.rds"
)

boundary_file <- file.path(
  processed_dir,
  "cabo_delgado_admin2_canonical.rds"
)

output_dir <- here::here(
  "outputs",
  "acled_spatial_review"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Basic helpers -----------------------------------------------------------
write_qa <- function(x, filename) {
  readr::write_csv(
    x,
    file.path(output_dir, filename),
    na = ""
  )
  invisible(x)
}

fmt_pct <- function(x) {
  round(100 * x, 3)
}

# 3. Load inputs -------------------------------------------------------------
if (!file.exists(event_file)) {
  stop(
    "Harmonized ACLED file not found:\n",
    event_file,
    "\nRun script 03 first."
  )
}

if (!file.exists(boundary_file)) {
  stop(
    "Canonical boundary file not found:\n",
    boundary_file,
    "\nRun script 03 first."
  )
}

events <- readRDS(event_file) |>
  tibble::as_tibble()

districts <- readRDS(boundary_file)

if (!inherits(districts, "sf")) {
  stop("Canonical boundary object is not an sf object.")
}

required_event_fields <- c(
  "event_id_cnty",
  "event_date",
  "district_id",
  "canonical_name",
  "longitude",
  "latitude",
  "geo_precision",
  "spatial_status"
)

missing_event_fields <- setdiff(required_event_fields, names(events))

if (length(missing_event_fields) > 0L) {
  stop(
    "Harmonized ACLED file is missing required fields: ",
    paste(missing_event_fields, collapse = ", ")
  )
}

required_boundary_fields <- c(
  "district_id",
  "canonical_name"
)

missing_boundary_fields <- setdiff(
  required_boundary_fields,
  names(districts)
)

if (length(missing_boundary_fields) > 0L) {
  stop(
    "Canonical boundary is missing required fields: ",
    paste(missing_boundary_fields, collapse = ", ")
  )
}

events <- events |>
  mutate(
    event_date = as.Date(.data$event_date),
    longitude = as.numeric(.data$longitude),
    latitude = as.numeric(.data$latitude),
    geo_precision = as.integer(.data$geo_precision),
    in_primary_window =
      .data$event_date >= analysis_start &
      .data$event_date <= analysis_end
  )

if (anyDuplicated(events$event_id_cnty) > 0L) {
  stop("Duplicate event_id_cnty detected in harmonized ACLED file.")
}

if (any(is.na(events$district_id))) {
  stop(
    "At least one event has missing district_id. ",
    "Resolve script 03 harmonization before running this review."
  )
}

if (is.na(sf::st_crs(districts))) {
  stop("Canonical district geometry has no CRS.")
}

# 4. Identify events requiring spatial review --------------------------------
review <- events |>
  filter(
    .data$spatial_status != "text_spatial_agree"
  )

if (nrow(review) == 0L) {
  message("No spatial mismatches/outside events found. Nothing to review.")
  quit(save = "no", status = 0)
}

if (any(is.na(review$longitude) | is.na(review$latitude))) {
  missing_xy <- review |>
    filter(is.na(.data$longitude) | is.na(.data$latitude))

  write_qa(
    missing_xy,
    "ERROR_missing_coordinates_in_review_cases.csv"
  )

  stop(
    nrow(missing_xy),
    " review events have missing longitude/latitude. ",
    "Inspect ERROR_missing_coordinates_in_review_cases.csv."
  )
}

# 5. Convert review events and polygons to projected CRS for distances -------
review_sf <- review |>
  sf::st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

review_utm <- sf::st_transform(
  review_sf,
  distance_crs
)

districts_utm <- districts |>
  sf::st_make_valid() |>
  sf::st_transform(distance_crs)

# Ensure exactly one polygon per district_id.
district_id_counts <- districts_utm |>
  sf::st_drop_geometry() |>
  count(.data$district_id, name = "n_features")

if (any(district_id_counts$n_features != 1L)) {
  write_qa(
    district_id_counts,
    "ERROR_boundary_feature_count_by_district.csv"
  )

  stop(
    "Canonical boundary must contain exactly one feature per district_id."
  )
}

# 6. Distance to ACLED-assigned canonical district polygon -------------------
assigned_index <- match(
  review_utm$district_id,
  districts_utm$district_id
)

if (any(is.na(assigned_index))) {
  bad <- unique(review_utm$district_id[is.na(assigned_index)])

  stop(
    "Some event district IDs do not exist in canonical polygons: ",
    paste(bad, collapse = ", ")
  )
}

assigned_geom <- sf::st_geometry(
  districts_utm[assigned_index, ]
)

distance_to_assigned_m <- as.numeric(
  sf::st_distance(
    sf::st_geometry(review_utm),
    assigned_geom,
    by_element = TRUE
  )
)

# 7. Nearest canonical district polygon --------------------------------------
nearest_index <- sf::st_nearest_feature(
  review_utm,
  districts_utm
)

nearest_geom <- sf::st_geometry(
  districts_utm[nearest_index, ]
)

distance_to_nearest_m <- as.numeric(
  sf::st_distance(
    sf::st_geometry(review_utm),
    nearest_geom,
    by_element = TRUE
  )
)

nearest_district_id <- districts_utm$district_id[nearest_index]
nearest_district_name <- districts_utm$canonical_name[nearest_index]

# 8. Build diagnostic review table -------------------------------------------
review_table <- review |>
  mutate(
    distance_to_assigned_polygon_m =
      round(distance_to_assigned_m, 1),

    nearest_district_id =
      nearest_district_id,

    nearest_district_name =
      nearest_district_name,

    distance_to_nearest_polygon_m =
      round(distance_to_nearest_m, 1),

    nearest_is_assigned =
      .data$nearest_district_id == .data$district_id,

    review_priority = case_when(
      .data$geo_precision == 1L ~ "HIGH",
      .data$geo_precision == 2L ~ "MEDIUM",
      .data$geo_precision == 3L ~ "LOW",
      TRUE ~ "UNKNOWN"
    ),

    # Purely descriptive bins. These do not imply reassignment.
    distance_band = case_when(
      .data$distance_to_nearest_polygon_m == 0 ~
        "inside_another_polygon_or_on_boundary",

      .data$distance_to_nearest_polygon_m <= 100 ~
        "outside_<=100m",

      .data$distance_to_nearest_polygon_m <= 500 ~
        "outside_100-500m",

      .data$distance_to_nearest_polygon_m <= 1000 ~
        "outside_500m-1km",

      .data$distance_to_nearest_polygon_m <= 5000 ~
        "outside_1-5km",

      .data$distance_to_nearest_polygon_m <= 10000 ~
        "outside_5-10km",

      TRUE ~
        "outside_>10km"
    ),

    diagnostic_interpretation = case_when(
      .data$spatial_status == "text_spatial_disagree" ~
        "point falls inside a different canonical district polygon",

      .data$spatial_status == "outside_all_polygons" ~
        "point falls outside every canonical district polygon",

      .data$spatial_status == "boundary_ambiguous" ~
        "point intersects more than one polygon/boundary",

      .data$spatial_status == "text_unresolved" ~
        "textual district assignment unresolved",

      TRUE ~
        .data$spatial_status
    )
  ) |>
  arrange(
    factor(
      .data$review_priority,
      levels = c("HIGH", "MEDIUM", "LOW", "UNKNOWN")
    ),
    desc(.data$distance_to_assigned_polygon_m),
    .data$event_date
  )

# Put the most useful review columns first.
preferred_columns <- c(
  "review_priority",
  "event_id_cnty",
  "event_date",
  "in_primary_window",
  "admin2_original",
  "district_id",
  "canonical_name",
  "geo_precision",
  "location",
  "latitude",
  "longitude",
  "spatial_status",
  "spatial_district_id",
  "spatial_canonical_name",
  "nearest_district_id",
  "nearest_district_name",
  "nearest_is_assigned",
  "distance_to_assigned_polygon_m",
  "distance_to_nearest_polygon_m",
  "distance_band",
  "diagnostic_interpretation",
  "event_type",
  "sub_event_type",
  "fatalities",
  "source",
  "source_scale",
  "notes"
)

review_table <- review_table |>
  select(
    any_of(preferred_columns),
    everything()
  )

# 9. Summary tables ----------------------------------------------------------
n_events <- nrow(events)
n_review <- nrow(review_table)
n_agree <- sum(events$spatial_status == "text_spatial_agree", na.rm = TRUE)

summary_table <- tibble(
  metric = c(
    "total_ACLED_events",
    "text_spatial_agreements",
    "non_agree_events_total",
    "pct_non_agree_events",
    "text_spatial_disagreements",
    "outside_all_polygons",
    "boundary_ambiguous",
    "text_unresolved",
    "non_agree_geo_precision_1",
    "non_agree_geo_precision_2",
    "non_agree_geo_precision_3",
    "non_agree_in_primary_2017_2024",
    "non_agree_outside_primary_window",
    "max_distance_to_assigned_polygon_m",
    "median_distance_to_assigned_polygon_m",
    "max_distance_to_nearest_polygon_m",
    "median_distance_to_nearest_polygon_m"
  ),
  value = c(
    n_events,
    n_agree,
    n_review,
    fmt_pct(n_review / n_events),
    sum(review_table$spatial_status == "text_spatial_disagree", na.rm = TRUE),
    sum(review_table$spatial_status == "outside_all_polygons", na.rm = TRUE),
    sum(review_table$spatial_status == "boundary_ambiguous", na.rm = TRUE),
    sum(review_table$spatial_status == "text_unresolved", na.rm = TRUE),
    sum(review_table$geo_precision == 1L, na.rm = TRUE),
    sum(review_table$geo_precision == 2L, na.rm = TRUE),
    sum(review_table$geo_precision == 3L, na.rm = TRUE),
    sum(review_table$in_primary_window, na.rm = TRUE),
    sum(!review_table$in_primary_window, na.rm = TRUE),
    max(review_table$distance_to_assigned_polygon_m, na.rm = TRUE),
    median(review_table$distance_to_assigned_polygon_m, na.rm = TRUE),
    max(review_table$distance_to_nearest_polygon_m, na.rm = TRUE),
    median(review_table$distance_to_nearest_polygon_m, na.rm = TRUE)
  )
)

by_status_precision <- review_table |>
  count(
    .data$spatial_status,
    .data$geo_precision,
    .data$review_priority,
    name = "n_events",
    sort = TRUE
  )

by_assigned_district <- review_table |>
  group_by(
    .data$district_id,
    .data$canonical_name
  ) |>
  summarise(
    n_non_agree = n(),
    n_geo_precision_1 =
      sum(.data$geo_precision == 1L, na.rm = TRUE),
    n_text_spatial_disagree =
      sum(.data$spatial_status == "text_spatial_disagree", na.rm = TRUE),
    n_outside_all_polygons =
      sum(.data$spatial_status == "outside_all_polygons", na.rm = TRUE),
    median_distance_to_assigned_m =
      median(.data$distance_to_assigned_polygon_m, na.rm = TRUE),
    max_distance_to_assigned_m =
      max(.data$distance_to_assigned_polygon_m, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(
    desc(.data$n_geo_precision_1),
    desc(.data$n_non_agree),
    .data$district_id
  )

distance_bins <- review_table |>
  count(
    .data$spatial_status,
    .data$geo_precision,
    .data$distance_band,
    name = "n_events"
  ) |>
  arrange(
    .data$spatial_status,
    .data$geo_precision,
    .data$distance_band
  )

# 10. Write diagnostic tables ------------------------------------------------
write_qa(
  summary_table,
  "00_spatial_review_summary.csv"
)

write_qa(
  review_table,
  "01_all_non_agree_events.csv"
)

write_qa(
  review_table |>
    filter(.data$geo_precision == 1L),
  "02_priority_geo_precision_1.csv"
)

write_qa(
  review_table |>
    filter(.data$spatial_status == "text_spatial_disagree"),
  "03_text_spatial_disagreements.csv"
)

write_qa(
  review_table |>
    filter(.data$spatial_status == "outside_all_polygons"),
  "04_outside_all_polygons.csv"
)

write_qa(
  by_assigned_district,
  "05_non_agree_by_assigned_district.csv"
)

write_qa(
  by_status_precision,
  "06_non_agree_by_status_and_precision.csv"
)

write_qa(
  distance_bins,
  "07_distance_bins.csv"
)

# 11. Maps -------------------------------------------------------------------
# Join enriched diagnostic fields back to review geometry.
review_map_sf <- review_sf |>
  select(
    .data$event_id_cnty,
    geometry
  ) |>
  left_join(
    review_table |>
      select(
        .data$event_id_cnty,
        .data$review_priority,
        .data$geo_precision,
        .data$spatial_status,
        .data$canonical_name,
        .data$nearest_district_name,
        .data$distance_to_nearest_polygon_m
      ),
    by = "event_id_cnty"
  )

district_labels <- districts |>
  sf::st_transform(4326) |>
  sf::st_point_on_surface() |>
  mutate(
    label = .data$canonical_name
  )

p_all <- ggplot() +
  geom_sf(
    data = sf::st_transform(districts, 4326),
    fill = NA,
    linewidth = 0.45
  ) +
  geom_sf(
    data = review_map_sf,
    aes(
      shape = .data$spatial_status
    ),
    size = 2.2,
    alpha = 0.8
  ) +
  geom_sf_text(
    data = district_labels,
    aes(label = .data$label),
    size = 2.5,
    check_overlap = TRUE
  ) +
  labs(
    title = "ACLED spatial review: all non-agree events",
    subtitle = paste0(
      "Diagnostic only; district_id is not reassigned. n = ",
      nrow(review_map_sf)
    ),
    shape = "Spatial status",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11)

ggsave(
  filename = file.path(
    output_dir,
    "08_map_all_non_agree.png"
  ),
  plot = p_all,
  width = 10,
  height = 10,
  dpi = 300
)

review_geo1_sf <- review_map_sf |>
  filter(.data$geo_precision == 1L)

if (nrow(review_geo1_sf) > 0L) {
  p_geo1 <- ggplot() +
    geom_sf(
      data = sf::st_transform(districts, 4326),
      fill = NA,
      linewidth = 0.45
    ) +
    geom_sf(
      data = review_geo1_sf,
      aes(
        shape = .data$spatial_status
      ),
      size = 2.8,
      alpha = 0.9
    ) +
    geom_sf_text(
      data = district_labels,
      aes(label = .data$label),
      size = 2.5,
      check_overlap = TRUE
    ) +
    labs(
      title = "Priority review: ACLED geo_precision = 1",
      subtitle = paste0(
        "Only non-agree events; n = ",
        nrow(review_geo1_sf),
        ". No automatic reassignment."
      ),
      shape = "Spatial status",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    filename = file.path(
      output_dir,
      "09_map_geo_precision_1_priority.png"
    ),
    plot = p_geo1,
    width = 10,
    height = 10,
    dpi = 300
  )
}

# 12. Session information ----------------------------------------------------
capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "R_sessionInfo.txt"
  )
)

# 13. Console report ---------------------------------------------------------
message(
  "\nacled spatial review complete"
)

message(
  "Total ACLED events                    : ",
  n_events
)

message(
  "Text/spatial agreements               : ",
  n_agree
)

message(
  "Events requiring spatial review       : ",
  n_review,
  " (",
  fmt_pct(n_review / n_events),
  "%)"
)

message(
  "  text_spatial_disagree               : ",
  sum(
    review_table$spatial_status ==
      "text_spatial_disagree",
    na.rm = TRUE
  )
)

message(
  "  outside_all_polygons                : ",
  sum(
    review_table$spatial_status ==
      "outside_all_polygons",
    na.rm = TRUE
  )
)

message(
  "High priority (geo_precision = 1)     : ",
  sum(
    review_table$geo_precision == 1L,
    na.rm = TRUE
  )
)

message(
  "Medium priority (geo_precision = 2)   : ",
  sum(
    review_table$geo_precision == 2L,
    na.rm = TRUE
  )
)

message(
  "Low priority (geo_precision = 3)      : ",
  sum(
    review_table$geo_precision == 3L,
    na.rm = TRUE
  )
)

message(
  "Review cases inside 2017-2024 window  : ",
  sum(
    review_table$in_primary_window,
    na.rm = TRUE
  )
)

message(
  "Review cases outside 2017-2024 window : ",
  sum(
    !review_table$in_primary_window,
    na.rm = TRUE
  )
)

message(
  "Output directory                      : ",
  output_dir
)

message(
  "\n"
)

message(
  "IMPORTANT: This script is diagnostic only. ",
  "No ACLED district_id was changed.\n"
)

print(
  tibble::as_tibble(summary_table),
  n = Inf,
  width = Inf
)

message("\n--- Priority cases: geo_precision = 1 ---\n")

print(
  review_table |>
    filter(.data$geo_precision == 1L) |>
    select(
      any_of(c(
        "event_id_cnty",
        "event_date",
        "in_primary_window",
        "admin2_original",
        "canonical_name",
        "location",
        "spatial_status",
        "spatial_canonical_name",
        "nearest_district_name",
        "distance_to_assigned_polygon_m",
        "distance_to_nearest_polygon_m",
        "distance_band"
      ))
    ) |>
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)
