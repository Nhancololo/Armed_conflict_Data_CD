# Recompute the reporting covariates of 05b_conflict_covariates.R.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  dplyr,
  tidyr,
  tibble,
  readr,
  stringr,
  purrr,
  lubridate,
  here
)

# 1. Paths -------------------------------------------------------------------

panel_file <- here::here(
  "data", "processed",
  "cabo_delgado_district_month_panel_2017_2024.rds"
)

event_file <- here::here(
  "data", "processed",
  "cabo_delgado_event_level_2017_2024.rds"
)

out_file <- here::here(
  "data", "processed",
  "covariates_internal_conflict_2017_2024.csv"
)

qa_dir <- here::here(
  "outputs",
  "covariates",
  "05B_internal_conflict"
)

dir.create(
  qa_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (p in c(panel_file, event_file, out_file)) {
  if (!file.exists(p)) {
    stop("Missing required input: ", p)
  }
}

# 2. Read and validate frozen inputs -----------------------------------------

panel <- readRDS(panel_file) |>
  as_tibble() |>
  mutate(
    year_month = as.Date(year_month)
  ) |>
  arrange(
    district_id,
    time_id
  )

events <- readRDS(event_file) |>
  as_tibble()

old <- readr::read_csv(
  out_file,
  show_col_types = FALSE
) |>
  as_tibble() |>
  mutate(
    year_month = as.Date(year_month)
  ) |>
  arrange(
    district_id,
    time_id
  )

if (
  nrow(panel) != 1632L ||
  n_distinct(panel$district_id) != 17L ||
  n_distinct(panel$time_id) != 96L ||
  anyDuplicated(panel$cell_id)
) {
  stop("Frozen panel is not the expected 17 x 96 design.")
}

if (
  nrow(old) != 1632L ||
  anyDuplicated(old$cell_id)
) {
  stop("Existing 05B output is not the expected 17 x 96 design.")
}

required_event_cols <- c(
  "event_date",
  "district_id",
  "geo_precision",
  "time_precision"
)

missing_event_cols <- setdiff(
  required_event_cols,
  names(events)
)

if (length(missing_event_cols)) {
  stop(
    "Event-level file is missing required column(s): ",
    paste(
      missing_event_cols,
      collapse = ", "
    )
  )
}

has_source <-
  "source" %in%
    names(events)

has_source_scale <-
  "source_scale" %in%
    names(events)

events <- events |>
  mutate(
    event_date =
      as.Date(event_date),

    year_month =
      as.Date(
        lubridate::floor_date(
          event_date,
          "month"
        )
      ),

    source =
      if (has_source) {
        as.character(source)
      } else {
        NA_character_
      },

    source_scale =
      if (has_source_scale) {
        as.character(source_scale)
      } else {
        NA_character_
      },

    geo_precision =
      suppressWarnings(
        as.integer(geo_precision)
      ),

    time_precision =
      suppressWarnings(
        as.integer(time_precision)
      )
  )

unknown_event_ids <- setdiff(
  unique(
    events$district_id[
      !is.na(events$district_id)
    ]
  ),
  unique(panel$district_id)
)

if (length(unknown_event_ids)) {
  stop(
    "Event-level file contains unknown district_id value(s): ",
    paste(
      unknown_event_ids,
      collapse = ", "
    )
  )
}

# 3. Helpers -----------------------------------------------------------------

count_source_tokens <- function(x, n_events) {

  if (n_events == 0L) {
    return(0L)
  }

  x <- as.character(x)

  x <- x[
    !is.na(x) &
      nzchar(
        trimws(x)
      )
  ]

  if (!length(x)) {
    return(NA_integer_)
  }

  tokens <- unlist(
    strsplit(
      x,
      ";",
      fixed = TRUE
    ),
    use.names = FALSE
  )

  tokens <- trimws(tokens)

  tokens <- tokens[
    nzchar(tokens)
  ]

  if (!length(tokens)) {
    return(NA_integer_)
  }

  length(
    unique(tokens)
  )
}

safe_precision_share <- function(x, n_events) {

  if (n_events == 0L) {
    return(NA_real_)
  }

  x <- x[
    !is.na(x)
  ]

  if (!length(x)) {
    return(NA_real_)
  }

  mean(
    x == 1L
  )
}

same_numeric <- function(x, y, tolerance = 1e-12) {

  same_na <-
    identical(
      is.na(x),
      is.na(y)
    )

  if (!same_na) {
    return(FALSE)
  }

  ok <- !is.na(x)

  if (!any(ok)) {
    return(TRUE)
  }

  max(
    abs(
      x[ok] -
        y[ok]
    )
  ) <= tolerance
}

# 4. Rebuild leave-one-district-out reporting environment --------------------

cells <- panel |>
  select(
    district_id,
    year_month,
    time_id,
    cell_id
  )

loo <- purrr::pmap_dfr(
  cells,
  function(
      district_id,
      year_month,
      time_id,
      cell_id
  ) {

    # Rename pmap arguments before entering the dplyr data mask.
    focal_district_id <- district_id
    focal_year_month <- year_month

    e <- events |>
      filter(
        .data$year_month ==
          .env$focal_year_month,

        .data$district_id !=
          .env$focal_district_id
      )

    n_e <- nrow(e)

    tibble(
      district_id =
        focal_district_id,

      year_month =
        focal_year_month,

      time_id =
        time_id,

      cell_id =
        cell_id,

      reporting_events_loo =
        n_e,

      reporting_source_diversity_loo =
        count_source_tokens(
          e$source,
          n_e
        ),

      reporting_source_scale_diversity_loo =
        count_source_tokens(
          e$source_scale,
          n_e
        ),

      reporting_geo1_share_loo =
        safe_precision_share(
          e$geo_precision,
          n_e
        ),

      reporting_time1_share_loo =
        safe_precision_share(
          e$time_precision,
          n_e
        )
    )
  }
)

loo <- loo |>
  arrange(
    district_id,
    time_id
  ) |>
  group_by(district_id) |>
  mutate(
    reporting_events_loo_lag1 =
      lag(reporting_events_loo),

    reporting_source_diversity_loo_lag1 =
      lag(
        reporting_source_diversity_loo
      ),

    reporting_source_scale_diversity_loo_lag1 =
      lag(
        reporting_source_scale_diversity_loo
      ),

    reporting_geo1_share_loo_lag1 =
      lag(
        reporting_geo1_share_loo
      ),

    reporting_time1_share_loo_lag1 =
      lag(
        reporting_time1_share_loo
      )
  ) |>
  ungroup()

# 5. Independent validation against frozen panel counts ----------------------

panel_expected <- panel |>
  group_by(year_month) |>
  mutate(
    province_events_all =
      sum(
        n_events_all,
        na.rm = TRUE
      ),

    province_geo1 =
      sum(
        n_geo_precision_1,
        na.rm = TRUE
      ),

    province_time1 =
      sum(
        n_time_precision_1,
        na.rm = TRUE
      ),

    expected_reporting_events_loo =
      province_events_all -
        n_events_all,

    expected_reporting_geo1_share_loo =
      if_else(
        expected_reporting_events_loo > 0,
        (
          province_geo1 -
            n_geo_precision_1
        ) /
          expected_reporting_events_loo,
        NA_real_
      ),

    expected_reporting_time1_share_loo =
      if_else(
        expected_reporting_events_loo > 0,
        (
          province_time1 -
            n_time_precision_1
        ) /
          expected_reporting_events_loo,
        NA_real_
      )
  ) |>
  ungroup() |>
  arrange(
    district_id,
    time_id
  ) |>
  group_by(district_id) |>
  mutate(
    expected_reporting_events_loo_lag1 =
      lag(
        expected_reporting_events_loo
      ),

    expected_reporting_geo1_share_loo_lag1 =
      lag(
        expected_reporting_geo1_share_loo
      ),

    expected_reporting_time1_share_loo_lag1 =
      lag(
        expected_reporting_time1_share_loo
      )
  ) |>
  ungroup()

check <- loo |>
  select(
    district_id,
    year_month,
    reporting_events_loo_lag1,
    reporting_geo1_share_loo_lag1,
    reporting_time1_share_loo_lag1
  ) |>
  left_join(
    panel_expected |>
      select(
        district_id,
        year_month,
        expected_reporting_events_loo_lag1,
        expected_reporting_geo1_share_loo_lag1,
        expected_reporting_time1_share_loo_lag1
      ),
    by = c(
      "district_id",
      "year_month"
    )
  )

if (
  !same_numeric(
    check$reporting_events_loo_lag1,
    check$expected_reporting_events_loo_lag1
  )
) {
  stop(
    "Rebuilt reporting_events_loo_lag1 does not reproduce the ",
    "independent frozen-panel calculation."
  )
}

if (
  !same_numeric(
    check$reporting_geo1_share_loo_lag1,
    check$expected_reporting_geo1_share_loo_lag1
  )
) {
  stop(
    "Rebuilt reporting_geo1_share_loo_lag1 does not reproduce the ",
    "independent frozen-panel calculation."
  )
}

if (
  !same_numeric(
    check$reporting_time1_share_loo_lag1,
    check$expected_reporting_time1_share_loo_lag1
  )
) {
  stop(
    "Rebuilt reporting_time1_share_loo_lag1 does not reproduce the ",
    "independent frozen-panel calculation."
  )
}

if (
  n_distinct(
    loo$reporting_events_loo_lag1,
    na.rm = TRUE
  ) <= 1L
) {
  stop(
    "reporting_events_loo_lag1 is still degenerate after reconstruction."
  )
}

# 6. Replace only reporting-environment columns in the existing 05B output ---

reporting_cols <- c(
  "reporting_events_loo_lag1",
  "reporting_source_diversity_loo_lag1",
  "reporting_source_scale_diversity_loo_lag1",
  "reporting_geo1_share_loo_lag1",
  "reporting_time1_share_loo_lag1"
)

old_nonreporting <- old |>
  select(
    -any_of(
      reporting_cols
    )
  )

corrected <- old_nonreporting |>
  left_join(
    loo |>
      select(
        district_id,
        year_month,
        all_of(
          reporting_cols
        )
      ),
    by = c(
      "district_id",
      "year_month"
    )
  ) |>
  arrange(
    district_id,
    time_id
  )

if (
  nrow(corrected) != 1632L ||
  anyDuplicated(corrected$cell_id)
) {
  stop(
    "Corrected 05B output corrupted the frozen 17 x 96 key."
  )
}

# Verify every non-reporting column is unchanged.
common_nonreporting <- intersect(
  names(old_nonreporting),
  names(corrected)
)

for (v in common_nonreporting) {

  lhs <- old_nonreporting[[v]]
  rhs <- corrected[[v]]

  same <-
    isTRUE(
      all.equal(
        lhs,
        rhs,
        check.attributes = FALSE
      )
    )

  if (!same) {
    stop(
      "Unexpected change outside reporting block: ",
      v
    )
  }
}

# 7. QA and write ------------------------------------------------------------

qa_summary <- tibble(
  variable = reporting_cols,
  n_missing = vapply(
    corrected[
      reporting_cols
    ],
    function(x) {
      sum(
        is.na(x)
      )
    },
    integer(1)
  ),
  n_unique_nonmissing = vapply(
    corrected[
      reporting_cols
    ],
    function(x) {
      n_distinct(
        x,
        na.rm = TRUE
      )
    },
    integer(1)
  )
)

readr::write_csv(
  qa_summary,
  file.path(
    qa_dir,
    "01_reporting_environment_corrected_QA.csv"
  ),
  na = ""
)

readr::write_csv(
  check,
  file.path(
    qa_dir,
    "02_reporting_environment_panel_crosscheck.csv"
  ),
  na = ""
)

readr::write_csv(
  corrected,
  out_file,
  na = ""
)

message("\n05b reporting fix complete")
message("Output: ", out_file)
message(
  "reporting_events_loo_lag1 unique nonmissing values: ",
  n_distinct(
    corrected$reporting_events_loo_lag1,
    na.rm = TRUE
  )
)
message(
  "reporting_geo1_share_loo_lag1 nonmissing: ",
  sum(
    !is.na(
      corrected$reporting_geo1_share_loo_lag1
    )
  )
)
message(
  "reporting_time1_share_loo_lag1 nonmissing: ",
  sum(
    !is.na(
      corrected$reporting_time1_share_loo_lag1
    )
  )
)
message("Frozen-panel independent cross-check: PASS")
message("\n")
