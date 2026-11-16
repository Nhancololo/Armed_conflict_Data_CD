# Merge the covariate files and check their coverage.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  dplyr,
  tidyr,
  tibble,
  readr,
  stringr,
  lubridate,
  slider,
  purrr,
  here
)

# 1. Configuration -----------------------------------------------------------

require_core <- TRUE

panel_file <- here::here(
  "data", "processed",
  "cabo_delgado_district_month_panel_2017_2024.rds"
)

registry_file <- here::here(
  "data", "processed",
  "cabo_delgado_covariate_registry.csv"
)

gee_dir <- here::here(
  "data", "raw",
  "covariates",
  "gee"
)

out_raw_rds <- here::here(
  "data", "processed",
  "cabo_delgado_panel_with_covariates_raw.rds"
)

out_raw_csv <- here::here(
  "data", "processed",
  "cabo_delgado_panel_with_covariates_raw.csv"
)

out_z_rds <- here::here(
  "data", "processed",
  "cabo_delgado_panel_with_covariates_standardized.rds"
)

out_z_csv <- here::here(
  "data", "processed",
  "cabo_delgado_panel_with_covariates_standardized.csv"
)

qa_dir <- here::here(
  "outputs",
  "covariates",
  "05L_merge_QA"
)

dir.create(
  dirname(out_raw_rds),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  qa_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Frozen outcome panel ----------------------------------------------------

if (!file.exists(panel_file)) {
  stop("Missing frozen panel: ", panel_file)
}

panel <- readRDS(panel_file) |>
  as_tibble() |>
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
  stop(
    "Frozen panel failed the expected 17 x 96 domain/key checks."
  )
}

canonical_ids <- panel |>
  distinct(
    district_id,
    canonical_name
  ) |>
  arrange(district_id)

if (
  nrow(canonical_ids) != 17L ||
  anyDuplicated(canonical_ids$district_id)
) {
  stop("Frozen panel does not contain exactly 17 canonical district IDs.")
}

# 3. Generic helpers ---------------------------------------------------------

assert_known_district_ids <- function(x, label) {
  unknown <- setdiff(
    unique(
      x$district_id[
        !is.na(x$district_id)
      ]
    ),
    canonical_ids$district_id
  )
  
  if (length(unknown)) {
    stop(
      label,
      " contains unknown district_id value(s): ",
      paste(
        unknown,
        collapse = ", "
      )
    )
  }
}

join_dynamic <- function(master, x, label) {
  x <- as_tibble(x)
  
  required <- c(
    "district_id",
    "year_month"
  )
  
  missing_required <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing_required)) {
    stop(
      label,
      " lacks required key column(s): ",
      paste(
        missing_required,
        collapse = ", "
      )
    )
  }
  
  x <- x |>
    mutate(
      year_month =
        as.Date(year_month)
    )
  
  if (
    anyNA(x$district_id) ||
    anyNA(x$year_month)
  ) {
    stop(
      label,
      " contains missing district_id/year_month keys."
    )
  }
  
  assert_known_district_ids(
    x,
    label
  )
  
  if (
    anyDuplicated(
      x[
        c(
          "district_id",
          "year_month"
        )
      ]
    )
  ) {
    stop(
      label,
      " has duplicate district-month keys."
    )
  }
  
  # Frozen panel metadata are authoritative. External copies are discarded.
  drop_meta <- intersect(
    c(
      "canonical_name",
      "time_id",
      "cell_id",
      "month_label",
      "district_index",
      "year",
      "month"
    ),
    names(x)
  )
  
  x <- x |>
    select(
      -any_of(drop_meta)
    )
  
  vars <- setdiff(
    names(x),
    c(
      "district_id",
      "year_month"
    )
  )
  
  overlap <- intersect(
    vars,
    names(master)
  )
  
  if (length(overlap)) {
    stop(
      label,
      " would overwrite existing variable(s): ",
      paste(
        overlap,
        collapse = ", "
      )
    )
  }
  
  out <- master |>
    left_join(
      x,
      by = c(
        "district_id",
        "year_month"
      )
    )
  
  if (
    nrow(out) != nrow(master) ||
    anyDuplicated(out$cell_id)
  ) {
    stop(
      label,
      " corrupted the frozen panel key."
    )
  }
  
  out
}

join_static <- function(master, x, label) {
  x <- as_tibble(x)
  
  if (
    !"district_id" %in% names(x) ||
    anyNA(x$district_id) ||
    anyDuplicated(x$district_id)
  ) {
    stop(
      label,
      " has an invalid static district key."
    )
  }
  
  assert_known_district_ids(
    x,
    label
  )
  
  # --------------------------------------------------------------------------
  # Frozen-panel geometry metadata are authoritative.
  #
  # 05I intentionally carries district_area_km2 for transparency, but the
  # frozen outcome panel already contains the same canonical district area.
  # Before discarding the duplicate copy, verify that the two values agree.
  # This prevents a silent overwrite while also avoiding a false merge error.
  # --------------------------------------------------------------------------
  
  if (
    "district_area_km2" %in% names(x) &&
    "district_area_km2" %in% names(master)
  ) {
    panel_area <- master |>
      distinct(
        district_id,
        panel_district_area_km2 = district_area_km2
      )
    
    if (
      nrow(panel_area) !=
      n_distinct(master$district_id) ||
      anyDuplicated(panel_area$district_id)
    ) {
      stop(
        "Frozen panel has inconsistent district_area_km2 values ",
        "within at least one district."
      )
    }
    
    area_check <- x |>
      select(
        district_id,
        static_district_area_km2 = district_area_km2
      ) |>
      left_join(
        panel_area,
        by = "district_id"
      ) |>
      mutate(
        absolute_difference_km2 =
          abs(
            static_district_area_km2 -
              panel_district_area_km2
          ),
        
        tolerance_km2 =
          pmax(
            1e-6,
            abs(panel_district_area_km2) * 1e-8
          ),
        
        agrees =
          is.finite(static_district_area_km2) &
          is.finite(panel_district_area_km2) &
          absolute_difference_km2 <= tolerance_km2
      )
    
    readr::write_csv(
      area_check,
      file.path(
        qa_dir,
        "06_static_area_consistency.csv"
      ),
      na = ""
    )
    
    if (any(!area_check$agrees)) {
      stop(
        "Static-spatial district_area_km2 does not agree with the ",
        "frozen panel. Inspect 06_static_area_consistency.csv."
      )
    }
  }
  
  # canonical_name, district_index and district_area_km2 are metadata already
  # governed by the frozen panel. Keep only genuinely new static covariates.
  x <- x |>
    select(
      -any_of(
        c(
          "canonical_name",
          "district_index",
          "district_area_km2"
        )
      )
    )
  
  vars <- setdiff(
    names(x),
    "district_id"
  )
  
  overlap <- intersect(
    vars,
    names(master)
  )
  
  if (length(overlap)) {
    stop(
      label,
      " would overwrite existing variable(s): ",
      paste(
        overlap,
        collapse = ", "
      )
    )
  }
  
  out <- master |>
    left_join(
      x,
      by = "district_id"
    )
  
  if (
    nrow(out) != nrow(master) ||
    anyDuplicated(out$cell_id)
  ) {
    stop(
      label,
      " corrupted the frozen panel key."
    )
  }
  
  out
}

resolve_first_existing <- function(
    directory,
    candidates,
    label
) {
  paths <- file.path(
    directory,
    candidates
  )
  
  existing <- paths[
    file.exists(paths)
  ]
  
  if (!length(existing)) {
    return(NA_character_)
  }
  
  if (length(existing) > 1L) {
    message(
      label,
      ": multiple candidate files exist; using highest-priority file: ",
      existing[1]
    )
  }
  
  existing[1]
}

read_gee <- function(path) {
  x <- readr::read_csv(
    path,
    show_col_types = FALSE,
    na = c(
      "",
      "NA",
      "null"
    )
  ) |>
    as_tibble() |>
    select(
      -any_of(
        c(
          "system:index",
          ".geo"
        )
      )
    )
  
  if ("year_month" %in% names(x)) {
    x$year_month <-
      as.Date(
        x$year_month
      )
  }
  
  if ("district_id" %in% names(x)) {
    assert_known_district_ids(
      x,
      basename(path)
    )
  }
  
  x
}

assert_n_rows <- function(
    x,
    expected,
    label
) {
  if (nrow(x) != expected) {
    stop(
      label,
      " has ",
      nrow(x),
      " rows; expected ",
      expected,
      "."
    )
  }
}

# 4. Merge manifest ----------------------------------------------------------

manifest <- tibble(
  block = character(),
  file = character(),
  status = character(),
  n_rows = integer()
)

add_manifest <- function(
    block,
    file,
    status,
    n_rows = NA_integer_
) {
  manifest <<- bind_rows(
    manifest,
    tibble(
      block = block,
      file = file,
      status = status,
      n_rows = n_rows
    )
  )
}

# 5. Start master table with descriptive outcome categories ------------------

master <- panel |>
  mutate(
    # Descriptive response categories only. Never eligible as predictors.
    fatality_threat_level_all =
      case_when(
        fatalities_all < 5 ~ "Low",
        fatalities_all <= 10 ~ "Moderate",
        fatalities_all <= 20 ~ "Substantial",
        fatalities_all <= 50 ~ "Severe",
        TRUE ~ "Critical"
      ),
    
    fatality_threat_level_pv =
      case_when(
        fatalities_political_violence < 5 ~ "Low",
        fatalities_political_violence <= 10 ~ "Moderate",
        fatalities_political_violence <= 20 ~ "Substantial",
        fatalities_political_violence <= 50 ~ "Severe",
        TRUE ~ "Critical"
      )
  )

# 6. R-built dynamic covariate blocks ----------------------------------------

dynamic_blocks <- c(
  internal_conflict =
    here::here(
      "data", "processed",
      "covariates_internal_conflict_2017_2024.csv"
    ),
  
  chirps =
    here::here(
      "data", "processed",
      "covariates_chirps_2017_2024.csv"
    ),
  
  ibtracs =
    here::here(
      "data", "processed",
      "covariates_ibtracs_2017_2024.csv"
    ),
  
  idb =
    here::here(
      "data", "processed",
      "covariates_idb_demography_2017_2024.csv"
    ),
  
  worldpop =
    here::here(
      "data", "processed",
      "covariates_worldpop_population_2017_2024.csv"
    ),
  
  mining =
    here::here(
      "data", "processed",
      "covariates_mining_2017_2024.csv"
    ),
  
  interventions =
    here::here(
      "data", "processed",
      "covariates_interventions_financing_2017_2024.csv"
    ),
  
  iom_dtm =
    here::here(
      "data", "processed",
      "covariates_iom_dtm_2017_2024.csv"
    ),
  
  food_prices =
    here::here(
      "data", "processed",
      "covariates_food_prices_2017_2024.csv"
    )
)

for (block_name in names(dynamic_blocks)) {
  path <- dynamic_blocks[[block_name]]
  
  if (file.exists(path)) {
    x <- readr::read_csv(
      path,
      show_col_types = FALSE
    )
    
    master <- join_dynamic(
      master,
      x,
      block_name
    )
    
    add_manifest(
      block_name,
      path,
      "loaded",
      nrow(x)
    )
  } else {
    add_manifest(
      block_name,
      path,
      "missing_optional_or_not_yet_run"
    )
  }
}

# 7. Static spatial block ----------------------------------------------------

static_file <- here::here(
  "data", "processed",
  "covariates_static_spatial.csv"
)

if (file.exists(static_file)) {
  static_covariates <- readr::read_csv(
    static_file,
    show_col_types = FALSE
  )
  
  master <- join_static(
    master,
    static_covariates,
    "static_spatial"
  )
  
  add_manifest(
    "static_spatial",
    static_file,
    "loaded",
    nrow(static_covariates)
  )
} else {
  add_manifest(
    "static_spatial",
    static_file,
    "missing_not_yet_run"
  )
}

# 8. Google Earth Engine file resolution -------------------------------------

viirs_file <- resolve_first_existing(
  gee_dir,
  c(
    "05H_viirs_monthly_2012_2024.csv"
  ),
  "VIIRS"
)

modis_file <- resolve_first_existing(
  gee_dir,
  c(
    "05H_modis_vi_monthly_2001_2024.csv"
  ),
  "MODIS"
)

era5_file <- resolve_first_existing(
  gee_dir,
  c(
    "05H_era5land_monthly_2017_2024.csv"
  ),
  "ERA5-Land"
)

# Corrected 05H2 outputs have priority. Old 05H files are intentionally not
# accepted as a fallback because QA showed that they were scientifically
# incomplete / incorrectly normalized.
dynamic_world_file <- resolve_first_existing(
  gee_dir,
  c(
    "05H2_dynamicworld_annual_2016_2024.csv",
    "05H2_dynamicworld_annual_2016_2024_corrected.csv"
  ),
  "Dynamic World corrected"
)

hansen_file <- resolve_first_existing(
  gee_dir,
  c(
    "05H2_hansen_annual_2016_2024.csv",
    "05H2_hansen_annual_2016_2024_corrected.csv"
  ),
  "Hansen corrected"
)

# 9. VIIRS nighttime lights --------------------------------------------------

if (!is.na(viirs_file)) {
  viirs <- read_gee(
    viirs_file
  )
  
  assert_n_rows(
    viirs,
    2601L,
    "GEE VIIRS"
  )
  
  viirs <- viirs |>
    mutate(
      calendar_month =
        lubridate::month(
          year_month
        ),
      
      ntl_avg_rad_raw =
        as.numeric(
          ntl_avg_rad
        ),
      
      ntl_cf_cvg =
        as.numeric(
          ntl_cf_cvg
        ),
      
      ntl_negative_raw =
        ntl_avg_rad_raw < 0,
      
      # Radiance used for substantive activity measures.
      ntl_avg_rad_nonneg =
        pmax(
          ntl_avg_rad_raw,
          0
        ),
      
      ntl_log1p =
        log1p(
          ntl_avg_rad_nonneg
        )
    ) |>
    select(
      -ntl_avg_rad
    )
  
  # Pre-study seasonal baseline. Jan-Mar have 2013-2016 observations; Apr-Dec
  # additionally include 2012.
  ntl_baseline <- viirs |>
    filter(
      year_month <
        as.Date("2017-01-01")
    ) |>
    group_by(
      district_id,
      calendar_month
    ) |>
    summarise(
      ntl_log1p_clim_mean =
        mean(
          ntl_log1p,
          na.rm = TRUE
        ),
      
      ntl_log1p_clim_sd =
        sd(
          ntl_log1p,
          na.rm = TRUE
        ),
      
      ntl_baseline_n =
        sum(
          is.finite(
            ntl_log1p
          )
        ),
      
      .groups = "drop"
    )
  
  viirs_study <- viirs |>
    filter(
      year_month >=
        as.Date("2017-01-01"),
      
      year_month <=
        as.Date("2024-12-01")
    ) |>
    left_join(
      ntl_baseline,
      by = c(
        "district_id",
        "calendar_month"
      )
    ) |>
    arrange(
      district_id,
      year_month
    ) |>
    group_by(district_id) |>
    mutate(
      ntl_log1p_anom =
        ntl_log1p -
        ntl_log1p_clim_mean,
      
      ntl_log1p_z =
        if_else(
          is.finite(
            ntl_log1p_clim_sd
          ) &
            ntl_log1p_clim_sd > 0,
          ntl_log1p_anom /
            ntl_log1p_clim_sd,
          NA_real_
        ),
      
      ntl_log1p_lag1 =
        lag(ntl_log1p),
      
      ntl_log1p_anom_lag1 =
        lag(
          ntl_log1p_anom
        ),
      
      ntl_log1p_z_lag1 =
        lag(
          ntl_log1p_z
        ),
      
      ntl_yoy_change =
        ntl_log1p -
        lag(
          ntl_log1p,
          12
        ),
      
      ntl_yoy_change_lag1 =
        lag(
          ntl_yoy_change
        ),
      
      ntl_cf_cvg_lag1 =
        lag(
          ntl_cf_cvg
        )
    ) |>
    ungroup() |>
    select(
      district_id,
      year_month,
      starts_with("ntl_")
    )
  
  if (
    nrow(viirs_study) != 1632L ||
    anyDuplicated(
      viirs_study[
        c(
          "district_id",
          "year_month"
        )
      ]
    )
  ) {
    stop("VIIRS study-period output is not exactly 17 x 96.")
  }
  
  viirs_qa <- viirs |>
    summarise(
      n_rows = n(),
      n_negative_raw =
        sum(
          ntl_negative_raw,
          na.rm = TRUE
        ),
      min_raw =
        min(
          ntl_avg_rad_raw,
          na.rm = TRUE
        ),
      min_cf_cvg =
        min(
          ntl_cf_cvg,
          na.rm = TRUE
        )
    )
  
  readr::write_csv(
    viirs_qa,
    file.path(
      qa_dir,
      "10_viirs_source_QA.csv"
    ),
    na = ""
  )
  
  master <- join_dynamic(
    master,
    viirs_study,
    "GEE_VIIRS"
  )
  
  add_manifest(
    "GEE_VIIRS",
    viirs_file,
    "loaded",
    nrow(viirs)
  )
} else {
  add_manifest(
    "GEE_VIIRS",
    "",
    "missing_not_yet_exported"
  )
}

# 10. MODIS NDVI / EVI -------------------------------------------------------

if (!is.na(modis_file)) {
  modis <- read_gee(
    modis_file
  )
  
  assert_n_rows(
    modis,
    4896L,
    "GEE MODIS"
  )
  
  modis <- modis |>
    mutate(
      calendar_month =
        lubridate::month(
          year_month
        ),
      
      ndvi =
        as.numeric(ndvi),
      
      evi =
        as.numeric(evi)
    )
  
  modis_climatology <- modis |>
    filter(
      year_month >=
        as.Date("2001-01-01"),
      year_month <=
        as.Date("2016-12-01")
    ) |>
    group_by(
      district_id,
      calendar_month
    ) |>
    summarise(
      ndvi_clim_mean =
        mean(
          ndvi,
          na.rm = TRUE
        ),
      
      ndvi_clim_sd =
        sd(
          ndvi,
          na.rm = TRUE
        ),
      
      evi_clim_mean =
        mean(
          evi,
          na.rm = TRUE
        ),
      
      evi_clim_sd =
        sd(
          evi,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  modis_study <- modis |>
    filter(
      year_month >=
        as.Date("2017-01-01"),
      year_month <=
        as.Date("2024-12-01")
    ) |>
    left_join(
      modis_climatology,
      by = c(
        "district_id",
        "calendar_month"
      )
    ) |>
    arrange(
      district_id,
      year_month
    ) |>
    group_by(district_id) |>
    mutate(
      ndvi_anom =
        ndvi -
        ndvi_clim_mean,
      
      ndvi_z =
        if_else(
          is.finite(
            ndvi_clim_sd
          ) &
            ndvi_clim_sd > 0,
          ndvi_anom /
            ndvi_clim_sd,
          NA_real_
        ),
      
      evi_anom =
        evi -
        evi_clim_mean,
      
      evi_z =
        if_else(
          is.finite(
            evi_clim_sd
          ) &
            evi_clim_sd > 0,
          evi_anom /
            evi_clim_sd,
          NA_real_
        ),
      
      ndvi_lag1 =
        lag(ndvi),
      
      ndvi_anom_lag1 =
        lag(
          ndvi_anom
        ),
      
      ndvi_z_lag1 =
        lag(
          ndvi_z
        ),
      
      evi_lag1 =
        lag(evi),
      
      evi_anom_lag1 =
        lag(
          evi_anom
        ),
      
      evi_z_lag1 =
        lag(
          evi_z
        ),
      
      ndvi_prev3_mean =
        slider::slide_dbl(
          lag(ndvi),
          ~mean(.x, na.rm = TRUE),
          .before = 2,
          .complete = TRUE
        ),
      
      ndvi_anom_prev3_mean =
        slider::slide_dbl(
          lag(ndvi_anom),
          ~mean(.x, na.rm = TRUE),
          .before = 2,
          .complete = TRUE
        )
    ) |>
    ungroup() |>
    select(
      district_id,
      year_month,
      ndvi,
      evi,
      starts_with("ndvi_"),
      starts_with("evi_")
    )
  
  master <- join_dynamic(
    master,
    modis_study,
    "GEE_MODIS"
  )
  
  add_manifest(
    "GEE_MODIS",
    modis_file,
    "loaded",
    nrow(modis)
  )
} else {
  add_manifest(
    "GEE_MODIS",
    "",
    "missing_not_yet_exported"
  )
}

# 11. ERA5-Land (secondary climate sensitivity) ------------------------------

if (!is.na(era5_file)) {
  era5 <- read_gee(
    era5_file
  )
  
  assert_n_rows(
    era5,
    1632L,
    "GEE ERA5-Land"
  )
  
  era5_study <- era5 |>
    filter(
      year_month >=
        as.Date("2017-01-01"),
      year_month <=
        as.Date("2024-12-01")
    ) |>
    arrange(
      district_id,
      year_month
    ) |>
    group_by(district_id) |>
    mutate(
      era5_available =
        is.finite(
          era5_temperature_c
        ) &
        is.finite(
          era5_soil_water_l1
        ),
      
      era5_temperature_c_lag1 =
        lag(
          era5_temperature_c
        ),
      
      era5_soil_water_l1_lag1 =
        lag(
          era5_soil_water_l1
        )
    ) |>
    ungroup() |>
    select(
      district_id,
      year_month,
      starts_with("era5_")
    )
  
  era5_missing_by_district <- era5_study |>
    group_by(district_id) |>
    summarise(
      missing_temperature =
        sum(
          is.na(
            era5_temperature_c
          )
        ),
      missing_soil_water =
        sum(
          is.na(
            era5_soil_water_l1
          )
        ),
      .groups = "drop"
    )
  
  readr::write_csv(
    era5_missing_by_district,
    file.path(
      qa_dir,
      "11_era5_missingness_by_district.csv"
    ),
    na = ""
  )
  
  master <- join_dynamic(
    master,
    era5_study,
    "GEE_ERA5"
  )
  
  add_manifest(
    "GEE_ERA5",
    era5_file,
    "loaded_with_explicit_missingness",
    nrow(era5)
  )
} else {
  add_manifest(
    "GEE_ERA5",
    "",
    "missing_optional"
  )
}

# 12. Annual Earth Engine helper ---------------------------------------------

annual_with_lag_to_monthly <- function(
    x,
    value_cols,
    label
) {
  if (
    !"district_id" %in% names(x)
  ) {
    stop(
      label,
      " lacks district_id."
    )
  }
  
  if (!"year" %in% names(x)) {
    if (!"year_month" %in% names(x)) {
      stop(
        label,
        " lacks year and year_month."
      )
    }
    
    x$year <-
      lubridate::year(
        x$year_month
      )
  }
  
  x <- x |>
    mutate(
      year =
        as.integer(year)
    )
  
  annual <- x |>
    select(
      district_id,
      year,
      all_of(value_cols)
    ) |>
    distinct()
  
  if (
    anyDuplicated(
      annual[
        c(
          "district_id",
          "year"
        )
      ]
    )
  ) {
    stop(
      label,
      " has duplicate district-year keys."
    )
  }
  
  # Compute the lag before filtering 2017-2024 so 2016 correctly supplies the
  # one-year lag for 2017.
  annual <- annual |>
    arrange(
      district_id,
      year
    ) |>
    group_by(district_id) |>
    mutate(
      across(
        all_of(value_cols),
        ~lag(.x),
        .names = "{.col}_lag1y"
      )
    ) |>
    ungroup() |>
    filter(
      year >= 2017,
      year <= 2024
    )
  
  domain <- master |>
    distinct(
      district_id,
      year_month
    ) |>
    mutate(
      year =
        lubridate::year(
          year_month
        )
    )
  
  domain |>
    left_join(
      annual,
      by = c(
        "district_id",
        "year"
      )
    ) |>
    select(-year)
}

# 13. Corrected Dynamic World ------------------------------------------------

if (!is.na(dynamic_world_file)) {
  dynamic_world <- read_gee(
    dynamic_world_file
  )
  
  assert_n_rows(
    dynamic_world,
    153L,
    "GEE corrected Dynamic World"
  )
  
  dw_class_vars <- c(
    "water_share",
    "forest_share",
    "grass_share",
    "flooded_vegetation_share",
    "cropland_share",
    "shrub_share",
    "built_share",
    "bare_share",
    "snow_ice_share"
  )
  
  missing_dw <- setdiff(
    dw_class_vars,
    names(dynamic_world)
  )
  
  if (length(missing_dw)) {
    stop(
      "Corrected Dynamic World file is missing class share(s): ",
      paste(
        missing_dw,
        collapse = ", "
      )
    )
  }
  
  dynamic_world <- dynamic_world |>
    mutate(
      dw_recomputed_share_sum =
        rowSums(
          pick(
            all_of(
              dw_class_vars
            )
          ),
          na.rm = FALSE
        )
    )
  
  if (
    any(
      abs(
        dynamic_world$dw_recomputed_share_sum -
        1
      ) > 1e-6,
      na.rm = TRUE
    )
  ) {
    stop(
      "Corrected Dynamic World class shares do not sum to 1 within tolerance."
    )
  }
  
  readr::write_csv(
    dynamic_world |>
      select(
        district_id,
        year,
        any_of(
          c(
            "dw_class_share_sum",
            "dw_recomputed_share_sum",
            "dw_valid_area_fraction_of_geometry",
            "dw_source_image_count"
          )
        )
      ),
    file.path(
      qa_dir,
      "12_dynamicworld_source_QA.csv"
    ),
    na = ""
  )
  
  dw_monthly <- annual_with_lag_to_monthly(
    dynamic_world,
    dw_class_vars,
    "GEE_DynamicWorld"
  )
  
  master <- join_dynamic(
    master,
    dw_monthly,
    "GEE_DynamicWorld"
  )
  
  add_manifest(
    "GEE_DynamicWorld",
    dynamic_world_file,
    "loaded",
    nrow(dynamic_world)
  )
} else {
  add_manifest(
    "GEE_DynamicWorld",
    "",
    "missing_corrected_05H2_export"
  )
}

# 14. Corrected Hansen Global Forest Change ----------------------------------

if (!is.na(hansen_file)) {
  hansen <- read_gee(
    hansen_file
  )
  
  assert_n_rows(
    hansen,
    153L,
    "GEE corrected Hansen"
  )
  
  # Registry-compatible aliases retain the old variable names, but now use the
  # corrected mapped-land denominator explicitly.
  hansen <- hansen |>
    mutate(
      forest_loss_share_annual =
        forest_loss_share_land_annual,
      forest_loss_share_cumulative =
        forest_loss_share_land_cumulative
    )
  
  hansen_vars <- c(
    "treecover2000_pct",
    "forest2000_share_30",
    "forest_loss_share_land_annual",
    "forest_loss_share_land_cumulative",
    "forest_loss_share_forest2000_30_annual",
    "forest_loss_share_forest2000_30_cumulative",
    "forest_loss_share_annual",
    "forest_loss_share_cumulative"
  )
  
  missing_hansen <- setdiff(
    hansen_vars,
    names(hansen)
  )
  
  if (length(missing_hansen)) {
    stop(
      "Corrected Hansen file is missing variable(s): ",
      paste(
        missing_hansen,
        collapse = ", "
      )
    )
  }
  
  hansen_ordered <- hansen |>
    arrange(
      district_id,
      year
    ) |>
    group_by(district_id) |>
    mutate(
      land_cum_diff =
        forest_loss_share_land_cumulative -
        lag(
          forest_loss_share_land_cumulative
        ),
      
      forest30_cum_diff =
        forest_loss_share_forest2000_30_cumulative -
        lag(
          forest_loss_share_forest2000_30_cumulative
        )
    ) |>
    ungroup()
  
  if (
    any(
      hansen_ordered$land_cum_diff < -1e-10,
      na.rm = TRUE
    ) ||
    any(
      hansen_ordered$forest30_cum_diff < -1e-10,
      na.rm = TRUE
    )
  ) {
    stop(
      "Corrected Hansen cumulative forest-loss series is not monotone."
    )
  }
  
  readr::write_csv(
    hansen_ordered |>
      summarise(
        max_land_annual =
          max(
            forest_loss_share_land_annual,
            na.rm = TRUE
          ),
        max_land_cumulative =
          max(
            forest_loss_share_land_cumulative,
            na.rm = TRUE
          ),
        max_forest30_annual =
          max(
            forest_loss_share_forest2000_30_annual,
            na.rm = TRUE
          ),
        max_forest30_cumulative =
          max(
            forest_loss_share_forest2000_30_cumulative,
            na.rm = TRUE
          )
      ),
    file.path(
      qa_dir,
      "13_hansen_source_QA.csv"
    ),
    na = ""
  )
  
  hansen_monthly <- annual_with_lag_to_monthly(
    hansen,
    hansen_vars,
    "GEE_Hansen"
  )
  
  master <- join_dynamic(
    master,
    hansen_monthly,
    "GEE_Hansen"
  )
  
  add_manifest(
    "GEE_Hansen",
    hansen_file,
    "loaded",
    nrow(hansen)
  )
} else {
  add_manifest(
    "GEE_Hansen",
    "",
    "missing_corrected_05H2_export"
  )
}

# 15. Cross-source derived variables -----------------------------------------

# Prefer IDB population for continuity. WorldPop remains an
# independent population sensitivity; no source is silently substituted.
if (
  all(
    c(
      "idp_stock",
      "population"
    ) %in% names(master)
  )
) {
  master <- master |>
    mutate(
      idp_rate_per_1000 =
        1000 *
        idp_stock /
        population,
      
      returnee_rate_per_1000 =
        1000 *
        returnee_stock /
        population
    )
} else if (
  all(
    c(
      "idp_stock",
      "population_worldpop"
    ) %in% names(master)
  )
) {
  master <- master |>
    mutate(
      idp_rate_per_1000 =
        1000 *
        idp_stock /
        population_worldpop,
      
      returnee_rate_per_1000 =
        1000 *
        returnee_stock /
        population_worldpop
    )
}

if ("population" %in% names(master)) {
  master$log_population <-
    log1p(
      master$population
    )
}

if (
  "population_worldpop" %in%
  names(master)
) {
  master$log_population_worldpop <-
    log1p(
      master$population_worldpop
    )
}

# 16. Standardization --------------------------------------------------------

# Only columns added after the frozen panel are eligible for standardization.
# This prevents outcome/event-count fields from becoming candidate predictors.
covariate_cols <- setdiff(
  names(master),
  names(panel)
)

numeric_covariate_cols <- covariate_cols[
  vapply(
    master[covariate_cols],
    is.numeric,
    logical(1)
  )
]

continuous <- numeric_covariate_cols[
  vapply(
    master[
      numeric_covariate_cols
    ],
    function(x) {
      finite_values <-
        unique(
          x[
            is.finite(x)
          ]
        )
      
      length(finite_values) > 2L &&
        is.finite(
          sd(
            x,
            na.rm = TRUE
          )
        ) &&
        sd(
          x,
          na.rm = TRUE
        ) > 0
    },
    logical(1)
  )
]

standardized <- master

# Generic full-panel standardization must never overwrite a scientifically
# meaningful pre-existing column such as ntl_log1p_z, ndvi_z, or evi_z.
#
# Default generic name:
#   <variable>_z
#
# If that name already exists in the RAW merged panel, preserve the existing
# column and write the generic full-panel z-score as:
#   <variable>_global_z
standardization_map <- tibble(
  variable = continuous,
  preferred_name =
    paste0(
      continuous,
      "_z"
    )
) |>
  mutate(
    name_collision =
      preferred_name %in%
      names(master),
    
    standardized_name =
      if_else(
        name_collision,
        paste0(
          variable,
          "_global_z"
        ),
        preferred_name
      )
  )

if (
  anyDuplicated(
    standardization_map$standardized_name
  )
) {
  stop(
    "Generic standardization generated duplicate output names."
  )
}

for (i in seq_len(nrow(standardization_map))) {
  variable <-
    standardization_map$variable[i]
  
  output_name <-
    standardization_map$standardized_name[i]
  
  mu <- mean(
    master[[variable]],
    na.rm = TRUE
  )
  
  sigma <- sd(
    master[[variable]],
    na.rm = TRUE
  )
  
  standardized[[output_name]] <- (
    master[[variable]] -
      mu
  ) / sigma
}

readr::write_csv(
  standardization_map,
  file.path(
    qa_dir,
    "07_standardization_name_map.csv"
  ),
  na = ""
)

# Sensitivity-only state-presence index from post-onset/current snapshots.
state_components <- intersect(
  c(
    "school_density_z",
    "hospital_density_z",
    "cell_tower_density_z"
  ),
  names(standardized)
)

if (length(state_components) >= 2L) {
  standardized$state_presence_index_z <-
    rowMeans(
      as.data.frame(
        standardized[
          state_components
        ]
      ),
      na.rm = TRUE
    )
  
  standardized$state_presence_index_z[
    !is.finite(
      standardized$state_presence_index_z
    )
  ] <- NA_real_
}

# Protected source-specific z/anomaly variables must be bit-for-bit preserved.
protected_preexisting_z <- intersect(
  c(
    "ntl_log1p_z",
    "ndvi_z",
    "evi_z"
  ),
  names(master)
)

for (v in protected_preexisting_z) {
  same_values <-
    isTRUE(
      all.equal(
        master[[v]],
        standardized[[v]],
        check.attributes = FALSE
      )
    )
  
  if (!same_values) {
    stop(
      "Standardization overwrote protected scientific column: ",
      v
    )
  }
}

# 17. Final key QA -----------------------------------------------------------

if (
  nrow(master) != 1632L ||
  anyDuplicated(master$cell_id)
) {
  stop(
    "Merged panel corrupted the frozen cell key."
  )
}

months_per_district <- master |>
  count(
    district_id,
    name = "n_months"
  )

if (
  any(
    months_per_district$n_months != 96L
  )
) {
  stop(
    "Merged panel lost months for at least one district."
  )
}

# 18. Missingness / correlations / registry availability ---------------------

missingness <- tibble(
  variable = covariate_cols,
  
  n_missing =
    vapply(
      master[
        covariate_cols
      ],
      function(x) {
        sum(
          is.na(x)
        )
      },
      integer(1)
    ),
  
  pct_missing =
    100 *
    n_missing /
    nrow(master),
  
  n_unique_nonmissing =
    vapply(
      master[
        covariate_cols
      ],
      function(x) {
        length(
          unique(
            x[
              !is.na(x)
            ]
          )
        )
      },
      integer(1)
    )
) |>
  arrange(
    desc(pct_missing),
    variable
  )

numeric_covariates <- intersect(
  covariate_cols,
  names(master)[
    vapply(
      master,
      is.numeric,
      logical(1)
    )
  ]
)

correlation_pairs <- tibble(
  var1 = character(),
  var2 = character(),
  correlation = numeric(),
  n_complete = integer()
)

if (
  length(
    numeric_covariates
  ) >= 2L
) {
  pairs <- combn(
    numeric_covariates,
    2,
    simplify = FALSE
  )
  
  correlation_pairs <- purrr::map_dfr(
    pairs,
    function(pair) {
      ok <-
        is.finite(
          master[[pair[1]]]
        ) &
        is.finite(
          master[[pair[2]]]
        )
      
      if (sum(ok) < 30L) {
        return(NULL)
      }
      
      r <- suppressWarnings(
        cor(
          master[[pair[1]]][ok],
          master[[pair[2]]][ok]
        )
      )
      
      if (
        !is.finite(r) ||
        abs(r) < 0.90
      ) {
        return(NULL)
      }
      
      tibble(
        var1 = pair[1],
        var2 = pair[2],
        correlation = r,
        n_complete = sum(ok)
      )
    }
  ) |>
    arrange(
      desc(
        abs(correlation)
      )
    )
}

readr::write_csv(
  manifest,
  file.path(
    qa_dir,
    "00_merge_manifest.csv"
  ),
  na = ""
)

readr::write_csv(
  missingness,
  file.path(
    qa_dir,
    "01_covariate_missingness.csv"
  ),
  na = ""
)

readr::write_csv(
  correlation_pairs,
  file.path(
    qa_dir,
    "02_extreme_correlations_abs_ge_0.90.csv"
  ),
  na = ""
)

readr::write_csv(
  standardization_map |>
    select(
      variable,
      standardized_name,
      name_collision
    ),
  file.path(
    qa_dir,
    "03_standardized_continuous_variables.csv"
  ),
  na = ""
)

if (file.exists(registry_file)) {
  registry <- readr::read_csv(
    registry_file,
    show_col_types = FALSE
  )
  
  availability <- registry |>
    transmute(
      candidate_id,
      variable,
      block,
      priority,
      role,
      available_exact =
        variable %in%
        names(master),
      available_z =
        vapply(
          variable,
          function(v) {
            hit <- standardization_map |>
              filter(
                .data$variable == v
              )
            
            if (!nrow(hit)) {
              return(FALSE)
            }
            
            hit$standardized_name[1] %in%
              names(standardized)
          },
          logical(1)
        ),
      causal_note
    )
  
  readr::write_csv(
    availability,
    file.path(
      qa_dir,
      "04_registry_availability.csv"
    ),
    na = ""
  )
  
  unregistered <- tibble(
    variable = setdiff(
      covariate_cols,
      registry$variable
    )
  )
  
  readr::write_csv(
    unregistered,
    file.path(
      qa_dir,
      "05_covariates_not_yet_in_registry.csv"
    ),
    na = ""
  )
}

# 19. Core completeness gate -------------------------------------------------

if (require_core) {
  core_blocks <- c(
    "internal_conflict",
    "chirps",
    "ibtracs",
    "interventions",
    "static_spatial",
    "GEE_VIIRS",
    "GEE_MODIS",
    "GEE_DynamicWorld",
    "GEE_Hansen"
  )
  
  missing_core <- manifest |>
    filter(
      block %in%
        core_blocks,
      status != "loaded"
    )
  
  population_available <-
    any(
      manifest$block %in%
        c(
          "idb",
          "worldpop"
        ) &
        manifest$status ==
        "loaded"
    )
  
  if (
    nrow(missing_core) ||
    !population_available
  ) {
    readr::write_csv(
      missing_core,
      file.path(
        qa_dir,
        "STOP_missing_core_blocks.csv"
      ),
      na = ""
    )
    
    stop(
      "Core covariate universe is incomplete. Inspect 00_merge_manifest.csv ",
      "and STOP_missing_core_blocks.csv. At least IDB or WorldPop must also ",
      "be loaded."
    )
  }
}

# 20. Save final analysis tables ---------------------------------------------

saveRDS(
  master,
  out_raw_rds,
  compress = "xz"
)

readr::write_csv(
  master,
  out_raw_csv,
  na = ""
)

saveRDS(
  standardized,
  out_z_rds,
  compress = "xz"
)

readr::write_csv(
  standardized,
  out_z_csv,
  na = ""
)

message(
  "Merged raw panel: ",
  out_raw_rds
)

message(
  "Standardized panel: ",
  out_z_rds
)

message(
  "Rows: ",
  nrow(master),
  "; external/derived covariate columns: ",
  length(covariate_cols)
)

message(
  "No external missing values were imputed. ",
  "Missingness remains explicit in QA."
)
