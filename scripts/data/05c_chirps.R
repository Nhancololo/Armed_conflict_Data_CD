# Download CHIRPS rainfall and build the precipitation anomaly.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(
  terra, sf, dplyr, tidyr, tibble, readr, lubridate, slider, here, curl
)

baseline_start <- as.Date("1981-01-01")
baseline_end   <- as.Date("2016-12-01")
study_start    <- as.Date("2017-01-01")
study_end      <- as.Date("2024-12-01")
allow_download_fallback <- TRUE

boundary_file <- here::here("data", "processed",
  "cabo_delgado_admin2_canonical.rds")
raw_dir <- here::here("data", "raw", "covariates", "chirps")
cache_dir <- file.path(raw_dir, "district_month_cache")
out <- here::here("data", "processed", "covariates_chirps_2017_2024.csv")
qa_dir <- here::here("outputs", "covariates", "05C_chirps")
for (d in c(raw_dir, cache_dir, dirname(out), qa_dir)) dir.create(d,
  recursive = TRUE, showWarnings = FALSE)

if (!file.exists(boundary_file)) stop("Missing canonical boundary.")
d <- readRDS(boundary_file) |> sf::st_make_valid() |> sf::st_transform(4326)
d <- d |> arrange(district_id)
v <- terra::vect(d)

month_seq <- seq.Date(baseline_start, study_end, by = "month")
base_url <- "https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/cogs"

extract_month <- function(date) {
  ym <- format(date, "%Y.%m")
  ym_dash <- format(date, "%Y-%m")
  cache <- file.path(cache_dir, paste0("chirps_", ym_dash, ".csv"))
  if (file.exists(cache)) return(readr::read_csv(cache, show_col_types = FALSE))

  url <- paste0(base_url, "/chirps-v2.0.", ym, ".cog")
  remote <- paste0("/vsicurl/", url)

  get_values <- function(rpath) {
    r <- terra::rast(rpath)
    # Crop first to reduce I/O. A small geographic padding prevents edge
    # effects.
    e <- terra::ext(v)
    e <- terra::ext(
      terra::xmin(e) - 0.1, terra::xmax(e) + 0.1,
      terra::ymin(e) - 0.1, terra::ymax(e) + 0.1
    )
    rc <- terra::crop(r, e)

    # Exact polygon-cell fractions are used explicitly so boundary pixels are
    # area-weighted rather than treated as fully inside a district.
    z <- terra::extract(rc, v, exact = TRUE)
    value_col <- setdiff(names(z), c("ID", "fraction", "weight"))
    if (length(value_col) != 1L) stop("Unexpected CHIRPS extraction columns for ", ym_dash)
    wt_col <- intersect(c("fraction", "weight"), names(z))
    if (!length(wt_col)) {
      z$.w <- 1
    } else {
      z$.w <- z[[wt_col[1]]]
    }
    z$.value <- as.numeric(z[[value_col]])
    agg <- z |>
      as_tibble() |>
      group_by(ID) |>
      summarise(
        precip_mm = if (all(is.na(.data$.value))) NA_real_ else
          weighted.mean(.data$.value, w = .data$.w, na.rm = TRUE),
        .groups = "drop"
      )

    tibble(
      district_id = d$district_id[agg$ID],
      canonical_name = d$canonical_name[agg$ID],
      year_month = as.Date(paste0(ym_dash, "-01")),
      precip_mm = agg$precip_mm
    )
  }

  ans <- try(get_values(remote), silent = TRUE)
  if (inherits(ans, "try-error") && allow_download_fallback) {
    message("Remote COG stream failed for ", ym_dash,
      "; downloading temporary full COG.")
    tmp <- tempfile(fileext = ".cog")
    curl::curl_download(url, tmp, quiet = FALSE, mode = "wb")
    on.exit(unlink(tmp), add = TRUE)
    ans <- get_values(tmp)
  }
  if (inherits(ans, "try-error")) stop("CHIRPS failed for ", ym_dash, ": ",
    as.character(ans))
  if (nrow(ans) != 17L || anyNA(ans$precip_mm)) stop("CHIRPS extraction incomplete for ", ym_dash)
  readr::write_csv(ans, cache, na = "")
  ans
}

message("Extracting/caching ", length(month_seq), " CHIRPS months (1981-2024).")
all_monthly <- bind_rows(lapply(month_seq, extract_month))

clim <- all_monthly |>
  filter(year_month >= baseline_start, year_month <= baseline_end) |>
  mutate(calendar_month = month(year_month)) |>
  group_by(district_id, calendar_month) |>
  summarise(
    precip_clim_mean = mean(precip_mm, na.rm = TRUE),
    precip_clim_sd = sd(precip_mm, na.rm = TRUE),
    n_clim_years = n(),
    .groups = "drop"
  )

if (any(clim$n_clim_years != 36L)) stop("Expected 36 baseline years for every district/calendar-month.")

study <- all_monthly |>
  filter(year_month >= study_start, year_month <= study_end) |>
  mutate(calendar_month = month(year_month)) |>
  left_join(clim, by = c("district_id", "calendar_month")) |>
  mutate(
    chirps_precip_mm = precip_mm,
    chirps_precip_anom_mm = precip_mm - precip_clim_mean,
    chirps_precip_z = if_else(precip_clim_sd > 0,
      chirps_precip_anom_mm / precip_clim_sd, NA_real_),
    chirps_drought_z_le_m1 = chirps_precip_z <= -1,
    chirps_severe_drought_z_le_m1_5 = chirps_precip_z <= -1.5,
    chirps_wet_z_ge_1 = chirps_precip_z >= 1
  ) |>
  arrange(district_id, year_month) |>
  group_by(district_id) |>
  mutate(
    chirps_precip_lag1_mm = lag(chirps_precip_mm),
    chirps_precip_z_lag1 = lag(chirps_precip_z),
    chirps_prev3_mm_sum = slider::slide_dbl(lag(chirps_precip_mm), ~ sum(.x,
      na.rm = TRUE), .before = 2, .complete = TRUE),
    chirps_prev3_z_mean = slider::slide_dbl(lag(chirps_precip_z), ~ mean(.x,
      na.rm = TRUE), .before = 2, .complete = TRUE),
    chirps_prev3_drought = chirps_prev3_z_mean <= -1
  ) |>
  ungroup() |>
  select(
    district_id, canonical_name, year_month,
    chirps_precip_mm, chirps_precip_anom_mm, chirps_precip_z,
    chirps_drought_z_le_m1, chirps_severe_drought_z_le_m1_5, chirps_wet_z_ge_1,
    chirps_precip_lag1_mm, chirps_precip_z_lag1,
    chirps_prev3_mm_sum, chirps_prev3_z_mean, chirps_prev3_drought,
    precip_clim_mean, precip_clim_sd
  )

if (nrow(study) != 1632L) stop("CHIRPS output is not 17 x 96.")
readr::write_csv(study, out, na = "")
readr::write_csv(clim, file.path(qa_dir, "00_fixed_climatology_1981_2016.csv"),
  na = "")
readr::write_csv(
  study |> summarise(across(starts_with("chirps_"),
    list(missing = ~sum(is.na(.x)), min = ~min(.x, na.rm=TRUE), max = ~max(.x,
    na.rm=TRUE)))),
  file.path(qa_dir, "01_summary.csv"), na = ""
)
message("CHIRPS covariates written: ", out)
