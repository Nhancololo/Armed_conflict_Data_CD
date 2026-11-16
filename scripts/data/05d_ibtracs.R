# Download IBTrACS tracks and build the cyclone exposure covariate.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(
  dplyr, tidyr, tibble, readr, lubridate, sf, purrr, slider, here, curl, units
)

start <- as.Date("2017-01-01")
end   <- as.Date("2024-12-31")
crs_m <- 32737

raw_dir <- here::here("data", "raw", "covariates", "ibtracs")
qa_dir <- here::here("outputs", "covariates", "05D_ibtracs")
out_cfg <- here::here("data", "processed", "covariates_ibtracs_2017_2024.csv")
boundary_file <- here::here("data", "processed",
  "cabo_delgado_admin2_canonical.rds")
for (d in c(raw_dir, qa_dir, dirname(out_cfg))) dir.create(d, recursive = TRUE,
  showWarnings = FALSE)

ibtracs_urls <- c(
  "https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.SI.list.v04r01.csv",
  "https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.SI.list.v04r01.csv.gz"
)
raw_csv <- file.path(raw_dir, "ibtracs.SI.list.v04r01.csv")

if (!file.exists(raw_csv)) {
  ok <- FALSE
  for (u in ibtracs_urls) {
    message("Trying IBTrACS: ", u)
    tmp <- tempfile(fileext = if (grepl("gz$", u)) ".gz" else ".csv")
    z <- try(curl::curl_download(u, tmp, quiet = FALSE, mode = "wb"),
      silent = TRUE)
    if (inherits(z, "try-error")) next
    if (grepl("gz$", u)) {
      con_in <- gzfile(tmp, "rb"); con_out <- file(raw_csv, "wb")
      repeat {
        buf <- readBin(con_in, "raw", n = 1024 * 1024)
        if (!length(buf)) break
        writeBin(buf, con_out)
      }
      close(con_in); close(con_out)
    } else file.copy(tmp, raw_csv, overwrite = TRUE)
    unlink(tmp)
    ok <- file.exists(raw_csv) && file.info(raw_csv)$size > 1000
    if (ok) break
  }
  if (!ok) stop("Could not download IBTrACS SI CSV.")
}

# IBTrACS CSV has a units row under the header. Parsing ISO_TIME naturally drops
# it.
x <- readr::read_csv(raw_csv, show_col_types = FALSE, progress = FALSE)
needed <- c("SID", "SEASON", "NAME", "ISO_TIME", "LAT", "LON")
if (!all(needed %in% names(x))) stop("IBTrACS schema changed; missing: ",
  paste(setdiff(needed, names(x)), collapse=", "))

num <- function(v) suppressWarnings(as.numeric(v))
has_usa_wind <- "USA_WIND" %in% names(x)
has_wmo_wind <- "WMO_WIND" %in% names(x)
x <- x |>
  mutate(
    time = suppressWarnings(lubridate::ymd_hms(ISO_TIME, tz = "UTC")),
    date = as.Date(time),
    LAT = num(LAT), LON = num(LON),
    USA_WIND_NUM = if (has_usa_wind) num(USA_WIND) else NA_real_,
    WMO_WIND_NUM = if (has_wmo_wind) num(WMO_WIND) else NA_real_,
    wind_knots = dplyr::coalesce(USA_WIND_NUM, WMO_WIND_NUM),
    year_month = as.Date(lubridate::floor_date(date, "month"))
  ) |>
  filter(!is.na(time), !is.na(LAT), !is.na(LON), date >= start, date <= end)

if (!nrow(x)) stop("No IBTrACS SI observations in study period.")

storms_sf <- sf::st_as_sf(x, coords = c("LON", "LAT"), crs = 4326,
  remove = FALSE) |>
  sf::st_transform(crs_m)
d <- readRDS(boundary_file) |> sf::st_make_valid() |> sf::st_transform(crs_m) |> arrange(district_id)

months <- seq.Date(as.Date("2017-01-01"), as.Date("2024-12-01"), by = "month")

one_month <- function(m) {
  pts <- storms_sf[storms_sf$year_month == m, , drop = FALSE]
  if (!nrow(pts)) {
    return(tibble(
      district_id = d$district_id,
      canonical_name = d$canonical_name,
      year_month = m,
      cyclone_track_points = 0L,
      cyclone_storm_count_100km = 0L,
      cyclone_min_track_distance_km = NA_real_,
      cyclone_max_wind_100km = NA_real_,
      cyclone_within_50km = FALSE,
      cyclone_within_100km = FALSE,
      cyclone_within_200km = FALSE
    ))
  }

  # 17 x n track-point distance matrix, metres.
  dist_m <- units::drop_units(sf::st_distance(d, pts))
  purrr::map_dfr(seq_len(nrow(d)), function(i) {
    di <- as.numeric(dist_m[i, ])
    within100 <- which(di <= 100000)
    within50 <- any(di <= 50000)
    within200 <- any(di <= 200000)
    tibble(
      district_id = d$district_id[i],
      canonical_name = d$canonical_name[i],
      year_month = m,
      cyclone_track_points = nrow(pts),
      cyclone_storm_count_100km = if (length(within100)) dplyr::n_distinct(pts$SID[within100]) else 0L,
      cyclone_min_track_distance_km = min(di, na.rm = TRUE) / 1000,
      cyclone_max_wind_100km = if (length(within100) && any(is.finite(pts$wind_knots[within100]))) max(pts$wind_knots[within100], na.rm = TRUE) else NA_real_,
      cyclone_within_50km = within50,
      cyclone_within_100km = length(within100) > 0L,
      cyclone_within_200km = within200
    )
  })
}

out <- bind_rows(lapply(months, one_month)) |>
  arrange(district_id, year_month) |>
  group_by(district_id) |>
  mutate(
    cyclone_exposure_lag1 = lag(cyclone_within_100km),
    cyclone_post_3m = slider::slide_lgl(lag(cyclone_within_100km), ~ any(.x,
      na.rm = TRUE), .before = 2, .complete = FALSE),
    cyclone_post_6m = slider::slide_lgl(lag(cyclone_within_100km), ~ any(.x,
      na.rm = TRUE), .before = 5, .complete = FALSE)
  ) |>
  ungroup()

if (nrow(out) != 1632L) stop("IBTrACS output is not 17 x 96.")
readr::write_csv(out, out_cfg, na = "")

storm_manifest <- x |>
  sf::st_drop_geometry() |>
  group_by(SID, NAME, SEASON) |>
  summarise(first_date = min(date), last_date = max(date),
    max_wind_knots = suppressWarnings(max(wind_knots,
    na.rm=TRUE)), .groups="drop") |>
  mutate(max_wind_knots = ifelse(is.infinite(max_wind_knots), NA,
    max_wind_knots))
readr::write_csv(storm_manifest, file.path(qa_dir,
  "00_storm_manifest_2017_2024.csv"), na = "")
readr::write_csv(out |> filter(cyclone_within_200km), file.path(qa_dir,
  "01_exposed_cells_within_200km.csv"), na = "")

message("IBTrACS cyclone covariates written: ", out_cfg)
