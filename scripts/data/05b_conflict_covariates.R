# Build lagged conflict and reporting covariates.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(
  dplyr, tidyr, tibble, readr, stringr, purrr, slider, spdep, lubridate, here
)
`%||%` <- function(x, y) if (is.null(x)) y else x

panel_file <- here::here("data", "processed",
  "cabo_delgado_district_month_panel_2017_2024.rds")
event_file <- here::here("data", "processed",
  "cabo_delgado_event_level_2017_2024.rds")
adj_file <- here::here("data", "processed",
  "cabo_delgado_adjacency_topological.rds")
out_cfg <- here::here("data", "processed",
  "covariates_internal_conflict_2017_2024.csv")
qa_dir <- here::here("outputs", "covariates", "05B_internal_conflict")
dir.create(dirname(out_cfg), recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

for (p in c(panel_file, event_file,
  adj_file)) if (!file.exists(p)) stop("Missing input: ", p)
panel <- readRDS(panel_file) |> as_tibble()
events <- readRDS(event_file) |> as_tibble()
adj <- readRDS(adj_file)

if (nrow(panel) != 1632L || n_distinct(panel$district_id) != 17L || n_distinct(panel$time_id) != 96L) {
  stop("Panel domain is not the frozen 17 x 96 design.")
}
if (anyDuplicated(panel$cell_id)) stop("Duplicate cell_id in panel.")

# Helper: sum of prior k months, excluding current month.
past_sum <- function(x, k) {
  slider::slide_dbl(dplyr::lag(x), ~ sum(.x, na.rm = TRUE), .before = k - 1L,
    .complete = TRUE)
}
past_mean <- function(x, k) {
  slider::slide_dbl(dplyr::lag(x), ~ mean(.x, na.rm = TRUE), .before = k - 1L,
    .complete = TRUE)
}

# Number of months since the most recent PRIOR month with x > 0. Current-month
# information is never used to define this predictor. NA means no prior positive
# month has yet been observed in the study window.
months_since_last_positive <- function(x) {
  out <- rep(NA_integer_, length(x))
  last <- NA_integer_
  for (i in seq_along(x)) {
    if (!is.na(last)) out[i] <- i - last
    if (is.finite(x[i]) && x[i] > 0) last <- i
  }
  out
}

hist <- panel |>
  arrange(district_id, time_id) |>
  group_by(district_id) |>
  mutate(
    events_all_lag1 = lag(n_events_all, 1),
    fatalities_all_lag1 = lag(fatalities_all, 1),
    events_pv_lag1 = lag(n_events_political_violence, 1),
    fatalities_pv_lag1 = lag(fatalities_political_violence, 1),
    events_all_prev3_sum = past_sum(n_events_all, 3),
    fatalities_all_prev3_sum = past_sum(fatalities_all, 3),
    events_all_prev6_sum = past_sum(n_events_all, 6),
    fatalities_all_prev6_sum = past_sum(fatalities_all, 6),
    events_pv_prev3_sum = past_sum(n_events_political_violence, 3),
    fatalities_pv_prev3_sum = past_sum(fatalities_political_violence, 3),

    # Duration dependence / conflict recency.
    months_since_last_event = months_since_last_positive(n_events_all),
    months_since_last_pv_event = months_since_last_positive(n_events_political_violence),

    # Lagged composition of prior violence/unrest. These are endogenous process
    # controls, not exogenous causal exposures.
    battles_prev3_sum = past_sum(n_battles, 3),
    violence_against_civilians_prev3_sum = past_sum(n_violence_against_civilians, 3),
    explosions_remote_prev3_sum = past_sum(n_explosions_remote_violence, 3),
    riots_prev3_sum = past_sum(n_riots, 3),
    protests_prev3_sum = past_sum(n_protests, 3),
    strategic_events_prev3_sum = past_sum(n_strategic_events, 3),
    civilian_targeting_prev3_sum = past_sum(n_civilian_targeting_events, 3)
  ) |>
  ungroup()

# Neighbor spillovers. Use connected model graph for primary, raw Queen as
# sensitivity.
nb_model <- adj$nb_queen_connected %||% adj$nb_queen
nb_raw <- adj$nb_queen_raw %||% nb_model
if (is.null(nb_model)) stop("Adjacency RDS has no nb object.")

id_order <- adj$district_index$district_id
if (is.null(id_order)) stop("Adjacency RDS lacks district_index.")

nb_sum <- function(df, nb, value_col) {
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  stopifnot(nrow(W) == length(id_order))
  out <- rep(NA_real_, nrow(df))
  for (tt in sort(unique(df$time_id))) {
    idx <- which(df$time_id == tt)
    block <- df[idx, ] |> arrange(match(district_id, id_order))
    if (!identical(as.character(block$district_id),
      as.character(id_order))) stop("District order mismatch in neighbor calculation.")
    v <- block[[value_col]]
    z <- as.numeric(W %*% v)
    out[idx] <- z[match(df$district_id[idx], block$district_id)]
  }
  out
}

hist <- hist |>
  arrange(time_id, match(district_id, id_order))

hist$neighbor_events_lag1 <- nb_sum(hist, nb_model, "events_all_lag1")
hist$neighbor_fatalities_lag1 <- nb_sum(hist, nb_model, "fatalities_all_lag1")
hist$neighbor_events_prev3_sum <- nb_sum(hist, nb_model, "events_all_prev3_sum")
hist$neighbor_fatalities_prev3_sum <- nb_sum(hist, nb_model,
  "fatalities_all_prev3_sum")
hist$neighbor_events_pv_lag1 <- nb_sum(hist, nb_model, "events_pv_lag1")
hist$neighbor_fatalities_pv_lag1 <- nb_sum(hist, nb_model, "fatalities_pv_lag1")
hist$neighbor_events_pv_prev3_sum <- nb_sum(hist, nb_model,
  "events_pv_prev3_sum")
hist$neighbor_fatalities_pv_prev3_sum <- nb_sum(hist, nb_model,
  "fatalities_pv_prev3_sum")
hist$neighbor_events_lag1_rawqueen <- nb_sum(hist, nb_raw, "events_all_lag1")
hist$neighbor_fatalities_lag1_rawqueen <- nb_sum(hist, nb_raw,
  "fatalities_all_lag1")

# Reporting environment ---------------------------------------------------------
# For each district-month, aggregate ACLED evidence OUTSIDE the focal district.
# The resulting environment is lagged by one month before use.
has_source <- "source" %in% names(events)
has_source_scale <- "source_scale" %in% names(events)

events <- events |>
  mutate(
    event_date = as.Date(event_date),
    year_month = as.Date(lubridate::floor_date(event_date, "month")),
    source = if (has_source) as.character(source) else NA_character_,
    source_scale = if (has_source_scale) as.character(source_scale) else NA_character_
  )

count_source_tokens <- function(x) {
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(0L)
  tok <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  tok <- trimws(tok)
  length(unique(tok[nzchar(tok)]))
}

cells <- panel |> select(district_id, year_month, time_id, cell_id)
loo <- purrr::pmap_dfr(cells, function(district_id, year_month, time_id,
  cell_id) {
  e <- events |> filter(.data$year_month == year_month,
    .data$district_id != district_id)
  tibble(
    district_id = district_id,
    year_month = year_month,
    time_id = time_id,
    cell_id = cell_id,
    reporting_events_loo = nrow(e),
    reporting_source_diversity_loo = count_source_tokens(e$source),
    reporting_source_scale_diversity_loo = count_source_tokens(e$source_scale),
    reporting_geo1_share_loo = if (nrow(e)) mean(e$geo_precision == 1L,
      na.rm = TRUE) else NA_real_,
    reporting_time1_share_loo = if (nrow(e)) mean(e$time_precision == 1L,
      na.rm = TRUE) else NA_real_
  )
})

loo <- loo |>
  arrange(district_id, time_id) |>
  group_by(district_id) |>
  mutate(
    reporting_events_loo_lag1 = lag(reporting_events_loo),
    reporting_source_diversity_loo_lag1 = lag(reporting_source_diversity_loo),
    reporting_source_scale_diversity_loo_lag1 = lag(reporting_source_scale_diversity_loo),
    reporting_geo1_share_loo_lag1 = lag(reporting_geo1_share_loo),
    reporting_time1_share_loo_lag1 = lag(reporting_time1_share_loo)
  ) |>
  ungroup()

out <- hist |>
  select(
    district_id, canonical_name, year_month, time_id, cell_id,
    ends_with("lag1"), contains("prev3"), contains("prev6"), contains("rawqueen"),
    starts_with("months_since_last")
  ) |>
  left_join(
    loo |> select(district_id, year_month, matches("^reporting_.*_lag1$")),
    by = c("district_id", "year_month")
  ) |>
  arrange(district_id, time_id)

if (nrow(out) != 1632L || anyDuplicated(out$cell_id)) stop("Internal covariate output domain corrupted.")

readr::write_csv(out, out_cfg, na = "")

missing <- out |>
  summarise(across(where(is.numeric), ~ sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  arrange(desc(n_missing))
readr::write_csv(missing, file.path(qa_dir, "00_missingness.csv"), na = "")

message("Internal conflict/reporting covariates written: ", out_cfg)
message(
  "First-month lag NAs are intentional. Rolling windows require complete",
  " past history.")
print(missing, n = Inf)
