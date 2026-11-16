# Build force deployment and mining covariates.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  dplyr,
  tidyr,
  tibble,
  readr,
  lubridate,
  here
)

# 1. Files and frozen panel --------------------------------------------------

panel_file <- here::here(
  "data", "processed",
  "cabo_delgado_district_month_panel_2017_2024.rds"
)

out_file <- here::here(
  "data", "processed",
  "covariates_interventions_financing_2017_2024.csv"
)

qa_dir <- here::here(
  "outputs",
  "covariates",
  "05J_interventions_financing"
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(panel_file)) {
  stop("Missing frozen panel: ", panel_file)
}

panel <- readRDS(panel_file) |>
  as_tibble() |>
  select(
    district_id,
    canonical_name,
    year_month,
    time_id,
    cell_id
  ) |>
  arrange(district_id, time_id)

if (
  nrow(panel) != 1632L ||
  dplyr::n_distinct(panel$district_id) != 17L ||
  dplyr::n_distinct(panel$time_id) != 96L ||
  anyDuplicated(panel$cell_id)
) {
  stop("Frozen panel failed the expected 17 x 96 key checks.")
}

# 2. Prespecified chronology -------------------------------------------------

# Rwanda Security Forces (RSF)
rsf_start <- as.Date("2021-07-09")
rsf_end_study <- as.Date("2024-12-31")  # study censoring, not withdrawal

# Official retrospective descriptions identify Palma and Mocimboa da Praia as
# core RSF areas of responsibility since the 2021 deployment.
rsf_core_aor <- c("CD_009", "CD_015")

# Initial Rwanda National Police unit destinations reported on 9 July 2021:
# Mweda/Mueda and Afungi (Palma). This is a sector-designation sensitivity,
# not proof of continuous police presence through the whole study period.
rnp_initial_sectors <- c("CD_011", "CD_015")

# District-specific expansion / joint-operation milestones.
rsf_macomia_milestone <- as.Date("2022-03-30")
rsf_ancuabe_milestone <- as.Date("2022-12-21")

# Approximate RSF force-size milestones. These are sensitivity-only.
rsf_2500_milestone <- as.Date("2022-12-01")
rsf_4000_milestone <- as.Date("2024-05-01")

# SAMIM
samim_start <- as.Date("2021-07-15")
samim_drawdown_start <- as.Date("2023-12-15")
samim_end <- as.Date("2024-07-15")

# SADC documented these as SAMIM operational districts in Nov 2022. This is
# evidence of presence by that date, not proof that Nov 2022 was the start date.
samim_documented_by_2022_11 <- c(
  "CD_011", # Mueda
  "CD_005", # Macomia
  "CD_014", # Nangade
  "CD_016"  # Pemba
)
samim_documentation_date <- as.Date("2022-11-14")

# Private military / security-provider chronology (sensitivity only).
dag_start_approx <- as.Date("2020-04-01")
dag_end <- as.Date("2021-04-06")
wagner_start_approx <- as.Date("2019-09-01")
wagner_end_approx <- as.Date("2019-11-30")
paramount_burnham_milestone <- as.Date("2021-02-01")

# Islamic State organizational / external-support milestones.
# These are not measured financing amounts.
external_support_milestone <- as.Date("2020-09-01")
al_karrar_milestone <- as.Date("2021-07-01")
moz_province_milestone <- as.Date("2022-05-01")

# 3. Helpers -----------------------------------------------------------------

month_fraction <- function(month_start, exposure_start, exposure_end) {
  month_start <- as.Date(month_start)
  month_end <- lubridate::ceiling_date(month_start, "month") -
    lubridate::days(1)

  lo <- pmax(month_start, exposure_start)
  hi <- pmin(month_end, exposure_end)

  overlap_days <- pmax(
    0,
    as.numeric(hi - lo) + 1
  )

  month_days <- as.numeric(month_end - month_start) + 1

  overlap_days / month_days
}

active_in_month <- function(month_start, exposure_start, exposure_end) {
  month_end <- lubridate::ceiling_date(month_start, "month") -
    lubridate::days(1)

  month_end >= exposure_start &
    month_start <= exposure_end
}

# 4. District-month covariates -----------------------------------------------

out <- panel |>
  mutate(
    # --------------------------------------------------------------------------
    # Rwanda / RSF
    # --------------------------------------------------------------------------
    rsf_active = active_in_month(
      year_month,
      rsf_start,
      rsf_end_study
    ),

    rsf_exposure_fraction = month_fraction(
      year_month,
      rsf_start,
      rsf_end_study
    ),

    # District-specific core area-of-responsibility indicator.
    rsf_core_aor_palma_mocimboa =
      rsf_active &
      district_id %in% rsf_core_aor,

    rnp_initial_sector_mueda_afungi =
      rsf_active &
      district_id %in% rnp_initial_sectors,

    # Legacy registry-compatible name. Interpret only as an initial sector
    # designation sensitivity, not as observed continuous troop presence.
    rdf_initial_sector =
      rnp_initial_sector_mueda_afungi,

    # Milestones: evidence of deployment / operations from these dates onward.
    # They are not interpreted as continuously observed troop counts.
    rsf_macomia_post_2022_03 =
      district_id == "CD_005" &
      year_month >= floor_date(rsf_macomia_milestone, "month"),

    rsf_ancuabe_post_2022_12 =
      district_id == "CD_001" &
      year_month >= floor_date(rsf_ancuabe_milestone, "month"),

    # Approximate province-wide force-size history; sensitivity only.
    rsf_force_size_approx = case_when(
      !rsf_active ~ 0,
      year_month >= rsf_4000_milestone ~ 4000,
      year_month >= rsf_2500_milestone ~ 2500,
      TRUE ~ 1000
    ),

    # Registry-compatible aliases retained so the existing eligibility matrix
    # does not silently lose previously declared candidates.
    rdf_active = rsf_active,
    rdf_exposure_fraction = rsf_exposure_fraction,
    rdf_troop_strength_approx = rsf_force_size_approx,
    rdf_ancuabe_post_2022_12 = rsf_ancuabe_post_2022_12,

    # --------------------------------------------------------------------------
    # SAMIM
    # --------------------------------------------------------------------------
    samim_active = active_in_month(
      year_month,
      samim_start,
      samim_end
    ),

    samim_exposure_fraction = month_fraction(
      year_month,
      samim_start,
      samim_end
    ),

    samim_drawdown =
      active_in_month(
        year_month,
        samim_drawdown_start,
        samim_end
      ),

    samim_documented_area_by_2022_11 =
      samim_active &
      district_id %in% samim_documented_by_2022_11 &
      year_month >= floor_date(samim_documentation_date, "month"),

    # --------------------------------------------------------------------------
    # Other security providers
    # --------------------------------------------------------------------------
    dag_active_approx = active_in_month(
      year_month,
      dag_start_approx,
      dag_end
    ),

    dag_exposure_fraction_approx = month_fraction(
      year_month,
      dag_start_approx,
      dag_end
    ),

    wagner_active_approx = active_in_month(
      year_month,
      wagner_start_approx,
      wagner_end_approx
    ),

    dag_active = dag_active_approx,
    wagner_active = wagner_active_approx,

    post_paramount_burnham_2021_02 =
      year_month >= paramount_burnham_milestone,

    # --------------------------------------------------------------------------
    # External-support / organizational milestones
    # --------------------------------------------------------------------------
    post_external_support_2020_09 =
      year_month >= external_support_milestone,

    post_al_karrar_2021_07 =
      year_month >= al_karrar_milestone,

    post_independent_moz_province_2022_05 =
      year_month >= moz_province_milestone,

    financing_history_stage = case_when(
      year_month < external_support_milestone ~
        "pre_documented_external_support_milestone",

      year_month < al_karrar_milestone ~
        "documented_external_support",

      year_month < moz_province_milestone ~
        "al_karrar_regional_coordination",

      TRUE ~
        "independent_mozambique_province"
    )
  ) |>
  arrange(district_id, time_id)

if (
  nrow(out) != 1632L ||
  anyDuplicated(out$cell_id)
) {
  stop("Intervention output corrupted the 17 x 96 panel key.")
}

# 5. Reproducibility metadata ------------------------------------------------

source_registry <- tibble::tribble(
  ~construct, ~date_or_period, ~precision, ~interpretation, ~source,
  "RSF deployment", "2021-07-09 onward", "day",
  "Province-wide deployment calendar",
  "Rwanda Ministry of Defence, deployment announcement, 2021-07-09/10",
  "RSF core AOR", "from 2021 deployment", "retrospective official description",
  "Palma and Mocimboa da Praia core areas of responsibility",
  "Rwanda Ministry of Defence, security-review statements",
  "RSF Macomia milestone", "2022-03-30", "day",
  "Joint operations/deployment milestone; not continuous-presence proof",
  "Rwanda Ministry of Defence, 2022-05-01 report",
  "RSF Ancuabe milestone", "2022-12-21", "day",
  "New deployment / operations milestone",
  "Rwanda Ministry of Defence, 2022-12-23 and 2023-02-05 reports",
  "RSF force size", "1000 -> ~2500 -> ~4000", "approximate milestone",
  "Sensitivity-only force-size history",
  "Rwanda official source for ~2500 (2022-12); ACLED/Cabo Ligado for ~4000 (2024-05)",
  "SAMIM deployment", "2021-07-15 to 2024-07-15", "day",
  "Official mission calendar",
  "SADC official mission summaries",
  "SAMIM drawdown", "2023-12-15 to 2024-07-15", "day",
  "Official phased drawdown window",
  "SADC official 2024 mission update",
  "SAMIM documented operational districts", "documented by 2022-11", "documentation date",
  "Mueda, Macomia, Nangade, Pemba; not claimed as start date",
  "SADC Executive Secretary visit, 2022-11",
  "DAG", "2020-04 to 2021-04-06", "approximate start / day end",
  "Private military contractor period; sensitivity only",
  "Public reporting / human-rights documentation",
  "Wagner", "2019-09 to 2019-11", "approximate month",
  "Short private-military episode; sensitivity only",
  "Public reporting; treat cautiously",
  "External-support milestone", "2020-09", "month",
  "Historical milestone, not financing amount",
  "UN/CTC reporting cited in the study literature",
  "Al-Karrar milestone", "2021-07", "month",
  "Organizational milestone, not financing amount",
  "Weiss et al. / CTC literature",
  "Independent IS Mozambique province", "2022-05", "month",
  "Organizational-status milestone",
  "Weiss et al. / CTC literature"
)

readr::write_csv(
  out,
  out_file,
  na = ""
)

readr::write_csv(
  source_registry,
  file.path(qa_dir, "00_intervention_source_registry.csv"),
  na = ""
)

readr::write_csv(
  out |>
    distinct(
      year_month,
      rsf_active,
      rsf_force_size_approx,
      samim_active,
      samim_drawdown,
      dag_active_approx,
      wagner_active_approx,
      post_external_support_2020_09,
      post_al_karrar_2021_07,
      post_independent_moz_province_2022_05,
      financing_history_stage
    ) |>
    arrange(year_month),
  file.path(qa_dir, "01_monthly_timeline.csv"),
  na = ""
)

readr::write_csv(
  out |>
    group_by(district_id, canonical_name) |>
    summarise(
      months_rsf_core_aor =
        sum(rsf_core_aor_palma_mocimboa),
      months_rsf_macomia_post =
        sum(rsf_macomia_post_2022_03),
      months_rsf_ancuabe_post =
        sum(rsf_ancuabe_post_2022_12),
      months_samim_documented_area =
        sum(samim_documented_area_by_2022_11),
      .groups = "drop"
    ),
  file.path(qa_dir, "02_district_specific_intervention_months.csv"),
  na = ""
)

message(
  "Intervention / organizational-history covariates written: ",
  out_file
)
message(
  "Financing and deployment milestones are sensitivity variables; ",
  "they are not measured capital flows or clean causal treatments."
)
