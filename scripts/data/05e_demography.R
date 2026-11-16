# Build demographic covariates from the Census Bureau IDB.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  readxl,
  dplyr,
  tidyr,
  tibble,
  readr,
  stringr,
  stringi,
  lubridate,
  sf,
  here
)

# 1. Configuration -----------------------------------------------------------
years_cfg <- 2017:2024
growth_base_year <- min(years_cfg) - 1L
idb_parent_code <- "MZ002"       # Cabo Delgado ADM1 in the supplied IDB export
idb_district_level <- 2L
distance_crs <- 32737L            # UTM 37S

out_cfg <- here::here(
  "data", "processed",
  "covariates_idb_demography_2017_2024.csv"
)

qa_dir <- here::here(
  "outputs",
  "covariates",
  "05E_idb"
)

canonical_table_file <- here::here(
  "data", "processed",
  "cabo_delgado_admin2_canonical.csv"
)

canonical_sf_file <- here::here(
  "data", "processed",
  "cabo_delgado_admin2_canonical.rds"
)

dir.create(dirname(out_cfg), recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(
  here::here("data", "raw", "covariates", "idb"),
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Helpers -----------------------------------------------------------------
norm <- function(x) {
  x |>
    as.character() |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", " ") |>
    stringr::str_squish()
}

write_qa <- function(x, filename) {
  readr::write_csv(
    x,
    file.path(qa_dir, filename),
    na = ""
  )
  invisible(x)
}

first_existing <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  ok <- paths[file.exists(paths)]
  if (!length(ok)) return(NA_character_)
  normalizePath(ok[1], winslash = "/", mustWork = TRUE)
}

# 3. Resolve official manual IDB workbook ------------------------------------
env_file <- Sys.getenv("IDB_MOZAMBIQUE_XLSX", unset = "")

manual_candidates <- c(
  env_file,
  here::here("data", "raw", "covariates", "idb", "Mozambique.xlsx"),
  here::here("Mozambique.xlsx"),
  here::here("data", "raw", "Mozambique.xlsx")
)

manual_file <- first_existing(manual_candidates)

if (is.na(manual_file)) {
  stop(
    "Official IDB Mozambique Excel workbook was not found.\n",
    "Place the manually downloaded file at:\n  ",
    here::here("data", "raw", "covariates", "idb", "Mozambique.xlsx"),
    "\n(or set environment variable IDB_MOZAMBIQUE_XLSX to its exact path).\n",
    "The v5 pipeline intentionally does not depend on the fragile Census API."
  )
}

message("IDB manual workbook: ", manual_file)

sheets <- readxl::excel_sheets(manual_file)

if (!"Final" %in% sheets) {
  stop(
    "IDB workbook does not contain sheet 'Final'. Available sheets: ",
    paste(sheets, collapse = ", ")
  )
}

# 4. Read and validate workbook schema ---------------------------------------
idb_all <- readxl::read_excel(
  manual_file,
  sheet = "Final",
  .name_repair = "unique"
) |>
  tibble::as_tibble()

required_columns <- c(
  "CODE", "BASENAME", "ADM_LEVEL", "YR", "GENC",
  "SEXRATIO", "DEPND", "DEPND0_14", "DEPND65_", "POP",
  "POP15_19", "POP20_24", "POP25_29",
  "MPOP15_19", "MPOP20_24", "MPOP25_29"
)

missing_columns <- setdiff(required_columns, names(idb_all))

if (length(missing_columns) > 0L) {
  stop(
    "Manual IDB workbook is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

idb_all <- idb_all |>
  mutate(
    CODE = as.character(.data$CODE),
    BASENAME = as.character(.data$BASENAME),
    GENC = as.character(.data$GENC),
    ADM_LEVEL = suppressWarnings(as.integer(.data$ADM_LEVEL)),
    YR = suppressWarnings(as.integer(.data$YR))
  )

# Verify that the expected parent exists in the workbook.
parent_check <- idb_all |>
  filter(
    .data$CODE == idb_parent_code,
    .data$ADM_LEVEL == 1L
  ) |>
  distinct(.data$CODE, .data$BASENAME, .data$ADM_LEVEL)

if (
  nrow(parent_check) != 1L ||
  norm(parent_check$BASENAME[1]) != "cabo delgado"
) {
  write_qa(parent_check, "ERROR_parent_code_check.csv")
  stop(
    "Could not verify IDB parent code MZ002 as Cabo Delgado ADM1. ",
    "Inspect ERROR_parent_code_check.csv."
  )
}

# 5. Canonical 17-district domain --------------------------------------------
if (!file.exists(canonical_table_file)) {
  stop("Canonical district table not found: ", canonical_table_file)
}

if (!file.exists(canonical_sf_file)) {
  stop("Canonical district geometry not found: ", canonical_sf_file)
}

canonical <- readr::read_csv(
  canonical_table_file,
  show_col_types = FALSE
) |>
  select(district_id, canonical_name) |>
  distinct() |>
  mutate(norm_name = norm(.data$canonical_name))

if (
  nrow(canonical) != 17L ||
  n_distinct(canonical$district_id) != 17L ||
  n_distinct(canonical$canonical_name) != 17L
) {
  stop("Canonical district table is not exactly 17 unique districts.")
}

# Source-specific name exception. Do not make this a global alias because
# "Pemba" has different meanings in some legacy administrative sources.
idb_alias <- tibble::tribble(
  ~source_name,       ~canonical_name_alias,
  "Cidade de Pemba", "Pemba"
) |>
  mutate(
    source_norm = norm(.data$source_name),
    alias_norm = norm(.data$canonical_name_alias)
  )

# 6. Extract Cabo Delgado ADM2 rows, including 2016 for growth calculation ---
needed_years <- growth_base_year:max(years_cfg)

raw <- idb_all |>
  filter(
    .data$GENC == "MZ",
    .data$ADM_LEVEL == idb_district_level,
    stringr::str_starts(.data$CODE, idb_parent_code),
    .data$YR %in% needed_years
  ) |>
  transmute(
    source_code = .data$CODE,
    source_name = .data$BASENAME,
    year = .data$YR,
    sex_ratio = suppressWarnings(as.numeric(.data$SEXRATIO)),
    dependency_ratio = suppressWarnings(as.numeric(.data$DEPND)),
    youth_dependency_ratio = suppressWarnings(as.numeric(.data$DEPND0_14)),
    old_age_dependency_ratio = suppressWarnings(as.numeric(.data$DEPND65_)),
    population = suppressWarnings(as.numeric(.data$POP)),
    pop15_19 = suppressWarnings(as.numeric(.data$POP15_19)),
    pop20_24 = suppressWarnings(as.numeric(.data$POP20_24)),
    pop25_29 = suppressWarnings(as.numeric(.data$POP25_29)),
    mpop15_19 = suppressWarnings(as.numeric(.data$MPOP15_19)),
    mpop20_24 = suppressWarnings(as.numeric(.data$MPOP20_24)),
    mpop25_29 = suppressWarnings(as.numeric(.data$MPOP25_29)),
    source_norm = norm(.data$BASENAME)
  )

# Workbook structure expected from the supplied official export:
# 17 ADM2 districts x 9 years (2016-2024) = 153 rows.
expected_raw_rows <- 17L * length(needed_years)

if (
  nrow(raw) != expected_raw_rows ||
  n_distinct(raw$source_code) != 17L ||
  n_distinct(raw$year) != length(needed_years)
) {
  write_qa(raw, "ERROR_manual_IDB_domain.csv")
  stop(
    "Manual IDB ADM2 extraction is not exactly 17 districts x ",
    length(needed_years),
    " years (2016-2024). Inspect ERROR_manual_IDB_domain.csv."
  )
}

# Each source code must map to a single name and each district-year must be
# unique.
if (anyDuplicated(raw[c("source_code", "year")]) > 0L) {
  dup <- raw |>
    count(.data$source_code, .data$year, name = "n") |>
    filter(.data$n > 1L)
  write_qa(dup, "ERROR_duplicate_IDB_district_year.csv")
  stop("Duplicate IDB district-year rows detected.")
}

# 7. Exact source-specific district reconciliation ---------------------------
raw <- raw |>
  left_join(
    idb_alias |>
      select(source_norm, alias_norm),
    by = "source_norm"
  ) |>
  mutate(
    target_norm = dplyr::coalesce(.data$alias_norm, .data$source_norm)
  ) |>
  left_join(
    canonical |>
      select(
        district_id,
        canonical_name,
        target_norm = norm_name
      ),
    by = "target_norm"
  )

crosswalk <- raw |>
  distinct(
    .data$source_code,
    .data$source_name,
    .data$source_norm,
    .data$target_norm,
    .data$district_id,
    .data$canonical_name
  ) |>
  arrange(.data$source_code)

write_qa(crosswalk, "01_idb_manual_district_crosswalk.csv")

if (
  any(is.na(crosswalk$district_id)) ||
  n_distinct(crosswalk$district_id) != 17L ||
  nrow(crosswalk) != 17L
) {
  stop(
    "IDB manual district names did not reconcile one-to-one to all 17 canonical ",
    "districts. Inspect 01_idb_manual_district_crosswalk.csv."
  )
}

# Explicit QA that the workbook contains the 17 source labels we expect.
write_qa(
  raw |>
    distinct(.data$source_code, .data$source_name) |>
    arrange(.data$source_code),
  "00_idb_manual_source_districts.csv"
)

# 8. District areas from frozen canonical geometry ---------------------------
district_sf <- readRDS(canonical_sf_file)

if (!inherits(district_sf, "sf")) {
  stop("Canonical district geometry is not an sf object.")
}

# Compute area outside dplyr::mutate(). Using st_area(.) inside mutate.sf is
# fragile because the magrittr/base-pipe placeholder is not available inside
# dplyr's data mask. Work explicitly with the projected sf object instead.
if (!"district_id" %in% names(district_sf)) {
  stop("Canonical district geometry is missing district_id.")
}

district_sf_metric <- district_sf |>
  sf::st_make_valid() |>
  sf::st_transform(distance_crs)

if (
  nrow(district_sf_metric) != 17L ||
  dplyr::n_distinct(district_sf_metric$district_id) != 17L
) {
  stop("Canonical district geometry is not exactly 17 unique district IDs.")
}

district_area_km2 <- as.numeric(
  sf::st_area(district_sf_metric)
) / 1e6

area_table <- tibble::tibble(
  district_id = as.character(district_sf_metric$district_id),
  district_area_km2 = district_area_km2
) |>
  arrange(.data$district_id)

if (
  nrow(area_table) != 17L ||
  any(!is.finite(area_table$district_area_km2)) ||
  any(area_table$district_area_km2 <= 0)
) {
  stop("Canonical district areas are invalid.")
}

if (!setequal(area_table$district_id, canonical$district_id)) {
  write_qa(
    dplyr::full_join(
      canonical |> select(district_id, canonical_name),
      area_table,
      by = "district_id"
    ),
    "ERROR_canonical_geometry_ID_mismatch.csv"
  )
  stop(
    "Canonical geometry district IDs do not match the canonical district table. ",
    "Inspect ERROR_canonical_geometry_ID_mismatch.csv."
  )
}

# 9. Construct annual demographic covariates ---------------------------------
required_numeric <- c(
  "population",
  "sex_ratio",
  "dependency_ratio",
  "youth_dependency_ratio",
  "old_age_dependency_ratio",
  "pop15_19", "pop20_24", "pop25_29",
  "mpop15_19", "mpop20_24", "mpop25_29"
)

bad_numeric <- vapply(
  raw[required_numeric],
  function(x) any(!is.finite(x)),
  logical(1)
)

if (any(bad_numeric)) {
  stop(
    "Non-finite values found in required IDB numeric fields: ",
    paste(names(bad_numeric)[bad_numeric], collapse = ", ")
  )
}

if (any(raw$population <= 0)) {
  stop("IDB population contains non-positive values.")
}

annual_with_base <- raw |>
  select(
    district_id,
    canonical_name,
    source_code,
    source_name,
    year,
    all_of(required_numeric)
  ) |>
  left_join(area_table, by = "district_id") |>
  arrange(.data$district_id, .data$year) |>
  group_by(.data$district_id) |>
  mutate(
    # Annual percent change, using 2016 only as the lag for 2017.
    population_growth_rate =
      100 * (.data$population / dplyr::lag(.data$population) - 1),
    
    population_density =
      .data$population / .data$district_area_km2,
    
    youth_share_15_24 =
      100 * (.data$pop15_19 + .data$pop20_24) / .data$population,
    
    youth_share_15_29 =
      100 * (
        .data$pop15_19 +
          .data$pop20_24 +
          .data$pop25_29
      ) / .data$population,
    
    male_youth_share_15_29 =
      100 * (
        .data$mpop15_19 +
          .data$mpop20_24 +
          .data$mpop25_29
      ) / .data$population
  ) |>
  ungroup()

annual <- annual_with_base |>
  filter(.data$year %in% years_cfg) |>
  transmute(
    .data$district_id,
    .data$canonical_name,
    .data$year,
    .data$population,
    .data$population_density,
    .data$population_growth_rate,
    .data$sex_ratio,
    .data$dependency_ratio,
    .data$youth_dependency_ratio,
    .data$old_age_dependency_ratio,
    .data$youth_share_15_24,
    .data$youth_share_15_29,
    .data$male_youth_share_15_29
  ) |>
  arrange(.data$district_id, .data$year)

if (nrow(annual) != 17L * length(years_cfg)) {
  write_qa(annual, "ERROR_annual_domain.csv")
  stop("IDB annual data are not exactly 17 districts x 8 years.")
}

if (
  n_distinct(annual$district_id) != 17L ||
  n_distinct(annual$year) != 8L
) {
  stop("IDB annual district/year domain failed final QA.")
}

# 2017 growth must be available because 2016 was loaded as the baseline.
if (any(!is.finite(annual$population_growth_rate))) {
  stop(
    "population_growth_rate contains missing/non-finite values. ",
    "The manual workbook must include 2016-2024 for all 17 districts."
  )
}

# 10. Expand annual estimates to the complete district-month panel -----------
monthly_domain <- tidyr::crossing(
  district_id = canonical$district_id,
  year_month = seq.Date(
    as.Date("2017-01-01"),
    as.Date("2024-12-01"),
    by = "month"
  )
) |>
  mutate(year = lubridate::year(.data$year_month)) |>
  left_join(
    canonical |>
      select(district_id, canonical_name),
    by = "district_id"
  )

out <- monthly_domain |>
  left_join(
    annual,
    by = c("district_id", "canonical_name", "year")
  ) |>
  arrange(.data$district_id, .data$year_month)

if (nrow(out) != 1632L) {
  stop("IDB monthly expansion is not exactly 17 x 96 = 1,632 rows.")
}

covariate_cols <- setdiff(
  names(annual),
  c("district_id", "canonical_name", "year")
)

if (
  any(
    vapply(
      out[covariate_cols],
      function(x) any(!is.finite(x)),
      logical(1)
    )
  )
) {
  stop("At least one final IDB covariate contains non-finite values.")
}

# 11. Save outputs and provenance --------------------------------------------
readr::write_csv(out, out_cfg, na = "")
write_qa(area_table, "02_canonical_district_areas_km2.csv")
write_qa(annual, "03_idb_annual_final.csv")

source_info <- tibble(
  field = c(
    "source",
    "manual_file",
    "manual_file_md5",
    "sheet",
    "parent_code",
    "parent_name",
    "district_adm_level",
    "input_years",
    "analysis_years",
    "districts",
    "annual_rows",
    "monthly_rows",
    "district_area_method",
    "population_density_definition",
    "population_growth_definition",
    "monthly_expansion"
  ),
  value = c(
    "U.S. Census Bureau International Data Base manual Excel export",
    manual_file,
    unname(tools::md5sum(manual_file)),
    "Final",
    idb_parent_code,
    "Cabo Delgado",
    as.character(idb_district_level),
    paste(range(needed_years), collapse = "-"),
    paste(range(years_cfg), collapse = "-"),
    "17",
    as.character(nrow(annual)),
    as.character(nrow(out)),
    paste0("st_area after st_make_valid + st_transform(EPSG:", distance_crs,
      ")"),
    "IDB POP divided by frozen canonical district polygon area (km^2)",
    "100 * (POP_t / POP_{t-1} - 1); 2016 used only as baseline for 2017",
    "Annual values repeated within calendar year; no monthly interpolation"
  )
)

write_qa(source_info, "04_idb_source_and_definitions.csv")

variable_definitions <- tibble::tribble(
  ~variable, ~manual_IDB_field_or_formula, ~unit_or_scale,
  "population", "POP", "persons",
  "population_density", "POP / canonical district area", "persons per km2",
  "population_growth_rate", "100 * (POP_t / POP_t-1 - 1)", "percent per year",
  "sex_ratio", "SEXRATIO", "as supplied by IDB export",
  "dependency_ratio", "DEPND", "IDB ratio scale",
  "youth_dependency_ratio", "DEPND0_14", "IDB ratio scale",
  "old_age_dependency_ratio", "DEPND65_", "IDB ratio scale",
  "youth_share_15_24", "100*(POP15_19+POP20_24)/POP", "percent of total population",
  "youth_share_15_29", "100*(POP15_19+POP20_24+POP25_29)/POP", "percent of total population",
  "male_youth_share_15_29", "100*(MPOP15_19+MPOP20_24+MPOP25_29)/POP", "percent of total population"
)

write_qa(variable_definitions, "05_idb_variable_definitions.csv")

message("\nidb manual-excel build complete")
message("Source workbook      : ", manual_file)
message("Cabo Delgado code    : ", idb_parent_code)
message("District level       : ADM_LEVEL = ", idb_district_level)
message("Canonical districts  : ", n_distinct(annual$district_id))
message("Annual rows          : ", nrow(annual), " (17 x 8)")
message("Monthly rows         : ", nrow(out), " (17 x 96)")
message("Output               : ", out_cfg)
message("QA directory         : ", qa_dir)
message("\n")

print(
  crosswalk |>
    select(
      source_code,
      source_name,
      district_id,
      canonical_name
    ),
  n = Inf,
  width = Inf
)

message(
  "\nIMPORTANT: Annual IDB estimates are repeated within calendar year; ",
  "no monthly interpolation was performed."
)
