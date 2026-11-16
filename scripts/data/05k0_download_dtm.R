# Download displacement records from the IOM DTM API.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  dtmapi,
  dplyr,
  tidyr,
  tibble,
  readr,
  stringr,
  stringi,
  lubridate,
  here
)

# 1. Configuration -----------------------------------------------------------

study_start <- as.Date("2017-01-01")
study_end   <- as.Date("2024-12-31")

country_name   <- "Mozambique"
admin1_target  <- "Cabo Delgado"
operation_name <- "Cabo Delgado Conflict"

force_api_key_prompt <- FALSE
refresh_operation_catalog <- FALSE

api_max_attempts <- 4L
api_retry_wait_seconds <- c(2, 5, 10)

# Strong QA: every Admin-2 date must have an Admin-1 comparator and match it.
require_admin1_for_every_admin2_date <- TRUE
require_exact_admin1_reconciliation <- TRUE

raw_dir <- here::here("data", "raw", "covariates", "iom_dtm")
qa_dir  <- here::here("outputs", "covariates", "05K0_iom_dtm_api")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir,  recursive = TRUE, showWarnings = FALSE)

operations_cache <- file.path(qa_dir, "00_operations_mozambique.csv")

raw_admin2_file <- file.path(
  raw_dir,
  "api_raw_admin2_cabo_delgado_2017_2024.csv"
)

raw_admin1_file <- file.path(
  raw_dir,
  "api_raw_admin1_cabo_delgado_2017_2024.csv"
)

manifest_file <- file.path(raw_dir, "dtm_manifest.csv")

# 2. Helpers -----------------------------------------------------------------

normalize_name <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", " ") |>
    stringr::str_squish()
}

resolve_column <- function(data, candidates, required = TRUE) {
  nms <- names(data)
  hit <- match(tolower(candidates), tolower(nms))
  hit <- hit[!is.na(hit)]
  
  if (length(hit)) {
    return(nms[hit[1]])
  }
  
  if (required) {
    stop(
      "Required API column not found. Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable: ",
      paste(nms, collapse = ", ")
    )
  }
  
  NA_character_
}

safe_character <- function(data, column) {
  if (is.na(column) || !column %in% names(data)) {
    return(rep(NA_character_, nrow(data)))
  }
  as.character(data[[column]])
}

safe_integer <- function(data, column) {
  if (is.na(column) || !column %in% names(data)) {
    return(rep(NA_integer_, nrow(data)))
  }
  suppressWarnings(as.integer(data[[column]]))
}

safe_numeric <- function(data, column) {
  if (is.na(column) || !column %in% names(data)) {
    return(rep(NA_real_, nrow(data)))
  }
  suppressWarnings(as.numeric(data[[column]]))
}

parse_api_date <- function(x) {
  as.Date(substr(as.character(x), 1, 10))
}

collapse_text <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  
  if (!length(x)) {
    return(NA_character_)
  }
  
  paste(x, collapse = " | ")
}

unique_integer_or_stop <- function(x, label, group_label) {
  x <- sort(unique(as.integer(x[!is.na(x)])))
  
  if (!length(x)) {
    return(NA_integer_)
  }
  
  if (length(x) > 1L) {
    stop(
      "Multiple ", label, " values within ", group_label, ": ",
      paste(x, collapse = ", ")
    )
  }
  
  x[1]
}

make_origin_key <- function(pcode, name) {
  pcode <- stringr::str_squish(as.character(pcode))
  name  <- stringr::str_squish(as.character(name))
  
  dplyr::case_when(
    !is.na(pcode) &
      nzchar(pcode) &
      normalize_name(pcode) != "not available" ~ paste0("PCODE:", pcode),
    
    !is.na(name) &
      nzchar(name) ~ paste0("NAME:", normalize_name(name)),
    
    TRUE ~ "ORIGIN_NOT_AVAILABLE"
  )
}

make_reason_key <- function(reason) {
  reason <- stringr::str_squish(as.character(reason))
  
  dplyr::case_when(
    is.na(reason) | !nzchar(reason) ~ "REASON_NOT_AVAILABLE",
    TRUE ~ paste0("REASON:", normalize_name(reason))
  )
}

write_schema <- function(data, path) {
  readr::write_csv(
    tibble(
      column = names(data),
      class = vapply(
        data,
        function(x) paste(class(x), collapse = "/"),
        character(1)
      )
    ),
    path,
    na = ""
  )
}

api_request <- function(label, request_fun) {
  last_error <- NULL
  
  for (attempt in seq_len(api_max_attempts)) {
    result <- tryCatch(request_fun(), error = function(e) e)
    
    if (!inherits(result, "error")) {
      if (attempt > 1L) {
        message(label, " succeeded on attempt ", attempt, ".")
      }
      return(result)
    }
    
    last_error <- result
    message(
      label, " failed on attempt ", attempt, "/", api_max_attempts,
      ": ", conditionMessage(result)
    )
    
    if (attempt < api_max_attempts) {
      wait <- api_retry_wait_seconds[
        min(attempt, length(api_retry_wait_seconds))
      ]
      message("Retrying in ", wait, " seconds...")
      Sys.sleep(wait)
    }
  }
  
  stop(
    label, " failed after ", api_max_attempts, " attempts.\n",
    "Last error: ", conditionMessage(last_error)
  )
}

# 3. Authentication ----------------------------------------------------------

key_loaded <- nzchar(Sys.getenv("DTM_SUBSCRIPTION_KEY"))

if (force_api_key_prompt || !key_loaded) {
  if (force_api_key_prompt) {
    message("Replacing the current DTM API key.")
  } else {
    message("DTM_SUBSCRIPTION_KEY is not loaded.")
  }
  
  dtmapi::set_subscription_key()
  key_loaded <- nzchar(Sys.getenv("DTM_SUBSCRIPTION_KEY"))
}

if (!key_loaded) {
  stop("DTM API subscription key is unavailable.")
}

message("DTM API key status: loaded (value not printed).")

# 4. Validate operation ------------------------------------------------------

if (file.exists(operations_cache) && !refresh_operation_catalog) {
  operations_mozambique <- readr::read_csv(
    operations_cache,
    show_col_types = FALSE
  )
} else {
  operations_all <- api_request(
    "DTM operation catalogue",
    function() dtmapi::get_all_operations() |> tibble::as_tibble()
  )
  
  op_country_col <- resolve_column(
    operations_all,
    c("admin0Name", "countryName", "country")
  )
  
  operations_mozambique <- operations_all |>
    filter(
      normalize_name(.data[[op_country_col]]) ==
        normalize_name(country_name)
    )
  
  readr::write_csv(
    operations_mozambique,
    operations_cache,
    na = ""
  )
}

op_name_col <- resolve_column(
  operations_mozambique,
  c("operation", "operationName")
)

operation_hit <- operations_mozambique |>
  filter(
    normalize_name(.data[[op_name_col]]) ==
      normalize_name(operation_name)
  )

if (nrow(operation_hit) != 1L) {
  stop(
    "Could not uniquely validate operation '", operation_name,
    "'. Inspect: ", operations_cache
  )
}

message("Using DTM operation: ", operation_name)

# 5. Download Admin-2 and Admin-1 API-v3 data --------------------------------

message("Downloading Admin-2 IDP history...")

admin2_raw <- api_request(
  "DTM Admin-2 request",
  function() {
    dtmapi::get_idp_admin2_data(
      Operation = operation_name,
      CountryName = country_name,
      Admin1Name = admin1_target,
      FromReportingDate = format(study_start, "%Y-%m-%d"),
      ToReportingDate = format(study_end, "%Y-%m-%d")
    ) |>
      tibble::as_tibble()
  }
)

message("Downloading Admin-1 IDP history for reconciliation...")

admin1_raw <- api_request(
  "DTM Admin-1 request",
  function() {
    dtmapi::get_idp_admin1_data(
      Operation = operation_name,
      CountryName = country_name,
      Admin1Name = admin1_target,
      FromReportingDate = format(study_start, "%Y-%m-%d"),
      ToReportingDate = format(study_end, "%Y-%m-%d")
    ) |>
      tibble::as_tibble()
  }
)

if (!nrow(admin2_raw)) stop("Admin-2 API returned zero rows.")
if (!nrow(admin1_raw)) stop("Admin-1 API returned zero rows.")

readr::write_csv(admin2_raw, raw_admin2_file, na = "")
readr::write_csv(admin1_raw, raw_admin1_file, na = "")

write_schema(
  admin2_raw,
  file.path(qa_dir, "01_admin2_api_schema.csv")
)

write_schema(
  admin1_raw,
  file.path(qa_dir, "02_admin1_api_schema.csv")
)

# 6. Standardize Admin-2 rows ------------------------------------------------

a2_id       <- resolve_column(admin2_raw, c("id"), FALSE)
a2_a1_name  <- resolve_column(admin2_raw, c("admin1Name", "admin1_name"))
a2_a1_code  <- resolve_column(admin2_raw, c("admin1Pcode", "admin1_pcode"),
  FALSE)
a2_a2_name  <- resolve_column(admin2_raw, c("admin2Name", "admin2_name"))
a2_a2_code  <- resolve_column(admin2_raw, c("admin2Pcode", "admin2_pcode"),
  FALSE)
a2_idp      <- resolve_column(admin2_raw, c("numPresentIdpInd",
  "num_present_idp_ind"))
a2_date     <- resolve_column(admin2_raw, c("reportingDate", "reporting_date"))
a2_round    <- resolve_column(admin2_raw, c("roundNumber", "round_number"),
  FALSE)
a2_op       <- resolve_column(admin2_raw, c("operation", "operationName"),
  FALSE)
a2_reason   <- resolve_column(admin2_raw, c("displacementReason",
  "displacement_reason"), FALSE)
a2_orig_nm  <- resolve_column(admin2_raw, c("idpOriginAdmin1Name",
  "idp_origin_admin1_name"), FALSE)
a2_orig_pc  <- resolve_column(admin2_raw, c("idpOriginAdmin1Pcode",
  "idp_origin_admin1_pcode"), FALSE)
a2_male     <- resolve_column(admin2_raw, c("numberMales", "number_males"),
  FALSE)
a2_female   <- resolve_column(admin2_raw, c("numberFemales", "number_females"),
  FALSE)
a2_type     <- resolve_column(admin2_raw, c("assessmentType",
  "assessment_type"), FALSE)

admin2 <- tibble(
  source_id = safe_integer(admin2_raw, a2_id),
  Province = safe_character(admin2_raw, a2_a1_name),
  province_pcode = safe_character(admin2_raw, a2_a1_code),
  District = safe_character(admin2_raw, a2_a2_name),
  district_pcode = safe_character(admin2_raw, a2_a2_code),
  IDPs = safe_numeric(admin2_raw, a2_idp),
  assessment_date = parse_api_date(admin2_raw[[a2_date]]),
  operation = safe_character(admin2_raw, a2_op),
  round_number = safe_integer(admin2_raw, a2_round),
  displacement_reason = safe_character(admin2_raw, a2_reason),
  origin_admin1_name = safe_character(admin2_raw, a2_orig_nm),
  origin_admin1_pcode = safe_character(admin2_raw, a2_orig_pc),
  number_males = safe_numeric(admin2_raw, a2_male),
  number_females = safe_numeric(admin2_raw, a2_female),
  assessment_type = safe_character(admin2_raw, a2_type)
) |>
  filter(
    !is.na(assessment_date),
    assessment_date >= study_start,
    assessment_date <= study_end,
    normalize_name(Province) == normalize_name(admin1_target)
  ) |>
  mutate(
    origin_key = make_origin_key(
      origin_admin1_pcode,
      origin_admin1_name
    ),
    reason_key = make_reason_key(displacement_reason)
  ) |>
  distinct()

if (!nrow(admin2)) stop("No standardized Admin-2 rows remain.")
if (any(admin2$IDPs < 0,
  na.rm = TRUE)) stop("Negative Admin-2 IDP count found.")

# 7. Admin-2 QA and atomic-stratum validation --------------------------------

# 7.1 Sex totals, where available
bad_sex <- admin2 |>
  filter(
    !is.na(IDPs),
    !is.na(number_males),
    !is.na(number_females),
    number_males + number_females != IDPs
  )

if (nrow(bad_sex)) {
  readr::write_csv(
    bad_sex,
    file.path(qa_dir, "ERROR_admin2_sex_counts.csv"),
    na = ""
  )
  stop("Admin-2 sex counts do not reproduce IDP totals.")
}

# 7.2 One measurement frame per destination district/date
frame_check <- admin2 |>
  group_by(assessment_date, district_pcode, District) |>
  summarise(
    n_rounds = n_distinct(round_number[!is.na(round_number)]),
    rounds = collapse_text(round_number),
    n_assessment_types = n_distinct(
      assessment_type[
        !is.na(assessment_type) &
          nzchar(assessment_type)
      ]
    ),
    assessment_types = collapse_text(assessment_type),
    .groups = "drop"
  )

bad_frames <- frame_check |>
  filter(n_rounds > 1L | n_assessment_types > 1L)

if (nrow(bad_frames)) {
  readr::write_csv(
    bad_frames,
    file.path(
      qa_dir,
      "ERROR_admin2_multiple_frames_same_district_date.csv"
    ),
    na = ""
  )
  stop(
    "Multiple rounds/assessment types occur for the same Admin-2 date."
  )
}

# 7.3 Atomic API-v3 record must be unique
atomic_duplicates <- admin2 |>
  count(
    assessment_date,
    district_pcode,
    District,
    origin_key,
    reason_key,
    name = "n_rows"
  ) |>
  filter(n_rows > 1L)

if (nrow(atomic_duplicates)) {
  detail <- admin2 |>
    semi_join(
      atomic_duplicates,
      by = c(
        "assessment_date",
        "district_pcode",
        "District",
        "origin_key",
        "reason_key"
      )
    ) |>
    arrange(
      assessment_date,
      District,
      origin_key,
      reason_key,
      source_id
    )
  
  readr::write_csv(
    detail,
    file.path(qa_dir, "ERROR_admin2_repeated_atomic_stratum.csv"),
    na = ""
  )
  
  stop(
    "Repeated non-identical Admin-2 origin/reason strata remain. ",
    "No summation was performed."
  )
}

# 7.4 Show cases where reason resolves the earlier origin-only duplicates
origin_reason_qa <- admin2 |>
  group_by(
    assessment_date,
    district_pcode,
    District,
    origin_key
  ) |>
  summarise(
    n_rows = n(),
    n_reason_strata = n_distinct(reason_key),
    reported_reasons = collapse_text(displacement_reason),
    idp_values = paste(sort(unique(IDPs)), collapse = " | "),
    .groups = "drop"
  ) |>
  filter(n_rows > 1L) |>
  arrange(assessment_date, District, origin_key)

readr::write_csv(
  origin_reason_qa,
  file.path(
    qa_dir,
    "03a_admin2_origin_repetitions_resolved_by_reason.csv"
  ),
  na = ""
)

readr::write_csv(
  admin2 |>
    arrange(
      assessment_date,
      District,
      origin_key,
      reason_key
    ),
  file.path(qa_dir, "03_admin2_atomic_strata.csv"),
  na = ""
)

# 8. Aggregate validated Admin-2 strata to district totals -------------------

admin2_totals <- admin2 |>
  group_by(
    assessment_date,
    Province,
    province_pcode,
    District,
    district_pcode
  ) |>
  summarise(
    IDPs = if (all(is.na(IDPs))) NA_real_ else sum(IDPs, na.rm = TRUE),
    n_origin_strata = n_distinct(origin_key),
    n_reason_strata = n_distinct(reason_key),
    n_atomic_strata = n_distinct(paste(origin_key, reason_key, sep = "||")),
    origin_admin1_names = collapse_text(origin_admin1_name),
    displacement_reasons = collapse_text(displacement_reason),
    operation = collapse_text(operation),
    round_number = unique_integer_or_stop(
      round_number,
      "round_number",
      paste(District[1], assessment_date[1])
    ),
    assessment_type = collapse_text(assessment_type),
    .groups = "drop"
  ) |>
  arrange(assessment_date, District)

if (
  anyDuplicated(
    admin2_totals[c("assessment_date", "district_pcode")]
  )
) {
  stop("Admin-2 district/date key remains duplicated after aggregation.")
}

readr::write_csv(
  admin2_totals,
  file.path(qa_dir, "04_admin2_district_totals.csv"),
  na = ""
)

# 9. Standardize Admin-1 rows using the same atomic definition ---------------

a1_name    <- resolve_column(admin1_raw, c("admin1Name", "admin1_name"))
a1_code    <- resolve_column(admin1_raw, c("admin1Pcode", "admin1_pcode"),
  FALSE)
a1_idp     <- resolve_column(admin1_raw, c("numPresentIdpInd",
  "num_present_idp_ind"))
a1_date    <- resolve_column(admin1_raw, c("reportingDate", "reporting_date"))
a1_round   <- resolve_column(admin1_raw, c("roundNumber", "round_number"),
  FALSE)
a1_reason  <- resolve_column(admin1_raw, c("displacementReason",
  "displacement_reason"), FALSE)
a1_orig_nm <- resolve_column(admin1_raw, c("idpOriginAdmin1Name",
  "idp_origin_admin1_name"), FALSE)
a1_orig_pc <- resolve_column(admin1_raw, c("idpOriginAdmin1Pcode",
  "idp_origin_admin1_pcode"), FALSE)
a1_type    <- resolve_column(admin1_raw, c("assessmentType", "assessment_type"),
  FALSE)

admin1 <- tibble(
  Province = safe_character(admin1_raw, a1_name),
  province_pcode = safe_character(admin1_raw, a1_code),
  IDPs = safe_numeric(admin1_raw, a1_idp),
  assessment_date = parse_api_date(admin1_raw[[a1_date]]),
  round_number = safe_integer(admin1_raw, a1_round),
  displacement_reason = safe_character(admin1_raw, a1_reason),
  origin_admin1_name = safe_character(admin1_raw, a1_orig_nm),
  origin_admin1_pcode = safe_character(admin1_raw, a1_orig_pc),
  assessment_type = safe_character(admin1_raw, a1_type)
) |>
  filter(
    !is.na(assessment_date),
    assessment_date >= study_start,
    assessment_date <= study_end,
    normalize_name(Province) == normalize_name(admin1_target)
  ) |>
  mutate(
    origin_key = make_origin_key(
      origin_admin1_pcode,
      origin_admin1_name
    ),
    reason_key = make_reason_key(displacement_reason)
  ) |>
  distinct()

if (!nrow(admin1)) stop("No standardized Admin-1 rows remain.")

# One measurement frame per Admin-1 date
a1_frames <- admin1 |>
  group_by(assessment_date, province_pcode, Province) |>
  summarise(
    n_rounds = n_distinct(round_number[!is.na(round_number)]),
    n_assessment_types = n_distinct(
      assessment_type[
        !is.na(assessment_type) &
          nzchar(assessment_type)
      ]
    ),
    rounds = collapse_text(round_number),
    assessment_types = collapse_text(assessment_type),
    .groups = "drop"
  )

bad_a1_frames <- a1_frames |>
  filter(n_rounds > 1L | n_assessment_types > 1L)

if (nrow(bad_a1_frames)) {
  readr::write_csv(
    bad_a1_frames,
    file.path(qa_dir, "ERROR_admin1_multiple_frames_same_date.csv"),
    na = ""
  )
  stop("Multiple rounds/assessment types occur for one Admin-1 date.")
}

# Atomic Admin-1 key must also be unique
a1_duplicates <- admin1 |>
  count(
    assessment_date,
    province_pcode,
    origin_key,
    reason_key,
    name = "n_rows"
  ) |>
  filter(n_rows > 1L)

if (nrow(a1_duplicates)) {
  detail <- admin1 |>
    semi_join(
      a1_duplicates,
      by = c(
        "assessment_date",
        "province_pcode",
        "origin_key",
        "reason_key"
      )
    ) |>
    arrange(
      assessment_date,
      origin_key,
      reason_key
    )
  
  readr::write_csv(
    detail,
    file.path(qa_dir, "ERROR_admin1_repeated_atomic_stratum.csv"),
    na = ""
  )
  
  stop("Repeated non-identical Admin-1 origin/reason strata remain.")
}

readr::write_csv(
  admin1 |>
    arrange(
      assessment_date,
      origin_key,
      reason_key
    ),
  file.path(qa_dir, "05_admin1_atomic_strata.csv"),
  na = ""
)

admin1_totals <- admin1 |>
  group_by(
    assessment_date,
    Province,
    province_pcode
  ) |>
  summarise(
    admin1_IDPs = if (all(is.na(IDPs))) NA_real_ else sum(IDPs, na.rm = TRUE),
    admin1_n_origin_strata = n_distinct(origin_key),
    admin1_n_reason_strata = n_distinct(reason_key),
    admin1_n_atomic_strata = n_distinct(
      paste(origin_key, reason_key, sep = "||")
    ),
    admin1_round = unique_integer_or_stop(
      round_number,
      "round_number",
      paste(Province[1], assessment_date[1])
    ),
    .groups = "drop"
  ) |>
  arrange(assessment_date)

readr::write_csv(
  admin1_totals,
  file.path(qa_dir, "06_admin1_totals.csv"),
  na = ""
)

# 10. Reconcile Admin-2 and Admin-1 totals -----------------------------------

admin2_by_date <- admin2_totals |>
  group_by(assessment_date) |>
  summarise(
    admin2_IDPs = if (all(is.na(IDPs))) NA_real_ else sum(IDPs, na.rm = TRUE),
    n_admin2_districts = sum(!is.na(IDPs)),
    n_admin2_atomic_strata = sum(n_atomic_strata),
    .groups = "drop"
  )

reconciliation <- full_join(
  admin1_totals |>
    select(
      assessment_date,
      admin1_IDPs,
      admin1_n_atomic_strata,
      admin1_round
    ),
  admin2_by_date,
  by = "assessment_date"
) |>
  mutate(
    difference = admin2_IDPs - admin1_IDPs,
    exact_match = case_when(
      is.na(admin1_IDPs) | is.na(admin2_IDPs) ~ NA,
      TRUE ~ difference == 0
    ),
    coverage_status = case_when(
      !is.na(admin1_IDPs) & !is.na(admin2_IDPs) ~ "both",
      is.na(admin1_IDPs) & !is.na(admin2_IDPs) ~ "admin2_only",
      !is.na(admin1_IDPs) & is.na(admin2_IDPs) ~ "admin1_only",
      TRUE ~ "neither"
    )
  ) |>
  arrange(assessment_date)

readr::write_csv(
  reconciliation,
  file.path(qa_dir, "07_admin1_admin2_reconciliation.csv"),
  na = ""
)

if (
  require_admin1_for_every_admin2_date &&
  any(reconciliation$coverage_status == "admin2_only")
) {
  readr::write_csv(
    reconciliation |>
      filter(coverage_status == "admin2_only"),
    file.path(qa_dir, "ERROR_admin2_date_without_admin1.csv"),
    na = ""
  )
  
  stop(
    "At least one Admin-2 date has no Admin-1 comparator. ",
    "No final manifest was written."
  )
}

bad_reconciliation <- reconciliation |>
  filter(
    coverage_status == "both",
    !exact_match
  )

if (
  require_exact_admin1_reconciliation &&
  nrow(bad_reconciliation)
) {
  readr::write_csv(
    bad_reconciliation,
    file.path(qa_dir, "ERROR_admin1_admin2_totals_do_not_match.csv"),
    na = ""
  )
  
  stop(
    "Admin-2 district totals do not exactly reproduce the Admin-1 ",
    "Cabo Delgado total for some dates. No final manifest was written."
  )
}

# 11. Write one district-total CSV per Admin-2 assessment date ---------------

assessment_summary <- admin2_totals |>
  group_by(assessment_date) |>
  summarise(
    n_districts = n_distinct(District),
    total_idps = if (all(is.na(IDPs))) NA_real_ else sum(IDPs, na.rm = TRUE),
    n_atomic_strata = sum(n_atomic_strata),
    round_number = unique_integer_or_stop(
      round_number,
      "round_number",
      paste("assessment", assessment_date[1])
    ),
    .groups = "drop"
  ) |>
  arrange(assessment_date)

readr::write_csv(
  assessment_summary,
  file.path(qa_dir, "08_assessment_dates.csv"),
  na = ""
)

assessment_dates <- assessment_summary$assessment_date
written_files <- character(length(assessment_dates))

for (i in seq_along(assessment_dates)) {
  date_i <- assessment_dates[i]
  
  out_i <- admin2_totals |>
    filter(assessment_date == date_i) |>
    transmute(
      Province,
      District,
      IDPs,
      Returnees = NA_real_,
      district_pcode,
      n_origin_strata,
      n_reason_strata,
      n_atomic_strata,
      origin_admin1_names,
      displacement_reasons,
      operation,
      round_number,
      assessment_type
    ) |>
    arrange(District)
  
  file_i <- paste0(
    "dtm_api_",
    format(date_i, "%Y%m%d"),
    ".csv"
  )
  
  readr::write_csv(
    out_i,
    file.path(raw_dir, file_i),
    na = ""
  )
  
  written_files[i] <- file_i
}

# 12. Build dtm_manifest.csv for 05K -----------------------------------------

manifest <- assessment_summary |>
  transmute(
    assessment_date,
    file = paste0(
      "dtm_api_",
      format(assessment_date, "%Y%m%d"),
      ".csv"
    ),
    url = NA_character_,
    sheet = NA_character_,
    province_col = "Province",
    district_col = "District",
    idp_col = "IDPs",
    returnee_col = NA_character_,
    aggregation_mode = "district_total",
    notes = paste0(
      "Official IOM-DTM API v3; operation=",
      operation_name,
      "; district totals reconstructed from unique origin x ",
      "displacement-reason strata and reconciled against Admin-1; round=",
      ifelse(is.na(round_number), "not supplied", round_number),
      "; returnees unavailable from this acquisition route"
    )
  )

if (anyDuplicated(manifest$assessment_date)) {
  stop("Generated manifest contains duplicate assessment dates.")
}

readr::write_csv(
  manifest,
  manifest_file,
  na = ""
)

# 13. Final report -----------------------------------------------------------

if (!all(file.exists(file.path(raw_dir, written_files)))) {
  stop("At least one assessment CSV was not written.")
}

message("\n05k0 complete")
message("Operation                  : ", operation_name)
message("Admin-2 raw rows           : ", nrow(admin2_raw))
message("Admin-2 atomic rows        : ", nrow(admin2))
message("District-date totals       : ", nrow(admin2_totals))
message("Assessment dates           : ", length(assessment_dates))
message("Earliest assessment        : ", min(assessment_dates))
message("Latest assessment          : ", max(assessment_dates))
message(
  "Admin1/Admin2 common dates : ",
  sum(reconciliation$coverage_status == "both")
)
message(
  "Exact reconciliation       : ",
  if (
    all(
      reconciliation$exact_match[
        reconciliation$coverage_status == "both"
      ]
    )
  ) "PASS" else "CHECK QA"
)
message("Manifest                   : ", manifest_file)
message("Returnees                   : unavailable -> NA")
message("\nNext:")
message(
  'source(here::here("scripts", ',
  '"05k_ingest_dtm.R"))'
)
message("\n")
