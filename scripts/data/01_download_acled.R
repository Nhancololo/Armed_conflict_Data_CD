# Download ACLED event records for Mozambique through the API.

# 0. Packages ----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  httr2,
  jsonlite,
  dplyr,
  tibble,
  readr,
  askpass,
  here,
  rstudioapi
)

# 1. Configuration -----------------------------------------------------------
acled_token_url <- "https://acleddata.com/oauth/token"
acled_api_url   <- "https://acleddata.com/api/acled/read"

start_date_cfg <- as.Date("2017-01-01")
end_date_cfg   <- Sys.Date()
country_cfg    <- "Mozambique"
province   <- "Cabo Delgado"

# Request ACLED population-exposure fields:
# population_1km, population_2km, population_5km, population_best.
population_mode <- "full"

# One event per row; actor1 and actor2 remain separate columns.
export_type_cfg <- "dyadic"

# {here} automatically anchors paths at the R Project root.
output_dir <- here::here(
  "data", "raw",
  "acled"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("\nacled download configuration")
message("Project root  : ", here::here())
message("Country       : ", country_cfg)
message("Province copy : ", province)
message("Start date    : ", start_date_cfg)
message("End date      : ", end_date_cfg)
message("Export type   : ", export_type_cfg)
message("Total row cap : NONE (cursor pagination until next_cursor = NULL)")
message("Output dir    : ", output_dir)
message("\n")

stopifnot(start_date_cfg <= end_date_cfg)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

# 2. Interactive credentials: always username -> password --------------------
get_acled_credentials <- function() {
  
  # In RStudio, use a dialog box. This avoids source() accidentally feeding the
  # next line of the script into readline(). Outside RStudio, fall back to
  # console.
  if (interactive() && rstudioapi::isAvailable()) {
    username <- rstudioapi::showPrompt(
      title   = "myACLED login",
      message = "Email / username:",
      default = ""
    )
  } else if (interactive()) {
    username <- readline("myACLED email/username: ")
  } else {
    stop(
      "This script is configured for interactive login. Run it in an ",
      "interactive R/RStudio session.")
  }
  
  username <- trimws(username)
  
  if (!nzchar(username)) {
    stop("Username/email was left blank. Download cancelled.")
  }
  
  password <- askpass::askpass("myACLED password:")
  
  if (is.null(password) || !nzchar(password)) {
    stop("Password was left blank. Download cancelled.")
  }
  
  list(username = username, password = password)
}

# 3. OAuth authentication ----------------------------------------------------
get_acled_token <- function(username, password) {
  message("Authenticating with myACLED OAuth...")
  
  response <- httr2::request(acled_token_url) |>
    httr2::req_body_form(
      username   = username,
      password   = password,
      grant_type = "password",
      client_id  = "acled",
      scope      = "authenticated"
    ) |>
    httr2::req_retry(
      max_tries = 5,
      retry_on_failure = TRUE
    ) |>
    httr2::req_timeout(seconds = 60) |>
    httr2::req_perform()
  
  token_body <- httr2::resp_body_json(response, simplifyVector = TRUE)
  
  if (is.null(token_body$access_token) || !nzchar(token_body$access_token)) {
    stop("ACLED authentication returned no access_token.")
  }
  
  message(
    "Authentication successful. Token lifetime reported by ACLED: ",
    token_body$expires_in %||% NA,
    " seconds."
  )
  
  token_body
}

credentials <- get_acled_credentials()
token_info <- get_acled_token(
  username = credentials$username,
  password = credentials$password
)
access_token <- token_info$access_token

# Remove plaintext credentials immediately after the token is obtained.
rm(credentials)

# 4. Request one API page ----------------------------------------------------
request_acled_page <- function(cursor, access_token) {
  
  # IMPORTANT:
  # There is deliberately no "limit" parameter here.
  # ACLED may still serve the response internally in pages (currently default
  # up to 5000 rows for dyadic calls). The cursor loop below follows every page
  # until next_cursor = NULL, so this does not limit the total dataset.
  params <- list(
    `_format`        = "json",
    country          = country_cfg,
    event_date       = paste0(start_date_cfg, "|", end_date_cfg),
    event_date_where = "BETWEEN",
    export_type      = export_type_cfg,
    population       = population_mode,
    with_total       = "true",
    cursor           = cursor
  )
  
  response <- httr2::request(acled_api_url) |>
    httr2::req_url_query(!!!params) |>
    httr2::req_headers(
      Authorization = paste("Bearer", access_token),
      Accept = "application/json"
    ) |>
    httr2::req_retry(
      max_tries = 6,
      retry_on_failure = TRUE
    ) |>
    httr2::req_timeout(seconds = 120) |>
    httr2::req_perform()
  
  body <- httr2::resp_body_json(response, simplifyVector = TRUE)
  
  if (!is.null(body$status) && !identical(as.integer(body$status), 200L)) {
    stop("ACLED API returned status ", body$status, ".")
  }
  
  if (!is.null(body$success) && !isTRUE(body$success)) {
    stop("ACLED API returned success = FALSE.")
  }
  
  data <- body$data
  
  if (is.null(data) || length(data) == 0L) {
    data <- tibble::tibble()
  } else {
    data <- tibble::as_tibble(data)
  }
  
  list(
    data                    = data,
    next_cursor             = body$next_cursor %||% NULL,
    count                   = body$count %||% nrow(data),
    total_count             = body$total_count %||% NA,
    last_update             = body$last_update %||% NA,
    messages                = body$messages %||% NULL,
    data_query_restrictions = body$data_query_restrictions %||% NULL,
    filename                = body$filename %||% NULL
  )
}

# 5. Download all pages with cursor pagination -------------------------------
message("\nDownloading ALL ACLED events matching the query...")

cursor <- "0"
pages <- list()
page_log <- list()
page_number <- 1L
first_metadata <- NULL
seen_cursors <- character(0)

repeat {
  message("  Request ", page_number, " | cursor = ", cursor)
  
  if (cursor %in% seen_cursors) {
    stop("Cursor loop detected before request; download aborted for safety.")
  }
  seen_cursors <- c(seen_cursors, cursor)
  
  page <- request_acled_page(
    cursor = cursor,
    access_token = access_token
  )
  
  if (is.null(first_metadata)) {
    first_metadata <- list(
      total_count = page$total_count,
      last_update = page$last_update,
      messages = page$messages,
      data_query_restrictions = page$data_query_restrictions,
      api_filename = page$filename
    )
    
    message("  ACLED total_count reported: ",
      first_metadata$total_count %||% "not supplied")
  }
  
  pages[[page_number]] <- page$data
  
  next_cursor <- page$next_cursor
  
  page_log[[page_number]] <- tibble::tibble(
    request_number = page_number,
    cursor_used = as.character(cursor),
    rows_returned = nrow(page$data),
    cumulative_rows = sum(vapply(pages, nrow, integer(1))),
    next_cursor = if (
      is.null(next_cursor) || length(next_cursor) == 0L || all(is.na(next_cursor))
    ) NA_character_ else as.character(next_cursor)[1]
  )
  
  message(
    "    rows this request = ", nrow(page$data),
    " | cumulative = ", sum(vapply(pages, nrow, integer(1)))
  )
  
  if (
    is.null(next_cursor) ||
    length(next_cursor) == 0L ||
    all(is.na(next_cursor)) ||
    !nzchar(as.character(next_cursor)[1])
  ) {
    message("  next_cursor = NULL -> final page reached.")
    break
  }
  
  cursor <- as.character(next_cursor)[1]
  page_number <- page_number + 1L
}

acled_moz_raw <- dplyr::bind_rows(pages)
page_log_df <- dplyr::bind_rows(page_log)

message("Download completed: ", format(nrow(acled_moz_raw), big.mark = ","),
  " rows.")

if (nrow(acled_moz_raw) == 0L) {
  stop(
    "The ACLED API returned zero rows. Inspect account restrictions, ",
    "authentication, and date filters.")
}

# 6. Verify completeness before saving analytical derivatives ----------------
api_total <- suppressWarnings(as.numeric(first_metadata$total_count %||% NA))

if (is.finite(api_total) && nrow(acled_moz_raw) != api_total) {
  stop(
    "INCOMPLETE DOWNLOAD: local rows = ", nrow(acled_moz_raw),
    " but ACLED total_count = ", api_total,
    ". Raw files will not be treated as complete."
  )
}

if (is.finite(api_total) && nrow(acled_moz_raw) == api_total) {
  message("Completeness check PASSED: downloaded rows == ACLED total_count == ",
    api_total)
}

# 7. Preserve RAW download exactly as returned -------------------------------
raw_rds <- file.path(
  output_dir,
  "acled_mozambique_raw.rds"
)

raw_csv <- file.path(
  output_dir,
  "acled_mozambique_raw.csv.gz"
)

saveRDS(acled_moz_raw, raw_rds, compress = "xz")
readr::write_csv(acled_moz_raw, raw_csv, na = "")

# 8. Typed derivative (RAW remains untouched) --------------------------------
acled_moz <- acled_moz_raw

if ("event_date" %in% names(acled_moz)) {
  acled_moz$event_date <- as.Date(acled_moz$event_date)
}

integer_fields <- intersect(
  c("year", "time_precision", "geo_precision", "fatalities"),
  names(acled_moz)
)

for (v in integer_fields) {
  acled_moz[[v]] <- suppressWarnings(as.integer(acled_moz[[v]]))
}

numeric_fields <- intersect(
  c(
    "latitude", "longitude", "timestamp",
    "population_1km", "population_2km", "population_5km", "population_best"
  ),
  names(acled_moz)
)

for (v in numeric_fields) {
  acled_moz[[v]] <- suppressWarnings(as.numeric(acled_moz[[v]]))
}

# Keep textual/categorical variables as character at this stage.

acled_cabo_delgado <- acled_moz |>
  dplyr::filter(.data$admin1 == province)

# 9. Quality assurance -------------------------------------------------------
id_field <- if ("event_id_cnty" %in% names(acled_moz)) "event_id_cnty" else NA_character_

n_duplicate_ids <- if (!is.na(id_field)) {
  sum(duplicated(acled_moz[[id_field]]))
} else {
  NA_integer_
}

qa_summary <- tibble::tibble(
  metric = c(
    "download_timestamp",
    "requested_start_date",
    "requested_end_date",
    "mozambique_rows",
    "cabo_delgado_rows",
    "unique_event_ids",
    "duplicate_event_ids",
    "api_total_count",
    "download_equals_api_total",
    "api_requests",
    "min_event_date",
    "max_event_date",
    "zero_fatality_events",
    "positive_fatality_events"
  ),
  value = c(
    as.character(Sys.time()),
    as.character(start_date_cfg),
    as.character(end_date_cfg),
    as.character(nrow(acled_moz)),
    as.character(nrow(acled_cabo_delgado)),
    if (!is.na(id_field)) as.character(dplyr::n_distinct(acled_moz[[id_field]])) else NA_character_,
    as.character(n_duplicate_ids),
    as.character(api_total),
    if (is.finite(api_total)) as.character(nrow(acled_moz) == api_total) else NA_character_,
    as.character(nrow(page_log_df)),
    if ("event_date" %in% names(acled_moz)) as.character(min(acled_moz$event_date, na.rm = TRUE)) else NA_character_,
    if ("event_date" %in% names(acled_moz)) as.character(max(acled_moz$event_date, na.rm = TRUE)) else NA_character_,
    if ("fatalities" %in% names(acled_moz)) as.character(sum(acled_moz$fatalities == 0, na.rm = TRUE)) else NA_character_,
    if ("fatalities" %in% names(acled_moz)) as.character(sum(acled_moz$fatalities > 0, na.rm = TRUE)) else NA_character_
  )
)

if ("country" %in% names(acled_moz)) {
  unexpected_country <- setdiff(unique(acled_moz$country), country_cfg)
  if (length(unexpected_country) > 0L) {
    warning("Unexpected country values returned: ", paste(unexpected_country,
      collapse = ", "))
  }
}

if (!is.na(n_duplicate_ids) && n_duplicate_ids > 0L) {
  warning(
    "Duplicate event_id_cnty values detected in a dyadic download: ",
    n_duplicate_ids,
    ". Investigate before aggregation."
  )
}

qa_by_year <- acled_moz |>
  dplyr::mutate(
    year_qa = if ("year" %in% names(acled_moz)) {
      .data$year
    } else {
      as.integer(format(.data$event_date, "%Y"))
    }
  ) |>
  dplyr::group_by(.data$year_qa) |>
  dplyr::summarise(
    events = dplyr::n(),
    zero_fatality_events = if ("fatalities" %in% names(acled_moz)) sum(.data$fatalities == 0, na.rm = TRUE) else NA_integer_,
    fatalities = if ("fatalities" %in% names(acled_moz)) sum(.data$fatalities,
      na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) |>
  dplyr::rename(year = year_qa)

qa_cabo_by_year <- acled_cabo_delgado |>
  dplyr::mutate(
    year_qa = if ("year" %in% names(acled_cabo_delgado)) {
      .data$year
    } else {
      as.integer(format(.data$event_date, "%Y"))
    }
  ) |>
  dplyr::group_by(.data$year_qa) |>
  dplyr::summarise(
    events = dplyr::n(),
    zero_fatality_events = if ("fatalities" %in% names(acled_cabo_delgado)) sum(.data$fatalities == 0, na.rm = TRUE) else NA_integer_,
    fatalities = if ("fatalities" %in% names(acled_cabo_delgado)) sum(.data$fatalities, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) |>
  dplyr::rename(year = year_qa)

qa_event_types <- if ("event_type" %in% names(acled_cabo_delgado)) {
  acled_cabo_delgado |>
    dplyr::count(.data$event_type, name = "n_events", sort = TRUE)
} else {
  tibble::tibble()
}

# 10. Save typed derivatives, QA, pagination log, metadata -------------------
typed_moz_rds <- file.path(
  output_dir,
  "acled_mozambique_typed.rds"
)

typed_cd_rds <- file.path(
  output_dir,
  "acled_cabo_delgado_typed.rds"
)

typed_cd_csv <- file.path(
  output_dir,
  "acled_cabo_delgado_typed.csv.gz"
)

saveRDS(acled_moz, typed_moz_rds, compress = "xz")
saveRDS(acled_cabo_delgado, typed_cd_rds, compress = "xz")
readr::write_csv(acled_cabo_delgado, typed_cd_csv, na = "")

readr::write_csv(
  page_log_df,
  file.path(output_dir, "pagination_log.csv")
)
readr::write_csv(
  qa_summary,
  file.path(output_dir, "qa_summary.csv")
)
readr::write_csv(
  qa_by_year,
  file.path(output_dir, "qa_mozambique_by_year.csv")
)
readr::write_csv(
  qa_cabo_by_year,
  file.path(output_dir, "qa_cabodelgado_by_year.csv")
)

if (nrow(qa_event_types) > 0L) {
  readr::write_csv(
    qa_event_types,
    file.path(output_dir, "qa_cabodelgado_event_types.csv")
  )
}

metadata <- list(
  downloaded_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  project_root = here::here(),
  api_endpoint = acled_api_url,
  authentication = "OAuth Bearer token; username/password prompted interactively and not stored",
  request = list(
    country = country_cfg,
    start_date = as.character(start_date_cfg),
    end_date = as.character(end_date_cfg),
    export_type = export_type_cfg,
    population = population_mode,
    explicit_limit_parameter = NULL,
    pagination = "cursor until next_cursor = NULL"
  ),
  api = first_metadata,
  local = list(
    rows_mozambique = nrow(acled_moz),
    rows_cabo_delgado = nrow(acled_cabo_delgado),
    duplicate_event_ids = n_duplicate_ids,
    requests = nrow(page_log_df),
    complete_vs_total_count = if (is.finite(api_total)) nrow(acled_moz) == api_total else NA
  )
)

jsonlite::write_json(
  metadata,
  path = file.path(output_dir, "acled_download_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null"
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "r_sessioninfo.txt")
)

# 11. Final console report ---------------------------------------------------
message("\nacled download complete")
message("Mozambique events : ", nrow(acled_moz))
message("Cabo Delgado      : ", nrow(acled_cabo_delgado))
message("ACLED total_count : ", api_total)
message("Complete          : ",
  if (is.finite(api_total)) nrow(acled_moz) == api_total else "total_count unavailable")
message("Duplicate IDs     : ", n_duplicate_ids)
message("API requests      : ", nrow(page_log_df))
message(
  "Date range in data: ",
  if ("event_date" %in% names(acled_moz)) {
    paste0(min(acled_moz$event_date, na.rm = TRUE), " to ",
      max(acled_moz$event_date, na.rm = TRUE))
  } else {
    "event_date unavailable"
  }
)
message("Output directory  : ", output_dir)
message("\n")

print(qa_summary, n = Inf)
print(qa_cabo_by_year, n = Inf)
if (nrow(qa_event_types) > 0L) print(qa_event_types, n = Inf)
