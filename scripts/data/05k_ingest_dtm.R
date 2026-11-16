# Build displacement covariates from the IOM DTM records.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  sf,
  dplyr,
  tidyr,
  tibble,
  readr,
  readxl,
  stringr,
  stringi,
  lubridate,
  curl,
  here
)

# 1. Paths -------------------------------------------------------------------

boundary_file <- here::here(
  "data", "processed",
  "cabo_delgado_admin2_canonical.rds"
)

raw_dir <- here::here(
  "data", "raw",
  "covariates",
  "iom_dtm"
)

manifest_file <- file.path(
  raw_dir,
  "dtm_manifest.csv"
)

out_file <- here::here(
  "data", "processed",
  "covariates_iom_dtm_2017_2024.csv"
)

qa_dir <- here::here(
  "outputs",
  "covariates",
  "05K_iom_dtm"
)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(boundary_file)) {
  stop("Missing canonical boundary: ", boundary_file)
}

# 2. Canonical 17-district key -----------------------------------------------

district_sf <- readRDS(boundary_file) |>
  sf::st_make_valid() |>
  arrange(district_id)

canonical <- district_sf |>
  sf::st_drop_geometry() |>
  select(district_id, canonical_name)

if (
  nrow(canonical) != 17L ||
  anyDuplicated(canonical$district_id)
) {
  stop("Canonical district key is not exactly 17 unique districts.")
}

normalize_name <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", " ") |>
    stringr::str_squish()
}

canonical <- canonical |>
  mutate(normalized_name = normalize_name(canonical_name))

# Conservative source-specific aliases only.
aliases <- tibble::tribble(
  ~normalized_name, ~canonical_name,
  "chiure", "Chiúre",
  "mecufi", "Mecúfi",
  "mocimboa da praia", "Mocímboa da Praia",
  "cidade de pemba", "Pemba",
  "pemba cidade", "Pemba",
  "pemba metuge", "Metuge",
  "pemba metugi", "Metuge"
) |>
  left_join(
    canonical |>
      select(district_id, canonical_name),
    by = "canonical_name"
  )

# 3. Manifest contract -------------------------------------------------------

if (!file.exists(manifest_file)) {
  template <- tibble(
    assessment_date = c(
      "2023-04-30",
      "2024-06-30"
    ),
    file = c(
      "REPLACE_round_18.xlsx",
      "REPLACE_round_XX.xlsx"
    ),
    url = c("", ""),
    sheet = c("1", "1"),
    province_col = c("Province", "Province"),
    district_col = c("District", "District"),
    idp_col = c("IDPs", "IDPs"),
    returnee_col = c("Returnees", "Returnees"),
    aggregation_mode = c(
      "district_total",
      "sum_rows"
    ),
    notes = c(
      paste(
        "Use district_total only when there is exactly one district total row",
        "per district after filtering Cabo Delgado."
      ),
      paste(
        "Use sum_rows only when rows are non-overlapping lower-level/site",
        "records that should legitimately be summed to district level."
      )
    )
  )

  readr::write_csv(
    template,
    manifest_file,
    na = ""
  )

  stop(
    "Created IOM-DTM manifest template: ",
    manifest_file,
    "\nReview every column, replace example rows, declare aggregation_mode, ",
    "then rerun. Do not infer the schema from the examples."
  )
}

manifest <- readr::read_csv(
  manifest_file,
  show_col_types = FALSE,
  na = c("", "NA")
) |>
  mutate(
    assessment_date = as.Date(assessment_date),
    aggregation_mode = stringr::str_to_lower(
      stringr::str_squish(aggregation_mode)
    )
  ) |>
  filter(
    !is.na(assessment_date),
    !is.na(file),
    nzchar(file)
  ) |>
  arrange(assessment_date)

required_manifest_cols <- c(
  "assessment_date",
  "file",
  "url",
  "sheet",
  "province_col",
  "district_col",
  "idp_col",
  "returnee_col",
  "aggregation_mode",
  "notes"
)

missing_manifest_cols <- setdiff(
  required_manifest_cols,
  names(manifest)
)

if (length(missing_manifest_cols)) {
  stop(
    "DTM manifest is missing required columns: ",
    paste(missing_manifest_cols, collapse = ", ")
  )
}

if (!nrow(manifest)) {
  stop("DTM manifest contains no valid assessment rows.")
}

if (anyDuplicated(manifest$assessment_date)) {
  stop(
    "Manifest has duplicate assessment_date values. ",
    "Use one consolidated source specification per assessment date."
  )
}

allowed_modes <- c(
  "district_total",
  "sum_rows"
)

bad_modes <- setdiff(
  unique(manifest$aggregation_mode),
  allowed_modes
)

if (length(bad_modes)) {
  stop(
    "Unsupported aggregation_mode: ",
    paste(bad_modes, collapse = ", "),
    ". Allowed values: ",
    paste(allowed_modes, collapse = ", ")
  )
}

# 4. Source helpers ----------------------------------------------------------

download_if_needed <- function(local_file, url) {
  if (file.exists(local_file)) {
    return(invisible(local_file))
  }

  if (is.na(url) || !nzchar(url)) {
    stop(
      "Missing DTM source file and no URL supplied: ",
      local_file
    )
  }

  message("Downloading IOM-DTM source: ", url)

  last_error <- NULL

  for (attempt in 1:4) {
    result <- tryCatch(
      {
        curl::curl_download(
          url,
          local_file,
          mode = "wb",
          quiet = FALSE
        )
        NULL
      },
      error = function(e) e
    )

    if (is.null(result) && file.exists(local_file)) {
      return(invisible(local_file))
    }

    last_error <- result

    if (file.exists(local_file)) {
      unlink(local_file)
    }

    if (attempt < 4) {
      Sys.sleep(2^(attempt - 1))
    }
  }

  stop(
    "Failed to download DTM source after 4 attempts: ",
    url,
    "\nLast error: ",
    conditionMessage(last_error)
  )
}

read_tabular_source <- function(path, sheet) {
  ext <- tolower(tools::file_ext(path))

  if (ext == "csv") {
    return(
      readr::read_csv(
        path,
        show_col_types = FALSE
      ) |>
        as_tibble()
    )
  }

  if (ext %in% c("xlsx", "xls")) {
    sheet_arg <- suppressWarnings(
      as.integer(sheet)
    )

    if (is.na(sheet_arg)) {
      sheet_arg <- sheet
    }

    return(
      readxl::read_excel(
        path,
        sheet = sheet_arg
      ) |>
        as_tibble()
    )
  }

  stop(
    "Unsupported IOM-DTM file extension: ",
    ext
  )
}

parse_count <- function(x) {
  readr::parse_number(
    as.character(x),
    na = c("", "NA", "N/A", "-", "..")
  )
}

# 5. Harmonize one assessment ------------------------------------------------

read_assessment <- function(row) {
  local_file <- file.path(
    raw_dir,
    row$file
  )

  download_if_needed(
    local_file,
    row$url
  )

  x <- read_tabular_source(
    local_file,
    row$sheet
  )

  required_source_cols <- c(
    row$district_col,
    row$idp_col
  )

  if (
    !is.na(row$province_col) &&
    nzchar(row$province_col)
  ) {
    required_source_cols <- c(
      required_source_cols,
      row$province_col
    )
  }

  if (
    !is.na(row$returnee_col) &&
    nzchar(row$returnee_col)
  ) {
    required_source_cols <- c(
      required_source_cols,
      row$returnee_col
    )
  }

  missing_source_cols <- setdiff(
    required_source_cols,
    names(x)
  )

  if (length(missing_source_cols)) {
    stop(
      "Manifest/source mismatch for ",
      row$file,
      ". Missing columns: ",
      paste(missing_source_cols, collapse = ", ")
    )
  }

  if (
    !is.na(row$province_col) &&
    nzchar(row$province_col)
  ) {
    x <- x |>
      filter(
        stringr::str_detect(
          normalize_name(
            .data[[row$province_col]]
          ),
          "^cabo delgado$"
        )
      )
  }

  returnee <- if (
    !is.na(row$returnee_col) &&
    nzchar(row$returnee_col)
  ) {
    parse_count(
      x[[row$returnee_col]]
    )
  } else {
    rep(NA_real_, nrow(x))
  }

  z <- x |>
    transmute(
      source_district =
        as.character(.data[[row$district_col]]),

      normalized_name =
        normalize_name(source_district),

      idp_stock =
        parse_count(.data[[row$idp_col]]),

      returnee_stock =
        returnee
    ) |>
    left_join(
      canonical |>
        select(
          district_id,
          canonical_name,
          normalized_name
        ),
      by = "normalized_name"
    ) |>
    left_join(
      aliases |>
        select(
          normalized_name,
          district_id_alias = district_id,
          canonical_name_alias = canonical_name
        ),
      by = "normalized_name"
    ) |>
    mutate(
      district_id =
        coalesce(
          district_id,
          district_id_alias
        ),

      canonical_name =
        coalesce(
          canonical_name,
          canonical_name_alias
        )
    ) |>
    select(
      -district_id_alias,
      -canonical_name_alias
    )

  unresolved <- z |>
    filter(is.na(district_id)) |>
    distinct(
      source_district,
      normalized_name
    )

  if (nrow(unresolved)) {
    qa_file <- file.path(
      qa_dir,
      paste0(
        "UNRESOLVED_",
        format(
          row$assessment_date,
          "%Y%m%d"
        ),
        ".csv"
      )
    )

    readr::write_csv(
      unresolved,
      qa_file,
      na = ""
    )

    stop(
      "Unresolved IOM-DTM district names for assessment ",
      row$assessment_date,
      ". Inspect ",
      qa_file,
      " and add only validated explicit aliases."
    )
  }

  if (
    any(
      z$idp_stock < 0,
      na.rm = TRUE
    ) ||
    any(
      z$returnee_stock < 0,
      na.rm = TRUE
    )
  ) {
    stop(
      "Negative displacement stock found in ",
      row$file,
      "."
    )
  }

  if (row$aggregation_mode == "district_total") {
    duplicates <- z |>
      count(
        district_id,
        name = "n_rows"
      ) |>
      filter(n_rows > 1L)

    if (nrow(duplicates)) {
      qa_file <- file.path(
        qa_dir,
        paste0(
          "ERROR_duplicate_district_totals_",
          format(
            row$assessment_date,
            "%Y%m%d"
          ),
          ".csv"
        )
      )

      readr::write_csv(
        duplicates,
        qa_file,
        na = ""
      )

      stop(
        "aggregation_mode='district_total' but multiple rows exist ",
        "for at least one district. Inspect ",
        qa_file,
        "."
      )
    }

    aggregated <- z |>
      select(
        district_id,
        canonical_name,
        idp_stock,
        returnee_stock
      )
  } else {
    # sum_rows is allowed only because the manifest explicitly declares that
    # the source rows are non-overlapping lower-level/site records.
    aggregated <- z |>
      group_by(
        district_id,
        canonical_name
      ) |>
      summarise(
        idp_stock =
          if (all(is.na(idp_stock))) {
            NA_real_
          } else {
            sum(idp_stock, na.rm = TRUE)
          },

        returnee_stock =
          if (all(is.na(returnee_stock))) {
            NA_real_
          } else {
            sum(returnee_stock, na.rm = TRUE)
          },

        .groups = "drop"
      )
  }

  canonical |>
    select(
      district_id,
      canonical_name
    ) |>
    left_join(
      aggregated,
      by = c(
        "district_id",
        "canonical_name"
      )
    ) |>
    mutate(
      assessment_date = row$assessment_date,
      source_file = row$file,
      aggregation_mode =
        row$aggregation_mode
    )
}

# 6. Read all declared assessments -------------------------------------------

rounds <- bind_rows(
  lapply(
    seq_len(nrow(manifest)),
    function(i) {
      read_assessment(
        manifest[i, ]
      )
    }
  )
) |>
  arrange(
    assessment_date,
    district_id
  )

if (
  anyDuplicated(
    rounds[
      c(
        "district_id",
        "assessment_date"
      )
    ]
  )
) {
  stop(
    "Harmonized IOM-DTM assessments contain duplicate district/date keys."
  )
}

readr::write_csv(
  rounds,
  file.path(
    qa_dir,
    "00_all_assessments_harmonized.csv"
  ),
  na = ""
)

# 7. Causal-order monthly carry-forward --------------------------------------

months <- seq.Date(
  as.Date("2017-01-01"),
  as.Date("2024-12-01"),
  by = "month"
)

round_dates <- sort(
  unique(rounds$assessment_date)
)

monthly <- tidyr::crossing(
  district_id = canonical$district_id,
  year_month = months
) |>
  left_join(
    canonical |>
      select(
        district_id,
        canonical_name
      ),
    by = "district_id"
  ) |>
  rowwise() |>
  mutate(
    dtm_assessment_date_used = {
      candidate_dates <-
        round_dates[
          round_dates < year_month
        ]

      if (length(candidate_dates)) {
        max(candidate_dates)
      } else {
        as.Date(NA)
      }
    }
  ) |>
  ungroup() |>
  left_join(
    rounds |>
      select(
        district_id,
        dtm_assessment_date_used =
          assessment_date,
        idp_stock,
        returnee_stock,
        source_file,
        aggregation_mode
      ),
    by = c(
      "district_id",
      "dtm_assessment_date_used"
    )
  ) |>
  mutate(
    months_since_dtm_assessment =
      if_else(
        !is.na(
          dtm_assessment_date_used
        ),
        interval(
          dtm_assessment_date_used,
          year_month
        ) %/% months(1),
        NA_integer_
      ),

    dtm_stale_gt_6m =
      months_since_dtm_assessment > 6L
  ) |>
  arrange(
    district_id,
    year_month
  )

if (
  nrow(monthly) != 1632L ||
  anyDuplicated(
    monthly[
      c(
        "district_id",
        "year_month"
      )
    ]
  )
) {
  stop("IOM-DTM monthly output is not exactly 17 x 96.")
}

# 8. QA and provenance -------------------------------------------------------

source_provenance <- purrr::map_dfr(
  seq_len(nrow(manifest)),
  function(i) {
    path <- file.path(
      raw_dir,
      manifest$file[i]
    )

    exists <- file.exists(path)

    tibble(
      assessment_date =
        manifest$assessment_date[i],
      file =
        manifest$file[i],
      aggregation_mode =
        manifest$aggregation_mode[i],
      notes =
        manifest$notes[i],
      exists = exists,
      bytes =
        if (exists) {
          as.numeric(file.info(path)$size)
        } else {
          NA_real_
        },
      md5 =
        if (exists) {
          unname(tools::md5sum(path))
        } else {
          NA_character_
        }
    )
  }
)

coverage_by_assessment <- rounds |>
  group_by(assessment_date) |>
  summarise(
    districts_with_idp =
      sum(!is.na(idp_stock)),
    districts_with_returnees =
      sum(!is.na(returnee_stock)),
    .groups = "drop"
  )

staleness <- monthly |>
  count(
    dtm_assessment_date_used,
    months_since_dtm_assessment,
    sort = TRUE
  )

readr::write_csv(
  monthly,
  out_file,
  na = ""
)

readr::write_csv(
  source_provenance,
  file.path(
    qa_dir,
    "01_source_provenance.csv"
  ),
  na = ""
)

readr::write_csv(
  coverage_by_assessment,
  file.path(
    qa_dir,
    "02_coverage_by_assessment.csv"
  ),
  na = ""
)

readr::write_csv(
  staleness,
  file.path(
    qa_dir,
    "03_assessment_staleness_distribution.csv"
  ),
  na = ""
)

message(
  "IOM-DTM sensitivity covariates written: ",
  out_file
)
message(
  "No interpolation was used. Month t receives only the latest assessment ",
  "strictly before month t, and staleness is retained explicitly."
)
