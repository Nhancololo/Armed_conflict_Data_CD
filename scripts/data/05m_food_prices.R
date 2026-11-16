# Build food price covariates.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  sf,
  dplyr,
  tidyr,
  tibble,
  readr,
  stringr,
  stringi,
  lubridate,
  curl,
  rvest,
  xml2,
  purrr,
  here
)

# 1. Paths -------------------------------------------------------------------

raw_dir <- here::here(
  "data", "raw",
  "covariates",
  "food_prices"
)

boundary_file <- here::here(
  "data", "processed",
  "cabo_delgado_admin2_canonical.rds"
)

out_file <- here::here(
  "data", "processed",
  "covariates_food_prices_2017_2024.csv"
)

qa_dir <- here::here(
  "outputs",
  "covariates",
  "05M_food_prices"
)

# Automatic World Bank RTFP acquisition --------------------------------------
#
# The World Bank/NADA catalog uses a terms page before exposing the direct
# download link. The downloader below:
#   1) opens the official "Get Microdata" page;
#   2) asks you once in the R console to confirm that YOU accept the terms;
#   3) submits the NADA "accept" form;
#   4) discovers the current MOZ_RTFP_mkt_*.csv link dynamically;
#   5) downloads it into raw_dir and continues the 05M pipeline.
#
# No login is required for this direct/open-access workflow, but the terms
# acceptance is intentionally not silently fabricated.
auto_download_rtfp <- TRUE

world_bank_catalog_url_cfg <-
  "https://microdata.worldbank.org/catalog/4500"

world_bank_get_microdata_url <-
  paste0(
    world_bank_catalog_url_cfg,
    "/get-microdata"
  )

# Terms acceptance ------------------------------------------------------------
#
# When this script is run with source(), readline() can return an empty string
# immediately in some RStudio/source contexts. To avoid that fragile console
# interaction, acceptance is read from an environment variable that survives
# rm(list = ls()).
#
# If you have read and accept the World Bank dataset terms, run ONCE in the
# R console before sourcing this script:
#
#   Sys.setenv(WB_RTFP_TERMS_ACCEPTED = "yes")
#
# Then the acquisition proceeds without a prompt.
#
# To revoke that session-level consent later:
#
#   Sys.unsetenv("WB_RTFP_TERMS_ACCEPTED")
#
# The script never prints or stores any account credential.
accept_world_bank_terms <- identical(
  tolower(
    stringr::str_squish(
      Sys.getenv(
        "WB_RTFP_TERMS_ACCEPTED",
        unset = ""
      )
    )
  ),
  "yes"
)

download_max_attempts <- 4L
download_retry_seconds <- c(2, 5, 10)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(boundary_file)) {
  stop("Missing canonical boundary: ", boundary_file)
}

# 2. Locate or automatically download the official World Bank source ---------

looks_like_rtfp_csv <- function(path) {
  
  if (
    !file.exists(path) ||
    file.info(path)$size <= 0
  ) {
    return(FALSE)
  }
  
  # zip magic bytes "PK" mean this is not a CSV even if it was accidentally
  # saved with a .csv extension.
  con <- file(
    path,
    open = "rb"
  )
  
  on.exit(
    close(con),
    add = TRUE
  )
  
  magic <- readBin(
    con,
    what = "raw",
    n = 4L
  )
  
  if (
    length(magic) >= 2L &&
    identical(
      as.integer(magic[1:2]),
      c(80L, 75L)
    )
  ) {
    return(FALSE)
  }
  
  header <- tryCatch(
    readLines(
      path,
      n = 1L,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) character()
  )
  
  if (!length(header)) {
    return(FALSE)
  }
  
  required_tokens <- c(
    "adm1_name",
    "adm2_name",
    "mkt_name",
    "year",
    "month",
    "c_food_price_index"
  )
  
  all(
    vapply(
      required_tokens,
      function(token) {
        stringr::str_detect(
          header,
          fixed(token)
        )
      },
      logical(1)
    )
  )
}

source_files_all <- list.files(
  raw_dir,
  pattern = "^MOZ_RTFP_mkt_.*\\.csv$",
  full.names = TRUE
)

source_files <- source_files_all[
  vapply(
    source_files_all,
    looks_like_rtfp_csv,
    logical(1)
  )
]

ignored_local_files <- setdiff(
  source_files_all,
  source_files
)

if (length(ignored_local_files)) {
  
  readr::write_csv(
    tibble(
      file = ignored_local_files,
      bytes = as.numeric(
        file.info(
          ignored_local_files
        )$size
      ),
      reason =
        "Filename matched MOZ_RTFP_mkt_*.csv but content is not the expected RTFP CSV"
    ),
    file.path(
      qa_dir,
      "00_ignored_invalid_local_downloads.csv"
    ),
    na = ""
  )
  
  message(
    "Ignoring ",
    length(ignored_local_files),
    " invalid local file(s) from previous failed download attempt(s)."
  )
}

auto_downloaded <- FALSE
download_url_used <- NA_character_
download_archive_used <- NA_character_
terms_action_used <- NA_character_

http_fetch_text <- function(
    url,
    handle,
    method = c("GET", "POST"),
    postfields = NULL,
    label = "HTTP request"
) {
  method <- match.arg(method)
  last_error <- NULL
  
  for (attempt in seq_len(download_max_attempts)) {
    
    result <- tryCatch(
      {
        if (method == "GET") {
          curl::handle_setopt(
            handle,
            httpget = TRUE,
            post = FALSE,
            followlocation = TRUE,
            useragent = paste(
              "R/05M Cabo-Delgado research downloader;",
              "World Bank RTFP"
            )
          )
        } else {
          curl::handle_setopt(
            handle,
            post = TRUE,
            postfields = postfields,
            followlocation = TRUE,
            useragent = paste(
              "R/05M Cabo-Delgado research downloader;",
              "World Bank RTFP"
            ),
            httpheader = c(
              "Content-Type: application/x-www-form-urlencoded"
            )
          )
        }
        
        response <- curl::curl_fetch_memory(
          url,
          handle = handle
        )
        
        if (
          response$status_code < 200L ||
          response$status_code >= 400L
        ) {
          stop(
            "HTTP status ",
            response$status_code
          )
        }
        
        rawToChar(response$content)
      },
      error = function(e) e
    )
    
    if (!inherits(result, "error")) {
      return(result)
    }
    
    last_error <- result
    
    message(
      label,
      " failed on attempt ",
      attempt,
      "/",
      download_max_attempts,
      ": ",
      conditionMessage(result)
    )
    
    if (attempt < download_max_attempts) {
      wait_i <- download_retry_seconds[
        min(
          attempt,
          length(download_retry_seconds)
        )
      ]
      
      message(
        "Retrying in ",
        wait_i,
        " seconds..."
      )
      
      Sys.sleep(wait_i)
    }
  }
  
  stop(
    label,
    " failed after ",
    download_max_attempts,
    " attempts. Last error: ",
    conditionMessage(last_error)
  )
}

confirm_world_bank_terms <- function() {
  
  if (isTRUE(accept_world_bank_terms)) {
    message(
      "World Bank terms acceptance: confirmed through ",
      "WB_RTFP_TERMS_ACCEPTED=yes."
    )
    return(TRUE)
  }
  
  stop(
    "World Bank dataset terms have not been explicitly accepted for this ",
    "R session. If you have read and accept them, run:\n\n",
    '  Sys.setenv(WB_RTFP_TERMS_ACCEPTED = "yes")\n\n',
    "and then source 05M again."
  )
}

# Extract every named field that a browser would submit with the NADA terms
# form.
# This is important because current NADA versions may include hidden CSRF
# fields;
# posting only accept=1 can therefore return the terms page again.
extract_nada_form_payload <- function(form_node) {
  
  inputs <- rvest::html_elements(
    form_node,
    "input[name]"
  )
  
  input_payload <- if (length(inputs)) {
    
    tibble(
      name = rvest::html_attr(inputs, "name"),
      type = tolower(
        dplyr::coalesce(
          rvest::html_attr(inputs, "type"),
          "text"
        )
      ),
      value = dplyr::coalesce(
        rvest::html_attr(inputs, "value"),
        ""
      ),
      checked = !is.na(
        rvest::html_attr(inputs, "checked")
      )
    ) |>
      filter(
        !is.na(name),
        nzchar(name),
        !type %in% c(
          "button",
          "submit",
          "reset",
          "file",
          "image"
        )
      ) |>
      filter(
        !type %in% c("checkbox", "radio") |
          checked |
          tolower(name) == "accept"
      ) |>
      mutate(
        value = if_else(
          tolower(name) == "accept",
          "1",
          value
        )
      ) |>
      select(
        name,
        value
      )
    
  } else {
    tibble(
      name = character(),
      value = character()
    )
  }
  
  
  selects <- rvest::html_elements(
    form_node,
    "select[name]"
  )
  
  select_payload <- purrr::map_dfr(
    selects,
    function(node) {
      
      field_name <- rvest::html_attr(
        node,
        "name"
      )
      
      options <- rvest::html_elements(
        node,
        "option"
      )
      
      if (
        is.na(field_name) ||
        !nzchar(field_name) ||
        !length(options)
      ) {
        return(
          tibble(
            name = character(),
            value = character()
          )
        )
      }
      
      selected <- options[
        !is.na(
          rvest::html_attr(
            options,
            "selected"
          )
        )
      ]
      
      if (!length(selected)) {
        selected <- options[1]
      }
      
      tibble(
        name = field_name,
        value = dplyr::coalesce(
          rvest::html_attr(
            selected,
            "value"
          ),
          rvest::html_text2(selected)
        )
      )
    }
  )
  
  
  textareas <- rvest::html_elements(
    form_node,
    "textarea[name]"
  )
  
  textarea_payload <- if (length(textareas)) {
    tibble(
      name = rvest::html_attr(
        textareas,
        "name"
      ),
      value = rvest::html_text2(
        textareas
      )
    ) |>
      filter(
        !is.na(name),
        nzchar(name)
      )
  } else {
    tibble(
      name = character(),
      value = character()
    )
  }
  
  
  payload <- bind_rows(
    input_payload,
    select_payload,
    textarea_payload
  )
  
  # NADA's direct-access controller specifically checks POST["accept"].
  if (!any(tolower(payload$name) == "accept")) {
    payload <- bind_rows(
      payload,
      tibble(
        name = "accept",
        value = "1"
      )
    )
  } else {
    payload$value[
      tolower(payload$name) == "accept"
    ] <- "1"
  }
  
  payload
}

encode_form_payload <- function(payload) {
  
  paste0(
    curl::curl_escape(
      payload$name
    ),
    "=",
    curl::curl_escape(
      payload$value
    ),
    collapse = "&"
  )
}

find_accept_form <- function(html) {
  
  doc <- rvest::read_html(html)
  
  forms <- rvest::html_elements(
    doc,
    "form"
  )
  
  if (!length(forms)) {
    return(NULL)
  }
  
  has_accept <- vapply(
    forms,
    function(form) {
      length(
        rvest::html_elements(
          form,
          '[name="accept"]'
        )
      ) > 0L
    },
    logical(1)
  )
  
  if (!any(has_accept)) {
    return(NULL)
  }
  
  forms[
    which(has_accept)[1]
  ]
}

submit_nada_accept_form <- function(
    page_html,
    page_url,
    handle
) {
  
  form <- find_accept_form(
    page_html
  )
  
  if (is.null(form)) {
    stop(
      "Could not identify the NADA direct-access acceptance form at ",
      page_url
    )
  }
  
  confirm_world_bank_terms()
  
  action <- rvest::html_attr(
    form,
    "action"
  )
  
  if (
    is.na(action) ||
    !nzchar(action)
  ) {
    action_url <- page_url
  } else {
    action_url <- xml2::url_absolute(
      action,
      page_url
    )
  }
  
  method <- toupper(
    dplyr::coalesce(
      rvest::html_attr(
        form,
        "method"
      ),
      "GET"
    )
  )
  
  payload <- extract_nada_form_payload(
    form
  )
  
  # Save field NAMES only; values can include CSRF tokens and should not be
  # written to QA output.
  readr::write_csv(
    tibble(
      submitted_field_name =
        unique(payload$name)
    ),
    file.path(
      qa_dir,
      "01_terms_form_field_names.csv"
    ),
    na = ""
  )
  
  encoded <- encode_form_payload(
    payload
  )
  
  message(
    "Submitting the official NADA terms form with ",
    nrow(payload),
    " browser-form field(s), including hidden fields."
  )
  
  if (method == "POST") {
    
    html_out <- http_fetch_text(
      action_url,
      handle = handle,
      method = "POST",
      postfields = encoded,
      label = "World Bank NADA terms acceptance"
    )
    
  } else {
    
    separator <- if (
      stringr::str_detect(
        action_url,
        fixed("?")
      )
    ) {
      "&"
    } else {
      "?"
    }
    
    html_out <- http_fetch_text(
      paste0(
        action_url,
        separator,
        encoded
      ),
      handle = handle,
      method = "GET",
      label = "World Bank NADA terms acceptance"
    )
  }
  
  list(
    html = html_out,
    action_url = action_url,
    method = method
  )
}

discover_download_link <- function(
    html,
    base_url
) {
  
  doc <- rvest::read_html(
    html
  )
  
  links <- rvest::html_elements(
    doc,
    "a[href]"
  )
  
  link_candidates <- if (length(links)) {
    
    tibble(
      text = stringr::str_squish(
        rvest::html_text2(
          links
        )
      ),
      
      raw_url = rvest::html_attr(
        links,
        "href"
      ),
      
      title = rvest::html_attr(
        links,
        "title"
      ),
      
      data_filename =
        rvest::html_attr(
          links,
          "data-filename"
        ),
      
      data_extension =
        rvest::html_attr(
          links,
          "data-extension"
        )
    )
    
  } else {
    tibble(
      text = character(),
      raw_url = character(),
      title = character(),
      data_filename = character(),
      data_extension = character()
    )
  }
  
  
  # Some NADA themes may store URLs in data-* attributes.
  attr_nodes <- rvest::html_elements(
    doc,
    "[data-url], [data-href], [data-download-url]"
  )
  
  attr_candidates <- purrr::map_dfr(
    attr_nodes,
    function(node) {
      
      values <- c(
        rvest::html_attr(
          node,
          "data-url"
        ),
        rvest::html_attr(
          node,
          "data-href"
        ),
        rvest::html_attr(
          node,
          "data-download-url"
        )
      )
      
      values <- values[
        !is.na(values) &
          nzchar(values)
      ]
      
      if (!length(values)) {
        return(
          tibble(
            text = character(),
            raw_url = character(),
            title = character(),
            data_filename = character(),
            data_extension = character()
          )
        )
      }
      
      tibble(
        text = stringr::str_squish(
          rvest::html_text2(
            node
          )
        ),
        raw_url = values,
        title = rvest::html_attr(
          node,
          "title"
        ),
        data_filename =
          rvest::html_attr(
            node,
            "data-filename"
          ),
        data_extension =
          rvest::html_attr(
            node,
            "data-extension"
          )
      )
    }
  )
  
  
  bind_rows(
    link_candidates,
    attr_candidates
  ) |>
    filter(
      !is.na(raw_url),
      nzchar(raw_url)
    ) |>
    mutate(
      absolute_url =
        xml2::url_absolute(
          raw_url,
          base_url
        ),
      
      resolved_filename =
        dplyr::coalesce(
          data_filename,
          title,
          NA_character_
        ),
      
      # This is the only family wanted by 05M.
      # Excludes MOZ_RTP_details_* and MOZ_RTP_ticker_info_*.
      is_market_rtfp_archive =
        stringr::str_detect(
          resolved_filename,
          regex(
            "^MOZ_RTFP_mkt_.*\\.zip$",
            ignore_case = TRUE
          )
        ),
      
      version_date_text =
        stringr::str_extract(
          resolved_filename,
          "\\d{4}-\\d{2}-\\d{2}(?=\\.zip$)"
        ),
      
      version_date =
        suppressWarnings(
          as.Date(
            version_date_text
          )
        ),
      
      nada_download_route =
        stringr::str_detect(
          raw_url,
          regex(
            "/catalog/4500/download/\\d+|access_direct.*/download|/download/",
            ignore_case = TRUE
          )
        )
    ) |>
    filter(
      is_market_rtfp_archive
    ) |>
    distinct(
      absolute_url,
      .keep_all = TRUE
    ) |>
    arrange(
      desc(version_date),
      desc(resolved_filename)
    )
}

download_binary_with_retry <- function(
    url,
    destfile,
    handle,
    label = "RTFP data download"
) {
  
  last_error <- NULL
  
  for (attempt in seq_len(download_max_attempts)) {
    
    result <- tryCatch(
      {
        curl::handle_setopt(
          handle,
          httpget = TRUE,
          post = FALSE,
          followlocation = TRUE,
          useragent = paste(
            "R/05M Cabo-Delgado research downloader;",
            "World Bank RTFP"
          )
        )
        
        curl::curl_download(
          url,
          destfile = destfile,
          quiet = FALSE,
          mode = "wb",
          handle = handle
        )
        
        if (
          !file.exists(destfile) ||
          file.info(destfile)$size <= 0
        ) {
          stop(
            "Downloaded file is empty or missing."
          )
        }
        
        TRUE
      },
      error = function(e) e
    )
    
    if (identical(result, TRUE)) {
      return(
        invisible(destfile)
      )
    }
    
    last_error <- result
    
    if (file.exists(destfile)) {
      unlink(destfile)
    }
    
    message(
      label,
      " failed on attempt ",
      attempt,
      "/",
      download_max_attempts,
      ": ",
      conditionMessage(result)
    )
    
    if (attempt < download_max_attempts) {
      wait_i <- download_retry_seconds[
        min(
          attempt,
          length(download_retry_seconds)
        )
      ]
      
      Sys.sleep(wait_i)
    }
  }
  
  stop(
    label,
    " failed after ",
    download_max_attempts,
    " attempts. Last error: ",
    conditionMessage(last_error)
  )
}

try_nada_direct_endpoint <- function(
    handle
) {
  
  # The open-source NADA direct-access controller is Access_direct::index($sid).
  # Try the standard CodeIgniter route variants used by NADA deployments.
  base <- "https://microdata.worldbank.org"
  
  candidates <- c(
    paste0(
      base,
      "/index.php/access_direct/index/4500"
    ),
    paste0(
      base,
      "/index.php/access_direct/4500"
    ),
    paste0(
      base,
      "/access_direct/index/4500"
    ),
    paste0(
      base,
      "/access_direct/4500"
    )
  )
  
  for (candidate in candidates) {
    
    message(
      "Trying NADA direct-access endpoint: ",
      candidate
    )
    
    page <- tryCatch(
      http_fetch_text(
        candidate,
        handle = handle,
        method = "GET",
        label = "NADA direct-access page"
      ),
      error = function(e) NULL
    )
    
    if (is.null(page)) {
      next
    }
    
    direct_links <- discover_download_link(
      page,
      candidate
    )
    
    if (nrow(direct_links)) {
      return(
        list(
          html = page,
          page_url = candidate,
          links = direct_links,
          terms_action = NA_character_
        )
      )
    }
    
    form <- find_accept_form(
      page
    )
    
    if (!is.null(form)) {
      
      accepted <- tryCatch(
        submit_nada_accept_form(
          page,
          candidate,
          handle
        ),
        error = function(e) {
          message(
            "Direct-endpoint form attempt failed: ",
            conditionMessage(e)
          )
          NULL
        }
      )
      
      if (is.null(accepted)) {
        next
      }
      
      links_after <- discover_download_link(
        accepted$html,
        accepted$action_url
      )
      
      if (nrow(links_after)) {
        return(
          list(
            html = accepted$html,
            page_url = accepted$action_url,
            links = links_after,
            terms_action = accepted$action_url
          )
        )
      }
    }
  }
  
  NULL
}

auto_download_world_bank_rtfp <- function() {
  
  message(
    "\nNo local MOZ_RTFP_mkt_*.csv found."
  )
  
  message(
    "Attempting automatic acquisition from the official World Bank ",
    "Microdata Library..."
  )
  
  h <- curl::new_handle()
  
  terms_html <- http_fetch_text(
    world_bank_get_microdata_url,
    handle = h,
    method = "GET",
    label = "World Bank Get Microdata page"
  )
  
  writeLines(
    terms_html,
    file.path(
      qa_dir,
      "00_world_bank_terms_page.html"
    ),
    useBytes = TRUE
  )
  
  
  # First inspect the wrapper page. If it already exposes data links, use them.
  links <- discover_download_link(
    terms_html,
    world_bank_get_microdata_url
  )
  
  accepted_html <- terms_html
  accepted_base_url <-
    world_bank_get_microdata_url
  form_action <- NA_character_
  
  
  if (!nrow(links)) {
    
    form <- find_accept_form(
      terms_html
    )
    
    if (!is.null(form)) {
      
      accepted <- submit_nada_accept_form(
        terms_html,
        world_bank_get_microdata_url,
        h
      )
      
      accepted_html <- accepted$html
      accepted_base_url <-
        accepted$action_url
      form_action <- accepted$action_url
      
      writeLines(
        accepted_html,
        file.path(
          qa_dir,
          "02_world_bank_download_page.html"
        ),
        useBytes = TRUE
      )
      
      links <- discover_download_link(
        accepted_html,
        accepted_base_url
      )
    }
  }
  
  
  # Fallback: talk directly to NADA's documented Access_direct controller.
  # This avoids depending on World Bank theme/wrapper markup.
  if (!nrow(links)) {
    
    direct_result <- try_nada_direct_endpoint(
      h
    )
    
    if (!is.null(direct_result)) {
      
      accepted_html <- direct_result$html
      accepted_base_url <-
        direct_result$page_url
      links <- direct_result$links
      
      if (
        is.na(form_action) &&
        !is.na(
          direct_result$terms_action
        )
      ) {
        form_action <-
          direct_result$terms_action
      }
      
      writeLines(
        accepted_html,
        file.path(
          qa_dir,
          "03_nada_direct_download_page.html"
        ),
        useBytes = TRUE
      )
    }
  }
  
  
  if (!nrow(links)) {
    stop(
      "World Bank terms were accepted, but no NADA download URL could be ",
      "resolved from either the Get Microdata page or the direct-access ",
      "controller. QA HTML pages were saved in ",
      qa_dir,
      "."
    )
  }
  
  
  # discover_download_link() already restricts the candidates to the
  # MOZ_RTFP_mkt_*.zip family and sorts newest version first.
  chosen <- links[1, ]
  
  archive_name <- chosen$resolved_filename
  
  if (
    is.na(archive_name) ||
    !stringr::str_detect(
      archive_name,
      regex(
        "^MOZ_RTFP_mkt_.*\\.zip$",
        ignore_case = TRUE
      )
    )
  ) {
    stop(
      "The selected NADA resource is not the expected market-level RTFP archive."
    )
  }
  
  archive_path <- file.path(
    raw_dir,
    archive_name
  )
  
  message(
    "Selected current market-level RTFP archive: ",
    archive_name
  )
  
  message(
    "Downloading RTFP archive from NADA: ",
    chosen$absolute_url
  )
  
  download_binary_with_retry(
    chosen$absolute_url,
    destfile = archive_path,
    handle = h
  )
  
  
  # Verify zip magic.
  con_zip <- file(
    archive_path,
    open = "rb"
  )
  
  zip_magic <- readBin(
    con_zip,
    what = "raw",
    n = 4L
  )
  
  close(con_zip)
  
  if (
    length(zip_magic) < 2L ||
    !identical(
      as.integer(
        zip_magic[1:2]
      ),
      c(80L, 75L)
    )
  ) {
    stop(
      "The selected World Bank resource is not a ZIP archive. Retained at ",
      archive_path
    )
  }
  
  
  zip_listing <- tryCatch(
    utils::unzip(
      archive_path,
      list = TRUE
    ),
    error = function(e) e
  )
  
  if (inherits(zip_listing, "error")) {
    stop(
      "Downloaded RTFP archive cannot be listed as ZIP: ",
      conditionMessage(zip_listing)
    )
  }
  
  readr::write_csv(
    as_tibble(
      zip_listing
    ),
    file.path(
      qa_dir,
      "04_rtfp_zip_contents.csv"
    ),
    na = ""
  )
  
  
  csv_inside <- zip_listing$Name[
    stringr::str_detect(
      basename(
        zip_listing$Name
      ),
      regex(
        "^MOZ_RTFP_mkt_.*\\.csv$",
        ignore_case = TRUE
      )
    )
  ]
  
  if (!length(csv_inside)) {
    stop(
      "The downloaded market-level RTFP ZIP contains no ",
      "MOZ_RTFP_mkt_*.csv file. Inspect 04_rtfp_zip_contents.csv."
    )
  }
  
  if (length(csv_inside) > 1L) {
    
    # Prefer a CSV whose date matches the selected archive version.
    archive_date <- stringr::str_extract(
      archive_name,
      "\\d{4}-\\d{2}-\\d{2}(?=\\.zip$)"
    )
    
    same_version <- csv_inside[
      stringr::str_detect(
        basename(csv_inside),
        fixed(
          archive_date
        )
      )
    ]
    
    if (length(same_version) == 1L) {
      csv_inside <- same_version
    } else {
      stop(
        "RTFP ZIP contains multiple candidate market CSVs and no unique ",
        "same-version file could be identified. Inspect 04_rtfp_zip_contents.csv."
      )
    }
  }
  
  
  extract_dir <- tempfile(
    "rtfp_extract_"
  )
  
  dir.create(
    extract_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  extracted_paths <- utils::unzip(
    archive_path,
    files = csv_inside,
    exdir = extract_dir
  )
  
  if (
    length(extracted_paths) != 1L ||
    !file.exists(
      extracted_paths
    )
  ) {
    stop(
      "Could not extract the unique RTFP market CSV from ",
      archive_name
    )
  }
  
  
  csv_name <- basename(
    extracted_paths
  )
  
  destination <- file.path(
    raw_dir,
    csv_name
  )
  
  copied <- file.copy(
    extracted_paths,
    destination,
    overwrite = TRUE
  )
  
  if (!isTRUE(copied)) {
    stop(
      "Could not copy extracted RTFP CSV to ",
      destination
    )
  }
  
  
  if (!looks_like_rtfp_csv(destination)) {
    stop(
      "Extracted file does not have the expected RTFP market CSV schema: ",
      destination
    )
  }
  
  
  message(
    "Automatic World Bank RTFP download + extraction completed."
  )
  
  message(
    "Archive: ",
    archive_path
  )
  
  message(
    "CSV:     ",
    destination
  )
  
  list(
    file = destination,
    archive = archive_path,
    url = chosen$absolute_url,
    terms_action = form_action
  )
}

if (
  !length(source_files) &&
  isTRUE(auto_download_rtfp)
) {
  
  downloaded <-
    auto_download_world_bank_rtfp()
  
  source_files <-
    downloaded$file
  
  auto_downloaded <- TRUE
  
  download_url_used <-
    downloaded$url
  
  download_archive_used <-
    downloaded$archive
  
  terms_action_used <-
    downloaded$terms_action
}

if (!length(source_files)) {
  stop(
    "No World Bank RTFP market CSV found in ",
    raw_dir,
    " and automatic download is disabled or did not produce a source file."
  )
}

# If several downloads exist, use the most recently modified copy and record it.
source_file <- source_files[
  which.max(
    file.info(source_files)$mtime
  )
]

raw <- readr::read_csv(
  source_file,
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) |>
  as_tibble()

required_cols <- c(
  "adm1_name",
  "adm2_name",
  "mkt_name",
  "lat",
  "lon",
  "year",
  "month",
  "c_food_price_index"
)

missing_cols <- setdiff(
  required_cols,
  names(raw)
)

if (length(missing_cols)) {
  stop(
    "RTFP schema is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# 3. Canonical district key --------------------------------------------------

normalize_name <- function(x) {
  x |>
    stringi::stri_trans_general("Latin-ASCII") |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", " ") |>
    stringr::str_squish()
}

districts <- readRDS(boundary_file) |>
  sf::st_make_valid() |>
  sf::st_transform(4326) |>
  arrange(district_id)

canonical <- districts |>
  sf::st_drop_geometry() |>
  transmute(
    district_id,
    canonical_name,
    normalized_name =
      normalize_name(canonical_name)
  )

aliases <- tibble::tribble(
  ~normalized_name, ~canonical_name,
  "cidade de pemba", "Pemba",
  "pemba cidade", "Pemba",
  "pemba metuge", "Metuge",
  "mocimboa da praia", "Mocímboa da Praia",
  "chiure", "Chiúre",
  "mecufi", "Mecúfi"
) |>
  left_join(
    canonical |>
      select(
        district_id,
        canonical_name
      ),
    by = "canonical_name"
  )

# 4. Source fields and conservative harmonization ----------------------------

has_inflation <-
  "inflation_food_price_index" %in%
  names(raw)

has_trust <-
  "trust_food_price_index" %in%
  names(raw)

has_spatial_interp <-
  "spatially_interpolated" %in%
  names(raw)

has_coverage <-
  "data_coverage" %in%
  names(raw)

z <- raw |>
  filter(
    stringr::str_detect(
      normalize_name(adm1_name),
      "^cabo delgado$"
    ),
    year >= 2017,
    year <= 2024
  ) |>
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    
    year_month = as.Date(
      sprintf(
        "%04d-%02d-01",
        year,
        month
      )
    ),
    
    normalized_adm2 =
      normalize_name(adm2_name),
    
    lat =
      suppressWarnings(
        as.numeric(lat)
      ),
    
    lon =
      suppressWarnings(
        as.numeric(lon)
      ),
    
    food_price_index =
      suppressWarnings(
        as.numeric(
          c_food_price_index
        )
      ),
    
    food_price_inflation =
      if (has_inflation) {
        suppressWarnings(
          as.numeric(
            inflation_food_price_index
          )
        )
      } else {
        NA_real_
      },
    
    food_price_trust =
      if (has_trust) {
        suppressWarnings(
          as.numeric(
            trust_food_price_index
          )
        )
      } else {
        NA_real_
      },
    
    food_price_spatially_interpolated =
      if (has_spatial_interp) {
        suppressWarnings(
          as.numeric(
            spatially_interpolated
          )
        )
      } else {
        NA_real_
      },
    
    food_price_data_coverage =
      if (has_coverage) {
        suppressWarnings(
          as.numeric(
            data_coverage
          )
        )
      } else {
        NA_real_
      }
  ) |>
  left_join(
    canonical,
    by = c(
      "normalized_adm2" =
        "normalized_name"
    )
  ) |>
  left_join(
    aliases |>
      select(
        normalized_name,
        district_id_alias = district_id,
        canonical_name_alias = canonical_name
      ),
    by = c(
      "normalized_adm2" =
        "normalized_name"
    )
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

if (
  any(
    !is.na(z$month) &
    !z$month %in% 1:12
  )
) {
  stop("RTFP source contains an invalid month outside 1:12.")
}

# 5. Coordinate fallback for unmatched markets -------------------------------

unmatched_idx <- which(
  is.na(z$district_id) &
    is.finite(z$lon) &
    is.finite(z$lat)
)

if (length(unmatched_idx)) {
  points <- sf::st_as_sf(
    z[unmatched_idx, ],
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  within <- sf::st_within(
    points,
    districts
  )
  
  for (k in seq_along(unmatched_idx)) {
    if (length(within[[k]]) == 1L) {
      row_id <- unmatched_idx[k]
      district_row <- within[[k]]
      
      z$district_id[row_id] <-
        districts$district_id[
          district_row
        ]
      
      z$canonical_name[row_id] <-
        districts$canonical_name[
          district_row
        ]
    }
  }
}

unresolved <- z |>
  filter(is.na(district_id)) |>
  distinct(
    adm2_name,
    mkt_name,
    lat,
    lon
  )

readr::write_csv(
  unresolved,
  file.path(
    qa_dir,
    "00_unresolved_markets.csv"
  ),
  na = ""
)

# 6. District-month aggregation ----------------------------------------------

resolved <- z |>
  filter(!is.na(district_id))

if (!nrow(resolved)) {
  stop(
    "No Cabo Delgado RTFP market rows were resolved to the 17 districts."
  )
}

# Keep the market-month structure visible in QA before district averaging.
market_month_duplicates <- resolved |>
  count(
    district_id,
    year_month,
    mkt_name,
    name = "n_rows"
  ) |>
  filter(n_rows > 1L)

readr::write_csv(
  market_month_duplicates,
  file.path(
    qa_dir,
    "01_market_month_duplicate_rows.csv"
  ),
  na = ""
)

agg <- resolved |>
  group_by(
    district_id,
    canonical_name,
    year_month
  ) |>
  summarise(
    food_price_index =
      if (all(is.na(food_price_index))) {
        NA_real_
      } else {
        mean(
          food_price_index,
          na.rm = TRUE
        )
      },
    
    food_price_inflation =
      if (all(is.na(food_price_inflation))) {
        NA_real_
      } else {
        mean(
          food_price_inflation,
          na.rm = TRUE
        )
      },
    
    food_price_trust =
      if (all(is.na(food_price_trust))) {
        NA_real_
      } else {
        mean(
          food_price_trust,
          na.rm = TRUE
        )
      },
    
    food_price_spatial_interp_share =
      if (
        all(
          is.na(
            food_price_spatially_interpolated
          )
        )
      ) {
        NA_real_
      } else {
        mean(
          food_price_spatially_interpolated,
          na.rm = TRUE
        )
      },
    
    food_price_data_coverage =
      if (all(is.na(food_price_data_coverage))) {
        NA_real_
      } else {
        mean(
          food_price_data_coverage,
          na.rm = TRUE
        )
      },
    
    # Count markets only when the primary index is actually observed.
    food_price_markets_n =
      n_distinct(
        mkt_name[
          !is.na(food_price_index)
        ]
      ),
    
    .groups = "drop"
  )

# 7. Complete 17 x 96 domain and lagged predictors ---------------------------

domain <- tidyr::crossing(
  district_id = districts$district_id,
  year_month = seq.Date(
    as.Date("2017-01-01"),
    as.Date("2024-12-01"),
    by = "month"
  )
) |>
  left_join(
    canonical |>
      select(
        district_id,
        canonical_name
      ),
    by = "district_id"
  )

out <- domain |>
  left_join(
    agg,
    by = c(
      "district_id",
      "canonical_name",
      "year_month"
    )
  ) |>
  arrange(
    district_id,
    year_month
  ) |>
  group_by(district_id) |>
  mutate(
    food_price_index_lag1 =
      lag(food_price_index),
    
    food_price_inflation_lag1 =
      lag(food_price_inflation),
    
    food_price_trust_lag1 =
      lag(food_price_trust)
  ) |>
  ungroup()

if (
  nrow(out) != 1632L ||
  anyDuplicated(
    out[
      c(
        "district_id",
        "year_month"
      )
    ]
  )
) {
  stop("Food-price output is not exactly 17 x 96.")
}

# 8. QA and provenance -------------------------------------------------------

coverage <- out |>
  group_by(
    district_id,
    canonical_name
  ) |>
  summarise(
    months_with_price =
      sum(
        !is.na(food_price_index)
      ),
    
    months_with_lag1 =
      sum(
        !is.na(food_price_index_lag1)
      ),
    
    .groups = "drop"
  )

source_provenance <- tibble(
  source_file = source_file,
  bytes = as.numeric(
    file.info(source_file)$size
  ),
  md5 = unname(
    tools::md5sum(source_file)
  ),
  auto_downloaded = auto_downloaded,
  terms_acceptance_via_env =
    accept_world_bank_terms,
  world_bank_catalog_url =
    world_bank_catalog_url_cfg,
  download_url_used =
    download_url_used,
  download_archive_used =
    download_archive_used,
  terms_action_used =
    terms_action_used,
  acquisition_timestamp_utc =
    format(
      Sys.time(),
      tz = "UTC",
      usetz = TRUE
    )
)

readr::write_csv(
  out,
  out_file,
  na = ""
)

readr::write_csv(
  coverage,
  file.path(
    qa_dir,
    "02_coverage_by_district.csv"
  ),
  na = ""
)

readr::write_csv(
  source_provenance,
  file.path(
    qa_dir,
    "03_source_file.csv"
  ),
  na = ""
)

message(
  "Food-price sensitivity covariates written: ",
  out_file
)
message(
  "No district-month price gaps were imputed. ",
  "Prefer lag-1 metrics and inspect trust/interpolation/coverage fields."
)
