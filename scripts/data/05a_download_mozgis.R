# Download vector layers from MozGIS.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  httr2,
  jsonlite,
  sf,
  dplyr,
  purrr,
  tibble,
  readr,
  here,
  fs
)

# 1. Configuration -----------------------------------------------------------
raw_dir <- here::here("data", "raw", "covariates", "mozgis")
qa_dir  <- here::here("outputs", "covariates", "05A_mozgis")

boundary_file <- here::here(
  "data", "processed",
  "cabo_delgado_admin2_canonical.rds"
)

buffer_m <- 50000
default_chunk_size <- 100L
max_retries <- 4L
fail_if_any_layer_fails <- TRUE

for (d in c(raw_dir, qa_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(boundary_file)) {
  stop("Missing canonical boundary: ", boundary_file)
}

districts <- readRDS(boundary_file)

if (!inherits(districts, "sf")) {
  stop("Canonical boundary is not sf.")
}

if (is.na(sf::st_crs(districts))) {
  stop("Canonical boundary has no CRS.")
}

# Canonical geometry is used only to define the download/cropping extent.
districts4326 <- districts |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

study_buffer4326 <- districts |>
  sf::st_make_valid() |>
  sf::st_transform(32737) |>
  sf::st_union() |>
  sf::st_buffer(buffer_m) |>
  sf::st_transform(4326)

study_bbox <- sf::st_bbox(study_buffer4326)
study_bbox_sfc <- sf::st_as_sfc(study_bbox)

bbox_esri_json <- jsonlite::toJSON(
  list(
    xmin = unname(study_bbox[["xmin"]]),
    ymin = unname(study_bbox[["ymin"]]),
    xmax = unname(study_bbox[["xmax"]]),
    ymax = unname(study_bbox[["ymax"]]),
    spatialReference = list(wkid = 4326)
  ),
  auto_unbox = TRUE,
  digits = 16
)

# 2. HTTP helpers ------------------------------------------------------------
# Local null-coalescing helper to avoid another package dependency.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

perform_with_retry <- function(req, max_tries = max_retries) {
  last_error <- NULL

  for (attempt in seq_len(max_tries)) {
    ans <- try(
      httr2::req_perform(req),
      silent = TRUE
    )

    if (!inherits(ans, "try-error")) {
      return(ans)
    }

    last_error <- ans

    if (attempt < max_tries) {
      Sys.sleep(min(2^(attempt - 1), 8))
    }
  }

  stop(
    "HTTP request failed after ",
    max_tries,
    " attempts: ",
    as.character(last_error)
  )
}

parse_json_response <- function(resp, url) {
  # MozGIS/ArcGIS can return valid JSON while declaring Content-Type:
  # text/plain.
  # Decode the raw payload explicitly instead of MIME-gating the body parser.
  raw <- httr2::resp_body_raw(resp)

  if (length(raw) == 0L) {
    stop("Empty response body from ", url)
  }

  txt <- rawToChar(raw)
  txt <- sub("^\ufeff", "", txt)
  txt <- trimws(txt)

  out <- tryCatch(
    jsonlite::fromJSON(txt, simplifyVector = TRUE),
    error = function(e) {
      content_type <- tryCatch(
        httr2::resp_header(resp, "content-type"),
        error = function(...) NA_character_
      )

      preview <- substr(
        gsub("[\r\n\t]+", " ", txt),
        1L,
        500L
      )

      stop(
        "Could not parse ArcGIS response as JSON from ",
        url,
        ". HTTP status=",
        httr2::resp_status(resp),
        "; content-type=",
        content_type %||% NA_character_,
        "; body preview='",
        preview,
        "'. Parser error: ",
        conditionMessage(e)
      )
    }
  )

  if (!is.null(out$error)) {
    code <- out$error$code %||% NA
    message_txt <- out$error$message %||% "Unknown ArcGIS REST error"
    details <- out$error$details %||% character()

    details_txt <- if (length(details) > 0L) {
      paste(as.character(details), collapse = " | ")
    } else {
      ""
    }

    stop(
      "ArcGIS REST error from ",
      url,
      ": code=",
      code,
      "; message=",
      message_txt,
      if (nzchar(details_txt)) paste0("; details=", details_txt) else ""
    )
  }

  out
}

safe_get_json <- function(url, query = list()) {
  req <- httr2::request(url) |>
    httr2::req_timeout(120)

  if (length(query) > 0L) {
    req <- do.call(
      httr2::req_url_query,
      c(list(req), query)
    )
  }

  resp <- perform_with_retry(req)
  parse_json_response(resp, url)
}

safe_post <- function(url, body) {
  req <- httr2::request(url) |>
    httr2::req_timeout(180)

  req <- do.call(
    httr2::req_body_form,
    c(list(req), body)
  ) |>
    httr2::req_method("POST")

  perform_with_retry(req)
}

safe_post_json <- function(url, body) {
  resp <- safe_post(url, body)
  parse_json_response(resp, url)
}

# 3. ArcGIS REST helpers -----------------------------------------------------
arcgis_layer_metadata <- function(layer_url) {
  m <- safe_get_json(
    layer_url,
    list(f = "pjson")
  )

  if (is.null(m$type)) {
    stop("Layer metadata did not return a layer type: ", layer_url)
  }

  if (!identical(as.character(m$type), "Feature Layer")) {
    stop(
      "Expected Feature Layer but service reports '",
      as.character(m$type),
      "': ",
      layer_url
    )
  }

  fields <- m$fields
  oid_field <- NA_character_

  if (!is.null(fields) && NROW(fields) > 0L) {
    if (is.data.frame(fields)) {
      hit <- which(fields$type == "esriFieldTypeOID")
      if (length(hit) > 0L) oid_field <- as.character(fields$name[hit[1]])
    } else if (is.list(fields)) {
      hit <- which(vapply(
        fields,
        function(z) identical(z$type, "esriFieldTypeOID"),
        logical(1)
      ))
      if (length(hit) > 0L) oid_field <- as.character(fields[[hit[1]]]$name)
    }
  }

  formats <- paste(
    as.character(m$supportedQueryFormats %||% ""),
    collapse = ","
  )

  use_format <- if (
    grepl("geojson", formats, ignore.case = TRUE)
  ) {
    "geojson"
  } else {
    "json"
  }

  list(
    name = as.character(m$name %||% NA_character_),
    type = as.character(m$type %||% NA_character_),
    geometry_type = as.character(m$geometryType %||% NA_character_),
    max_record_count = as.integer(m$maxRecordCount %||% default_chunk_size),
    supported_query_formats = formats,
    query_format = use_format,
    oid_field = oid_field
  )
}

arcgis_ids <- function(layer_url, spatial_filter = TRUE) {
  body <- list(
    where = "1=1",
    returnIdsOnly = "true",
    f = "json"
  )

  if (isTRUE(spatial_filter)) {
    body <- c(
      body,
      list(
        geometry = bbox_esri_json,
        geometryType = "esriGeometryEnvelope",
        inSR = "4326",
        spatialRel = "esriSpatialRelIntersects"
      )
    )
  }

  x <- safe_post_json(
    paste0(layer_url, "/query"),
    body
  )

  ids <- x$objectIds

  if (is.null(ids)) {
    # Empty spatial query is allowed; absent objectIds plus no ArcGIS error is
    # interpreted as zero returned features.
    return(integer())
  }

  sort(unique(as.integer(ids)))
}

parse_arcgis_feature_response <- function(resp, query_format, layer_url) {
  raw <- httr2::resp_body_raw(resp)

  # Before giving the payload to GDAL, detect an ArcGIS JSON error body.
  txt_head <- rawToChar(raw[seq_len(min(length(raw), 10000L))])

  if (grepl('"error"\\s*:', txt_head)) {
    parsed <- try(
      jsonlite::fromJSON(rawToChar(raw)),
      silent = TRUE
    )

    if (!inherits(parsed, "try-error") && !is.null(parsed$error)) {
      stop(
        "ArcGIS REST feature-query error from ",
        layer_url,
        ": ",
        paste(c(parsed$error$code, parsed$error$message), collapse = " | ")
      )
    }
  }

  ext <- if (identical(query_format, "geojson")) ".geojson" else ".json"
  tmp <- tempfile(fileext = ext)
  writeBin(raw, tmp)
  on.exit(unlink(tmp), add = TRUE)

  out <- try(
    sf::st_read(
      tmp,
      quiet = TRUE,
      stringsAsFactors = FALSE
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    stop(
      "Could not parse ArcGIS ",
      query_format,
      " feature response from: ",
      layer_url
    )
  }

  out
}

arcgis_feature_chunk_once <- function(
    layer_url,
    ids,
    query_format = "geojson") {

  if (length(ids) == 0L) {
    return(
      sf::st_sf(
        geometry = sf::st_sfc(crs = 4326)
      )
    )
  }

  body <- list(
    where = "1=1",
    objectIds = paste(ids, collapse = ","),
    outFields = "*",
    returnGeometry = "true",
    outSR = "4326",
    f = query_format
  )

  resp <- safe_post(
    paste0(layer_url, "/query"),
    body
  )

  parse_arcgis_feature_response(
    resp,
    query_format,
    layer_url
  )
}

# If a geometry-heavy chunk fails, split it recursively rather than silently
# reducing the requested dataset.
arcgis_feature_chunk <- function(
    layer_url,
    ids,
    query_format = "geojson",
    depth = 0L,
    max_depth = 12L) {

  ans <- try(
    arcgis_feature_chunk_once(
      layer_url,
      ids,
      query_format
    ),
    silent = TRUE
  )

  if (!inherits(ans, "try-error")) {
    return(ans)
  }

  if (length(ids) <= 1L || depth >= max_depth) {
    stop(
      "ArcGIS feature chunk failed even after adaptive splitting. ",
      "Layer: ", layer_url,
      "; OBJECTID(s): ", paste(ids, collapse = ","),
      "; last error: ", as.character(ans)
    )
  }

  cut <- floor(length(ids) / 2)

  left <- arcgis_feature_chunk(
    layer_url,
    ids[seq_len(cut)],
    query_format,
    depth + 1L,
    max_depth
  )

  right <- arcgis_feature_chunk(
    layer_url,
    ids[(cut + 1L):length(ids)],
    query_format,
    depth + 1L,
    max_depth
  )

  dplyr::bind_rows(left, right)
}

arcgis_download <- function(
    name,
    layer_url,
    crop = TRUE,
    chunk_size = default_chunk_size) {

  message("\n--- ", name, " ---")
  message(layer_url)

  meta <- arcgis_layer_metadata(layer_url)

  message("Layer name                : ", meta$name)
  message("Geometry                  : ", meta$geometry_type)
  message("Supported query formats   : ", meta$supported_query_formats)
  message("Selected query format     : ", meta$query_format)
  message("MaxRecordCount            : ", meta$max_record_count)
  message("OID field                 : ", meta$oid_field)

  # National/service count for provenance only; geometry is not downloaded here.
  ids_service <- arcgis_ids(
    layer_url,
    spatial_filter = FALSE
  )

  # Download only features intersecting the buffered Cabo Delgado extent.
  ids_study <- arcgis_ids(
    layer_url,
    spatial_filter = TRUE
  )

  message("Feature IDs in service    : ", length(ids_service))
  message("Feature IDs in 50-km bbox : ", length(ids_study))

  effective_chunk <- max(
    1L,
    min(
      as.integer(chunk_size),
      as.integer(meta$max_record_count %||% chunk_size)
    )
  )

  if (length(ids_study) == 0L) {
    obj <- sf::st_sf(
      geometry = sf::st_sfc(crs = 4326)
    )
  } else {
    groups <- split(
      ids_study,
      ceiling(seq_along(ids_study) / effective_chunk)
    )

    pieces <- purrr::map(
      groups,
      ~ arcgis_feature_chunk(
        layer_url = layer_url,
        ids = .x,
        query_format = meta$query_format
      )
    )

    obj <- dplyr::bind_rows(pieces)
  }

  if (!inherits(obj, "sf")) {
    obj <- sf::st_as_sf(obj)
  }

  if (is.na(sf::st_crs(obj))) {
    sf::st_crs(obj) <- 4326
  }

  if (nrow(obj) > 0L) {
    obj <- obj |>
      sf::st_make_valid() |>
      sf::st_transform(4326)
  }

  n_downloaded <- nrow(obj)

  # Exact local crop is retained as a second independent spatial check.
  if (crop && nrow(obj) > 0L) {
    hit <- lengths(
      sf::st_intersects(
        obj,
        study_bbox_sfc
      )
    ) > 0L

    obj <- obj[hit, , drop = FALSE]
  }

  out_rds <- file.path(
    raw_dir,
    paste0(name, ".rds")
  )

  saveRDS(
    obj,
    out_rds,
    compress = "xz"
  )

  tibble(
    layer = name,
    status = "ok",
    error = NA_character_,
    source_url = layer_url,
    layer_name = meta$name,
    source_geometry_type = meta$geometry_type,
    supported_query_formats = meta$supported_query_formats,
    selected_query_format = meta$query_format,
    max_record_count = meta$max_record_count,
    oid_field = meta$oid_field,
    n_service_ids = length(ids_service),
    n_bbox_ids = length(ids_study),
    n_downloaded = n_downloaded,
    n_study_extent = nrow(obj),
    output_geometry_type = paste(
      unique(as.character(sf::st_geometry_type(obj))),
      collapse = ";"
    ),
    output = out_rds
  )
}

arcgis_download_safe <- function(name, layer_url, crop = TRUE) {
  tryCatch(
    arcgis_download(
      name = name,
      layer_url = layer_url,
      crop = crop
    ),
    error = function(e) {
      message("FAILED: ", name, " -> ", conditionMessage(e))

      tibble(
        layer = name,
        status = "error",
        error = conditionMessage(e),
        source_url = layer_url,
        layer_name = NA_character_,
        source_geometry_type = NA_character_,
        supported_query_formats = NA_character_,
        selected_query_format = NA_character_,
        max_record_count = NA_integer_,
        oid_field = NA_character_,
        n_service_ids = NA_integer_,
        n_bbox_ids = NA_integer_,
        n_downloaded = NA_integer_,
        n_study_extent = NA_integer_,
        output_geometry_type = NA_character_,
        output = NA_character_
      )
    }
  )
}

# 4. Current public MozGIS sources -------------------------------------------
# These explicit URLs are intentional: changes in upstream services should be
# visible in version control and QA rather than silently discovered/replaced.
S <- list(
  roads_primary =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Rede_de_Transportes/FeatureServer/2",

  roads_secondary =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Rede_de_Transportes/FeatureServer/3",

  roads_tertiary =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Rede_de_Transportes/FeatureServer/4",

  roads_vicinal =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Rede_de_Transportes/FeatureServer/5",

  roads_other =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Rede_de_Transportes/FeatureServer/6",

  border_posts =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Rede_de_Transportes/FeatureServer/8",

  aerodromes_transport =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Rede_de_Transportes/FeatureServer/10",

  mining_concessions =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Recursos_Geol%C3%B3gicos/FeatureServer/1",

  hospitals =
    "https://www.mozgis.gov.mz/server/rest/services/MISAU/EDUSUnidadesSanitarias/FeatureServer/1",

  cell_towers_2020 =
    "https://www.mozgis.gov.mz/server/rest/services/Portaldb/ServiTores2020/FeatureServer/0",

  airports_2022 =
    "https://www.mozgis.gov.mz/server/rest/services/MEF/MEFAAPInfraestruturasAeroportuarias/MapServer/1",

  ports =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Mundo/MapServer/26",

  district_capitals =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Topon%C3%ADmia/FeatureServer/1",

  province_capitals =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Topon%C3%ADmia/FeatureServer/0",

  coastline =
    "https://www.mozgis.gov.mz/server/rest/services/Data/Bases/MapServer/11",

  mozambique_2017 =
    "https://www.mozgis.gov.mz/server/rest/services/Data/UnidadesAdministrativasCenso2017/FeatureServer/0",

  provinces_2017 =
    "https://www.mozgis.gov.mz/server/rest/services/Data/UnidadesAdministrativasCenso2017/FeatureServer/1"
)

# Education layers 1..9 are combined and tagged with layer type.
edu_urls <- setNames(
  paste0(
    "https://www.mozgis.gov.mz/server/rest/services/Data/Educa%C3%A7%C3%A3o/FeatureServer/",
    1:9
  ),
  c(
    "creche",
    "primary",
    "secondary",
    "technical",
    "other",
    "orphanage",
    "education_centres",
    "higher_institute",
    "university"
  )
)

# 5. Download all main layers ------------------------------------------------
log_main <- purrr::imap_dfr(
  S,
  ~ arcgis_download_safe(
    name = .y,
    layer_url = .x,
    crop = TRUE
  )
)

# 6. Download and combine education layers -----------------------------------
edu_log <- purrr::imap_dfr(
  edu_urls,
  ~ arcgis_download_safe(
    name = paste0("education_", .y),
    layer_url = .x,
    crop = TRUE
  )
)

edu_ok <- edu_log |>
  filter(
    .data$status == "ok",
    !is.na(.data$output),
    file.exists(.data$output)
  )

edu_pieces <- purrr::map2(
  edu_ok$output,
  sub("^education_", "", edu_ok$layer),
  function(path, nm) {
    x <- readRDS(path)
    if (nrow(x) > 0L) {
      x$education_type <- nm
    }
    x
  }
)

if (length(edu_pieces) > 0L) {
  education_all <- dplyr::bind_rows(edu_pieces)
} else {
  education_all <- sf::st_sf(
    education_type = character(),
    geometry = sf::st_sfc(crs = 4326)
  )
}

education_all_file <- file.path(
  raw_dir,
  "education_all.rds"
)

saveRDS(
  education_all,
  education_all_file,
  compress = "xz"
)

education_combined_log <- tibble(
  layer = "education_all",
  status = if (nrow(edu_ok) > 0L) "ok" else "error",
  error = if (nrow(edu_ok) > 0L) NA_character_ else "No education sublayer downloaded successfully.",
  source_url = paste(edu_urls, collapse = " | "),
  layer_name = "Education combined",
  source_geometry_type = NA_character_,
  supported_query_formats = NA_character_,
  selected_query_format = NA_character_,
  max_record_count = NA_integer_,
  oid_field = NA_character_,
  n_service_ids = NA_integer_,
  n_bbox_ids = NA_integer_,
  n_downloaded = NA_integer_,
  n_study_extent = nrow(education_all),
  output_geometry_type = paste(
    unique(as.character(sf::st_geometry_type(education_all))),
    collapse = ";"
  ),
  output = education_all_file
)

# 7. QA / provenance ---------------------------------------------------------
prov <- bind_rows(
  log_main,
  edu_log,
  education_combined_log
)

readr::write_csv(
  prov,
  file.path(
    qa_dir,
    "00_mozgis_download_manifest.csv"
  ),
  na = ""
)

successful_outputs <- prov |>
  filter(
    .data$status == "ok",
    !is.na(.data$output),
    file.exists(.data$output)
  ) |>
  pull(.data$output) |>
  unique()

md5 <- if (length(successful_outputs) > 0L) {
  tools::md5sum(successful_outputs)
} else {
  character()
}

readr::write_csv(
  tibble(
    file = names(md5),
    md5 = unname(md5)
  ),
  file.path(
    qa_dir,
    "01_download_md5.csv"
  ),
  na = ""
)

schema <- purrr::map_dfr(
  successful_outputs,
  function(p) {
    x <- readRDS(p)

    tibble(
      file = basename(p),
      column = names(sf::st_drop_geometry(x))
    )
  }
)

readr::write_csv(
  schema,
  file.path(
    qa_dir,
    "02_schema.csv"
  ),
  na = ""
)

failed <- prov |>
  filter(.data$status == "error")

readr::write_csv(
  failed,
  file.path(
    qa_dir,
    "03_failed_layers.csv"
  ),
  na = ""
)

message("\nmozgis download summary")
message("Successful rows/layers : ", sum(prov$status == "ok"))
message("Failed rows/layers     : ", sum(prov$status == "error"))
message("QA directory           : ", qa_dir)
message("\n")

print(
  prov |>
    select(
      layer,
      status,
      n_service_ids,
      n_bbox_ids,
      n_study_extent,
      selected_query_format,
      error
    ),
  n = Inf,
  width = Inf
)

if (
  isTRUE(fail_if_any_layer_fails) &&
  nrow(failed) > 0L
) {
  stop(
    nrow(failed),
    " MozGIS layer download(s) failed. Successful layers were preserved. ",
    "Inspect outputs/covariates/05A_mozgis/03_failed_layers.csv."
  )
}

message("\nMozGIS vector download complete with no layer failures.")
