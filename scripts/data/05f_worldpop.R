# Download WorldPop rasters and build population covariates.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(
  terra, sf, dplyr, tidyr, tibble, readr, stringr, rvest, xml2, curl,
  lubridate, here
)

years <- 2017:2024
base_url <- "https://data.worldpop.org/GIS/Population/Global_2015_2030/R2024B"

# Robust download settings. Large WorldPop GeoTIFF transfers occasionally end
# with a transient "Broken pipe". Downloads therefore go to a temporary file,
# are validated, and are promoted to the final .tif only after completion.
download_max_attempts <- 5L
download_connect_timeout <- 60L
download_retry_base_seconds <- 5
boundary_file <- here::here("data", "processed",
  "cabo_delgado_admin2_canonical.rds")
raw_dir <- here::here("data", "raw", "covariates", "worldpop_R2024B")
qa_dir <- here::here("outputs", "covariates", "05F_worldpop")
out <- here::here("data", "processed",
  "covariates_worldpop_population_2017_2024.csv")
for (d in c(raw_dir, qa_dir, dirname(out))) dir.create(d, recursive = TRUE,
  showWarnings = FALSE)

if (!file.exists(boundary_file)) stop("Missing canonical boundary: ",
  boundary_file)
districts <- readRDS(boundary_file) |> sf::st_make_valid() |> sf::st_transform(4326) |> arrange(district_id)
if (nrow(districts) != 17L) stop("Expected 17 canonical districts.")
v <- terra::vect(districts)

# Recursive HTML-directory discovery. WorldPop R2024B stores country rasters
# below YEAR/ISO/v1/100m/{constrained,unconstrained}/ rather than directly in
# YEAR/ISO. We deliberately prefer the constrained 100 m total-population grid.
list_hrefs <- function(url) {
  html <- try(xml2::read_html(url), silent = TRUE)
  if (inherits(html, "try-error")) return(character())
  href <- rvest::html_attr(rvest::html_elements(html, "a"), "href")
  href[!is.na(href) & !href %in% c("../", "./")]
}

join_url <- function(base, href) {
  if (grepl("^https?://", href)) return(href)
  paste0(sub("/*$", "/", base), href)
}

discover_tif <- function(year) {
  root <- paste0(base_url, "/", year, "/MOZ/")
  preferred_dirs <- c(
    paste0(root, "v1/100m/constrained/"),
    paste0(root, "v1/100m/unconstrained/")
  )
  
  # Prefer the known R2024B directory structure. If it changes, fall back to a
  # shallow recursive traversal (maximum depth 4) from the MOZ directory.
  candidate_urls <- character()
  for (u in preferred_dirs) {
    h <- list_hrefs(u)
    tif <- h[stringr::str_detect(h, regex("\\.tif$", ignore_case = TRUE))]
    if (length(tif)) candidate_urls <- c(candidate_urls, vapply(tif,
      function(z) join_url(u, z), character(1)))
  }
  
  if (!length(candidate_urls)) {
    frontier <- root
    seen <- character()
    for (depth in 0:4) {
      next_frontier <- character()
      for (u in frontier) {
        if (u %in% seen) next
        seen <- c(seen, u)
        h <- list_hrefs(u)
        if (!length(h)) next
        full <- vapply(h, function(z) join_url(u, z), character(1))
        candidate_urls <- c(candidate_urls, full[stringr::str_detect(full,
          regex("\\.tif$", ignore_case = TRUE))])
        dirs <- full[stringr::str_detect(h, "/$")]
        next_frontier <- c(next_frontier, dirs)
      }
      frontier <- unique(next_frontier)
      if (!length(frontier)) break
    }
  }
  
  candidate_urls <- unique(candidate_urls)
  candidate_urls <- candidate_urls[!stringr::str_detect(candidate_urls,
    regex("(age|sex|male|female|women|men|birth|preg|dependency)",
    ignore_case = TRUE))]
  if (!length(candidate_urls)) stop("No plausible WorldPop total-population tif discovered for ", year, ". Root: ", root)
  
  fname <- basename(candidate_urls)
  score <- rep(0, length(candidate_urls))
  score <- score + 200 * stringr::str_detect(candidate_urls,
    regex("/constrained/", ignore_case = TRUE))
  score <- score + 100 * stringr::str_detect(fname, regex("_CN_",
    ignore_case = TRUE))
  score <- score + 50  * stringr::str_detect(fname, regex("100m",
    ignore_case = TRUE))
  score <- score + 20  * stringr::str_detect(fname, regex("pop",
    ignore_case = TRUE))
  score <- score + 10  * stringr::str_detect(fname, fixed(as.character(year)))
  chosen_url <- candidate_urls[order(score, decreasing = TRUE)[1]]
  
  tibble(
    year = year,
    root = root,
    candidates = paste(candidate_urls, collapse = " | "),
    chosen = basename(chosen_url),
    url = chosen_url
  )
}

# Robust large-file download helpers -----------------------------------------

remote_content_length <- function(url) {
  h <- curl::new_handle()
  curl::handle_setopt(
    h,
    nobody = TRUE,
    followlocation = TRUE,
    connecttimeout = download_connect_timeout,
    timeout = 60L
  )
  
  res <- try(
    curl::curl_fetch_memory(url, handle = h),
    silent = TRUE
  )
  
  if (inherits(res, "try-error")) {
    return(NA_real_)
  }
  
  hdr <- rawToChar(res$headers)
  
  m <- stringr::str_match_all(
    hdr,
    stringr::regex(
      "(?im)^content-length:\\s*([0-9]+)\\s*$"
    )
  )[[1]]
  
  if (nrow(m) == 0L) {
    return(NA_real_)
  }
  
  suppressWarnings(
    as.numeric(m[nrow(m), 2])
  )
}

validate_worldpop_tif <- function(path, expected_bytes = NA_real_) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  
  size <- as.numeric(file.info(path)$size)
  
  if (!is.finite(size) || size < 1e6) {
    return(FALSE)
  }
  
  # If the server supplied Content-Length, require an exact byte match.
  # This detects the important case where a failed transfer left a large but
  # truncated file that terra might still be able to open partially.
  if (
    is.finite(expected_bytes) &&
    expected_bytes > 0 &&
    size != expected_bytes
  ) {
    return(FALSE)
  }
  
  rr <- try(
    terra::rast(path),
    silent = TRUE
  )
  
  if (inherits(rr, "try-error")) {
    return(FALSE)
  }
  
  if (
    terra::nlyr(rr) != 1L ||
    terra::nrow(rr) <= 0L ||
    terra::ncol(rr) <= 0L
  ) {
    return(FALSE)
  }
  
  # If the server did not provide Content-Length, force a complete raster read
  # once. This is slower, but protects against reusing a large truncated TIFF
  # left by an interrupted download.
  if (!is.finite(expected_bytes)) {
    full_check <- try(
      terra::global(
        rr,
        fun = "sum",
        na.rm = TRUE
      ),
      silent = TRUE
    )
    
    if (inherits(full_check, "try-error")) {
      return(FALSE)
    }
    
    vals <- suppressWarnings(
      as.numeric(
        unlist(full_check, use.names = FALSE)
      )
    )
    
    if (!length(vals) || !any(is.finite(vals))) {
      return(FALSE)
    }
  }
  
  TRUE
}

download_worldpop_tif <- function(
    url,
    dest,
    max_attempts = download_max_attempts
) {
  expected_bytes <- remote_content_length(url)
  
  # Reuse a pre-existing complete raster. This means the successful 2017-2023
  # files from a previous run are not downloaded again.
  if (validate_worldpop_tif(dest, expected_bytes)) {
    return(
      tibble::tibble(
        status = "reused_valid",
        attempts = 0L,
        expected_bytes = expected_bytes,
        local_bytes = as.numeric(file.info(dest)$size),
        file = dest
      )
    )
  }
  
  # A previous version of the script downloaded directly to `dest`, so a
  # Broken-pipe failure can leave a truncated .tif there. Remove it explicitly.
  if (file.exists(dest)) {
    message(
      "Existing WorldPop file is incomplete/invalid; removing: ",
      basename(dest),
      " [",
      format(file.info(dest)$size, big.mark = ","),
      " bytes]"
    )
    unlink(dest)
  }
  
  part <- paste0(dest, ".part.tif")
  
  last_error <- NULL
  
  for (attempt in seq_len(max_attempts)) {
    if (file.exists(part)) {
      unlink(part)
    }
    
    message(
      "Download attempt ",
      attempt,
      "/",
      max_attempts,
      ": ",
      basename(dest)
    )
    
    h <- curl::new_handle()
    
    curl::handle_setopt(
      h,
      followlocation = TRUE,
      connecttimeout = download_connect_timeout,
      timeout = 0L
    )
    
    ans <- tryCatch(
      {
        curl::curl_download(
          url,
          part,
          mode = "wb",
          quiet = FALSE,
          handle = h
        )
        NULL
      },
      error = function(e) e
    )
    
    if (is.null(ans)) {
      if (validate_worldpop_tif(part, expected_bytes)) {
        if (file.exists(dest)) {
          unlink(dest)
        }
        
        promoted <- file.rename(
          part,
          dest
        )
        
        if (!isTRUE(promoted)) {
          stop(
            "Downloaded WorldPop file validated but could not be renamed from ",
            part,
            " to ",
            dest,
            "."
          )
        }
        
        return(
          tibble::tibble(
            status = "downloaded",
            attempts = attempt,
            expected_bytes = expected_bytes,
            local_bytes = as.numeric(file.info(dest)$size),
            file = dest
          )
        )
      }
      
      last_error <- simpleError(
        paste0(
          "Transfer completed but file validation failed. local_bytes=",
          if (file.exists(part)) file.info(part)$size else NA,
          "; expected_bytes=",
          expected_bytes
        )
      )
    } else {
      last_error <- ans
    }
    
    if (file.exists(part)) {
      unlink(part)
    }
    
    if (attempt < max_attempts) {
      wait <- download_retry_base_seconds * 2^(attempt - 1L)
      
      message(
        "WorldPop transfer failed: ",
        conditionMessage(last_error),
        "\nRetrying in ",
        wait,
        " seconds..."
      )
      
      Sys.sleep(wait)
    }
  }
  
  stop(
    "WorldPop download failed after ",
    max_attempts,
    " attempts for:\n",
    url,
    "\nLast error: ",
    conditionMessage(last_error)
  )
}

manifest <- bind_rows(lapply(years, discover_tif))
readr::write_csv(manifest, file.path(qa_dir,
  "00_worldpop_discovery_manifest.csv"), na = "")

# Exact/fractional raster aggregation. WorldPop total-population rasters contain
# estimated persons per pixel. When terra returns polygon-cell fractions, use
# them explicitly; otherwise sum the extracted population values.
weighted_sum_extract <- function(r, polygons) {
  z <- terra::extract(r, polygons, exact = TRUE)
  val_cols <- setdiff(names(z), c("ID", "weight", "fraction"))
  if (length(val_cols) != 1L) stop("Unexpected raster extraction schema: ",
    paste(names(z), collapse = ", "))
  val <- z[[val_cols]]
  wt_name <- intersect(c("weight", "fraction"), names(z))
  z$weighted_value <- if (length(wt_name)) val * z[[wt_name[1]]] else val
  z |>
    as_tibble() |>
    group_by(ID) |>
    summarise(population_worldpop = sum(weighted_value, na.rm = TRUE),
      .groups = "drop")
}

annual_list <- vector("list", length(years))
download_log <- vector("list", length(years))

for (ii in seq_along(years)) {
  yy <- years[ii]
  row <- manifest |> filter(year == yy)
  local <- file.path(raw_dir, basename(row$chosen))
  
  message("\nworldpop ", yy, "")
  message("URL  : ", row$url)
  message("File : ", local)
  
  dl <- download_worldpop_tif(
    url = row$url,
    dest = local
  )
  
  download_log[[ii]] <- dl |>
    mutate(
      year = yy,
      url = row$url,
      .before = 1
    )
  
  r <- terra::rast(local)
  if (terra::nlyr(r) != 1L) stop("Expected one-band WorldPop raster for ", yy)
  
  # Crop first for speed; keep source CRS and reproject vector, not raster.
  vv <- terra::project(v, terra::crs(r))
  rc <- terra::crop(r, terra::ext(vv))
  agg <- weighted_sum_extract(rc, vv)
  if (nrow(agg) != 17L) stop("WorldPop did not return 17 district aggregates for ", yy)
  
  # Area is calculated geodesically/projected independently of raster cells.
  area_km2 <- as.numeric(sf::st_area(sf::st_transform(districts, 32737))) / 1e6
  annual_list[[ii]] <- tibble(
    district_id = districts$district_id[agg$ID],
    canonical_name = districts$canonical_name[agg$ID],
    year = yy,
    population_worldpop = agg$population_worldpop,
    district_area_km2 = area_km2[agg$ID],
    population_density_worldpop = agg$population_worldpop / area_km2[agg$ID],
    worldpop_file = basename(local)
  )
}

download_log <- bind_rows(download_log)

readr::write_csv(
  download_log,
  file.path(qa_dir, "01_worldpop_download_log.csv"),
  na = ""
)

annual <- bind_rows(annual_list) |> arrange(district_id, year)
if (nrow(annual) != 17L * length(years) || anyNA(annual$population_worldpop)) stop("WorldPop annual domain/values incomplete.")

monthly <- tidyr::crossing(
  district_id = districts$district_id,
  year_month = seq.Date(as.Date("2017-01-01"), as.Date("2024-12-01"),
    by = "month")
) |>
  mutate(year = lubridate::year(year_month)) |>
  left_join(districts |> sf::st_drop_geometry() |> select(district_id,
    canonical_name), by = "district_id") |>
  left_join(annual |> select(-canonical_name), by = c("district_id", "year")) |>
  arrange(district_id, year_month)

if (nrow(monthly) != 1632L) stop("WorldPop monthly expansion is not 17 x 96.")
readr::write_csv(monthly, out, na = "")
readr::write_csv(annual, file.path(qa_dir, "02_worldpop_annual_17x8.csv"),
  na = "")
raw_tifs <- list.files(
  raw_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

readr::write_csv(
  tibble::tibble(
    file = raw_tifs,
    bytes = as.numeric(file.info(raw_tifs)$size),
    md5 = unname(tools::md5sum(raw_tifs))
  ),
  file.path(qa_dir, "03_raw_md5.csv"),
  na = ""
)
message("\nworldpop build complete")
message("Years                    : ", paste(range(years), collapse = "-"))
message("Canonical districts      : ", dplyr::n_distinct(annual$district_id))
message("Annual rows              : ", nrow(annual), " (17 x 8)")
message("Monthly rows             : ", nrow(monthly), " (17 x 96)")
message("Reused valid downloads   : ",
  sum(download_log$status == "reused_valid"))
message("New successful downloads : ", sum(download_log$status == "downloaded"))
message("Output                   : ", out)
message("QA directory             : ", qa_dir)
message("")
message(
  "Annual values are repeated within year only; no monthly interpolation",
  " is performed.")
