# Fatality modelling: model search, diagnostics and validation

# 1. Environment setup and package initialisation -----------------------------

if (!requireNamespace("pacman",
  quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  bamlss, gamlss.dist, dplyr, tidyr, ggplot2, ggrepel, ggspatial,
  lubridate, sf, spdep, scales, patchwork, ncvreg, cowplot,
  goftest, ddst, coda, stringr, forcats, here, gt, gtsummary, grDevices,
  RColorBrewer, readr
)

# Load core spatial functions and custom families
source(here("scripts", "analysis", "functions.R"))

theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40", size = 11),
      text = element_text(family = "serif")
    )
)

# RBST_SEED_SHIFT and RBST_OUT_SUFFIX give an independent search chain whose
# outputs go to a separate folder (e.g. 1 and "_chain2").
seed_shift <- as.integer(Sys.getenv("RBST_SEED_SHIFT", "0"))
out_suffix <- Sys.getenv("RBST_OUT_SUFFIX", "")
set.seed(20260920L + seed_shift)

# FAMILIES set: this process fits those families and writes their outputs;
# RBST_DESCRIPTIVE=only runs the descriptive part alone. Unset, all families
# are loaded from the checkpoints and the summary tables are written.
results_root <- here(Sys.getenv("RBST_RESULTS_ROOT", "outputs"))
subset_run <- nzchar(Sys.getenv("FAMILIES"))
desc_only <- identical(Sys.getenv("RBST_DESCRIPTIVE"), "only")
dir_desc <- if (subset_run && !desc_only) {
  file.path(tempdir(), "descriptive")
} else {
  file.path(results_root, paste0("descriptive", out_suffix))
}
dir_mod  <- file.path(results_root, paste0("models", out_suffix))
dir.create(dir_desc, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_mod,  recursive = TRUE, showWarnings = FALSE)

# RBST_TABLES_ONLY=1 rebuilds the LaTeX tables and summary figures from the
# CSV files in dir_mod, without fitting.
if (nzchar(Sys.getenv("RBST_TABLES_ONLY"))) {
  write_result_tables(dir_mod)
  quit(save = "no")
}

# Utility wrappers ------------------------------------------------------------

timeline_x_scale <- function(limits) {
  list(
    scale_x_date(
      limits = limits,
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = c(0.01, 0)
    )
  )
}

save_summary_docx <- function(tbl, path, font_size = 8) {
  if (inherits(tbl, "gtsummary")) {
    tbl <- gtsummary::as_gt(tbl)
  }
  
  tbl |>
    gt::opt_table_font(font = "Times New Roman") |>
    gt::tab_options(table.font.size = font_size) |>
    gt::gtsave(filename = path)
}

# 2. Data pre-processing and spatial weight matrix ----------------------------

analysis <- readRDS(here("data","processed" ,"cabo_delgado_data.rds"))

df_raw <- analysis$panel
districts_sf <- analysis$districts
adjacency <- analysis$adjacency
event_sf <- analysis$events
strategic_geography <- analysis$strategic_geography
mozambique_admin1 <- analysis$mozambique_admin1

# Predictor-to-parameter admissibility table
admissibility_table <- if ("admissibility" %in% names(analysis)) {
  analysis$admissibility
} else if ("eligibility" %in% names(analysis)) {
  analysis[["eligibility"]]
} else {
  stop(
    "Analysis object missing predictor-to-parameter admissibility table.",
    call. = FALSE
  )
}

if (!is.data.frame(admissibility_table)) {
  stop("Admissibility table must be a data frame.", call. = FALSE)
}

event_data <- sf::st_drop_geometry(event_sf)

# Spatial penalty matrix
id_order <- adjacency$district_index$district_id
districts_sf <- districts_sf[match(id_order, districts_sf$district_id), ]

nb_map <- adjacency$nb_queen_connected
adj_mat <- nb2mat(nb_map, style = "B", zero.policy = TRUE)
dimnames(adj_mat) <- list(id_order, id_order)

penalty_mat <- diag(rowSums(adj_mat)) - adj_mat
rownames(penalty_mat) <- colnames(penalty_mat) <- id_order
listw_obj <- nb2listw(nb_map, style = "W", zero.policy = TRUE)

# UTM projection of district centroids
shape_ll <- sf::st_transform(districts_sf, 4326)
bb_ll <- sf::st_bbox(shape_ll)
lon0 <- as.numeric((bb_ll["xmin"] + bb_ll["xmax"]) / 2)
lat0 <- as.numeric((bb_ll["ymin"] + bb_ll["ymax"]) / 2)
utm_zone <- floor((lon0 + 180) / 6) + 1L
utm_epsg <- if (lat0 >= 0) 32600L + utm_zone else 32700L + utm_zone

shape_metric <- sf::st_transform(districts_sf, utm_epsg)
centroids_metric <- sf::st_centroid(sf::st_geometry(shape_metric))
coords_mat <- sf::st_coordinates(centroids_metric)

df_coords <- data.frame(
  district_id = districts_sf$district_id,
  x_km = coords_mat[, 1] / 1000,
  y_km = coords_mat[, 2] / 1000
)

df_panel <- df_raw %>%
  left_join(df_coords, by = "district_id") %>%
  mutate(
    year_month = as.Date(year_month),
    date_obj = year_month,
    fatalities = fatalities_all,
    district_id = factor(district_id, levels = id_order),
    log_pop = log(population)
  )

message(
  "Data preparation complete. Full district-month panel N: ", nrow(df_panel)
)

# 3. Response variable distribution -------------------------------------------

grDevices::cairo_pdf(
  file.path(dir_desc, "Datadistribution1.pdf"), width = 14, height = 3.5
)

op <- par(
  mfrow = c(1, 6), oma = c(2.5, 2.5, 1.4, 0),
  mar = c(1, 1, 0, 0.2), mgp = c(1, 0, 0), tcl = 0.1
)

hist(df_panel$fatalities, main = "", xlab = "")
box()
mtext('(a) Raw monthly fatality counts', side=3, line=0)
mtext('Fatalities', side=1, line=1.2)
mtext('Frequency', side=2, line=1.2)

hist(df_panel$fatalities[df_panel$fatalities>0], main = "", xlab = "")
box()
mtext('(b) Only fatalities > 0', side=3, line=0)
mtext('Fatalities', side=1, line=1.2)

hist(df_panel$fatalities[df_panel$fatalities>1], main = "", xlab = "")
box()
mtext('(c) Only fatalities > 1', side=3, line=0)
mtext('Fatalities', side=1, line=1.2)

hist(
  df_panel$fatalities[df_panel$fatalities > 1 & df_panel$fatalities <= 50],
  main = "", xlab = ""
)
box()
mtext('(d) 1 < Fatalities <= 50', side = 3, line = 0)
mtext('Fatalities', side=1, line=1.2)

hist(log(df_panel$fatalities +1), main = "", xlab = "")
box()
mtext('(e) log transformation', side=3, line=0)
mtext('log(Fatality + 1)', side=1, line=1.2)

p_seq <- seq(0, 1, by = 0.01)
plot(
  p_seq, quantile(df_panel$fatalities, probs = p_seq), type = "s", main = "",
  xlab = "", ylab = ""
)
box()
mtext('(f) Empirical Quantiles', side=3, line=0)
mtext(expression('Probability'~tau), side=1, line=1.2)

par(op)
invisible(dev.off())

# 4. Descriptive statistics tables --------------------------------------------

df_table <- df_panel %>%
  mutate(
    threat_level = factor(
      fatality_threat_level_all,
      levels = c("Low", "Moderate", "Substantial", "Severe", "Critical")
    ),
    event_state_label = factor(
      event_state_all,
      levels = c(
        "no_recorded_event", "recorded_event_zero_reported_fatalities",
        "positive_reported_fatalities"
      ),
      labels = c(
        "No recorded event", "Recorded event, zero reported fatalities",
        "Positive reported fatalities"
      )
    ),
    geographic_position = factor(
      case_when(
        touches_tanzania & coastal ~ "Tanzania border + coast",
        touches_tanzania ~ "Tanzania border",
        touches_nampula & coastal ~ "Nampula border + coast",
        touches_nampula ~ "Nampula border",
        coastal ~ "Coast",
        TRUE ~ "Interior"
      )
    ),
    cyclone_100km = factor(
      cyclone_within_100km, levels = c(FALSE, TRUE), labels = c("No", "Yes")
    ),
    wagner_presence = factor(
      wagner_active, levels = c(FALSE, TRUE), labels = c("No", "Yes")
    ),
    dag_presence = factor(
      dag_active, levels = c(FALSE, TRUE), labels = c("No", "Yes")
    ),
    rwanda_presence = factor(
      rdf_active, levels = c(FALSE, TRUE), labels = c("No", "Yes")
    ),
    samim_presence = factor(
      samim_active, levels = c(FALSE, TRUE), labels = c("No", "Yes")
    ),
    samim_drawdown_label = factor(
      samim_drawdown, levels = c(FALSE, TRUE), labels = c("No", "Yes")
    ),
    financing_stage = factor(
      financing_history_stage,
      levels = c(
        "pre_documented_external_support_milestone",
        "documented_external_support", "al_karrar_regional_coordination",
        "independent_mozambique_province"
      ),
      labels = c(
        "Before external-support milestone", "External support",
        "Al-Karrar regional coordination", "Independent Mozambique province"
      )
    )
  )

table_1 <- df_table %>%
  select(
    event_state_label, geographic_position, cyclone_100km, wagner_presence,
    dag_presence, rwanda_presence, samim_presence, samim_drawdown_label,
    financing_stage, population_density, sex_ratio, dependency_ratio,
    chirps_precip_anom_mm, ntl_log1p, ndvi_anom, road_density_primary_km_km2,
    distance_tanzania_km, distance_nampula_km, distance_coast_km,
    distance_afungi_lng_km, distance_airport_km, distance_port_km,
    mining_active_area_share, threat_level
  ) %>%
  gtsummary::tbl_summary(
    by = threat_level, percent = "column", missing = "ifany",
    label = list(
      event_state_label ~ "Event State",
      geographic_position ~ "Geographic Position",
      cyclone_100km ~ "Cyclone Exposure (<100km)",
      wagner_presence ~ "Wagner Group Active",
      dag_presence ~ "Dyck Advisory Group Active",
      rwanda_presence ~ "Rwanda Security Forces Active",
      samim_presence ~ "SAMIM Active",
      samim_drawdown_label ~ "SAMIM Drawdown Phase",
      financing_stage ~ "Insurgency Financing Stage",
      population_density ~ "Population Density",
      sex_ratio ~ "Sex Ratio",
      dependency_ratio ~ "Dependency Ratio",
      chirps_precip_anom_mm ~ "Precipitation Anomaly (mm)",
      ntl_log1p ~ "Nighttime Lights (Log)",
      ndvi_anom ~ "NDVI Anomaly",
      road_density_primary_km_km2 ~ "Primary Road Density",
      distance_tanzania_km ~ "Distance to Tanzania (km)",
      distance_nampula_km ~ "Distance to Nampula (km)",
      distance_coast_km ~ "Distance to Coast (km)",
      distance_afungi_lng_km ~ "Distance to Afungi LNG (km)",
      distance_airport_km ~ "Distance to Airport (km)",
      distance_port_km ~ "Distance to Port (km)",
      mining_active_area_share ~ "Active Mining Area Share"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2
  ) %>%
  add_overall() %>%
  modify_caption(paste(
    "**Table 1. District-month characteristics by reported fatality",
    "threat level**"
  ))

save_summary_docx(
  table_1, file.path(dir_desc, "Table_1_DistrictMonth_Descriptives.docx")
)

district_table <- df_table %>%
  tbl_cross(
    row = canonical_name, col = threat_level, percent = "column",
    margin = "row", label = list(canonical_name ~ "District")
  ) %>%
  modify_caption(
    "**Table S1. All 17 districts by reported fatality threat level**"
  )

save_summary_docx(
  district_table, file.path(dir_desc, "Table_S1_All_Districts.docx")
)

event_desc <- event_data %>%
  mutate(
    event_threat_level = factor(
      event_threat_level,
      levels = c("Low", "Moderate", "Substantial", "Severe", "Critical")
    ),
    event_type = factor(event_type),
    sub_event_type = factor(sub_event_type),
    event_geographic_position = factor(event_geographic_position)
  )

table_s2 <- event_desc %>%
  select(
    event_type, sub_event_type, event_geographic_position,
    distance_event_airport_km,
    distance_event_port_km, distance_event_border_post_km,
    distance_event_afungi_lng_km,
    distance_event_mining_concession_km, event_threat_level
  ) %>%
  tbl_summary(
    by = event_threat_level, percent = "column", missing = "ifany",
    label = list(
      event_type ~ "Event Type",
      sub_event_type ~ "Sub-Event Type",
      event_geographic_position ~ "Event Geographic Position",
      distance_event_airport_km ~ "Distance to Airport (km)",
      distance_event_port_km ~ "Distance to Port (km)",
      distance_event_border_post_km ~ "Distance to Border Post (km)",
      distance_event_afungi_lng_km ~ "Distance to Afungi LNG (km)",
      distance_event_mining_concession_km ~ "Distance to Mining Concession (km)"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2
  ) %>%
  add_overall() %>%
  modify_caption(
    "**Table S2. Event-level characteristics and strategic proximity**"
  )

save_summary_docx(
  table_s2, file.path(dir_desc, "Table_S2_Event_Descriptives.docx")
)

actor_columns <- intersect(
  c("actor1", "actor2", "assoc_actor_1", "assoc_actor_2"), names(event_data)
)
if (!length(actor_columns)) {
  stop("Event-level data contain no ACLED actor fields.")
}

event_id_col <- intersect(c("event_id_cnty", "event_id"), names(event_data))[1]
if (is.na(event_id_col)) event_id_col <- "event_record_id"

actor_source <- event_data
if (event_id_col == "event_record_id") {
  actor_source$event_record_id <- seq_len(nrow(actor_source))
}

actor_inventory <- actor_source %>%
  select(
    all_of(event_id_col), event_date, fatalities, all_of(actor_columns)
  ) %>%
  pivot_longer(
    cols = all_of(actor_columns), names_to = "actor_field", values_to = "actor"
  ) %>%
  filter(!is.na(actor), nzchar(stringr::str_squish(actor))) %>%
  separate_rows(actor, sep = ";") %>%
  mutate(actor = stringr::str_squish(actor)) %>%
  filter(nzchar(actor)) %>%
  group_by(actor) %>%
  summarise(
    first_recorded = min(event_date, na.rm = TRUE),
    last_recorded = max(event_date, na.rm = TRUE),
    n_events = n_distinct(.data[[event_id_col]]),
    fatalities_in_associated_events = sum(fatalities, na.rm = TRUE),
    actor_fields = paste(sort(unique(actor_field)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(first_recorded, actor)

readr::write_csv(
  actor_inventory, file.path(dir_desc, "Table_S3_ACLED_Actor_Inventory.csv"),
  na = ""
)

actor_gt <- gt::gt(actor_inventory) %>%
  gt::cols_label(
    actor = "ACLED Actor",
    first_recorded = "First Recorded",
    last_recorded = "Last Recorded",
    n_events = "Events",
    fatalities_in_associated_events = "Fatalities in Associated Events",
    actor_fields = "Actor Field(s)"
  ) %>%
  gt::opt_table_font(font = "Times New Roman") %>%
  gt::tab_options(table.font.size = 8)

gt::gtsave(actor_gt, file.path(dir_desc, "Table_S3_ACLED_Actor_Inventory.docx"))

# 5. Temporal analysis, interventions and financing stages --------------------

df_timeline <- df_panel %>%
  group_by(date_obj) %>%
  summarise(total_fatalities = sum(fatalities, na.rm = TRUE), .groups = "drop")

df_year <- df_panel %>%
  group_by(year = lubridate::year(date_obj)) %>%
  summarise(total_fatalities = sum(fatalities, na.rm = TRUE), .groups = "drop")

df_month_season <- df_panel %>%
  mutate(month_lbl = month(date_obj, label = TRUE, abbr = TRUE)) %>%
  group_by(month_lbl) %>%
  summarise(total_fatalities = sum(fatalities, na.rm = TRUE), .groups = "drop")

threshold <-
  mean(df_timeline$total_fatalities) + 1.5 * sd(df_timeline$total_fatalities)
peak_months <- df_timeline %>%
  filter(total_fatalities > threshold) %>%
  arrange(desc(total_fatalities)) %>%
  slice_head(n = 7)

get_headline_desc <- function(target_date, events) {
  top_event <- events %>%
    filter(year_month == target_date) %>%
    arrange(desc(fatalities)) %>%
    slice(1)
  
  if (!nrow(top_event)) return("Violence peak")
  
  event_label <- if ("sub_event_type" %in% names(top_event) &&
                     !is.na(top_event$sub_event_type[1])) {
    as.character(top_event$sub_event_type[1])
  } else {
    as.character(top_event$event_type[1])
  }
  
  district_label <- df_panel %>%
    filter(district_id == top_event$district_id[1]) %>%
    distinct(canonical_name) %>%
    pull(canonical_name)
  
  paste0(
    district_label[1], ": ", event_label, " (", top_event$fatalities[1], "†)"
  )
}

annotations <- peak_months %>%
  rowwise() %>%
  mutate(
    legend_label = get_headline_desc(date_obj, event_data),
    plot_label = legend_label
  ) %>%
  ungroup()

unique_events <- unique(annotations$legend_label)
palette_pts <- brewer.pal(
  n = max(3, length(unique_events)), name = "Dark2"
)[seq_along(unique_events)]
names(palette_pts) <- unique_events

g_year <- ggplot(df_year, aes(x = factor(year), y = total_fatalities)) +
  geom_col(fill = "#2c3e50", width = 0.7, alpha = 0.9) +
  geom_text(
    aes(label = comma(total_fatalities)), vjust = -0.5, size = 3,
    fontface = "bold", family = "serif"
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Fatalities") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, face = "bold", family = "serif"
    )
  )

g_season <- ggplot(df_month_season, aes(x = month_lbl, y = total_fatalities)) +
  geom_col(fill = "#e67e22", width = 0.7, alpha = 0.8) +
  geom_text(
    aes(label = comma(total_fatalities)), vjust = -0.5, size = 2.5,
    family = "serif"
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "serif"))

global_date_limits <- c(min(df_panel$date_obj), max(df_panel$date_obj))

g_timeline <- ggplot() +
  geom_col(
    data = df_timeline, aes(x = date_obj, y = total_fatalities), width = 20,
    alpha = 0.9
  ) +
  geom_smooth(
    data = df_timeline,
    aes(x = date_obj, y = total_fatalities, color = "Trend (Loess)"),
    method = "loess", span = 0.2, fill = "gray70", alpha = 0.3, linewidth = 0.8,
    se = TRUE
  ) +
  geom_point(
    data = annotations,
    aes(x = date_obj, y = total_fatalities, color = legend_label), size = 3,
    stroke = 1, shape = 21, fill = "white", show.legend = FALSE
  ) +
  geom_label_repel(
    data = annotations,
    aes(x = date_obj, y = total_fatalities, label = plot_label), size = 3,
    family = "serif", box.padding = 1.2, point.padding = 0.5, force = 15,
    max.overlaps = Inf, nudge_y = max(df_timeline$total_fatalities) * 0.1,
    direction = "both", segment.color = "gray40", segment.size = 0.3,
    min.segment.length = 0
  ) +
  timeline_x_scale(global_date_limits) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(
    name = NULL, values = c(palette_pts, "Trend (Loess)" = "#c0392b")
  ) +
  labs(
    subtitle = "Observed district-month fatalities and peak events",
    y = "Monthly Fatalities (Total)", x = NULL
  ) +
  theme(
    legend.position = "none", axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), plot.margin = margin(b = 0, r = 10, l = 10)
  )

force_month <- df_panel %>%
  group_by(date_obj) %>%
  summarise(
    `Wagner Group` = any(wagner_active %in% TRUE),
    `Dyck Advisory Group` = any(dag_active %in% TRUE),
    `Rwanda Security Forces (RDF + RNP)` = any(rdf_active %in% TRUE),
    `SAMIM` = any(samim_active %in% TRUE),
    `SAMIM drawdown` = any(samim_drawdown %in% TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-date_obj, names_to = "force", values_to = "active") %>%
  filter(active) %>%
  mutate(evidence_type = "Independent deployment / intervention timeline")

local_force_pattern <- paste(
  c("local force", "local forces", "communal militia", "naparama", "aclln",
    "aclin"),
  collapse = "|"
)
state_force_pattern <- paste(
  c("military forces of mozambique", "police forces of mozambique",
    "mozambique defence", "mozambican military", "rapid intervention",
    "unidade de interv"),
  collapse = "|"
)

event_actor_text <- event_data %>%
  transmute(
    year_month,
    actor_text = apply(
      select(., all_of(actor_columns)), 1,
      function(z) paste(z[!is.na(z)], collapse = "; ")
    )
  ) %>%
  mutate(
    actor_text_lower = stringr::str_to_lower(actor_text),
    local_force_observed = stringr::str_detect(
      actor_text_lower, local_force_pattern
    ),
    moz_state_force_observed = stringr::str_detect(
      actor_text_lower, state_force_pattern
    )
  )

event_force_month <- event_actor_text %>%
  group_by(year_month) %>%
  summarise(
    `Local Forces / Communal Militias` = any(
      local_force_observed, na.rm = TRUE
    ),
    `Mozambican State Forces` = any(moz_state_force_observed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(date_obj = year_month) %>%
  pivot_longer(-date_obj, names_to = "force", values_to = "active") %>%
  filter(active) %>%
  mutate(evidence_type = "Recorded actor involvement")

security_timeline <- bind_rows(force_month, event_force_month) %>%
  mutate(
    force = factor(
      force,
      levels = rev(
        c(
          "Wagner Group", "Dyck Advisory Group", "Mozambican State Forces",
          "Local Forces / Communal Militias",
          "Rwanda Security Forces (RDF + RNP)", "SAMIM", "SAMIM drawdown"
        )
      )
    )
  )

g_security <- ggplot(
  security_timeline, aes(x = date_obj, y = force, shape = evidence_type)
) +
  geom_point(size = 1.7) +
  timeline_x_scale(global_date_limits) +
  labs(x = NULL, y = NULL, shape = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_blank(), axis.title = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "bottom",
    text = element_text(family = "serif")
  )

financing_month <- df_panel %>%
  distinct(date_obj, financing_history_stage) %>%
  mutate(
    financing_stage = factor(
      financing_history_stage,
      levels = c(
        "pre_documented_external_support_milestone",
        "documented_external_support", "al_karrar_regional_coordination",
        "independent_mozambique_province"
      ),
      labels = c(
        "Pre-external support", "External support",
        "Al-Karrar regional coordination", "Independent Mozambique province"
      )
    )
  )

g_financing <- ggplot(
  financing_month,
  aes(x = date_obj, y = "Financing / organisation", fill = financing_stage)
) +
  geom_tile(height = 0.7, width = 31) +
  timeline_x_scale(global_date_limits) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(), legend.position = "bottom",
    text = element_text(family = "serif")
  )

temporal_context_figure <-
  (g_year + g_season) / g_timeline / g_security / g_financing +
  plot_layout(heights = c(1, 2.1, 1.05, 0.65)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E")))

print(temporal_context_figure)

ggsave(
  file.path(dir_desc, "Figure_Temporal_Context.pdf"),
  g_timeline / g_security / g_financing + plot_layout(heights = c(7, 2.6, 0.5)),
  width = 12, height = 10, device = cairo_pdf
)

# 6. Spatial distribution of observed fatalities ------------------------------

df_risk_total <- df_panel %>%
  group_by(district_id) %>%
  summarise(total_risk = sum(fatalities, na.rm = TRUE), .groups = "drop")

map_data_cd <- districts_sf %>%
  left_join(df_risk_total, by = "district_id") %>%
  mutate(
    risk_cat = cut(
      total_risk, breaks = c(0, 50, 100, 500, 1000, 1500, Inf),
      labels = c("<50", "51-100", "101-500", "501-1000", "1001-1500", ">1500"),
      include.lowest = TRUE
    )
  )

map_data_cd_plot <- sf::st_transform(map_data_cd, 4326)
map_label_pts <- sf_label_points(map_data_cd_plot)

as_display_points <- function(x) {
  y <- sf::st_transform(x, 32737)
  sf::st_geometry(y) <- suppressWarnings(
    sf::st_point_on_surface(sf::st_geometry(y))
  )
  sf::st_transform(y, 4326)
}

strategic_points <- bind_rows(
  as_display_points(strategic_geography$airports) |>
    select(asset_type, asset_label, geometry),
  as_display_points(strategic_geography$ports) |>
    select(asset_type, asset_label, geometry),
  as_display_points(strategic_geography$border_posts) |>
    select(asset_type, asset_label, geometry),
  as_display_points(strategic_geography$afungi_lng) |>
    select(asset_type, asset_label, geometry)
)

cd_metric <- sf::st_union(sf::st_transform(districts_sf, 32737))
strategic_metric <- sf::st_transform(strategic_points, 32737)
near_cabo_delgado <- lengths(
  sf::st_is_within_distance(strategic_metric, cd_metric, dist = 5000)
) > 0L
strategic_points <- strategic_points[near_cabo_delgado, ]

strategic_points$asset_type <- factor(
  strategic_points$asset_type,
  levels = c("Afungi LNG", "Airport / aerodrome", "Border post", "Port")
)

p_main <- ggplot(map_data_cd_plot) +
  geom_sf(aes(fill = risk_cat), color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "<50" = "#ffffcc", "51-100" = "#c2e699", "101-500" = "#78c679",
      "501-1000" = "#31a354", "1001-1500" = "#006837", ">1500" = "#003311"
    ),
    na.value = "white", name = "Fatalities"
  ) +
  geom_sf_text(
    data = map_label_pts, aes(label = canonical_name), inherit.aes = FALSE,
    size = 3, angle = -10, color = "black", fontface = "bold",
    check_overlap = TRUE, family = "serif"
  ) +
  geom_sf(
    data = strategic_points, aes(shape = asset_type), inherit.aes = FALSE,
    color = "#6A00A8", fill = "white", size = 3.1, stroke = 1.15
  ) +
  scale_shape_manual(
    values = c(
      "Afungi LNG" = 21, "Airport / aerodrome" = 24, "Border post" = 22,
      "Port" = 23
    ),
    drop = FALSE, name = "Strategic Assets"
  ) +
  annotation_scale(location = "br", width_hint = 0.3, pad_x = unit(0.2, "in")) +
  annotation_north_arrow(
    location = "br", which_north = "true", pad_x = unit(0.2, "in"),
    pad_y = unit(0.3, "in"), style = north_arrow_fancy_orienteering
  ) +
  guides(
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom", legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8), text = element_text(family = "serif")
  )

mozambique_admin1_plot <- sf::st_transform(mozambique_admin1, 4326)
p_inset <- ggplot(mozambique_admin1_plot) +
  geom_sf(fill = "white", color = "gray60", linewidth = 0.18) +
  geom_sf(
    data = mozambique_admin1_plot[mozambique_admin1_plot$is_cabo_delgado, ],
    fill = "black", color = "black", linewidth = 0.18
  ) +
  theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 90, hjust = 1, family = "serif", size = 6
    ),
    axis.text.y = element_text(family = "serif", size = 6),
    axis.title = element_blank()
  )

p_left_panel <- ggdraw(p_main) +
  draw_plot(p_inset, x = 0.005, y = 0.705, width = 0.285, height = 0.285)

ggsave(
  file.path(dir_desc, "Figure_map.pdf"), p_left_panel, width = 12, height = 10,
  device = cairo_pdf
)
unlink(file.path(dir_desc, "Figure_Strategic_Geography.pdf"), force = TRUE)

# 7. Spatial and temporal summaries -------------------------------------------
p_right_panel <-
  g_timeline / g_security / g_financing + plot_layout(heights = c(7, 2.6, 0.5))
ggsave(
  file.path(dir_desc, "Figure_tem.pdf"), p_right_panel, width = 12, height = 10,
  device = cairo_pdf
)

left_col <- plot_grid(
  p_left_panel + theme(plot.margin = margin(0, 0, 0, 0)),
  (g_year | g_season) + theme(plot.margin = margin(10, 0, 0, 0)),
  ncol = 1, rel_heights = c(2.5, 1), labels = c("A)", "B)"), label_size = 14,
  label_fontfamily = "serif"
)

descriptive_figure <- plot_grid(
  left_col, p_right_panel + theme(plot.margin = margin(0, 0, 0, 2)),
  ncol = 2, rel_widths = c(1.5, 2), labels = c("", "C)"), label_size = 14,
  align = "h", axis = "tb", label_fontfamily = "serif"
)

print(descriptive_figure)
ggsave(
  file.path(dir_desc, "Figure_Descriptive.pdf"), descriptive_figure, width = 12,
  height = 8, device = cairo_pdf
)

if (desc_only) quit(save = "no")

# 8. Model estimation ---------------------------------------------------------

# Candidate scope switches and the correlation level flagged for static
# (district-constant) candidates
include_contextual_between <- FALSE
include_sensitivity_candidates <- FALSE
static_correlation_limit <- 0.80

spatial_representation <- c(
  road_access = "distance",
  coast_access = "distance",
  tanzania_access = "distance",
  nampula_access = "distance",
  palma_lng_access = "afungi_distance"
)

# Model-search and posterior sampling configurations
model_search_iter <- as.integer(Sys.getenv("MODEL_SEARCH_ITER", "25000"))
model_search_burn <- as.integer(Sys.getenv("MODEL_SEARCH_BURN", "5000"))
parameter_thin <- 20L
model_thin_cfg <- 1L
frame_cache_limit <- 8L
representative_state_limit <- 80L

posterior_iter <- as.integer(Sys.getenv("POSTERIOR_ITER", "40000"))
posterior_burn <- as.integer(Sys.getenv("POSTERIOR_BURN", "5000"))
posterior_thin_cfg <- 10L
posterior_chains_cfg <- 4L
posterior_cores_cfg <- as.integer(Sys.getenv("POSTERIOR_CORES", "4"))
posterior_segment_cfg <- as.integer(Sys.getenv("POSTERIOR_SEGMENT", "2000"))
chain_seed_cfg <- 20260918L + seed_shift

# Primary empirical specification
use_population_offset <- FALSE
include_spatiotemporal_tensor <- FALSE

# Structural priors: N(0, 5^2) for intercepts and year effects of every
# parameter; exponential (PC) prior for the spatial standard deviation with
# P(sqrt(v_S) > 3) = 0.05, i.e. the bamlss "sd" prior with theta below
structural_fixed_sd <- 5
mrf_sd_upper <- 3
mrf_sd_prob <- 0.05
mrf_pc_theta <- (mrf_sd_upper / -log(mrf_sd_prob))^2
options(rbst.structural_sd = structural_fixed_sd, rbst.mrf_theta = mrf_pc_theta)

# Priors for predictors
candidate_linear_sd <- 1.0
candidate_smooth_null_sd <- 2.5
candidate_smooth_penalized_sd <- 1.0

diagnostic_seed <- 20260913L
watson_mc_cfg <- 1999L
tb2_mc_cfg <- 1999L
tb2_m_cfg <- NULL

z_score <- function(x, label) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  sigma <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
    stop(
      "Cannot standardise '", label, "': non-finite or zero variance.",
      call. = FALSE
    )
  }
  (x - mu) / sigma
}

# Transform lagged count variables and decompose time-varying covariates
df_panel <- df_panel %>%
  group_by(district_id) %>%
  mutate(across(
    c(
      population_density, dependency_ratio, chirps_precip_anom_mm,
      ndvi_anom_lag1, ntl_log1p_lag1
    ),
    list(between = ~mean(.x, na.rm = TRUE),
         within = ~.x - mean(.x, na.rm = TRUE))
  )) %>%
  ungroup() %>%
  mutate(
    log_population_z = z_score(log(population), "log_population"),
    across(
      c(
        fatalities_all_lag1, neighbor_fatalities_lag1,
        reporting_events_loo_lag1, reporting_source_diversity_loo_lag1
      ),
      ~z_score(log1p(pmax(.x, 0)), cur_column()), .names = "{.col}_log1p_z"
    ),
    across(
      ends_with("_within") | ends_with("_between"),
      ~z_score(.x, cur_column()), .names = "{.col}_z"
    )
  )

# Force-deployment phases
intervention_fields <- c("rdf_active", "samim_active", "samim_drawdown")
if (anyNA(df_panel[intervention_fields])) {
  stop("Intervention indicators contain missing values.", call. = FALSE)
}
samim_without_rdf <- df_panel$samim_active %in% TRUE &
  !(df_panel$rdf_active %in% TRUE)
if (any(samim_without_rdf)) {
  stop("SAMIM-active months without Rwanda-active months found.", call. = FALSE)
}

df_panel <- df_panel %>%
  mutate(
    foreign_force_phase = case_when(
      rdf_active %in% TRUE &
        samim_active %in% TRUE &
        samim_drawdown %in% TRUE ~ "rwanda_samim_drawdown",
      rdf_active %in% TRUE & samim_active %in% TRUE ~ "rwanda_samim",
      rdf_active %in% TRUE & !(samim_active %in% TRUE) ~ "rwanda_only",
      TRUE ~ "pre_joint_intervention"
    ),
    phase_rwanda_samim = as.integer(foreign_force_phase == "rwanda_samim"),
    phase_rwanda_samim_drawdown = as.integer(
      foreign_force_phase == "rwanda_samim_drawdown"
    ),
    phase_rwanda_only = as.integer(foreign_force_phase == "rwanda_only"),
    mining_concession_present = as.integer(mining_active_license_count > 0)
  )

# Model variable = concept used in the predictor-to-parameter
# admissibility table
candidate_concepts <- c(
  fatalities_all_lag1_log1p_z = "fatalities_all_lag1",
  neighbor_fatalities_lag1_log1p_z = "neighbor_fatalities_lag1",
  log_population_z = "population_density",
  dependency_ratio_within_z = "dependency_ratio",
  chirps_precip_anom_mm_within_z = "chirps_precip_anom_mm",
  ndvi_anom_lag1_within_z = "ndvi_anom",
  ntl_log1p_lag1_within_z = "ntl_log1p",
  cyclone_within_100km = "cyclone_within_100km",
  terrain_ruggedness_z = "terrain_ruggedness",
  distance_primary_road_km_z = "distance_primary_road_km",
  road_density_primary_km_km2_z = "road_density_primary_km_km2",
  distance_coast_km_z = "distance_coast_km",
  coastal = "coastal",
  distance_tanzania_km_z = "distance_tanzania_km",
  touches_tanzania = "touches_tanzania",
  distance_afungi_lng_km_z = "distance_afungi_lng_km",
  distance_palma_km_z = "distance_palma_km",
  palma_neighbor_raw = "palma_neighbor_raw",
  distance_port_km_z = "distance_port_km",
  phase_rwanda_samim = "rdf_active",
  phase_rwanda_samim_drawdown = "rdf_active",
  phase_rwanda_only = "rdf_active",
  dag_active = "dag_active",
  reporting_events_loo_lag1_log1p_z = "reporting_events_loo_lag1",
  reporting_source_diversity_loo_lag1_log1p_z =
    "reporting_source_diversity_loo_lag1",
  population_density_between_z = "population_density",
  dependency_ratio_between_z = "dependency_ratio",
  chirps_precip_anom_mm_between_z = "chirps_precip_anom_mm",
  ndvi_anom_lag1_between_z = "ndvi_anom",
  ntl_log1p_lag1_between_z = "ntl_log1p",
  population_density_within_z = "population_density",
  touches_nampula = "touches_nampula",
  distance_nampula_km_z = "distance_nampula_km",
  distance_airport_km_z = "distance_airport_km",
  distance_border_post_km_z = "distance_border_post_km",
  mining_concession_present = "mining_active_license_count",
  wagner_active = "wagner_active"
)

# Alternative representations of one mechanism and candidates outside the
# primary set (contextual district means, sensitivity distances)
candidate_roles <- tibble::tribble(
  ~variable, ~mechanism, ~representation, ~scope,
  "distance_primary_road_km_z", "road_access", "distance", "primary",
  "road_density_primary_km_km2_z", "road_access", "density", "primary",
  "distance_coast_km_z", "coast_access", "distance", "primary",
  "coastal", "coast_access", "indicator", "primary",
  "distance_tanzania_km_z", "tanzania_access", "distance", "primary",
  "touches_tanzania", "tanzania_access", "indicator", "primary",
  "distance_nampula_km_z", "nampula_access", "distance", "primary",
  "touches_nampula", "nampula_access", "indicator", "primary",
  "distance_afungi_lng_km_z", "palma_lng_access", "afungi_distance",
  "primary",
  "distance_palma_km_z", "palma_lng_access", "palma_distance", "primary",
  "palma_neighbor_raw", "palma_lng_access", "palma_neighbor", "primary",
  "distance_airport_km_z", NA, NA, "sensitivity",
  "distance_border_post_km_z", NA, NA, "sensitivity",
  "population_density_between_z", NA, NA, "contextual",
  "dependency_ratio_between_z", NA, NA, "contextual",
  "chirps_precip_anom_mm_between_z", NA, NA, "contextual",
  "ndvi_anom_lag1_between_z", NA, NA, "contextual",
  "ntl_log1p_lag1_between_z", NA, NA, "contextual"
)

candidate_spec <- tibble::tibble(
  concept = sub("(_log1p)?_z$", "", names(candidate_concepts)),
  variable = names(candidate_concepts),
  admissibility_concept = unname(candidate_concepts)
) %>%
  left_join(candidate_roles, by = "variable") %>%
  mutate(scope = dplyr::coalesce(scope, "primary"))

representation_choice <- tibble(
  mechanism = names(spatial_representation),
  selected_representation = unname(spatial_representation)
)

candidate_spec <- candidate_spec %>%
  left_join(representation_choice, by = "mechanism") %>%
  filter(is.na(mechanism) | representation == selected_representation) %>%
  filter(scope == "primary" |
           (scope == "contextual" & include_contextual_between) |
           (scope == "sensitivity" & include_sensitivity_candidates)) %>%
  select(
    concept, variable, admissibility_concept, scope, mechanism, representation
  )

if (anyDuplicated(candidate_spec$variable)) {
  stop(
    "Candidate specification contains duplicated model variables.",
    call. = FALSE
  )
}
candidate_pool <- candidate_spec$variable

missing_candidate_variables <- setdiff(candidate_pool, names(df_panel))
if (length(missing_candidate_variables)) {
  stop(
    "Model data are missing candidate variable(s): ",
    paste(missing_candidate_variables, collapse = ", "),
    call. = FALSE
  )
}

preferred_admissibility_columns <- c(
  "variable", "admissible_mu", "admissible_sigma", "admissible_nu",
  "admissible_tau"
)
legacy_admissibility_columns <- c(
  "variable", "elig_mu", "elig_sigma", "elig_nu", "elig_tau"
)

if (all(preferred_admissibility_columns %in% names(admissibility_table))) {
  admissibility_lookup <- admissibility_table %>%
    select(
      admissibility_concept = variable,
      mu = admissible_mu, sigma = admissible_sigma,
      nu = admissible_nu, tau = admissible_tau
    )
} else if (all(legacy_admissibility_columns %in% names(admissibility_table))) {
  admissibility_lookup <- admissibility_table %>%
    select(
      admissibility_concept = variable,
      mu = elig_mu, sigma = elig_sigma, nu = elig_nu, tau = elig_tau
    )
} else {
  stop(
    "The admissibility table lacks the admissible_* (or elig_*) columns.",
    call. = FALSE
  )
}

candidate_admissibility_df <- candidate_spec %>%
  left_join(admissibility_lookup, by = "admissibility_concept")

admissibility_flags <- candidate_admissibility_df %>% select(mu, sigma, nu, tau)
if (any(!stats::complete.cases(admissibility_flags))) {
  incomplete_concepts <- candidate_admissibility_df %>%
    filter(!if_all(c(mu, sigma, nu, tau), ~ !is.na(.x))) %>%
    pull(admissibility_concept) %>%
    unique()
  stop(
    "Predictor-to-parameter admissibility matrix is incomplete for: ",
    paste(incomplete_concepts, collapse = ", "),
    call. = FALSE
  )
}

candidate_admissibility <- as.matrix(
  candidate_admissibility_df %>% select(mu, sigma, nu, tau)
)
rownames(candidate_admissibility) <- candidate_admissibility_df$variable

# The admissibility table is read and checked above, but every parameter is
# then opened to every candidate. Remove this line to apply the table.
candidate_admissibility[, c("mu", "sigma", "nu", "tau")] <- TRUE

is_time_invariant <- vapply(
  candidate_pool,
  function(v) {
    by_district <- split(df_panel[[v]], df_panel$district_id)
    all(
      vapply(
        by_district, function(z) length(unique(z[!is.na(z)])) <= 1L, logical(1)
      )
    )
  },
  logical(1)
)

# Static candidates are informed by the districts only. Of each pair with
# |r| >= static_correlation_limit across districts, the one listed first in
# candidate_concepts is kept and the other leaves the candidate set
static_pruning <- data.frame(
  dropped = character(0), kept = character(0), correlation = numeric(0)
)
static_first <- candidate_pool[is_time_invariant]
if (length(static_first) >= 2L) {
  static_data <- df_panel %>%
    group_by(district_id) %>%
    summarise(
      across(all_of(static_first), ~ as.numeric(dplyr::first(.x))),
      .groups = "drop"
    )
  static_cor <- stats::cor(
    as.data.frame(static_data[static_first]), use = "pairwise.complete.obs"
  )
  kept <- character(0)
  for (v in static_first) {
    r <- static_cor[v, kept]
    hit <- kept[is.finite(r) & abs(r) >= static_correlation_limit]
    if (length(hit)) {
      static_pruning <- rbind(static_pruning, data.frame(
        dropped = v, kept = hit[1L], correlation = static_cor[v, hit[1L]]
      ))
    } else {
      kept <- c(kept, v)
    }
  }
}
write.csv(
  static_pruning, file.path(dir_mod, "static_candidate_pruning.csv"),
  row.names = FALSE
)
if (nrow(static_pruning)) {
  message(
    "Static candidates removed for |r| >= ", static_correlation_limit, ": ",
    paste0(static_pruning$dropped, " (kept ", static_pruning$kept, ")",
           collapse = "; ")
  )
  keep_candidate <- !candidate_pool %in% static_pruning$dropped
  candidate_pool <- candidate_pool[keep_candidate]
  is_time_invariant <- is_time_invariant[keep_candidate]
  candidate_spec <- candidate_spec[
    candidate_spec$variable %in% candidate_pool, , drop = FALSE
  ]
  candidate_admissibility <- candidate_admissibility[
    candidate_pool, , drop = FALSE
  ]
}

# A candidate that the base predictors of a parameter already reproduce
# (R2 >= base_r2_limit) has no information of its own there. The base
# predictors are district and year effects in mu, year effects in nu and the
# intercept in sigma and tau. A static covariate has R2 = 1 in mu
base_r2_limit <- as.numeric(Sys.getenv("BASE_R2_LIMIT", "0.95"))
base_dat <- df_panel
base_dat$year_factor <- factor(lubridate::year(base_dat$year_month))
base_dat$district_id <- factor(base_dat$district_id)
base_r2 <- function(x, form) {
  ok <- is.finite(x)
  M <- stats::model.matrix(form, base_dat[ok, , drop = FALSE])
  tss <- sum((x[ok] - mean(x[ok]))^2)
  if (tss > 0) 1 - sum(stats::lm.fit(M, x[ok])$residuals^2) / tss else NA_real_
}
base_redundancy <- do.call(rbind, lapply(candidate_pool, function(v) {
  x <- base_dat[[v]]
  num <- is.numeric(x) || is.logical(x)
  data.frame(
    variable = v,
    R2_mu = if (num) base_r2(as.numeric(x), ~ district_id + year_factor)
    else NA_real_,
    R2_nu = if (num) base_r2(as.numeric(x), ~ year_factor) else NA_real_
  )
}))
for (k in c("mu", "nu")) {
  out <- base_redundancy[[paste0("R2_", k)]] >= base_r2_limit
  out[is.na(out)] <- FALSE
  candidate_admissibility[base_redundancy$variable[out], k] <- FALSE
}
write.csv(
  base_redundancy, file.path(dir_mod, "base_redundancy.csv"),
  row.names = FALSE
)
print(base_redundancy, digits = 3, row.names = FALSE)

write.csv(
  candidate_spec, file.path(dir_mod, "candidate_specification.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    variable = rownames(candidate_admissibility), candidate_admissibility,
    check.names = FALSE
  ),
  file.path(dir_mod, "candidate_admissibility.csv"),
  row.names = FALSE
)

static_candidates <- candidate_pool[is_time_invariant]

model_candidate_pool <- ifelse(
  candidate_pool %in% static_candidates,
  paste0("L::", candidate_pool),
  candidate_pool
)

candidate_component_policy <- data.frame(
  variable = candidate_pool,
  time_invariant = is_time_invariant,
  search_component = ifelse(
    is_time_invariant, "linear_only", "linear_or_nonlinear"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  candidate_component_policy,
  file.path(dir_mod, "candidate_component_policy.csv"),
  row.names = FALSE
)

df_model <- df_panel %>%
  filter(
    if_all(all_of(candidate_pool), ~ !is.na(.x)),
    is.finite(fatalities),
    is.finite(population),
    population > 0
  ) %>%
  mutate(
    t_index = as.numeric(factor(year_month, levels = sort(unique(year_month)))),
    year_factor = factor(lubridate::year(year_month)),
    district_id = factor(district_id, levels = id_order),
    log_pop = log(population)
  ) %>%
  droplevels()

message(
  "Model sample N: ", nrow(df_model), "; months: ",
  n_distinct(df_model$year_month), "; candidates: ", length(candidate_pool)
)

# Share of each time-varying candidate explained by districts, by years and
# by years given districts. Diagnostic only; the candidate set is unchanged
r2_gain <- function(x, base, full) {
  sse <- function(M) sum(stats::lm.fit(M, x)$residuals^2)
  s0 <- sse(base)
  if (s0 > 0) 1 - sse(full) / s0 else NA_real_
}
X_one <- matrix(1, nrow(df_model), 1L)
X_d <- stats::model.matrix(~ district_id, df_model)
X_y <- stats::model.matrix(~ year_factor, df_model)
X_dy <- stats::model.matrix(~ district_id + year_factor, df_model)
tv_numeric <- candidate_pool[!is_time_invariant & vapply(
  candidate_pool, function(v) is.numeric(df_model[[v]]) ||
    is.logical(df_model[[v]]), logical(1)
)]
time_collinearity <- do.call(rbind, lapply(tv_numeric, function(v) {
  x <- as.numeric(df_model[[v]])
  z <- (x - mean(x)) / stats::sd(x)
  Z <- data.matrix(df_model[setdiff(tv_numeric, v)])
  data.frame(
    variable = v,
    R2_district = r2_gain(x, X_one, X_d),
    R2_year = r2_gain(x, X_one, X_y),
    R2_year_given_district = r2_gain(x, X_d, X_dy),
    R2_others_given_dy = r2_gain(x, X_dy, cbind(X_dy, Z)),
    skewness = mean(z^3),
    max_abs_z = max(abs(z))
  )
}))
write.csv(
  time_collinearity, file.path(dir_mod, "time_collinearity.csv"),
  row.names = FALSE
)
print(time_collinearity, digits = 3, row.names = FALSE)
if (isTRUE(as.logical(Sys.getenv("DIAGNOSE_ONLY", "FALSE")))) {
  quit(save = "no")
}

# Base predictors for all models ----------------------------------------------
# xt is kept as a short object: bamlss drops spatial terms by their label, and
# a label longer than 60 characters breaks across lines and is not removed
mrf_xt <- list(penalty = penalty_mat, prior = "sd", theta = mrf_pc_theta)
base_mu_formula <- fatalities ~
  year_factor +
  s(district_id, bs = "mrf", xt = mrf_xt)

if (isTRUE(include_spatiotemporal_tensor)) {
  base_mu_formula <- stats::update(
    base_mu_formula,
    . ~ . +
      ti(t_index, x_km, y_km, bs = c("ps", "tp"), d = c(1, 2), k = c(5, 5))
  )
}

f_sigma <- ~ 1
f_nu <- ~ year_factor
f_tau <- ~ 1

auxiliary_formulas <- list(sigma = f_sigma, nu = f_nu, tau = f_tau, s = f_nu)

model_offsets <- if (isTRUE(use_population_offset)) {
  list(
    ZINB = list(mu = df_model$log_pop),
    ZANBI = list(mu = df_model$log_pop),
    PIG = list(mu = df_model$log_pop),
    ZIPIG = list(mu = df_model$log_pop),
    DEL = list(mu = df_model$log_pop),
    SICHEL = list(mu = df_model$log_pop),
    ZOINB = list(mu = df_model$log_pop),
    RSB = list(lambda = df_model$log_pop)
  )
} else {
  NULL
}

all_families <- c(
  "DEL", "SICHEL", "ZOINB", "RSB", "ZIPIG", "PIG", "ZANBI", "ZINB"
)

# FAMILIES="DEL,SICHEL" fits only those families in this process; unset,
# all of them are fitted and assembled.
run_families <- Sys.getenv("FAMILIES", "")
run_families <- if (nzchar(run_families)) {
  intersect(all_families, trimws(strsplit(run_families, ",")[[1]]))
} else {
  all_families
}

checkpoint_path <- Sys.getenv(
  "CHECKPOINT_DIR", file.path(dir_mod, "checkpoints")
)

# RBST_POOL_DIRS lists checkpoint folders of independent searches (other
# seeds) whose allocation traces are pooled with this one.
pool_dirs <- trimws(strsplit(Sys.getenv("RBST_POOL_DIRS", ""), ",")[[1]])
pool_dirs <- pool_dirs[nzchar(pool_dirs)]

# RBST_LINEARIZE names covariates whose selected smooth effect is replaced by
# a linear one, as a sensitivity analysis (use with RBST_OUT_SUFFIX).
linearize_vars <- trimws(strsplit(Sys.getenv("RBST_LINEARIZE", ""), ",")[[1]])
linearize_vars <- linearize_vars[nzchar(linearize_vars)]

fit_results <- fit_all_spatial_families(
  families = run_families, base_formula = base_mu_formula,
  candidate_pool = model_candidate_pool,
  data = df_model, penalty_mat = penalty_mat,
  auxiliary_formulas = auxiliary_formulas,
  offsets = model_offsets, candidate_admissibility = candidate_admissibility,
  n_iter = model_search_iter, n_burn = model_search_burn,
  thin = parameter_thin, model_thin = model_thin_cfg,
  optimizer_maxit = 400L, optimizer_eps = .Machine$double.eps^0.25,
  optimizer_nu = TRUE,
  within_sweeps = 5L, a_V = 0.5, a_A = 0.2,
  c_linear = candidate_linear_sd,
  c0 = candidate_smooth_null_sd,
  cplus = candidate_smooth_penalized_sd,
  a_lambda = 1, b_lambda = 1, epsilon_birth = 0.25, omega_birth = 0.25,
  proposal_ridge = 1e-6,
  full_cross_information = TRUE, fd_eps = 1e-4, swap_weight = 1,
  lambda_prior_mix = 0.10, beta_prior_mix = 0.15, proposal_link_rms = 1.0,
  lambda_opt_maxit = 25L,
  cache_frames = TRUE, cache_blocks = TRUE,
  max_frame_cache = frame_cache_limit,
  max_representative_states = representative_state_limit,
  store_theta = TRUE,
  pip_threshold = 0.5, support_threshold = 0.30,
  continue_on_family_error = FALSE,
  run_posterior_sampling = !nzchar(Sys.getenv("RBST_SKIP_POSTERIOR")),
  posterior_n_iter = posterior_iter, posterior_n_burn = posterior_burn,
  posterior_thin = posterior_thin_cfg, posterior_chains = posterior_chains_cfg,
  posterior_cores = posterior_cores_cfg, chain_seed = chain_seed_cfg,
  reference_min_chains = posterior_chains_cfg, reference_min_draws = 400L,
  ess_reference = 400, split_drift_reference = 0.50, rhat_reference = 1.01,
  checkpoint_dir = checkpoint_path,
  posterior_segment = posterior_segment_cfg,
  pooled_search_dirs = if (length(pool_dirs)) pool_dirs,
  linearize = linearize_vars
)

# Summary files of a partial run go to a separate folder
dir_out <- dir_mod
full_run <- setequal(run_families, all_families)
if (!full_run) {
  dir_out <- file.path(dir_mod, "partial", paste(run_families, collapse = "_"))
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  message("Fitted: ", paste(run_families, collapse = ", "), ".")
}

fitted_models <- fit_results$models
selected_vars_list <- fit_results$selections
supported_vars_list <- fit_results$supported_selections
selection_types_list <- fit_results$selection_types
effect_evidence_list <- fit_results$effect_evidence
model_formulas <- fit_results$formulas
search_results <- fit_results$model_search
family_status <- fit_results$family_status

saveRDS(search_results, file.path(dir_out, "model_search_results.rds"))
write.csv(
  family_status, file.path(dir_out, "family_status.csv"), row.names = FALSE
)
saveRDS(
  fit_results$diagnostics, file.path(dir_out, "model_search_diagnostics.rds")
)

mixing_summary <- dplyr::bind_rows(lapply(names(search_results), function(fam) {
  d <- search_results[[fam]]$diagnostics
  ar <- d$move_acceptance_rate
  
  data.frame(
    Family = fam,
    Birth_accept = unname(ar["birth"]),
    Death_accept = unname(ar["death"]),
    Swap_accept = unname(ar["swap"]),
    Birth_attempted = unname(d$move_attempted["birth"]),
    Death_attempted = unname(d$move_attempted["death"]),
    Swap_attempted = unname(d$move_attempted["swap"]),
    Proposal_failure_rate = d$proposal_failures / max(1, sum(d$move_attempted)),
    Retained_model_draws = d$retained_model_draws,
    Unique_models_postburn = d$unique_models_visited_postburn,
    Postburn_model_switches = d$postburn_model_switches,
    ESS_model_size = d$ess_model_size,
    Max_association_PIP_MCSE = d$max_association_pip_mcse,
    Max_nonlinear_PIP_MCSE = d$max_nonlinear_pip_mcse,
    Selection_threshold_stable = isTRUE(d$selection_threshold_stable),
    Allocation_diagnostics_within_reference = isTRUE(
      d$allocation_diagnostics_within_reference
    ),
    Within_success_rate = d$within_model_success_rate,
    Elapsed_sec = d$elapsed_seconds,
    Frame_cache_limit = d$frame_cache_limit,
    Frame_cache_peak = d$frame_cache_peak,
    Frame_cache_evictions = d$frame_cache_evictions,
    Retained_theta_states = d$retained_theta_states,
    Representative_states = d$representative_states_retained,
    Diagnostics_outside_reference = paste(unique(c(
      d$kernel_diagnostics_outside_reference,
      d$mixing_diagnostics_outside_reference,
      if (!isTRUE(d$selection_threshold_stable)) {
        "selection_threshold_instability"
      }
    )), collapse = "; "),
    stringsAsFactors = FALSE
  )
}))

if (nrow(mixing_summary) > 0L) {
  print(mixing_summary)
  write.csv(
    mixing_summary, file.path(dir_out, "model_space_mixing_summary.csv"),
    row.names = FALSE
  )
}

pip_mcse_table <- dplyr::bind_rows(lapply(names(search_results), function(fam) {
  fit <- search_results[[fam]]
  pip <- fit$posterior_inclusion_probabilities
  pip_nl <- fit$posterior_nonlinear_probabilities
  pip_nl_cond <- fit$posterior_nonlinear_given_association
  
  out <- as.data.frame(as.table(pip), stringsAsFactors = FALSE)
  names(out) <- c("Covariate", "Parameter", "PIP_Association")
  key <- paste(out$Covariate, out$Parameter, sep = "::")
  
  out$MCSE_Association <- unname(fit$diagnostics$association_pip_mcse[key])
  out$PIP_Nonlinear <- mapply(function(v, p) {
    if (!is.null(pip_nl) &&
        v %in% rownames(pip_nl) &&
        p %in% colnames(pip_nl)) pip_nl[v, p] else NA_real_
  }, out$Covariate, out$Parameter)
  
  out$MCSE_Nonlinear <- unname(fit$diagnostics$nonlinear_pip_mcse[key])
  out$Pr_Nonlinear_Given_Association <- mapply(function(v, p) {
    if (!is.null(pip_nl_cond) &&
        v %in% rownames(pip_nl_cond) &&
        p %in% colnames(pip_nl_cond)) pip_nl_cond[v, p] else NA_real_
  }, out$Covariate, out$Parameter)
  
  out$Family <- fam
  out$Model_Space_Diagnostics <- ifelse(
    isTRUE(fit$diagnostics$allocation_diagnostics_within_reference),
    "within reference values", "outside reference values"
  )
  
  cols <- c(
    "Family", "Model_Space_Diagnostics", "Parameter", "Covariate",
    "PIP_Association", "MCSE_Association", "PIP_Nonlinear", "MCSE_Nonlinear",
    "Pr_Nonlinear_Given_Association"
  )
  out[, cols]
}))

if (nrow(pip_mcse_table) > 0L) {
  pip_mcse_table <- pip_mcse_table |>
    dplyr::mutate(Association_Evidence = dplyr::case_when(
      !is.finite(PIP_Association) ~
        "allocation not admissible or not estimable",
      PIP_Association > 0.5 ~ "median probability model",
      PIP_Association >= 0.30 ~ "secondary posterior support",
      TRUE ~ "low posterior inclusion support"
    ))
  
  write.csv(
    pip_mcse_table, file.path(dir_out, "pip_mcse.csv"), row.names = FALSE
  )
  
  supported_association_table <- pip_mcse_table |>
    dplyr::filter(is.finite(PIP_Association), PIP_Association >= 0.30) |>
    dplyr::arrange(Family, Parameter, dplyr::desc(PIP_Association))
  
  write.csv(
    supported_association_table,
    file.path(dir_out, "supported_associations_pip30.csv"), row.names = FALSE
  )
}

effect_evidence_table <-
  dplyr::bind_rows(lapply(names(effect_evidence_list), function(fam) {
    z <- effect_evidence_list[[fam]]
    if (is.null(z) || !nrow(z)) return(NULL)
    z$Family <- fam
    z[, c("Family", setdiff(names(z), "Family")), drop = FALSE]
  }))

if (nrow(effect_evidence_table) > 0L) {
  write.csv(
    effect_evidence_table, file.path(dir_out, "effect_evidence.csv"),
    row.names = FALSE
  )
}

selection_type_table <-
  dplyr::bind_rows(lapply(names(selection_types_list), function(fam) {
    tp <- selection_types_list[[fam]]
    dplyr::bind_rows(lapply(names(tp), function(par) {
      z <- tp[[par]]
      if (is.null(z) || !length(z)) return(NULL)
      data.frame(
        Family = fam, Parameter = par, Covariate = names(z),
        Selected_Form = unname(z), stringsAsFactors = FALSE
      )
    }))
  }))

if (nrow(selection_type_table) > 0L) {
  write.csv(
    selection_type_table, file.path(dir_out, "selected_effect_forms.csv"),
    row.names = FALSE
  )
}

run_manifest <- data.frame(
  Analysis_Unit = "district_month", N_Analysis_Rows = nrow(df_model),
  N_Districts = dplyr::n_distinct(df_model$district_id),
  N_Months = dplyr::n_distinct(df_model$year_month),
  Candidate_Count = length(candidate_pool),
  Families = paste(run_families, collapse = ";"),
  Model_Search_Iterations = model_search_iter,
  Model_Search_Burnin = model_search_burn,
  Parameter_Thin = parameter_thin,
  Model_Thin = model_thin_cfg,
  Frame_Cache_Limit = frame_cache_limit,
  Representative_State_Limit = representative_state_limit,
  Population_Offset = use_population_offset,
  Spatiotemporal_Tensor = include_spatiotemporal_tensor,
  Candidate_Linear_SD = candidate_linear_sd,
  Candidate_Smooth_Null_SD = candidate_smooth_null_sd,
  Candidate_Smooth_Penalized_SD = candidate_smooth_penalized_sd,
  Posterior_Iterations = posterior_iter,
  Posterior_Burnin = posterior_burn,
  Posterior_Thin = posterior_thin_cfg,
  Posterior_Chains = posterior_chains_cfg,
  Posterior_Concurrent_Chains = min(posterior_cores_cfg, posterior_chains_cfg),
  Chain_Seed = chain_seed_cfg,
  Posterior_Families = length(fitted_models), Diagnostic_Seed = diagnostic_seed,
  Watson_MC = watson_mc_cfg, tb2_mc_cfg = tb2_mc_cfg, stringsAsFactors = FALSE
)
write.csv(
  run_manifest, file.path(dir_out, "run_manifest.csv"), row.names = FALSE
)

model_space_within_reference <- vapply(run_families, function(fam) {
  !is.null(search_results[[fam]]) &&
    isTRUE(
      search_results[[fam]]$diagnostics$allocation_diagnostics_within_reference
    )
}, logical(1))

audit_list <- lapply(run_families, function(fam) {
  a <- search_results[[fam]]$diagnostics$selection_audit
  if (is.null(a) && !is.null(search_results[[fam]])) {
    a <- selection_audit(search_results[[fam]])
  }
  if (is.null(a)) return(NULL)
  cbind(Family = fam, a, stringsAsFactors = FALSE)
})
selection_audit_all <- dplyr::bind_rows(audit_list)
if (nrow(selection_audit_all)) {
  write.csv(
    selection_audit_all, file.path(dir_out, "selection_audit.csv"),
    row.names = FALSE
  )
  write.csv(
    selection_concordance(selection_audit_all),
    file.path(dir_out, "selection_concordance.csv"), row.names = FALSE
  )
}
agree_list <- lapply(run_families, function(fam) {
  a <- search_results[[fam]]$diagnostics$replicate_agreement
  if (is.null(a) || !nrow(a)) return(NULL)
  cbind(Family = fam, a, stringsAsFactors = FALSE)
})
replicate_agreement_all <- dplyr::bind_rows(agree_list)
if (nrow(replicate_agreement_all)) {
  write.csv(
    replicate_agreement_all, file.path(dir_out, "selection_replicates.csv"),
    row.names = FALSE
  )
}

model_space_diagnostics <- data.frame(
  Family = run_families,
  Diagnostics_within_reference = unname(model_space_within_reference),
  stringsAsFactors = FALSE
)
write.csv(
  model_space_diagnostics, file.path(dir_out, "model_space_diagnostics.csv"),
  row.names = FALSE
)

if (length(model_space_within_reference) &&
    any(!model_space_within_reference)) {
  warning(
    "Model-space diagnostics outside reference values for: ",
    paste(
      names(model_space_within_reference)[!model_space_within_reference],
      collapse = ", "
    ),
    ". Treat their selections and posterior summaries as provisional.",
    call. = FALSE
  )
}

if (!nzchar(Sys.getenv("RBST_SKIP_POSTERIOR")) &&
    !identical(sort(names(fitted_models)), sort(run_families))) {
  stop(
    "Not all requested families produced post-selection posterior samples.",
    call. = FALSE
  )
}

sampler_diagnostics_summary <-
  dplyr::bind_rows(lapply(names(fitted_models), function(fam) {
    d <- attr(fitted_models[[fam]], "sampler_diagnostics")
    
    data.frame(
      Family = fam,
      Diagnostics_within_reference = !is.null(d) &&
        isTRUE(d$diagnostics_within_reference),
      Chains = if (is.null(d)) NA_integer_ else d$n_chains,
      Minimum_draws_per_chain = if (is.null(d)) {
        NA_integer_
      } else {
        d$min_draws_per_chain
      },
      Minimum_ESS_observed = if (is.null(d)) {
        NA_real_
      } else {
        d$minimum_ess_observed
      },
      ESS_reference = if (is.null(d)) NA_real_ else d$ess_reference,
      Maximum_split_drift = if (is.null(d)) {
        NA_real_
      } else {
        d$max_split_drift_observed
      },
      Split_drift_reference = if (is.null(d)) {
        NA_real_
      } else {
        d$split_drift_reference
      },
      Maximum_split_Rhat = if (is.null(d)) NA_real_ else d$max_rhat_observed,
      Rhat_reference = if (is.null(d)) NA_real_ else d$rhat_reference,
      Diagnostics_outside_reference = paste(
        d$diagnostics_outside_reference, collapse = "; "
      ),
      stringsAsFactors = FALSE
    )
  }))
write.csv(
  sampler_diagnostics_summary, file.path(dir_out, "sampler_diagnostics.csv"),
  row.names = FALSE
)

coefficient_diagnostics <-
  dplyr::bind_rows(lapply(names(fitted_models), function(fam) {
    d <- attr(fitted_models[[fam]], "sampler_diagnostics")
    if (is.null(d) || is.null(d$coefficient_diagnostics)) return(NULL)
    
    z <- d$coefficient_diagnostics
    z$Family <- fam
    z[, c("Family", setdiff(names(z), "Family")), drop = FALSE]
  }))

if (nrow(coefficient_diagnostics) > 0L) {
  write.csv(
    coefficient_diagnostics, file.path(dir_out, "coefficient_diagnostics.csv"),
    row.names = FALSE
  )
}

# 9. Goodness-of-fit diagnostics and comparisons (in-sample) ------------------
results_list <- list()
for (name in names(fitted_models)) {
  if (!is.null(fitted_models[[name]])) {
    results_list[[name]] <- compute_all_metrics(
      model_obj = fitted_models[[name]],
      label = name,
      data = df_model,
      y_col = "fatalities",
      space_col = "district_id",
      listw_obj = listw_obj,
      residual_seed = diagnostic_seed,
      watson_mc = watson_mc_cfg,
      tb2_mc = tb2_mc_cfg,
      tb2_m = tb2_m_cfg
    )
  }
}

if (length(results_list) > 0) {
  model_metrics <- do.call(rbind, results_list) %>% arrange(WAIC)
  
  model_diagnostics_table <- model_metrics %>%
    rename(
      `Distribution Family` = Model,
      `WAIC` = WAIC,
      `DIC` = DIC,
      `Log-Likelihood` = LogLik,
      `Finite RQR (%)` = RQR_Finite_Pct,
      `Extreme RQR (|r|>3) (%)` = Extreme_RQR_Pct,
      `Moran's I` = Moran_I,
      `Moran's p-value` = Moran_P,
      `KS statistic` = KS,
      `KS p-value` = KS_P,
      `Anderson-Darling statistic` = AD,
      `Anderson-Darling p-value` = AD_P,
      `Watson U2` = Watson,
      `Watson p-value` = Watson_P,
      `Neyman smooth statistic` = NS,
      `Neyman smooth p-value` = NS_P,
      `Zamanzade TB2` = TB2,
      `TB2 p-value` = TB2_P,
      `TB2 spacing m` = TB2_m,
      `Quadratic Score` = QS,
      `Logarithmic Score` = LS,
      `Spherical Score` = SPS,
      `CRPS` = CRPS,
      `RPS (K = max y)` = RPS_K,
      `Max omitted tail mass` = TailMass,
      `Max |PMF sum - CDF|` = PmfCdfGap
    )
  
  print(model_diagnostics_table)
  write.csv(
    model_diagnostics_table, file.path(dir_out, "model_fit_diagnostics.csv"),
    row.names = FALSE
  )
  
  metrics_gt <- gt::gt(model_diagnostics_table) |>
    gt::opt_table_font(font = "Times New Roman") |>
    gt::tab_options(table.font.size = 8) |>
    gt::fmt_number(decimals = 3)
  
  gt::gtsave(
    metrics_gt,
    filename = file.path(dir_out, "Table_Model_Fit_Diagnostics.docx")
  )
}

# Exceedance Brier scores and central 90% interval scores
threshold_scores <- lapply(names(fitted_models), function(name) {
  if (is.null(fitted_models[[name]])) return(NULL)
  tryCatch(
    compute_threshold_scores(fitted_models[[name]], df_model, "fatalities"),
    error = function(e) NULL
  )
})
names(threshold_scores) <- names(fitted_models)
for (part in c("exceedance", "interval")) {
  tab <- do.call(rbind, lapply(names(threshold_scores), function(name) {
    z <- threshold_scores[[name]][[part]]
    if (is.null(z)) NULL else cbind(Family = name, z)
  }))
  if (!is.null(tab)) {
    write.csv(
      tab, file.path(dir_out, paste0("predictive_", part, "_scores.csv")),
      row.names = FALSE
    )
  }
}

# Functional-level R-hat and ESS (location parameter, exceedance events)
functional_diag <- lapply(names(fitted_models), function(name) {
  if (is.null(fitted_models[[name]])) return(NULL)
  tryCatch(
    compute_functional_diagnostics(
      fitted_models[[name]], df_model, "fatalities"
    ),
    error = function(e) {
      message("Functional diagnostics failed for ", name, ": ",
              conditionMessage(e))
      NULL
    }
  )
})
names(functional_diag) <- names(fitted_models)
for (part in c("scalar", "curve", "location")) {
  tab <- do.call(rbind, lapply(names(functional_diag), function(name) {
    z <- functional_diag[[name]][[part]]
    if (is.null(z)) NULL else cbind(Family = name, z)
  }))
  if (!is.null(tab)) {
    write.csv(
      tab, file.path(dir_out, paste0("functional_diagnostics_", part, ".csv")),
      row.names = FALSE
    )
  }
}

# 10. Diagnostics and spatio-temporal outputs (all models) --------------------
pit_store <- list()
for (model_name in names(fitted_models)) {
  if (is.null(fitted_models[[model_name]])) next
  
  message(sprintf("\n>>> Processing outputs for model: %s", model_name))
  current_model_obj <- fitted_models[[model_name]]
  current_covariates <- selected_vars_list[[model_name]]
  current_selection_types <- selection_types_list[[model_name]]
  
  model_dir <- file.path(dir_mod, model_name)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_diag <- plot_rqr_diagnostics(
    model = current_model_obj,
    data = df_model,
    y_col = "fatalities",
    title = paste("Diagnostics:", model_name),
    seed = diagnostic_seed
  )
  ggsave(
    file.path(model_dir, "Figure_Quantile_Residuals.pdf"), p_diag, width = 14,
    height = 4.5, device = cairo_pdf
  )
  
  pit_res <- get_pit_plots(
    model = current_model_obj,
    data = df_model,
    y_col = "fatalities",
    label = paste("PIT:", model_name),
    seed = diagnostic_seed
  )
  ggsave(
    file.path(model_dir, "Figure_PIT_Histogram.pdf"), pit_res$hist, width = 6,
    height = 5, device = cairo_pdf
  )
  ggsave(
    file.path(model_dir, "Figure_PIT_ECDF.pdf"), pit_res$ecdf, width = 6,
    height = 5, device = cairo_pdf
  )
  
  rqr_map <- quantile_residuals(
    model = current_model_obj,
    data = df_model,
    y_col = "fatalities",
    seed = diagnostic_seed
  )
  
  pit_store[[model_name]] <- stats::pnorm(rqr_map[is.finite(rqr_map)])
  tryCatch(
    ggsave(
      file.path(model_dir, "Figure_PIT_Panels.pdf"),
      plot_pit_panels(pit_store[[model_name]], model_name),
      width = 14, height = 3.6, device = cairo_pdf
    ),
    error = function(e) message("PIT panels skipped: ", conditionMessage(e))
  )
  
  p_sp_resid <- plot_spatial_risk(
    data = df_model %>% mutate(Resid = rqr_map),
    shape = districts_sf,
    space_col = "district_id",
    y_col = "Resid",
    shape_id = "district_id",
    label_col = "canonical_name",
    aggregation = "mean"
  ) +
    geom_sf(data = districts_sf, fill = NA, color = "gray95") +
    labs(title = paste("Spatial residual map:", model_name), subtitle = "")
  ggsave(
    file.path(model_dir, "Figure_Spatial_Residuals.pdf"), p_sp_resid, width = 8,
    height = 8, device = cairo_pdf
  )
  
  tryCatch(
    save_trace_pdf(current_model_obj,
                   file.path(model_dir, "Figure_Trace_Plots.pdf")),
    error = function(e) message("Trace plots skipped: ", conditionMessage(e))
  )
  
  fixed_res <- plot_fixed_effects(current_model_obj)
  ggsave(
    file.path(model_dir, "Figure_Fixed_Effects.pdf"), fixed_res$plot,
    width = 12, height = 7, device = cairo_pdf
  )
  
  smooth_plots_by_param <- list()
  for (param_key in names(current_covariates)) {
    selected_in_param <- current_covariates[[param_key]]
    active_smooths_param <- list()
    for (v in selected_in_param) {
      type_v <- if (!is.null(current_selection_types[[param_key]]) &&
                    v %in% names(current_selection_types[[param_key]])) {
        unname(current_selection_types[[param_key]][v])
      } else { NA_character_ }
      
      if (identical(type_v, "nonlinear") &&
          v %in% names(df_model) &&
          is.numeric(df_model[[v]]) &&
          length(unique(df_model[[v]][is.finite(df_model[[v]])])) > 7) {
        active_smooths_param[[sprintf("s(%s)", v)]] <- v
      }
    }
    if (length(active_smooths_param)) {
      subtitle <- bquote("Distributional parameter:" ~ .(as.name(param_key)))
      gp <- tryCatch(
        plot_smooth_effects(
          current_model_obj, df_model, active_smooths_param,
          parameter = param_key
        ) + ggplot2::labs(subtitle = subtitle),
        error = function(e) {
          message(
            sprintf(
              "  [Output warning] Smooth effects for %s/%s failed: %s",
              model_name, param_key, conditionMessage(e)
            )
          )
          NULL
        }
      )
      if (!is.null(gp)) smooth_plots_by_param[[param_key]] <- gp
    }
  }
  if (length(smooth_plots_by_param)) {
    ggsave(
      file.path(model_dir, "Figure_NonLinear_Smooth_Effects.pdf"),
      patchwork::wrap_plots(smooth_plots_by_param, ncol = 1),
      width = 14, height = max(5, 4.5 * length(smooth_plots_by_param)),
      device = cairo_pdf
    )
  }
  
  plot_func <- if (model_name == "RSB") "location" else "mean"
  
  if (!identical(model_name, "RSB")) {
    fitted_mean <- compute_expected_value(current_model_obj)
    
    if (
      length(fitted_mean) != nrow(df_model) ||
      any(!is.finite(fitted_mean)) ||
      any(fitted_mean < 0)
    ) {
      stop(
        "Fitted means are not finite, non-negative and aligned for ",
        model_name,
        ".",
        call. = FALSE
      )
    }
    
    prediction_audit <- df_model %>%
      transmute(
        analysis_row_id = dplyr::row_number(),
        district_id = as.character(district_id),
        canonical_name,
        year_month,
        year = lubridate::year(year_month),
        observed_fatalities = fatalities,
        fitted_mean = fitted_mean
      )
    
    readr::write_csv(
      prediction_audit,
      file.path(model_dir, "fitted_mean_by_district_month.csv"),
      na = ""
    )
    
    annual_calibration <- prediction_audit %>%
      group_by(district_id, canonical_name, year) %>%
      summarise(
        observed_fatalities = sum(observed_fatalities),
        fitted_fatalities = sum(fitted_mean),
        .groups = "drop"
      ) %>%
      mutate(
        prediction_error = fitted_fatalities - observed_fatalities,
        observed_to_fitted_ratio = dplyr::if_else(
          fitted_fatalities > 0,
          observed_fatalities / fitted_fatalities,
          NA_real_
        )
      )
    
    readr::write_csv(
      annual_calibration,
      file.path(model_dir, "annual_observed_vs_fitted.csv"),
      na = ""
    )
    
    calibration_plot <- ggplot(
      annual_calibration,
      aes(x = observed_fatalities, y = fitted_fatalities)
    ) +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      geom_point(alpha = 0.75) +
      facet_wrap(~ year, scales = "free") +
      scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
      labs(
        x = "Observed annual fatalities",
        y = "Fitted annual fatalities",
        subtitle = paste("Distribution family:", model_name)
      ) +
      theme_bw(base_size = 10)
    
    ggsave(
      file.path(model_dir, "Figure_Observed_vs_Fitted_Annual.pdf"),
      calibration_plot,
      width = 12,
      height = 7,
      device = cairo_pdf
    )
  }
  
  # Use RBST_MAP_SCALE_MAX to lock the legend across independently fitted
  # families. The separate replot utility derives a shared maximum from
  # completed family CSVs, without repeating posterior sampling.
  map_scale_max <- Sys.getenv("RBST_MAP_SCALE_MAX", "")
  if (identical(plot_func, "location")) map_scale_max <- ""

  g_evolution <- plot_expected_spatiotemporal_risk(
    model = current_model_obj,
    data = df_model,
    shape = districts_sf,
    space_col = "district_id",
    time_col = "date_obj",
    shape_id = "district_id",
    label_col = "canonical_name",
    time_unit = "year",
    functional = plot_func,
    scale_max = if (nzchar(map_scale_max)) map_scale_max else NULL
  ) +
    geom_sf(data = districts_sf, fill = NA, color = "gray95")
  
  ggsave(
    file.path(model_dir, "Figure_SpatioTemporal_Evolution.pdf"), g_evolution,
    width = 14, height = 9, bg = "white", device = cairo_pdf
  )
  
  summary_tbl <- tryCatch(
    create_model_summary_table(current_model_obj, include_intercept = TRUE),
    error = function(e) {
      message(
        sprintf(
          "  [Output warning] Summary table for %s failed: %s", model_name,
          conditionMessage(e)
        )
      )
      NULL
    }
  )
  if (!is.null(summary_tbl)) {
    gt::gtsave(
      summary_tbl,
      file.path(
        model_dir, paste0("Table_Summary_Estimates_", model_name, ".docx")
      )
    )
  }
}

if (full_run && length(pit_store)) {
  tryCatch(
    ggsave(
      file.path(dir_out, "Figure_PIT_ECDF_Difference_AllFamilies.pdf"),
      plot_ecdf_difference_grid(pit_store), width = 12, height = 6.5,
      device = cairo_pdf
    ),
    error = function(e) message("PIT grid skipped: ", conditionMessage(e))
  )
}

# 11. Print and export model covariate allocations ----------------------------
formula_summary_list <- list()
for (m in names(model_formulas)) {
  cat(sprintf("\nModel: %s\n", m))
  for (param in names(model_formulas[[m]])) {
    cat(sprintf("  - %s: %s\n", param, model_formulas[[m]][[param]]))
    formula_summary_list[[length(formula_summary_list) + 1]] <- data.frame(
      Model = m,
      Parameter = param,
      Formula = model_formulas[[m]][[param]]
    )
  }
}

if (length(formula_summary_list)) {
  write.csv(
    do.call(rbind, formula_summary_list),
    file.path(dir_out, "model_formulas_summary.csv"),
    row.names = FALSE
  )
}

# 12. Result tables and summary figures ---------------------------------------
if (full_run) {
  tryCatch(
    write_result_tables(dir_mod),
    error = function(e) message("Result tables failed: ", conditionMessage(e))
  )
}

# 13. Out-of-sample validation over consecutive later periods ----------------
# RBST_OOS=1 refits each family on an initial window, scores the next blocks of
# months and writes the scores of the model-averaged predictive distribution.
if (nzchar(Sys.getenv("RBST_OOS"))) {
  oos_iter <- as.integer(Sys.getenv("OOS_ITER", "8000"))
  oos_burn <- as.integer(Sys.getenv("OOS_BURN", "3000"))
  oos_scores <- validate_prospective_bma(
    families = run_families, base_formula = base_mu_formula,
    candidate_pool = model_candidate_pool, data = df_model,
    auxiliary_formulas = auxiliary_formulas, offsets = model_offsets,
    candidate_admissibility = candidate_admissibility,
    time_col = "date_obj", n_blocks = 3L, min_train_frac = 0.60,
    n_iter = oos_iter, n_burn = oos_burn, thin = parameter_thin,
    c_linear = candidate_linear_sd, c0 = candidate_smooth_null_sd,
    cplus = candidate_smooth_penalized_sd,
    within_sweeps = 5L, a_V = 0.5, a_A = 0.2,
    cache_frames = TRUE, cache_blocks = TRUE,
    max_frame_cache = frame_cache_limit,
    max_representative_states = representative_state_limit
  )
  write.csv(
    oos_scores, file.path(dir_out, "prospective_validation.csv"),
    row.names = FALSE
  )
}
