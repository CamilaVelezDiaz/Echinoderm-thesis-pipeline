# ==============================================================================
# RAREFACTION AND RICHNESS ESTIMATION  (v2 - incidence / coverage-based)
# ==============================================================================
# Camila Velez | JCU Marine Biology MSc | Thesis Chapter - Analysis
#
# Question: given current sampling effort, how complete is the species
# inventory of each depth zone, and how does richness compare among zones
# once sampling completeness is standardised?
#
# WHY THIS VERSION (response to examiner comment on Methods pp. 16-17,
# Results pp. 32-33, Discussion p. 48, Appendix S p. 96):
#   v1 treated the number of occurrence RECORDS per species as "abundance".
#   Records reflect heterogeneous collecting, database and digitisation
#   effort, not comparable individuals, so abundance-based rarefaction to a
#   common number of records (and extrapolating Abyssal from 3 to 234
#   records) was not defensible as formal inference.
#
#   v2 changes:
#   1. SAMPLING UNIT = COLLECTING EVENT, defined as a unique location
#      (latitude/longitude rounded to 0.01 deg, ~1 km) on a unique date.
#      Each species is scored present/absent per event (incidence data), so
#      repeated records of the same species from one haul/dive/station -
#      including duplicates across databases - count once.
#   2. iNEXT is run with datatype = "incidence_freq"; Chao2 replaces Chao1.
#   3. Zones are compared at a COMMON SAMPLE COVERAGE (Chao & Jost 2012),
#      not a common number of records. The standard coverage is iNEXT's
#      default: the lowest coverage any zone reaches at DOUBLE its own
#      number of events, so no zone is extrapolated beyond 2x its sample.
#   4. Zones with fewer than MIN_EVENTS_FOR_INFERENCE events (Abyssal) are
#      NOT rarefied, extrapolated or tested. They are reported descriptively
#      as evidence of insufficient sampling.
#   5. The pairwise z-test / compact letters are removed. Zones are compared
#      by whether their 95% CIs at the common coverage overlap
#      (non-overlap = conservative evidence of a difference; Chao et al. 2014).
#   6. Sensitivity analysis: the species-level analysis is repeated with
#      0.5-degree grid cells as the sampling unit, which retains undated
#      records that the event definition has to drop.
#
#   Interpretation caveat (state in Methods): even event-based curves
#   describe the specimen/observation RECORD, which was not collected under
#   a standardised design. Estimates are measures of inventory completeness
#   and relative richness of the recorded fauna, not unbiased estimates of
#   true biological richness.
#
# Input:  echino_wide.csv
# Outputs (all in OUT_DIR):
#   sampling_summary_by_zone.csv            descriptive table, ALL zones
#   chao2_richness_by_zone.csv              observed, Chao2, completeness
#   coverage_standardised_richness.csv      richness at common coverage
#   coverage_ci_overlap.csv                 pairwise CI-overlap comparison
#   chao2_genus_richness_by_zone.csv        genus-level equivalents
#   coverage_standardised_richness_genus.csv
#   coverage_ci_overlap_genus.csv
#   sensitivity_gridcell_species.csv        0.5-deg grid-cell sensitivity
#   fig_rarefaction_events_species.png/.pdf richness vs collecting events
#   fig_coverage_species.png/.pdf           completeness + richness vs coverage
#   fig_rarefaction_events_genus.png        genus-level (supplementary)
# ==============================================================================

# =============================================================================
# SETUP
# =============================================================================
# setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(ggplot2)
library(iNEXT)
library(patchwork)

set.seed(42)
NBOOT <- 999                     # bootstrap replicates for CIs
MIN_EVENTS_FOR_INFERENCE <- 20   # zones with fewer events: descriptive only
COORD_DECIMALS <- 2              # event = location rounded to 0.01 deg + date
GRID_SIZE_DEG  <- 0.5            # sensitivity-analysis grid cell size
OUT_DIR <- "."

zone_order <- c("Continental Shelf", "Upper Slope", "Lower Slope", "Abyssal")
zone_colors <- c(
  "Continental Shelf" = "#1B9E77",
  "Upper Slope"       = "#D95F02",
  "Lower Slope"       = "#7570B3",
  "Abyssal"           = "#E7298A"
)
zone_labels_with_range <- c(
  "Continental Shelf" = "Continental Shelf\n(0 - 200 m)",
  "Upper Slope"       = "Upper Slope\n(200 - 1000 m)",
  "Lower Slope"       = "Lower Slope\n(1000 - 4000 m)",
  "Abyssal"           = "Abyssal\n(>4000 m)"
)

# =============================================================================
# LOAD DATA
# =============================================================================
echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE, guess_max = 50000)

if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c("Species"="Species","Subspecies"="Species","Genus"="Genus","Subgenus"="Genus",
                "Family"="Family","Order"="Order","Class"="Class","Phylum"="Phylum")
  echino_wide <- echino_wide %>% mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}
if (!"crosses_zone_boundary" %in% names(echino_wide)) {
  echino_wide <- echino_wide %>%
    mutate(crosses_zone_boundary = str_starts(coalesce(depth_zone, ""), "Spans"))
}
echino_wide <- echino_wide %>%
  mutate(crosses_zone_boundary = as.logical(toupper(as.character(crosses_zone_boundary))),
         event_date = na_if(substr(as.character(best_eventDate), 1, 10), ""))

# -----------------------------------------------------------------------------
# GEOGRAPHIC SCOPE SENSITIVITY (PNG + New Caledonia)
# Run once with "all" and once with "australia_only"; outputs go to separate
# folders (richness_all/, richness_australia_only/) so they can be compared.
# -----------------------------------------------------------------------------
RICH_SCOPE <- "all"   # <- "all" or "australia_only"

echino_wide <- echino_wide %>%
  mutate(
    loc_text   = str_to_lower(paste(verbatimLocality, `Locality Name`, `Locality Remarks`)),
    outside_au = coalesce(country != "Australia", FALSE) |
      str_detect(loc_text, "new caledonia|nouvelle.cal|noum[eé]a")
  )
cat(sprintf("Records flagged outside Australia: %d\n", sum(echino_wide$outside_au)))
if (RICH_SCOPE == "australia_only") echino_wide <- filter(echino_wide, !outside_au)
cat(sprintf("Richness scope: %s | %d records used\n", RICH_SCOPE, nrow(echino_wide)))

OUT_DIR <- paste0("richness_", RICH_SCOPE)
dir.create(OUT_DIR, showWarnings = FALSE)

# =============================================================================
# BUILD ANALYSIS TABLES
# =============================================================================
# taxon_level = "species" uses species-level records; "genus" uses species-
# and genus-level records collapsed to genus.
prep_records <- function(level = c("species", "genus")) {
  level <- match.arg(level)
  ranks <- if (level == "species") "Species" else c("Species", "Genus")
  echino_wide %>%
    filter(rank6 %in% ranks, !is.na(accepted_name_final),
           !coalesce(crosses_zone_boundary, FALSE), depth_zone %in% zone_order) %>%
    mutate(
      taxon = if (level == "species") accepted_name_final else word(accepted_name_final, 1),
      has_coords = !is.na(best_latitude) & !is.na(best_longitude),
      event_id = ifelse(has_coords & !is.na(event_date),
                        paste(round(best_latitude, COORD_DECIMALS),
                              round(best_longitude, COORD_DECIMALS), event_date, sep = "_"),
                        NA_character_),
      cell_id  = ifelse(has_coords,
                        paste(floor(best_latitude / GRID_SIZE_DEG) * GRID_SIZE_DEG,
                              floor(best_longitude / GRID_SIZE_DEG) * GRID_SIZE_DEG, sep = "_"),
                        NA_character_)
    ) %>%
    filter(!is.na(taxon), taxon != "")
}

# iNEXT incidence_freq format: c(number of sampling units, incidence freq of each taxon)
make_incidence <- function(d, unit_col) {
  d <- d[!is.na(d[[unit_col]]), ]
  n_units <- n_distinct(d[[unit_col]])
  freqs <- d %>% distinct(.data[[unit_col]], taxon) %>% count(taxon) %>% pull(n)
  c(n_units, sort(freqs, decreasing = TRUE))
}

incidence_by_zone <- function(d, unit_col, zones) {
  out <- lapply(zones, function(z) make_incidence(filter(d, depth_zone == z), unit_col))
  setNames(out, zones)
}

# Descriptive summary for EVERY zone (this is where Abyssal is reported)
sampling_summary <- function(d, level_label) {
  ev <- d %>% filter(!is.na(event_id))
  inc <- ev %>% distinct(depth_zone, event_id, taxon) %>% count(depth_zone, taxon, name = "n_events")
  inc_sum <- inc %>% group_by(depth_zone) %>%
    summarise(uniques_Q1 = sum(n_events == 1), duplicates_Q2 = sum(n_events == 2), .groups = "drop")
  d %>% group_by(depth_zone) %>%
    summarise(
      n_records          = n(),
      n_taxa_observed    = n_distinct(taxon),
      records_undated_or_no_coords = sum(is.na(event_id)),
      n_collecting_events = n_distinct(event_id, na.rm = TRUE),
      n_grid_cells_0.5deg = n_distinct(cell_id, na.rm = TRUE),
      .groups = "drop") %>%
    left_join(inc_sum, by = "depth_zone") %>%
    mutate(level = level_label,
           pct_taxa_unique = round(100 * uniques_Q1 / n_taxa_observed, 1),
           analysis = ifelse(n_collecting_events >= MIN_EVENTS_FOR_INFERENCE,
                             "rarefaction / coverage-standardised",
                             "descriptive only (insufficient sampling)"),
           depth_zone = factor(depth_zone, levels = zone_order)) %>%
    arrange(depth_zone) %>%
    relocate(level)
}

# Pairwise CI-overlap comparison (replaces z-test + compact letters)
ci_overlap_table <- function(est) {
  z <- as.character(est$depth_zone)
  pr <- combn(seq_along(z), 2)
  bind_rows(lapply(seq_len(ncol(pr)), function(k) {
    a <- pr[1, k]; b <- pr[2, k]
    overlap <- est$ci_lower[a] <= est$ci_upper[b] && est$ci_lower[b] <= est$ci_upper[a]
    tibble(zone_a = z[a], zone_b = z[b],
           richness_a = est$richness[a], richness_b = est$richness[b],
           ci_a = sprintf("%.1f-%.1f", est$ci_lower[a], est$ci_upper[a]),
           ci_b = sprintf("%.1f-%.1f", est$ci_lower[b], est$ci_upper[b]),
           ci_overlap = overlap,
           interpretation = ifelse(overlap, "no clear difference",
                                   "non-overlapping 95% CIs (difference)"))
  }))
}

# Coverage-standardised estimates at iNEXT's default common coverage
coverage_standardise <- function(inc_list) {
  w <- character(0)
  est <- withCallingHandlers(
    estimateD(inc_list, q = 0, datatype = "incidence_freq", base = "coverage",
              level = NULL, nboot = NBOOT, conf = 0.95),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  if (length(w)) cat("  estimateD warnings:", unique(w), sep = "\n   ")
  est %>%
    transmute(depth_zone = factor(Assemblage, levels = zone_order),
              standard_coverage = round(SC, 3),
              events_needed = round(t, 0),
              method = Method,
              richness = round(qD, 1),
              ci_lower = round(qD.LCL, 1),
              ci_upper = round(qD.UCL, 1)) %>%
    arrange(depth_zone)
}

asymptotic_table <- function(inext_obj) {
  info <- inext_obj$DataInfo %>% transmute(depth_zone = Assemblage, n_events = T,
                                           sample_coverage = round(SC, 3))
  inext_obj$AsyEst %>%
    filter(Diversity == "Species richness") %>%
    transmute(depth_zone = Assemblage, observed = Observed,
              chao2 = round(Estimator, 1), chao2_se = round(s.e., 1),
              chao2_ci_lower = round(LCL, 1), chao2_ci_upper = round(UCL, 1)) %>%
    left_join(info, by = "depth_zone") %>%
    mutate(pct_of_chao2_observed = round(100 * observed / chao2, 1),
           depth_zone = factor(depth_zone, levels = zone_order)) %>%
    relocate(n_events, sample_coverage, .after = depth_zone) %>%
    arrange(depth_zone)
}

# Run iNEXT separately for each zone with LOG-SPACED knots from 1 to 2x that
# zone's own number of events (plus the observed point), then merge into one
# iNEXT object. With a single shared knot set, the Continental Shelf (thousands
# of events) gets almost no knots in its steep early section, which shows up as
# a jagged coverage curve; per-zone knots fix that without changing estimates.
run_inext <- function(inc_list, nknots = 80) {
  parts <- lapply(names(inc_list), function(z) {
    T_obs <- inc_list[[z]][1]
    sizes <- sort(unique(c(round(exp(seq(0, log(2 * T_obs), length.out = nknots))), T_obs)))
    iNEXT(inc_list[z], q = 0, datatype = "incidence_freq", size = sizes,
          se = TRUE, conf = 0.95, nboot = NBOOT)
  })
  out <- list(
    DataInfo = bind_rows(lapply(parts, `[[`, "DataInfo")),
    iNextEst = list(
      size_based     = bind_rows(lapply(parts, function(o) o$iNextEst$size_based)),
      coverage_based = bind_rows(lapply(parts, function(o) o$iNextEst$coverage_based))),
    AsyEst = bind_rows(lapply(parts, `[[`, "AsyEst")))
  class(out) <- "iNEXT"
  out
}

subset_inext <- function(obj, zones) {
  obj$DataInfo <- filter(obj$DataInfo, Assemblage %in% zones)
  obj$iNextEst$size_based     <- filter(obj$iNextEst$size_based, Assemblage %in% zones)
  obj$iNextEst$coverage_based <- filter(obj$iNextEst$coverage_based, Assemblage %in% zones)
  obj$AsyEst <- filter(obj$AsyEst, Assemblage %in% zones)
  obj
}

style_zone_scales <- function(p, zones) {
  suppressMessages(p +
                     scale_colour_manual(values = zone_colors, limits = zones, labels = zone_labels_with_range[zones]) +
                     scale_fill_manual(values = zone_colors, limits = zones, labels = zone_labels_with_range[zones]) +
                     scale_shape_discrete(limits = zones, labels = zone_labels_with_range[zones]))
}

# =============================================================================
# DESCRIPTIVE SAMPLING SUMMARY (all zones, species and genus level)
# =============================================================================
sp <- prep_records("species")
gn <- prep_records("genus")

summary_tbl <- bind_rows(sampling_summary(sp, "species"), sampling_summary(gn, "genus"))
cat("\n=== Sampling summary by depth zone (all zones) ===\n")
print(as.data.frame(summary_tbl))
write_csv(summary_tbl, file.path(OUT_DIR, "sampling_summary_by_zone.csv"))

zones_inf_sp <- summary_tbl %>% filter(level == "species", n_collecting_events >= MIN_EVENTS_FOR_INFERENCE) %>%
  pull(depth_zone) %>% as.character()
zones_inf_gn <- summary_tbl %>% filter(level == "genus", n_collecting_events >= MIN_EVENTS_FOR_INFERENCE) %>%
  pull(depth_zone) %>% as.character()
cat("\nZones used for rarefaction (species):", paste(zones_inf_sp, collapse = ", "), "\n")
cat("Zones reported descriptively only (species):",
    paste(setdiff(zone_order, zones_inf_sp), collapse = ", "), "\n")

# =============================================================================
# PART 1: SPECIES-LEVEL INCIDENCE RAREFACTION / EXTRAPOLATION + CHAO2
# =============================================================================
inc_sp <- incidence_by_zone(sp, "event_id", zones_inf_sp)
cat("\n=== Running iNEXT (species, incidence, nboot =", NBOOT, ") ===\n")
inext_sp <- run_inext(inc_sp)

chao2_sp <- asymptotic_table(inext_sp)
cat("\n-- Chao2 asymptotic species richness (collecting events as sampling units) --\n")
print(as.data.frame(chao2_sp))
write_csv(chao2_sp, file.path(OUT_DIR, "chao2_richness_by_zone.csv"))

# Figure 1: richness vs number of collecting events, with inset for slope zones
p_events <- style_zone_scales(ggiNEXT(inext_sp, type = 1, se = TRUE), zones_inf_sp) +
  labs(x = "Number of collecting events", y = "Species richness",
       title = "Incidence-based rarefaction and extrapolation by depth zone",
       subtitle = "Solid = rarefaction; dashed = extrapolation (to 2x observed events); shading = 95% CI") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

small_sp <- setdiff(zones_inf_sp, "Continental Shelf")
if (length(small_sp) >= 1) {
  inext_sp_small <- subset_inext(inext_sp, small_sp)
  p_inset <- style_zone_scales(ggiNEXT(inext_sp_small, type = 1, se = TRUE), small_sp) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none", plot.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = "grey40"),
          plot.margin = margin(4, 8, 4, 4)) +
    labs(x = NULL, y = NULL)
  p_events <- p_events + inset_element(p_inset, left = 0.50, bottom = 0.06, right = 0.97, top = 0.50)
}
ggsave(file.path(OUT_DIR, "fig_rarefaction_events_species.png"), p_events, width = 9, height = 6, dpi = 300)
ggsave(file.path(OUT_DIR, "fig_rarefaction_events_species.pdf"), p_events, width = 9, height = 6,
       device = cairo_pdf)

# Figure 2: (a) sample completeness vs events, (b) richness vs sample coverage
p_cov_a <- style_zone_scales(ggiNEXT(inext_sp, type = 2, se = TRUE), zones_inf_sp) +
  scale_x_log10() +
  labs(x = "Number of collecting events (log scale)", y = "Sample coverage", title = "(a) Sample completeness") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
p_cov_b <- style_zone_scales(ggiNEXT(inext_sp, type = 3, se = TRUE), zones_inf_sp) +
  labs(x = "Sample coverage", y = "Species richness", title = "(b) Coverage-based rarefaction") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
p_cov <- (p_cov_a | p_cov_b) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "fig_coverage_species.png"), p_cov, width = 11, height = 5.5, dpi = 300)
ggsave(file.path(OUT_DIR, "fig_coverage_species.pdf"), p_cov, width = 11, height = 5.5,
       device = cairo_pdf)

# =============================================================================
# PART 2: COVERAGE-STANDARDISED RICHNESS + CI-OVERLAP COMPARISON
# =============================================================================
cat("\n=== Coverage-standardised species richness ===\n")
cov_sp <- coverage_standardise(inc_sp)
print(as.data.frame(cov_sp))
write_csv(cov_sp, file.path(OUT_DIR, "coverage_standardised_richness.csv"))

ovl_sp <- ci_overlap_table(cov_sp)
cat("\n-- Pairwise comparison by 95% CI overlap (species) --\n")
print(as.data.frame(ovl_sp))
write_csv(ovl_sp, file.path(OUT_DIR, "coverage_ci_overlap.csv"))

# =============================================================================
# PART 3 (SUPPLEMENTARY): GENUS LEVEL
# =============================================================================
inc_gn <- incidence_by_zone(gn, "event_id", zones_inf_gn)
cat("\n=== Running iNEXT (genus, incidence) ===\n")
inext_gn <- run_inext(inc_gn)
chao2_gn <- asymptotic_table(inext_gn)
print(as.data.frame(chao2_gn))
write_csv(chao2_gn, file.path(OUT_DIR, "chao2_genus_richness_by_zone.csv"))

cov_gn <- coverage_standardise(inc_gn)
print(as.data.frame(cov_gn))
write_csv(cov_gn, file.path(OUT_DIR, "coverage_standardised_richness_genus.csv"))
ovl_gn <- ci_overlap_table(cov_gn)
print(as.data.frame(ovl_gn))
write_csv(ovl_gn, file.path(OUT_DIR, "coverage_ci_overlap_genus.csv"))

p_gn <- style_zone_scales(ggiNEXT(inext_gn, type = 1, se = TRUE), zones_inf_gn) +
  labs(x = "Number of collecting events", y = "Genus richness",
       title = "Genus-level incidence-based rarefaction (supplementary)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "fig_rarefaction_events_genus.png"), p_gn, width = 9, height = 6, dpi = 300)

# =============================================================================
# PART 4: SENSITIVITY - 0.5-DEGREE GRID CELLS AS SAMPLING UNIT (species)
# =============================================================================
zones_cell <- summary_tbl %>% filter(level == "species", n_grid_cells_0.5deg >= MIN_EVENTS_FOR_INFERENCE) %>%
  pull(depth_zone) %>% as.character()
cat("\n=== Sensitivity: 0.5-degree grid cells (zones:", paste(zones_cell, collapse = ", "), ") ===\n")
inc_cell <- incidence_by_zone(sp, "cell_id", zones_cell)
inext_cell <- run_inext(inc_cell)
sens <- asymptotic_table(inext_cell) %>%
  rename(n_cells = n_events) %>%
  left_join(coverage_standardise(inc_cell) %>%
              select(depth_zone, standard_coverage, cov_richness = richness,
                     cov_ci_lower = ci_lower, cov_ci_upper = ci_upper, cov_method = method),
            by = "depth_zone")
print(as.data.frame(sens))
write_csv(sens, file.path(OUT_DIR, "sensitivity_gridcell_species.csv"))

cat("\n=== Done ===\n")