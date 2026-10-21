# =============================================================================
# SPECIMEN PROVENANCE ("BRAIN DRAIN") ANALYSIS
# =============================================================================
# Addresses: of the specimens COLLECTED in the NE Australia/QLD study
# region, where are they now physically HELD? A QLD museum, an interstate
# Australian museum, or an overseas institution?
#
# IMPORTANT SCOPE DECISION: this analysis is restricted to records with a
# genuine physical specimen (obs_quality_flag 1 - preserved specimens only), 
# EXCLUDING Flag 2: Human observation and citizen-science/observation-only
# records (Flag 3: iNaturalist, Diveboard, etc.) and known non-physical
# aggregator/platform entries (Atlas of Living Australia as a platform
# name, rather than a specimen-holding institution). "Where is this
# specimen held" is a meaningless question for a photo observation with
# no voucher - including these would inflate the apparent number of
# "institutions" without representing genuine specimen custody.
#
# CASE-DUPLICATE FIX: institutionCode values include case-inconsistent
# duplicates (QM/qm, AM/am, MCZ/mcz) from different source datasets -
# standardised to uppercase before classification, mirroring the same
# case-standardisation already applied to catalogNumber/occurrenceID
# elsewhere in this pipeline.
#
# CLASSIFICATION CONFIDENCE: every institution code below is tagged with
# a confidence level. HIGH = well-known major museum/institution, location
# certain. LOW/UNVERIFIED = code could not be confidently identified from
# available information - flagged explicitly rather than guessed, and
# these records are EXCLUDED from the location breakdown until manually
# verified (their raw counts are still reported separately).
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================
setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")
# =============================================================================
# Library
# =============================================================================
library(tidyverse)

echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)

# Defensive check (matching every other script in this pipeline) - rank6
# is needed for Analysis 1 below (species-level identification rate by
# holding location), even though it wasn't needed for the earlier parts
# of this script.
if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c(
    "Species" = "Species", "Subspecies" = "Species",
    "Genus"   = "Genus",   "Subgenus"   = "Genus",
    "Family"  = "Family", "Order" = "Order",
    "Class"   = "Class",  "Phylum" = "Phylum"
  )
  echino_wide <- echino_wide %>%
    mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}

inst_cols <- names(echino_wide)[str_detect(names(echino_wide), "^institutionCode")]

echino_wide <- echino_wide %>%
  mutate(
    institution_raw = coalesce(!!!syms(inst_cols)),
    institution_std = str_to_upper(str_trim(institution_raw))
  )

# -----------------------------------------------------------------------------
# CLASSIFICATION LOOKUP
# -----------------------------------------------------------------------------
# location_category: "QLD" | "Other Australian state/territory" | "International"
# confidence: "high" | "low_unverified"
# state_or_country: specific state (Australia) or country (international),
#   for a more granular breakdown if useful
# -----------------------------------------------------------------------------
institution_lookup <- tribble(
  ~institution_std,                                    ~full_name,                                              ~location_category,                  ~state_or_country,        ~confidence,
  # --- Queensland ---
  "QM",                                                 "Queensland Museum",                                     "QLD",                                "Queensland",              "high",
  "MTQ",                                                 "Museum of Tropical Queensland",                         "QLD",                                "Queensland",              "high",
  "AIMS",                                                "Australian Institute of Marine Science",                "QLD",                                "Queensland",              "high",
  "UQ",                                                   "University of Queensland",                              "QLD",                                "Queensland",              "high",
  "TORRES STRAIT CRC",                                    "Torres Strait Cooperative Research Centre",             "QLD",                                "Queensland",              "high",
  # --- Other Australian states/territories ---
  "AM",                                                   "Australian Museum",                                     "Other Australian state/territory",  "New South Wales",        "high",
  "CSIRO, AUSTRALIA",                                     "CSIRO National Research Collections Australia",         "Other Australian state/territory",  "Tasmania",                "high",
  "NMV",                                                  "Museums Victoria",                                      "Other Australian state/territory",  "Victoria",                "high",
  "MUSEUMS VICTORIA, AUSTRALIA",                           "Museums Victoria",                                      "Other Australian state/territory",  "Victoria",                "high",
  "MAGNT",                                                 "Museum and Art Gallery of the Northern Territory",      "Other Australian state/territory",  "Northern Territory",     "high",
  "WAM",                                                   "Western Australian Museum",                             "Other Australian state/territory",  "Western Australia",      "high",
  "SAMA",                                                  "South Australian Museum",                               "Other Australian state/territory",  "South Australia",        "high",
  "TMAG",                                                  "Tasmanian Museum and Art Gallery",                      "Other Australian state/territory",  "Tasmania",                "high",
  "INSTITUTE FOR MARINE AND ANTARCTIC STUDIES (IMAS)",     "Institute for Marine and Antarctic Studies (UTAS)",     "Other Australian state/territory",  "Tasmania",                "high",
  "DEW",                                                   "Department for Environment and Water (SA)",             "Other Australian state/territory",  "South Australia",        "high",
  "BRS",                                                   "Bureau of Rural Sciences - Commercial Fisheries Presence, Australia, 2000-2002", "Other Australian state/territory", "ACT (federal)", "high",
  # --- International ---
  "NORWEGIAN INSTITUTE FOR WATER RESEARCH",                "Norwegian Institute for Water Research (NIVA)",         "International",                     "Norway",                  "high",
  "UF",                                                    "Florida Museum of Natural History (Univ. of Florida)",  "International",                     "United States",           "high",
  "MCZ",                                                   "Museum of Comparative Zoology, Harvard University",     "International",                     "United States",           "high",
  "USNM",                                                  "Smithsonian National Museum of Natural History",       "International",                     "United States",           "high",
  "NHMUK",                                                 "Natural History Museum, London",                        "International",                     "United Kingdom",          "high",
  "SPC",                                                   "Pacific Community (Secretariat of the Pacific Community)", "International",                  "New Caledonia/Pacific",   "high",
  "CAS",                                                   "California Academy of Sciences",                        "International",                     "United States",           "high",
  "OGL",                                                   "Ocean Genome Legacy Center, Northeastern University",   "International",                     "United States",           "high",
  "UNIVERSITY OF GLASGOW",                                 "University of Glasgow",                                 "International",                     "United Kingdom",          "high",
  "NIWA",                                                  "National Institute of Water and Atmospheric Research",  "International",                     "New Zealand",             "high",
  "NIWA INVERTEBRATE COLLECTION",                          "National Institute of Water and Atmospheric Research",  "International",                     "New Zealand",             "high",
  "SIO",                                                   "Scripps Institution of Oceanography",                   "International",                     "United States",           "high",
  "RBINS",                                                 "Royal Belgian Institute of Natural Sciences",           "International",                     "Belgium",                 "high",
  "RBINS-SCIENTIFIC HERITAGE",                             "Royal Belgian Institute of Natural Sciences",           "International",                     "Belgium",                 "high",
  "LACM",                                                  "Natural History Museum of Los Angeles County",         "International",                     "United States",           "high",
  "SBMNH",                                                 "Santa Barbara Museum of Natural History",               "International",                     "United States",           "high",
  "STRI",                                                  "Smithsonian Tropical Research Institute",               "International",                     "Panama",                  "high",
  "YPM",                                                   "Yale Peabody Museum",                                   "International",                     "United States",           "high",
  "NRM",                                                   "Swedish Museum of Natural History",                     "International",                     "Sweden",                  "high",
  "ROM",                                                   "Royal Ontario Museum",                                  "International",                     "Canada",                  "high",
  "UM-RSMAS",                                              "Univ. of Miami Rosenstiel School of Marine Science",   "International",                     "United States",           "high",
  "BBM",                                                   "Bernice Pauahi Bishop Museum",                          "International",                     "United States (Hawaii)",  "high",
  "DMNS",                                                  "Denver Museum of Nature and Science",                   "International",                     "United States",           "high",
  "TCWC",                                                  "Texas A&M Biodiversity Research and Teaching Collections", "International",                   "United States",           "high",
  "EARTH SCIENCES NEW ZEALAND",                            "Earth Sciences New Zealand",                            "International",                     "New Zealand",             "high",
  "FMNH",                                                  "Field Museum of Natural History",                       "International",                     "United States",           "high",
  "MACN",                                                  "Museo Argentino de Ciencias Naturales",                 "International",                     "Argentina",               "high",
  "UIO",                                                   "University of Oslo",                                    "International",                     "Norway",                  "high",
  "ZMB",                                                   "Museum f\u00fcr Naturkunde Berlin",                     "International",                     "Germany",                 "high",
  "AK",                                                    "Auckland War Memorial Museum, NZ Marine Collection",    "International",                     "New Zealand",             "high",
  "LMU",                                                    "Ludwig-Maximilians-Universit\u00e4t M\u00fcnchen (SNSB - Bavarian State Zoological Collection); voucher specimen linked to a GEOME/DIPnet genomic sample record", "International", "Germany", "high",
  # --- Non-physical-specimen platforms (excluded from location analysis) ---
  "INATURALIST",                                            "iNaturalist (observation platform, no voucher)",        "Non-physical/excluded",              NA,                        "n/a",
  "DIVEBOARD",                                              "Diveboard (observation platform, no voucher)",          "Non-physical/excluded",              NA,                        "n/a",
  "ATLAS OF LIVING AUSTRALIA",                              "ALA (aggregator platform, not a holding institution)",  "Non-physical/excluded",              NA,                        "n/a",
  "OCEAN CONNECT",                                          "Rocky Shore Citizen Science project (observation only, no voucher)", "Non-physical/excluded",   NA,                        "high",
  "S9-5",                                                   "Asia-Pacific literature-compiled occurrence dataset (JAMSTEC/Japan MoE) - not physical specimens", "Non-physical/excluded", NA, "high",
  "ENTANGLED WILDLIFE AUSTRALIA",                           "Marine & Freshwater Fishing Entanglements project (iNaturalist human observation, no voucher)", "Non-physical/excluded", NA, "high"
)

cat("=== Institution classification lookup ===\n\n")
cat(sprintf("Total institution codes classified: %d\n", nrow(institution_lookup)))
cat(sprintf("  High confidence:  %d\n", sum(institution_lookup$confidence == "high")))
cat(sprintf("  Low/unverified:   %d\n", sum(institution_lookup$confidence == "low_unverified")))
cat(sprintf("  Non-physical (n/a): %d\n\n", sum(institution_lookup$confidence == "n/a")))

echino_wide <- echino_wide %>%
  left_join(institution_lookup, by = "institution_std")

n_unmatched <- sum(!is.na(echino_wide$institution_std) &
                     !echino_wide$institution_std %in% institution_lookup$institution_std)
if (n_unmatched > 0) {
  cat(sprintf("\u26A0\uFE0F  %d records have an institution_std value NOT in the lookup table -\n", n_unmatched))
  cat("   these are genuinely new/unclassified codes. Check:\n")
  print(echino_wide %>%
          filter(!is.na(institution_std), !institution_std %in% institution_lookup$institution_std) %>%
          count(institution_std, sort = TRUE))
}

# -----------------------------------------------------------------------------
# RESTRICT TO PHYSICAL-SPECIMEN RECORDS (Flag 1), excluding
# citizen-science/observation records and non-physical platforms
# -----------------------------------------------------------------------------
# Two separate exclusion steps happen here, and were previously being
# collapsed into a single, misleadingly-labeled total (nrow(echino_wide) -
# nrow(provenance_data)), which actually spans BOTH steps below combined -
# not just the "no institution code" step it appeared to describe. Split
# out explicitly so each number means exactly what it says.
# Flag 1 only: preserved, physically held specimens
flag1_data <- echino_wide %>% filter(obs_quality_flag == 1)

provenance_data <- flag1_data %>%
  filter(location_category %in% c("QLD", "Other Australian state/territory", "International"))

n_excluded_step1 <- nrow(echino_wide) - nrow(flag1_data)
n_excluded_step2 <- nrow(flag1_data) - nrow(provenance_data)

cat(sprintf("\n=== Records included in provenance analysis ===\n"))
cat(sprintf("Full dataset: %d records\n", nrow(echino_wide)))
cat(sprintf("Excluded (Flag 2–3: survey/observation records and citizen science): %d\n", n_excluded_step1))
cat(sprintf("Flag 1 (specimen-based) records: %d\n", nrow(flag1_data)))
cat(sprintf("Of those, excluded (no institution code, or non-physical platform code: %d\n", n_excluded_step2))
cat(sprintf("Of those, with a classifiable (high or low-confidence) institution: %d\n",
            nrow(provenance_data)))

# -----------------------------------------------------------------------------
# SUMMARY: where are QLD-collected specimens now held?
# -----------------------------------------------------------------------------
cat("\n\n=== Specimen holding location summary ===\n\n")
location_summary <- provenance_data %>%
  count(location_category, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(location_summary)
write_csv(location_summary, "specimen_provenance_location_summary.csv")

cat("\n── Detail: which specific institutions/states/countries ──\n\n")
detail_summary <- provenance_data %>%
  count(location_category, state_or_country, full_name, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 2))
print(detail_summary, n = 100)
write_csv(detail_summary, "specimen_provenance_detail.csv")

# -----------------------------------------------------------------------------
# COUNTRY-LEVEL ROLLUP: Australia treated as ONE country here (collapsing
# the QLD / other-Australian-state distinction, which is already covered
# by location_summary/detail_summary above), so international countries
# can be ranked directly against Australia as a whole - e.g. "USA holds
# X%, driven primarily by [museum]" rather than a flat 44-row list.
# -----------------------------------------------------------------------------
cat("\n\n=== Country-level summary ===\n\n")

provenance_country <- provenance_data %>%
  mutate(country = if_else(location_category == "International", state_or_country, "Australia"))

country_summary <- provenance_country %>%
  count(country, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1))
cat("── Records by country (Australia = all states/territories combined) ──\n\n")
print(country_summary, n = Inf)
write_csv(country_summary, "specimen_provenance_country_summary.csv")

cat("\n── Within each country, the leading holding institution ──\n\n")
country_top_institution <- provenance_country %>%
  count(country, full_name, sort = TRUE) %>%
  group_by(country) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(country_summary %>% dplyr::select(country, country_total = n, country_pct = pct), by = "country") %>%
  mutate(pct_of_country = round(100 * n / country_total, 1)) %>%
  arrange(desc(country_total)) %>%
  dplyr::select(country, country_total, country_pct, leading_institution = full_name,
                n_at_leading_institution = n, pct_of_country)
print(country_top_institution, n = Inf)
write_csv(country_top_institution, "specimen_provenance_country_leading_institution.csv")

cat("\n── Full institution breakdown within each country (top 3 per country) ──\n\n")
country_institution_detail <- provenance_country %>%
  count(country, full_name, sort = TRUE) %>%
  group_by(country) %>%
  mutate(pct_of_country = round(100 * n / sum(n), 1)) %>%
  slice_max(n, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(n))
print(country_institution_detail, n = 100)
write_csv(country_institution_detail, "specimen_provenance_country_institution_detail.csv")

# --- Figure: country totals, ranked ---
p_country <- country_summary %>%
  mutate(country = fct_reorder(country, n)) %>%
  ggplot(aes(x = country, y = n)) +
  geom_col(fill = "#2a78d6") +
  geom_text(aes(label = if_else(n / sum(n) * 100 < 0.1, "<0.1%", sprintf("%.1f%%", pct))),
            hjust = -0.15, size = 3) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of specimen-based records",
       title = "Specimen holdings by country",
       subtitle = "Australia = all states/territories combined; specimen-based records only (Flag 1)") +
  theme_minimal(base_size = 12)
ggsave("figure_specimen_provenance_by_country.png", p_country, width = 8, height = 7, dpi = 300)
cat("\nFigure saved: figure_specimen_provenance_by_country.png\n")

cat("\n=== Country-level summary complete ===\n")
cat("Outputs: specimen_provenance_country_summary.csv,\n")
cat("specimen_provenance_country_leading_institution.csv,\n")
cat("specimen_provenance_country_institution_detail.csv,\n")
cat("figure_specimen_provenance_by_country.png\n")

# -----------------------------------------------------------------------------
# FIGURE: proportion of specimens held in QLD vs other Australia vs international
# -----------------------------------------------------------------------------
p_provenance <- location_summary %>%
  mutate(location_category = factor(location_category,
                                    levels = c("QLD", "Other Australian state/territory", "International"))) %>%
  ggplot(aes(x = location_category, y = pct)) +
  geom_col(fill = "#2a78d6") +
geom_text(aes(label = if_else(n / sum(n) * 100 < 0.1, "<0.1%", sprintf("%.1f%%", pct))),
          hjust = -0.15, size = 3) +  scale_y_continuous(limits = c(0, max(location_summary$pct) * 1.2),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "% of specimen-based records",
       title = "Where are Queensland-collected echinoderm specimens now held?",
       subtitle = sprintf("n = %d specimen-based records (Flag 1) with a classifiable institution",
                          nrow(provenance_data))) +
  theme_minimal(base_size = 12)

ggsave("figure_specimen_provenance.png", p_provenance, width = 8, height = 6, dpi = 300)
cat("\nFigure saved: figure_specimen_provenance.png\n")

# -----------------------------------------------------------------------------
# International/interstate breakdown figure (which specific institutions
# hold the most QLD material outside Queensland)
# -----------------------------------------------------------------------------
top_external <- detail_summary %>%
  filter(location_category != "QLD") %>%
  arrange(desc(n)) %>%
  head(15)

p_external <- top_external %>%
  mutate(full_name = fct_reorder(full_name, n)) %>%
  ggplot(aes(x = full_name, y = n, fill = location_category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Other Australian state/territory" = "#eda100",
                               "International" = "#e34948")) +
  labs(x = NULL, y = "Number of records", fill = "Location",
       title = "Top institutions holding QLD-collected specimens outside Queensland",
       subtitle = "Specimen-based records only (Flag 1)") +
  theme_minimal(base_size = 12)

ggsave("figure_specimen_provenance_top_institutions.png", p_external, width = 11, height = 7, dpi = 300)
cat("Figure saved: figure_specimen_provenance_top_institutions.png\n")

cat("\n=== Specimen provenance analysis complete ===\n")
cat("Outputs: specimen_provenance_location_summary.csv, specimen_provenance_detail.csv,\n")
cat("figure_specimen_provenance.png, figure_specimen_provenance_top_institutions.png\n")
cat("\nAll institution codes are now classified with 'high' confidence, following\n")
cat("manual verification of BRS, AK, LMU, Ocean Connect, S9-5, and Entangled\n")
cat("Wildlife Australia (see comments in the lookup table above for sources).\n")

# =============================================================================
# TEMPORAL ANALYSIS: has specimen provenance changed over the collection's
# history? Connects directly to Temporal Trends (best_year/decade already
# resolved there) - tests whether local retention has improved or worsened
# across collecting eras, rather than only reporting a static snapshot.
# =============================================================================
cat("\n\n=== TEMPORAL ANALYSIS: specimen provenance by collection decade ===\n\n")

if (!"best_year" %in% names(provenance_data)) {
  stop("best_year not found - run the Section 3 addition in ",
       "ECHINODERM_POST-PROCESSING.R (majority-vote year resolution) before this section.")
}

provenance_decade <- provenance_data %>%
  filter(!is.na(best_year)) %>%
  mutate(decade = floor(best_year / 10) * 10)

n_no_year <- nrow(provenance_data) - nrow(provenance_decade)
cat(sprintf("Provenance records with a resolvable collection year: %d of %d (%.1f%%)\n",
            nrow(provenance_decade), nrow(provenance_data),
            100 * nrow(provenance_decade) / nrow(provenance_data)))
cat(sprintf("Excluded (no resolvable year): %d\n\n", n_no_year))

# --- Descriptive: % held locally (QLD), by decade ---
provenance_by_decade <- provenance_decade %>%
  count(decade, location_category) %>%
  group_by(decade) %>%
  mutate(n_decade = sum(n), pct = round(100 * n / n_decade, 1)) %>%
  ungroup()

cat("── % of specimens by holding location, per collection decade ──\n\n")
provenance_by_decade_wide <- provenance_by_decade %>%
  dplyr::select(decade, location_category, pct) %>%
  pivot_wider(names_from = location_category, values_from = pct, values_fill = 0) %>%
  left_join(provenance_decade %>% count(decade, name = "n_records"), by = "decade") %>%
  arrange(decade)
print(provenance_by_decade_wide, n = Inf)
write_csv(provenance_by_decade_wide, "specimen_provenance_by_decade.csv")

# --- Figure: stacked proportional bar, provenance composition by decade ---
p_provenance_decade <- provenance_by_decade %>%
  mutate(location_category = factor(location_category,
                                    levels = c("QLD", "Other Australian state/territory", "International"))) %>%
  ggplot(aes(x = factor(decade), y = pct, fill = location_category)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("QLD" = "#2a78d6", "Other Australian state/territory" = "#eda100",
                               "International" = "#e34948")) +
  labs(x = "Collection decade", y = "% of specimen-based records", fill = "Holding location",
       title = "Specimen holding location by collection decade",
       subtitle = "Specimen-based records only (Flag 1), with a resolvable collection year") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figure_specimen_provenance_by_decade.png", p_provenance_decade, width = 9, height = 6, dpi = 300)
cat("\nFigure saved: figure_specimen_provenance_by_decade.png\n")

# --- Formal test: chi-square, decade x location_category ---
cat("\n── Chi-square test: holding location vs. collection decade ──\n\n")
cat("H0: the distribution of holding location (QLD / other Australia /\n")
cat("    international) is the same across collection decades\n\n")

decade_loc_table <- table(provenance_decade$decade, provenance_decade$location_category)
print(decade_loc_table)

chisq_prov_standard <- chisq.test(decade_loc_table)
n_low_expected_prov <- sum(chisq_prov_standard$expected < 5)

cat(sprintf("\nMinimum expected cell count: %.2f\n", min(chisq_prov_standard$expected)))
cat(sprintf("Cells below 5: %d of %d\n", n_low_expected_prov, length(chisq_prov_standard$expected)))

if (n_low_expected_prov > 0) {
  low_rows <- rownames(chisq_prov_standard$expected)[
    apply(chisq_prov_standard$expected, 1, function(r) any(r < 5))
  ]
  cat(sprintf("\n\u26A0\uFE0F  Standard chi-square not valid (%d cells <5, driven by decade(s):\n", n_low_expected_prov))
  cat("   ", paste(low_rows, collapse = ", "), ") - using Monte Carlo simulated\n")
  cat("   p-value instead (matching the identical correction already applied\n")
  cat("   to Taxonomic Completeness and Temporal Trends' chi-square tests).\n\n")
  
  set.seed(42)
  chisq_prov_result <- chisq.test(decade_loc_table, simulate.p.value = TRUE, B = 10000)
  print(chisq_prov_result)
  cat("\nUSE THIS p-value in the thesis (cite: Monte Carlo simulation, B=10000).\n")
} else {
  chisq_prov_result <- chisq_prov_standard
  print(chisq_prov_result)
  cat("\nAll expected counts >= 5 - standard chi-square is valid as-is.\n")
}

write_csv(
  tibble(
    statistic = "X-squared",
    value = unname(chisq_prov_result$statistic),
    df_note = if (n_low_expected_prov > 0) "NA (Monte Carlo simulated)" else as.character(unname(chisq_prov_standard$parameter)),
    p_value = chisq_prov_result$p.value,
    method = if (n_low_expected_prov > 0) "Monte Carlo simulation (B=10000)" else "Pearson chi-square",
    min_expected_count = min(chisq_prov_standard$expected),
    cells_below_5 = n_low_expected_prov
  ),
  "specimen_provenance_by_decade_chisq_test.csv"
)
cat("\nSaved: specimen_provenance_by_decade_chisq_test.csv\n")

# --- Simple trend check: % held in QLD, earliest vs most recent well-
# sampled decades (descriptive support for interpreting the chi-square
# result's direction, not a separate formal test) ---
qld_trend <- provenance_by_decade_wide %>% dplyr::select(decade, QLD, n_records)
cat("\n── QLD-retention trend (for interpreting the direction of the chi-square result) ──\n\n")
print(qld_trend, n = Inf)

cat("\n=== Temporal provenance analysis complete ===\n")
cat("Outputs: specimen_provenance_by_decade.csv, figure_specimen_provenance_by_decade.png,\n")
cat("specimen_provenance_by_decade_chisq_test.csv\n")

# =============================================================================
# ANALYSIS 2: DOES PROVENANCE DIFFER BY ECHINODERM CLASS?
# =============================================================================
# Were certain classes disproportionately retained locally vs. shipped
# elsewhere? Mirrors the same by-class chi-square approach already used
# for TRi (Taxonomic Representativeness Index, Section 2.x.3).
# =============================================================================
cat("\n\n=== ANALYSIS 2: holding location by echinoderm class ===\n\n")

provenance_class <- provenance_data %>% filter(!is.na(echino_class))
n_no_class <- nrow(provenance_data) - nrow(provenance_class)
cat(sprintf("Provenance records with a known echino_class: %d of %d (%.1f%%)\n",
            nrow(provenance_class), nrow(provenance_data),
            100 * nrow(provenance_class) / nrow(provenance_data)))
cat(sprintf("Excluded (unknown class): %d\n\n", n_no_class))

class_provenance_table <- table(provenance_class$echino_class, provenance_class$location_category)
cat("── Holding location by echinoderm class (counts) ──\n\n")
print(class_provenance_table)

class_provenance_pct <- provenance_class %>%
  count(echino_class, location_category) %>%
  group_by(echino_class) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()
cat("\n── Holding location by echinoderm class (% within class) ──\n\n")
print(class_provenance_pct %>%
        dplyr::select(echino_class, location_category, pct) %>%
        pivot_wider(names_from = location_category, values_from = pct, values_fill = 0),
      n = Inf)
write_csv(class_provenance_pct, "provenance_by_class.csv")

# --- Chi-square test ---
cat("\n── Chi-square test: holding location vs. echinoderm class ──\n\n")
cat("H0: the distribution of holding location is the same across classes\n\n")

chisq_class_standard <- chisq.test(class_provenance_table)
n_low_class <- sum(chisq_class_standard$expected < 5)
cat(sprintf("Minimum expected cell count: %.2f\n", min(chisq_class_standard$expected)))
cat(sprintf("Cells below 5: %d of %d\n", n_low_class, length(chisq_class_standard$expected)))

if (n_low_class > 0) {
  cat("\n\u26A0\uFE0F  Standard chi-square not valid (sparse cells) - using Monte Carlo\n")
  cat("   simulated p-value instead:\n\n")
  set.seed(42)
  chisq_class_result <- chisq.test(class_provenance_table, simulate.p.value = TRUE, B = 10000)
  print(chisq_class_result)
  cat("\nUSE THIS p-value in the thesis (cite: Monte Carlo simulation, B=10000).\n")
} else {
  chisq_class_result <- chisq_class_standard
  print(chisq_class_result)
  cat("\nAll expected counts >= 5 - standard chi-square is valid as-is.\n")
}

write_csv(
  tibble(
    statistic = "X-squared",
    value = unname(chisq_class_result$statistic),
    df_note = if (n_low_class > 0) "NA (Monte Carlo simulated)" else as.character(unname(chisq_class_standard$parameter)),
    p_value = chisq_class_result$p.value,
    method = if (n_low_class > 0) "Monte Carlo simulation (B=10000)" else "Pearson chi-square",
    min_expected_count = min(chisq_class_standard$expected),
    cells_below_5 = n_low_class
  ),
  "provenance_by_class_chisq_test.csv"
)

# --- Figure: stacked proportional bar, provenance composition by class ---
p_class_provenance <- class_provenance_pct %>%
  mutate(location_category = factor(location_category,
                                    levels = c("QLD", "Other Australian state/territory", "International"))) %>%
  ggplot(aes(x = echino_class, y = pct, fill = location_category)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("QLD" = "#2a78d6", "Other Australian state/territory" = "#eda100",
                               "International" = "#e34948")) +
  labs(x = "Echinoderm class", y = "% of specimen-based records", fill = "Holding location",
       title = "Specimen holding location by echinoderm class",
       subtitle = "Specimen-based records only (Flag 1)") +
  theme_minimal(base_size = 12)
ggsave("figure_provenance_by_class.png", p_class_provenance, width = 8, height = 6, dpi = 300)
cat("\nFigure saved: figure_provenance_by_class.png\n")

cat("\n=== Provenance x class analysis complete ===\n")
cat("Outputs: provenance_by_class.csv, provenance_by_class_chisq_test.csv,\n")
cat("figure_provenance_by_class.png\n")