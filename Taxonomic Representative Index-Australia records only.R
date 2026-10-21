# ==============================================================================
# TAXONOMIC REPRESENTATIVENESS INDEX (TRi) - Australia Records only
# ==============================================================================
# TRi = CSi / CSe x 100
#   CSe = species documented by AFD as occurring in the study region
#         (Queensland/Torres Strait/Gulf of Carpentaria/Coral Sea),
#         WoRMS-resolved to current accepted names
#   CSi = distinct CSe species also present in the occurrence dataset
#
# This is the SINGLE, FINAL version of the whole TRi pipeline, consolidating
# everything built and debugged across multiple sessions:
#   PART A: retrieve species lists per class from ALA, scrape AFD for
#           distribution data, apply known manual corrections
#   PART B: build CSe using reviewed IMCRA regions
#   PART C: WoRMS-resolve every CSe species to its current accepted name
#   PART D: compute CSi and TRi (overall, by class, depth-stratified
#           secondary), plus gap-species outputs in both directions,
#           plus figures for both the by-class and by-depth-zone results
#
# Run top to bottom. Each part is checkpointed to disk, so re-running after
# an interruption skips work already done. To force a completely fresh
# pull, delete: afd_echinoderm_species_list.csv, afd_echinoderm_master.csv,
# cse_afd_checklist.csv, cse_worms_resolved.csv
# ==============================================================================
# =============================================================================
# SETUP
# =============================================================================
setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")
# =============================================================================
# Library
# =============================================================================
library(httr)
library(jsonlite)
library(rvest)
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
CLASSES            <- c("Asteroidea", "Ophiuroidea", "Echinoidea", "Holothuroidea", "Crinoidea")
SPECIES_LIST_CSV   <- "afd_echinoderm_species_list.csv"
MASTER_CSV         <- "afd_echinoderm_master.csv"
CSE_AFD_CSV        <- "cse_afd_checklist.csv"
CSE_RESOLVED_CSV   <- "cse_worms_resolved.csv"
REQUEST_DELAY      <- 1.0
# ==============================================================================
# PART A: RETRIEVE SPECIES LISTS + SCRAPE AFD DISTRIBUTION DATA
# ==============================================================================
resolve_class_guid <- function(class_name) {
  url <- paste0("https://bie-ws.ala.org.au/ws/search.json?q=", class_name, "&fq=rank:class")
  resp <- try(GET(url, timeout(15)), silent = TRUE)
  if (inherits(resp, "try-error") || status_code(resp) != 200) {
    warning(paste("GUID lookup failed for", class_name)); return(NA_character_)
  }
  parsed <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
  results <- parsed$searchResults$results
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("No GUID found for", class_name)); return(NA_character_)
  }
  hit <- results %>% filter(str_to_lower(name) == str_to_lower(class_name))
  if (nrow(hit) == 0) hit <- results[1, ]
  hit$guid[1]
}
EMPTY_SPECIES_TBL <- tibble(class = character(), species = character())
fetch_class_species_list <- function(class_name, guid, n_attempts = 5) {
  if (is.na(guid)) return(EMPTY_SPECIES_TBL)
  
  url <- paste0("https://bie-ws.ala.org.au/ws/download?q=",
                URLencode(paste0('rkid_class:"', guid, '"')),
                "&fq=rank:species")
  
  all_results <- vector("list", n_attempts)
  for (i in seq_len(n_attempts)) {
    resp <- try(GET(url, timeout(60)), silent = TRUE)
    if (inherits(resp, "try-error") || status_code(resp) != 200) next
    
    raw_text <- content(resp, as = "text", encoding = "UTF-8")
    df <- try(read_csv(raw_text, show_col_types = FALSE), silent = TRUE)
    if (inherits(df, "try-error") || !"scientificName" %in% names(df)) next
    
    all_results[[i]] <- df %>%
      filter(taxonRank == "species") %>%
      transmute(class = class_name, species = scientificName)
    
    if (i < n_attempts) Sys.sleep(1.0)
  }
  
  combined <- bind_rows(all_results) %>% distinct(class, species)
  message(sprintf("    %s: union of %d attempts = %d distinct species",
                  class_name, n_attempts, nrow(combined)))
  if (nrow(combined) == 0) warning(paste("All attempts failed for", class_name))
  combined
}
message("=== PART A: species lists + AFD distribution scrape ===\n")
if (file.exists(SPECIES_LIST_CSV)) {
  message("Found existing species list, loading from disk")
  species_master <- read_csv(SPECIES_LIST_CSV, show_col_types = FALSE)
} else {
  species_master <- map_dfr(CLASSES, function(cls) {
    guid <- resolve_class_guid(cls)
    message(sprintf("  %-15s guid: %s", cls, ifelse(is.na(guid), "NOT FOUND", guid)))
    Sys.sleep(REQUEST_DELAY)
    fetch_class_species_list(cls, guid)
  })
  
  if (nrow(species_master) == 0 || !"species" %in% names(species_master)) {
    stop("All class downloads failed - check ALA API status before re-running.")
  }
  
  species_master <- species_master %>%
    filter(
      !is.na(species),
      str_detect(species, "^[A-Z][a-z]+"),
      !str_detect(species, regex("\\bsp\\.?\\b|\\bspp\\.?\\b|\\bcf\\.|\\baff\\.|\\bnr\\.|\\d|\u00d7", ignore_case = TRUE))
    ) %>%
    distinct(class, species)
  
  write_csv(species_master, SPECIES_LIST_CSV)
}
message(sprintf("\n-> %d species across %d classes", nrow(species_master), n_distinct(species_master$class)))
print(count(species_master, class))
# --- Scraping helpers ---
extract_afd_field <- function(txt, label) {
  m <- str_match(txt, paste0(label, "\\n+([^\\n]+)"))
  m[1, 2]
}
parse_depth_range <- function(depth_str) {
  if (is.na(depth_str)) return(c(min = NA_real_, max = NA_real_))
  rng <- str_match(depth_str, "([\\d.]+)\\s*-\\s*([\\d.]+)\\s*m")
  if (!is.na(rng[1, 1])) return(c(min = as.numeric(rng[1, 2]), max = as.numeric(rng[1, 3])))
  single <- str_match(depth_str, "([\\d.]+)\\s*m")
  if (!is.na(single[1, 1])) return(c(min = NA_real_, max = as.numeric(single[1, 2])))
  c(min = NA_real_, max = NA_real_)
}
scrape_afd_species <- function(binomial) {
  slug <- str_replace_all(binomial, " ", "_")
  url  <- paste0("https://biodiversity.org.au/afd/taxa/", slug)
  
  resp <- try(GET(url, timeout(15), user_agent("JCU thesis research")), silent = TRUE)
  if (inherits(resp, "try-error") || status_code(resp) != 200) {
    return(tibble(species = binomial, afd_url = url, scrape_status = "fetch_failed"))
  }
  txt <- try(read_html(resp) %>% html_text2(), silent = TRUE)
  if (inherits(txt, "try-error") || is.na(txt)) {
    return(tibble(species = binomial, afd_url = url, scrape_status = "parse_failed"))
  }
  
  states_block <- extract_afd_field(txt, "States")
  extra_block  <- extract_afd_field(txt, "Extra Distribution Information")
  imcra_block  <- extract_afd_field(txt, "IMCRA")
  eco_block    <- extract_afd_field(txt, "Ecological Descriptors")
  depth        <- str_extract(extra_block %||% "", "Depth range[^.]*m")
  depth_range  <- parse_depth_range(depth)
  
  tibble(
    species        = binomial,
    afd_url        = url,
    states         = states_block,
    qld_flag       = str_detect(states_block %||% "", regex("Queensland", ignore_case = TRUE)),
    imcra_regions  = imcra_block,
    depth_text     = depth,
    depth_min_m    = depth_range["min"],
    depth_max_m    = depth_range["max"],
    ecological_descriptors = eco_block,
    scrape_status  = "ok"
  )
}
message("\nScraping AFD species pages...")
already_done <- if (file.exists(MASTER_CSV)) read_csv(MASTER_CSV, show_col_types = FALSE) else tibble(species = character())
to_scrape <- species_master %>% anti_join(already_done, by = "species") %>% pull(species)
message(sprintf("%d already scraped, %d remaining", nrow(already_done), length(to_scrape)))
if (length(to_scrape) > 0) {
  results_new <- vector("list", length(to_scrape))
  for (i in seq_along(to_scrape)) {
    results_new[[i]] <- scrape_afd_species(to_scrape[i])
    Sys.sleep(REQUEST_DELAY)
    if (i %% 50 == 0 || i == length(to_scrape)) {
      message(sprintf("  ...%d / %d scraped", i, length(to_scrape)))
      write_csv(bind_rows(already_done, bind_rows(results_new)), MASTER_CSV)
    }
  }
  species_detail <- bind_rows(already_done, bind_rows(results_new))
} else {
  species_detail <- already_done
}
afd_master <- species_master %>%
  dplyr::select(class, species) %>%
  left_join(species_detail %>% dplyr::select(-any_of("class")), by = "species") %>%
  dplyr::select(class, species, afd_url, states, qld_flag, imcra_regions,
                depth_text, depth_min_m, depth_max_m, ecological_descriptors, scrape_status)
known_fixes <- tribble(
  ~species,                                              ~correct_url,
  "Asterodiscides soelae",     "https://biodiversity.org.au/afd/taxa/Asterodiscides_soleae",
  "Stephanometra tenuispina",  "https://biodiversity.org.au/afd/taxa/Stephanometra_tenuipinna",
  "Cheiraster teres",          "https://biodiversity.org.au/afd/taxa/Cheiraster_(Luidiaster)_teres",
  "Antedon loveni",            "https://biodiversity.org.au/afd/taxa/Colobometra_perspinosa",
  "Astropecten zebra",         "https://biodiversity.org.au/afd/taxa/Astropecten_vappa",
  "Macrophiothrix microplax",  "https://biodiversity.org.au/afd/taxa/Macrophiothrix_longipeda",
  "Deima validum",             "https://biodiversity.org.au/afd/taxa/Deima_validum_validum",
  "Oneirophanta mutabilis",    "https://biodiversity.org.au/afd/taxa/Oneirophanta_mutabilis_mutabilis",
  "Actinopyga mauritiana",     "https://biodiversity.org.au/afd/taxa/Actinopyga_mauritiana;Actinopyga",
  "Araeosoma coriaceum",       "https://biodiversity.org.au/afd/taxa/Araeosoma_coriaceum;Araeosoma",
  "Nanometra johnstoni",       "https://biodiversity.org.au/afd/taxa/Nanometra_johnstoni;Nanometra",
  "Ophiocreas sibogae",        "https://biodiversity.org.au/afd/taxa/Ophiocreas_sibogae;Ophiocreas",
  "Ophiomastix elegans",       "https://biodiversity.org.au/afd/taxa/Ophiomastix_elegans;Ophiomastix"
)
excluded_from_afd <- tribble(
  ~species,                   ~reason,
  "Pseudoceramaster glasbyi", "Valid WoRMS name, no matching record in AFD - likely not yet reflected in AFD's database"
)
to_fix <- afd_master %>% filter(species %in% known_fixes$species, scrape_status != "ok")
if (nrow(to_fix) > 0) {
  message(sprintf("\nApplying %d known manual correction(s)...", nrow(to_fix)))
  fix_results <- vector("list", nrow(to_fix))
  for (i in seq_len(nrow(to_fix))) {
    sp  <- to_fix$species[i]
    url <- known_fixes$correct_url[known_fixes$species == sp]
    resp <- try(GET(url, timeout(20), user_agent("JCU thesis research")), silent = TRUE)
    if (inherits(resp, "try-error") || status_code(resp) != 200) {
      fix_results[[i]] <- tibble(species = sp, afd_url = url, scrape_status = "fetch_failed")
      next
    }
    txt <- read_html(resp) %>% html_text2()
    states_block <- extract_afd_field(txt, "States")
    extra_block  <- extract_afd_field(txt, "Extra Distribution Information")
    imcra_block  <- extract_afd_field(txt, "IMCRA")
    eco_block    <- extract_afd_field(txt, "Ecological Descriptors")
    depth        <- str_extract(extra_block %||% "", "Depth range[^.]*m")
    depth_range  <- parse_depth_range(depth)
    fix_results[[i]] <- tibble(
      species = sp, afd_url = url, states = states_block,
      qld_flag = str_detect(states_block %||% "", regex("Queensland", ignore_case = TRUE)),
      imcra_regions = imcra_block, depth_text = depth,
      depth_min_m = depth_range["min"], depth_max_m = depth_range["max"],
      ecological_descriptors = eco_block, scrape_status = "ok"
    )
    message(sprintf("  [%d/%d] %s -> %s", i, nrow(to_fix), sp, fix_results[[i]]$scrape_status))
    Sys.sleep(1.5)
  }
  fix_results <- bind_rows(fix_results)
  afd_master <- afd_master %>%
    filter(!(species %in% fix_results$species)) %>%
    bind_rows(fix_results %>% dplyr::select(any_of(names(afd_master)))) %>%
    arrange(class, species)
}
write_csv(afd_master, MASTER_CSV)
message("\n=== Scrape summary ===")
message(sprintf("Total species: %d | ok: %d | still failing: %d",
                nrow(afd_master), sum(afd_master$scrape_status == "ok"),
                sum(afd_master$scrape_status != "ok")))
still_failing <- afd_master %>% filter(scrape_status != "ok") %>% pull(species)
if (length(still_failing) > 0) {
  message("Still failing (expected: only genuinely AFD-absent species):")
  print(still_failing)
}
# ==============================================================================
# PART B: BUILD CSe USING REVIEWED IMCRA REGIONS
# ==============================================================================
message("\n\n=== PART B: building CSe ===\n")
QLD_RELEVANT_IMCRA <- c(
  "Northeast Shelf Province", "Northeast Shelf Transition",
  "Northeast Province", "Northeast Transition",
  "Central Eastern Shelf Province", "Central Eastern Shelf Transition",
  "Central Eastern Province", "Central Eastern Transition",
  "Northern Shelf Province",
  "Kenn Province", "Kenn Transition"
)
imcra_inventory <- afd_master %>%
  filter(!is.na(imcra_regions)) %>%
  separate_rows(imcra_regions, sep = ",\\s*") %>%
  mutate(imcra_label = str_trim(str_remove(imcra_regions, "\\s*\\(\\d+\\)$"))) %>%
  count(imcra_label, sort = TRUE)
message("IMCRA labels currently in use (top 10, for reference):")
print(head(imcra_inventory, 10))
species_imcra_flag <- afd_master %>%
  filter(!is.na(imcra_regions)) %>%
  separate_rows(imcra_regions, sep = ",\\s*") %>%
  mutate(imcra_label = str_trim(str_remove(imcra_regions, "\\s*\\(\\d+\\)$"))) %>%
  group_by(species) %>%
  summarise(imcra_qld_match = any(imcra_label %in% QLD_RELEVANT_IMCRA), .groups = "drop")
species_summary <- afd_master %>%
  group_by(species) %>%
  summarise(
    class         = first(na.omit(class)),
    states        = first(na.omit(states)),
    qld_flag      = any(qld_flag, na.rm = TRUE),
    depth_min_m   = suppressWarnings(min(depth_min_m, na.rm = TRUE)),
    depth_max_m   = suppressWarnings(max(depth_max_m, na.rm = TRUE)),
    scrape_status = first(na.omit(scrape_status)),
    .groups = "drop"
  ) %>%
  mutate(
    depth_min_m = if_else(is.infinite(depth_min_m), NA_real_, depth_min_m),
    depth_max_m = if_else(is.infinite(depth_max_m), NA_real_, depth_max_m)
  )
cse_afd <- species_summary %>%
  left_join(species_imcra_flag, by = "species") %>%
  mutate(
    imcra_qld_match = coalesce(imcra_qld_match, FALSE),
    in_cse = qld_flag | imcra_qld_match
  ) %>%
  dplyr::select(class, species, states, qld_flag, imcra_qld_match, in_cse,
                depth_min_m, depth_max_m, scrape_status)
write_csv(cse_afd, CSE_AFD_CSV)
message(sprintf("\nCSe (AFD, pre-WoRMS-resolution): %d species", sum(cse_afd$in_cse)))
print(cse_afd %>% filter(in_cse) %>% count(class, name = "n_species"))
# ==============================================================================
# PART C: WoRMS-RESOLVE EVERY CSe SPECIES TO ITS CURRENT ACCEPTED NAME
# ==============================================================================
message("\n\n=== PART C: WoRMS resolution ===\n")
lookup_worms_accepted_name <- function(binomial) {
  url <- paste0("https://www.marinespecies.org/rest/AphiaRecordsByName/",
                URLencode(binomial), "?marine_only=false")
  resp <- try(GET(url, timeout(15)), silent = TRUE)
  if (inherits(resp, "try-error") || status_code(resp) != 200) {
    return(list(status = "worms_fetch_failed", accepted_name = NA_character_))
  }
  raw_text <- content(resp, as = "text", encoding = "UTF-8")
  parsed <- try(fromJSON(raw_text, flatten = TRUE), silent = TRUE)
  if (inherits(parsed, "try-error") || length(parsed) == 0 ||
      (is.data.frame(parsed) && nrow(parsed) == 0)) {
    return(list(status = "worms_no_match", accepted_name = NA_character_))
  }
  df <- if (is.data.frame(parsed)) parsed else bind_rows(parsed)
  if ("status" %in% names(df)) {
    accepted_row <- df %>% filter(status == "accepted")
    if (nrow(accepted_row) >= 1) df <- accepted_row
  }
  rec <- df[1, ]
  list(status = "ok", accepted_name = rec$valid_name %||% rec$scientificname %||% NA_character_)
}
cse <- read_csv(CSE_AFD_CSV, show_col_types = FALSE) %>% filter(in_cse)
already_resolved <- if (file.exists(CSE_RESOLVED_CSV)) {
  read_csv(CSE_RESOLVED_CSV, show_col_types = FALSE) %>%
    dplyr::select(species, worms_accepted_name, worms_match_status)
} else {
  tibble(species = character())
}
to_resolve <- cse %>% anti_join(already_resolved, by = "species") %>% pull(species)
message(sprintf("%d already resolved, %d remaining", nrow(already_resolved), length(to_resolve)))
if (length(to_resolve) == 0 && nrow(already_resolved) != nrow(cse)) {
  warning("Cached WoRMS resolution size doesn't match current CSe size but 0 species ",
          "need resolving - the cache is likely stale. Delete ", CSE_RESOLVED_CSV, " and re-run.")
}
if (length(to_resolve) > 0) {
  results <- vector("list", length(to_resolve))
  for (i in seq_along(to_resolve)) {
    lookup <- lookup_worms_accepted_name(to_resolve[i])
    results[[i]] <- tibble(
      species = to_resolve[i],
      worms_accepted_name = if (lookup$status == "ok") lookup$accepted_name else NA_character_,
      worms_match_status = lookup$status
    )
    Sys.sleep(REQUEST_DELAY)
    if (i %% 50 == 0 || i == length(to_resolve)) {
      message(sprintf("  ...%d / %d resolved", i, length(to_resolve)))
      write_csv(bind_rows(already_resolved, bind_rows(results)), CSE_RESOLVED_CSV)
    }
  }
  resolved_lookup <- bind_rows(already_resolved, bind_rows(results))
} else {
  resolved_lookup <- already_resolved
}
cse_resolved <- cse %>%
  left_join(resolved_lookup, by = "species") %>%
  mutate(
    match_name = coalesce(worms_accepted_name, species),
    name_was_updated = !is.na(worms_accepted_name) & worms_accepted_name != species
  )
write_csv(cse_resolved, CSE_RESOLVED_CSV)
message(sprintf("\nWoRMS resolution: %d succeeded, %d name(s) updated, %d fell back to AFD name",
                sum(cse_resolved$worms_match_status == "ok", na.rm = TRUE),
                sum(cse_resolved$name_was_updated, na.rm = TRUE),
                sum(cse_resolved$worms_match_status != "ok" | is.na(cse_resolved$worms_match_status))))
dup_check <- cse_resolved %>% count(match_name) %>% filter(n > 1)
if (nrow(dup_check) > 0) {
  message(sprintf(
    "Note: %d match_name(s) appear twice post-resolution (two AFD entries resolving to the same current species) - CSe uses DISTINCT match_name.",
    nrow(dup_check)
  ))
}
cse_final_n <- n_distinct(cse_resolved$match_name)
message(sprintf("\nFINAL CSe (distinct species, WoRMS-resolved): %d", cse_final_n))
# ==============================================================================
# PART D: COMPUTE CSi AND TRi
# ==============================================================================
message("\n\n=== PART D: computing CSi and TRi ===\n")
normalize_name <- function(x) {
  x %>% str_replace_all("\\([A-Za-z]+\\)", " ") %>% str_squish() %>% str_to_lower()
}
normalize_name_binomial_only <- function(x) {
  x %>% str_replace_all("\\([A-Za-z]+\\)", " ") %>% str_squish() %>% str_to_lower() %>% word(1, 2)
}
known_spelling_aliases <- tribble(
  ~variant_a,                  ~variant_b,
  "stichopus hermanni",        "stichopus herrmanni",
  "stephanometra tenuispina",  "stephanometra tenuipinna",
  "leiaster leachii",          "leiaster leachi",
  "ophiocentrus asperus",      "ophiocentrus aspera"
)
resolve_known_aliases <- function(key) {
  match_idx <- match(key, known_spelling_aliases$variant_b)
  if_else(!is.na(match_idx), known_spelling_aliases$variant_a[match_idx], key)
}
echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)
if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c("Species"="Species","Subspecies"="Species","Genus"="Genus","Subgenus"="Genus",
                "Family"="Family","Order"="Order","Class"="Class","Phylum"="Phylum")
  echino_wide <- echino_wide %>% mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}
if (!"crosses_zone_boundary" %in% names(echino_wide)) {
  echino_wide <- echino_wide %>% mutate(crosses_zone_boundary = str_starts(coalesce(depth_zone, ""), "Spans"))
}

# ------------------------------------------------------------------------------
# GEOGRAPHIC SCOPE SENSITIVITY (PNG + New Caledonia)
# ------------------------------------------------------------------------------
TRI_SCOPE <- "australia_only"            # <- run once with "all", then once with "australia_only"

echino_wide <- echino_wide %>%
  mutate(
    loc_text   = str_to_lower(paste(verbatimLocality, `Locality Name`, `Locality Remarks`)),
    outside_au = coalesce(country != "Australia", FALSE) |
      str_detect(loc_text, "new caledonia|nouvelle.cal|noum[e\u00e9]a")
  )

cat(sprintf("Records flagged outside Australia: %d (PNG %d, other %d)\n",
            sum(echino_wide$outside_au),
            sum(echino_wide$country == "Papua New Guinea", na.rm = TRUE),
            sum(echino_wide$outside_au) - sum(echino_wide$country == "Papua New Guinea", na.rm = TRUE)))
# expected: 1023 flagged (PNG 1017, other 6)

echino_tri <- if (TRI_SCOPE == "australia_only") filter(echino_wide, !outside_au) else echino_wide
cat(sprintf("TRi scope: %s | %d records used\n", TRI_SCOPE, nrow(echino_tri)))

# every output from here on goes into its own folder (reads above are already done)
base_dir <- getwd()
out_dir  <- file.path(base_dir, paste0("TRi_", TRI_SCOPE))
dir.create(out_dir, showWarnings = FALSE)
setwd(out_dir)
occurrence_species <- echino_tri %>%
  filter(rank6 == "Species", !is.na(accepted_name_final)) %>%
  distinct(accepted_name_final) %>%
  mutate(match_key = resolve_known_aliases(normalize_name(accepted_name_final)),
         match_key_binomial = resolve_known_aliases(normalize_name_binomial_only(accepted_name_final)))
cse_normalized <- cse_resolved %>%
  distinct(match_name, class) %>%
  mutate(match_key = resolve_known_aliases(normalize_name(match_name)),
         match_key_binomial = resolve_known_aliases(normalize_name_binomial_only(match_name)))
cse_total <- n_distinct(cse_normalized$match_key)
matched_full <- cse_normalized %>% inner_join(occurrence_species, by = "match_key")
cse_unmatched_full <- cse_normalized %>% anti_join(occurrence_species, by = "match_key")
matched_binomial_fallback <- cse_unmatched_full %>%
  dplyr::select(-match_key) %>%
  inner_join(occurrence_species %>% dplyr::select(-match_key), by = "match_key_binomial") %>%
  distinct(match_key_binomial, class, .keep_all = TRUE)
matched_species <- bind_rows(
  matched_full %>% mutate(matched_via = "full_name"),
  matched_binomial_fallback %>% mutate(matched_via = "binomial_fallback")
) %>% mutate(csi_id = coalesce(match_key, match_key_binomial))
csi_total <- n_distinct(matched_species$csi_id)
tri_overall <- round(csi_total / cse_total * 100, 1)
message(sprintf("CSe = %d | CSi = %d | TRi = %.1f%%", cse_total, csi_total, tri_overall))
write_csv(tibble(metric = c("CSe", "CSi", "TRi (%)"),
                 value = c(cse_total, csi_total, tri_overall)),
          "TRi_overall_summary.csv")
# By class
tri_by_class <- cse_normalized %>%
  group_by(class) %>% summarise(CSe = n_distinct(match_key), .groups = "drop") %>%
  left_join(matched_species %>% group_by(class) %>% summarise(CSi = n_distinct(csi_id), .groups = "drop"),
            by = "class") %>%
  mutate(CSi = coalesce(CSi, 0L), TRi = round(CSi / CSe * 100, 1)) %>%
  arrange(desc(TRi))
print(tri_by_class, n = Inf)
write_csv(tri_by_class, "TRi_by_class.csv")
# -----------------------------------------------------------------------------
# CHI-SQUARE TEST: does the proportion of CSe species recovered differ
# significantly among echinoderm classes?
# -----------------------------------------------------------------------------
message("\n── Chi-square test: TRi across echinoderm classes ──\n")
message("H0: the proportion of CSe species recovered does not differ\n")
message("    significantly among classes\n")
class_matrix <- tri_by_class %>%
  transmute(collected = CSi, not_collected = CSe - CSi) %>%
  as.matrix()
rownames(class_matrix) <- tri_by_class$class
chisq_by_class <- chisq.test(class_matrix)
print(chisq_by_class)
message("\n── Expected cell count check ──\n")
expected_table <- as_tibble(round(chisq_by_class$expected, 2)) %>%
  mutate(class = tri_by_class$class, .before = 1)
print(expected_table)
min_expected <- min(chisq_by_class$expected)
n_below_5 <- sum(chisq_by_class$expected < 5)
message(sprintf("\nMinimum expected cell count: %.2f", min_expected))
message(sprintf("Cells below 5: %d of %d", n_below_5, length(chisq_by_class$expected)))
if (n_below_5 == 0) {
  message("\nVERIFIED: all expected cell counts exceed 5 - standard Pearson")
  message("chi-square is valid as run, no correction needed. This differs from")
  message("the Taxonomic Completeness chi-square (8/24 cells <5, Monte Carlo")
  message("simulation used) and TIa's pairwise tests (2/6 comparisons <5,")
  message("corrected to Fisher's exact test) - this by-class comparison has")
  message("larger per-cell sample sizes and does not share their violation.")
  chi_pvalue_for_reporting <- chisq_by_class$p.value
} else {
  message("\nWARNING: chi-square approximation may not be reliable here -")
  message("consider Monte Carlo simulation or Fisher's exact test instead.")
  chi_pvalue_for_reporting <- NA_real_
}
write_csv(expected_table, "TRi_by_class_expected_cell_counts.csv")
write_csv(
  tibble(statistic = "X-squared", value = unname(chisq_by_class$statistic),
         df = unname(chisq_by_class$parameter), p_value = chisq_by_class$p.value,
         min_expected_count = min_expected, cells_below_5 = n_below_5),
  "TRi_by_class_chisq_test.csv"
)
# Gap species (both directions)
matched_cse_keys <- unique(c(matched_full$match_key, matched_binomial_fallback$match_key_binomial))
gap_species <- cse_normalized %>%
  filter(!(match_key %in% matched_cse_keys), !(match_key_binomial %in% matched_cse_keys)) %>%
  dplyr::select(class, species = match_name) %>% arrange(class, species)
write_csv(gap_species, "csi_gap_species.csv")
occurrence_record_counts <- echino_tri %>%
  filter(rank6 == "Species", !is.na(accepted_name_final)) %>%
  count(accepted_name_final, name = "n_records")
occurrence_only_species <- occurrence_species %>%
  filter(!(match_key %in% matched_cse_keys), !(match_key_binomial %in% matched_cse_keys)) %>%
  left_join(occurrence_record_counts, by = "accepted_name_final") %>%
  dplyr::select(species = accepted_name_final, n_records) %>% arrange(desc(n_records))
write_csv(occurrence_only_species, "occurrence_only_species.csv")
message(sprintf("\nGap species: %d in CSe not collected | %d collected not in CSe",
                nrow(gap_species), nrow(occurrence_only_species)))
# Secondary: depth-stratified TRi
zone_bounds <- tribble(
  ~depth_bin,           ~zone_min, ~zone_max,
  "Continental Shelf",  0,         200,
  "Upper Slope",        200,       1000,
  "Lower Slope",        1000,      4000,
  "Abyssal",    4000,      Inf
)
cse_with_depth <- cse_normalized %>%
  left_join(cse_resolved %>% distinct(match_name, depth_min_m, depth_max_m) %>%
              mutate(match_key = resolve_known_aliases(normalize_name(match_name))), by = "match_key") %>%
  filter(!is.na(depth_min_m) | !is.na(depth_max_m))
cse_zone_expected <- cse_with_depth %>%
  mutate(d_min = coalesce(depth_min_m, 0), d_max = coalesce(depth_max_m, 6000)) %>%
  dplyr::cross_join(zone_bounds) %>%
  filter(d_min <= zone_max, d_max >= zone_min) %>%
  distinct(match_key, class, depth_bin)
occurrence_species_by_zone <- echino_tri %>%
  filter(rank6 == "Species", !is.na(accepted_name_final), !crosses_zone_boundary,
         depth_zone %in% zone_bounds$depth_bin) %>%
  distinct(depth_zone, accepted_name_final) %>%
  rename(depth_bin = depth_zone) %>%
  mutate(match_key = resolve_known_aliases(normalize_name(accepted_name_final)))
tri_by_zone <- cse_zone_expected %>%
  group_by(depth_bin) %>% summarise(CSe_zone = n_distinct(match_key), .groups = "drop") %>%
  left_join(
    cse_zone_expected %>% inner_join(occurrence_species_by_zone, by = c("depth_bin", "match_key")) %>%
      group_by(depth_bin) %>% summarise(CSi_zone = n_distinct(match_key), .groups = "drop"),
    by = "depth_bin"
  ) %>%
  mutate(CSi_zone = coalesce(CSi_zone, 0L), TRi_zone = round(CSi_zone / CSe_zone * 100, 1),
         depth_bin = factor(depth_bin, levels = zone_bounds$depth_bin)) %>%
  arrange(depth_bin)
depth_coverage_pct <- round(100 * n_distinct(cse_with_depth$match_key) / cse_total, 1)
message(sprintf("\nDepth-stratified TRi (SECONDARY - %d/%d CSe species with depth data, %s%%):",
                n_distinct(cse_with_depth$match_key), cse_total, depth_coverage_pct))
print(tri_by_zone, n = Inf)
write_csv(tri_by_zone, "TRi_depth_stratified_secondary.csv")

# ==============================================================================
# TRi BY CLASS x DEPTH ZONE (combined matrix) - EXPLORATORY/SECONDARY
# ==============================================================================
# Reframes the separate TRi-by-class and TRi-by-depth-zone results above into
# a single combined cross-tabulation: for each echinoderm class, does
# regional representativeness collapse specifically at certain depths,
# rather than uniformly across its whole depth range? This is a genuinely
# new cut of already-computed data - cse_zone_expected already carries
# class alongside depth_bin, so this simply groups by both simultaneously
# rather than each independently.
#
# SAME "exploratory/secondary" caveat as the depth-stratified TRi above
# applies here (restricted to the subset of CSe species with usable AFD
# depth data) - AND is compounded further here, since splitting by class
# on top of depth produces some genuinely small CSe cells. Cells with
# CSe < 5 are flagged explicitly (matching the expected-cell-count
# convention used for chi-square validity throughout this thesis) rather
# than presented as a precise, trustworthy percentage - e.g. TRi=100% from
# a CSe of 1 species is not a meaningful representativeness estimate.
# ==============================================================================
cat("\n\n=== TRi BY CLASS x DEPTH ZONE (combined matrix, exploratory/secondary) ===\n\n")

# Depth-range legend labels (name on top, range below), matching the
# convention used across every other depth-zone figure in this thesis -
# defined here (rather than at its original location further below) so
# both this new heatmap and the existing by-zone bar figure can use it
zone_labels_with_range <- c(
  "Continental Shelf" = "Continental Shelf\n(0 - 200 m)",
  "Upper Slope"        = "Upper Slope\n(200 - 1000 m)",
  "Lower Slope"         = "Lower Slope\n(1000 - 4000 m)",
  "Abyssal"     = "Abyssal\n(>4000 m)"
)

tri_by_class_zone <- cse_zone_expected %>%
  group_by(class, depth_bin) %>%
  summarise(CSe_cz = n_distinct(match_key), .groups = "drop") %>%
  left_join(
    cse_zone_expected %>%
      inner_join(occurrence_species_by_zone, by = c("depth_bin", "match_key")) %>%
      group_by(class, depth_bin) %>%
      summarise(CSi_cz = n_distinct(match_key), .groups = "drop"),
    by = c("class", "depth_bin")
  ) %>%
  mutate(
    CSi_cz = coalesce(CSi_cz, 0L),
    TRi_cz = round(CSi_cz / CSe_cz * 100, 1),
    reliable = CSe_cz >= 5,
    depth_bin = factor(depth_bin, levels = zone_bounds$depth_bin)
  ) %>%
  arrange(class, depth_bin)

cat("── Full class x depth-zone TRi matrix (CSi/CSe, TRi%) ──\n\n")
print(tri_by_class_zone, n = Inf)
n_unreliable <- sum(!tri_by_class_zone$reliable)
cat(sprintf("\n%d of %d class x depth-zone cells have CSe < 5 - flagged as unreliable\n",
            n_unreliable, nrow(tri_by_class_zone)))
cat("(not excluded from the table, but should not be quoted as a precise\n")
cat("representativeness estimate; reported alongside its raw CSi/CSe count)\n\n")
write_csv(tri_by_class_zone, "TRi_by_class_and_depth_zone.csv")

cat("── Wide format: TRi% by class (rows) x depth zone (columns) ──\n\n")
tri_class_zone_wide <- tri_by_class_zone %>%
  dplyr::select(class, depth_bin, TRi_cz) %>%
  pivot_wider(names_from = depth_bin, values_from = TRi_cz)
print(tri_class_zone_wide, n = Inf)
write_csv(tri_class_zone_wide, "TRi_by_class_and_depth_zone_wide.csv")

# --- Figure: heatmap, class x depth zone, TRi% as fill, unreliable cells
# marked distinctly rather than shown as a solid, trustworthy colour ---

p_tri_heatmap <- tri_by_class_zone %>%
  mutate(
    label_text = ifelse(reliable,
                        sprintf("%.0f%%\n(%d/%d)", TRi_cz, CSi_cz, CSe_cz),
                        sprintf("(%d/%d)\nsmall n", CSi_cz, CSe_cz)),
    # Text colour: white on dark/low-TRi cells (viridis maps low values to
    # dark purple), dark grey on light/high-TRi cells. Faded (unreliable)
    # cells always get dark text, since the alpha blending with the white
    # background already lightens them regardless of their underlying
    # TRi value - a low-TRi unreliable cell is NOT rendered as dark purple.
    text_colour = ifelse(!reliable | TRi_cz >= 35, "grey15", "white")
  ) %>%
  ggplot(aes(x = depth_bin, y = class, fill = TRi_cz)) +
  geom_tile(aes(alpha = reliable), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label_text, colour = text_colour),
            size = 3, lineheight = 0.9, show.legend = FALSE) +
  scale_fill_viridis_c(option = "viridis", name = "TRi (%)", na.value = "grey90") +
  scale_colour_identity() +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
  scale_x_discrete(labels = zone_labels_with_range) +
  labs(x = "Depth zone", y = "Echinoderm class",
       title = "Taxonomic Representativeness Index by class and depth zone",
       subtitle = sprintf(
         "Exploratory/secondary: CSe restricted to %d of %d species (%.1f%%) with usable AFD depth data.\nFaded cells (CSe<5) are unreliable and shown for completeness only, not as precise estimates.",
         n_distinct(cse_with_depth$match_key), cse_total, depth_coverage_pct
       )) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())
ggsave("figure_TRi_by_class_and_depth_zone.png", p_tri_heatmap, width = 10, height = 6.5, dpi = 300)
cat("Figure saved: figure_TRi_by_class_and_depth_zone.png\n")
# Vector version for poster
ggsave("figure_TRi_by_class_and_depth_zone.pdf",  p_tri_heatmap, width = 10, height = 6.5,
       device = cairo_pdf)
cat("\nFigure saved: figure_TRi_by_class_and_depth_zone.pdf (vector)\n")


cat("Outputs: TRi_by_class_and_depth_zone.csv, TRi_by_class_and_depth_zone_wide.csv,\n")
cat("figure_TRi_by_class_and_depth_zone.png\n")

cat("\n=== TRi by class x depth zone complete ===\n")


# ==============================================================================
# TRi BY DEPTH ZONE - FIGURE (exploratory/secondary analysis)
# ==============================================================================
# Consumes tri_by_zone directly (not a hardcoded duplicate) - the figure
# always matches the live computation above.
#
# CAVEAT carried into the figure's own subtitle: this is the EXPLORATORY/
# SECONDARY depth-stratified TRi, restricted to the subset of CSe species
# with usable AFD-documented depth ranges - not the primary TRi result.
#
# Two-line zone labels (name on top, depth range below) match the
# convention used across the other depth-zone figures in this thesis
# (Taxonomic Completeness Figure 1, TIa Figures 4-5). A small point marker
# is added at each bar's actual value, since Abyssal's TRi is 0% here -
# a zero-height bar would otherwise be invisible with nothing marking it.
# ==============================================================================

tri_by_zone_plot <- tri_by_zone %>%
  mutate(depth_bin_plot = factor(depth_bin, levels = rev(zone_bounds$depth_bin)))

axis_max <- max(tri_by_zone_plot$TRi_zone) * 1.35  # just enough headroom for the label

p_tri_zone <- ggplot(tri_by_zone_plot, aes(x = depth_bin_plot, y = TRi_zone)) +
  geom_col(fill = "#2a78d6", width = 0.65) +
  geom_point(size = 1.8, colour = "#0c447c") +
  geom_text(aes(label = sprintf("%.1f%%  (%d/%d)", TRi_zone, CSi_zone, CSe_zone)),
            hjust = -0.08, size = 3.4, colour = "grey20") +
  coord_flip(ylim = c(0, axis_max), clip = "off") +
  scale_x_discrete(labels = zone_labels_with_range) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Taxonomic Representativeness Index (TRi, %)",
       title = "Taxonomic Representativeness Index by depth zone",
       subtitle = "Exploratory: CSe restricted to species with usable AFD depth data (56.9%)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

ggsave("figure_TRi_by_depth_zone.png", p_tri_zone, width = 8, height = 4, dpi = 300)
cat("Figure saved: figure_TRi_by_depth_zone.png\n")

message("\n=== TRi pipeline complete ===")
message("Outputs: TRi_overall_summary.csv, TRi_by_class.csv, csi_gap_species.csv,")
message("occurrence_only_species.csv, TRi_depth_stratified_secondary.csv,")
message("TRi_by_class_chisq_test.csv, TRi_by_class_expected_cell_counts.csv,")
message("figure_TRi_by_depth_zone.png")
# ==============================================================================
# TRi BY ECHINODERM CLASS - FIGURE
# ==============================================================================
tri_by_class_plot <- tri_by_class %>%
  mutate(
    class = factor(class, levels = class[order(-TRi)]),
    label_y = TRi + 6
  )
overall_TRi <- round(100 * sum(tri_by_class_plot$CSi) / sum(tri_by_class_plot$CSe), 1)
cat(sprintf("Overall TRi (reference line): %.1f%% (CSi=%d, CSe=%d)\n",
            overall_TRi, sum(tri_by_class_plot$CSi), sum(tri_by_class_plot$CSe)))
chisq_subtitle_text <- if (!is.na(chi_pvalue_for_reporting)) {
  sprintf("chi-square test: %s (p = %.3f)",
          if (chi_pvalue_for_reporting < 0.05) "significant" else "not significant",
          chi_pvalue_for_reporting)
} else {
  "chi-square test: expected cell count assumption violated - see console output"
}
p_tri_class <- ggplot(tri_by_class_plot, aes(x = class, y = TRi)) +
  geom_col(fill = "#2a78d6") +
  geom_hline(yintercept = overall_TRi, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_text(aes(y = label_y, label = sprintf("%.1f%%\n(%d/%d)", TRi, CSi, CSe)),
            size = 3, colour = "grey20", lineheight = 0.9) +
  scale_y_continuous(limits = c(0, max(tri_by_class_plot$label_y) * 1.15),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Echinoderm class", y = "Taxonomic Representativeness Index (TRi, %)",
       title = "Taxonomic Representativeness Index by echinoderm class",
       subtitle = sprintf("Dashed line = overall TRi (%.1f%%, CSi=%d of CSe=%d); %s",
                          overall_TRi, sum(tri_by_class_plot$CSi), sum(tri_by_class_plot$CSe),
                          chisq_subtitle_text)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave("figure_TRi_by_class.png", p_tri_class, width = 9, height = 6, dpi = 300)
cat("Figure saved: figure_TRi_by_class.png\n")
print(tri_by_class_plot %>% dplyr::select(class, CSe, CSi, TRi) %>% arrange(desc(TRi)))

setwd(base_dir)
