# ==============================================================================
# DISTRIBUTION MAPS
# ==============================================================================
# Map 1: all occurrence records, coloured by depth zone
# Map 2: species richness AND sampling density per grid cell
#
# Input:  echino_wide.csv
# Output: map1_records_by_depth_zone.png
#         map2_species_richness_grid.png
#         map2_sampling_density_grid.png
#         grid_richness_density_summary.csv
# ==============================================================================
# =============================================================================
# SETUP
# =============================================================================
setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")
# =============================================================================
# Library
# =============================================================================
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggspatial)
echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)
if (!"rank6" %in% names(echino_wide)) {
  rank_map <- c("Species"="Species","Subspecies"="Species","Genus"="Genus","Subgenus"="Genus",
                "Family"="Family","Order"="Order","Class"="Class","Phylum"="Phylum")
  echino_wide <- echino_wide %>% mutate(rank6 = recode(taxonomic_resolution_level, !!!rank_map))
}
zone_order <- c("Continental Shelf", "Upper Slope", "Lower Slope", "Abyssal")
zone_labels_with_range <- c(
  "Continental Shelf" = "Continental Shelf\n(0 - 200 m)",
  "Upper Slope"        = "Upper Slope\n(200 - 1000 m)",
  "Lower Slope"         = "Lower Slope\n(1000 - 4000 m)",
  "Abyssal"     = "Abyssal\n(>4000 m)"
)
lat_candidates <- c("best_latitude", "decimalLatitude", "latitude", "lat")
lon_candidates <- c("best_longitude", "decimalLongitude", "longitude", "lon", "lng")
lat_col <- lat_candidates[lat_candidates %in% names(echino_wide)][1]
lon_col <- lon_candidates[lon_candidates %in% names(echino_wide)][1]
if (is.na(lat_col) || is.na(lon_col)) {
  stop("Could not find latitude/longitude columns automatically.")
}
message(sprintf("Using coordinate columns: %s (latitude), %s (longitude)", lat_col, lon_col))
echino_geo <- echino_wide %>%
  rename(lat = all_of(lat_col), lon = all_of(lon_col)) %>%
  filter(!is.na(lat), !is.na(lon))
message(sprintf("%d of %d records have usable coordinates (%.1f%%)",
                nrow(echino_geo), nrow(echino_wide), 100 * nrow(echino_geo) / nrow(echino_wide)))
australia <- ne_countries(scale = "large", country = "Australia", returnclass = "sf")
bbox <- echino_geo %>%
  summarise(xmin = min(lon, na.rm = TRUE) - 1, xmax = max(lon, na.rm = TRUE) + 1,
            ymin = min(lat, na.rm = TRUE) - 1, ymax = max(lat, na.rm = TRUE) + 1)
message(sprintf("\nMap extent: lon %.1f to %.1f, lat %.1f to %.1f",
                bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax))
map_aspect <- (bbox$xmax - bbox$xmin) * cos(mean(c(bbox$ymin, bbox$ymax)) * pi / 180) /
  (bbox$ymax - bbox$ymin)
fig_height <- 9.5
fig_width  <- fig_height * map_aspect
message(sprintf("Figure size: %.1f x %.1f in (aspect %.2f)", fig_width, fig_height, map_aspect))
world_context <- ne_countries(scale = "large", continent = "Oceania", returnclass = "sf") %>%
  filter(name != "Australia")
reference_cities <- tribble(
  ~city,             ~lon,   ~lat,    ~label_lon, ~label_lat, ~hjust, ~vjust,
  "Cairns",           145.45, -16.95, 145.30,     -16.95,     1,      0.5,
  "Townsville",       146.45, -19.35, 146.30,     -19.35,     1,      0.5,
  "Mackay",           148.85, -21.20, 148.70,     -21.20,     1,      0.5,
  "Rockhampton",      150.05, -23.40, 149.90,     -23.40,     1,      0.5,
  "Brisbane",         152.55, -27.50, 152.40,     -27.50,     1,      0.5
)
if (!"coord_qc_exclude_recommended" %in% names(echino_geo)) {
  stop("coord_qc_exclude_recommended not found - run Section 1d of ECHINODERM_POST-PROCESSING.R first.")
}
n_flagged <- sum(echino_geo$coord_qc_exclude_recommended, na.rm = TRUE)
cat(sprintf("\n%d of %d georeferenced records flagged coord_qc_exclude_recommended = TRUE (%.2f%%)\n",
            n_flagged, nrow(echino_geo), 100 * n_flagged / nrow(echino_geo)))

echino_geo_qc <- echino_geo %>% filter(!coord_qc_exclude_recommended)

echino_geo_zoned <- echino_geo_qc %>%
  filter(depth_zone %in% zone_order) %>%
  mutate(depth_zone = factor(depth_zone, levels = zone_order))

n_no_zone <- nrow(echino_geo_qc) - nrow(echino_geo_qc %>% filter(depth_zone %in% zone_order))
if (n_no_zone > 0) {
  message(sprintf("\nNote: %d georeferenced records have no depth zone assigned and are excluded from Map 1.", n_no_zone))
}

map1 <- ggplot() +
  geom_sf(data = world_context, fill = "grey93", colour = "grey70", linewidth = 0.25) +
  geom_sf(data = australia, fill = "grey90", colour = "grey60", linewidth = 0.3) +
  geom_point(data = echino_geo_zoned, aes(x = lon, y = lat, colour = depth_zone),
             size = 0.6, alpha = 0.5) +
  geom_point(data = reference_cities, aes(x = lon, y = lat),
             shape = 18, size = 2.2, colour = "grey20") +
  geom_text(data = reference_cities,
            aes(x = label_lon, y = label_lat, label = city, hjust = hjust, vjust = vjust),
            size = 3, colour = "grey20", fontface = "italic") +
  annotate("text", x = 144.5, y = -23, label = "AUSTRALIA",
           size = 5, colour = "grey45", fontface = "bold") +
  annotate("text", x = 157, y = -19.5, label = "CORAL SEA",
           size = 4, colour = "grey40", fontface = "italic") +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  scale_colour_manual(
    values = c("Continental Shelf" = "#2a78d6", "Upper Slope" = "#eda100",
               "Lower Slope" = "#e34948", "Abyssal" = "#4a3aa7"),
    limits = zone_order,
    labels = zone_labels_with_range[zone_order]
  ) +
  labs(x = NULL, y = NULL, colour = "Depth zone",
       title = "Echinoderm occurrence records across the study region",
       subtitle = sprintf("n = %d georeferenced, depth-zoned records", nrow(echino_geo_zoned)),
       caption = "Map shown in geographic coordinates (WGS84); scale bar is approximate\nand most accurate at the map's central latitude.") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1.5, "cm"), width = unit(1.5, "cm"),
                         pad_x = unit(0.3, "in"), pad_y = unit(0.85, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", panel.grid = element_line(colour = "grey92"))
ggsave("map1_records_by_depth_zone.png", map1, width = fig_width, height = fig_height, dpi = 300)
cat("Figure saved: map1_records_by_depth_zone.png\n")

# ==============================================================================
# MAP 1b: ALL GEOREFERENCED RECORDS, INCLUDING THOSE WITHOUT DEPTH
# ==============================================================================
# Map 1 above only plots records with a resolved depth zone (depth_zone %in%
# zone_order), silently excluding the ~43% of records with depth_zone =
# "No depth data". This matters: checking the Coral Sea cluster specifically
# found over half of it (4,446 of ~8,002 candidate records) has no depth
# data and would be entirely invisible on Map 1. This companion map plots
# every georeferenced record regardless of depth status, with "No depth
# data" as its own explicit category rather than an invisible gap.
# ==============================================================================
cat("\n=== MAP 1b: all georeferenced records, including no-depth-data ===\n\n")

zone_order_with_nodata <- c(zone_order, "No depth data")
zone_labels_with_nodata <- c(zone_labels_with_range, "No depth data" = "No depth\ndata recorded")

echino_geo_all <- echino_geo %>%
  filter(!coord_qc_exclude_recommended) %>%
  mutate(depth_zone_plot = if_else(depth_zone %in% zone_order, depth_zone, "No depth data"),
         depth_zone_plot = factor(depth_zone_plot, levels = zone_order_with_nodata))

n_nodata <- sum(echino_geo_all$depth_zone_plot == "No depth data")
cat(sprintf("Records plotted: %d (%d depth-zoned + %d with no depth data)\n",
            nrow(echino_geo_all), nrow(echino_geo_all) - n_nodata, n_nodata))

map_scale <- 2  

map1b <- ggplot() +
  geom_sf(data = world_context, fill = "grey93", colour = "grey70", linewidth = 0.25) +
  geom_sf(data = australia, fill = "grey90", colour = "grey60", linewidth = 0.3) +
  geom_point(data = echino_geo_all %>% filter(depth_zone_plot != "No depth data"),
             aes(x = lon, y = lat, colour = depth_zone_plot,
                 shape = depth_zone_plot, size = depth_zone_plot),
             alpha = 0.5) +
  geom_point(data = echino_geo_all %>% filter(depth_zone_plot != "No depth data"),
             aes(x = lon, y = lat, colour = depth_zone_plot), size = 1.1, alpha = 0.5) +
  geom_point(data = reference_cities, aes(x = lon, y = lat),
             shape = 18, size = 4, colour = "grey20") +
  geom_text(data = reference_cities,
            aes(x = label_lon, y = label_lat, label = city, hjust = hjust, vjust = vjust),
            size = 6, colour = "grey20", fontface = "italic") +
  annotate("text", x = 144.5, y = -23, label = "AUSTRALIA",
           size = 9, colour = "grey45", fontface = "bold") +
  annotate("text", x = 157, y = -19.5, label = "CORAL SEA",
           size = 7, colour = "grey40", fontface = "italic") +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  scale_colour_manual(
    values = c("Continental Shelf" = "#2a78d6", "Upper Slope" = "#eda100",
               "Lower Slope" = "#e34948", "Abyssal" = "#4a3aa7", "No depth data" = "grey60"),
    limits = zone_order_with_nodata,
    labels = zone_labels_with_nodata[zone_order_with_nodata]
  ) +
  scale_shape_manual(
    values = c("Continental Shelf" = 16, "Upper Slope" = 16, "Lower Slope" = 16,
                                             "Abyssal" = 17, "No depth data" = 16),
    limits = zone_order_with_nodata, guide = "none"
  ) +
  scale_size_manual(
    values = c("Continental Shelf" = 1.1, "Upper Slope" = 1.1, "Lower Slope" = 1.1,
                                             "Abyssal" = 2, "No depth data" = 1.1),
    limits = zone_order_with_nodata, guide = "none"
  ) +
  labs(x = NULL, y = NULL, colour = "Depth zone",
       title = "Echinoderm occurrence records across the study region",
       subtitle = sprintf("n = %d georeferenced records, including %d (%.1f%%) with no resolved depth",
                          nrow(echino_geo_all), n_nodata, 100 * n_nodata / nrow(echino_geo_all)),
       caption = "Map shown in geographic coordinates (WGS84); scale bar is approximate\nand most accurate at the map's central latitude.") +
  guides(colour = guide_legend(override.aes = list(size = c(5, 5, 5, 6.5, 5),
                                                   shape = c(16, 16, 16, 17, 16),
                                                   alpha = 1))) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(3, "cm"), width = unit(3, "cm"),
                         pad_x = unit(0.3, "in"), pad_y = unit(0.85, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   text_cex = 1.8,
                   pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")) +
  theme_minimal(base_size = 24) +   # was 12 -> bump this
  theme(
    legend.position = "right",
    panel.grid = element_line(colour = "grey92"),
    plot.title = element_text(size = rel(1.25), face = "bold"),
    plot.subtitle = element_text(size = rel(0.75), colour = "grey30"),
    plot.caption = element_text(size = rel(0.55), colour = "grey40"),
    legend.title = element_text(size = rel(1.0))
  )

ggsave("map1b_all_records_including_no_depth.png", map1b,
       width = fig_width * map_scale, height = fig_height * map_scale, dpi = 300)
cat("Figure saved: map1b_all_records_including_no_depth.png\n")

# Sampling density map

GRID_SIZE <- 0.5
echino_geo_for_grid <- echino_geo %>% filter(!coord_qc_exclude_recommended)
cat(sprintf("\nMap 2 grid built from %d records (%d excluded)\n",
            nrow(echino_geo_for_grid), nrow(echino_geo) - nrow(echino_geo_for_grid)))
echino_sf <- st_as_sf(echino_geo_for_grid, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
grid <- st_make_grid(
  st_as_sfc(st_bbox(c(xmin = bbox$xmin, xmax = bbox$xmax, ymin = bbox$ymin, ymax = bbox$ymax), crs = 4326)),
  cellsize = GRID_SIZE, square = TRUE
) %>%
  st_sf(grid_id = seq_along(.), geometry = .)
echino_gridded <- st_join(echino_sf, grid, join = st_within)
density_by_cell <- echino_gridded %>%
  st_drop_geometry() %>%
  filter(!is.na(grid_id)) %>%
  count(grid_id, name = "n_records")
richness_by_cell <- echino_gridded %>%
  st_drop_geometry() %>%
  filter(!is.na(grid_id), rank6 == "Species", !is.na(accepted_name_final)) %>%
  group_by(grid_id) %>%
  summarise(n_species = n_distinct(accepted_name_final), .groups = "drop")
grid_summary <- grid %>%
  left_join(density_by_cell, by = "grid_id") %>%
  left_join(richness_by_cell, by = "grid_id") %>%
  mutate(n_records = coalesce(n_records, 0L), n_species = coalesce(n_species, 0L))
write_csv(st_drop_geometry(grid_summary), "grid_richness_density_summary.csv")
message(sprintf("\nGrid: %d cells total, %d with at least 1 record, %d with at least 1 species-level record",
                nrow(grid_summary), sum(grid_summary$n_records > 0), sum(grid_summary$n_species > 0)))
map2_density <- ggplot() +
  geom_sf(data = grid_summary %>% filter(n_records > 0),
          aes(fill = n_records), colour = NA) +
  geom_sf(data = world_context, fill = NA, colour = "grey60", linewidth = 0.25) +
  geom_sf(data = australia, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_point(data = reference_cities, aes(x = lon, y = lat),
             shape = 18, size = 2.2, colour = "grey20") +
  geom_text(data = reference_cities,
            aes(x = label_lon, y = label_lat, label = city, hjust = hjust, vjust = vjust),
            size = 3, colour = "grey20", fontface = "italic") +
  annotate("text", x = 144.5, y = -23, label = "AUSTRALIA",
           size = 5, colour = "grey45", fontface = "bold") +
  annotate("text", x = 157, y = -19.5, label = "CORAL SEA",
           size = 4, colour = "grey40", fontface = "italic") +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  scale_fill_viridis_c(option = "plasma", trans = "log10",
                       name = "Records\n(log scale)", na.value = "grey95") +
  labs(x = NULL, y = NULL,
       title = "Sampling density",
       subtitle = sprintf("%.1f\u00b0 x %.1f\u00b0 grid cells", GRID_SIZE, GRID_SIZE),
       caption = "Map shown in geographic coordinates (WGS84).") +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1.4, "cm"), width = unit(1.4, "cm"),
                         pad_x = unit(0.3, "in"), pad_y = unit(0.8, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_line(colour = "grey92"))

ggsave("map2_sampling_density_grid.png", map2_density, width = fig_width, height = fig_height, dpi = 300)
cat("Figure saved: map2_sampling_density_grid.png\n")

# Richness map

map2_richness <- ggplot() +
  geom_sf(data = grid_summary %>% filter(n_species > 0),
          aes(fill = n_species), colour = NA) +
  geom_sf(data = world_context, fill = NA, colour = "grey60", linewidth = 0.25) +
  geom_sf(data = australia, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_point(data = reference_cities, aes(x = lon, y = lat),
             shape = 18, size = 2.2, colour = "grey20") +
  geom_text(data = reference_cities,
            aes(x = label_lon, y = label_lat, label = city, hjust = hjust, vjust = vjust),
            size = 3, colour = "grey20", fontface = "italic") +
  annotate("text", x = 144.5, y = -23, label = "AUSTRALIA",
           size = 5, colour = "grey45", fontface = "bold") +
  annotate("text", x = 157, y = -19.5, label = "CORAL SEA",
           size = 4, colour = "grey40", fontface = "italic") +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  scale_fill_viridis_c(option = "viridis", name = "Species\nrichness", na.value = "grey95") +
  labs(x = NULL, y = NULL,
       title = "Species richness",
       subtitle = sprintf("%.1f\u00b0 x %.1f\u00b0 grid cells, species-level records only", GRID_SIZE, GRID_SIZE),
       caption = "Map shown in geographic coordinates (WGS84).") +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1.4, "cm"), width = unit(1.4, "cm"),
                         pad_x = unit(0.3, "in"), pad_y = unit(0.8, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.25,
                   pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_line(colour = "grey92"))
ggsave("map2_species_richness_grid.png", map2_richness, width = fig_width, height = fig_height, dpi = 300)
cat("Figure saved: map2_species_richness_grid.png\n")
check_data <- grid_summary %>% st_drop_geometry() %>% filter(n_records > 0)
cor_density_richness <- cor(check_data$n_records, check_data$n_species, method = "spearman")
cat(sprintf("\nSpearman correlation, sampling density vs species richness per cell: %.2f\n", cor_density_richness))
cat("\n=== Distribution maps complete ===\n")
cat("Outputs: map1_records_by_depth_zone.png, map2_species_richness_grid.png,\n")
cat("map2_sampling_density_grid.png, grid_richness_density_summary.csv\n")