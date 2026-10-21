# =============================================================================
# TAXONOMIC COMPLETENESS - full analysis
# Camila Velez | JCU Marine Biology MSc | Thesis Chapter — Analysis
#
# Scope: the SIX-LEVEL RESOLUTION HIERARCHY (Species/Genus/Family/Order/
# Class/Phylum).
#
# Order (matches the Methods subsection structure):
#   1. Descriptive: rank composition by depth zone (table + figure)
#   2. Chi-square test on the full rank x depth-zone table (global test,
#      no assumptions about depth's functional form)
#   3. GAM check of the depth relationship's shape - justifies how depth
#      enters the regression models below - + breakpoint estimation if
#      non-linearity is found (pooled across all classes)
#   3b. Class-specific version of the same check: does the curve SHAPE differ by
#     echinoderm class? k = 4 chosen a priori to limit overfitting in sparse classes.
#     Sensitivity checks (k = 5, 8) improved AIC but destabilised sparse classes
#     (Crinoidea class coefficient +0.19 -> -0.80 -> -2.04; oscillating curves),
#     so k = 4 is retained. Per-class breakpoints use single-breakpoint segmented models.
#   4. Ordinal logistic regression (proportional odds model), using the
#      depth specification justified in step 3
#   4b. Full pairwise class comparisons (extends 4's vs-Ophiuroidea-only
#      comparisons to all C(5,2) pairs)
#   5. Proportional odds assumption check (Brant test); given a violation,
#      cutpoint-specific binary models reported as the corrected result,
#      with separation detection + Firth's correction applied automatically
#   5b. Pairwise class comparisons, at each cutpoint separately
#   5c. Depth x class interaction, at each cutpoint (LRT) - reported at the
#      omnibus level only; see the note at the end of Section 5c for why a
#      depth-stratified pairwise breakdown was not pursued further.
#  6: Deep-only GAM Analysis 
#
# NOTE: the GAM/breakpoint result from step 3 also justifies the log-depth
# specification used in the TIa logistic regression (section 2.x.2) - it
# is deliberately only run ONCE, here, and referenced (not repeated) there.
#
# Packages needed beyond tidyverse: mgcv, segmented, brant, logistf, emmeans
# (all optional - script degrades gracefully with clear messages if any
# are missing).
# Install if needed:
# install.packages(c("mgcv","segmented","brant","logistf","emmeans"))
#
# CHANGE LOG (this version):
#   - Section 3b: class_breakpoints for the class-specific figure is now
#     built FROM the per-class breakpoint results computed in this same
#     run, instead of being hardcoded. The per-class loop was moved to run
#     BEFORE the figure is built (previously it ran after, leaving the
#     figure permanently out of sync with whatever the loop found).
#   - New shared helper fit_glm_flag_separation(): centralises the
#     separation-detection logic (R's own "fitted probabilities numerically
#     0 or 1" warning + SE > 10x estimate check) so Section 5 and 5c apply
#     IDENTICAL criteria instead of drifting copies.
#   - Section 5c: the interaction model itself is now checked for
#     separation (sparse depth x class cells, e.g. Crinoidea/Echinoidea at
#     Abyssal, are exactly the setup that produces it) and flagged in the
#     output rather than silently trusting the LRT p-value.
#   - segmented()'s breakpoint search is now seeded (BREAKPOINT_SEED) at
#     every call site, since its internal bootstrap-restart is otherwise
#     unseeded and can converge to different local optima - including
#     different accept/reject verdicts on the 3-breakpoint degeneracy
#     check - on IDENTICAL data across reruns.
#   - Section 5d (depth-stratified pairwise class comparisons at
#     representative depths) has been REMOVED. It required splitting the
#     data into depth x class x cutpoint cells thin enough (e.g. ~49
#     Crinoidea records near 1000m) that roughly half its estimates needed
#     a separation correction to be trustworthy, and even then rested on
#     very little data per cell. Section 5c's omnibus LRT result (depth's
#     effect differs by class, at every cutpoint) is kept as the
#     appropriate level of claim for this dataset - see the note at the
#     end of Section 5c.
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================
setwd("C:\\Users\\Camilita\\Desktop\\JCU\\Thesis\\Phase 1")

# =============================================================================
# Library
# =============================================================================
library(tidyverse)

# -----------------------------------------------------------------------------
# SETUP - always loads fresh from disk (see earlier scripts for why this
# matters: reusing an in-memory echino_wide can silently be a filtered copy
# left behind by a previous script in the same session)
# -----------------------------------------------------------------------------
echino_wide <- read_csv("echino_wide.csv", show_col_types = FALSE)

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

if (!"crosses_zone_boundary" %in% names(echino_wide)) {
  echino_wide <- echino_wide %>%
    mutate(crosses_zone_boundary = str_starts(coalesce(depth_zone, ""), "Spans"))
}

# Fixed seed for segmented()'s internal bootstrap-restart algorithm. Without
# this, segmented() can converge to different local optima on IDENTICAL data
# across reruns - not just slightly different breakpoint values, but
# different structural conclusions (e.g. whether the 3-breakpoint degeneracy
# check passes or fails at all). Report this seed value in Methods so the
# breakpoint results are exactly reproducible.
BREAKPOINT_SEED <- 42

zone_order_shallow_first <- c("Continental Shelf", "Upper Slope", "Lower Slope", "Abyssal")

# Depth-range labels for figure axes/legends
zone_labels_with_range <- c(
  "Continental Shelf" = "Continental Shelf\n(0 - 200 m)",
  "Upper Slope"        = "Upper Slope\n(200 - 1000 m)",
  "Lower Slope"         = "Lower Slope\n(1000 - 4000 m)",
  "Abyssal"     = "Abyssal\n(>4000 m)"
)

rank_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")

echino_wide_depth <- echino_wide %>%
  filter(!crosses_zone_boundary, depth_zone %in% zone_order_shallow_first) %>%
  mutate(depth_bin_shallow_first = factor(depth_zone, levels = zone_order_shallow_first))


# =============================================================================
# SECTION 1: DESCRIPTIVE - rank composition by depth zone (table + figure)
# =============================================================================
cat("=== SECTION 1: Rank composition by depth zone (descriptive) ===\n\n")

table1_rank_by_depthbin <- echino_wide_depth %>%
  count(depth_bin_shallow_first, rank6) %>%
  group_by(depth_bin_shallow_first) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  dplyr::select(depth_bin_shallow_first, rank6, pct) %>%
  pivot_wider(names_from = rank6, values_from = pct, values_fill = 0) %>%
  dplyr::select(depth_bin_shallow_first, any_of(rank_order)) %>%
  left_join(echino_wide_depth %>% count(depth_bin_shallow_first, name = "n_records"),
            by = "depth_bin_shallow_first") %>%
  rename(depth_bin = depth_bin_shallow_first)

cat("── Table 1: % records by taxonomic rank, per depth zone ──\n")
print(table1_rank_by_depthbin, n = Inf)
write_csv(table1_rank_by_depthbin, "table1_rank_by_depthbin.csv")

# Figure: 100% stacked horizontal bar, depth-ordered shelf-at-top
rank_by_zone <- echino_wide_depth %>%
  count(depth_bin_shallow_first, rank6) %>%
  group_by(depth_bin_shallow_first) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    rank6 = factor(rank6, levels = rev(rank_order)),
    depth_bin_plot = factor(depth_bin_shallow_first, levels = rev(zone_order_shallow_first))
  )

p_rank <- rank_by_zone %>%
  ggplot(aes(x = depth_bin_plot, y = pct, fill = rank6)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_x_discrete(labels = zone_labels_with_range) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "Phylum"  = "#999999",
      "Class"   = "#E69F00",
      "Order"   = "#CC79A7",
      "Family"  = "#F0E442",
      "Genus"   = "#009E73",
      "Species" = "#0072B2"
    ), breaks = rank_order) +
  labs(x = NULL, y = "% of records", fill = "Rank",
       title = "Taxonomic rank composition across the depth gradient") +
  theme_minimal(base_size = 24) +   # was 12 -> bump this
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = rel(1.3), face = "bold"),
    legend.title = element_text(size = rel(1.1))
  )

# Render at ~2x the original physical size, same aspect ratio (2:1)
ggsave("figure_rank_composition_by_depth.png", p_rank, width = 16, height = 8, dpi = 300)
cat("\nFigure saved: figure_rank_composition_by_depth.png\n")

# Vector version for poster (preferred — scales losslessly when placed)
ggsave("figure_rank_composition_by_depth.pdf", p_rank, width = 16, height = 8,
       device = cairo_pdf)
cat("\nFigure saved: figure_rank_composition_by_depth.pdf (vector)\n")

# =============================================================================
# SECTION 2: CHI-SQUARE TEST on the full rank x depth-zone table
# =============================================================================
cat("\n\n=== SECTION 2: Chi-square test (full rank x depth-zone table) ===\n")
cat("H0: the distribution across all 6 ranks is the same across depth zones\n\n")

zone_data <- echino_wide %>%
  filter(has_depth, !crosses_zone_boundary, depth_zone %in% zone_order_shallow_first) %>%
  mutate(depth_bin = factor(depth_zone, levels = zone_order_shallow_first))

full_table <- table(zone_data$depth_bin, zone_data$rank6)
cat("── Full contingency table ──\n")
print(full_table)

chisq_standard <- chisq.test(full_table)
n_low_expected <- sum(chisq_standard$expected < 5)

cat("\nCells with expected count < 5:", n_low_expected, "of", length(chisq_standard$expected), "\n")

if (n_low_expected > 0) {
  cat("\n\u26A0\uFE0F  Standard chi-square not valid (sparse cells) - using Monte Carlo\n")
  cat("   simulated p-value instead:\n\n")
  set.seed(42)
  chisq_result <- chisq.test(full_table, simulate.p.value = TRUE, B = 10000)
  print(chisq_result)
  cat("\nUSE THIS p-value in the thesis (cite: Monte Carlo simulation, B=10000).\n")
} else {
  chisq_result <- chisq_standard
  print(chisq_result)
  cat("\nAll expected counts >= 5 - standard chi-square is valid as-is.\n")
}


# =============================================================================
# SECTION 3: GAM CHECK of the depth relationship's shape + breakpoint
# estimation. This justifies the log-depth specification used in BOTH the
# ordinal model below AND the TIa logistic regression in section 2.x.2 -
# run once here, referenced (not repeated) there.
# =============================================================================
cat("\n\n=== SECTION 3: GAM check of depth relationship shape ===\n\n")

model_data <- echino_wide %>%
  filter(has_depth, !crosses_zone_boundary, !is.na(echino_class)) %>%
  mutate(
    species_level = as.integer(rank6 == "Species"),
    log_depth = log10(depth_median + 1),
    echino_class = factor(echino_class),
    echino_class = relevel(echino_class, ref = "Ophiuroidea"),
    obs_quality_flag = factor(obs_quality_flag)
  )

final_breakpoints_m <- NULL  # populated once the 2-vs-3 decision is made below

if (requireNamespace("mgcv", quietly = TRUE)) {
  library(mgcv)
  
  model_linear <- glm(
    species_level ~ log_depth + echino_class + obs_quality_flag,
    data = model_data, family = binomial
  )
  model_gam <- gam(
    species_level ~ s(log_depth) + echino_class + obs_quality_flag,
    data = model_data, family = binomial, method = "REML"
  )
  
  cat("── Smooth term summary (edf close to 1 = linear; higher = curved) ──\n\n")
  print(summary(model_gam))
  
  cat("\n── Model comparison ──\n")
  cat("Linear model AIC:", round(AIC(model_linear), 1), "\n")
  cat("GAM model AIC:   ", round(AIC(model_gam), 1), "\n")
  
  edf <- summary(model_gam)$s.table[1, "edf"]
  delta_aic <- AIC(model_linear) - AIC(model_gam)
  
  if (edf > 2 && delta_aic > 2) {
    cat("\nRESULT: relationship is NOT well-approximated by a straight line\n")
    cat("(edf =", round(edf, 1), ", \u0394AIC =", round(delta_aic, 1), "in favour of the GAM).\n")
    cat("Proceeding to estimate the breakpoint location.\n")
  } else {
    cat("\nRESULT: linear log-depth specification appears reasonable\n")
    cat("(edf close to 1 and/or AIC difference small) - no breakpoint needed.\n")
  }
  
  # Save the fitted curve for the figure (figure itself is built further
  # below, AFTER the breakpoint decision - see note there)
  depth_range <- seq(min(model_data$log_depth), max(model_data$log_depth), length.out = 200)
  pred_data <- tibble(
    log_depth = depth_range,
    echino_class = factor("Ophiuroidea", levels = levels(model_data$echino_class)),
    obs_quality_flag = factor("1", levels = levels(model_data$obs_quality_flag))
  )
  pred <- predict(model_gam, newdata = pred_data, type = "response", se.fit = TRUE)
  pred_data <- pred_data %>%
    filter(is.finite(log_depth)) %>%
    mutate(
      depth_m = 10^log_depth - 1,
      fit = pred$fit[is.finite(log_depth)],
      lower = pmax(0, pred$fit[is.finite(log_depth)] - 1.96 * pred$se.fit[is.finite(log_depth)]),
      upper = pmin(1, pred$fit[is.finite(log_depth)] + 1.96 * pred$se.fit[is.finite(log_depth)])
    )
  
  # --- Breakpoint estimation, only if non-linearity was found ---
  if (edf > 2 && delta_aic > 2 && requireNamespace("segmented", quietly = TRUE)) {
    library(segmented)
    
    cat("\n── Breakpoint estimation (segmented regression) ──\n\n")
    cat("(seed =", BREAKPOINT_SEED, "- fixed for reproducibility of segmented()'s\n")
    cat("internal bootstrap restarts)\n\n")
    
    set.seed(BREAKPOINT_SEED)
    seg_model <- tryCatch(
      segmented(model_linear, seg.Z = ~log_depth, npsi = 2),
      error = function(e) {
        cat("2-breakpoint model failed to converge, trying 1 breakpoint\n")
        set.seed(BREAKPOINT_SEED)
        segmented(model_linear, seg.Z = ~log_depth, npsi = 1)
      }
    )
    print(summary(seg_model))
    
    breakpoints_log <- seg_model$psi[, "Est."]
    breakpoints_depth <- round(10^breakpoints_log - 1, 1)
    breakpoint_ci <- confint(seg_model)
    
    cat("\nRaw breakpoint CI matrix (check column order matches Est./lower/upper):\n")
    print(breakpoint_ci)
    
    cat("\n── Estimated breakpoint(s), back-transformed to metres ──\n")
    breakpoint_table <- tibble()
    for (i in seq_along(breakpoints_depth)) {
      ci_log <- breakpoint_ci[i, ]
      ci_depth_lo <- round(10^ci_log[2] - 1, 1)
      ci_depth_hi <- round(10^ci_log[3] - 1, 1)
      cat(sprintf("Breakpoint %d: %.1f m (95%% CI: %.1f - %.1f m)\n",
                  i, breakpoints_depth[i], ci_depth_lo, ci_depth_hi))
      breakpoint_table <- bind_rows(breakpoint_table, tibble(
        breakpoint_number = i, depth_m = breakpoints_depth[i],
        ci_lower_m = ci_depth_lo, ci_upper_m = ci_depth_hi
      ))
    }
    write_csv(breakpoint_table, "depth_breakpoints.csv")
    cat("\nCompare to the ecological zone boundaries (200m, 1000m, 4000m).\n")
    final_breakpoints_m <- breakpoints_depth  # default: 2-breakpoint result, may be overridden below
    
    # --- Check whether a 3rd breakpoint further out improves the fit ---
    cat("\n── Checking whether a 3rd breakpoint improves the fit ──\n\n")
    cat("(seed =", BREAKPOINT_SEED, "again, for the same reproducibility reason)\n\n")
    
    set.seed(BREAKPOINT_SEED)
    seg_model_3 <- tryCatch(
      segmented(model_linear, seg.Z = ~log_depth, npsi = 3),
      error = function(e) {
        cat("3-breakpoint model failed to converge:", conditionMessage(e), "\n")
        NULL
      }
    )
    
    if (!is.null(seg_model_3)) {
      aic_2 <- AIC(seg_model)
      aic_3 <- AIC(seg_model_3)
      cat("2-breakpoint model AIC:", round(aic_2, 1), "\n")
      cat("3-breakpoint model AIC:", round(aic_3, 1), "\n")
      
      if (aic_3 < aic_2 - 2) {
        cat("\nAIC suggests the 3-breakpoint model fits better (\u0394AIC =",
            round(aic_2 - aic_3, 1), ") - but AIC alone is not enough here.\n")
        cat("Checking for degeneracy (breakpoints collapsing onto each other,\n")
        cat("which produces a spurious AIC improvement from local noise-fitting\n")
        cat("rather than a genuine third structural change):\n\n")
        
        breakpoints_log_3 <- seg_model_3$psi[, "Est."]
        
        min_gap <- min(diff(sort(breakpoints_log_3)))
        full_range <- diff(range(model_data$log_depth))
        gap_is_tiny <- min_gap < 0.05 * full_range
        
        seg_coefs <- summary(seg_model_3)$coefficients
        u_rows <- grepl("^U[0-9]+\\.log_depth$", rownames(seg_coefs))
        u_unstable <- any(abs(seg_coefs[u_rows, "Std. Error"]) >
                            10 * abs(seg_coefs[u_rows, "Estimate"]))
        
        if (gap_is_tiny || u_unstable) {
          cat("\u26A0\uFE0F  DEGENERACY DETECTED - the 3-breakpoint model is NOT reliable:\n")
          if (gap_is_tiny) {
            cat("  - Two breakpoints have collapsed onto nearly the same location\n")
            cat("    (smallest gap:", round(min_gap, 3), "on the log-depth scale,\n")
            cat("    versus a full data range of", round(full_range, 2), ").\n")
          }
          if (u_unstable) {
            cat("  - At least one jump-size coefficient has a standard error\n")
            cat("    more than 10x its estimate - the fit is numerically unstable.\n")
          }
          cat("\nREJECTING the 3-breakpoint model despite its lower AIC - this is\n")
          cat("overfitting to local noise, not a genuine third regime change.\n")
          cat("STICK WITH the 2-breakpoint result saved above (depth_breakpoints.csv)\n")
          cat("as the final answer.\n")
        } else {
          cat("No degeneracy detected - the 3-breakpoint model appears genuine.\n\n")
          
          print(summary(seg_model_3))
          
          breakpoints_depth_3 <- round(10^breakpoints_log_3 - 1, 1)
          breakpoint_ci_3 <- confint(seg_model_3)
          
          cat("\n── 3-breakpoint estimates, back-transformed to metres ──\n")
          breakpoint_table_3 <- tibble()
          for (i in seq_along(breakpoints_depth_3)) {
            ci_log <- breakpoint_ci_3[i, ]
            ci_depth_lo <- round(10^ci_log[2] - 1, 1)
            ci_depth_hi <- round(10^ci_log[3] - 1, 1)
            cat(sprintf("Breakpoint %d: %.1f m (95%% CI: %.1f - %.1f m)\n",
                        i, breakpoints_depth_3[i], ci_depth_lo, ci_depth_hi))
            breakpoint_table_3 <- bind_rows(breakpoint_table_3, tibble(
              breakpoint_number = i, depth_m = breakpoints_depth_3[i],
              ci_lower_m = ci_depth_lo, ci_upper_m = ci_depth_hi
            ))
          }
          write_csv(breakpoint_table_3, "depth_breakpoints_3psi.csv")
          cat("\nUse THESE 3 breakpoints instead of the 2-breakpoint set above -\n")
          cat("saved to depth_breakpoints_3psi.csv.\n")
          final_breakpoints_m <- breakpoints_depth_3  # override: 3-breakpoint result confirmed genuine
        }
        
      } else {
        cat("\nRESULT: the 3-breakpoint model does NOT fit meaningfully better\n")
        cat("(\u0394AIC =", round(aic_2 - aic_3, 1), ", below the >2 threshold) - the\n")
        cat("2-breakpoint model (already saved above) remains the best choice.\n")
      }
    }
  } else if (!requireNamespace("segmented", quietly = TRUE)) {
    cat("\nPackage 'segmented' not installed - run install.packages('segmented')\n")
    cat("to estimate the breakpoint location precisely.\n")
  }
  
} else {
  cat("Package 'mgcv' not installed - run install.packages('mgcv').\n")
}

# --- GAM figure, built HERE (after the breakpoint decision above is
# finalised) so the plotted vertical lines always match whatever was
# actually determined to be the correct number of breakpoints (2 or 3). ---
if (exists("pred_data") && !is.null(final_breakpoints_m)) {
  cat(sprintf("\n── Building GAM figure with %d confirmed breakpoint(s): %s ──\n",
              length(final_breakpoints_m), paste(final_breakpoints_m, collapse = ", ")))
  
  breakpoint_caption <- if (length(final_breakpoints_m) == 1) {
    sprintf("Dashed line marks the estimated breakpoint (%.1f m)", final_breakpoints_m)
  } else {
    sprintf("Dashed lines mark the %d estimated breakpoints (%s m)",
            length(final_breakpoints_m), paste(final_breakpoints_m, collapse = ", "))
  }
  
  p_gam <- pred_data %>%
    filter(depth_m > 0) %>%
    ggplot(aes(x = depth_m, y = fit * 100)) +
    geom_vline(xintercept = final_breakpoints_m, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100), fill = "#2a78d6", alpha = 0.2) +
    geom_line(color = "#2a78d6", linewidth = 1) +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 4000), labels = scales::comma) +
    labs(
      x = "Depth (m, log scale)", y = "Predicted probability of species-level ID (%)",
      title = "GAM-fitted depth relationship (Ophiuroidea, Flag 1 reference)",
      caption = breakpoint_caption
    ) +
    theme_minimal(base_size = 12)
  
  ggsave("figure_gam_depth_shape.png", p_gam, width = 8, height = 5, dpi = 300)
  cat("Figure saved: figure_gam_depth_shape.png\n")
} else if (exists("pred_data")) {
  p_gam <- pred_data %>%
    filter(depth_m > 0) %>%
    ggplot(aes(x = depth_m, y = fit * 100)) +
    geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100), fill = "#2a78d6", alpha = 0.2) +
    geom_line(color = "#2a78d6", linewidth = 1) +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 4000), labels = scales::comma) +
    labs(
      x = "Depth (m, log scale)", y = "Predicted probability of species-level ID (%)",
      title = "GAM-fitted depth relationship (Ophiuroidea, Flag 1 reference)"
    ) +
    theme_minimal(base_size = 12)
  
  ggsave("figure_gam_depth_shape.png", p_gam, width = 8, height = 5, dpi = 300)
  cat("Figure saved: figure_gam_depth_shape.png (no breakpoints estimated)\n")
}


# =============================================================================
# SECTION 3b: DOES THE DEPTH RELATIONSHIP'S SHAPE DIFFER BY CLASS?
# =============================================================================
cat("\n\n=== SECTION 3b: Class-specific depth relationship shape ===\n\n")

cat("── Sample size per class x depth zone ──\n")
print(table(model_data$echino_class, model_data$depth_zone))
cat("\nNote: Crinoidea and Echinoidea have ZERO records at Abyssal-\n")
cat("their per-class models below may fail or be unstable for this reason.\n\n")

# breakpoint_by_class is populated by the per-class loop below and then
# used to build class_breakpoints for the figure - initialise here so the
# figure code has something sensible to fall back on even if the loop
# produces nothing (e.g. package missing, or every class fails to fit).
breakpoint_by_class <- tibble()

if (requireNamespace("mgcv", quietly = TRUE)) {
  
  model_gam_byclass <- gam(
    species_level ~ echino_class + s(log_depth, by = echino_class, k = 4) + obs_quality_flag,
    data = model_data, family = binomial, method = "REML"
  )
  
  cat("── Per-class smooth term summary (k=4 constrained) ──\n\n")
  print(summary(model_gam_byclass))
  
  aic_shared_gam <- AIC(model_gam)
  aic_byclass_gam <- AIC(model_gam_byclass)
  cat("\nShared-shape GAM AIC:  ", round(aic_shared_gam, 1), "\n")
  cat("Class-specific GAM AIC:", round(aic_byclass_gam, 1), "\n")
  cat("\u0394AIC:", round(aic_shared_gam - aic_byclass_gam, 1), "in favour of class-specific model\n")
  cat("\nNOTE: before trusting this AIC comparison, inspect\n")
  cat("figure_gam_depth_shape_by_class.png for any class whose curve still\n")
  cat("oscillates implausibly even at k=4 - that class's AIC contribution is\n")
  cat("likely inflated by overfitting to sparse/heterogeneous data, not\n")
  cat("evidence of a genuine shape difference. Cross-check any such class\n")
  cat("against the per-class breakpoint table below, which uses a much\n")
  cat("simpler, less overfit-prone model.\n")
  
  depth_range_bc <- seq(min(model_data$log_depth), max(model_data$log_depth), length.out = 200)
  pred_grid <- expand_grid(
    log_depth = depth_range_bc, echino_class = levels(model_data$echino_class)
  ) %>%
    mutate(
      echino_class = factor(echino_class, levels = levels(model_data$echino_class)),
      obs_quality_flag = factor("1", levels = levels(model_data$obs_quality_flag))
    )
  pred_bc <- predict(model_gam_byclass, newdata = pred_grid, type = "response", se.fit = TRUE)
  pred_grid <- pred_grid %>%
    mutate(
      depth_m = 10^log_depth - 1, fit = pred_bc$fit,
      lower = pmax(0, pred_bc$fit - 1.96 * pred_bc$se.fit),
      upper = pmin(1, pred_bc$fit + 1.96 * pred_bc$se.fit)
    ) %>%
    filter(is.finite(depth_m), depth_m > 0)
  
  # --- Per-class breakpoint estimation (simpler, single-breakpoint model) ---
  # MOVED to run BEFORE the class-specific figure below, so the figure's
  # dashed reference lines are built FROM these results directly rather
  # than a hardcoded vector that goes stale every time the dataset changes.
  if (requireNamespace("segmented", quietly = TRUE)) {
    cat("\n── Per-class breakpoint estimation ──\n\n")
    
    classes <- levels(model_data$echino_class)
    
    for (cls in classes) {
      cat("── ", cls, " ──\n")
      d <- model_data %>% filter(echino_class == cls)
      
      m_linear_cls <- tryCatch(
        glm(species_level ~ log_depth + obs_quality_flag, data = d, family = binomial),
        error = function(e) NULL
      )
      if (is.null(m_linear_cls)) { cat("  Linear model failed to fit - skipping.\n\n"); next }
      
      set.seed(BREAKPOINT_SEED)
      seg_m <- tryCatch(
        segmented(m_linear_cls, seg.Z = ~log_depth, npsi = 1),
        error = function(e) { cat("  Segmented model failed to converge.\n"); NULL }
      )
      if (is.null(seg_m)) {
        cat("  No stable breakpoint found (likely sparse deep records).\n\n")
        next
      }
      
      bp_log <- seg_m$psi[, "Est."]
      bp_depth <- round(10^bp_log - 1, 1)
      ci_cls <- tryCatch(confint(seg_m), error = function(e) NULL)
      
      seg_coefs_cls <- summary(seg_m)$coefficients
      u_row_cls <- grepl("^U1\\.log_depth$", rownames(seg_coefs_cls))
      u_unstable_cls <- any(abs(seg_coefs_cls[u_row_cls, "Std. Error"]) >
                              10 * abs(seg_coefs_cls[u_row_cls, "Estimate"]))
      
      if (u_unstable_cls) {
        cat("  \u26A0\uFE0F Unstable jump coefficient - NOT reporting for", cls, ".\n\n")
        next
      }
      
      if (!is.null(ci_cls)) {
        ci_lo <- round(10^ci_cls[1, 2] - 1, 1)
        ci_hi <- round(10^ci_cls[1, 3] - 1, 1)
        cat(sprintf("  Breakpoint: %.1f m (95%% CI: %.1f - %.1f m)\n\n", bp_depth, ci_lo, ci_hi))
        breakpoint_by_class <- bind_rows(breakpoint_by_class, tibble(
          echino_class = cls, n = nrow(d), breakpoint_m = bp_depth,
          ci_lower_m = ci_lo, ci_upper_m = ci_hi
        ))
      }
    }
    
    cat("── Summary: per-class breakpoints ──\n")
    if (nrow(breakpoint_by_class) > 0) {
      print(breakpoint_by_class, n = Inf)
      write_csv(breakpoint_by_class, "depth_breakpoints_by_class.csv")
      cat("\nCompare to the pooled estimate. Overlapping CIs across classes\n")
      cat("indicate where the decline first steepens\n")
      cat("though the earlier interaction test showed the RATE of decline\n")
      cat("differs significantly by class.\n")
    } else {
      cat("No class produced a stable, trustworthy breakpoint estimate.\n")
    }
  } else {
    cat("Package 'segmented' not installed - run install.packages('segmented').\n")
  }
  
  # --- Class-specific figure, built AFTER the per-class breakpoints above,
  # using class_breakpoints DERIVED from breakpoint_by_class rather than a
  # hardcoded vector (this was the source of a stale-figure bug in an
  # earlier version of this script). ---
  class_breakpoints <- breakpoint_by_class %>%
    dplyr::select(echino_class, breakpoint_m)
  
  p_byclass <- pred_grid %>%
    ggplot(aes(x = depth_m, y = fit * 100, color = echino_class, fill = echino_class)) +
    { if (nrow(class_breakpoints) > 0)
      geom_vline(data = class_breakpoints,
                 aes(xintercept = breakpoint_m, colour = echino_class),
                 linetype = "dashed", linewidth = 0.4, show.legend = FALSE)
    } +
    geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 4000), labels = scales::comma) +
    scale_color_manual(values = c(
      "Asteroidea" = "#648FFF", "Crinoidea" = "#117733", "Echinoidea" = "#FE6100",
      "Holothuroidea" = "#785EF0", "Ophiuroidea" = "#DC267F"
    )) +
    scale_fill_manual(values = c(
      "Asteroidea" = "#648FFF", "Crinoidea" = "#117733", "Echinoidea" = "#FE6100",
      "Holothuroidea" = "#785EF0", "Ophiuroidea" = "#DC267F"
    )) +
    labs(x = "Depth (m, log scale)", y = "Predicted probability of species-level ID (%)",
         color = "Class", fill = "Class",
         title = "Class-specific GAM-fitted depth relationships (Flag 1 reference)",
         caption = if (nrow(class_breakpoints) > 0)
           "Dashed lines mark each class's own estimated breakpoint"
         else NULL) +
    theme_minimal(base_size = 12)
  
  ggsave("figure_gam_depth_shape_by_class.png", p_byclass, width = 9, height = 5.5, dpi = 300)
  cat("\nFigure saved: figure_gam_depth_shape_by_class.png\n")
  
} else {
  cat("Package 'mgcv' not installed - skipping class-specific shape check.\n")
}


# =============================================================================
# SECTION 4: ORDINAL LOGISTIC REGRESSION (full resolution hierarchy)
# Uses the log-depth specification justified in Section 3.
# =============================================================================
cat("\n\n=== SECTION 4: Ordinal logistic regression ===\n\n")

model_data_ord <- model_data %>%
  filter(rank6 %in% c("Species", "Genus", "Family", "Order", "Class")) %>%
  mutate(rank_ordered = factor(rank6, levels = c("Species", "Genus", "Family", "Order", "Class"),
                               ordered = TRUE))

cat("rank6 distribution in this subset (Phylum excluded):\n")
print(table(model_data_ord$rank6))
cat("\nRecords in ordinal model:", nrow(model_data_ord), "\n\n")

ordinal_model <- MASS::polr(
  rank_ordered ~ log_depth + echino_class + obs_quality_flag,
  data = model_data_ord, Hess = TRUE
)
print(summary(ordinal_model))

coef_table <- coef(summary(ordinal_model))
p_values <- pnorm(abs(coef_table[, "t value"]), lower.tail = FALSE) * 2
coef_table <- cbind(coef_table, p_value = p_values)
cat("\n── Coefficients with p-values ──\n")
print(round(coef_table, 4))

or_ci <- exp(cbind(OR = coef(ordinal_model), confint(ordinal_model)))
cat("\n── Odds ratios (95% CI) - odds of a WORSE resolution category ──\n")
print(round(or_ci, 3))

ordinal_results <- as_tibble(round(or_ci, 3), rownames = "term") %>%
  rename(ci_lower = `2.5 %`, ci_upper = `97.5 %`) %>%
  left_join(
    as_tibble(round(coef_table, 4), rownames = "term") %>% dplyr::select(term, p_value),
    by = "term"
  )
write_csv(ordinal_results, "ordinal_regression_taxonomic_completeness.csv")


# =============================================================================
# SECTION 4b: FULL PAIRWISE CLASS COMPARISONS
# =============================================================================
# Section 4's ordinal model only shows each class's comparison to
# Ophiuroidea (an arbitrary reference level) - e.g. it cannot directly
# answer "does Crinoidea differ from Holothuroidea?" without refitting.
# This extends the SAME fitted model to all C(5,2)=10 pairwise class
# comparisons directly, using estimated marginal means on the model's
# linear predictor scale (equivalent to comparing log-odds of worse
# resolution between every pair of classes, the same quantity already
# reported vs Ophiuroidea in Section 4), with Holm correction.
# =============================================================================
cat("\n\n=== SECTION 4b: Full pairwise class comparisons ===\n\n")
cat("Extends Section 4's ordinal model (each class vs Ophiuroidea only) to\n")
cat("all pairwise class comparisons - does the model support statements\n")
cat("like 'Crinoidea differs from Holothuroidea', not just each class vs\n")
cat("the reference?\n\n")

if (requireNamespace("emmeans", quietly = TRUE)) {
  library(emmeans)
  
  class_emm <- emmeans(ordinal_model, ~ echino_class)
  
  cat("── Estimated marginal means (linear predictor scale) ──\n\n")
  print(class_emm)
  cat("\nNOTE ON DIRECTION: in Section 4's model, a HIGHER coefficient means\n")
  cat("HIGHER odds of worse resolution. By the same convention, a HIGHER\n")
  cat("emmean above means WORSE resolution, a LOWER emmean means BETTER\n")
  cat("resolution. Verify this against the known direction before trusting\n")
  cat("any interpretation drawn from it.\n\n")
  
  cat("── All pairwise comparisons (Holm-corrected) ──\n\n")
  pairwise_contrasts <- pairs(class_emm, adjust = "holm")
  print(pairwise_contrasts)
  
  pairwise_df <- as_tibble(pairwise_contrasts) %>%
    mutate(
      OR = round(exp(estimate), 3),
      ci_lower = round(exp(estimate - 1.96 * SE), 3),
      ci_upper = round(exp(estimate + 1.96 * SE), 3),
      p_holm = round(p.value, 4)
    ) %>%
    dplyr::select(contrast, OR, ci_lower, ci_upper, p_holm)
  
  cat("\n── Pairwise odds ratios (95% CI), Holm-corrected ──\n\n")
  print(pairwise_df, n = Inf)
  write_csv(pairwise_df, "ordinal_pairwise_class_comparisons.csv")
  
  cat("\nRead the table above directly for which pairs are significant after\n")
  cat("Holm correction - do not assume a pattern without checking p_holm.\n")
  
} else {
  cat("Package 'emmeans' not installed - run install.packages('emmeans')\n")
  cat("to compute full pairwise class comparisons.\n")
}


# =============================================================================
# SHARED HELPER: separation-detection, used by Section 5 and 5c so both
# three apply IDENTICAL criteria instead of three independent copies that
# could quietly drift out of sync with each other over time.
#
# Flags a fitted binomial GLM as showing (quasi-)complete separation if
# EITHER:
#   (a) R itself raised its own "fitted probabilities numerically 0 or 1"
#       warning during fitting, or
#   (b) any coefficient's standard error is more than 10x its estimate
#       (the same threshold already used for the segmented-regression
#       degeneracy check in Section 3, kept consistent here)
# =============================================================================
fit_glm_flag_separation <- function(formula, data) {
  separation_warning <- FALSE
  m <- withCallingHandlers(
    glm(formula, data = data, family = binomial),
    warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w))) {
        separation_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  coefs <- summary(m)$coefficients
  se_unstable <- any(abs(coefs[, "Std. Error"]) > 10 * abs(coefs[, "Estimate"]), na.rm = TRUE)
  list(model = m, separation_detected = separation_warning || se_unstable)
}


# =============================================================================
# SECTION 5: PROPORTIONAL ODDS CHECK + cutpoint-specific models
# (the corrected result, given a violation is expected)
# =============================================================================
cat("\n\n=== SECTION 5: Proportional odds assumption check ===\n\n")
cat("H0: each predictor's effect is the SAME at every cutpoint\n\n")

if (requireNamespace("brant", quietly = TRUE)) {
  brant_result <- brant::brant(ordinal_model)
  print(brant_result)
  cat("\nA significant (p < .05) result means that predictor's effect genuinely\n")
  cat("differs across cutpoints - the single ordinal OR is an oversimplification\n")
  cat("for that term. See cutpoint-specific models below for the corrected view.\n")
} else {
  cat("Package 'brant' not installed - run install.packages('brant').\n")
  cat("Proceeding directly to cutpoint-specific models below, which are\n")
  cat("informative regardless of whether the formal test is run.\n")
}

cat("\n\n── Cutpoint-specific binary models ──\n")
cat("(eport these ALONGSIDE the overall ordinal OR if\n")
cat("proportional odds is violated - depth's effect size at each cutpoint\n")
cat("is itself the finding, not noise to average away)\n\n")

cutpoint_defs <- list(
  list(better = c("Species"), label = "Species | worse"),
  list(better = c("Species", "Genus"), label = "\u2264Genus | worse"),
  list(better = c("Species", "Genus", "Family"), label = "\u2264Family | worse"),
  list(better = c("Species", "Genus", "Family", "Order"), label = "\u2264Order | Class")
)

cutpoint_results <- map_dfr(cutpoint_defs, function(cp) {
  d <- model_data_ord %>% mutate(y = as.integer(rank_ordered %in% cp$better))
  
  fit <- fit_glm_flag_separation(y ~ log_depth + echino_class + obs_quality_flag, d)
  m <- fit$model
  coefs <- summary(m)$coefficients
  
  or_depth <- exp(-coef(m)["log_depth"])
  ci_std <- confint(m, "log_depth", trace = FALSE)
  
  result <- tibble(
    cutpoint = cp$label,
    n = nrow(d),
    OR_worse_resolution = round(or_depth, 3),
    ci_lower = round(exp(-ci_std[2]), 3),
    ci_upper = round(exp(-ci_std[1]), 3),
    p_value = signif(coefs["log_depth", "Pr(>|z|)"], 3),
    separation_flag = fit$separation_detected,
    method = "standard glm"
  )
  
  if (fit$separation_detected && requireNamespace("logistf", quietly = TRUE)) {
    cat(sprintf("\n\u26A0\uFE0F  Separation detected at cutpoint '%s' - fitting Firth's penalized\n", cp$label))
    cat("   logistic regression for comparison...\n")
    
    m_firth <- logistf::logistf(y ~ log_depth + echino_class + obs_quality_flag, data = d)
    or_depth_firth <- exp(-coef(m_firth)["log_depth"])
    ci_firth <- confint(m_firth)
    
    cat(sprintf("   Standard glm:  log_depth OR = %.3f\n", or_depth))
    cat(sprintf("   Firth's:       log_depth OR = %.3f\n", or_depth_firth))
    
    result_firth <- tibble(
      cutpoint = cp$label,
      n = nrow(d),
      OR_worse_resolution = round(or_depth_firth, 3),
      ci_lower = round(exp(-ci_firth["log_depth", 2]), 3),
      ci_upper = round(exp(-ci_firth["log_depth", 1]), 3),
      p_value = signif(unname(m_firth$prob["log_depth"]), 3),
      separation_flag = TRUE,
      method = "Firth's penalized (recommended)"
    )
    result <- bind_rows(result, result_firth)
  } else if (fit$separation_detected) {
    cat(sprintf("\n\u26A0\uFE0F  Separation detected at cutpoint '%s' but package 'logistf' not\n", cp$label))
    cat("   installed - run install.packages('logistf') for the corrected estimate.\n")
  }
  
  result
})

cat("── OR (per 10-fold depth increase) of being on the WORSE side of each cutpoint ──\n")
print(cutpoint_results, n = Inf)
write_csv(cutpoint_results, "cutpoint_specific_odds_ratios.csv")

n_separation <- sum(cutpoint_results$separation_flag & cutpoint_results$method == "standard glm")
if (n_separation > 0) {
  cat(sprintf("\n%d of 4 cutpoints showed signs of quasi-complete separation with\n", n_separation))
  cat("this dataset - use the 'Firth's penalized' row (not 'standard glm') for\n")
  cat("those cutpoints when reporting. Compare the two OR values above to\n")
  cat("confirm whether depth's estimate is materially affected or stable.\n")
} else {
  cat("\nNo separation detected at any cutpoint with this dataset - standard\n")
  cat("glm estimates are valid as-is, no correction needed.\n")
}
cat("\nIf these ORs increase from top to bottom, depth has a progressively\n")
cat("STRONGER effect on the worst resolution outcomes - the concrete finding\n")
cat("behind the proportional-odds violation.\n")


# =============================================================================
# SECTION 5b: PAIRWISE CLASS COMPARISONS, AT EACH CUTPOINT SEPARATELY
# =============================================================================
cat("\n\n=== SECTION 5b: pairwise class comparisons, by cutpoint ===\n\n")
cat("Checks WHERE (at which cutpoint) the class differences found in\n")
cat("Section 4b are concentrated, rather than relying on a single value\n")
cat("averaged across all four cutpoints.\n\n")

if (requireNamespace("emmeans", quietly = TRUE)) {
  library(emmeans)
  
  cutpoint_class_pairwise <- map_dfr(cutpoint_defs, function(cp) {
    d <- model_data_ord %>% mutate(y = as.integer(rank_ordered %in% cp$better))
    m <- glm(y ~ log_depth + echino_class + obs_quality_flag, data = d, family = binomial)
    
    class_emm_cp <- emmeans(m, ~ echino_class)
    pairwise_cp <- pairs(class_emm_cp, adjust = "holm")
    
    as_tibble(pairwise_cp) %>%
      mutate(
        cutpoint = cp$label,
        OR = round(exp(estimate), 3),
        ci_lower = round(exp(estimate - 1.96 * SE), 3),
        ci_upper = round(exp(estimate + 1.96 * SE), 3),
        p_holm = round(p.value, 4)
      ) %>%
      dplyr::select(cutpoint, contrast, OR, ci_lower, ci_upper, p_holm)
  })
  
  cat("── Pairwise class comparisons at each cutpoint (Holm-corrected within each cutpoint) ──\n\n")
  print(cutpoint_class_pairwise, n = Inf)
  write_csv(cutpoint_class_pairwise, "cutpoint_specific_pairwise_class_comparisons.csv")
  
  cat("\n── Focus: Holothuroidea vs each other class, across cutpoints ──\n\n")
  holothuroidea_by_cutpoint <- cutpoint_class_pairwise %>%
    filter(str_detect(contrast, "Holothuroidea"))
  print(holothuroidea_by_cutpoint, n = Inf)
  
  cat("\n── Focus: Asteroidea vs Echinoidea, across cutpoints ──\n\n")
  ast_echi_by_cutpoint <- cutpoint_class_pairwise %>%
    filter(str_detect(contrast, "Asteroidea") & str_detect(contrast, "Echinoidea"))
  print(ast_echi_by_cutpoint, n = Inf)
  
  cat("\nINTERPRETATION: read the actual numbers above before writing any\n")
  cat("specific claim - do not assume a pattern without checking it.\n")
  
  cutpoint_order <- c("Species | worse", "\u2264Genus | worse", "\u2264Family | worse", "\u2264Order | Class")
  
  p_cutpoint_class <- cutpoint_class_pairwise %>%
    filter(str_starts(contrast, "Ophiuroidea")) %>%
    mutate(
      other_class = str_remove(contrast, "^Ophiuroidea - "),
      cutpoint = factor(cutpoint, levels = cutpoint_order)
    ) %>%
    ggplot(aes(x = cutpoint, y = OR, color = other_class, group = other_class)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    scale_y_log10() +
    scale_color_manual(values = c(
      "Asteroidea" = "#648FFF", "Crinoidea" = "#117733", "Echinoidea" = "#FE6100",
      "Holothuroidea" = "#785EF0"
    )) +
    labs(x = "Cutpoint", y = "OR vs Ophiuroidea (log scale)", color = "Class",
         title = "Class-specific odds ratios across taxonomic cutpoints.",
         subtitle = "OR > 1 = worse resolution than Ophiuroidea at that cutpoint") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  
  ggsave("figure_cutpoint_class_comparisons.png", p_cutpoint_class, width = 9, height = 6, dpi = 300)
  cat("\nFigure saved: figure_cutpoint_class_comparisons.png\n")
  
  cat("\n=== Section 5b complete ===\n")
  
} else {
  cat("Package 'emmeans' not installed - skipping cutpoint-specific pairwise\n")
  cat("class comparisons.\n")
}


# =============================================================================
# SECTION 5c: DOES DEPTH x CLASS INTERACT, AT EACH CUTPOINT?
# =============================================================================
# Section 5b's class comparisons all CONTROL for depth - they show which
# classes differ overall, but not whether the SIZE of a class's
# disadvantage itself depends on depth. This adds a log_depth:echino_class
# interaction to each of the four cutpoint-specific models and tests it
# via LRT. Both the main-effects model AND the interaction model are now
# checked for separation using the shared helper above: sparse depth x
# class cells (e.g. Crinoidea/Echinoidea have ZERO Abyssal records - see
# Section 3b) are exactly the setup that produces quasi-complete
# separation once an interaction term is added, and this would otherwise
# make the LRT p-value unreliable without any visible warning.
# =============================================================================
cat("\n\n=== SECTION 5c: depth x class interaction, at each cutpoint ===\n\n")
cat("H0: depth's effect on reaching/failing this cutpoint is the SAME\n")
cat("    across all echinoderm classes (no interaction)\n\n")

interaction_by_cutpoint <- map_dfr(cutpoint_defs, function(cp) {
  d <- model_data_ord %>% mutate(y = as.integer(rank_ordered %in% cp$better))
  
  fit_main <- fit_glm_flag_separation(y ~ log_depth + echino_class + obs_quality_flag, d)
  fit_interact <- fit_glm_flag_separation(y ~ log_depth * echino_class + obs_quality_flag, d)
  
  if (fit_interact$separation_detected) {
    cat(sprintf("\u26A0\uFE0F  Separation detected in the INTERACTION model at cutpoint '%s' -\n", cp$label))
    cat("   likely a sparse depth x class cell (e.g. Crinoidea/Echinoidea at\n")
    cat("   Abyssal). The LRT p-value below may be unreliable for this\n")
    cat("   cutpoint - treat with caution; this is why the depth-stratified\n")
    cat("   breakdown (previously Section 5d) was dropped rather than pursued.\n\n")
  }
  
  lrt <- anova(fit_main$model, fit_interact$model, test = "LRT")
  
  tibble(
    cutpoint = cp$label,
    n = nrow(d),
    deviance_diff = round(lrt$Deviance[2], 2),
    df = lrt$Df[2],
    p_value = signif(lrt$`Pr(>Chi)`[2], 3),
    main_model_separation_flag = fit_main$separation_detected,
    interaction_model_separation_flag = fit_interact$separation_detected
  )
})

cat("── Likelihood ratio test: does adding depth x class improve the fit, at each cutpoint? ──\n\n")
print(interaction_by_cutpoint, n = Inf)
write_csv(interaction_by_cutpoint, "depth_class_interaction_by_cutpoint.csv")

n_sig_interaction <- sum(interaction_by_cutpoint$p_value < 0.05)
n_interaction_separation <- sum(interaction_by_cutpoint$interaction_model_separation_flag)

cat(sprintf("\n%d of %d cutpoints show a significant depth x class interaction (p<.05)\n",
            n_sig_interaction, nrow(interaction_by_cutpoint)))
if (n_interaction_separation > 0) {
  cat(sprintf("%d of %d cutpoints' interaction models show signs of separation -\n",
              n_interaction_separation, nrow(interaction_by_cutpoint)))
  cat("treat the LRT result for those specific cutpoints as provisional.\n")
}

if (n_sig_interaction > 0) {
  cat("\nAt cutpoint(s) with a significant interaction, the SIZE of each\n")
  cat("class's disadvantage genuinely changes with depth - the Section 5b\n")
  cat("comparisons (which hold depth constant) may understate or overstate\n")
  cat("the gap at specific depths.\n")
} else {
  cat("\nNo cutpoint shows a significant depth x class interaction - the\n")
  cat("class differences reported in Section 5b are consistent across\n")
  cat("depth, so those comparisons can be reported as depth-independent.\n")
}

cat("\nSTOPPING HERE rather than unpacking this interaction at specific\n")
cat("representative depths (an earlier version of this script did so in a\n")
cat("Section 5d). That breakdown required splitting the data into\n")
cat("depth x class x cutpoint cells, and several classes are too sparse in\n")
cat("deep water to support it reliably - see Section 3b (Crinoidea and\n")
cat("Echinoidea have ZERO Abyssal records) and this section's own\n")
cat("separation flags above. The resulting depth-stratified pairwise ORs\n")
cat("were technically computable but rested on very thin per-cell sample\n")
cat("sizes (e.g. as few as ~49 Crinoidea records near 1000m) and required\n")
cat("separation corrections for roughly half the cutpoints - more\n")
cat("machinery than the data could cleanly support. The omnibus finding\n")
cat("above (depth x class interaction is significant at every cutpoint) is\n")
cat("the appropriate level of claim for this dataset: report THAT a\n")
cat("class's resolution disadvantage changes with depth, not precise ORs\n")
cat("for specific depth x class x cutpoint combinations.\n")

cat("\n=== Section 5c complete ===\n")

# =============================================================================
# SECTION 6: Deep-only GAM Analysis 
# =============================================================================

wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n
  c0 <- (p + z^2 / (2 * n)) / d
  h  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  tibble(lo = c0 - h, hi = c0 + h)
}

# ---- 1. Where are the deep records? Finer bands, with n and CIs ----------
deep_bins <- model_data %>%
  filter(depth_median >= 200) %>%
  mutate(bin = cut(depth_median, c(200, 500, 1000, 1500, 2000, 2500, 3000, 4000, Inf),
                   right = FALSE,
                   labels = c("200–500", "500–1000", "1000–1500", "1500–2000",
                              "2000–2500", "2500–3000", "3000–4000", ">4000")))
deep_bins %>%
  group_by(bin) %>% summarise(n = n(), x = sum(species_level), .groups = "drop") %>%
  mutate(pct_species = round(100 * x / n, 1), wilson(x, n)) %>%
  print()

# ---- 2. GAM on deep records only (≥ 200 m), so the smooth focuses on the slope ----
deep <- model_data %>% filter(depth_median >= 200)
m_deep <- gam(species_level ~ s(log_depth, k = 5) + echino_class + obs_quality_flag,
              data = deep, family = binomial, method = "REML")
summary(m_deep)

nd <- tibble(log_depth = seq(min(deep$log_depth), max(deep$log_depth), length.out = 200),
             echino_class = factor("Ophiuroidea", levels = levels(model_data$echino_class)),
             obs_quality_flag = factor("1", levels = levels(model_data$obs_quality_flag)))
p <- predict(m_deep, nd, type = "link", se.fit = TRUE)
curve_deep <- nd %>% mutate(depth_m = 10^log_depth - 1, fit = plogis(p$fit),
                            lo = plogis(p$fit - 1.96 * p$se.fit), hi = plogis(p$fit + 1.96 * p$se.fit))

ggplot(curve_deep, aes(depth_m, fit * 100)) +
  geom_ribbon(aes(ymin = lo * 100, ymax = hi * 100), fill = "#2a78d6", alpha = 0.2) +
  geom_line(colour = "#2a78d6", linewidth = 1) +
  geom_rug(data = deep, aes(x = depth_median, y = NULL), sides = "b", alpha = 0.1) +
  geom_vline(xintercept = c(1000, 4000), linetype = "dashed", colour = "grey50") +
  labs(x = "Depth (m)", y = "Predicted probability of species-level ID (%)",
       title = "GAM-fitted depth relationship, records ≥200 m (Ophiuroidea, Flag 1 reference)",
       caption = "Dashed lines mark the 1,000 m and 4,000 m zone boundaries") +
  theme_minimal(base_size = 14)
ggsave("figure_gam_deep_only.png", width = 9, height = 5.5, dpi = 300)

# ---- 3. Do the zones capture the pattern, or does the smooth fit better? ----
m_zone   <- gam(species_level ~ factor(depth_zone) + echino_class + obs_quality_flag,
                data = model_data, family = binomial, method = "REML")
m_smooth <- gam(species_level ~ s(log_depth) + echino_class + obs_quality_flag,
                data = model_data, family = binomial, method = "REML")
AIC(m_zone, m_smooth)

# Lower Slope point test
m_zone_deep <- gam(species_level ~ factor(depth_zone) + echino_class + obs_quality_flag,
                   data = deep, family = binomial, method = "REML")
AIC(m_zone_deep, m_deep)

# =============================================================================
# SECTION 7: SENSITIVITY TO SOURCE / INSTITUTION STRUCTURE
# Records are nested within contributing datasets and institutions. Source is
# treated as a NUISANCE term (random intercept / sensitivity), not as a research
# question. Run once with GROUP_VAR = "primary_source", then again with
# GROUP_VAR = "institutionCode", and compare.
# =============================================================================
# install.packages(c("lme4", "ordinal", "sandwich", "lmtest", "forcats"))  # if needed
library(lme4); library(ordinal); library(sandwich); library(lmtest)

GROUP_VAR <- "merged"      # options: "primary_source", "institutionCode", "merged"

add_group <- function(df) {
  df %>% mutate(grp = case_when(
    GROUP_VAR == "institutionCode" ~ as.character(forcats::fct_lump_min(
      factor(coalesce(institutionCode, "NA")), min = 50, other_level = "Other")),
    GROUP_VAR == "merged" ~ as.character(forcats::fct_collapse(factor(primary_source),
                                                               GBRSBD = c("ALA_GBRSBD", "CSIRO_GBRSBD"),
                                                               Other  = c("iDigBio", "GBIF_echinodermata"))),
    TRUE ~ primary_source),
    grp = factor(grp))
}
model_data_src     <- add_group(model_data)
model_data_ord_src <- add_group(model_data_ord)

cat("\n=== SECTION 7: grouping variable =", GROUP_VAR, "===\n")
cat("Groups in models:", nlevels(droplevels(model_data_src$grp)), "\n")
print(sort(table(model_data_src$grp), decreasing = TRUE))

# ---- 7A. Pooled GAM with a random intercept for group ------------------------
m_gam_re <- gam(species_level ~ s(log_depth, k = 10) + echino_class + obs_quality_flag +
                  s(grp, bs = "re"),
                data = model_data_src, family = binomial, method = "REML")
m_lin_re <- gam(species_level ~ log_depth + echino_class + obs_quality_flag +
                  s(grp, bs = "re"),
                data = model_data_src, family = binomial, method = "REML")
cat("\n── 7A. Pooled GAM with random intercept ──\n")
print(summary(m_gam_re))
cat("Delta AIC vs linear (with random effect):",
    round(AIC(m_lin_re) - AIC(m_gam_re), 1), "\n")

nd <- tibble(
  log_depth = seq(min(model_data_src$log_depth), max(model_data_src$log_depth), length.out = 300),
  echino_class = factor("Ophiuroidea", levels = levels(model_data_src$echino_class)),
  obs_quality_flag = factor("1", levels = levels(model_data_src$obs_quality_flag)),
  grp = factor(levels(model_data_src$grp)[1], levels = levels(model_data_src$grp)))
pr <- predict(m_gam_re, nd, type = "link", se.fit = TRUE, exclude = "s(grp)")
curve_re <- nd %>% mutate(depth_m = 10^log_depth - 1, fit = plogis(pr$fit),
                          lo = plogis(pr$fit - 1.96 * pr$se.fit),
                          hi = plogis(pr$fit + 1.96 * pr$se.fit)) %>% filter(depth_m > 0)
p_re <- ggplot(curve_re, aes(depth_m, fit * 100)) +
  geom_vline(xintercept = c(9.0, 22.2, 77.1), linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo * 100, ymax = hi * 100), fill = "#2a78d6", alpha = 0.2) +
  geom_line(colour = "#2a78d6", linewidth = 1) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 4000), labels = scales::comma) +
  labs(x = "Depth (m, log scale)", y = "Predicted probability of species-level ID (%)",
       title = "Pooled GAM with random intercept for source (Ophiuroidea, Flag 1)") +
  theme_minimal(base_size = 12)
ggsave(paste0("figure_gam_depth_shape_re_", GROUP_VAR, ".png"), p_re, width = 8, height = 5, dpi = 300)

# ---- 7B. Ordinal model with random intercept ---------------------------------
cat("\n── 7B. Ordinal model (clmm) with random intercept ──\n")
m_clmm <- clmm(rank_ordered ~ log_depth + echino_class + obs_quality_flag + (1 | grp),
               data = model_data_ord_src, Hess = TRUE)
print(summary(m_clmm))
cf <- coef(summary(m_clmm))
cf <- cf[!grepl("\\|", rownames(cf)), , drop = FALSE]      # drop threshold rows
clmm_or <- tibble(term = rownames(cf), OR = exp(cf[, "Estimate"]),
                  ci_lower = exp(cf[, "Estimate"] - 1.96 * cf[, "Std. Error"]),
                  ci_upper = exp(cf[, "Estimate"] + 1.96 * cf[, "Std. Error"]),
                  p_value = cf[, "Pr(>|z|)"]) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
print(clmm_or)
write_csv(clmm_or, paste0("ordinal_clmm_random_intercept_", GROUP_VAR, ".csv"))

# ---- 7C. Cutpoint-specific models with random intercept ----------------------
cat("\n── 7C. Cutpoint-specific models with random intercept ──\n")
cutpoint_re <- map_dfr(cutpoint_defs, function(cp) {
  d <- model_data_ord_src %>% mutate(y = as.integer(rank_ordered %in% cp$better))
  m <- glmer(y ~ log_depth + echino_class + obs_quality_flag + (1 | grp), data = d,
             family = binomial,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  s <- summary(m)$coefficients
  b <- s["log_depth", "Estimate"]; se <- s["log_depth", "Std. Error"]
  tibble(cutpoint = cp$label, n = nrow(d),
         OR_worse_resolution = round(exp(-b), 2),
         ci_lower = round(exp(-(b + 1.96 * se)), 2),
         ci_upper = round(exp(-(b - 1.96 * se)), 2),
         group_sd = round(as.data.frame(VarCorr(m))$sdcor[1], 2))
})
print(cutpoint_re)
write_csv(cutpoint_re, paste0("cutpoint_random_intercept_", GROUP_VAR, ".csv"))

# ---- 7D. Cluster-robust SEs (cluster = institution), species-level binary -----
cat("\n── 7D. Cluster-robust SEs, clustered by institution ──\n")
m_glm <- glm(species_level ~ log_depth + echino_class + obs_quality_flag,
             data = model_data_src, family = binomial)
inst_vec <- coalesce(model_data_src$institutionCode, "NA")
vc <- vcovCL(m_glm, cluster = inst_vec)
print(coeftest(m_glm, vcov. = vc))
b <- coef(m_glm)["log_depth"]; se <- sqrt(vc["log_depth", "log_depth"])
cat(sprintf("Depth OR (per 10-fold, worse ID): %.2f (95%% CI %.2f-%.2f), %d institutions\n",
            exp(-b), exp(-(b + 1.96 * se)), exp(-(b - 1.96 * se)), n_distinct(inst_vec)))

# ---- 7E. Leave-one-source-out -------------------------------------------------
cat("\n── 7E. Leave-one-source-out ──\n")
loso_fit <- function(dd, label) {
  dd <- droplevels(dd)
  m  <- tryCatch(glm(species_level ~ log_depth + echino_class + obs_quality_flag,
                     data = dd, family = binomial), error = function(e) NULL)
  mo <- tryCatch(MASS::polr(rank_ordered ~ log_depth + echino_class + obs_quality_flag,
                            data = dd, Hess = TRUE), error = function(e) NULL)
  g <- function(model, nm, f) if (is.null(model) || !nm %in% names(coef(model))) NA_real_ else f(coef(model)[[nm]])
  tibble(dropped = label, n = nrow(dd),
         species_depth_OR = round(g(m,  "log_depth", function(x) exp(-x)), 2),
         ordinal_depth_OR = round(g(mo, "log_depth", exp), 2),
         holothuroidea_species_coef = round(g(m, "echino_classHolothuroidea", identity), 2),
         holothuroidea_ordinal_OR   = round(g(mo, "echino_classHolothuroidea", exp), 2))
}
big <- names(which(table(model_data_ord_src$grp) > 200))
loso <- bind_rows(
  loso_fit(model_data_ord_src, "none"),
  map_dfr(big, ~ loso_fit(filter(model_data_ord_src, as.character(grp) != .x), .x)),
  loso_fit(filter(model_data_ord_src, !primary_source %in% c("ALA_GBRSBD", "CSIRO_GBRSBD")),
           "both GBRSBD sources"))
print(loso, n = Inf)
write_csv(loso, paste0("leave_one_source_out_", GROUP_VAR, ".csv"))

# ---- 7F. Is the shallow shape (breakpoints) driven by GBRSBD? -----------------
cat("\n── 7F. Shallow-water shape with and without GBRSBD ──\n")
no_gb <- model_data_src %>% filter(!primary_source %in% c("ALA_GBRSBD", "CSIRO_GBRSBD"))
bins <- c(-0.01, 5, 9, 22, 40, 77, 120, 200)
bin_tab <- bind_rows(all_sources = model_data_src, without_GBRSBD = no_gb, .id = "dataset") %>%
  filter(depth_median <= 200) %>%
  mutate(bin = cut(depth_median, bins)) %>%
  group_by(dataset, bin) %>%
  summarise(n = n(), pct_species = round(100 * mean(species_level), 1), .groups = "drop") %>%
  pivot_wider(names_from = dataset, values_from = c(n, pct_species))
print(bin_tab, n = Inf)

lin_nogb <- glm(species_level ~ log_depth + echino_class + obs_quality_flag,
                data = no_gb, family = binomial)
for (k in 1:3) {
  set.seed(BREAKPOINT_SEED)
  sm <- tryCatch(segmented(lin_nogb, seg.Z = ~log_depth, npsi = k), error = function(e) NULL)
  if (is.null(sm)) { cat("npsi =", k, "without GBRSBD: did not converge\n"); next }
  cat(sprintf("npsi = %d without GBRSBD: AIC = %.1f; breakpoints (m): %s\n", k, AIC(sm),
              paste(round(10^sm$psi[, "Est."] - 1, 1), collapse = ", ")))
}

# ---- 7H. Model specification checks requested in examiner comments ----------
# (a) package versions, (b) zero-depth records,
# (d) zero-depth sensitivity of the breakpoints, (e) starting-value
# sensitivity of the breakpoints. Uses model_data / BREAKPOINT_SEED from
# Section 3.
cat("\n── 7H. Model specification checks ──\n")

# (a) Versions to report in Methods
cat("mgcv", as.character(packageVersion("mgcv")),
    "| segmented", as.character(packageVersion("segmented")), "\n")

# (b) Zero-depth records: how many, where from, how well identified
n_zero <- sum(model_data$depth_median == 0)
cat(sprintf("Records at 0 m: %d of %d (%.1f%%); %.1f%% species-level\n",
            n_zero, nrow(model_data), 100 * n_zero / nrow(model_data),
            100 * mean(model_data$species_level[model_data$depth_median == 0])))
cat("Records with 0 < depth < 1 m:",
    sum(model_data$depth_median > 0 & model_data$depth_median < 1), "\n")
print(table(model_data$primary_source[model_data$depth_median == 0]))

# Helper: AIC recomputed independently of segmented(). On some systems
# segmented() reports a deviance/AIC that does not match its own coefficients
# (e.g. a 3-breakpoint AIC higher than the linear model, which is impossible
# for a nested model). Here the breakpoints are held at segmented()'s estimates
# and an ordinary binomial GLM is refitted with the hinge terms (x - psi)+;
# 2 df are added per breakpoint for estimating psi.
refit_aic <- function(seg, d) {
  psi <- seg$psi[, "Est."]
  H <- sapply(psi, function(p) pmax(d$log_depth - p, 0))
  m <- glm(species_level ~ log_depth + H + echino_class + obs_quality_flag,
           data = d, family = binomial, control = glm.control(maxit = 100))
  AIC(m) + 2 * length(psi)
}

# Helper: fit 2- and 3-breakpoint models and return back-transformed results
seg_summary <- function(d, label, psi_start = NULL) {
  lin <- glm(species_level ~ log_depth + echino_class + obs_quality_flag,
             data = d, family = binomial)
  out <- tibble()
  for (k in 2:3) {
    set.seed(BREAKPOINT_SEED)
    s <- tryCatch(
      if (is.null(psi_start)) segmented(lin, seg.Z = ~log_depth, npsi = k)
      else segmented(lin, seg.Z = ~log_depth, psi = log10(psi_start + 1)),
      error = function(e) NULL)
    if (is.null(s)) next
    ci <- confint(s)  # delta method, log10 scale
    out <- bind_rows(out, tibble(
      analysis = label, n = nrow(d), n_breakpoints = nrow(s$psi),
      AIC_linear = round(AIC(lin), 1),
      AIC_segmented_reported = round(AIC(s), 1),
      AIC_refit = round(refit_aic(s, d), 1),
      breakpoints_m = paste(sprintf("%.1f (%.1f-%.1f)", 10^ci[, 1] - 1, 10^ci[, 2] - 1,
                                    10^ci[, 3] - 1), collapse = "; ")))
    if (!is.null(psi_start)) break
  }
  out
}

# Check of the main Section 3 models (seg_model = 2 bp, seg_model_3 = 3 bp):
# if AIC_segmented_reported > AIC of the linear model, the reported value is wrong
# and AIC_refit should be used for the 2- vs 3-breakpoint comparison (Appendix K).
main_aic_check <- tibble(
  model = c("linear", "2 breakpoints", "3 breakpoints"),
  AIC_reported = round(c(AIC(model_linear), AIC(seg_model), AIC(seg_model_3)), 1),
  AIC_refit = round(c(AIC(model_linear), refit_aic(seg_model, model_data),
                      refit_aic(seg_model_3, model_data)), 1)
)
cat("\nAIC check for the main breakpoint models:\n")
print(main_aic_check)

# (d) Zero-depth sensitivity
zero_sens <- bind_rows(
  seg_summary(model_data, "Main model"),
  seg_summary(filter(model_data, depth_median > 0), "Excluding 0 m records"),
  seg_summary(mutate(model_data, log_depth = if_else(depth_median == 0, log10(1.5), log_depth)),
              "0 m recoded as 0.5 m")
)
print(zero_sens, width = Inf)
write_csv(zero_sens, "breakpoint_zero_depth_sensitivity.csv")

# (e) Starting-value sensitivity of the 3-breakpoint model
starts <- list(c(5, 15, 40), c(9, 22, 77), c(15, 60, 150), c(20, 100, 500),
               c(10, 50, 200), c(30, 200, 1000))
start_sens <- map_dfr(starts, function(st) {
  r <- seg_summary(model_data, paste("start:", paste(st, collapse = "/")), psi_start = st)
  r
})
print(arrange(start_sens, AIC_refit), width = Inf)
write_csv(start_sens, "breakpoint_starting_value_sensitivity.csv")

# ---- 7G. Descriptive: % species-level by source x class ----------------------
cat("\n── 7G. % species-level (n) by source and class ──\n")
src_class <- model_data_src %>%
  group_by(primary_source, echino_class) %>%
  summarise(n = n(), pct = round(100 * mean(species_level), 1), .groups = "drop") %>%
  mutate(cell = paste0(pct, " (", n, ")")) %>%
  dplyr::select(primary_source, echino_class, cell) %>%
  pivot_wider(names_from = echino_class, values_from = cell, values_fill = "-")
print(src_class, n = Inf)
write_csv(src_class, "species_level_by_source_and_class.csv")

# Species-level: source as a fixed effect instead of random
m_fe <- glm(species_level ~ log_depth + echino_class + obs_quality_flag + primary_source,
            data = model_data, family = binomial)
b  <- coef(m_fe)["log_depth"]
se <- summary(m_fe)$coefficients["log_depth", "Std. Error"]
c(OR = exp(-b), lower = exp(-(b + 1.96 * se)), upper = exp(-(b - 1.96 * se)))

cat("\n=== Section 7 complete (", GROUP_VAR, ") ===\n")

cat("\n\n=== Taxonomic completeness analysis complete ===\n")
cat("Outputs: table1_rank_by_depthbin.csv, figure_rank_composition_by_depth.png,\n")
cat("figure_gam_depth_shape.png, depth_breakpoints.csv,\n")
cat("figure_gam_depth_shape_by_class.png, depth_breakpoints_by_class.csv,\n")
cat("ordinal_regression_taxonomic_completeness.csv, cutpoint_specific_odds_ratios.csv,\n")
cat("ordinal_pairwise_class_comparisons.csv, cutpoint_specific_pairwise_class_comparisons.csv,\n")
cat("figure_cutpoint_class_comparisons.png, depth_class_interaction_by_cutpoint.csv\n")