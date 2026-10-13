# =====================================================================
# 02  Figure 2 (Bray–Curtis PCoA + day-to-day dissimilarity), PERMANOVA,
#     per-run and windowed PERMANOVA, Figure S3 (Jaccard), stabilised average
# Chapter 3 — "Stability as a dynamic process" (J Appl Microbiol submission)
# Inputs : data/microbiome/TR*/alldaysTR*.csv, alldaysTR*m.csv
# Outputs: figures/Figure2.pdf, figures/FigureS3_Jaccard.pdf;
#          PERMANOVA R2 (Results 3.2), windowed PERMANOVA p, 0.238 +/- 0.151
# NOTE   : input files hold read counts; converted to relative abundance
#          before Bray–Curtis (Panel A, Panel B, all PERMANOVAs).
# CHECK  : Methods say PERMANOVA included compartment, timepoint and run;
#          the model below is Stage + Run.
# Run from the repository root:  source("R/<this file>")
# =====================================================================

# ---- 1. Figure 2 -----------------------------------------------------
base_dir <- "data/microbiome"
out_dir  <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(dplyr)
library(tidyr)
library(readr)
library(vegan)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

runs <- c("TR1", "TR2", "TR3", "TR4")

# ---------------------------------------------------------------------
# 1. LOAD RAW MICROBIOME DATA (genus-level relative abundance + metadata)
# ---------------------------------------------------------------------

microbiome_list <- list()

for (run in runs) {
  taxa_file <- file.path(base_dir, run, paste0("alldays", run, ".csv"))
  meta_file <- file.path(base_dir, run, paste0("alldays", run, "m.csv"))

  taxa <- read.csv(taxa_file, stringsAsFactors = FALSE)
  meta <- read_csv(meta_file, show_col_types = FALSE)

  merged <- meta %>%
    left_join(taxa, by = "id") %>%
    separate(Sample.name,
             into = c("Run2", "Day", "Compartment"),
             sep = "_") %>%
    mutate(
      Day = as.numeric(Day),
      Run = run,
      SampleID = paste(run, Day, Compartment, sep = "_")
    ) %>%
    dplyr::select(-Run2)

  microbiome_list[[run]] <- merged
}

microbiome_all <- bind_rows(microbiome_list)

# ---------------------------------------------------------------------
# 2. FILTER TO DAYS 1-10 (all sampling days)
# ---------------------------------------------------------------------

microbiome_filtered <- microbiome_all %>%
  filter(Day %in% 1:10)

# ---------------------------------------------------------------------
# 3. BUILD RAW RELATIVE-ABUNDANCE MATRIX (genus-level, all taxa)
# ---------------------------------------------------------------------

taxa_cols <- setdiff(
  colnames(microbiome_filtered),
  c("id", "Run", "Day", "Compartment", "Sample.name", "SampleID")
)

microbiome_matrix <- microbiome_filtered[, taxa_cols]
row_names_temp <- rownames(microbiome_matrix)
microbiome_matrix <- as.data.frame(lapply(microbiome_matrix, as.numeric))
rownames(microbiome_matrix) <- row_names_temp
microbiome_matrix[is.na(microbiome_matrix)] <- 0
rownames(microbiome_matrix) <- microbiome_filtered$SampleID

# Convert read counts to relative abundance (each sample sums to 1)
microbiome_matrix <- microbiome_matrix / rowSums(microbiome_matrix)

# =======================================================================
# PANEL A — Bray-Curtis PCoA, colored by Day, shaped by Run
# (computed on the raw relative-abundance matrix, per Methods 2.5)
# =======================================================================

dist_micro <- vegdist(microbiome_matrix, method = "bray")
pcoa <- cmdscale(dist_micro, k = 2, eig = TRUE)

eig_vals     <- pcoa$eig
eig_vals_pos <- eig_vals[eig_vals > 0]
pc1_var <- round(100 * eig_vals_pos[1] / sum(eig_vals_pos), 1)
pc2_var <- round(100 * eig_vals_pos[2] / sum(eig_vals_pos), 1)

# metadata aligned to microbiome_matrix row order
meta_df <- dplyr::select(microbiome_filtered, SampleID, Run, Day, Compartment)
match_idx <- match(rownames(microbiome_matrix), meta_df$SampleID)
meta_df <- dplyr::slice(meta_df, match_idx)

points_df <- as.data.frame(pcoa$points)
colnames(points_df) <- c("PC1", "PC2")
pcoa_df <- cbind(points_df, meta_df)
pcoa_df$Day_group <- as.factor(pcoa_df$Day)
pcoa_df$Run <- as.factor(pcoa_df$Run)

day_colors <- setNames(
  RColorBrewer::brewer.pal(10, "Spectral"),
  sort(unique(pcoa_df$Day_group))
)

panelA <- ggplot(pcoa_df, aes(PC1, PC2, color = Day_group, shape = Run)) +
  geom_point(size = 3) +
  scale_color_manual(values = day_colors, name = "Day_group") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Run") +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PCoA1 (", pc1_var, "%)"),
    y = paste0("PCoA2 (", pc2_var, "%)")
  )

# =======================================================================
# PANEL B — Temporal stability: mean +/- SD Bray-Curtis, day-to-day,
# within each run, faceted by compartment
# =======================================================================

bray_by_run <- microbiome_filtered %>%
  group_split(Run) %>%
  lapply(function(df) {

    cols <- setdiff(colnames(df),
                    c("id", "Run", "Day", "Compartment", "Sample.name", "SampleID"))

    mat <- as.data.frame(df[, cols])
    mat[] <- lapply(mat, as.numeric)
    mat[is.na(mat)] <- 0
    mat <- mat / rowSums(mat)   # relative abundance

    rownames(mat) <- paste0(df$Day, "_", df$Compartment)

    dist <- vegdist(mat, method = "bray")
    dist_df <- as.data.frame(as.matrix(dist))
    dist_df$Sample1 <- rownames(dist_df)

    dist_long <- tidyr::pivot_longer(
      dist_df, -Sample1,
      names_to = "Sample2", values_to = "Bray"
    )
    dist_long$Run <- unique(df$Run)
    dist_long
  }) %>%
  bind_rows()

bray_consecutive <- bray_by_run %>%
  separate(Sample1, into = c("Day1", "Comp1"), sep = "_") %>%
  separate(Sample2, into = c("Day2", "Comp2"), sep = "_") %>%
  mutate(Day1 = as.numeric(Day1), Day2 = as.numeric(Day2)) %>%
  filter(Comp1 == Comp2, Day2 == Day1 + 1)

bray_summary <- bray_consecutive %>%
  group_by(Comp1, Day1) %>%
  summarise(
    mean_bray = mean(Bray, na.rm = TRUE),
    sd_bray   = sd(Bray, na.rm = TRUE),
    .groups = "drop"
  )

panelB <- ggplot(bray_summary, aes(x = Day1, y = mean_bray)) +
  geom_ribbon(aes(ymin = mean_bray - sd_bray, ymax = mean_bray + sd_bray),
              alpha = 0.2, fill = "steelblue") +
  geom_line(linewidth = 1, color = "steelblue4") +
  geom_point(size = 2, color = "steelblue4") +
  facet_wrap(~Comp1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = length(unique(bray_summary$Day1)))) +
  theme_classic(base_size = 14) +
  labs(
    x = "Day",
    y = "Bray-Curtis distance (Day-to-Day)"
  )

# =======================================================================
# COMBINE — patchwork auto A/B tags, upper-left, bold
# =======================================================================

figure2 <- panelA + panelB +
  plot_layout(ncol = 2, widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

figure2

# =======================================================================
# EXPORT — vector PDF (no DPI concerns; both panels are data plots)
# =======================================================================

ggsave(
  filename = file.path(out_dir, "Figure2.pdf"),
  plot = figure2,
  width = 12,
  height = 5.5,
  units = "in"
)

cat("Saved:", file.path(out_dir, "Figure2.pdf"), "\n")
cat("PC1 variance explained:", pc1_var, "%\n")
cat("PC2 variance explained:", pc2_var, "%\n")


# ---- 2. PERMANOVA (main, per run, windowed) and 3. Figure S3 (Jaccard) ----
library(dplyr)
library(vegan)

# ---------------------------------------------------------------------
# Shared metadata, aligned to microbiome_matrix row order, with Stage
# ---------------------------------------------------------------------

meta_perm <- dplyr::select(microbiome_filtered, SampleID, Run, Day, Compartment)
match_idx <- match(rownames(microbiome_matrix), meta_perm$SampleID)
meta_perm <- dplyr::slice(meta_perm, match_idx)

meta_perm$Stage <- ifelse(meta_perm$Day <= 5, "Early-stage", "Late-stage")
meta_perm$Stage <- factor(meta_perm$Stage, levels = c("Early-stage", "Late-stage"))
meta_perm$Run <- factor(meta_perm$Run)
meta_perm$Compartment <- factor(meta_perm$Compartment)

# dist_micro should already exist from Figure2_complete.R;
# recompute here in case you're running this fresh
dist_micro <- vegdist(microbiome_matrix, method = "bray")

# =======================================================================
# 1. MAIN PERMANOVA — Stage + Run
# =======================================================================

set.seed(123)
perm_stage_run <- adonis2(dist_micro ~ Stage + Run, data = meta_perm,
                          permutations = 999, by = "terms")
cat("\n===== MAIN PERMANOVA: Stage + Run =====\n")
print(perm_stage_run)

# If your original model also included Compartment as a covariate,
# uncomment and use this instead:
# perm_stage_run_full <- adonis2(dist_micro ~ Stage + Compartment + Run,
#                                 data = meta_perm, permutations = 999)
# print(perm_stage_run_full)

# =======================================================================
# 2. PER-RUN PERMANOVA — "temporal effects significant across all runs"
# Tests whether Day is a significant predictor of composition WITHIN
# each run separately
# =======================================================================

cat("\n===== PER-RUN PERMANOVA: Day within each Run =====\n")

run_ids <- unique(meta_perm$Run)
perm_by_run <- list()

for (r in run_ids) {
  ids <- meta_perm$SampleID[meta_perm$Run == r]
  mat_r <- microbiome_matrix[rownames(microbiome_matrix) %in% ids, , drop = FALSE]

  meta_r <- meta_perm[meta_perm$Run == r, ]
  idx_r  <- match(rownames(mat_r), meta_r$SampleID)
  meta_r <- meta_r[idx_r, ]

  dist_r <- vegdist(mat_r, method = "bray")
  res_r  <- adonis2(dist_r ~ Day, data = meta_r, permutations = 999)

  perm_by_run[[as.character(r)]] <- res_r
  cat("\n--- Run:", as.character(r), "---\n")
  print(res_r)
}

# =======================================================================
# 3. WINDOWED PERMANOVA — Day effect within expanding windows anchored
# at the end of the operational period (Day 10), per compartment
# =======================================================================

cat("\n===== WINDOWED PERMANOVA: Day within expanding windows, by Compartment =====\n")

compartments  <- unique(meta_perm$Compartment)
window_starts <- 7:1   # produces windows: Days 7-10, 6-10, 5-10, ..., 1-10

windowed_results <- data.frame()

for (comp in compartments) {
  for (start_day in window_starts) {

    ids <- meta_perm$SampleID[
      meta_perm$Compartment == comp & meta_perm$Day >= start_day
    ]

    # need at least 2 distinct Day values in the window to test an effect
    days_in_window <- unique(meta_perm$Day[meta_perm$SampleID %in% ids])
    if (length(days_in_window) < 2) next

    mat_w <- microbiome_matrix[rownames(microbiome_matrix) %in% ids, , drop = FALSE]
    meta_w <- meta_perm[meta_perm$SampleID %in% ids, ]
    idx_w  <- match(rownames(mat_w), meta_w$SampleID)
    meta_w <- meta_w[idx_w, ]

    dist_w <- vegdist(mat_w, method = "bray")
    res_w  <- adonis2(dist_w ~ Day, data = meta_w, permutations = 999)

    windowed_results <- rbind(windowed_results, data.frame(
      Compartment = as.character(comp),
      Window = paste0("Days ", start_day, "-10"),
      n_samples = nrow(mat_w),
      R2 = round(res_w$R2[1], 3),
      p  = res_w$`Pr(>F)`[1]
    ))
  }
}

cat("\n")
print(windowed_results)

# Quick check against your original numbers:
cat("\nAC window R2 values:\n")
print(windowed_results[windowed_results$Compartment == "AC", ])  

# =====================================================================
# JACCARD PCoA + PERMANOVA — rebuilt from RAW presence/absence
# (previous version thresholded the CLR-transformed matrix at >0,
#  which is NOT valid presence/absence — CLR values are centered on
#  each sample's own geometric mean, not on "taxon detected or not")
#
# Run this AFTER Figure2_complete.R (needs microbiome_matrix,
# microbiome_filtered already built)
# =====================================================================

library(dplyr)
library(vegan)
library(ggplot2)
library(RColorBrewer)

# ---------------------------------------------------------------------
# 1. RAW presence/absence: any taxon with count/abundance > 0 is present
# ---------------------------------------------------------------------

microbiome_pa <- microbiome_matrix
microbiome_pa <- ifelse(microbiome_pa > 0, 1, 0)

# ---------------------------------------------------------------------
# 2. JACCARD DISTANCE + PCoA
# ---------------------------------------------------------------------

dist_jaccard <- vegdist(microbiome_pa, method = "jaccard")
pcoa_j <- cmdscale(dist_jaccard, k = 2, eig = TRUE)

eig_vals_j     <- pcoa_j$eig
eig_vals_j_pos <- eig_vals_j[eig_vals_j > 0]
pc1_var_j <- round(100 * eig_vals_j_pos[1] / sum(eig_vals_j_pos), 1)
pc2_var_j <- round(100 * eig_vals_j_pos[2] / sum(eig_vals_j_pos), 1)

# metadata aligned to microbiome_pa row order
meta_j <- dplyr::select(microbiome_filtered, SampleID, Run, Day, Compartment)
match_idx_j <- match(rownames(microbiome_pa), meta_j$SampleID)
meta_j <- dplyr::slice(meta_j, match_idx_j)

points_df_j <- as.data.frame(pcoa_j$points)
colnames(points_df_j) <- c("PC1", "PC2")
pcoa_j_df <- cbind(points_df_j, meta_j)

pcoa_j_df$Day_group <- as.factor(pcoa_j_df$Day)
pcoa_j_df$Run <- as.factor(pcoa_j_df$Run)

day_colors_j <- setNames(
  RColorBrewer::brewer.pal(10, "Spectral"),
  sort(unique(pcoa_j_df$Day_group))
)

panelJaccard <- ggplot(pcoa_j_df, aes(PC1, PC2, color = Day_group, shape = Run)) +
  geom_point(size = 3) +
  scale_color_manual(values = day_colors_j, name = "Day_group") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Run") +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PCoA1 (", pc1_var_j, "%)"),
    y = paste0("PCoA2 (", pc2_var_j, "%)")
  )

panelJaccard

ggsave(
  filename = file.path(out_dir, "FigureS3_Jaccard.pdf"),
  plot = panelJaccard,
  width = 7,
  height = 5.5,
  units = "in"
)

cat("Jaccard PC1 variance explained:", pc1_var_j, "%\n")
cat("Jaccard PC2 variance explained:", pc2_var_j, "%\n")

# ---------------------------------------------------------------------
# 3. JACCARD PERMANOVA (Stage + Run, matching your main Bray-Curtis test)
# ---------------------------------------------------------------------

meta_j$Stage <- ifelse(meta_j$Day <= 5, "Early-stage", "Late-stage")
meta_j$Stage <- factor(meta_j$Stage, levels = c("Early-stage", "Late-stage"))
meta_j$Run <- factor(meta_j$Run)
meta_j$Compartment <- factor(meta_j$Compartment)

set.seed(123)
perm_jaccard <- adonis2(dist_jaccard ~ Stage + Run, data = meta_j,
                        permutations = 999, by = "terms")
cat("\n===== JACCARD PERMANOVA: Stage + Run =====\n")
print(perm_jaccard)

# =======================================================================
# PANEL B (Jaccard) — day-to-day Jaccard dissimilarity, mean +/- SD,
# faceted by compartment. Mirrors the Bray-Curtis Panel B structure.
# =======================================================================

library(tidyr)
library(patchwork)
library(scales)

jaccard_by_run <- microbiome_filtered %>%
  group_split(Run) %>%
  lapply(function(df) {

    cols <- setdiff(colnames(df),
                    c("id", "Run", "Day", "Compartment", "Sample.name", "SampleID"))

    mat <- as.data.frame(df[, cols])
    mat[] <- lapply(mat, as.numeric)
    mat[is.na(mat)] <- 0
    mat <- (mat > 0) * 1   # raw presence/absence, same rule as Panel A

    rownames(mat) <- paste0(df$Day, "_", df$Compartment)

    dist <- vegdist(mat, method = "jaccard")
    dist_df <- as.data.frame(as.matrix(dist))
    dist_df$Sample1 <- rownames(dist_df)

    dist_long <- tidyr::pivot_longer(
      dist_df, -Sample1,
      names_to = "Sample2", values_to = "Jaccard"
    )
    dist_long$Run <- unique(df$Run)
    dist_long
  }) %>%
  bind_rows()

jaccard_consecutive <- jaccard_by_run %>%
  separate(Sample1, into = c("Day1", "Comp1"), sep = "_") %>%
  separate(Sample2, into = c("Day2", "Comp2"), sep = "_") %>%
  mutate(Day1 = as.numeric(Day1), Day2 = as.numeric(Day2)) %>%
  filter(Comp1 == Comp2, Day2 == Day1 + 1)

jaccard_summary <- jaccard_consecutive %>%
  group_by(Comp1, Day1) %>%
  summarise(
    mean_jaccard = mean(Jaccard, na.rm = TRUE),
    sd_jaccard   = sd(Jaccard, na.rm = TRUE),
    .groups = "drop"
  )

panelJaccardB <- ggplot(jaccard_summary, aes(x = Day1, y = mean_jaccard)) +
  geom_ribbon(aes(ymin = mean_jaccard - sd_jaccard, ymax = mean_jaccard + sd_jaccard),
              alpha = 0.2, fill = "steelblue") +
  geom_line(linewidth = 1, color = "steelblue4") +
  geom_point(size = 2, color = "steelblue4") +
  facet_wrap(~Comp1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = length(unique(jaccard_summary$Day1)))) +
  theme_classic(base_size = 14) +
  labs(
    x = "Day",
    y = "Jaccard distance (Day-to-Day)"
  )

panelJaccardB

# ---------------------------------------------------------------------
# COMBINE Panel A + Panel B into FigureS2 — same style as Figure 2
# ---------------------------------------------------------------------

figureS2 <- panelJaccard + panelJaccardB +
  plot_layout(ncol = 2, widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

figureS2

ggsave(
  filename = file.path(out_dir, "FigureS3_Jaccard.pdf"),
  plot = figureS2,
  width = 12,
  height = 5.5,
  units = "in"
)

cat("\nJaccard day-to-day stabilization value (mean +/- SD across all compartments, after Day 5):\n")
jaccard_late <- jaccard_summary %>% filter(Day1 >= 5)
cat(round(mean(jaccard_late$mean_jaccard), 3), "+/-", round(mean(jaccard_late$sd_jaccard), 3), "\n")


# ---- 4. Stabilised day-to-day Bray–Curtis after day 5 (Results 3.2) ----

library(vegan)
library(dplyr)

# Convert counts to relative abundance (Bray-Curtis should use
# proportions, not raw counts or CLR-transformed data)
rel_abund <- microbiome_filtered %>%
  dplyr::select(all_of(taxa_cols)) %>%
  as.matrix()
rel_abund[is.na(rel_abund)] <- 0
rel_abund <- rel_abund / rowSums(rel_abund)

meta_seq <- microbiome_filtered %>%
  dplyr::select(SampleID, Run, Day, Compartment)

# Function: for one Run x Compartment, compute Bray-Curtis
# between consecutive days (Day t vs Day t+1)
compute_sequential_bc <- function(run_name, comp_name) {
  idx <- meta_seq$Run == run_name & meta_seq$Compartment == comp_name
  meta_sub <- meta_seq[idx, ]
  mat_sub  <- rel_abund[idx, , drop = FALSE]

  # Order by Day
  ord <- order(meta_sub$Day)
  meta_sub <- meta_sub[ord, ]
  mat_sub  <- mat_sub[ord, , drop = FALSE]

  days <- meta_sub$Day
  n <- nrow(mat_sub)
  if (n < 2) return(NULL)

  out <- data.frame(
    Run = run_name,
    Compartment = comp_name,
    Day = days[-1],           # dissimilarity assigned to the later day
    BrayCurtis = NA_real_
  )

  for (i in 2:n) {
    # Only compute for truly consecutive days (skip gaps)
    if (days[i] - days[i-1] == 1) {
      d <- vegdist(rbind(mat_sub[i-1, ], mat_sub[i, ]), method = "bray")
      out$BrayCurtis[i-1] <- as.numeric(d)
    }
  }

  out
}

# Run across all Run x Compartment combinations
combos <- expand.grid(
  Run = unique(meta_seq$Run),
  Compartment = unique(meta_seq$Compartment),
  stringsAsFactors = FALSE
)

sequential_bc <- bind_rows(
  lapply(1:nrow(combos), function(i) {
    compute_sequential_bc(combos$Run[i], combos$Compartment[i])
  })
) %>%
  filter(!is.na(BrayCurtis))

# --- Sanity check: reproduce the shape of Fig 2b ---
sequential_bc %>%
  group_by(Compartment, Day) %>%
  summarise(mean_BC = mean(BrayCurtis, na.rm = TRUE),
            sd_BC   = sd(BrayCurtis, na.rm = TRUE), .groups = "drop") %>%
  arrange(Compartment, Day) %>%
  print(n = Inf)

# --- The number you need: stabilized average after Day 5 ---
stabilized_avg <- sequential_bc %>%
  filter(Day >= 6) %>%
  summarise(
    mean_BC = round(mean(BrayCurtis, na.rm = TRUE), 3),
    sd_BC   = round(sd(BrayCurtis, na.rm = TRUE), 3)
  )

print(stabilized_avg)

# Per-compartment breakdown, in case you want to report it that way instead
stabilized_avg_by_compartment <- sequential_bc %>%
  filter(Day >= 6) %>%
  group_by(Compartment) %>%
  summarise(
    mean_BC = round(mean(BrayCurtis, na.rm = TRUE), 3),
    sd_BC   = round(sd(BrayCurtis, na.rm = TRUE), 3)
  )

print(stabilized_avg_by_compartment)

