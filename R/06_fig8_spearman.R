# =====================================================================
# 06  Figure 8 — Spearman correlations, DIABLO-selected genera x metabolites
# Chapter 3 — "Stability as a dynamic process" (J Appl Microbiol submission)
# Inputs : objects from 05 (diablo_res, microbiome_final, metab_final) — run 05 first
# Outputs: figures/Figure8_Spearman_heatmap.*; table of FDR-significant pairs
# NOTE   : correlations use the CLR-transformed genus abundances from 05 — say so
#          in Methods. Asterisks are BH-adjusted (FDR) — legend should say FDR.
# CHANGED (presentation only): metabolite labels (Heptanoic acid, sentence case,
#          subscripts); Figure 8 file names.
# Run from the repository root:  source("R/<this file>")
# =====================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# PDF device that renders subscripts (NH3, H2S): the built-in Quartz PDF
# device on a Mac (no XQuartz needed); cairo_pdf elsewhere; otherwise pdf()
pdf_device <- function(filename, width, height, ...) {
  if (isTRUE(capabilities("aqua"))) {
    grDevices::quartz(type = "pdf", file = filename, width = width, height = height, ...)
  } else if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(filename, width = width, height = height, ...)
  } else {
    grDevices::pdf(filename, width = width, height = height, ...)
  }
}

############################
# STEP 1 — EXTRACT DIABLO-SELECTED GENERA
############################

micro_loadings <- diablo_res$loadings$microbiome[, 1]
selected_genera <- names(micro_loadings[micro_loadings != 0])

cat("DIABLO-selected genera:\n")
print(selected_genera)
cat("N =", length(selected_genera), "\n")

microbiome_spearman <- dplyr::select(microbiome_final, dplyr::all_of(selected_genera))

############################
# STEP 2 — CLEAN GENUS NAMES FOR DISPLAY
############################

clean_names <- colnames(microbiome_spearman)
clean_names <- gsub("\\.", " ", clean_names)
clean_names <- gsub("\\s+", " ", clean_names)
clean_names <- gsub(" unc(\\s|$)", " (unc.)", clean_names)
clean_names <- gsub(" IS(\\s|$)", " (Incertae Sedis)", clean_names)
clean_names <- trimws(clean_names)

# Strip the disambiguating suffix (make.names artefact) from any label
# that appears only once among the genera actually shown in this figure —
# no real ambiguity to preserve, so the arbitrary number is just confusing.
# (Same fix applied to Figure 7's loading plots.)
strip_singleton_suffix <- function(labels) {
  base_labels <- sub("\\(unc\\.\\)\\d+$", "unc.", labels)
  base_labels <- sub("\\(Incertae Sedis\\)\\d+$", "Incertae Sedis", base_labels)
  base_labels <- ifelse(base_labels == labels,
                        sub("\\d+$", "", labels),
                        base_labels)
  base_labels <- trimws(base_labels)

  counts <- table(base_labels)
  ifelse(counts[base_labels] > 1, labels, base_labels)
}

clean_names <- strip_singleton_suffix(clean_names)
colnames(microbiome_spearman) <- clean_names

############################
# STEP 3 — METABOLITE MATRIX
############################

metab_spearman <- metab_final %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::filter(SampleID %in% rownames(microbiome_spearman)) %>%
  dplyr::arrange(match(SampleID, rownames(microbiome_spearman))) %>%
  tibble::column_to_rownames("SampleID")

met_display <- colnames(metab_spearman)
met_display <- gsub("\\.", " ", met_display)
met_display <- gsub("\\s+", " ", met_display)
met_display <- trimws(met_display)
colnames(metab_spearman) <- met_display

############################
# STEP 4 — SPEARMAN CORRELATIONS WITH FDR
############################

results <- list()

for (tax in colnames(microbiome_spearman)) {
  for (met in colnames(metab_spearman)) {
    ct <- cor.test(
      microbiome_spearman[[tax]],
      metab_spearman[[met]],
      method = "spearman",
      exact  = FALSE
    )
    results[[length(results) + 1]] <- data.frame(
      Taxa       = tax,
      Metabolite = met,
      rho        = as.numeric(ct$estimate),
      p          = ct$p.value,
      stringsAsFactors = FALSE
    )
  }
}

res_df <- dplyr::bind_rows(results)
res_df$p_adj <- p.adjust(res_df$p, method = "BH")
res_df$star  <- cut(
  res_df$p_adj,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)
res_df$star_col <- ifelse(res_df$rho > 0, "black", "black")

cat("\nSignificant correlations (FDR < 0.05):", sum(res_df$p_adj < 0.05), "\n")
cat("Total tested:", nrow(res_df), "\n")

############################
# STEP 5 — ORDER ROWS BY DIABLO LOADING
# Early-stage (negative) at top, late-stage (positive) at bottom
############################

loading_order <- sort(micro_loadings[micro_loadings != 0])
name_map      <- setNames(clean_names, selected_genera)
clean_ordered <- name_map[names(loading_order)]

res_df$Taxa <- factor(res_df$Taxa, levels = clean_ordered)

############################
# STEP 6 — ORDER METABOLITE COLUMNS
# SCFA = straight-chain (acetic, propionic, butyric, valeric, caproic, heptanoic)
# BCFA = branched-chain only (isobutyric, isovaleric, isocaproic)
############################

metab_order <- c(
  # SCFA (straight-chain)
  "Acetic acid", "Propionic Acid", "Butyric Acid",
  "Valeric Acid", "Caproic Acid", "Heptapentoic Acid",   # column name in data; displayed as Heptanoic acid
  # BCFA (branched-chain only)
  "IsoButyric Acid", "Isovaleric Acid", "Isocaproic Acid",
  # Fermentation environment markers
  "NH3", "H2S", "FRAP"
)

# NOTE: check your actual column name — confirm "Heptapentoic Acid" matches
# your real column exactly, or it will silently fall into "remaining" below.

metab_order <- metab_order[metab_order %in% unique(res_df$Metabolite)]
remaining   <- setdiff(unique(res_df$Metabolite), metab_order)
metab_order <- c(metab_order, remaining)

res_df$Metabolite <- factor(res_df$Metabolite, levels = metab_order)

# Display labels (presentation only): sentence case, heptanoic acid, subscripts
metab_display_labels <- c(
  "Acetic acid" = "Acetic acid", "Propionic Acid" = "Propionic acid",
  "Butyric Acid" = "Butyric acid", "Valeric Acid" = "Valeric acid",
  "Caproic Acid" = "Caproic acid", "Heptapentoic Acid" = "Heptanoic acid",
  "IsoButyric Acid" = "Isobutyric acid", "Isovaleric Acid" = "Isovaleric acid",
  "Isocaproic Acid" = "Isocaproic acid",
  "NH3" = "NH<sub>3</sub>", "H2S" = "H<sub>2</sub>S", "FRAP" = "FRAP")

############################
# STEP 6B — EARLY/LATE STAGE SIDEBAR (y-axis annotation)
############################

stage_df <- data.frame(
  Taxa  = clean_ordered,
  Stage = ifelse(loading_order < 0, "Early-stage", "Late-stage")
)
stage_df$Taxa <- factor(stage_df$Taxa, levels = clean_ordered)

stage_colours <- c("Early-stage" = "#E69F00", "Late-stage" = "#009E73")

############################
# STEP 6C — METABOLITE GROUP LABELS FOR FACET HEADER
############################

metabolite_group_map <- c(
  "Acetic acid"        = "SCFA",
  "Propionic Acid"     = "SCFA",
  "Butyric Acid"       = "SCFA",
  "Valeric Acid"       = "SCFA",
  "Caproic Acid"       = "SCFA",
  "Heptapentoic Acid"  = "SCFA",
  "IsoButyric Acid"    = "BCFA",
  "Isovaleric Acid"    = "BCFA",
  "Isocaproic Acid"    = "BCFA",
  "NH3"                = "Other metabolites",
  "H2S"                = "Other metabolites",
  "FRAP"               = "Other metabolites"
)

res_df$MetaboliteGroup <- metabolite_group_map[as.character(res_df$Metabolite)]
res_df$MetaboliteGroup <- factor(res_df$MetaboliteGroup,
                                 levels = c("SCFA", "BCFA", "Other metabolites"))

############################
# STEP 7 — HEATMAP WITH FACETED COLUMN GROUPS
############################

p_heat <- ggplot(res_df, aes(x = Metabolite, y = Taxa, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = star, colour = star_col),
            size = 3.5, vjust = 0.75) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#D6604D",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Spearman ρ",
    guide    = guide_colourbar(display = "rectangles")   # raster legends can render flipped on the Quartz PDF device
  ) +
  scale_colour_identity() +
  facet_grid(~ MetaboliteGroup, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(metab_display_labels), metab_display_labels[x], x)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x       = ggtext::element_markdown(angle = 45, hjust = 1, size = 10),  # renders NH3/H2S subscripts on any device
    axis.text.y       = element_text(size = 9, face = "italic"),
    legend.position   = "right",
    legend.key.height = unit(1.2, "cm"),
    panel.border      = element_rect(colour = "grey70", fill = NA),
    strip.background  = element_rect(fill = "grey90", colour = "grey70"),
    strip.text        = element_text(size = 10, face = "bold")
  )

############################
# STEP 6D — SIDEBAR, HEIGHT-MATCHED TO HEATMAP
############################

p_stage <- ggplot(stage_df, aes(x = 1, y = Taxa, fill = Stage)) +
  geom_tile(width = 1, height = 0.9) +
  scale_fill_manual(values = stage_colours, name = "DIABLO stage") +
  facet_grid(~ "") +   # blank facet, matches strip row height of p_heat
  theme_void() +
  theme(
    legend.position  = "left",
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10),
    strip.text       = element_text(size = 10, colour = "white"),  # invisible but occupies space
    strip.background = element_rect(fill = "white", colour = "white")
  )

############################
# STEP 7B — COMBINE SIDEBAR + HEATMAP
############################

p_combined <- p_stage + p_heat +
  plot_layout(widths = c(0.06, 1), guides = "collect")   # both legends on the right, none clipped

p_combined

############################
# STEP 8 — SAVE
############################

out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tiff(
  file.path(out_dir, "Figure8_Spearman_heatmap.tiff"),
  width = 11.5, height = 7.3,
  units = "in", res = 600, compression = "lzw"
)
print(p_combined)
dev.off()

pdf_device(file.path(out_dir, "Figure8_Spearman_heatmap.pdf"), width = 11.5, height = 7.3)
print(p_combined)
dev.off()

cat("Saved:", file.path(out_dir, "Figure8_Spearman_heatmap.pdf"), "\n")

############################
# STEP 9 — PRINT SIGNIFICANT CORRELATIONS
############################

cat("\n=== Significant correlations (FDR < 0.05) ===\n")
sig <- res_df[res_df$p_adj < 0.05, ]
print(sig)  
