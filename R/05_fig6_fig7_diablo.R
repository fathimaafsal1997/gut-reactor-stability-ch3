# =====================================================================
# 05  Figures 6–7 — DIABLO (microbiome + metabolites), CV error, permutation test
# Chapter 3 — "Stability as a dynamic process" (J Appl Microbiol submission)
# Inputs : data/microbiome/TR*/..., data/Metabolites.csv
# Outputs: figures/Figure6_DIABLO_score.*, figures/Figure7_DIABLO_loadings.*;
#          BER (Results 3.4: 0.281), permutation p
# CHECK  : keepX is fixed (20 genera, 10 metabolites); Methods say it was
#          'selected based on cross-validation performance' — no tune.block.splsda here.
# CHECK  : permutation test now 999 permutations (the earlier run used 99, which
#          cannot give p < 0.001). Takes several hours.
# NOTE   : Leave this session open and run 06 next (it uses diablo_res).
# Run from the repository root:  source("R/<this file>")
# =====================================================================

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(mixOmics)
library(ggplot2)
library(patchwork)
library(igraph)

#####################
# STEP 1 — LOAD MICROBIOME DATA
#####################

base_dir <- "data/microbiome"
runs <- c("TR1","TR2","TR3","TR4")

microbiome_list <- list()

for(run in runs){
  taxa_file <- file.path(base_dir, run, paste0("alldays", run, ".csv"))
  meta_file <- file.path(base_dir, run, paste0("alldays", run, "m.csv"))

  taxa <- read.csv(taxa_file, stringsAsFactors = FALSE)
  meta <- read_csv(meta_file, show_col_types = FALSE)

  merged <- meta %>%
    left_join(taxa, by = "id") %>%
    separate(Sample.name,
             into = c("Run2","Day","Compartment"),
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

#####################
# STEP 2 — FILTER DAYS 1-10
#####################
microbiome_filtered <- microbiome_all %>%
  filter(Day %in% 1:10)

#####################
# STEP 3 — MICROBIOME MATRIX
# (no comma-bracket indexing — dplyr::select + lapply pattern only)
#####################

taxa_cols <- setdiff(colnames(microbiome_filtered),
                     c("id","Run","Day","Compartment","Sample.name","SampleID"))

microbiome_matrix <- dplyr::select(microbiome_filtered, dplyr::all_of(taxa_cols))
microbiome_matrix <- as.data.frame(lapply(microbiome_matrix, as.numeric))
rownames(microbiome_matrix) <- microbiome_filtered$SampleID

#####################
# STEP 3.5 — CLEAN GENUS NAMES
# For uncultured genus: report family level instead
# For uncultured family too: report order
#####################

clean_genus_smart <- function(x) {
  result <- character(length(x))

  for(i in seq_along(x)) {
    tag <- x[i]

    get_level <- function(prefix, filter_bad = TRUE) {
      pattern <- paste0(prefix, "([^.]*(?:\\.[^_][^.]*)*?)(?=\\.[a-z]__|$)")
      m <- regmatches(tag, regexpr(pattern, tag, perl = TRUE))
      if(length(m) == 0 || m == "") return("")
      val <- sub(paste0("^", prefix), "", m)
      val <- trimws(val)
      if(filter_bad && (tolower(val) == "uncultured" ||
                        grepl("incertae_sedis", tolower(val)) ||
                        val == "")) return("")
      return(val)
    }

    # Filtered levels — empty if uncultured or Incertae Sedis
    g  <- get_level("g__")
    f  <- get_level("f__")
    o  <- get_level("o__")
    cl <- get_level("c__")
    p  <- get_level("p__")

    # Raw genus — unfiltered, to distinguish IS vs uncultured for labelling
    g_raw <- get_level("g__", filter_bad = FALSE)

    if(g != "") {
      result[i] <- g
    } else if(f != "") {
      if(grepl("incertae_sedis", tolower(g_raw))) {
        result[i] <- paste0(f, " IS")
      } else {
        result[i] <- paste0(f, " unc")
      }
    } else if(o != "") {
      if(grepl("incertae_sedis", tolower(g_raw))) {
        result[i] <- paste0(o, " IS")
      } else {
        result[i] <- paste0(o, " unc")
      }
    } else if(cl != "") {
      result[i] <- paste0(cl, " unc")
    } else if(p != "") {
      result[i] <- paste0(p, " unc")
    } else {
      result[i] <- paste0("Unclassified_", i)
    }
  }

  make.names(result, unique = TRUE)
}

colnames(microbiome_matrix) <- clean_genus_smart(colnames(microbiome_matrix))

#####################
# STEP 4 — LOAD METABOLITES
#####################

metab <- read.csv(
  "data/Metabolites.csv",
  stringsAsFactors = FALSE
)
metab$SampleID <- paste(metab$Run, metab$Day, metab$Compartment, sep = "_")

#####################
# STEP 5 — ALIGN SAMPLES
#####################

metab <- metab %>%
  filter(SampleID %in% microbiome_filtered$SampleID) %>%
  arrange(match(SampleID, microbiome_filtered$SampleID))

#####################
# STEP 6 — METABOLITE MATRIX
# (no comma-bracket indexing — dplyr::select(-any_of(...)) instead)
#####################

metab_matrix <- metab %>%
  dplyr::select(-Run, -Day, -Compartment, -SampleID)

metab_matrix <- as.data.frame(lapply(metab_matrix, as.numeric))
rownames(metab_matrix) <- metab$SampleID

# Remove derived columns
drop_cols <- c("Total", "total.SCFA", "total.BCFA", "Others")
metab_matrix <- dplyr::select(metab_matrix, -dplyr::any_of(drop_cols))

# Standardise metabolite names — use "Acetic.acid" throughout
colnames(metab_matrix) <- gsub("Acetic.Acid|acetic.acid|Acetate|acetate",
                               "Acetic.acid",
                               colnames(metab_matrix))

#####################
# STEP 7 — CLR TRANSFORM
#####################

microbiome_counts <- as.matrix(microbiome_matrix)
microbiome_counts[is.na(microbiome_counts)] <- 0
microbiome_counts <- microbiome_counts + 1  # pseudocount

microbiome_log <- log(microbiome_counts)
microbiome_clr <- microbiome_log - rowMeans(microbiome_log)
microbiome_rel <- as.data.frame(microbiome_clr)
rownames(microbiome_rel) <- rownames(microbiome_matrix)

#####################
# STEP 8 — FILTER TAXA (top 50% variance)
#####################

taxa_var <- apply(microbiome_rel, 2, var)
keep_taxa <- taxa_var > quantile(taxa_var, 0.5)
microbiome_filtered_matrix <- dplyr::select(microbiome_rel, dplyr::all_of(names(taxa_var)[keep_taxa]))

#####################
# STEP 9 — ALIGN FINAL MATRICES
# (no comma-bracket indexing — rownames_to_column/arrange/column_to_rownames)
#####################

common_samples <- intersect(
  rownames(microbiome_filtered_matrix),
  rownames(metab_matrix)
)

microbiome_final <- microbiome_filtered_matrix %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::filter(SampleID %in% common_samples) %>%
  dplyr::arrange(SampleID) %>%
  tibble::column_to_rownames("SampleID")

metab_final <- metab_matrix %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::filter(SampleID %in% common_samples) %>%
  dplyr::arrange(SampleID) %>%
  tibble::column_to_rownames("SampleID")

stopifnot(identical(rownames(microbiome_final), rownames(metab_final)))

#####################
# STEP 10 — METADATA
# (arranged the same way as microbiome_final/metab_final, by SampleID)
#####################

meta_df <- microbiome_filtered %>%
  dplyr::filter(SampleID %in% common_samples) %>%
  dplyr::select(SampleID, Run, Day, Compartment) %>%
  dplyr::arrange(SampleID)

stopifnot(identical(meta_df$SampleID, rownames(microbiome_final)))

#####################
# STEP 11 — Y VARIABLE
# Early-stage = days 1-5, late-stage = days 6-10
#####################

day_info <- as.numeric(gsub(".*_(\\d+)_.*", "\\1", common_samples))

Y <- ifelse(day_info <= 5, "Early-stage", "Late-stage")
Y <- factor(Y, levels = c("Early-stage", "Late-stage"))

#####################
# STEP 12 — DIABLO INPUT
#####################

X <- list(
  microbiome = microbiome_final,
  metabolite = metab_final
)

design <- matrix(c(
  0, 0.5,
  0.5, 0
), nrow = 2, byrow = TRUE)
rownames(design) <- colnames(design) <- c("microbiome","metabolite")

#####################
# STEP 13 — TUNE NUMBER OF COMPONENTS (justifies ncomp)
# Run perf() on a 3-component model to check if comp 3 adds value
#####################

set.seed(1)
diablo_tune <- block.splsda(
  X = X,
  Y = Y,
  ncomp = 3,
  keepX = list(microbiome = 30, metabolite = 10),
  design = design
)

set.seed(1)
perf_res <- perf(
  diablo_tune,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)

# Inspect: if comp 3 error rate < comp 2, justify 3 components
print(perf_res$MajorityVote.error.rate)

# Based on perf output, set ncomp (2 or 3) — adjust after inspecting above
# Default: use ncomp = 2 unless perf shows clear improvement at comp 3
ncomp_final <- 2

#####################
# STEP 14 — FINAL DIABLO MODEL
# (sign-flip wrapped in a function so the bracket assignments run as
#  one block rather than as bare top-level console statements)
#####################

set.seed(1)
diablo_res <- block.splsda(
  X = X,
  Y = Y,
  ncomp = ncomp_final,
  keepX = list(microbiome = 20, metabolite = 10),  # top 20 genera
  design = design
)

orient_component <- function(diablo_obj, comp = 1) {
  diablo_obj$variates$microbiome[, comp]  <- -diablo_obj$variates$microbiome[, comp]
  diablo_obj$variates$metabolite[, comp]  <- -diablo_obj$variates$metabolite[, comp]
  diablo_obj$loadings$microbiome[, comp]  <- -diablo_obj$loadings$microbiome[, comp]
  diablo_obj$loadings$metabolite[, comp]  <- -diablo_obj$loadings$metabolite[, comp]
  diablo_obj
}

# Orient component so Late-stage (stable) is positive on Comp 1
late_mean  <- mean(diablo_res$variates$microbiome[Y == "Late-stage",  1])
early_mean <- mean(diablo_res$variates$microbiome[Y == "Early-stage", 1])

if(late_mean < early_mean){
  diablo_res <- orient_component(diablo_res, comp = 1)
}

# Orient Component 2 the same way as Component 1 — Late-stage positive
late_mean2  <- mean(diablo_res$variates$microbiome[Y == "Late-stage",  2])
early_mean2 <- mean(diablo_res$variates$microbiome[Y == "Early-stage", 2])

if (late_mean2 < early_mean2) {
  diablo_res <- orient_component(diablo_res, comp = 2)
}

#####################
# STEP 15 — SCORE PLOT (single panel, Fig 2 style)
# Shape = Run (TR1-TR4), Colour = stability state
# Same shape mapping as Fig 2 PCoA
# Fix: order points so early-stage plotted last (on top) — no orange hidden behind turquoise
#####################

scores_micro <- as.data.frame(diablo_res$variates$microbiome) %>%
  dplyr::select(1:2)
colnames(scores_micro) <- c("Comp1","Comp2")

# % explained variance for axis labels
expl_var <- diablo_res$prop_expl_var$microbiome * 100

plot_df <- data.frame(
  Comp1      = scores_micro$Comp1,
  Comp2      = scores_micro$Comp2,
  Stability  = Y,
  Run        = meta_df$Run,
  Day        = day_info,
  stringsAsFactors = FALSE
)

# Sort so Early-stage plotted last (on top) — fixes hidden-point issue
plot_df <- dplyr::arrange(plot_df, dplyr::desc(Stability))

# Colour palette — consistent with thesis (adjust hex to match Fig 2 exactly)
stability_colours <- c("Early-stage" = "#E69F00", "Late-stage" = "#009E73")

# Shape mapping — same as Fig 2 PCoA
run_shapes <- c("TR1" = 16, "TR2" = 17, "TR3" = 15, "TR4" = 18)

# Figure 6 is a single-panel figure — no tag/letter needed
p_score <- ggplot(plot_df, aes(x = Comp1, y = Comp2,
                               colour = Stability,
                               shape  = Run)) +
  geom_point(size = 3.5, alpha = 0.9, stroke = 0.3) +
  scale_colour_manual(values = stability_colours,
                      name   = "Stability state") +
  scale_shape_manual(values = run_shapes,
                     name   = "Experimental run") +
  labs(
    x = paste0("Component 1 (", round(expl_var[1], 1), "% expl. var.)"),
    y = paste0("Component 2 (", round(expl_var[2], 1), "% expl. var.)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.5, "cm")
  )

p_score


cat("late_mean2:", late_mean2, "\n")
cat("early_mean2:", early_mean2, "\n")
cat("Comp2 range:", range(diablo_res$variates$microbiome[, 2]), "\n")

#####################
# STEP 16 — LOADING PLOTS (top 20 genera, homogenised scale)
#####################

# --- Extract loadings manually for full control ---

get_loadings_df <- function(diablo_obj, block, comp = 1, n_top = 20) {
  ld <- diablo_obj$loadings[[block]][, comp]
  ld <- ld[ld != 0]
  ld <- sort(ld)

  if(length(ld) > n_top) {
    # keep top n by absolute value
    ld <- ld[order(abs(ld), decreasing = TRUE)[1:n_top]]
    ld <- sort(ld)
  }

  data.frame(
    Feature    = names(ld),
    Loading    = as.numeric(ld),
    Direction  = ifelse(as.numeric(ld) > 0, "Late-stage", "Early-stage"),
    stringsAsFactors = FALSE
  )
}

micro_ld  <- get_loadings_df(diablo_res, "microbiome", comp = 1, n_top = 20)
metab_ld  <- get_loadings_df(diablo_res, "metabolite", comp = 1, n_top = 10)

# Clean up metabolite names for display
metab_ld$Feature <- gsub("\\.", " ", metab_ld$Feature)

# Clean up genus names for display (remove make.names artefacts)
micro_ld$Feature <- gsub("\\.", " ", micro_ld$Feature)
micro_ld$Feature <- gsub("\\s+", " ", micro_ld$Feature)
micro_ld$Feature <- gsub(" unc(\\s|$)", " (unc.)", micro_ld$Feature)
micro_ld$Feature <- gsub(" IS(\\s|$)", " (Incertae Sedis)", micro_ld$Feature)
micro_ld$Feature <- trimws(micro_ld$Feature)
micro_ld$Feature <- factor(micro_ld$Feature, levels = micro_ld$Feature)

# Determine shared x-axis limits
x_max <- max(abs(c(micro_ld$Loading, metab_ld$Loading))) * 1.1
x_lim <- c(-x_max, x_max)

# Direction colours consistent with score plot
dir_colours <- c("Early-stage" = "#E69F00", "Late-stage" = "#009E73")

# --- Microbiome loading plot ---
metab_ld$Feature <- factor(metab_ld$Feature, levels = metab_ld$Feature)

# Figure 7 Panel A
p_micro_load <- ggplot(micro_ld, aes(x = Loading,
                                     y = Feature,
                                     fill = Direction)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = dir_colours, name = "Associated state") +
  scale_x_continuous(limits = x_lim,
                     name   = "Component 1 loading (contribution score)") +
  labs(y = "Genus (or family if uncultured)",
       tag = "A") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y     = element_text(size = 9, face = "italic"),
    legend.position = "none",
    plot.tag        = element_text(face = "bold", size = 14)
  )

# --- Metabolite loading plot ---

# Figure 7 Panel B
p_metab_load <- ggplot(metab_ld, aes(x = Loading,
                                     y = Feature,
                                     fill = Direction)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = dir_colours, name = "Associated state") +
  scale_x_continuous(limits = x_lim,
                     name   = "Component 1 loading (contribution score)") +
  labs(y = "Metabolite (mM)",
       tag = "B") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y     = element_text(size = 10),
    legend.position = "right",
    plot.tag        = element_text(face = "bold", size = 14)
  )

#####################
# STEP 17 — COMBINED FIGURE (score + loadings) — exploratory view only, not exported
#####################

fig_combined <- (p_micro_load | p_metab_load) +
  plot_layout(heights = c(1, 1.4))

fig_combined

#####################
# STEP 18 — EXPORT FIGURES
#####################

out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Figure 6: Score plot only (single panel as Dr. Frigon requested) — no tag letter
tiff(
  file.path(out_dir, "Figure6_DIABLO_score.tiff"),
  width = 6, height = 5, units = "in", res = 600, compression = "lzw"
)
print(p_score)
dev.off()

pdf(file.path(out_dir, "Figure6_DIABLO_score.pdf"), width = 6, height = 5)
print(p_score)
dev.off()

# Figure 7: Loading plots (microbiome + metabolite, homogenised scale) — tagged A / B
tiff(
  file.path(out_dir, "Figure7_DIABLO_loadings.tiff"),
  width = 12, height = 7, units = "in", res = 600, compression = "lzw"
)
print(p_micro_load | p_metab_load)
dev.off()

pdf(file.path(out_dir, "Figure7_DIABLO_loadings.pdf"), width = 12, height = 7)
print(p_micro_load | p_metab_load)
dev.off()

cat("Saved:", file.path(out_dir, "Figure6_DIABLO_score.pdf"), "\n")
cat("Saved:", file.path(out_dir, "Figure7_DIABLO_loadings.pdf"), "\n")

# ---- Cross-validation error and permutation test ----

# --- Cross-validation performance ---
set.seed(1)
perf_diablo <- perf(
  diablo_res,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 10,
  progressBar = TRUE
)

# Extract and print error rates
cat("\n=== Cross-validation error rates ===\n")
print(perf_diablo$MajorityVote.error.rate)

# Overall error rate on comp 1 and 2
# Fix CV error extraction
cv_error <- perf_diablo$MajorityVote.error.rate
cat("\nOverall BER (mahalanobis, comp1):", 
    round(cv_error$mahalanobis.dist["Overall.BER", "comp1"], 3), "\n")

real_ber <- cv_error$mahalanobis.dist["Overall.BER", "comp1"]

# --- Permutation test (run n=99 first for quick check) ---
set.seed(1)
n_perm <- 999   # Methods state 999; an earlier run used 99

perm_bers <- numeric(n_perm)

cat("\nRunning", n_perm, "permutations...\n")

for(p in 1:n_perm) {

  if(p %% 10 == 0) cat("Permutation", p, "of", n_perm, "\n")

  Y_perm <- sample(Y)

  perm_model <- tryCatch({
    block.splsda(
      X      = X,
      Y      = Y_perm,
      ncomp  = 2,
      keepX  = list(microbiome = 20, metabolite = 10),
      design = design
    )
  }, error = function(e) NULL)

  if(is.null(perm_model)) { perm_bers[p] <- NA; next }

  perm_perf <- tryCatch({
    perf(
      perm_model,
      validation  = "Mfold",
      folds       = 5,
      nrepeat     = 3,
      progressBar = FALSE
    )
  }, error = function(e) NULL)

  if(is.null(perm_perf)) { perm_bers[p] <- NA; next }

  perm_bers[p] <- perm_perf$MajorityVote.error.rate$mahalanobis.dist[
    "Overall.BER", "comp1"
  ]
}

perm_bers_clean <- perm_bers[!is.na(perm_bers)]

# p-value: proportion of permutations with BER <= real BER
perm_p <- (sum(perm_bers_clean <= real_ber) + 1) / (length(perm_bers_clean) + 1)  # +1: smallest possible p = 1/(n+1)

cat("\n=== Permutation Test Results ===\n")
cat("Real BER:       ", round(real_ber, 4), "\n")
cat("Permutation p:  ", round(perm_p, 4), "\n")
cat("N permutations: ", length(perm_bers_clean), "\n")

# Plot
perm_plot_df <- data.frame(BER = perm_bers_clean)

p_perm <- ggplot(perm_plot_df, aes(x = BER)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "grey40") +
  geom_vline(xintercept = real_ber, colour = "#E69F00",
             linewidth = 1.2, linetype = "dashed") +
  annotate("text",
           x     = real_ber - 0.01,
           y     = Inf,
           label = paste0("Observed BER = ", round(real_ber, 3),
                          "\np = ", round(perm_p, 4)),
           hjust = 1, vjust = 1.5, size = 3.5) +
  labs(
    x   = "Balanced error rate (permuted)",
    y   = "Count",
    tag = "D"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 14))

print(p_perm)


tiff(file.path(out_dir, "FigureS_DIABLO_permutation.tiff"),
     width = 5, height = 4, units = "in", res = 600, compression = "lzw")
print(p_perm)
dev.off()
