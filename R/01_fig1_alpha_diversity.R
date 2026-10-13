# =====================================================================
# 01  Figure 1 — alpha diversity (Hill number, q = 1) + Kruskal–Wallis tests
# Chapter 3 — "Stability as a dynamic process" (J Appl Microbiol submission)
# Inputs : data/microbiome/TR*/alldaysTR*.csv, alldaysTR*m.csv
# Outputs: figures/Figure1_Hill.pdf; Kruskal–Wallis p-values (Results 3.2)
# NOTE   : omnibus Kruskal–Wallis only (no post hoc tests).
# Run from the repository root:  source("R/<this file>")
# =====================================================================

library(tidyverse)
library(stringr)
library(vegan)
library(patchwork)

#############################################
# 1. BASE PATH + RUNS
#############################################

base_path <- "data/microbiome"
runs <- c("TR1", "TR2", "TR3", "TR4")

dir.create("figures", showWarnings = FALSE)

#############################################
# 2. FUNCTION TO LOAD ONE RUN
#############################################

load_run <- function(run_name) {

  metadata <- read.csv(file.path(base_path, run_name,
                                 paste0("alldays", run_name, "m.csv")))

  taxa <- read.csv(file.path(base_path, run_name,
                             paste0("alldays", run_name, ".csv")))

  data <- dplyr::left_join(metadata, taxa, by = "id")
  data$Run <- run_name

  data <- data %>%
    tidyr::separate(Sample.name,
                    into = c("Run","Day","Compartment"),
                    sep = "_",
                    remove = FALSE)

  return(data)
}

#############################################
# 3. LOAD ALL DATA
#############################################

all_data <- dplyr::bind_rows(lapply(runs, load_run))

#############################################
# 4. LONG FORMAT + GENUS EXTRACTION
#############################################

long_all <- all_data %>%
  tidyr::pivot_longer(
    cols = tidyselect::starts_with("d__"),
    names_to = "Taxonomy",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(
    Genus = stringr::str_extract(Taxonomy, "g__[^;]+"),
    Genus = gsub("g__", "", Genus),
    Genus = ifelse(is.na(Genus) | Genus == "" | Genus == "Unassigned",
                   "Unclassified", Genus)
  ) %>%
  dplyr::group_by(Sample.name) %>%
  dplyr::mutate(RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  dplyr::ungroup()

#############################################
# 5. OTU MATRIX (genus-level, per sample)
#############################################

otu_matrix <- long_all %>%
  dplyr::select(Sample.name, Genus, Abundance) %>%
  dplyr::group_by(Sample.name, Genus) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Genus,
    values_from = Abundance,
    values_fill = 0
  )

otu_matrix <- as.data.frame(otu_matrix)

rownames(otu_matrix) <- otu_matrix$Sample.name
otu_matrix$Sample.name <- NULL

otu_matrix <- data.matrix(otu_matrix)

meta <- long_all %>%
  dplyr::select(Sample.name, Run, Day, Compartment) %>%
  dplyr::distinct()

otu_matrix[is.na(otu_matrix)] <- 0
otu_matrix <- otu_matrix[rowSums(otu_matrix) > 0, ]

#############################################
# 6. HILL NUMBER (q = 1) — exponential of Shannon
#############################################

hill1 <- exp(vegan::diversity(otu_matrix, index = "shannon"))

alpha_df <- data.frame(
  Sample.name = rownames(otu_matrix),
  Hill1 = hill1
) %>%
  dplyr::left_join(meta, by = "Sample.name")

alpha_df <- alpha_df %>%
  dplyr::mutate(
    Day = factor(Day, levels = 1:10, ordered = TRUE),
    Compartment = factor(Compartment, levels = c("AC","TC","DC")),
    Run = factor(Run, levels = c("TR1","TR2","TR3","TR4"))
  )

#############################################
# 7. FUNCTION: ONE RUN PLOT
#############################################

plot_run <- function(run_id) {

  df <- alpha_df %>% dplyr::filter(Run == run_id)

  ggplot(df, aes(x = Compartment, y = Hill1)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = Day),
                width = 0.15,
                size = 2,
                alpha = 0.8) +
    scale_x_discrete(limits = c("AC", "TC", "DC")) +
    coord_cartesian(ylim = c(0, 25)) +
    theme_bw(base_size = 13) +
    labs(
      title = run_id,
      x = "Compartment",
      y = "Hill number (q = 1)",
      color = "Day"
    ) +
    scale_color_viridis_d(option = "C")
}

#############################################
# 8. BUILD INDIVIDUAL PANELS
#############################################

p_TR1 <- plot_run("TR1")
p_TR2 <- plot_run("TR2")
p_TR3 <- plot_run("TR3")
p_TR4 <- plot_run("TR4")

#############################################
# 9. COMBINE WITH A SINGLE SHARED LEGEND
#############################################

final_plot <- (p_TR1 | p_TR2) / (p_TR3 | p_TR4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

final_plot

#############################################
# 10. STATS (unchanged from original)
#############################################

kruskal.test(Hill1 ~ Compartment, data = alpha_df)
kruskal.test(Hill1 ~ Day, data = alpha_df)

#############################################
# 11. SAVE
#############################################

out_dir <- "figures"

ggsave(file.path(out_dir, "Figure1_Hill.png"), plot = final_plot,
       width = 9, height = 5, dpi = 600)
ggsave(file.path(out_dir, "Figure1_Hill.pdf"), plot = final_plot,
       width = 9, height = 5)

cat("Saved:", file.path(out_dir, "Figure1_Hill.pdf"), "\n")
