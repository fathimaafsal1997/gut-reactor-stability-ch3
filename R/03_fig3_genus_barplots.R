# =====================================================================
# 03  Figure 3 — relative abundance of the 20 most abundant genera (2 x 2 grid)
# Chapter 3 — "Stability as a dynamic process" (J Appl Microbiol submission)
# Inputs : data/microbiome/TR*/...
# Outputs: figures/Figure3_genus_barplots.pdf/.tiff
# CHECK  : taxonomy_lookup is typed by hand — verify every Family/Phylum.
#          Mitsuokella was listed as Bacteroidota; corrected below to
#          Firmicutes/Selenomonadaceae (SILVA 138). Phylum names mix old
#          (Firmicutes) and new (Pseudomonadota) nomenclature — make consistent.
# CHECK  : x position (SampleNum) follows row order in the file, not Day —
#          confirm files are sorted by day.
# Run from the repository root:  source("R/<this file>")
# =====================================================================

#install.packages(c("tidyverse", "stringr", "RColorBrewer", "ggplot2"))
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(ggplot2)

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

  data <- left_join(metadata, taxa, by = "id")
  data$Run <- run_name

  data <- data %>%
    separate(Sample.name,
             into = c("Run","Day","Compartment"),
             sep = "_",
             remove = FALSE)

  return(data)
}

#############################################
# 3. LOAD ALL DATA (FOR GLOBAL COLOR MAP)
#############################################

all_data <- bind_rows(lapply(runs, load_run))


#############################################
# 4. LONG FORMAT + GENUS EXTRACTION
#############################################

long_all <- all_data %>%
  pivot_longer(
    cols = starts_with("d__"),
    names_to = "Taxonomy",
    values_to = "Abundance"
  ) %>%
  mutate(
    Genus = str_extract(Taxonomy, "g__[^;]+"),
    Genus = gsub("g__", "", Genus),
    Genus = ifelse(is.na(Genus) | Genus == "" | Genus == "Unassigned",
                   "Unclassified", Genus)
  ) %>%
  group_by(Sample.name) %>%
  mutate(RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

#############################################
# 5. GLOBAL TOP 20 GENERA (CRITICAL FIX)
#############################################

top_genera <- long_all %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(RelAbundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abund)) %>%
  slice(1:20) %>%
  pull(Genus)

print(top_genera)

#############################################
# 5A. VALIDATION CHECK — run BEFORE 5B
# Confirms taxonomy_lookup actually covers your real top 20 genera
#############################################

# (Run this after defining taxonomy_lookup below, then re-check.)
# missing_genera <- setdiff(top_genera, taxonomy_lookup$Genus)
# print(missing_genera)
# If missing_genera is not empty, add those genera + their correct
# Phylum/Family to taxonomy_lookup before proceeding.

#############################################
# 5B. TAXONOMY LOOKUP TABLE (Phylum > Family > Genus)
#############################################

taxonomy_lookup <- tribble(
  ~Genus,                          ~Phylum,              ~Family,
  "Enterococcus",                  "Firmicutes",         "Enterococcaceae",
  "Bacteroides",                   "Bacteroidota",       "Bacteroidaceae",
  "Mitsuokella",                   "Firmicutes",         "Selenomonadaceae",   # CORRECTED (was Bacteroidota)
  "Bifidobacterium",               "Actinobacteriota",   "Bifidobacteriaceae",
  "Megasphaera",                   "Firmicutes",         "Veillonellaceae",
  "Clostridium_sensu_stricto_1",   "Firmicutes",         "Clostridiaceae",
  "Veillonella",                   "Firmicutes",         "Veillonellaceae",
  "Parabacteroides",               "Bacteroidota",       "Tannerellaceae",
  "Rhodobacter",                   "Pseudomonadota",     "Rhodobacteraceae",
  "Coprococcus",                   "Firmicutes",         "Lachnospiraceae",
  "Lactobacillus",                 "Firmicutes",         "Lactobacillaceae",
  "Lachnoclostridium",             "Firmicutes",         "Lachnospiraceae",
  ".Eubacterium._hallii_group",    "Firmicutes",         "Lachnospiraceae",
  "Escherichia.Shigella",          "Pseudomonadota",     "Enterobacteriaceae",
  "Blautia",                       "Firmicutes",         "Lachnospiraceae",
  "Faecalibacterium",              "Firmicutes",         "Oscillospiraceae",
  "Phascolarctobacterium",         "Firmicutes",         "Acidaminococcaceae",
  "Proteiniphilum",                "Bacteroidota",       "Dysgonomonadaceae",
  "Dorea",                         "Firmicutes",         "Lachnospiraceae",
  "Anaerostipes",                  "Firmicutes",         "Lachnospiraceae"
)
# NOTE: This is still a placeholder list from an earlier draft.
# CONFIRM it matches your actual top_genera before trusting the legend —
# run the check below right now.

missing_genera <- setdiff(top_genera, taxonomy_lookup$Genus)
cat("\nGenera in top_genera but MISSING from taxonomy_lookup:\n")
print(missing_genera)
cat("\nIf the above list is not empty, add these genera (with correct\n",
    "Phylum and Family) to taxonomy_lookup before proceeding further.\n")

# Clean up genus display names for legend
taxonomy_lookup <- taxonomy_lookup %>%
  mutate(
    Genus_display = Genus,
    Genus_display = gsub("_sensu_stricto_1", "", Genus_display),
    Genus_display = gsub("\\.", " ", Genus_display),
    Genus_display = gsub("^ ", "", Genus_display),
    Genus_display = gsub("Escherichia Shigella", "Escherichia/Shigella", Genus_display),
    Genus_display = gsub(" Eubacterium  hallii group", "Eubacterium hallii group", Genus_display)
  )

#############################################
# 6. FIXED COLOR MAP — SORTED BY PHYLUM > FAMILY > GENUS
#############################################

# Join top genera with taxonomy
top_genera_df <- data.frame(Genus = top_genera) %>%
  left_join(taxonomy_lookup, by = "Genus") %>%
  arrange(Phylum, Family, Genus)

# Create ordered genus levels for legend
ordered_genera <- top_genera_df$Genus

# Assign colors — grouped by Phylum using distinct palettes
phyla <- unique(top_genera_df$Phylum)

# Color palettes per phylum
phylum_palettes <- list(
  "Firmicutes"          = c(brewer.pal(9, "Blues")[3:9], 
                            brewer.pal(9, "RdPu")[3:8]),   # extends to 13 colors
  "Bacteroidota"        = brewer.pal(9, "Greens")[3:7],
  "Actinobacteriota"    = brewer.pal(9, "Purples")[4:6],
  "Pseudomonadota"      = brewer.pal(9, "Oranges")[4:6],
  "Thermodesulfobacteriota" = brewer.pal(9, "Reds")[4:6]
)

# Assign colors within each phylum group
genus_colors <- c()
for (phylum in phyla) {
  genera_in_phylum <- top_genera_df %>%
    filter(Phylum == phylum) %>%
    pull(Genus)
  palette <- phylum_palettes[[phylum]]
  if (is.null(palette)) palette <- brewer.pal(9, "Greys")[3:9]
  colors <- palette[1:length(genera_in_phylum)]
  genus_colors <- c(genus_colors, setNames(colors, genera_in_phylum))
}

genus_colors["Other"] <- "grey50"

# UPDATED: legend now includes Phylum alongside Family
legend_labels <- top_genera_df %>%
  mutate(label = paste0(Genus_display, " (", Family, "; ", Phylum, ")")) %>%
  pull(label, name = Genus)

legend_labels["Other"] <- "Other"


print(top_genera)
#############################################
# 6 (UPDATED). LEGEND LABELS — italic genus via ggtext markdown
#############################################

library(ggtext)

legend_labels <- top_genera_df %>%
  mutate(label = paste0("*", Genus_display, "*", " (", Family, "; ", Phylum, ")")) %>%
  pull(label, name = Genus)

legend_labels["Other"] <- "Other"

#############################################
# 7C. COMBINED PLOT — 2x2 GRID (TR1/TR2 top, TR3/TR4 bottom)
# Each panel keeps compartment facets; single shared legend on right,
# 1 column, genus names italicized via ggtext
#############################################

library(patchwork)

process_run_for_grid <- function(run_name) {

  data <- load_run(run_name)

  long <- data %>%
    pivot_longer(
      cols = starts_with("d__"),
      names_to = "Taxonomy",
      values_to = "Abundance"
    ) %>%
    mutate(
      Genus = str_extract(Taxonomy, "g__[^;]+"),
      Genus = gsub("g__", "", Genus),
      Genus = ifelse(is.na(Genus) | Genus == "" | Genus == "Unassigned",
                     "Unclassified", Genus)
    ) %>%
    group_by(Sample.name) %>%
    mutate(RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
    ungroup()

  long <- long %>%
    mutate(Genus2 = ifelse(Genus %in% top_genera, Genus, "Other"))

  long <- long %>%
    group_by(Compartment) %>%
    mutate(SampleNum = as.numeric(factor(Sample.name,
                                         levels = unique(Sample.name)))) %>%
    ungroup()

  long$Genus2 <- factor(long$Genus2, levels = c(ordered_genera, "Other"))
  long$Compartment <- factor(long$Compartment, levels = c("AC","TC","DC"))

  ggplot(long,
         aes(x = SampleNum, y = RelAbundance, fill = Genus2)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Compartment, scales = "free_x", nrow = 1) +
    scale_fill_manual(
      values = genus_colors,
      labels = legend_labels,
      name = "Genus (Family; Phylum)",
      guide = guide_legend(ncol = 1)
    ) +
    theme_bw(base_size = 11) +
    labs(title = run_name, x = "Day", y = "Relative abundance") +
    scale_x_continuous(breaks = 1:10) +
    theme(
      legend.text = element_markdown(size = 8),   # markdown-aware for italics
      legend.key.size = unit(0.4, "cm"),
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 7)
    )
}

p_TR1 <- process_run_for_grid("TR1")
p_TR2 <- process_run_for_grid("TR2")
p_TR3 <- process_run_for_grid("TR3")
p_TR4 <- process_run_for_grid("TR4")

p_grid <- (p_TR1 | p_TR2) / (p_TR3 | p_TR4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.text = element_markdown(size = 8))  # ensure collected legend also uses markdown

p_grid

#############################################
# SAVE
#############################################

ggsave(file.path("figures", "Figure3_genus_barplots.tiff"),
       p_grid, width = 14, height = 11,
       dpi = 600, compression = "lzw")

ggsave(file.path("figures", "Figure3_genus_barplots.pdf"),
       p_grid, width = 14, height = 11)

