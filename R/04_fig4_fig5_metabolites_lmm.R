# =====================================================================
# 04  Figures 4–5 — z-scored metabolite heatmaps; linear mixed-effects models
# Chapter 3 — "Stability as a dynamic process" (J Appl Microbiol submission)
# Inputs : data/Metabolites.csv (columns Run, Day, Compartment, metabolites)
# Outputs: figures/Figure4.pdf, figures/Figure5.pdf; Stage effects for NH3,
#          H2S, FRAP (Results 3.3)
# CHANGED (presentation only): metabolite labels not italic; NH3/H2S subscripts
#          via pdf_device (cairo_pdf or Quartz); Figure 5 panels share one colour scale.
# CHECK  : the LMM (Metabolite ~ Stage + Compartment + (1|Run)) is not yet
#          described in Methods.
# Run from the repository root:  source("R/<this file>")
# =====================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ---------------------------------------------------------------------
# 0. SETUP — adjust to your current machine's paths
# ---------------------------------------------------------------------
metab_file <- "data/Metabolites.csv"
out_dir    <- "figures"

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
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

metab <- read.csv(metab_file, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------
# 1. Metabolite groups
# ---------------------------------------------------------------------

env_metabs  <- c("NH3", "H2S", "FRAP")
scfa_metabs <- c("Acetic.Acid", "Propionic.Acid", "Butyric.Acid", "total.SCFA")
bcfa_metabs <- c("IsoButyric.Acid", "Isovaleric.Acid", "Isocaproic.Acid", "total.BCFA")

metabolite_labels <- c(
  "NH3"             = "NH\u2083",
  "H2S"             = "H\u2082S",
  "FRAP"            = "FRAP",
  "Acetic.Acid"     = "Acetate",
  "Propionic.Acid"  = "Propionate",
  "Butyric.Acid"    = "Butyrate",
  "total.SCFA"      = "Total SCFA",
  "IsoButyric.Acid" = "Isobutyrate",
  "Isovaleric.Acid" = "Isovalerate",
  "Isocaproic.Acid" = "Isocaproate",
  "total.BCFA"      = "Total BCFA"
)

# ---------------------------------------------------------------------
# 2. Long format + z-score (each metabolite standardised independently)
# ---------------------------------------------------------------------

metab_long <- metab %>%
  pivot_longer(cols = -c(Run, Day, Compartment),
               names_to = "Metabolite", values_to = "Value") %>%
  filter(!is.na(Value))

metab_long_z <- metab_long %>%
  group_by(Metabolite) %>%
  mutate(Value_z = as.numeric(scale(Value))) %>%
  ungroup() %>%
  mutate(
    Compartment = factor(Compartment, levels = c("AC", "TC", "DC")),
    Run         = factor(Run, levels = c("TR1", "TR2", "TR3", "TR4"))
  )

# ---------------------------------------------------------------------
# 3. Heatmap function (unchanged logic, no title baked in —
#    letting patchwork add the A/B tag instead, same fix as Figure 2)
# ---------------------------------------------------------------------

plot_heatmap <- function(data, metab_group, zlim_common = NULL) {

  ordered_metabs <- metab_group[metab_group %in% unique(data$Metabolite)]

  plot_data <- data %>%
    filter(Metabolite %in% ordered_metabs) %>%
    mutate(
      Metabolite = factor(Metabolite, levels = ordered_metabs),
      Label = metabolite_labels[as.character(Metabolite)]
    )

  zlim <- if (is.null(zlim_common)) max(abs(plot_data$Value_z), na.rm = TRUE) else zlim_common

  ggplot(plot_data, aes(x = Day, y = Compartment, fill = Value_z)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_vline(xintercept = 5.5, colour = "black", linewidth = 0.6, linetype = "dashed") +
    facet_grid(Label ~ Run, labeller = labeller(Label = label_value, Run = label_value)) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-zlim, zlim), name = "Z-score"
    ) +
    scale_x_continuous(breaks = 1:10, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text.x = element_text(face = "bold", size = 10),
      strip.text.y = element_text(face = "bold", size = 10, angle = 0, hjust = 0),
      strip.background = element_rect(fill = "grey92", colour = NA),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(margin = margin(t = 3)),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
      panel.spacing = unit(0.25, "lines"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.6, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    labs(x = "Day", y = NULL)
}

# ---------------------------------------------------------------------
# 4. FIGURE 4 — env metabolites (single figure, no A/B panels)
# ---------------------------------------------------------------------

figure4 <- plot_heatmap(metab_long_z, env_metabs)
figure4

ggsave(file.path(out_dir, "Figure4.pdf"), figure4, width = 10, height = 6.6, units = "in",
       device = pdf_device)   # renders NH3/H2S subscripts

# ---------------------------------------------------------------------
# 5. FIGURE 5 — SCFA (A) + BCFA (B), combined with patchwork
# ---------------------------------------------------------------------

# one shared colour scale for panels A and B, so patchwork shows a single key
zlim_5 <- max(abs(metab_long_z$Value_z[metab_long_z$Metabolite %in% c(scfa_metabs, bcfa_metabs)]), na.rm = TRUE)
p_scfa <- plot_heatmap(metab_long_z, scfa_metabs, zlim_5)
p_bcfa <- plot_heatmap(metab_long_z, bcfa_metabs, zlim_5)

figure5 <- p_scfa + p_bcfa +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

figure5

ggsave(file.path(out_dir, "Figure5.pdf"), figure5, width = 16, height = 8.8, units = "in",
       device = pdf_device)

cat("Saved:", file.path(out_dir, "Figure4.pdf"), "\n")
cat("Saved:", file.path(out_dir, "Figure5.pdf"), "\n")

# ---- Linear mixed-effects models: early vs late stage (Results 3.3) ----
# install.packages(c("lme4", "lmerTest"))  # if not already installed
library(lme4)
library(lmerTest)   # adds p-values to lmer output via Satterthwaite approximation
library(dplyr)

# Assumes you have a metabolite dataframe with columns:
# Run, Day, Compartment, FRAP, H2S, NH3 (adjust column names to match yours)

metab_stage <- metab %>%
  mutate(Stage = ifelse(Day <= 5, "Early", "Late"),
         Stage = factor(Stage, levels = c("Early", "Late")),
         Run = factor(Run),
         Compartment = factor(Compartment))

run_mixed_model <- function(response_var) {
  formula_str <- paste0(response_var, " ~ Stage + Compartment + (1|Run)")
  model <- lmer(as.formula(formula_str), data = metab_stage)
  cat("\n=== ", response_var, " ===\n")
  print(summary(model))
  return(model)
}

model_FRAP <- run_mixed_model("FRAP")
model_H2S  <- run_mixed_model("H2S")
model_NH3  <- run_mixed_model("NH3")

# Extract just the Stage effect p-values cleanly
extract_stage_result <- function(model, name) {
  s <- summary(model)$coefficients
  data.frame(
    Metabolite = name,
    Estimate   = round(s["StageLate", "Estimate"], 4),
    SE         = round(s["StageLate", "Std. Error"], 4),
    df         = round(s["StageLate", "df"], 1),
    t_value    = round(s["StageLate", "t value"], 3),
    p_value    = round(s["StageLate", "Pr(>|t|)"], 4)
  )
}

stage_results <- bind_rows(
  extract_stage_result(model_FRAP, "FRAP"),
  extract_stage_result(model_H2S, "H2S"),
  extract_stage_result(model_NH3, "NH3")
)

print(stage_results)

# ---- Descriptive values reported in Results 3.3 (mean +/- SD) ----
msd <- function(x) sprintf("%.2f +/- %.2f (n = %d)", mean(x, na.rm = TRUE),
                           sd(x, na.rm = TRUE), sum(!is.na(x)))
desc_vars <- c("NH3", "H2S", "FRAP", "Acetic.Acid", "Propionic.Acid",
               "Butyric.Acid", "total.SCFA", "total.BCFA")

cat("\n===== By stage (all runs and compartments) =====\n")
for (v in desc_vars) for (st in c("Early", "Late"))
  cat(sprintf("%-15s %-5s %s\n", v, st, msd(metab_stage[[v]][metab_stage$Stage == st])))

cat("\n===== By stage x compartment =====\n")
for (v in desc_vars) for (st in c("Early", "Late")) for (cp in c("AC", "TC", "DC"))
  cat(sprintf("%-15s %-5s %-2s %s\n", v, st, cp,
              msd(metab_stage[[v]][metab_stage$Stage == st & metab_stage$Compartment == cp])))

cat("\n===== Total BCFA by run x stage =====\n")
for (r in levels(metab_stage$Run)) for (st in c("Early", "Late"))
  cat(sprintf("%-4s %-5s %s\n", r, st,
              msd(metab_stage$total.BCFA[metab_stage$Run == r & metab_stage$Stage == st])))
