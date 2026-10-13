# =====================================================================
# run_all.R — runs scripts 01–06 in order and saves all printed output
# Usage: open this folder in RStudio (set it as working directory), then
#        source("run_all.R")
# Output: figures/  (all figures)
#         output/   (one .txt log per script + sessionInfo.txt)
# Note  : 06 uses objects created by 05, so they run in the same session.
#         05 (DIABLO, 999 permutations) is slow — leave it running.
# =====================================================================
install.packages(c("tidyverse","vegan","patchwork","RColorBrewer","scales",
                   "ggtext","igraph","lme4","lmerTest"))
install.packages("BiocManager"); BiocManager::install("mixOmics")
dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)

scripts <- c(
  "R/01_fig1_alpha_diversity.R",
  "R/02_fig2_figS3_beta_diversity_permanova.R",
  "R/03_fig3_genus_barplots.R",
  "R/04_fig4_fig5_metabolites_lmm.R",
  "R/05_fig6_fig7_diablo.R",
  "R/06_fig8_spearman.R"
)

for (s in scripts) {
  log_file <- file.path("output", sub("\\.R$", ".txt", basename(s)))
  cat("Running", s, "...\n")
  con <- file(log_file, open = "wt")
  sink(con); sink(con, type = "message")
  ok <- tryCatch({
    source(s, print.eval = TRUE)
    TRUE
  }, error = function(e) {
    cat("\n*** ERROR in", s, ":", conditionMessage(e), "\n")
    FALSE
  })
  sink(type = "message"); sink(); close(con)
  cat(if (ok) "  done ->" else "  ERROR, see", log_file, "\n")
  if (!ok) break
}

writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
cat("\nFinished. Figures in figures/, results in output/\n")
