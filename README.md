# Code for: Stability as a dynamic process: microbial and metabolic succession in a continuous in vitro gut reactor

Analysis code for the manuscript (submitted to *Journal of Applied Microbiology*).
Four reactor runs (TR1-TR4), ascending (AC), transverse (TC) and descending colon (DC), days 1–10.

## Folder layout

```
data/
  microbiome/TR1/alldaysTR1.csv     genus-level abundance table (one row per sample)
  microbiome/TR1/alldaysTR1m.csv    sample metadata (id, Sample.name = TR1_<day>_<compartment>)
  ...                               same for TR2, TR3, TR4
  Metabolites.csv                   Run, Day, Compartment + metabolite concentrations
R/                                  analysis scripts (run in order)
hydrodynamics/                      MATLAB scripts for HRT and turnover (Supplementary Text S1)
figures/                            created by the scripts
```

Raw 16S rRNA gene reads are in the NCBI Sequence Read Archive (BioProject PRJNA_______).

## Scripts and what they produce

| Script | Produces | Values reported in the paper |
|---|---|---|
| `R/01_fig1_alpha_diversity.R` | Figure 1 (Hill number, q = 1) | Kruskal–Wallis: compartment p = 2.04 × 10⁻⁸; day p = 1.66 × 10⁻⁵ |
| `R/02_fig2_figS3_beta_diversity_permanova.R` | Figure 2 (Bray–Curtis PCoA; day-to-day dissimilarity), Figure S3 (Jaccard) | PERMANOVA stage and run R²; Jaccard R²; windowed PERMANOVA (days 6–10); day-to-day Bray–Curtis after day 5 (0.238 ± 0.151) |
| `R/03_fig3_genus_barplots.R` | Figure 3 (20 most abundant genera) | — |
| `R/04_fig4_fig5_metabolites_lmm.R` | Figures 4–5 (z-scored metabolite heatmaps) | Linear mixed-effects models (stage effect on NH₃, H₂S, FRAP) |
| `R/05_fig6_fig7_diablo.R` | Figures 6–7 (DIABLO score and loading plots), DIABLO permutation histogram | Balanced error rate, permutation p |
| `R/06_fig8_spearman.R` | Figure 8 (Spearman heatmap) | FDR-significant genus–metabolite correlations |

### Hydrodynamics (MATLAB, Supplementary Text S1)

| Script | Produces |
|---|---|
| `hydrodynamics/HRT_feeding_starvation_schedule_annotated.m` | Figure S1 (compartment volumes over the feeding/starvation cycle, with HRT variables) |
| `hydrodynamics/HRT_washout_nominal_vs_effective_5V.m` | Figure S2 (washout: ideal CSTR vs real pulse feeding) |
| `hydrodynamics/HRT_calculation_RealOperation_5V.m` | Nominal and effective HRT of each compartment (Table S1) |

Run each script from the repository root. Script 06 uses objects created by 05, so run 05 and 06 in the same R session.

## Software

R 4.3.2; packages: tidyverse, vegan, patchwork, RColorBrewer, ggtext, lme4, lmerTest, mixOmics, igraph.
Exact versions: see `sessionInfo.txt`.
