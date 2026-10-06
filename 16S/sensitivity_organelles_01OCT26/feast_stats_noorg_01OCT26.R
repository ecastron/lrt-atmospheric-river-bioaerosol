# feast_stats_noorg_01OCT26.R
# Before versus During+After medians of the FEAST source proportions per pool with Benjamini-Hochberg-adjusted
# Wilcoxon tests.
# Inputs : FEAST proportion tables, 16S/ps_merged_withTree_07MAY26.RDS (sample dates)
# Run    : Rscript 16S/sensitivity_organelles_01OCT26/feast_stats_noorg_01OCT26.R

suppressPackageStartupMessages({library(phyloseq);library(dplyr);library(tidyr);library(readr)})
setwd("./16S")
ps <- readRDS("ps_merged_withTree_07MAY26.RDS"); sd <- as(sample_data(ps), "data.frame")
run <- function(f) {
  prop <- read_tsv(f, show_col_types = FALSE)
  prop$date <- as.Date(sd$date[match(prop$sample_id, rownames(sd))])
  prop$ar2 <- factor(ifelse(prop$date <= as.Date("2022-02-06"), "Before", "During+After"), c("Before","During+After"))
  prop |> pivot_longer(-c(sample_id, date, ar2), names_to = "source", values_to = "proportion") |>
    group_by(source) |>
    summarize(medB = round(100*median(proportion[ar2=="Before"]),1), medA = round(100*median(proportion[ar2=="During+After"]),1),
              p = suppressWarnings(wilcox.test(proportion ~ ar2, exact = FALSE)$p.value), .groups = "drop") |>
    mutate(p_BH = signif(p.adjust(p, "BH"), 3))
}
cat("Original (figures/feast_AR_event_12MAY26_proportions.tsv):\n"); print(as.data.frame(run("figures/feast_AR_event_12MAY26_proportions.tsv")))
cat("\nNo organelles:\n"); print(as.data.frame(run("sensitivity_organelles_01OCT26/feast_AR_event_noorg_01OCT26_proportions.tsv")))
