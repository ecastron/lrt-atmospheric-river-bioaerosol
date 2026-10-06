# feast_AR_event_postprocess_noorg_01OCT26.R
# Aggregates the per-source-sample FEAST output into per-pool proportions for each sink and writes the time-series
# figure and Wilcoxon statistics (FEAST is not rerun).
# Inputs : FEAST output of feast_AR_event_noorg_01OCT26.R, the organelle-free object
# Outputs: per-pool proportion table, figure and summary (written to the output folder set in the script)
# Run    : Rscript 16S/sensitivity_organelles_01OCT26/feast_AR_event_postprocess_noorg_01OCT26.R

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(patchwork); library(lubridate)
})

setwd("./16S")
date_tag <- "noorg_01OCT26"
fig_dir  <- "sensitivity_organelles_01OCT26"

feast_file <- file.path(fig_dir,
   paste0("feast_AR_event_", date_tag,
          "_source_contributions_matrix.txt"))
stopifnot(file.exists(feast_file))

## ---- Read per-source-sample matrix (sinks × source-samples) ---------------
raw <- read.table(feast_file, sep = "\t", header = TRUE,
                  check.names = FALSE, row.names = 1)
cat("FEAST output:", nrow(raw), "sinks ×", ncol(raw), "source-sample columns\n")

## Column names look like "PAA1_Patagonia_Air", "100_Malard_SO_Air", or
## simply "Unknown". Extract the pool name as everything after the LAST
## underscore for non-Unknown columns.
col_pool <- ifelse(colnames(raw) == "Unknown",
                   "Unknown",
                   sub("^.*?_", "", colnames(raw)))
cat("Pool inventory among source-sample columns:\n")
print(table(col_pool))

## Aggregate by pool: row-wise sum across source-samples belonging to the
## same pool.
pools <- unique(col_pool)
agg <- sapply(pools, function(p) rowSums(raw[, col_pool == p, drop = FALSE]))
agg <- as.data.frame(agg, stringsAsFactors = FALSE)
agg$sample_id <- sub("_RA_Air_sink$", "", rownames(raw))

cat("\nAggregated proportions (first 5 sinks):\n")
print(round(head(agg[, pools], 5), 3))

write_tsv(agg,
          file.path(fig_dir,
            paste0("feast_AR_event_", date_tag, "_proportions.tsv")))

## ---- Attach date + phase --------------------------------------------------
ps <- readRDS("sensitivity_organelles_01OCT26/ps_merged_withTree_noorg_01OCT26.RDS")
sd <- as(sample_data(ps), "data.frame")

agg$date <- as.Date(sd$date[match(agg$sample_id, rownames(sd))])
agg$phase <- with(agg, case_when(
  date <= as.Date("2022-02-06")                                ~ "Before",
  date >= as.Date("2022-02-07") & date <= as.Date("2022-02-08") ~ "During",
  TRUE                                                         ~ "After"))
agg$phase <- factor(agg$phase, c("Before","During","After"))
agg$ar2   <- ifelse(agg$phase == "Before", "Before", "During+After")
agg <- agg[order(agg$date), ]

pool_levels <- c("Patagonia_Air", "Patagonia_Terrestrial",
                 "Malard_SO_Air", "Risopatron_Local", "Unknown")
long <- agg |>
  pivot_longer(cols = all_of(pool_levels),
               names_to = "source", values_to = "proportion") |>
  mutate(source = factor(source, levels = pool_levels))

## ---- Visual ---------------------------------------------------------------
pal_src <- c(
  "Patagonia_Air"         = "#C44E52",
  "Patagonia_Terrestrial" = "#8C2A2A",
  "Malard_SO_Air"         = "#4C72B0",
  "Risopatron_Local"      = "#937860",
  "Unknown"               = "#BDBDBD")

ar_rect <- annotate("rect",
                    xmin = as.Date("2022-02-05"),
                    xmax = as.Date("2022-02-11"),
                    ymin = -Inf, ymax = Inf,
                    fill = "#C44E52", alpha = 0.10)

pA <- ggplot(long, aes(date, proportion, fill = source)) +
  ar_rect +
  geom_area(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal_src, name = "Source") +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "FEAST proportion",
       subtitle = paste0("FEAST per-sample source proportions across the ",
                         "24 daily Risopatrón Air samples (AR shaded Feb 5-11)")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

pB <- ggplot(long, aes(phase, proportion, fill = source)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  scale_fill_manual(values = pal_src, guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ source, ncol = 5, scales = "free_y") +
  labs(x = NULL, y = "FEAST proportion",
       subtitle = "By phase (Before / During / After)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9, face = "bold"))

fig <- (pA / pB) +
  plot_layout(heights = c(1.2, 0.9)) +
  plot_annotation(
    title    = "FEAST source apportionment of Risopatrón Air across the AR window",
    subtitle = paste0("Patagonia samples are not contemporaneous with the AR ",
                      "(lower-bound estimate of Patagonian contribution)")) &
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir,
       paste0("feast_AR_event_", date_tag, "_timeseries.pdf")),
       fig, width = 12, height = 9)
ggsave(file.path(fig_dir,
       paste0("feast_AR_event_", date_tag, "_timeseries.png")),
       fig, width = 12, height = 9, dpi = 300)

## ---- Wilcoxon Before vs During+After per source ---------------------------
cat("\n=== Wilcoxon Before vs During+After per source ===\n")
rows <- list()
for (src in pool_levels) {
  sub <- long |> filter(source == src) |>
    mutate(g = factor(ar2, levels = c("Before", "During+After")))
  w <- suppressWarnings(wilcox.test(proportion ~ g, data = sub,
                                    exact = FALSE))
  mb <- median(sub$proportion[sub$g == "Before"],       na.rm = TRUE)
  ma <- median(sub$proportion[sub$g == "During+After"], na.rm = TRUE)
  cat(sprintf("%-22s  median Before = %.4f  D+A = %.4f  W = %.1f  p = %.4f\n",
              src, mb, ma, w$statistic, w$p.value))
  rows[[src]] <- data.frame(source = src,
                            median_Before = round(mb, 4),
                            median_DurAft = round(ma, 4),
                            W = w$statistic,
                            p = signif(w$p.value, 3))
}
wilcox_df <- do.call(rbind, rows)
wilcox_df$p_BH <- signif(p.adjust(wilcox_df$p, "BH"), 3)
rownames(wilcox_df) <- NULL

sink(file.path(fig_dir,
               paste0("feast_AR_event_", date_tag, "_stats.txt")),
     split = FALSE)
cat("FEAST AR-event source apportionment: ", date_tag, "\n", sep = "")
cat("===========================================================\n")
cat("Sinks: 24 daily RA Coriolis_mu samples\n")
cat("Source pools (different_sources_flag = 0):\n")
cat("  Patagonia_Air         (n = 20 source samples)\n")
cat("  Patagonia_Terrestrial (n = 19)\n")
cat("  Risopatron_Local      (n = 63)\n")
cat("  Malard_SO_Air         (n = 135)\n")
cat("Per-source-sample contributions written by FEAST, aggregated to ")
cat("per-pool here by row-sum.\n")
cat("Rarefaction COVERAGE = 1000 reads, EM 1000 iters.\n")
cat("ASV prevalence pre-filter: ≥2 samples and ≥10 total reads.\n")
cat("\nWilcoxon Before vs During+After (BH-adjusted):\n")
print(wilcox_df, row.names = FALSE)
sink()

cat("\nDone. Outputs:\n")
for (s in c("_proportions.tsv","_timeseries.pdf",
            "_timeseries.png","_stats.txt"))
  cat("  ", file.path(fig_dir, paste0("feast_AR_event_", date_tag, s)), "\n")
