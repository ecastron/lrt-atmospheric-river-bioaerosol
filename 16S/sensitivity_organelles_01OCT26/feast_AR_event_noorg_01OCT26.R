# feast_AR_event_noorg_01OCT26.R
# FEAST (Shenhav et al. 2019) per-sample source tracking for the 24 daily Risopatron air samples (sinks) against four
# source pools: Patagonian air (Coriolis samples at PA_POR, PN_SG and Puerto Williams), Patagonian terrestrial (bulk
# soil at the same sites), Risopatron local soil and rhizosphere, and Southern Ocean air (Malard et al. 2022).
# FEAST is run with different_sources_flag = 0, 1,000 EM iterations and rarefaction coverage of 1,000 reads.
# The Patagonian air samples were collected in December 2021, before the event, so the pool is not contemporaneous.
# Input  : 16S/sensitivity_organelles_01OCT26/ps_merged_withTree_noorg_01OCT26.RDS (the organelle-free object from build_ps_noorg_01OCT26.R)
# Outputs: per-sink, per-source proportions, time-series plot and Wilcoxon summary (written to the output folder set
#          in the script)
# Run    : Rscript 16S/sensitivity_organelles_01OCT26/feast_AR_event_noorg_01OCT26.R

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(patchwork); library(FEAST)
  library(lubridate)
})

setwd("./16S")
date_tag <- "noorg_01OCT26"
fig_dir  <- "sensitivity_organelles_01OCT26"
dir.create(fig_dir, showWarnings = FALSE)

set.seed(100)
ps <- readRDS("sensitivity_organelles_01OCT26/ps_merged_withTree_noorg_01OCT26.RDS")
sd_all <- as(sample_data(ps), "data.frame")

## ---- Masks ---------------------------------------------------------------
pat_regions <- c("PA_POR", "PN_SG", "Puerto Williams")

sink_mask       <- with(sd_all,
  dataset == "LRT" & !is_control & region == "Risopatron" &
  sampler == "Coriolis_mu")

src_patair_mask <- with(sd_all,
  dataset == "LRT" & !is_control & region %in% pat_regions &
  sampler == "Coriolis_mu")

src_patland_mask <- with(sd_all,
  dataset == "LRT" & !is_control & region %in% pat_regions &
  sampler %in% c("soil", "rhizosphere"))

src_rasoil_mask <- with(sd_all,
  dataset == "LRT" & !is_control & region == "Risopatron" &
  sampler %in% c("soil", "rhizosphere"))

src_malard_mask <- with(sd_all,
  dataset == "Malard2022" & source2 == "SO Air")

cat("Sink set     :", sum(sink_mask),       "RA Coriolis_mu samples\n")
cat("Patagonia Air:", sum(src_patair_mask), "samples\n")
cat("Patagonia TS :", sum(src_patland_mask),"samples\n")
cat("RA soil/rhizo:", sum(src_rasoil_mask), "samples\n")
cat("Malard SO Air:", sum(src_malard_mask), "samples\n")

## ---- Subset phyloseq -----------------------------------------------------
keep <- sink_mask | src_patair_mask | src_patland_mask |
        src_rasoil_mask | src_malard_mask
ps_feast <- prune_samples(keep, ps)
ps_feast <- prune_samples(sample_sums(ps_feast) >= 1000, ps_feast)
ps_feast <- prune_taxa(taxa_sums(ps_feast) > 0, ps_feast)
cat("\nFEAST input before ASV filter: ", nsamples(ps_feast), " samples × ",
    ntaxa(ps_feast), " ASVs\n", sep = "")

## ASV prevalence filter: FEAST's EM scales poorly with very high ASV
## counts and singletons cannot be reliably traced to a source anyway.
## Keep only ASVs present in >=2 samples with >=10 total reads.
otu_chk <- as(otu_table(ps_feast), "matrix")
if (taxa_are_rows(ps_feast)) otu_chk <- t(otu_chk)
prev_per_asv <- colSums(otu_chk > 0)
sum_per_asv  <- colSums(otu_chk)
keep_asv     <- prev_per_asv >= 2 & sum_per_asv >= 10
ps_feast     <- prune_taxa(keep_asv, ps_feast)
cat("FEAST input after ASV filter:  ", nsamples(ps_feast), " samples × ",
    ntaxa(ps_feast), " ASVs (prevalence ≥2 samples, ≥10 total reads)\n",
    sep = "")

## ---- Build FEAST metadata ------------------------------------------------
sample_to_env <- function(sn) {
  case_when(
    sn %in% sample_names(ps)[sink_mask]        ~ "RA_Air_sink",
    sn %in% sample_names(ps)[src_patair_mask]  ~ "Patagonia_Air",
    sn %in% sample_names(ps)[src_patland_mask] ~ "Patagonia_Terrestrial",
    sn %in% sample_names(ps)[src_rasoil_mask]  ~ "Risopatron_Local",
    sn %in% sample_names(ps)[src_malard_mask]  ~ "Malard_SO_Air",
    TRUE                                       ~ NA_character_)
}

## FEAST convention when different_sources_flag = 0:
##   - sources share their pool across all sinks, but must have id = NA
##   - each sink must have its own unique integer id (1..n_sinks)
meta_feast <- data.frame(
  Env        = sample_to_env(sample_names(ps_feast)),
  SourceSink = ifelse(sample_names(ps_feast) %in% sample_names(ps)[sink_mask],
                      "Sink", "Source"),
  row.names  = sample_names(ps_feast),
  stringsAsFactors = FALSE)
sink_names <- rownames(meta_feast)[meta_feast$SourceSink == "Sink"]
meta_feast$id <- NA_integer_
meta_feast$id[match(sink_names, rownames(meta_feast))] <- seq_along(sink_names)

cat("\nFEAST metadata table:\n")
print(table(meta_feast$Env, meta_feast$SourceSink, useNA = "ifany"))

## ---- OTU matrix (samples × ASVs, integer counts) -------------------------
otu_mat <- as(otu_table(ps_feast), "matrix")
if (taxa_are_rows(ps_feast)) otu_mat <- t(otu_mat)
otu_mat <- round(otu_mat)
storage.mode(otu_mat) <- "integer"
stopifnot(identical(rownames(otu_mat), rownames(meta_feast)))

## ---- Run FEAST -----------------------------------------------------------
cat("\nRunning FEAST (this may take a few minutes) ...\n")
feast_out <- FEAST(
  C                       = otu_mat,
  metadata                = meta_feast,
  different_sources_flag  = 0,
  dir_path                = fig_dir,
  outfile                 = paste0("feast_AR_event_", date_tag),
  EM_iterations           = 1000,
  COVERAGE                = 1000
)
cat("FEAST done.\n")

## FEAST writes:
##   {fig_dir}/{outfile}_source_contributions_matrix.txt
##   {fig_dir}/{outfile}_source_contributions_matrix.tsv
## and returns the same matrix as an R object.

if (is.list(feast_out) && !is.null(feast_out$Proportions)) {
  props <- feast_out$Proportions
} else {
  props <- as.matrix(feast_out)
}
cat("\nFEAST proportions matrix:\n")
print(round(head(props, 5), 3))
cat("(rows = sinks, cols = sources + Unknown)\n")

write.table(props,
            file = file.path(fig_dir,
              paste0("feast_AR_event_", date_tag, "_proportions.tsv")),
            sep = "\t", quote = FALSE, col.names = NA)

## ---- Build per-sample data frame with phase ------------------------------
sd_sink <- sd_all[rownames(sd_all) %in% rownames(props), ]
prop_df <- as.data.frame(props) |>
  tibble::rownames_to_column("sample_id") |>
  mutate(date  = as.Date(sd_sink$date[match(sample_id, rownames(sd_sink))]),
         phase = case_when(
           date <= as.Date("2022-02-06")                                ~ "Before",
           date >= as.Date("2022-02-07") & date <= as.Date("2022-02-08") ~ "During",
           TRUE                                                         ~ "After"),
         phase = factor(phase, c("Before", "During", "After")),
         ar2   = if_else(phase == "Before", "Before", "During+After")) |>
  arrange(date)

long_df <- prop_df |>
  pivot_longer(cols = !c(sample_id, date, phase, ar2),
               names_to = "source", values_to = "proportion") |>
  mutate(source = factor(source,
    levels = c("Patagonia_Air", "Patagonia_Terrestrial",
               "Malard_SO_Air", "Risopatron_Local", "Unknown")))

## ---- Visual --------------------------------------------------------------
pal_src <- c(
  "Patagonia_Air"        = "#C44E52",  # red: terrestrial source
  "Patagonia_Terrestrial"= "#8C2A2A",  # darker red
  "Malard_SO_Air"        = "#4C72B0",  # blue: marine
  "Risopatron_Local"     = "#937860",  # brown: local
  "Unknown"              = "#BDBDBD")  # gray: unattributed

ar_rect <- annotate("rect",
                    xmin = as.Date("2022-02-05"),
                    xmax = as.Date("2022-02-11"),
                    ymin = -Inf, ymax = Inf,
                    fill = "#C44E52", alpha = 0.10)

pA <- ggplot(long_df,
             aes(x = date, y = proportion, fill = source)) +
  ar_rect +
  geom_area(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal_src, name = "Source") +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "FEAST proportion",
       subtitle = paste0("FEAST source proportions per daily Risopatron Air sample ",
                         "(AR window shaded, Feb 5-11)")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

pB <- ggplot(long_df,
             aes(x = phase, y = proportion, fill = source)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = pal_src, guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ source, ncol = 5, scales = "free_y") +
  labs(x = NULL, y = "FEAST proportion",
       subtitle = "Per-source distribution by phase") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9, face = "bold"))

fig <- (pA / pB) +
  plot_layout(heights = c(1.2, 0.9)) +
  plot_annotation(
    title = paste0("FEAST source apportionment of Risopatron Air across the AR window"),
    subtitle = paste0("Stegen-style time-resolved microbial source tracking; ",
                      "Patagonia samples not contemporaneous with AR (lower-bound estimate)")) &
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, paste0("feast_AR_event_", date_tag, "_timeseries.pdf")),
       fig, width = 12, height = 9)
ggsave(file.path(fig_dir, paste0("feast_AR_event_", date_tag, "_timeseries.png")),
       fig, width = 12, height = 9, dpi = 300)

## ---- Wilcoxon Before vs During+After per source --------------------------
cat("\n=== Wilcoxon Before vs During+After per source ===\n")
wilcox_rows <- list()
for (src in levels(long_df$source)) {
  sub <- long_df |> filter(source == src) |>
    mutate(g = factor(ar2, levels = c("Before", "During+After")))
  w <- suppressWarnings(wilcox.test(proportion ~ g, data = sub, exact = FALSE))
  mb <- median(sub$proportion[sub$g == "Before"],       na.rm = TRUE)
  ma <- median(sub$proportion[sub$g == "During+After"], na.rm = TRUE)
  cat(sprintf("%-22s  median Before = %.3f  D+A = %.3f  W = %.1f  p = %.4f\n",
              src, mb, ma, w$statistic, w$p.value))
  wilcox_rows[[src]] <- data.frame(source = src,
                                   median_Before = round(mb, 3),
                                   median_DurAft = round(ma, 3),
                                   W = w$statistic,
                                   p = signif(w$p.value, 3))
}
wilcox_df <- do.call(rbind, wilcox_rows)
wilcox_df$p_BH <- signif(p.adjust(wilcox_df$p, "BH"), 3)
row.names(wilcox_df) <- NULL

sink(file.path(fig_dir,
               paste0("feast_AR_event_", date_tag, "_stats.txt")), split = FALSE)
cat("FEAST AR-event source apportionment: ", date_tag, "\n", sep = "")
cat("===========================================================\n")
cat("Sinks: ", sum(sink_mask), " daily RA Coriolis_mu samples\n", sep = "")
cat("Source pools (different_sources_flag = 0): Patagonia_Air, ",
    "Patagonia_Terrestrial, Risopatron_Local, Malard_SO_Air\n", sep = "")
cat("Rarefaction COVERAGE = 1000 reads/sample, EM 1000 iters.\n")
cat("\nWilcoxon Before vs During+After (BH-adjusted):\n")
print(wilcox_df, row.names = FALSE)
sink()

cat("\nDone. Outputs:\n")
for (s in c("_proportions.tsv","_timeseries.pdf",
            "_timeseries.png","_stats.txt"))
  cat("  ", file.path(fig_dir, paste0("feast_AR_event_", date_tag, s)), "\n")
