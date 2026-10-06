# decontam_controls_11MAY26.R
# Contamination assessment following Fierer et al. 2025: summarizes the negative controls (read depth, composition,
# Bray-Curtis PCoA of buffer, Coriolis and deposition controls) and runs decontam::isContaminant (prevalence method)
# with the buffer controls as the negative reference, at thresholds 0.1 and 0.5.
# Inputs : 16S/ps_16S_LRT_withSeqs.RDS, 16S/metadata_LRT3.csv, metadata_LRT.csv,
#          16S/figures/ar_attribution_v3_summary_08MAY26.tsv (optional cross-check)
# Outputs: 16S/figures/decontam_summary_11MAY26.txt, decontam_contaminants_11MAY26.tsv,
#          decontam_control_diag_11MAY26.tsv, decontam_controls_pcoa_11MAY26.pdf
# Run    : Rscript 16S/decontam_controls_11MAY26.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(decontam)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(vegan)
})

setwd("./16S")
date_tag <- "11MAY26"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

## ---- 1. Load --------------------------------------------------------------
ps   <- readRDS("ps_16S_LRT_withSeqs.RDS")
meta <- read.csv("metadata_LRT3.csv", sep = ";", header = TRUE,
                 row.names = 1, check.names = FALSE)
sample_data(ps) <- sample_data(meta)

cat("Loaded:", nsamples(ps), "samples ×", ntaxa(ps), "ASVs\n")

sd <- as.data.frame(sample_data(ps))
cat("\nCategory2 table:\n"); print(table(sd$Category2, useNA = "ifany"))
cat("\nSource2 table:\n");   print(table(sd$Source2,   useNA = "ifany"))

## Identify controls by Source2: the metadata_LRT.csv uses lower-case
## tags (buffer_control etc.); metadata_LRT3 uses Category2 = "Control".
## Be defensive and try both.
ctrl_re <- "(?i)control|blank"
sd$is_control <- grepl(ctrl_re, sd$Source2) | grepl(ctrl_re, sd$Category2) |
                 grepl(ctrl_re, sd$Category)

## Pull a finer label so we can distinguish buffer / coriolis / deposition.
## Map from metadata_LRT.csv (lower-case sample_type) to provide the split.
meta_long <- read.csv("../metadata_LRT.csv", header = TRUE,
                      stringsAsFactors = FALSE)
sd$ctrl_type <- meta_long$sample_type[match(rownames(sd), meta_long$sample_name)]
sd$ctrl_type[!sd$is_control] <- NA_character_
sample_data(ps) <- sample_data(sd)

cat("\nControl type × region (sites listed in metadata_LRT.csv):\n")
print(table(sd$ctrl_type, sd$region, useNA = "ifany"))

## ---- 2. Control read-depth diagnostic -------------------------------------
depth_df <- data.frame(sample = sample_names(ps),
                       reads  = sample_sums(ps),
                       ctrl_type = sd$ctrl_type,
                       region    = sd$region,
                       is_control = sd$is_control,
                       Source    = sd$Source)

ctrl_depth <- depth_df |>
  filter(is_control) |>
  arrange(ctrl_type, region)
cat("\nRead depth per control sample:\n")
print(ctrl_depth, row.names = FALSE)

cat("\nMedian depth per control type:\n")
print(ctrl_depth |> group_by(ctrl_type) |>
        summarize(n = n(), median = median(reads), min = min(reads),
                  max = max(reads), .groups = "drop"))

## sea_water_control: real environmental seawater: EXCLUDE from any
## negative-reference role. Keep separately.
sd$is_neg_ref <- sd$is_control &
                 !is.na(sd$ctrl_type) &
                 sd$ctrl_type == "buffer_control"
sample_data(ps) <- sample_data(sd)

## ---- 3. Coriolis vs deposition: do they cluster? --------------------------
## If coriolis_controls (esp. Yelcho) compositionally resemble
## deposition_controls more than buffer_controls, that supports treating
## them as deposition controls (researchers' mis-labelling hypothesis).
ps_ctrl <- prune_samples(sd$is_control & !is.na(sd$ctrl_type) &
                         sd$ctrl_type != "sea_water_control", ps)
ps_ctrl <- prune_taxa(taxa_sums(ps_ctrl) > 0, ps_ctrl)
cat("\nControl-only subset:", nsamples(ps_ctrl), "samples ×",
    ntaxa(ps_ctrl), "ASVs\n")

if (nsamples(ps_ctrl) >= 6) {
  ps_ctrl_rel <- transform_sample_counts(ps_ctrl, function(x) {
    s <- sum(x); if (s == 0) x else x / s
  })
  bc <- phyloseq::distance(ps_ctrl_rel, method = "bray")
  meta_ctrl <- data.frame(sample_data(ps_ctrl_rel))
  meta_ctrl$ctrl_type <- factor(meta_ctrl$ctrl_type)
  meta_ctrl$region    <- factor(meta_ctrl$region)

  ## PERMANOVA: does control_type explain composition?
  set.seed(100)
  ad_type <- adonis2(bc ~ ctrl_type, data = meta_ctrl,
                     permutations = 9999, by = "margin")
  cat("\nPERMANOVA: composition ~ ctrl_type (all controls)\n")
  print(ad_type)

  ## Yelcho-only test (the hypothesis-of-interest site)
  yelcho_idx <- meta_ctrl$region == "Yelcho"
  if (sum(yelcho_idx) >= 4 &&
      length(unique(meta_ctrl$ctrl_type[yelcho_idx])) >= 2) {
    bc_y <- as.dist(as.matrix(bc)[yelcho_idx, yelcho_idx])
    set.seed(100)
    ad_y <- adonis2(bc_y ~ ctrl_type, data = meta_ctrl[yelcho_idx, ],
                    permutations = 9999, by = "margin")
    cat("\nPERMANOVA: composition ~ ctrl_type (Yelcho-only)\n")
    print(ad_y)
  } else {
    cat("\nYelcho-only PERMANOVA skipped (insufficient samples or 1 type only)\n")
  }

  ## PCoA plot
  pcoa <- cmdscale(bc, k = 2, eig = TRUE)
  pdat <- data.frame(PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2],
                     ctrl_type = meta_ctrl$ctrl_type,
                     region    = meta_ctrl$region,
                     sample    = rownames(meta_ctrl))
  eig <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)
  p <- ggplot(pdat, aes(PC1, PC2, color = ctrl_type, shape = region)) +
    geom_point(size = 3, alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = sample), size = 2.5,
                             max.overlaps = 30, show.legend = FALSE) +
    labs(x = paste0("PCoA1 (", eig[1], "%)"),
         y = paste0("PCoA2 (", eig[2], "%)"),
         title = "Negative controls: Bray-Curtis PCoA",
         subtitle = "Do coriolis_controls cluster with deposition_controls or with buffer_controls?") +
    theme_bw(base_size = 11)
  ggsave(file.path(fig_dir, paste0("decontam_controls_pcoa_", date_tag, ".pdf")),
         p, width = 7, height = 5)
  cat("\nWrote: ", file.path(fig_dir, paste0("decontam_controls_pcoa_",
                                             date_tag, ".pdf")), "\n")
}

## ---- 4. decontam (prevalence) ---------------------------------------------
## Negative reference = buffer_control only (true method blanks).
## Method = prevalence (no DNA concentration data available, so
## frequency method is not applicable).
##
## Threshold interpretation (decontam docs):
##   default 0.1: taxon is a contaminant if its prevalence-based
##                  score is < 0.1 (i.e. P-value of independence test)
##   0.5: recommended for "more aggressive" decontamination
##                  in high-contamination low-biomass settings
##
## We run both and report.
neg_vec <- sd$is_neg_ref
cat("\ndecontam negative reference = buffer_control (n =",
    sum(neg_vec), ")\n")

if (sum(neg_vec) < 2) stop("Need ≥2 buffer_control samples for prevalence method.")

contam_strict <- isContaminant(ps, method = "prevalence", neg = neg_vec,
                               threshold = 0.1, detailed = TRUE)
contam_aggr   <- isContaminant(ps, method = "prevalence", neg = neg_vec,
                               threshold = 0.5, detailed = TRUE)

cat("\nFlagged contaminant ASVs (threshold 0.1, default):",
    sum(contam_strict$contaminant), " / ", ntaxa(ps), "\n")
cat("Flagged contaminant ASVs (threshold 0.5, aggressive):",
    sum(contam_aggr$contaminant), " / ", ntaxa(ps), "\n")

## Reads removed per Source (air / soil / rhizosphere)
otu_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

sd_plain <- data.frame(sample = sample_names(ps),
                       Source = as.character(sd$Source),
                       region = as.character(sd$region),
                       is_control = as.logical(sd$is_control),
                       stringsAsFactors = FALSE)
rownames(sd_plain) <- sd_plain$sample
reads_removed <- function(flag_vec) {
  removed_reads_per_sample <- colSums(otu_mat[flag_vec, , drop = FALSE])
  total_reads_per_sample   <- colSums(otu_mat)
  pct <- 100 * removed_reads_per_sample /
              pmax(total_reads_per_sample, 1)
  tibble(sample = colnames(otu_mat),
         Source = sd_plain[colnames(otu_mat), "Source"],
         region = sd_plain[colnames(otu_mat), "region"],
         is_control = sd_plain[colnames(otu_mat), "is_control"],
         reads_removed = removed_reads_per_sample,
         reads_total   = total_reads_per_sample,
         pct_removed   = pct)
}

rr_strict <- reads_removed(contam_strict$contaminant)
rr_aggr   <- reads_removed(contam_aggr$contaminant)

cat("\nReads removed per Source × is_control (threshold 0.1):\n")
print(rr_strict |>
        group_by(Source, is_control) |>
        summarize(n = n(),
                  median_pct = round(median(pct_removed, na.rm = TRUE), 2),
                  max_pct    = round(max(pct_removed, na.rm = TRUE), 2),
                  .groups = "drop"))

cat("\nReads removed per Source × is_control (threshold 0.5):\n")
print(rr_aggr |>
        group_by(Source, is_control) |>
        summarize(n = n(),
                  median_pct = round(median(pct_removed, na.rm = TRUE), 2),
                  max_pct    = round(max(pct_removed, na.rm = TRUE), 2),
                  .groups = "drop"))

## ---- 5. Write contaminant table ------------------------------------------
tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
tax_df$asv <- rownames(tax_df)

contam_tbl <- contam_strict |>
  mutate(asv = rownames(contam_strict),
         flagged_strict = contaminant,
         p_strict       = p) |>
  select(asv, p_strict, flagged_strict, prev) |>
  left_join(contam_aggr |>
              mutate(asv = rownames(contam_aggr),
                     flagged_aggr = contaminant) |>
              select(asv, flagged_aggr),
            by = "asv") |>
  left_join(tax_df, by = "asv") |>
  arrange(p_strict)

write_tsv(contam_tbl, file.path(fig_dir,
          paste0("decontam_contaminants_", date_tag, ".tsv")))

## ---- 6. Cross-check vs AR-arrival taxa -----------------------------------
ar_path <- file.path(fig_dir, "ar_attribution_v3_summary_08MAY26.tsv")
if (file.exists(ar_path)) {
  ar_df <- read_tsv(ar_path, show_col_types = FALSE)
  cross <- ar_df |>
    left_join(contam_tbl |> select(asv, p_strict,
                                   flagged_strict, flagged_aggr),
              by = "asv") |>
    mutate(flagged_strict = ifelse(is.na(flagged_strict),
                                   FALSE, flagged_strict),
           flagged_aggr   = ifelse(is.na(flagged_aggr),
                                   FALSE, flagged_aggr))
  write_tsv(cross, file.path(fig_dir,
            paste0("decontam_AR_arrival_crosscheck_", date_tag, ".tsv")))

  cat("\n=== AR-arrival × contaminant cross-check ===\n")
  cat("AR-arrival ASVs total:        ", nrow(cross), "\n")
  cat("Flagged contaminant (0.1):    ", sum(cross$flagged_strict), "\n")
  cat("Flagged contaminant (0.5):    ", sum(cross$flagged_aggr),   "\n")

  cat("\nFlagged AR-arrival ASVs (strict 0.1):\n")
  print(cross |> filter(flagged_strict) |>
          select(asv, best_tax, attribution, p_strict, rel_durpost) |>
          arrange(desc(rel_durpost)), n = Inf)

  cat("\nFlagged AR-arrival ASVs by attribution (aggressive 0.5):\n")
  print(cross |> filter(flagged_aggr) |>
          count(attribution))
} else {
  cat("\nAR-arrival summary not found: ", ar_path, "\n")
}

## ---- 7. Write text summary -----------------------------------------------
sink(file.path(fig_dir, paste0("decontam_summary_", date_tag, ".txt")))
cat("decontam contamination assessment: ", date_tag, "\n")
cat("============================================\n")
cat("Object:        ps_16S_LRT_withSeqs.RDS\n")
cat("Method:        prevalence (decontam::isContaminant)\n")
cat("Negative ref:  buffer_control only (n = ", sum(neg_vec), ")\n",
    sep = "")
cat("Excluded:      sea_water_control (real environment)\n")
cat("\nFlagged ASVs:\n")
cat("  threshold 0.1 (default):    ", sum(contam_strict$contaminant), "\n")
cat("  threshold 0.5 (aggressive): ", sum(contam_aggr$contaminant),   "\n")
cat("\nTotal ASVs in object:        ", ntaxa(ps), "\n")
cat("\nFor reads-removed-per-Source breakdown see console output\n")
cat("Cross-check file: decontam_AR_arrival_crosscheck_", date_tag,
    ".tsv\n", sep = "")
sink()

cat("\nDone. Outputs in", fig_dir, "\n")
