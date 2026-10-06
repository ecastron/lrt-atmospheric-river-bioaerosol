# decontam_betaNTI_01OCT26.R
# Between-community phylogenetic turnover (betaNTI; Stegen et al. 2013) among the 24 daily Risopatron air samples,
# computed on the full and on the decontaminated object, to test whether the atmospheric-river signal depends on
# contaminant ASVs. betaNTI = (observed abundance-weighted betaMNTD - mean null) / sd null, with the null from
# tip-label shuffling on the cophenetic matrix; null permutations are drawn serially from set.seed(100) so results
# are reproducible. betaNTI > +2 indicates variable selection, < -2 homogeneous selection.
# Inputs : 16S/ps_merged_withTree_noorg_01OCT26.RDS, 16S/figures/decontam_contaminants_11MAY26.tsv
# Outputs: 16S/figures/decontam_betaNTI_01OCT26.{txt,tsv,pdf}
# Run    : Rscript 16S/decontam_betaNTI_01OCT26.R

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr); library(readr)
  library(picante); library(parallel); library(ggplot2); library(patchwork)
  library(lubridate)
})

setwd("./16S")
date_tag <- "01OCT26"
contam_tag <- "11MAY26"   # decontam flags (presence/absence based; unchanged by organelle removal)
fig_dir  <- "figures"
NULL_RUNS <- 999
NCORES    <- max(1, parallel::detectCores() - 2)
cat("Using", NCORES, "cores;", NULL_RUNS, "null reps\n")

set.seed(100)
ps <- readRDS("ps_merged_withTree_noorg_01OCT26.RDS")
contam <- read_tsv(file.path(fig_dir,
                  paste0("decontam_contaminants_", contam_tag, ".tsv")),
                  show_col_types = FALSE) |>
  filter(flagged_strict == TRUE)
refseq_chr <- as.character(refseq(ps))
asv_by_seq <- setNames(names(refseq_chr), refseq_chr)
contam_ids <- unname(asv_by_seq[contam$asv])
contam_ids <- contam_ids[!is.na(contam_ids)]
ps_decon <- prune_taxa(setdiff(taxa_names(ps), contam_ids), ps)

ra_subset <- function(p) {
  md <- as(sample_data(p), "data.frame")
  k <- md$dataset == "LRT" & !md$is_control &
       md$region == "Risopatron" & md$sampler == "Coriolis_mu"
  p2 <- prune_samples(k, p)
  p2 <- prune_taxa(taxa_sums(p2) > 0, p2)
  sd <- as(sample_data(p2), "data.frame") |>
    mutate(date  = as.Date(date),
           phase = case_when(
             date <= as.Date("2022-02-06") ~ "Before",
             date <= as.Date("2022-02-08") ~ "During",
             TRUE                          ~ "After"),
           phase = factor(phase, c("Before","During","After")),
           ar2   = if_else(phase == "Before", "Before", "During+After"))
  sample_data(p2) <- sample_data(sd)
  p2
}

ra_full  <- ra_subset(ps)
ra_decon <- ra_subset(ps_decon)

## ---- βNTI implementation -------------------------------------------------
## otu: samples × taxa (rows × cols)
## cophen: square cophenetic distance matrix with row/colnames = taxa
compute_bMNTD <- function(otu, cophen, weighted = TRUE) {
  as.matrix(picante::comdistnt(otu, cophen,
                               abundance.weighted = weighted,
                               exclude.conspecifics = FALSE))
}

compute_bNTI <- function(p, runs = NULL_RUNS, ncores = NCORES,
                         weighted = TRUE, label = "") {
  cat("  [", label, "] Preparing OTU + cophenetic ...\n", sep = "")
  otu <- as(otu_table(p), "matrix")
  if (taxa_are_rows(p)) otu <- t(otu)              # samples × taxa
  tr <- phy_tree(p)
  cophen <- cophenetic(tr)

  cat("  [", label, "] Observed βMNTD ...\n", sep = "")
  obs <- compute_bMNTD(otu, cophen, weighted = weighted)

  cat("  [", label, "] Null βMNTD (", runs, " reps, ", ncores,
      " cores) ...\n", sep = "")
  ## Permutations are drawn serially from one seed before mclapply, so the null
  ## does not depend on the number of cores (drawing inside forked workers is not reproducible)
  set.seed(100)
  perm_list <- replicate(runs, sample(colnames(cophen)), simplify = FALSE)
  null_list <- mclapply(perm_list, function(rand_names) {
    rc <- cophen
    rownames(rc) <- rand_names
    colnames(rc) <- rand_names
    compute_bMNTD(otu, rc, weighted = weighted)
  }, mc.cores = ncores, mc.preschedule = TRUE)

  null_arr <- simplify2array(null_list)            # n × n × runs
  mu  <- apply(null_arr, c(1, 2), mean)
  sg  <- apply(null_arr, c(1, 2), sd)
  bnti <- (obs - mu) / sg
  bnti[is.na(bnti)] <- 0                            # diagonal
  diag(bnti) <- 0
  list(obs = obs, bnti = bnti)
}

cat("\n--- FULL ---\n")
res_full  <- compute_bNTI(ra_full,  label = "FULL")
cat("\n--- DECON ---\n")
res_decon <- compute_bNTI(ra_decon, label = "DECON")

## ---- Build long-form pair table ------------------------------------------
long_from <- function(bnti, sdm, tag) {
  m <- bnti
  pairs <- which(upper.tri(m), arr.ind = TRUE)
  data.frame(
    A      = rownames(m)[pairs[, 1]],
    B      = colnames(m)[pairs[, 2]],
    bNTI   = m[upper.tri(m)],
    set    = tag,
    A_phase = sdm[rownames(m)[pairs[, 1]], "ar2"],
    B_phase = sdm[colnames(m)[pairs[, 2]], "ar2"],
    stringsAsFactors = FALSE)
}

sd_full  <- as(sample_data(ra_full),  "data.frame")
sd_decon <- as(sample_data(ra_decon), "data.frame")
df <- bind_rows(long_from(res_full$bnti,  sd_full,  "FULL"),
                long_from(res_decon$bnti, sd_decon, "DECON")) |>
  mutate(pair_type = case_when(
    A_phase == "Before" & B_phase == "Before"             ~ "Before-Before",
    A_phase != "Before" & B_phase != "Before"             ~ "D+A-D+A",
    TRUE                                                  ~ "Before-D+A"),
    pair_type = factor(pair_type,
                       c("Before-Before","Before-D+A","D+A-D+A")))

write_tsv(df, file.path(fig_dir,
          paste0("decontam_betaNTI_", date_tag, ".tsv")))

## ---- Summary by phase-pair × set -----------------------------------------
process_summary <- function(d) {
  d |> group_by(set, pair_type) |>
    summarize(n        = n(),
              median   = round(median(bNTI), 2),
              mean     = round(mean(bNTI),   2),
              q25      = round(quantile(bNTI, 0.25), 2),
              q75      = round(quantile(bNTI, 0.75), 2),
              hetero_pct = round(100 * mean(bNTI >   2), 1),  # variable sel.
              homo_pct   = round(100 * mean(bNTI <  -2), 1),  # homogeneous sel.
              stoch_pct  = round(100 * mean(abs(bNTI) <= 2), 1),
              .groups = "drop")
}

cat("\n=== βNTI by phase-pair and decontam set ===\n")
ps_summary <- process_summary(df)
print(ps_summary, n = Inf)

## ---- Wilcoxon: Before-Before vs Before-D+A ------------------------------
cat("\n=== Wilcoxon: Before-Before vs Before-D+A ===\n")
for (s in c("FULL","DECON")) {
  sub <- df |> filter(set == s, pair_type %in% c("Before-Before","Before-D+A"))
  w <- suppressWarnings(wilcox.test(bNTI ~ pair_type, data = sub,
                                    exact = FALSE))
  med_bb <- median(sub$bNTI[sub$pair_type == "Before-Before"])
  med_ba <- median(sub$bNTI[sub$pair_type == "Before-D+A"])
  cat(sprintf("  %-6s  W=%.0f  p=%.4f  median(BB)=%.2f  median(B-D+A)=%.2f\n",
              s, w$statistic, w$p.value, med_bb, med_ba))
}

## ---- Visual ---------------------------------------------------------------
df_plot <- df |> mutate(set = factor(set, c("FULL","DECON")))

pal <- c("Before-Before" = "#4C72B0",
         "Before-D+A"    = "#C44E52",
         "D+A-D+A"       = "#55A868")

p1 <- ggplot(df_plot, aes(pair_type, bNTI, fill = pair_type)) +
  geom_hline(yintercept = c(-2, 2), linetype = 2, color = "gray50") +
  geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
  geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
  facet_wrap(~ set) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = expression(beta*"NTI"),
       title = "βNTI by phase pair: Stegen process estimates",
       subtitle = "|βNTI| > 2 = selection; |βNTI| <= 2 = stochastic") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(fig_dir, paste0("decontam_betaNTI_", date_tag, ".pdf")),
       p1, width = 8, height = 5)

## ---- Text summary --------------------------------------------------------
sink(file.path(fig_dir, paste0("decontam_betaNTI_", date_tag, ".txt")),
     split = FALSE)
cat("Stegen βNTI at Risopatrón Air: ", date_tag, "\n")
cat("================================================\n")
cat("Abundance-weighted βMNTD; null = tip-label shuffle of cophenetic;\n")
cat("999 perms. Pair counts: Before-Before C(15,2) = 105;\n")
cat("Before-D+A 15×9 = 135; D+A-D+A C(9,2) = 36.\n\n")
print(ps_summary, n = Inf)
cat("\n--- Wilcoxon: Before-Before vs Before-D+A ---\n")
for (s in c("FULL","DECON")) {
  sub <- df |> filter(set == s, pair_type %in% c("Before-Before","Before-D+A"))
  w <- suppressWarnings(wilcox.test(bNTI ~ pair_type, data = sub,
                                    exact = FALSE))
  med_bb <- median(sub$bNTI[sub$pair_type == "Before-Before"])
  med_ba <- median(sub$bNTI[sub$pair_type == "Before-D+A"])
  cat(sprintf("  %-6s  W=%.0f  p=%.4g  median(BB)=%.2f  median(B-D+A)=%.2f\n",
              s, w$statistic, w$p.value, med_bb, med_ba))
}
sink()

cat("\nDone. Outputs:\n",
    " ", file.path(fig_dir, paste0("decontam_betaNTI_", date_tag, ".txt")), "\n",
    " ", file.path(fig_dir, paste0("decontam_betaNTI_", date_tag, ".tsv")), "\n",
    " ", file.path(fig_dir, paste0("decontam_betaNTI_", date_tag, ".pdf")), "\n",
    sep = "")
