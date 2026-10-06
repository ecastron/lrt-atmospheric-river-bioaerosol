# decontam_AR_arrival_crosscheck_01OCT26.R
# Checks whether any ASV in the atmospheric-river arrival groups (ar_patagonia_attribution_v3_01OCT26.R) was flagged
# by decontam. Arrival groups are keyed by ASV ID and decontam results by sequence, so IDs are mapped to sequences
# with refseq(); every member ASV of each best-taxonomy group is checked.
# Inputs : 16S/ps_merged_withTree_noorg_01OCT26.RDS, 16S/figures/decontam_contaminants_11MAY26.tsv,
#          16S/figures/ar_attribution_v3_summary_01OCT26.tsv
# Output : 16S/figures/decontam_AR_arrival_crosscheck_01OCT26.tsv
# Run    : Rscript 16S/decontam_AR_arrival_crosscheck_01OCT26.R

suppressPackageStartupMessages({library(phyloseq); library(readr); library(dplyr)})
setwd("./16S")
source("palettes.R")
ps     <- readRDS("ps_merged_withTree_noorg_01OCT26.RDS")
contam <- read_tsv("figures/decontam_contaminants_11MAY26.tsv", show_col_types = FALSE)
ar_df  <- read_tsv("figures/ar_attribution_v3_summary_01OCT26.tsv", show_col_types = FALSE)
seqs   <- setNames(as.character(refseq(ps)), taxa_names(ps))
stopifnot(all(nchar(contam$asv) > 100), all(ar_df$asv %in% names(seqs)))

# best-taxonomy label per ASV, identical rule to ar_patagonia_attribution_v3_01OCT26.R
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
ranks <- intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), colnames(tax))
clean <- function(x) { x <- as.character(x); x[is.na(x) | x == "" | x == "NA" | grepl("^uncultured", x, ignore.case = TRUE)] <- NA; x }
best <- apply(tax[, ranks], 1, function(v) { v <- clean(v); k <- max(c(0, which(!is.na(v))))
  if (k == 0) "Unassigned" else paste0(ranks[k], ":", v[k]) })

# Risopatron air reads in During+After, per ASV (for read-weighted contamination share)
md <- as(sample_data(ps), "data.frame")
ra <- md$dataset == "LRT" & md$region == "Risopatron" & md$sampler == "Coriolis_mu" & !md$is_control
da <- ra & as.Date(md$date) >= ar_event_start
otu <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)
da_reads <- colSums(otu[da, , drop = FALSE])

fl <- contam |> select(asv, p_strict, flagged_strict, flagged_aggr)
member <- tibble(asv_id = names(best), best_tax = best, seq = seqs[names(best)], da_reads = da_reads[names(best)]) |>
  filter(best_tax %in% ar_df$best_tax) |>
  left_join(fl, by = c("seq" = "asv")) |>
  mutate(in_table = !is.na(p_strict), flagged_strict = coalesce(flagged_strict, FALSE), flagged_aggr = coalesce(flagged_aggr, FALSE))
grp <- member |> group_by(best_tax) |>
  summarize(n_asv = n(), n_in_table = sum(in_table),
            any_flag_01 = any(flagged_strict), any_flag_05 = any(flagged_aggr),
            pct_reads_flag_01 = 100 * sum(da_reads[flagged_strict]) / max(1, sum(da_reads)),
            pct_reads_flag_05 = 100 * sum(da_reads[flagged_aggr]) / max(1, sum(da_reads)), .groups = "drop")
rep <- ar_df |> mutate(seq = seqs[asv]) |> left_join(fl, by = c("seq" = "asv")) |>
  transmute(best_tax, attribution, rel_durpost, rep_asv = asv, rep_in_table = !is.na(p_strict),
            rep_flag_01 = coalesce(flagged_strict, FALSE), rep_flag_05 = coalesce(flagged_aggr, FALSE))
cross <- rep |> left_join(grp, by = "best_tax")
write_tsv(cross, "figures/decontam_AR_arrival_crosscheck_01OCT26.tsv")

cat("AR-arrival best-taxonomy groups:", nrow(cross), "\n")
cat("Representative ASVs found in decontam table:", sum(cross$rep_in_table), "\n")
cat("Representative ASV flagged: 0.1 =", sum(cross$rep_flag_01), "| 0.5 =", sum(cross$rep_flag_05), "\n")
cat("Groups with any member flagged: 0.1 =", sum(cross$any_flag_01), "| 0.5 =", sum(cross$any_flag_05), "\n")
tot <- member |> summarize(r = sum(da_reads), f1 = sum(da_reads[flagged_strict]), f5 = sum(da_reads[flagged_aggr]))
cat(sprintf("Share of During+After arrival-group reads in flagged ASVs: 0.1 = %.2f%% | 0.5 = %.2f%%\n", 100*tot$f1/tot$r, 100*tot$f5/tot$r))
cat("\nGroups with flagged members (0.5), by During+After abundance:\n")
print(as.data.frame(cross |> filter(any_flag_05) |> arrange(desc(rel_durpost)) |>
  select(best_tax, attribution, rel_durpost, rep_flag_01, rep_flag_05, n_asv, pct_reads_flag_01, pct_reads_flag_05) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
