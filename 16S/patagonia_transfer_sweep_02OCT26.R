# patagonia_transfer_sweep_02OCT26.R
# Identifies bacterial ASVs shared between Patagonian sources and Antarctic Peninsula air but absent or rare in
# local Antarctic soil, Southern Ocean air and negative controls.
# Criteria: detected in >= 2 Patagonian air or >= 2 Patagonian bulk soil samples; detected in >= 2 Antarctic air
# samples (Risopatron and Yelcho); prevalence <= 0.05 in Risopatron soil and rhizosphere; prevalence <= 0.02 in
# Southern Ocean air (Malard et al. 2022, over sea and over land); absent from all buffer controls; not flagged by
# decontam at threshold 0.5. ASVs are grouped by occurrence at Risopatron before and after the event.
# Inputs : 16S/ps_merged_withTree_noorg_01OCT26.RDS, metadata_LRT.csv, 16S/figures/decontam_contaminants_11MAY26.tsv
# Outputs: 16S/figures/patagonia_transfer_sweep_02OCT26_{asvs,genera,supptable}.tsv and _summary.txt
# Run    : Rscript 16S/patagonia_transfer_sweep_02OCT26.R

suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(readr); library(tidyr)})
setwd("./16S")
source("palettes.R")

ps  <- readRDS("ps_merged_withTree_noorg_01OCT26.RDS")
md  <- as(sample_data(ps), "data.frame")
otu <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)
pa  <- otu > 0
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

# sample types for controls come from the long metadata (as in decontam_controls_11MAY26.R)
ml <- read.csv("../metadata_LRT.csv", check.names = FALSE); names(ml)[1:4] <- c("x", "sample_name", "sample_type", "site")
stype <- ml$sample_type[match(rownames(md), ml$sample_name)]

pat_regions <- c("PA_POR", "PN_SG", "Puerto Williams")
lrt  <- md$dataset == "LRT" & !md$is_control
pools <- list(
  pat_air  = lrt & md$region %in% pat_regions & md$sampler == "Coriolis_mu",
  pat_soil = lrt & md$region %in% pat_regions & md$sampler %in% c("soil", "rhizosphere"),
  ant_air  = lrt & md$region %in% c("Risopatron", "Yelcho") & md$sampler == "Coriolis_mu",
  ris_before = lrt & md$region == "Risopatron" & md$sampler == "Coriolis_mu" & as.Date(md$date) <  ar_event_start,
  ris_da     = lrt & md$region == "Risopatron" & md$sampler == "Coriolis_mu" & as.Date(md$date) >= ar_event_start,
  yelcho_air = lrt & md$region == "Yelcho" & md$sampler == "Coriolis_mu",
  local_soil = lrt & md$region == "Risopatron" & md$sampler %in% c("soil", "rhizosphere"),
  yelcho_soil = lrt & md$region == "Yelcho" & md$sampler %in% c("soil", "rhizosphere"),
  so_air   = md$dataset == "Malard2022" & md$source2 %in% c("SO Air", "SO Air over land"),
  buffer   = md$dataset == "LRT" & !is.na(stype) & stype == "buffer_control")
n <- sapply(pools, sum); print(n)
stopifnot(n["pat_air"] == 20, n["ant_air"] == 54, n["buffer"] >= 1, n["pat_soil"] == 19, n["yelcho_soil"] == 6)

cnt  <- sapply(pools, function(i) colSums(pa[i, , drop = FALSE]))
prev <- sweep(cnt, 2, n, "/")
rel  <- otu / pmax(rowSums(otu), 1)
ant_reads_pct <- 100 * colMeans(rel[pools$ant_air, , drop = FALSE])
da_pct  <- 100 * colMeans(rel[pools$ris_da, , drop = FALSE])
bef_pct <- 100 * colMeans(rel[pools$ris_before, , drop = FALSE])

contam <- read_tsv("figures/decontam_contaminants_11MAY26.tsv", show_col_types = FALSE)
seqs <- setNames(as.character(refseq(ps)), taxa_names(ps))
flag05 <- seqs %in% contam$asv[contam$flagged_aggr %in% TRUE]

base_ok <- cnt[, "ant_air"] >= 2 & prev[, "local_soil"] <= 0.05 & prev[, "so_air"] <= 0.02 &
           cnt[, "buffer"] == 0 & !flag05
tier1 <- base_ok & cnt[, "pat_air"] >= 2
tier2 <- base_ok & (cnt[, "pat_air"] >= 2 | cnt[, "pat_soil"] >= 2)

asv_tab <- tibble(asv = taxa_names(ps), Phylum = tax$Phylum, Family = tax$Family, Genus = tax$Genus,
                  tier = ifelse(tier1, "1 (Patagonian air)", ifelse(tier2, "2 (Patagonian soil only)", NA)),
                  n_pat_air = cnt[, "pat_air"], n_pat_soil = cnt[, "pat_soil"], n_ant_air = cnt[, "ant_air"],
                  n_ris_before = cnt[, "ris_before"], n_ris_da = cnt[, "ris_da"], n_yelcho = cnt[, "yelcho_air"],
                  n_local_soil = cnt[, "local_soil"], n_yelcho_soil = cnt[, "yelcho_soil"], n_so_air = cnt[, "so_air"],
                  ant_air_mean_pct = ant_reads_pct, ris_before_pct = bef_pct, ris_da_pct = da_pct) |>
  filter(!is.na(tier)) |> arrange(tier, desc(n_ant_air))
write_tsv(asv_tab, "figures/patagonia_transfer_sweep_02OCT26_asvs.tsv")

gen_tab <- asv_tab |> mutate(Genus = coalesce(Genus, paste0("(", coalesce(Family, Phylum), ")"))) |>
  group_by(Genus) |>
  summarize(Phylum = first(Phylum), n_asv = n(), n_asv_tier1 = sum(startsWith(tier, "1")),
            max_n_pat_air = max(n_pat_air), max_n_ant_air = max(n_ant_air),
            ant_air_mean_pct = sum(ant_air_mean_pct), ris_before_pct = sum(ris_before_pct), ris_da_pct = sum(ris_da_pct),
            .groups = "drop") |>
  arrange(desc(n_asv_tier1), desc(max_n_ant_air))
write_tsv(gen_tab, "figures/patagonia_transfer_sweep_02OCT26_genera.tsv")

# Supplementary Table 1 export: group, source, sequence and full taxonomy per ASV.
# Groups at Risopatron: baseline = detected before the AR (split by whether mean relative abundance
# was higher During+After); new = not detected in any Before sample, detected During+After;
# Yelcho only = not detected at Risopatron (Yelcho sampling starts 8 Feb, no local soil there).
ranks <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(tax))
supp <- asv_tab |>
  mutate(group = case_when(
           n_ris_before > 0 & ris_da_pct >  ris_before_pct ~ "Baseline, increased during/after AR",
           n_ris_before > 0                                ~ "Baseline, not increased",
           n_ris_da > 0                                    ~ "New during/after AR",
           TRUE                                            ~ "Yelcho only"),
         patagonian_source = case_when(n_pat_air >= 2 & n_pat_soil >= 2 ~ "air and soil",
                                       n_pat_air >= 2 ~ "air", TRUE ~ "soil")) |>
  left_join(tibble(asv = taxa_names(ps), sequence = seqs) |>
              bind_cols(tax[, setdiff(ranks, c("Phylum", "Family", "Genus")), drop = FALSE]), by = "asv") |>
  select(asv, group, patagonian_source, Kingdom, Phylum, Class, Order, Family, Genus, any_of("Species"),
         n_pat_air, n_pat_soil, n_ris_before, n_ris_da, n_yelcho, n_local_soil, n_yelcho_soil, n_so_air,
         ris_before_pct, ris_da_pct, ant_air_mean_pct, sequence) |>
  arrange(factor(group, c("New during/after AR", "Baseline, increased during/after AR",
                          "Baseline, not increased", "Yelcho only")), desc(n_ris_da + n_ris_before + n_yelcho))
write_tsv(supp, "figures/patagonia_transfer_sweep_02OCT26_supptable.tsv")

sink("figures/patagonia_transfer_sweep_02OCT26_summary.txt")
cat("Patagonia -> Antarctic air transfer sweep (organelle-free object)\n"); print(n)
cat("\nASVs: tier 1 =", sum(tier1), "| tier 2 (incl. tier 1) =", sum(tier2), "of", ntaxa(ps), "\n")
cat("Share of Antarctic air reads in tier-1 ASVs: ", round(sum(ant_reads_pct[tier1]), 2), "%; tier 2: ", round(sum(ant_reads_pct[tier2]), 2), "%\n", sep = "")
cat("\nGenera (top 40):\n"); print(as.data.frame(gen_tab |> mutate(across(where(is.numeric), ~ round(.x, 3)))) |> head(40), row.names = FALSE)
sink()
cat(readLines("figures/patagonia_transfer_sweep_02OCT26_summary.txt"), sep = "\n")
