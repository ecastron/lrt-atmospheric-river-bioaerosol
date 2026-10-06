# its_r1only_finalize_29MAY26.R
# Adds taxonomy and the ASV tree to the decontaminated ITS object, removes the reagent yeast Komagataella and keeps
# kingdom Fungi only.
# Inputs : ITS/ps_LRT_ITS_R1only_decontam_pretax_29MAY26.RDS, ITS/taxa_r1only_decontam_29MAY26.rds,
#          ITS/asvs_r1only_decontam_29MAY26.nwk (VeryFastTree)
# Output : ITS/ps_LRT_ITS_R1only_decontam_withTree_29MAY26.RDS
# Run    : Rscript ITS/its_r1only_finalize_29MAY26.R

suppressPackageStartupMessages({
  library(phyloseq); library(ape); library(phangorn); library(dplyr)
})
set.seed(100)
base <- "."
its_dir <- file.path(base, "ITS"); date_tag <- "29MAY26"

ps   <- readRDS(file.path(its_dir, paste0("ps_LRT_ITS_R1only_decontam_pretax_", date_tag, ".RDS")))
taxa <- readRDS(file.path(its_dir, paste0("taxa_r1only_decontam_", date_tag, ".rds")))
tr   <- ape::read.tree(file.path(its_dir, paste0("asvs_r1only_decontam_", date_tag, ".nwk")))

# attach taxonomy (already keyed by fASV id)
taxa <- taxa[taxa_names(ps), , drop = FALSE]
tax_table(ps) <- tax_table(taxa)
cat("Pretax ASVs:", ntaxa(ps), "\n")
cat("Kingdom breakdown:\n"); print(table(taxa[, "Kingdom"], useNA = "ifany"))

# --- taxonomy-driven cleaning ------------------------------------------
sd  <- as(sample_data(ps), "data.frame")
ra  <- sd$region == "Risopatron" & sd$sampler == "Coriolis_mu" & !sd$is_control
otu <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)

is_kom   <- !is.na(taxa[, "Genus"]) & taxa[, "Genus"] == "g__Komagataella"
is_fungi <- !is.na(taxa[, "Kingdom"]) & taxa[, "Kingdom"] == "k__Fungi"
kom_ids  <- taxa_names(ps)[is_kom]

cat(sprintf("\nKomagataella ASVs: %d (%.1f%% of all reads; %.1f%% of Risopatron-air reads)\n",
            length(kom_ids),
            100 * sum(otu[, kom_ids]) / sum(otu),
            100 * sum(otu[ra, kom_ids]) / sum(otu[ra, ])))
cat(sprintf("Non-Fungi / kingdom-NA ASVs dropped: %d\n", sum(!is_fungi)))

keep <- is_fungi & !is_kom
ps   <- prune_taxa(keep, ps)
ps   <- prune_taxa(taxa_sums(ps) > 0, ps)
cat("After Fungi-only + Komagataella removal:", ntaxa(ps), "ASVs\n")

# --- tree: midpoint-root, ladderize, prune to surviving taxa -----------
tr <- phangorn::midpoint(tr); tr <- ape::ladderize(tr)
tr <- ape::keep.tip(tr, intersect(tr$tip.label, taxa_names(ps)))
ps <- prune_taxa(tr$tip.label, ps)              # align taxa to tree tips
ps <- merge_phyloseq(ps, phy_tree(tr))
stopifnot(ape::is.rooted(phy_tree(ps)))

out <- file.path(its_dir, paste0("ps_LRT_ITS_R1only_decontam_withTree_", date_tag, ".RDS"))
saveRDS(ps, out)
cat("\nFinal object:\n"); print(ps); cat("Saved:", out, "\n")

# --- sanity: Risopatron-air richness by AR window ----------------------
source(file.path(base, "16S", "palettes.R"))
ps_ra <- prune_samples(ra, ps); ps_ra <- prune_taxa(taxa_sums(ps_ra) > 0, ps_ra)
sdra  <- as(sample_data(ps_ra), "data.frame")
sdra$d <- as.Date(sdra$date)
sdra$window <- cut(sdra$d, c(as.Date("2000-01-01"), ar_event_start - 1,
                             ar_event_end, as.Date("2100-01-01")),
                   labels = c("Before", "During", "After"))
sdra$reads <- sample_sums(ps_ra)
sdra$rich  <- estimate_richness(ps_ra, measures = "Observed")$Observed
cat("\nRisopatron air (cleaned R1-only) richness by AR window:\n")
print(sdra |> group_by(window) |>
        summarize(n = n(), median_reads = median(reads),
                  median_rich = median(rich), max_rich = max(rich),
                  n_singletons = sum(rich <= 1), .groups = "drop") |>
        as.data.frame())
cat("DONE finalize\n")
