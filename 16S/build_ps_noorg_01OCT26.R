# build_ps_noorg_01OCT26.R
# Removes chloroplast (Order == "Chloroplast") and mitochondrial (Family == "Mitochondria") ASVs from the merged
# 16S object. The result is the input for all bacterial analyses in the manuscript.
# Input  : 16S/ps_merged_withTree_07MAY26.RDS (build_tree.R)
# Output : 16S/ps_merged_withTree_noorg_01OCT26.RDS
# Run    : Rscript 16S/build_ps_noorg_01OCT26.R

suppressPackageStartupMessages({library(phyloseq); library(ape)})
base <- "./16S"
sd_dir <- file.path(base, "sensitivity_organelles_01OCT26")
ps <- readRDS(file.path(base, "ps_merged_withTree_07MAY26.RDS"))
tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
org <- (!is.na(tt$Order) & tt$Order == "Chloroplast") | (!is.na(tt$Family) & tt$Family == "Mitochondria")
cat("Organelle ASVs removed:", sum(org), "of", ntaxa(ps), "\n")
before <- sample_sums(ps)
ps_noorg <- prune_taxa(!org, ps)
after <- sample_sums(ps_noorg)
stopifnot(ape::is.rooted(phy_tree(ps_noorg)), ntaxa(ps_noorg) == ntaxa(ps) - sum(org))
md <- as(sample_data(ps), "data.frame")
cross <- before >= 1000 & after < 1000
cat("Samples crossing the 1000-read threshold after removal:", sum(cross), "\n")
if (any(cross)) print(data.frame(sample = names(before)[cross], group = paste(md$dataset, md$source2)[cross],
                                 before = before[cross], after = after[cross]))
saveRDS(ps_noorg, file.path(base, "ps_merged_withTree_noorg_01OCT26.RDS"))
