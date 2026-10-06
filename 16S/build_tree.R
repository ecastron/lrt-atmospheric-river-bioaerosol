# build_tree.R
# Aligns the merged ASV sequences with MAFFT, infers an approximately maximum-likelihood tree with FastTree (GTR),
# midpoint-roots it and attaches it to the merged phyloseq object.
# Inputs : 16S/ps_merged_<DATE>.RDS (merge_LRT_ACE.R); mafft on the PATH; FastTree binary at ./FastTree
# Outputs: 16S/asvs_merged_<DATE>.fasta, 16S/aln_merged_<DATE>.fasta, 16S/tree_merged_<DATE>.nwk,
#          16S/ps_merged_withTree_<DATE>.RDS
# Run    : Rscript 16S/build_tree.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(ape)
})

set.seed(100)
base     <- "./16S"
date_tag <- "07MAY26"   # tag of the outputs used in the manuscript; downstream scripts read *_07MAY26 files

ps_path     <- file.path(base, paste0("ps_merged_", date_tag, ".RDS"))
fasta_path  <- file.path(base, paste0("asvs_merged_", date_tag, ".fasta"))
aln_path    <- file.path(base, paste0("aln_merged_",  date_tag, ".fasta"))
tree_path   <- file.path(base, paste0("tree_merged_", date_tag, ".nwk"))
out_ps_path <- file.path(base, paste0("ps_merged_withTree_", date_tag, ".RDS"))

mafft_bin    <- "mafft"
fasttree_bin <- "./FastTree"

n_threads <- max(1, parallel::detectCores() - 1)

stopifnot(file.exists(ps_path),
          nchar(Sys.which(mafft_bin)) > 0,
          file.exists(fasttree_bin))

# --- 1) extract ASV sequences from merged phyloseq --------------------
ps <- readRDS(ps_path)
refs <- refseq(ps)
stopifnot(!is.null(refs), length(refs) == ntaxa(ps))
cat("ASVs to align: ", length(refs), "\n", sep = "")
writeXStringSet(refs, fasta_path)

# --- 2) align with MAFFT ----------------------------------------------
# --auto picks an algorithm based on input size; for short 16S V4 ASVs
# this typically lands on FFT-NS-2, which is fast and sufficient for
# tree-building at the genus/family resolution we care about.
cat("Running MAFFT...\n")
t0 <- Sys.time()
status <- system2(
  mafft_bin,
  args = c("--auto", "--thread", as.character(n_threads),
           "--quiet", shQuote(fasta_path)),
  stdout = aln_path, stderr = ""
)
stopifnot(status == 0, file.exists(aln_path), file.size(aln_path) > 0)
cat("MAFFT done in ", format(Sys.time() - t0), "\n", sep = "")

# --- 3) FastTree (GTR + CAT, double precision via -gamma omitted to
#         keep it fast; switch to -gamma if you want likelihoods that
#         are usable for model selection downstream) -------------------
cat("Running FastTree...\n")
t0 <- Sys.time()
status <- system2(
  fasttree_bin,
  args = c("-nt", "-gtr", "-quiet", shQuote(aln_path)),
  stdout = tree_path
)
stopifnot(status == 0, file.exists(tree_path), file.size(tree_path) > 0)
cat("FastTree done in ", format(Sys.time() - t0), "\n", sep = "")

# --- 4) read tree, midpoint-root, attach to phyloseq ------------------
tr <- read.tree(tree_path)
# Midpoint rooting (project convention): needs phangorn
# but it's a tiny operation. Fall back to ape::root at first tip if
# phangorn is unavailable.
if (requireNamespace("phangorn", quietly = TRUE)) {
  tr <- phangorn::midpoint(tr)
} else {
  warning("phangorn not installed: leaving tree unrooted")
}
tr <- ape::ladderize(tr)

# Sanity: tree tips must match ps taxa names
stopifnot(setequal(tr$tip.label, taxa_names(ps)))

phy_tree(ps) <- tr
saveRDS(ps, out_ps_path)

cat("\nWrote ", out_ps_path, "\n", sep = "")
cat("  taxa:    ", ntaxa(ps),    "\n", sep = "")
cat("  samples: ", nsamples(ps), "\n", sep = "")
cat("  tree tips:", Ntip(phy_tree(ps)), "\n", sep = "")
