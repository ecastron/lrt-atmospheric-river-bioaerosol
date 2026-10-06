# merge_LRT_ACE.R
# Merges the LRT and ACE sequence tables at the ASV level, removes chimeras on the union, assigns SILVA 138.1
# taxonomy once on the merged set and builds the merged phyloseq object.
# Inputs : 16S/seqtab_ACE_07MAY26.rds (ACE_dada2.R), 16S/ps_16S_LRT_withSeqs.RDS (LRT.R),
#          16S/metadata_merged_07MAY26.csv (harmonize_metadata.R), SILVA 138.1 training and species files
# Outputs: 16S/seqtab_merged_<DATE>.rds, 16S/taxa_merged_<DATE>.rds, 16S/ps_merged_<DATE>.RDS
# Run    : Rscript 16S/merge_LRT_ACE.R

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(Biostrings)
  library(dplyr)
  library(readr)
})

set.seed(100)
base     <- "./16S"
date_tag <- "07MAY26"   # tag of the outputs used in the manuscript; downstream scripts read *_07MAY26 files

seqtab_ace_path <- file.path(base, "seqtab_ACE_07MAY26.rds")
ps_lrt_path     <- file.path(base, "ps_16S_LRT_withSeqs.RDS")
metadata_path   <- file.path(base, "metadata_merged_07MAY26.csv")
silva_train     <- file.path(base, "silva_nr99_v138.1_train_set.fa")
silva_species   <- file.path(base, "silva_species_assignment_v138.1.fa")

stopifnot(file.exists(seqtab_ace_path),
          file.exists(ps_lrt_path),
          file.exists(metadata_path))

# --- LRT: phyloseq -> seqtab matrix ------------------------------------
# Need a samples-by-ASV count matrix where colnames are the actual
# ASV nucleotide sequences (the format DADA2 expects).
ps_lrt <- readRDS(ps_lrt_path)

otu <- otu_table(ps_lrt)
if (taxa_are_rows(otu)) otu <- t(otu)
otu <- as(otu, "matrix")

# If taxa_names are already DNA sequences, use them directly. If they
# are short labels (e.g. "ASV1") and sequences live in refseq(ps), swap.
looks_like_seq <- all(grepl("^[ACGTN]+$", taxa_names(ps_lrt)))
if (looks_like_seq) {
  seqtab_lrt <- otu
} else {
  refs <- refseq(ps_lrt)
  stopifnot(!is.null(refs),
            setequal(names(refs), taxa_names(ps_lrt)))
  colnames(otu) <- as.character(refs[colnames(otu)])
  seqtab_lrt <- otu
}
mode(seqtab_lrt) <- "integer"

# --- ACE: load seqtab from ACE_dada2.R ---------------------------------
seqtab_ace <- readRDS(seqtab_ace_path)

# Sanity: ASV sequence lengths should align between datasets. Tiny
# variation is OK (DADA2 sometimes produces +/-1 length variants),
# but a wholesale offset means primers / truncLen mismatch.
len_lrt <- nchar(getSequences(seqtab_lrt))
len_ace <- nchar(getSequences(seqtab_ace))
cat("LRT ASV length: median ", median(len_lrt),
    " range ", paste(range(len_lrt), collapse = "-"), "\n", sep = "")
cat("ACE ASV length: median ", median(len_ace),
    " range ", paste(range(len_ace), collapse = "-"), "\n", sep = "")
stopifnot(abs(median(len_lrt) - median(len_ace)) <= 2)

# Optional: prune to the dominant amplicon length to drop spurious
# variants. Keeps sequences within the expected V4 size class.
target_len <- median(c(len_lrt, len_ace))
seqtab_lrt <- seqtab_lrt[, abs(nchar(colnames(seqtab_lrt)) - target_len) <= 5, drop = FALSE]
seqtab_ace <- seqtab_ace[, abs(nchar(colnames(seqtab_ace)) - target_len) <= 5, drop = FALSE]

# --- merge -------------------------------------------------------------
# Same sample ID in both -> sum (shouldn't happen here, but explicit).
seqtab_merged <- mergeSequenceTables(seqtab_lrt, seqtab_ace,
                                     repeats = "sum")
cat("Merged: ", nrow(seqtab_merged), " samples x ",
    ncol(seqtab_merged), " ASVs\n", sep = "")

# --- chimera removal on the union --------------------------------------
seqtab_nochim <- removeBimeraDenovo(seqtab_merged,
                                    method      = "consensus",
                                    multithread = TRUE,
                                    verbose     = TRUE)
cat("Post-chimera: ", ncol(seqtab_nochim), " ASVs (",
    round(sum(seqtab_nochim) / sum(seqtab_merged), 3),
    " of reads retained)\n", sep = "")

saveRDS(seqtab_nochim,
        file.path(base, paste0("seqtab_merged_", date_tag, ".rds")))

# --- taxonomy ----------------------------------------------------------
taxa <- assignTaxonomy(seqtab_nochim, silva_train, multithread = TRUE)
if (!is.na(silva_species) && file.exists(silva_species)) {
  taxa <- addSpecies(taxa, silva_species)
}
saveRDS(taxa, file.path(base, paste0("taxa_merged_", date_tag, ".rds")))

# --- assemble phyloseq -------------------------------------------------
meta <- read_csv(metadata_path, show_col_types = FALSE) |>
  as.data.frame()
rownames(meta) <- meta$sample_id

# Keep only samples we have data for (intersection)
keep <- intersect(rownames(seqtab_nochim), rownames(meta))
cat("Samples in seqtab: ", nrow(seqtab_nochim),
    " | in metadata: ", nrow(meta),
    " | intersection: ", length(keep), "\n", sep = "")

# Re-key ASVs from full sequences to short labels (ASV0001, ...) for
# tractable downstream work; sequences preserved in refseq().
asv_seqs           <- colnames(seqtab_nochim)
asv_ids            <- sprintf("ASV%05d", seq_along(asv_seqs))
seqtab_short       <- seqtab_nochim[keep, , drop = FALSE]
colnames(seqtab_short) <- asv_ids
taxa_short         <- taxa
rownames(taxa_short) <- asv_ids
refs               <- DNAStringSet(asv_seqs)
names(refs)        <- asv_ids

ps_merged <- phyloseq(
  otu_table(seqtab_short, taxa_are_rows = FALSE),
  tax_table(taxa_short),
  sample_data(meta[keep, , drop = FALSE]),
  refs
)

saveRDS(ps_merged,
        file.path(base, paste0("ps_merged_", date_tag, ".RDS")))

message("Done.")
message("  seqtab_merged_", date_tag, ".rds")
message("  taxa_merged_",   date_tag, ".rds")
message("  ps_merged_",     date_tag, ".RDS")
