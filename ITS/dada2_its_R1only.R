# dada2_its_R1only.R
# DADA2 processing of the ITS1 amplicons using forward reads only. Many fungal ITS1 amplicons are too long for
# 2 x 251 bp reads to overlap, so pair merging loses them; R1 covers the first 250 bp of ITS1, enough for UNITE
# assignment and ASV-level analysis. The conserved small-subunit sequence at the start of R1 is trimmed with cutadapt
# before this step.
# Inputs : trimmed ITS R1 FASTQ files (BioProject PRJNA1540148), 16S/ps_merged_withTree_07MAY26.RDS (sample metadata)
# Outputs: ITS sequence table and ITS/ps_LRT_ITS_R1only_08MAY26b.RDS
# Run    : Rscript ITS/dada2_its_R1only.R

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(Biostrings)
  library(dplyr)
  library(readr)
})

set.seed(100)
base       <- "."
its_dir    <- file.path(base, "ITS")
trim_dir   <- file.path(its_dir, "trimmed")
filt_dir   <- file.path(its_dir, "filtered_R1")
date_tag   <- "08MAY26b"
ps_16s     <- file.path(base, "16S", "ps_merged_withTree_07MAY26.RDS")

if (!dir.exists(filt_dir)) dir.create(filt_dir, recursive = TRUE)

fnFs <- sort(list.files(trim_dir, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE))
sample_names <- sub("_R1\\.fastq\\.gz$", "", basename(fnFs))
filtFs <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
names(filtFs) <- sample_names
cat("Samples:", length(sample_names), "\n")

# --- 1) filterAndTrim (R1 only) -------------------------------------
cat("\n[1/4] filterAndTrim ...\n")
out <- dada2::filterAndTrim(
  fnFs, filtFs,
  maxN     = 0,
  maxEE    = 2,
  truncQ   = 2,
  rm.phix  = TRUE,
  minLen   = 50,
  multithread = TRUE,
  compress = TRUE,
  verbose  = TRUE
)
keep <- file.exists(filtFs) & file.size(filtFs) > 100
filtFs <- filtFs[keep]; sample_names <- sample_names[keep]
cat("Surviving:", length(sample_names), "\n")

# --- 2) learnErrors -------------------------------------------------
cat("\n[2/4] learnErrors ...\n")
errF <- dada2::learnErrors(filtFs, multithread = TRUE, verbose = TRUE)

# --- 3) dada (pseudo-pool) -----------------------------------------
cat("\n[3/4] dada ...\n")
dadaFs <- dada2::dada(filtFs, err = errF,
                      multithread = TRUE, pool = "pseudo")

# --- 4) sequence table + chimera ------------------------------------
seqtab <- dada2::makeSequenceTable(dadaFs)
cat("\nSequence table:", dim(seqtab),
    " ASV length min/median/max:",
    min(nchar(colnames(seqtab))),
    median(nchar(colnames(seqtab))),
    max(nchar(colnames(seqtab))), "\n")

seqlens <- nchar(colnames(seqtab))
keep_len <- seqlens >= 100 & seqlens <= 260
cat("Length-filter [100,260]:", sum(keep_len), "/", length(seqlens), "\n")
seqtab <- seqtab[, keep_len, drop = FALSE]

cat("\n[4/4] removeBimeraDenovo ...\n")
seqtab_nochim <- dada2::removeBimeraDenovo(
  seqtab, method = "consensus", multithread = TRUE, verbose = TRUE
)
cat("After chimera:", dim(seqtab_nochim),
    " | ", round(100 * sum(seqtab_nochim) / sum(seqtab), 1),
    "% reads retained\n")

# Save seqtab
seqtab_path <- file.path(its_dir,
                         paste0("seqtab_nochim_R1only_", date_tag, ".rds"))
saveRDS(seqtab_nochim, seqtab_path)
cat("Saved:", seqtab_path, "\n")

getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(out[keep, ],
               denoised = sapply(dadaFs, getN),
               nonchim  = rowSums(seqtab_nochim))
colnames(track)[1:2] <- c("input", "filtered")
track_df <- as.data.frame(track) |>
  tibble::rownames_to_column("sample_path") |>
  mutate(sample = sample_names)
write_tsv(track_df,
          file.path(its_dir,
                    paste0("track_reads_R1only_", date_tag, ".tsv")))

# ============================================================
# Build phyloseq with metadata
# ============================================================
ps_16s <- readRDS(ps_16s)
md_16s <- as(sample_data(ps_16s), "data.frame")
md_lrt <- md_16s[md_16s$dataset == "LRT", ]
keep_md <- intersect(rownames(md_lrt), rownames(seqtab_nochim))
seqtab_aligned <- seqtab_nochim[keep_md, , drop = FALSE]
md_aligned     <- md_lrt[keep_md, , drop = FALSE]
asv_seqs <- colnames(seqtab_aligned)
asv_ids  <- sprintf("fASV%05d", seq_along(asv_seqs))
colnames(seqtab_aligned) <- asv_ids

ps <- phyloseq(
  otu_table(seqtab_aligned, taxa_are_rows = FALSE),
  sample_data(md_aligned),
  refseq(setNames(DNAStringSet(asv_seqs), asv_ids))
)
print(ps)
ra_idx <- md_aligned$region == "Risopatron" &
          md_aligned$sampler == "Coriolis_mu" & !md_aligned$is_control
cat("Risopatron Air n =", sum(ra_idx),
    "  total reads in those =", sum(seqtab_aligned[ra_idx, ]), "\n")

ps_path <- file.path(its_dir, paste0("ps_LRT_ITS_R1only_", date_tag, ".RDS"))
saveRDS(ps, ps_path)
cat("Saved:", ps_path, "\n")
