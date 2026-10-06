# ACE_dada2.R
# DADA2 processing of the Antarctic Circumnavigation Expedition (ACE) air 16S V4 reads (Malard et al. 2022) with the
# same parameters as the LRT reads, so that the two sequence tables can be merged at the exact-sequence level.
# Error models are learned per sequencing run, so this dataset is processed separately; chimera removal is done
# after merging (merge_LRT_ACE.R) so that it sees the union of variants from both datasets.
# Both datasets start at the post-primer position (R1 begins TACG after 515F, R2 begins CCTG after 806R);
# reads are truncated to 150 + 150 bp as for LRT.
# Inputs : 16S/ACEDATA/Allfastq/*.fastq (BioProject PRJNA697829, uncompressed, paired)
# Outputs: 16S/seqtab_ACE_<DATE>.rds (sequence table, chimeras retained), 16S/track_ACE_<DATE>.csv
# Run    : Rscript 16S/ACE_dada2.R

suppressPackageStartupMessages({
  library(dada2)
})

set.seed(100)  # project convention

# --- paths -------------------------------------------------------------
base       <- "./16S"
ace_path   <- file.path(base, "ACEDATA", "Allfastq")
filt_dir   <- file.path(ace_path, "filtered")
out_dir    <- base
date_tag   <- "07MAY26"   # tag of the outputs used in the manuscript; downstream scripts read *_07MAY26 files

# --- find paired reads -------------------------------------------------
fnFs <- sort(list.files(ace_path, pattern = "_L001_R1_001\\.fastq$", full.names = TRUE))
fnRs <- sort(list.files(ace_path, pattern = "_L001_R2_001\\.fastq$", full.names = TRUE))
stopifnot(length(fnFs) == length(fnRs), length(fnFs) > 0)

# Sample name = leading token before first underscore (matches LRT)
sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Optional QC plots: uncomment to inspect before committing to truncLen
# plotQualityProfile(fnFs[1:4])
# plotQualityProfile(fnRs[1:4])

# --- filter & trim (params identical to LRT.R) ---------------------
if (!dir.exists(filt_dir)) dir.create(filt_dir)
filtFs <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen    = c(150, 150),   # MUST match LRT for mergeable ASVs
  maxN        = 0,
  maxEE       = c(2, 2),
  truncQ      = 2,
  rm.phix     = TRUE,
  compress    = TRUE,
  multithread = TRUE
)

# Drop samples whose filtered files are missing (zero reads passed filter)
keep         <- file.exists(filtFs) & file.exists(filtRs)
filtFs       <- filtFs[keep]
filtRs       <- filtRs[keep]
sample_names <- sample_names[keep]

# --- error model & DADA2 inference -------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = TRUE)

mergers    <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
seqtab_ace <- makeSequenceTable(mergers)

cat("ACE seqtab dims: ", dim(seqtab_ace)[1], " samples x ",
    dim(seqtab_ace)[2], " ASVs\n", sep = "")
cat("ASV length range: ",
    paste(range(nchar(getSequences(seqtab_ace))), collapse = "-"), "\n", sep = "")

# --- save (no chimera removal: done after merge with LRT) -------------
seqtab_path <- file.path(out_dir, paste0("seqtab_ACE_", date_tag, ".rds"))
saveRDS(seqtab_ace, seqtab_path)
message("Saved ", seqtab_path)

# --- read tracking -----------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- data.frame(
  sample    = sample_names,
  input     = out[keep, "reads.in"],
  filtered  = out[keep, "reads.out"],
  denoisedF = sapply(dadaFs,  getN),
  denoisedR = sapply(dadaRs,  getN),
  merged    = sapply(mergers, getN),
  row.names = NULL
)
track_path <- file.path(out_dir, paste0("track_ACE_", date_tag, ".csv"))
write.csv(track, track_path, row.names = FALSE)
message("Saved ", track_path)
