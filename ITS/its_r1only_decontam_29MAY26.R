# its_r1only_decontam_29MAY26.R
# Removes contaminant ITS ASVs with decontam (prevalence method; negative reference = buffer controls; seawater samples
# are environmental and excluded from the reference). Both thresholds (0.1 and 0.5) are reported; removal uses 0.1.
# Input  : ITS/ps_LRT_ITS_R1only_08MAY26b.RDS, metadata_LRT.csv
# Outputs: ITS/ps_LRT_ITS_R1only_decontam_pretax_29MAY26.RDS, ITS/asvs_r1only_decontam_29MAY26.fasta,
#          ITS/decontam_r1only_contaminants_29MAY26.tsv
# Run    : Rscript ITS/its_r1only_decontam_29MAY26.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(decontam)
  library(Biostrings)
  library(dplyr)
  library(readr)
})

set.seed(100)
base    <- "."
its_dir <- file.path(base, "ITS")
date_tag <- "29MAY26"

ps <- readRDS(file.path(its_dir, "ps_LRT_ITS_R1only_08MAY26b.RDS"))
cat("Loaded R1-only:", nsamples(ps), "samples x", ntaxa(ps), "ASVs\n")

sd <- as(sample_data(ps), "data.frame")

# --- map control type from metadata_LRT.csv (sample_type) --------------
meta_long <- read.csv(file.path(base, "metadata_LRT.csv"),
                      stringsAsFactors = FALSE)
sd$ctrl_type <- meta_long$sample_type[match(rownames(sd), meta_long$sample_name)]
cat("\nctrl_type for is_control samples:\n")
print(table(sd$ctrl_type[sd$is_control], useNA = "ifany"))

# Negative reference = buffer_control only (matches 16S convention)
sd$is_neg_ref <- !is.na(sd$ctrl_type) & sd$ctrl_type == "buffer_control"
sample_data(ps) <- sample_data(sd)
cat("\ndecontam negative reference = buffer_control (n =",
    sum(sd$is_neg_ref), ")\n")
stopifnot(sum(sd$is_neg_ref) >= 2)

# --- decontam prevalence -----------------------------------------------
contam_01 <- isContaminant(ps, method = "prevalence", neg = sd$is_neg_ref,
                           threshold = 0.1, detailed = TRUE)
contam_05 <- isContaminant(ps, method = "prevalence", neg = sd$is_neg_ref,
                           threshold = 0.5, detailed = TRUE)
cat("\nFlagged contaminant ASVs (threshold 0.1, default):   ",
    sum(contam_01$contaminant), "/", ntaxa(ps), "\n")
cat("Flagged contaminant ASVs (threshold 0.5, aggressive): ",
    sum(contam_05$contaminant), "/", ntaxa(ps), "\n")

# Reads removed in Risopatron air specifically (the focal subset)
ra <- sd$region == "Risopatron" & sd$sampler == "Coriolis_mu" & !sd$is_control
otu <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)
ra_reads_tot <- sum(otu[ra, ])
ra_reads_rm  <- sum(otu[ra, contam_01$contaminant])
cat(sprintf("\nRisopatron air: decontam(0.1) removes %d / %d reads (%.1f%%)\n",
            ra_reads_rm, ra_reads_tot, 100 * ra_reads_rm / ra_reads_tot))

# --- write contaminant table -------------------------------------------
contam_tbl <- contam_01 |>
  tibble::rownames_to_column("asv") |>
  transmute(asv, p_strict = p, prev,
            flagged_01 = contaminant) |>
  left_join(contam_05 |> tibble::rownames_to_column("asv") |>
              transmute(asv, flagged_05 = contaminant), by = "asv") |>
  arrange(p_strict)
write_tsv(contam_tbl,
          file.path(its_dir, paste0("decontam_r1only_contaminants_",
                                    date_tag, ".tsv")))

# --- remove decontam(0.1)-flagged, save pretax object ------------------
keep_asv <- !contam_01$contaminant
ps_clean <- prune_taxa(keep_asv, ps)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
cat("\nAfter decontam(0.1) removal:", ntaxa(ps_clean), "ASVs retained\n")

saveRDS(ps_clean,
        file.path(its_dir, paste0("ps_LRT_ITS_R1only_decontam_pretax_",
                                  date_tag, ".RDS")))

# survivor FASTA for taxonomy + tree
writeXStringSet(refseq(ps_clean),
                file.path(its_dir, paste0("asvs_r1only_decontam_",
                                          date_tag, ".fasta")))
cat("Wrote survivor FASTA (", ntaxa(ps_clean), "seqs) + pretax object.\n",
    sep = "")
