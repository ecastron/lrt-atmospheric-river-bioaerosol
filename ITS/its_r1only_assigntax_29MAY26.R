# its_r1only_assigntax_29MAY26.R
# Assigns UNITE taxonomy to the ITS ASVs retained after decontamination.
# Input  : ITS/asvs_r1only_decontam_29MAY26.fasta, UNITE general FASTA release for eukaryotes
# Output : ITS/taxa_r1only_decontam_29MAY26.rds
# Run    : Rscript ITS/its_r1only_assigntax_29MAY26.R

suppressPackageStartupMessages({ library(dada2); library(Biostrings) })
set.seed(100)
base    <- "."
its_dir <- file.path(base, "ITS")
date_tag <- "29MAY26"
unite <- file.path(its_dir, "sh_general_release_all_19.02.2025",
                   "sh_general_release_dynamic_all_19.02.2025.fasta")
stopifnot(file.exists(unite))

seqs <- readDNAStringSet(file.path(its_dir,
          paste0("asvs_r1only_decontam_", date_tag, ".fasta")))
ids  <- names(seqs); seqv <- as.character(seqs)
cat(length(seqv), "ASVs to classify\n")

taxa <- dada2::assignTaxonomy(seqv, refFasta = unite,
                              multithread = TRUE, tryRC = TRUE)
rownames(taxa) <- ids   # key by fASV id (not sequence)
cat("Kingdom Fungi:", sum(taxa[, "Kingdom"] == "k__Fungi", na.rm = TRUE),
    "| Genus-level:", sum(!is.na(taxa[, "Genus"])), "\n")
saveRDS(taxa, file.path(its_dir, paste0("taxa_r1only_decontam_",
                                        date_tag, ".rds")))
cat("DONE taxonomy\n")
