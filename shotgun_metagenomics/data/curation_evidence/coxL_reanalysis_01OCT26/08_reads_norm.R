# 08_reads_norm.R
# Read-level form I and form II CoxL per library, per million reads and per recA (copies per genome =
# (reads / median reference length in aa) / (recA reads / 353 aa)), for the QC-passing air libraries and Byers soil.
# recA counts come from the marker recruitment of 26_air_reads_func.sh and 28_byers_reads_func.sh on the same reads.
# Inputs : reads_out/*counts.tsv (04_reads_formI_formII.sh), ref_coxL.tsv, byers_bulksoil.tsv
# Output : reads_formI_formII_per_library.tsv
# Run    : Rscript 08_reads_norm.R (paths are set relative to the repository root)

suppressPackageStartupMessages({library(dplyr); library(readr); library(purrr); library(tidyr)})
setwd("./shotgun_metagenomics/data/curation_evidence/coxL_reanalysis_01OCT26")
loc <- "../../air_func_local"
ref <- read_tsv("ref_coxL.tsv", show_col_types = FALSE)
lenI <- median(ref$length[ref$form == "I"]); lenII <- median(ref$length[ref$form == "II"]); lenRecA <- 353
cnt <- list.files("reads_out", "counts.tsv$", full.names = TRUE) |> map_dfr(read_tsv, show_col_types = FALSE) |>
  mutate(set = ifelse(startsWith(sample, "air_"), "air", "Byers soil"),
         library = ifelse(set == "air", sub("^air_[0-9]+_[0-9]+_([A-Z]+[0-9]+)_.*$", "\\1", sample), sub("^byers_", "", sample)))
# Byers runs (SRR) -> sample ids (ASB) used by the marker recruitment files
runs <- read_tsv("byers_bulksoil.tsv", col_names = FALSE, show_col_types = FALSE)
cnt <- cnt |> mutate(library = ifelse(set == "Byers soil", runs$X2[match(library, runs$X1)], library))
rec <- function(f) { l <- readLines(f); h <- strsplit(l[1], "\t")[[1]]
  tibble(id = h[2], recA_nreads = as.numeric(h[length(h)]),
         recA = as.numeric(sub("^recA\t", "", grep("^recA\t", l, value = TRUE)) %||% NA)) }
recA <- bind_rows(list.files(file.path(loc, "reads_recruit"), "markercounts.tsv$", full.names = TRUE) |> map_dfr(rec) |>
                    mutate(library = sub("^[0-9]+_[0-9]+_([A-Z]+[0-9]+)_.*$", "\\1", id)),
                  list.files(file.path(loc, "byers_recruit"), "markercounts.tsv$", full.names = TRUE) |> map_dfr(rec) |>
                    mutate(library = id))
qc_pass <- c("PWA2", "RA3", "RA8", "RA14", "RA15", "RA16", "RA17", "RA18", "RA22", "RA23", "RA24")
tab <- cnt |> left_join(recA |> select(library, recA, recA_nreads), by = "library") |>
  mutate(qc_pass = set == "Byers soil" | library %in% qc_pass,
         nreads_match = abs(nreads - recA_nreads) / nreads < 0.01,
         formI_per_M = 1e6 * formI / nreads, formII_per_M = 1e6 * formII / nreads,
         formI_per_genome = (formI / lenI) / (recA / lenRecA), formII_per_genome = (formII / lenII) / (recA / lenRecA),
         frac_formI = formI / pmax(formI + formII, 1)) |>
  arrange(set, library)
write_tsv(tab, "reads_formI_formII_per_library.tsv")
cat("median ref length form I", lenI, "| form II", lenII, "\n")
print(as.data.frame(tab |> select(set, library, qc_pass, nreads, recA, nreads_match, formI, formII, outgroup, formI_per_genome, formII_per_genome, frac_formI) |>
  mutate(across(c(formI_per_genome, formII_per_genome, frac_formI), ~ signif(.x, 2)))))
s <- tab |> filter(qc_pass) |> group_by(set) |>
  summarize(n = n(), formI_reads = sum(formI), formII_reads = sum(formII),
            libs_with_formI = sum(formI > 0), median_formI_pg = median(formI_per_genome, na.rm = TRUE), median_formII_pg = median(formII_per_genome, na.rm = TRUE), n_recA = sum(!is.na(recA)),
            pooled_frac_formI = sum(formI) / sum(formI + formII), .groups = "drop")
print(as.data.frame(s))
w <- function(v) { a <- tab |> filter(qc_pass, set == "air") |> pull(v); b <- tab |> filter(set == "Byers soil") |> pull(v)
  cat(v, ": Wilcoxon air vs soil p =", signif(wilcox.test(a, b, exact = FALSE)$p.value, 2), "\n") }
w("formI_per_genome"); w("formII_per_genome"); w("frac_formI")
