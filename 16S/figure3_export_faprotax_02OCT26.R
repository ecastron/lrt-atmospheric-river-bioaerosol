# figure3_export_faprotax_02OCT26.R
# Exports the ASV table with taxonomy paths in the format expected by FAPROTAX collapse_table.py (LRT air, bulk soil
# and rhizosphere samples; controls and Union Glacier excluded). The FAPROTAX command is given at the end of the script.
# Input  : 16S/ps_merged_withTree_noorg_01OCT26.RDS
# Outputs: 16S/figure3_otu_with_taxonomy_02OCT26.tsv, 16S/figure3_lrt_metadata_02OCT26.csv
# Run    : Rscript 16S/figure3_export_faprotax_02OCT26.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(readr)
})

set.seed(100)
base     <- "./16S"
date_tag <- "02OCT26"

ps_path <- file.path(base, "ps_merged_withTree_noorg_01OCT26.RDS")
stopifnot(file.exists(ps_path))
ps <- readRDS(ps_path)

# --- LRT subset (matches Figure 2) ------------------------------------
md <- as(sample_data(ps), "data.frame")
keep <- md$dataset == "LRT" &
        !md$is_control &
        md$region  != "Union Glacier" &
        md$sampler %in% c("Coriolis_mu", "soil", "rhizosphere")
ps_lrt <- prune_samples(keep, ps)
ps_lrt <- prune_taxa(taxa_sums(ps_lrt) > 0, ps_lrt)
cat("LRT working set: ", nsamples(ps_lrt), " samples, ",
    ntaxa(ps_lrt), " ASVs\n", sep = "")

# --- abundance matrix (taxa as rows) ----------------------------------
otu <- as(otu_table(ps_lrt), "matrix")
if (!taxa_are_rows(ps_lrt)) otu <- t(otu)

# --- taxonomy path: key:value per rank, joined by ":" -----------------
# FAPROTAX matches by name within the path, so the rank labels make
# matching robust regardless of column order.
tax <- as.data.frame(tax_table(ps_lrt), stringsAsFactors = FALSE)
get_rank <- function(rank) {
  v <- as.character(tax[[rank]])
  v[is.na(v) | v == "" | v == "NA"] <- "unclassified"
  v
}
taxonomy_path <- paste(
  "Kingdom", get_rank("Kingdom"),
  "Phylum",  get_rank("Phylum"),
  "Class",   get_rank("Class"),
  "Order",   get_rank("Order"),
  "Family",  get_rank("Family"),
  "Genus",   get_rank("Genus"),
  "Species", get_rank("Species"),
  sep = ":"
)
names(taxonomy_path) <- rownames(tax)

# Align taxonomy_path to otu rows
taxonomy_path <- taxonomy_path[rownames(otu)]

# --- combined table ---------------------------------------------------
otu_tax <- data.frame(
  OTUID    = rownames(otu),
  taxonomy = taxonomy_path,
  otu,
  check.names = FALSE
)

otu_path <- file.path(base, paste0("figure3_otu_with_taxonomy_",
                                   date_tag, ".tsv"))
write.table(otu_tax, file = otu_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", otu_path)

# --- metadata for plotting --------------------------------------------
md_out <- as(sample_data(ps_lrt), "data.frame") |>
  mutate(sample_id = sample_names(ps_lrt),
         medium    = case_when(
           sampler == "Coriolis_mu" ~ "Air",
           sampler == "soil"        ~ "Bulk soil",
           sampler == "rhizosphere" ~ "Rhizosphere"
         )) |>
  select(sample_id, dataset, source2, region, sampler, medium,
         latitude, longitude, date, atm_river)

md_path <- file.path(base, paste0("figure3_lrt_metadata_",
                                  date_tag, ".csv"))
write_csv(md_out, md_path)
message("Wrote: ", md_path)

# Next step (run from the project 16S/ directory):
#   python3 FAPROTAX_1.2.12/collapse_table.py \
#     -i figure3_otu_with_taxonomy_<DATE>.tsv \
#     -o figure3_faprotax_<DATE>.tsv \
#     -g FAPROTAX_1.2.12/FAPROTAX.txt \
#     -d taxonomy \
#     -r figure3_faprotax_report_<DATE>.txt \
#     -n columns_after_collapsing -v -f
