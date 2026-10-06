# ar_patagonia_attribution_v3_01OCT26.R
# Genus-level (best available taxonomy) attribution of the taxa that increase at Risopatron during and after the
# atmospheric river: prevalence and abundance of each arrival group in Patagonian air and soil, Southern Ocean air
# (Malard et al. 2022) and Risopatron soil.
# Input  : 16S/ps_merged_withTree_noorg_01OCT26.RDS, 16S/palettes.R
# Outputs: 16S/figures/ar_attribution_v3_01OCT26.{txt,png,pdf}, 16S/figures/ar_attribution_v3_summary_01OCT26.tsv
# Run    : Rscript 16S/ar_patagonia_attribution_v3_01OCT26.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

set.seed(100)
base       <- "."
date_tag   <- "01OCT26"
ps_tag     <- "noorg_01OCT26"   # organelle-free object
fig_dir    <- file.path(base, "16S", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

source(file.path(base, "16S", "palettes.R"))

ps <- readRDS(file.path(base, "16S",
                        paste0("ps_merged_withTree_", ps_tag, ".RDS")))
md <- as(sample_data(ps), "data.frame")
cat("Loaded phyloseq:", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")

# ============================================================
# 1) Build best_tax column (deepest non-NA rank per ASV)
# ============================================================
tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
ranks <- intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                   colnames(tax))

clean_taxon <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA" |
    grepl("^uncultured", x, ignore.case = TRUE)] <- NA
  x
}

best_rank   <- character(nrow(tax))
best_taxon  <- character(nrow(tax))
for (i in seq_len(nrow(tax))) {
  found <- FALSE
  for (r in rev(ranks)) {
    v <- clean_taxon(tax[i, r])
    if (!is.na(v)) {
      best_rank[i]  <- r
      best_taxon[i] <- v
      found <- TRUE
      break
    }
  }
  if (!found) { best_rank[i] <- "Unassigned"; best_taxon[i] <- "Unassigned" }
}
best_tax_label <- ifelse(best_rank == "Unassigned",
                         "Unassigned",
                         paste0(best_rank, ":", best_taxon))

tt2 <- cbind(as.matrix(tax_table(ps)),
             best_rank = best_rank,
             best_tax  = best_tax_label)
tax_table(ps) <- tax_table(tt2)

# ============================================================
# 2) Agglomerate at best_tax level
# ============================================================
cat("\nAgglomerating at best_tax level ...\n")
ps_glom <- tax_glom(ps, taxrank = "best_tax", NArm = FALSE)
cat("Glommed:", ntaxa(ps_glom), "best_tax taxa\n")

# ============================================================
# 3) AR-arrival best_tax taxa at Risopatron Air
# ============================================================
md_glom <- as(sample_data(ps_glom), "data.frame")
ra_keep <- md_glom$dataset == "LRT" & md_glom$region == "Risopatron" &
           md_glom$sampler == "Coriolis_mu" & !md_glom$is_control
ps_ra <- prune_samples(ra_keep, ps_glom)
ps_ra <- prune_taxa(taxa_sums(ps_ra) > 0, ps_ra)

sd_ra <- as(sample_data(ps_ra), "data.frame") |>
  mutate(date  = as.Date(date),
         phase = case_when(
           date <  ar_event_start ~ "Before",
           date <= ar_event_end   ~ "During",
           TRUE                   ~ "After"),
         phase = factor(phase, c("Before","During","After"))) |>
  arrange(date)
ps_ra <- prune_samples(rownames(sd_ra), ps_ra)
sample_data(ps_ra) <- sample_data(sd_ra)

otu_ra <- as(otu_table(ps_ra), "matrix")
if (taxa_are_rows(ps_ra)) otu_ra <- t(otu_ra)
otu_ra <- otu_ra[sd_ra$sample_id, , drop = FALSE]
rel_ra <- sweep(otu_ra, 1, rowSums(otu_ra), "/")
rel_ra[is.na(rel_ra)] <- 0

before_idx  <- which(sd_ra$phase == "Before")
durpost_idx <- which(sd_ra$phase != "Before")
mean_before  <- colMeans(rel_ra[before_idx, ,  drop = FALSE])
mean_durpost <- colMeans(rel_ra[durpost_idx, , drop = FALSE])
log2fc       <- log2((mean_durpost + 1e-6) / (mean_before + 1e-6))
prev_after   <- colMeans(otu_ra[durpost_idx, , drop = FALSE] > 0)

ar_arrival <- which(prev_after >= 1/3 & log2fc >= 2)
cat("\nAR-arrival best_tax taxa: ", length(ar_arrival), " / ",
    ncol(rel_ra), "\n", sep = "")
ar_taxa_ids <- colnames(rel_ra)[ar_arrival]

# ============================================================
# 4) Source groups: BROADER Patagonia (all 3 aggregated regions)
#    + Malard SO Air, SO over land + RA soil control
# ============================================================
patagonian_regions <- c("PA_POR", "PN_SG", "Puerto Williams")

source_groups <- list(
  "RA Before"          = rownames(md_glom)[ra_keep & md_glom$atm_river == "no" &
                                            as.Date(md_glom$date) < ar_event_start],
  "Patagonia air"      = rownames(md_glom)[md_glom$dataset == "LRT" &
                                            md_glom$region %in% patagonian_regions &
                                            md_glom$sampler == "Coriolis_mu" &
                                            !md_glom$is_control],
  "Patagonia soil/rhizo" = rownames(md_glom)[md_glom$dataset == "LRT" &
                                              md_glom$region %in% patagonian_regions &
                                              md_glom$sampler %in% c("soil","rhizosphere") &
                                              !md_glom$is_control],
  "Malard SO Air"      = rownames(md_glom)[md_glom$dataset == "Malard2022" &
                                            md_glom$source2 == "SO Air"],
  "Malard SO over land" = rownames(md_glom)[md_glom$dataset == "Malard2022" &
                                              md_glom$source2 == "SO Air over land"],
  "RA soil (re-aero)"  = rownames(md_glom)[md_glom$dataset == "LRT" &
                                            md_glom$region == "Risopatron" &
                                            md_glom$sampler %in% c("soil","rhizosphere") &
                                            !md_glom$is_control]
)
cat("\nSource group sample counts:\n")
print(sapply(source_groups, length))

otu_full <- as(otu_table(ps_glom), "matrix")
if (taxa_are_rows(ps_glom)) otu_full <- t(otu_full)

per_source_summary <- function(asv, source_ids) {
  if (length(source_ids) == 0) return(c(prev = NA, mean_ra = NA))
  m_counts <- otu_full[source_ids, asv, drop = FALSE]
  m_full   <- otu_full[source_ids, , drop = FALSE]
  rs       <- rowSums(m_full)
  rs[rs == 0] <- 1
  m_rel    <- m_counts / rs
  c(prev = mean(m_counts > 0), mean_ra = mean(m_rel) * 100)
}

attribution <- lapply(ar_taxa_ids, function(asv) {
  rows <- lapply(names(source_groups), function(sg) {
    s <- per_source_summary(asv, source_groups[[sg]])
    data.frame(asv = asv, source = sg,
               prev = unname(s["prev"]),
               mean_ra = unname(s["mean_ra"]))
  }) |> bind_rows()
  rows
}) |> bind_rows()

tax_glom_df <- as.data.frame(tax_table(ps_glom), stringsAsFactors = FALSE) |>
  tibble::rownames_to_column("asv")

attr_full <- attribution |>
  left_join(data.frame(asv = ar_taxa_ids,
                       rel_durpost = mean_durpost[ar_taxa_ids] * 100,
                       log2fc      = log2fc[ar_taxa_ids]),
            by = "asv") |>
  left_join(tax_glom_df |> select(asv, best_rank, best_tax,
                                  Phylum, Class, Family, Genus),
            by = "asv")

# ============================================================
# 5) Attribution heuristic (broadened Patagonia, Malard contender)
# ============================================================
attr_wide <- attribution |>
  select(asv, source, prev) |>
  pivot_wider(names_from = source, values_from = prev) |>
  left_join(data.frame(asv = ar_taxa_ids,
                       rel_durpost = mean_durpost[ar_taxa_ids] * 100,
                       log2fc      = log2fc[ar_taxa_ids]),
            by = "asv") |>
  left_join(tax_glom_df |> select(asv, best_rank, best_tax,
                                  Phylum, Class, Family, Genus),
            by = "asv") |>
  arrange(desc(rel_durpost))

classify_source <- function(row) {
  pat   <- max(as.numeric(row[c("Patagonia air","Patagonia soil/rhizo")]),
               na.rm = TRUE)
  mar   <- max(as.numeric(row[c("Malard SO Air","Malard SO over land")]),
               na.rm = TRUE)
  loc   <- as.numeric(row["RA soil (re-aero)"])
  vals  <- c(Patagonian = pat, Marine = mar, "Local re-aerosol" = loc)
  if (max(vals, na.rm = TRUE) < 0.25) return("Unattributed")
  names(which.max(vals))
}
attr_wide$attribution <- apply(attr_wide, 1, classify_source)

cat("\n=== Attribution distribution (top 30 by post abundance) ===\n")
top30_w <- attr_wide |> slice_head(n = 30)
print(table(top30_w$attribution))

cat("\n=== Attribution distribution (all", nrow(attr_wide), "AR-arrival taxa) ===\n")
print(table(attr_wide$attribution))

# ============================================================
# 6) Heatmap of top 30
# ============================================================
top30_ids <- top30_w$asv
plotdat <- attr_full |>
  filter(asv %in% top30_ids) |>
  mutate(label = paste0(best_tax, " (", round(rel_durpost, 2), "%)"),
         label = factor(label,
                        levels = rev(unique(label[order(-rel_durpost)]))),
         source = factor(source,
                         levels = c("RA Before","Patagonia air","Patagonia soil/rhizo",
                                    "Malard SO Air","Malard SO over land",
                                    "RA soil (re-aero)")))

p <- ggplot(plotdat, aes(source, label, fill = prev)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(prev) | prev == 0, "",
                               sprintf("%.0f%%", prev * 100))),
            size = 2.6, color = "gray15") +
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603",
                      limits = c(0, 1),
                      labels = scales::percent_format(),
                      name = "Prevalence\nin source") +
  labs(x = NULL, y = NULL,
       title = "AR-arrival taxa at Risopatron Air: source attribution (v3, broad Patagonia)",
       subtitle = sprintf(
         "Top 30 of %d AR-arrival best-taxonomy taxa (post >= 4x pre, prev >= 1/3)",
         length(ar_arrival))) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid  = element_blank())

png_path <- file.path(fig_dir, paste0("ar_attribution_v3_", date_tag, ".png"))
pdf_path <- file.path(fig_dir, paste0("ar_attribution_v3_", date_tag, ".pdf"))
ggsave(png_path, p, width = 9, height = 9, dpi = 300)
ggsave(pdf_path, p, width = 9, height = 9)

# ============================================================
# 7) Text + tsv outputs
# ============================================================
out_path <- file.path(fig_dir, paste0("ar_attribution_v3_", date_tag, ".txt"))
sink(out_path)
cat("AR Patagonia attribution v3 - best-taxonomy level, broad Patagonia source\n")
cat("Date:", date_tag, "| input ps:", ps_tag, "\n\n")

cat("--- Source-group sample counts ---\n")
print(sapply(source_groups, length))
cat("\n(Patagonian aggregate: ",
    length(source_groups[["Patagonia air"]]) +
    length(source_groups[["Patagonia soil/rhizo"]]),
    " samples; Malard aggregate: ",
    length(source_groups[["Malard SO Air"]]) +
    length(source_groups[["Malard SO over land"]]),
    ")\n", sep = "")

cat("\n--- best_rank distribution (per ASV, pre-glom) ---\n")
print(table(best_rank))
cat("\nGlom level: best_tax: total taxa post-glom:", ntaxa(ps_glom), "\n")
cat("AR-arrival best_tax taxa (post>=4x pre, prev>=1/3 post):",
    length(ar_arrival), "\n")

cat("\n=== Attribution distribution (all", nrow(attr_wide), "AR-arrival taxa) ===\n")
print(table(attr_wide$attribution))
cat("\n=== Attribution distribution (top 30 by post abundance) ===\n")
print(table(top30_w$attribution))

cat("\n--- Top 30 AR-arrival taxa with source prevalences (%) ---\n")
print(top30_w |>
        mutate(across(c("RA Before","Patagonia air","Patagonia soil/rhizo",
                        "Malard SO Air","Malard SO over land",
                        "RA soil (re-aero)"),
                      ~ round(.x * 100, 1))) |>
        select(best_tax, rel_durpost, log2fc, attribution,
               `RA Before`,`Patagonia air`,`Patagonia soil/rhizo`,
               `Malard SO Air`,`Malard SO over land`,`RA soil (re-aero)`,
               Phylum, Class, Family, Genus),
      n = Inf, width = Inf)
sink()

write_tsv(attr_wide,
          file.path(fig_dir,
                    paste0("ar_attribution_v3_summary_", date_tag, ".tsv")))

cat("\nDone.\n")
cat("  ", out_path, "\n")
cat("  ", png_path, "\n")
