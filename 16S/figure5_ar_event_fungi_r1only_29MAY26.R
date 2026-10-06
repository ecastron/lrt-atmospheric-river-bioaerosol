# figure5_ar_event_fungi_r1only_29MAY26.R
# Fungal (ITS1) counterpart of the bacterial atmospheric-river analysis at Risopatron: weighted UniFrac PCoA trajectory,
# alpha diversity, genus heatmap and phylogenetic structure, with statistics reported in the text.
# Input  : ITS/ps_LRT_ITS_R1only_decontam_withTree_29MAY26.RDS (ITS chain), 16S/palettes.R
# Outputs: 16S/figures/figure5_ar_event_fungi_r1only_29MAY26.{pdf,png} and _stats_29MAY26.txt
# Run    : Rscript 16S/figure5_ar_event_fungi_r1only_29MAY26.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(picante)
  library(vegan)
  library(ape)
})

set.seed(100)
base       <- "."
its_dir    <- file.path(base, "ITS")
date_tag   <- "r1only_29MAY26"
ps_tag     <- "29MAY26"
fig_dir    <- file.path(base, "16S", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

source(file.path(base, "16S", "palettes.R"))

ps_path <- file.path(its_dir, paste0("ps_LRT_ITS_R1only_decontam_withTree_", ps_tag, ".RDS"))
stopifnot(file.exists(ps_path))
ps <- readRDS(ps_path)

log1p10 <- scales::trans_new(
  name      = "log1p10",
  transform = function(x) log10(x + 1),
  inverse   = function(x) 10^x - 1
)

# ============================================================
# 1) Subset Risopatron Air, assign phase
# ============================================================
md <- as(sample_data(ps), "data.frame")
keep <- md$dataset == "LRT" & md$region == "Risopatron" &
        md$sampler == "Coriolis_mu" & !md$is_control
ps_ra <- prune_samples(keep, ps)
ps_ra <- prune_taxa(taxa_sums(ps_ra) > 0, ps_ra)

# After removing the Komagataella reagent yeast + non-Fungi, samples that
# were essentially pure yeast can drop to zero non-yeast fungal reads.
# Weighted UniFrac / PCoA are undefined for all-zero samples, so drop
# them (and report) rather than feed NaN distances to eigen().
zero_reads <- sample_names(ps_ra)[sample_sums(ps_ra) == 0]
if (length(zero_reads)) {
  cat("Dropping zero-read Risopatron-air samples (pure yeast pre-cleaning): ",
      paste(zero_reads, collapse = ", "), "\n", sep = "")
  ps_ra <- prune_samples(sample_sums(ps_ra) > 0, ps_ra)
  ps_ra <- prune_taxa(taxa_sums(ps_ra) > 0, ps_ra)
}

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
stopifnot(nsamples(ps_ra) >= 20, ape::is.rooted(phy_tree(ps_ra)))
cat("Risopatron Air (fungi):", nsamples(ps_ra), "samples,",
    ntaxa(ps_ra), "ASVs (after pruning zeros)\n")

sample_levels <- sd_ra$sample_id
date_labels   <- format(sd_ra$date, "%b %d")
names(date_labels) <- sd_ra$sample_id

ar_rect <- annotate("rect",
                    xmin = ar_affected_start, xmax = ar_affected_end,
                    ymin = -Inf, ymax = Inf,
                    fill = "#C44E52", alpha = 0.10)

# ============================================================
# 2) Panel A: Weighted UniFrac PCoA trajectory
# ============================================================
ps_rel <- transform_sample_counts(ps_ra, function(x) {
  s <- sum(x); if (s == 0) x else x / s
})
D_wu  <- phyloseq::UniFrac(ps_rel, weighted = TRUE, normalized = TRUE)
ord   <- ape::pcoa(D_wu)
ax_pct <- 100 * ord$values$Relative_eig[1:2]

coords <- as.data.frame(ord$vectors[, 1:2]) |>
  tibble::rownames_to_column("sample_id") |>
  rename(Axis1 = Axis.1, Axis2 = Axis.2) |>
  left_join(sd_ra |> select(sample_id, date, phase), by = "sample_id") |>
  arrange(date)

pA <- ggplot(coords, aes(Axis1, Axis2)) +
  geom_path(arrow = arrow(length = unit(2.5, "pt"), type = "closed"),
            color = "gray40", linewidth = 0.4) +
  geom_point(aes(fill = phase), shape = 21, color = "black",
             size = 3, stroke = 0.4) +
  scale_fill_manual(values = pal_ar_phase, name = "Phase") +
  labs(x = sprintf("Axis 1 (%.1f%%)", ax_pct[1]),
       y = sprintf("Axis 2 (%.1f%%)", ax_pct[2]),
       subtitle = "Weighted UniFrac PCoA, time-ordered trajectory (fungi)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle    = element_text(size = 10))

# ============================================================
# 3) Panel B: Alpha diversity time series
# ============================================================
alpha <- phyloseq::estimate_richness(ps_ra,
                                     measures = c("Chao1", "Shannon"))
alpha$sample_id <- rownames(alpha)
alpha_long <- alpha |>
  select(sample_id, Chao1, Shannon) |>
  pivot_longer(c(Chao1, Shannon),
               names_to = "metric", values_to = "value") |>
  left_join(sd_ra |> select(sample_id, date, phase),
            by = "sample_id") |>
  mutate(metric = factor(metric, levels = c("Chao1", "Shannon")))

pB <- ggplot(alpha_long, aes(date, value)) +
  ar_rect +
  geom_line(aes(group = metric), color = "gray50", linewidth = 0.4) +
  geom_point(aes(fill = phase), shape = 21, color = "black",
             size = 2.6, stroke = 0.3) +
  scale_fill_manual(values = pal_ar_phase, name = "Phase") +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             strip.position = "left") +
  labs(x = NULL, y = NULL,
       subtitle = "Fungal alpha diversity (shaded: AR-affected, Feb 5-11)") +
  theme_bw(base_size = 11) +
  theme(strip.placement    = "outside",
        strip.background.y = element_rect(fill = "gray92", color = NA),
        strip.text.y.left  = element_text(angle = 90, face = "bold"),
        panel.grid.minor   = element_blank(),
        plot.subtitle      = element_text(size = 10))

# ============================================================
# 4) Panel C: Genus heatmap (top 20 + Other)
# ============================================================
ps_g <- tax_glom(ps_ra, taxrank = "Genus", NArm = FALSE)
ps_g <- transform_sample_counts(ps_g, function(x) {
  s <- sum(x); if (s == 0) x else x / s
})

genus_tax <- as.data.frame(tax_table(ps_g))
genus_tax$Genus_lbl <- ifelse(is.na(genus_tax$Genus) | genus_tax$Genus == "" |
                              grepl("^g__$", genus_tax$Genus),
                              "Unassigned",
                              gsub("^g__", "", genus_tax$Genus))
genus_sums <- taxa_sums(ps_g)
top_idx    <- order(genus_sums, decreasing = TRUE)[1:min(20, length(genus_sums))]
top_genera <- unique(genus_tax$Genus_lbl[top_idx])

heat_long <- psmelt(ps_g) |>
  mutate(Genus_lbl = ifelse(is.na(Genus) | Genus == "" |
                            grepl("^g__$", Genus),
                            "Unassigned",
                            gsub("^g__", "", Genus)),
         Taxon = ifelse(Genus_lbl %in% top_genera, Genus_lbl, "Other")) |>
  group_by(Sample, Taxon) |>
  summarize(Abundance = sum(Abundance), .groups = "drop") |>
  filter(Sample %in% sample_levels) |>
  mutate(Sample = factor(Sample, levels = sample_levels))

tax_order <- heat_long |>
  group_by(Taxon) |>
  summarize(tot = sum(Abundance), .groups = "drop") |>
  arrange(desc(tot)) |>
  pull(Taxon)
tax_order <- c(setdiff(tax_order, c("Other", "Unassigned")),
               intersect(c("Unassigned", "Other"), tax_order))
heat_long <- heat_long |>
  mutate(Taxon = factor(Taxon, levels = rev(tax_order)))

yes_idx <- which(sd_ra$atm_river == "yes")
ar_box  <- annotate("rect",
                    xmin = min(yes_idx) - 0.5,
                    xmax = max(yes_idx) + 0.5,
                    ymin = 0.5,
                    ymax = length(levels(heat_long$Taxon)) + 0.5,
                    color = "black", fill = NA, linewidth = 0.7)

pC <- ggplot(heat_long, aes(Sample, Taxon, fill = Abundance * 100)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#08306B", "#6BAED6", "#FFFFE5", "#FB6A4A", "#A50F15"),
    trans   = log1p10,
    breaks  = c(0.01, 0.1, 1, 10),
    labels  = function(x) paste0(x, "%"),
    limits  = c(0.01, 100),
    oob     = scales::squish,
    name    = "Relative\nabundance\n(log10+1)"
  ) +
  scale_x_discrete(labels = date_labels) +
  ar_box +
  labs(x = NULL, y = NULL,
       subtitle = "Fungal genus relative abundance (top 20 + Other; box = AR-affected)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y      = element_text(face = "italic", size = 9),
        panel.grid       = element_blank(),
        plot.subtitle    = element_text(size = 10))

# ============================================================
# 5) Panel D: Faith's PD only
# ============================================================
# NTI is not plotted in Panel D: pre-AR fungal samples are
# extremely sparse (4 of 15 Before samples carry only 1 ASV: RA2, RA4,
# RA6, RA11), and picante::ses.mntd requires >=2 taxa to define a
# nearest-neighbor distance. Those samples returned NaN, producing
# visible discontinuities in the NTI trace. Faith's PD is well-defined
# at n=1 (tip-to-root branch length) and the §2.5 statistics for fungal
# NTI were n.s. (p = 0.40), so NTI is left out of from the figure
# and reported in the stats file with a sparse-sample caveat.
otu_mat <- as(otu_table(ps_ra), "matrix")
if (taxa_are_rows(ps_ra)) otu_mat <- t(otu_mat)
otu_pa  <- (otu_mat > 0) + 0L
otu_pa  <- otu_pa[sample_levels, , drop = FALSE]
tr <- phy_tree(ps_ra)

cat("Computing Faith's PD ...\n")
pd_df <- picante::pd(otu_pa, tr, include.root = TRUE) |>
  tibble::rownames_to_column("sample_id") |>
  rename(Faith_PD = PD) |>
  select(sample_id, Faith_PD)

# NTI is still computed for the stats file (with the n>=2 caveat) but
# not plotted.
cat("Computing NTI via ses.mntd (999 perms; reported in stats only) ...\n")
set.seed(100)
co <- cophenetic(tr)
ses <- picante::ses.mntd(otu_pa, co,
                         null.model = "taxa.labels",
                         runs = 999,
                         abundance.weighted = FALSE)
nti_df <- data.frame(sample_id = rownames(ses),
                     NTI       = -1 * ses$mntd.obs.z,
                     n_ASV     = rowSums(otu_pa)[rownames(ses)])

phy_long <- pd_df |>
  inner_join(sd_ra |> select(sample_id, date, phase),
             by = "sample_id") |>
  mutate(metric = factor("Faith_PD",
                         levels = c("Faith_PD"),
                         labels = c("Faith's PD")),
         value  = Faith_PD)

pD <- ggplot(phy_long, aes(date, value)) +
  ar_rect +
  geom_line(aes(group = metric), color = "gray50", linewidth = 0.4) +
  geom_point(aes(fill = phase), shape = 21, color = "black",
             size = 2.6, stroke = 0.3) +
  scale_fill_manual(values = pal_ar_phase, name = "Phase") +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             strip.position = "left") +
  labs(x = NULL, y = NULL,
       subtitle = paste0("Phylogenetic richness (Faith's PD). NTI now ",
                         "well-defined (no singletons after R1-only ",
                         "recovery) but n.s.; reported in stats file.")) +
  theme_bw(base_size = 11) +
  theme(strip.placement    = "outside",
        strip.background.y = element_rect(fill = "gray92", color = NA),
        strip.text.y.left  = element_text(angle = 90, face = "bold"),
        panel.grid.minor   = element_blank(),
        plot.subtitle      = element_text(size = 9))

# ============================================================
# 6) Compose & save
# ============================================================
fig <- ((pA | pB) / pC / pD) +
  plot_layout(heights = c(1.0, 1.4, 0.55), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title    = "Figure 5. Fungal AR signal at Risopatron Air (ITS, R1-only; decontam + Komagataella removed)",
    subtitle = "Strict During (UTC): 2022-02-07/08. Visual shading/box: AR-affected Feb 5-11.") &
  theme(plot.tag        = element_text(face = "bold"),
        legend.position = "bottom")

pdf_path <- file.path(fig_dir, paste0("figure5_ar_event_fungi_", date_tag, ".pdf"))
png_path <- file.path(fig_dir, paste0("figure5_ar_event_fungi_", date_tag, ".png"))
ggsave(pdf_path, fig, width = 14, height = 13)
ggsave(png_path, fig, width = 14, height = 13, dpi = 300)
cat("Saved:\n  ", pdf_path, "\n  ", png_path, "\n", sep = "")

# ============================================================
# 7) Stats
# ============================================================
stats_path <- file.path(fig_dir,
                        paste0("figure5_ar_event_fungi_stats_", date_tag, ".txt"))
sink(stats_path)
cat("Figure 5 stats - Fungal AR signal at Risopatron Air (ITS, R1-only; decontam + Komagataella removed)\n")
cat("Date: ", date_tag, " | input: ", basename(ps_path), "\n", sep = "")
cat("\n=== Sample counts by phase ===\n")
print(table(sd_ra$phase))
cat("atm_river=='yes':", sum(sd_ra$atm_river == "yes"), "\n")

cat("\n=== PERMANOVA: wUniFrac ~ phase (9999 perms) ===\n")
# Align sd to D_wu's label order; sd_ra is chronological (after
# arrange(date)) but D_wu inherits the phyloseq's alphabetical sample
# order, and adonis2 position-matches when those disagree.
sd_perm <- sd_ra[attr(D_wu, "Labels"), ]
set.seed(100)
print(vegan::adonis2(D_wu ~ phase, data = sd_perm,
                     permutations = 9999, by = "margin"))
sd_perm$ar2 <- if_else(sd_perm$phase == "Before",
                       "Before", "During+After")
sd_ar2 <- sd_ra |> mutate(ar2 = if_else(phase == "Before",
                                        "Before", "During+After"))
cat("\n--- Pairwise: Before vs During+After ---\n")
set.seed(100)
print(vegan::adonis2(D_wu ~ ar2, data = sd_perm,
                     permutations = 9999, by = "margin"))

cat("\n=== Wilcoxon Before vs During+After ===\n")
metrics <- alpha |>
  select(sample_id, Chao1, Shannon) |>
  inner_join(pd_df,  by = "sample_id") |>
  inner_join(nti_df, by = "sample_id") |>
  inner_join(sd_ar2 |> select(sample_id, ar2), by = "sample_id")
for (m in c("Chao1", "Shannon", "Faith_PD", "NTI")) {
  v <- metrics[[m]]
  g <- factor(metrics$ar2, levels = c("Before", "During+After"))
  w  <- suppressWarnings(wilcox.test(v ~ g, exact = FALSE))
  cat(sprintf("%-12s  W=%.1f  p=%.4f  (median Before=%.3f, During+After=%.3f)\n",
              m, w$statistic, w$p.value,
              median(v[g == "Before"], na.rm = TRUE),
              median(v[g == "During+After"], na.rm = TRUE)))
}
n_na_nti <- sum(is.na(metrics$NTI))
n_singleton <- sum(metrics$n_ASV == 1, na.rm = TRUE)
cat(sprintf("\nNote: NTI is NA for %d samples (%d with n_ASV == 1; ses.mntd ",
            n_na_nti, n_singleton))
cat("requires >=2 taxa). NTI dropped from Figure 5 Panel D for this reason.\n")

# AR-arrival ASVs (concat)
otu_ra <- as(otu_table(ps_ra), "matrix")
if (taxa_are_rows(ps_ra)) otu_ra <- t(otu_ra)
otu_ra <- otu_ra[sample_levels, , drop = FALSE]
rel_ra <- sweep(otu_ra, 1, rowSums(otu_ra), "/")
rel_ra[is.na(rel_ra)] <- 0
before_idx  <- which(sd_ra$phase == "Before")
durpost_idx <- which(sd_ra$phase != "Before")
mean_before  <- colMeans(rel_ra[before_idx, ,  drop = FALSE])
mean_durpost <- colMeans(rel_ra[durpost_idx, , drop = FALSE])
log2fc       <- log2((mean_durpost + 1e-6) / (mean_before + 1e-6))
prev_after   <- colMeans(otu_ra[durpost_idx, , drop = FALSE] > 0)
ar_arrival_idx <- which(prev_after >= 1/3 & log2fc >= 2)
cat(sprintf("\nAR-arrival fungal ASVs (post >= 4x pre, prev >= 1/3 of post): %d / %d\n",
            length(ar_arrival_idx), ncol(rel_ra)))

# Top arrival ASVs with taxonomy
tt <- as.data.frame(tax_table(ps_ra))
arrival_df <- data.frame(
  asv          = colnames(rel_ra)[ar_arrival_idx],
  rel_durpost  = mean_durpost[ar_arrival_idx] * 100,
  log2fc       = log2fc[ar_arrival_idx]
) |>
  left_join(tt |> tibble::rownames_to_column("asv") |>
              select(asv, Phylum, Class, Order, Family, Genus),
            by = "asv") |>
  arrange(desc(rel_durpost))
cat("\n--- Top-30 AR-arrival fungal ASVs ---\n")
print(head(arrival_df, 30), row.names = FALSE)

# Phylum-level summary
cat("\n--- AR-arrival ASVs by Phylum ---\n")
print(table(arrival_df$Phylum, useNA = "always"))

# Per-genus shifts (top 20 by post abundance)
cat("\n--- Per-genus Wilcoxon Before vs During+After (top genera by post abundance, BH-adj) ---\n")
heat_with_phase <- heat_long |>
  rename(sample_id = Sample) |>
  mutate(sample_id = as.character(sample_id)) |>
  inner_join(sd_ar2 |> select(sample_id, ar2), by = "sample_id")
genus_tab <- heat_with_phase |>
  filter(Taxon != "Other") |>
  group_by(Taxon) |>
  summarize(
    mean_Before = mean(Abundance[ar2 == "Before"]) * 100,
    mean_DurAft = mean(Abundance[ar2 == "During+After"]) * 100,
    p_raw       = suppressWarnings(
      wilcox.test(Abundance ~ ar2, exact = FALSE)$p.value),
    .groups = "drop") |>
  mutate(p_BH = p.adjust(p_raw, method = "BH")) |>
  arrange(p_BH)
print(genus_tab, n = Inf)

sink()
cat("Stats:", stats_path, "\n")
