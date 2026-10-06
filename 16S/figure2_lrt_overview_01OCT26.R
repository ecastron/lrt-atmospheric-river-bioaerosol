# figure2_lrt_overview_01OCT26.R
# Figure 2. Phylum-level composition of air, bulk soil and rhizosphere samples per site (a) and Faith's
# phylogenetic diversity per site and medium with pairwise Wilcoxon tests (b). Union Glacier samples are excluded.
# Input  : 16S/ps_merged_withTree_noorg_01OCT26.RDS, 16S/palettes.R
# Outputs: 16S/figures/figure2_lrt_overview_01OCT26.{pdf,png}, 16S/figures/figure2_lrt_overview_stats_01OCT26.txt
# Run    : Rscript 16S/figure2_lrt_overview_01OCT26.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggh4x)
  library(patchwork)
  library(picante)
  library(rstatix)
  library(multcompView)
  library(scales)
})

set.seed(100)
base     <- "./16S"
in_tag   <- "noorg_01OCT26"   # input phyloseq object (organelle-free)
date_tag <- "01OCT26"   # output tag
fig_dir  <- file.path(base, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

source(file.path(base, "palettes.R"))

# Display labels for site strips (colors unchanged, from palettes.R).
site_labels_fig <- c(
  "PA_POR"          = "PA/POR",
  "PN_SG"           = "PN/SG",
  "Puerto Williams" = "PW",
  "Risopatron"      = "Risopatrón",
  "Yelcho"          = "Yelcho"
)

ps_path <- file.path(base, paste0("ps_merged_withTree_", in_tag, ".RDS"))
stopifnot(file.exists(ps_path))
ps <- readRDS(ps_path)

# --- 1) restrict to LRT, drop controls, drop Union Glacier ------------
md <- as(sample_data(ps), "data.frame")
keep <- md$dataset == "LRT" &
        !md$is_control &
        md$region != "Union Glacier" &
        md$sampler %in% c("Coriolis_mu", "soil", "rhizosphere")
ps_lrt <- prune_samples(keep, ps)
ps_lrt <- prune_samples(sample_sums(ps_lrt) >= 1000, ps_lrt)
ps_lrt <- prune_taxa(taxa_sums(ps_lrt) > 0, ps_lrt)
cat("LRT working set: ", nsamples(ps_lrt), " samples, ",
    ntaxa(ps_lrt), " ASVs\n", sep = "")

# Add medium (3-way) and ordered factors directly to sample_data so
# psmelt carries them through.
sd <- as(sample_data(ps_lrt), "data.frame") |>
  mutate(
    medium = case_when(
      sampler == "Coriolis_mu" ~ "Air",
      sampler == "soil"        ~ "Bulk soil",
      sampler == "rhizosphere" ~ "Rhizosphere",
      TRUE                     ~ NA_character_
    ),
    medium = factor(medium, levels = c("Air", "Bulk soil", "Rhizosphere")),
    region = factor(region, levels = site_order)
  )
sample_data(ps_lrt) <- sample_data(sd)

# ============================================================
# Panel A: phylum-level stacked bars
# ============================================================
ps_phy <- tax_glom(ps_lrt, taxrank = "Phylum", NArm = FALSE)
ps_phy <- transform_sample_counts(ps_phy, function(x) x / sum(x))

phy_df <- psmelt(ps_phy) |>
  mutate(Phylum = if_else(is.na(Phylum) | Phylum == "", "Unassigned",
                          as.character(Phylum)))

# Top 8 phyla by mean relative abundance, lump rest as "Other"
top_phyla <- phy_df |>
  group_by(Phylum) |>
  summarize(mean_ab = mean(Abundance), .groups = "drop") |>
  arrange(desc(mean_ab)) |>
  slice_head(n = 8) |>
  pull(Phylum)

# Order legend so well-known phyla in pal_phylum come first, then any
# extras (rare in the LRT data) at the end before "Other".
fill_levels <- c(intersect(names(pal_phylum), top_phyla),
                 setdiff(top_phyla, names(pal_phylum)),
                 "Other")
fill_colors <- ifelse(fill_levels %in% names(pal_phylum),
                      pal_phylum[fill_levels], "#999999")
names(fill_colors) <- fill_levels

phy_df <- phy_df |>
  mutate(Phylum_grp = if_else(Phylum %in% top_phyla, Phylum, "Other"))

# Within each (medium, region) panel, order samples by abs(latitude)
# south → north so the in-panel ordering is geographically meaningful.
sample_order <- as(sample_data(ps_lrt), "data.frame") |>
  mutate(sample_id = sample_names(ps_lrt),
         abs_lat   = abs(latitude)) |>
  arrange(medium, region, desc(abs_lat), sample_id) |>
  pull(sample_id)

phy_df <- phy_df |>
  mutate(
    Sample     = factor(Sample, levels = sample_order),
    Phylum_grp = factor(Phylum_grp, levels = fill_levels)
  )

pA <- ggplot(phy_df,
             aes(Sample, Abundance, fill = Phylum_grp)) +
  geom_col(width = 1) +
  scale_fill_manual(values = fill_colors, name = "Phylum") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  facet_nested(. ~ medium + region,
               scales = "free_x", space = "free_x",
               nest_line = element_line(linetype = 1),
               labeller = labeller(region = site_labels_fig),
               # site strips vertical so narrow panels (PW, n = 3) stay legible
               strip = strip_nested(
                 text_x = elem_list_text(angle = c(0, 90), size = c(7, 6)),
                 by_layer_x = TRUE)) +
  labs(x = NULL, y = "Relative abundance") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid         = element_blank(),
    panel.spacing.x    = unit(0.15, "lines"),
    strip.background   = element_rect(fill = "gray92", color = NA),
    strip.text         = element_text(size = 7),
    legend.position    = "right",
    legend.key.size    = unit(0.3, "cm")
  )

# ============================================================
# Panel B: Faith's PD per site, boxplots by medium with letters
# ============================================================
otu_mat <- as(otu_table(ps_lrt), "matrix")
if (taxa_are_rows(ps_lrt)) otu_mat <- t(otu_mat)
tr <- phy_tree(ps_lrt)

faith <- picante::pd(otu_mat, tr, include.root = TRUE)
faith$sample_id <- rownames(faith)

div_df <- as(sample_data(ps_lrt), "data.frame") |>
  mutate(sample_id = sample_names(ps_lrt)) |>
  left_join(faith, by = "sample_id")

# Per-site pairwise Wilcoxon → letters. Sites with only one medium
# (none here, but defensively) get "a".
make_letters <- function(sub) {
  # Drop factor levels with no observations in this site so
  # pairwise_wilcox_test doesn't try to test empty pairs.
  sub <- sub |> mutate(medium = droplevels(factor(medium)))
  mediums <- levels(sub$medium)
  if (length(mediums) < 2) {
    return(tibble(medium = mediums, cld = "a"))
  }
  pw <- rstatix::pairwise_wilcox_test(sub, PD ~ medium,
                                      p.adjust.method = "BH")
  pmat <- with(pw, setNames(p.adj, paste(group1, group2, sep = "-")))
  cld  <- multcompView::multcompLetters(pmat)$Letters
  tibble(medium = names(cld), cld = unname(cld))
}

letters_df <- div_df |>
  group_by(region) |>
  group_modify(~ make_letters(.x)) |>
  ungroup() |>
  mutate(medium = factor(medium,
                         levels = levels(div_df$medium)),
         region = factor(region, levels = site_order))

ypos <- div_df |>
  group_by(region, medium) |>
  summarize(y = max(PD, na.rm = TRUE) * 1.05, .groups = "drop") |>
  left_join(letters_df, by = c("region", "medium"))

pB <- ggplot(div_df, aes(medium, PD, fill = medium)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.65) +
  geom_jitter(width = 0.18, height = 0, size = 0.7,
              alpha = 0.55, shape = 21, color = "black") +
  geom_text(data = ypos,
            aes(medium, y, label = cld),
            inherit.aes = FALSE, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = pal_medium3, guide = "none") +
  facet_grid(. ~ region, scales = "free_x", space = "free_x",
             labeller = labeller(region = site_labels_fig)) +
  labs(x = NULL, y = "Faith's phylogenetic diversity") +
  theme_bw(base_size = 8) +
  theme(axis.text.x      = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray92", color = NA),
        strip.text       = element_text(size = 7))

# ============================================================
# Compose & save
# ============================================================
fig <- (pA / pB) +
  plot_layout(heights = c(1.0, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))

pdf_path <- file.path(fig_dir, paste0("figure2_lrt_overview_", date_tag, ".pdf"))
png_path <- file.path(fig_dir, paste0("figure2_lrt_overview_", date_tag, ".png"))
ggsave(pdf_path, fig, width = 183, height = 150, units = "mm")
ggsave(png_path, fig, width = 183, height = 150, units = "mm", dpi = 300)
cat("Saved:\n  ", pdf_path, "\n  ", png_path, "\n", sep = "")

# ============================================================
# Stats output
# ============================================================
stats_path <- file.path(fig_dir,
                        paste0("figure2_lrt_overview_stats_", date_tag, ".txt"))
sink(stats_path)
cat("Figure 2 stats - LRT composition + Faith's PD\n")
cat("Date: ", date_tag, "\n", sep = "")
cat("Input: ", basename(ps_path), "\n\n", sep = "")

cat("Sample counts (medium x region):\n")
print(table(div_df$medium, div_df$region))

cat("\nFaith's PD summary by region x medium:\n")
print(div_df |>
        group_by(region, medium) |>
        summarize(n = n(),
                  median_PD = median(PD, na.rm = TRUE),
                  mean_PD   = mean(PD,   na.rm = TRUE),
                  sd_PD     = sd(PD,     na.rm = TRUE),
                  .groups = "drop"))

# For the per-site tests, convert medium to character so empty factor
# levels don't trip pairwise_wilcox_test (e.g. PA_POR has no Rhizosphere).
div_df_test <- div_df |> mutate(medium = as.character(medium))

cat("\nKruskal-Wallis (PD ~ medium, within each site):\n")
print(div_df_test |>
        group_by(region) |>
        rstatix::kruskal_test(PD ~ medium))

cat("\nPairwise Wilcoxon within site (BH-adjusted):\n")
print(div_df_test |>
        group_by(region) |>
        rstatix::pairwise_wilcox_test(PD ~ medium,
                                      p.adjust.method = "BH"))

cat("\nCompact letter display (per site):\n")
print(letters_df)

cat("\nTop phyla used in Panel A (by mean relative abundance):\n")
print(top_phyla)
sink()
cat("Stats: ", stats_path, "\n", sep = "")
