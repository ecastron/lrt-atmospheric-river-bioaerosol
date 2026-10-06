# figureS8_decontam_checks_01OCT26.R
# Supplementary Figure 8. Contamination checks: (a) Bray-Curtis PCoA of the negative controls (buffer, Coriolis and
# deposition controls) with PERMANOVA; (b) betaNTI between Risopatron air sample pairs before and after removing
# decontam-flagged ASVs.
# Inputs : 16S/ps_16S_LRT_withSeqs.RDS, 16S/metadata_LRT3.csv, metadata_LRT.csv,
#          16S/figures/decontam_betaNTI_01OCT26.tsv (decontam_betaNTI_01OCT26.R)
# Outputs: 16S/figures/figureS8_decontam_checks_01OCT26.{pdf,png}
# Run    : Rscript 16S/figureS8_decontam_checks_01OCT26.R

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(readr); library(ggplot2)
  library(vegan); library(patchwork); library(ggrepel)
})

setwd("./16S")
source("palettes.R")
set.seed(100)
out_tag <- "01OCT26"
fig_dir <- "figures"

## ---- panel a: negative-control PCoA (as in decontam_controls_11MAY26.R) ---
ps   <- readRDS("ps_16S_LRT_withSeqs.RDS")
meta <- read.csv("metadata_LRT3.csv", sep = ";", header = TRUE,
                 row.names = 1, check.names = FALSE)
sample_data(ps) <- sample_data(meta)
sd <- as.data.frame(sample_data(ps))
ctrl_re <- "(?i)control|blank"
sd$is_control <- grepl(ctrl_re, sd$Source2) | grepl(ctrl_re, sd$Category2) |
                 grepl(ctrl_re, sd$Category)
meta_long <- read.csv("../metadata_LRT.csv", header = TRUE, stringsAsFactors = FALSE)
sd$ctrl_type <- meta_long$sample_type[match(rownames(sd), meta_long$sample_name)]
sd$ctrl_type[!sd$is_control] <- NA_character_
sample_data(ps) <- sample_data(sd)

ps_ctrl <- prune_samples(sd$is_control & !is.na(sd$ctrl_type) &
                         sd$ctrl_type != "sea_water_control", ps)
ps_ctrl <- prune_taxa(taxa_sums(ps_ctrl) > 0, ps_ctrl)
ps_ctrl_rel <- transform_sample_counts(ps_ctrl, function(x) {
  s <- sum(x); if (s == 0) x else x / s
})
bc <- phyloseq::distance(ps_ctrl_rel, method = "bray")
meta_ctrl <- data.frame(sample_data(ps_ctrl_rel))
## adonis2 position-matches rows: make sure metadata follows the distance labels
meta_ctrl <- meta_ctrl[attr(bc, "Labels"), , drop = FALSE]
stopifnot(identical(rownames(meta_ctrl), attr(bc, "Labels")))
meta_ctrl$ctrl_type <- factor(meta_ctrl$ctrl_type)

## Print the control-type PERMANOVA for checking only (not plotted)
set.seed(100)
print(adonis2(bc ~ ctrl_type, data = meta_ctrl, permutations = 9999, by = "margin"))
print(table(meta_ctrl$ctrl_type, meta_ctrl$region))

pcoa <- cmdscale(bc, k = 2, eig = TRUE)
eig  <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)
pdat <- data.frame(PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2],
                   ctrl_type = meta_ctrl$ctrl_type,
                   region    = factor(meta_ctrl$region, levels = site_order),
                   sample    = rownames(meta_ctrl))

ctrl_labels <- c("buffer_control"     = "Buffer control",
                 "coriolis_control"   = "Coriolis control",
                 "deposition_control" = "Deposition control")

pa <- ggplot(pdat, aes(PC1, PC2, color = ctrl_type, shape = region)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = 30,
                  show.legend = FALSE, seed = 100, segment.size = 0.2) +
  scale_color_manual(values = pal_ctrl_type, labels = ctrl_labels,
                      name = "Control type") +
  scale_shape_discrete(labels = site_labels, name = "Site") +
  labs(x = paste0("PCoA1 (", eig[1], "%)"),
       y = paste0("PCoA2 (", eig[2], "%)")) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        legend.key.size  = unit(3.5, "mm"),
        legend.text      = element_text(size = 6.5),
        legend.title     = element_text(size = 7))

## ---- panel b: betaNTI by phase pair, FULL vs DECON (output of decontam_betaNTI_01OCT26.R) --
df <- read_tsv(file.path(fig_dir, "decontam_betaNTI_01OCT26.tsv"),
               show_col_types = FALSE) |>
  mutate(set = factor(set, c("FULL", "DECON"),
                      labels = c("All ASVs", "decontam-filtered")),
         pair_type = factor(pair_type, c("Before-Before", "Before-D+A", "D+A-D+A")))
print(df |> group_by(set, pair_type) |>
        summarize(n = n(), median = round(median(bNTI), 2), .groups = "drop"))

pb <- ggplot(df, aes(pair_type, bNTI, fill = pair_type)) +
  geom_hline(yintercept = c(-2, 2), linetype = 2, color = "gray50") +
  geom_hline(yintercept = 0, linetype = 1, color = "gray80") +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85, linewidth = 0.3) +
  facet_wrap(~ set) +
  scale_fill_manual(values = pal_bnti_pair, guide = "none") +
  labs(x = NULL, y = expression(beta*"NTI")) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray92", color = NA),
        axis.text.x = element_text(angle = 30, hjust = 1))

fig <- (pa | pb) + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))

ggsave(file.path(fig_dir, paste0("figureS8_decontam_checks_", out_tag, ".pdf")),
       fig, width = 183, height = 85, units = "mm", device = cairo_pdf)
ggsave(file.path(fig_dir, paste0("figureS8_decontam_checks_", out_tag, ".png")),
       fig, width = 183, height = 85, units = "mm", dpi = 300)
cat("Saved figureS8_decontam_checks_", out_tag, ".{pdf,png}\n", sep = "")
