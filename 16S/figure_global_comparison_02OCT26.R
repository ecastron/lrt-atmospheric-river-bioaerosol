# figure_global_comparison_02OCT26.R
# Figure 1. Places the LRT air communities among Southern Ocean air (Malard et al. 2022) and soils: weighted UniFrac
# PCoA of all samples (a) and of air samples only (b), with PERMANOVA (adonis2) and a Mantel test of weighted
# UniFrac against geographic distance. Controls and the Union Glacier samples are excluded.
# Input  : 16S/ps_merged_withTree_noorg_01OCT26.RDS, 16S/palettes.R
# Outputs: 16S/figures/global_comparison_02OCT26.{pdf,png}, 16S/figures/global_comparison_stats_02OCT26.txt
# Run    : Rscript 16S/figure_global_comparison_02OCT26.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(geosphere)
  library(scales)
})

set.seed(100)
base     <- "./16S"
in_tag   <- "noorg_01OCT26"   # input phyloseq object (organelle-free)
date_tag <- "02OCT26"   # output tag
fig_dir  <- file.path(base, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

source(file.path(base, "palettes.R"))

# Display labels only; colors and shapes still come from palettes.R.
lab_source2 <- c(
  "Antarctic Air"    = "Antarctic Peninsula air",
  "Patagonian Air"   = "Patagonian air",
  "SO Air"           = "Southern Ocean air, over sea",
  "SO Air over land" = "Southern Ocean air, over land",
  "SO Precipitation" = "Southern Ocean precipitation",
  "Antarctic Soil"   = "Antarctic Peninsula soil",
  "Patagonian Soil"  = "Patagonian soil"
)
lab_dataset <- c("LRT" = "This study", "Malard2022" = "Malard et al. 2022")

ps_path <- file.path(base, paste0("ps_merged_withTree_", in_tag, ".RDS"))
stopifnot(file.exists(ps_path))
ps <- readRDS(ps_path)

# --- 1) clean: drop controls and samples without coordinates ----------
md <- as(sample_data(ps), "data.frame")
keep <- !md$is_control &
        !is.na(md$latitude) & !is.na(md$longitude) &
        md$source2 != "Control" &
        !(md$region %in% "Union Glacier")   # separate campaign, not on the transect; %in% keeps NA regions
stopifnot(sum(md$region %in% "Union Glacier" & !md$is_control) == 6)
ps <- prune_samples(keep, ps)
cat("After dropping controls / no-coords: ",
    nsamples(ps), " samples\n", sep = "")

# --- 2) drop low-depth samples (no rarefaction) -----------------------
# weighted UniFrac in phyloseq normalizes counts to fractions internally,
# so library-size differences are handled. We just remove samples too
# shallow to give stable proportion estimates.
depths <- sample_sums(ps)
cat("Read depth quantiles:\n")
print(round(quantile(depths, c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.95, 1))))

min_depth <- 1000
ps_r <- prune_samples(sample_sums(ps) >= min_depth, ps)
ps_r <- prune_taxa(taxa_sums(ps_r) > 0, ps_r)
cat("After min-depth filter (", min_depth, " reads): ",
    nsamples(ps_r), " samples, ", ntaxa(ps_r), " ASVs\n", sep = "")

# --- 3) weighted UniFrac on the full set ------------------------------
cat("Computing weighted UniFrac (this can take a few minutes)...\n")
unif_w <- UniFrac(ps_r, weighted = TRUE, parallel = FALSE)

# --- 4) PCoA + ggplot data --------------------------------------------
ord_all <- ordinate(ps_r, method = "PCoA", distance = unif_w)
eig <- ord_all$values$Relative_eig
ax_lab <- function(i) sprintf("PCoA%d (%.1f%%)", i, 100 * eig[i])

md_r <- as(sample_data(ps_r), "data.frame") |>
  mutate(sample_id = sample_names(ps_r))
scores_all <- ord_all$vectors[, 1:2] |>
  as.data.frame() |>
  tibble::rownames_to_column("sample_id") |>
  rename(Axis1 = Axis.1, Axis2 = Axis.2) |>
  left_join(md_r, by = "sample_id")

# --- 5) Panel A: all samples, color by source2, shape by dataset ------
pA <- ggplot(scores_all,
             aes(Axis1, Axis2, color = source2, shape = dataset)) +
  geom_point(size = 1.1, alpha = 0.8) +
  stat_ellipse(aes(group = source2), linewidth = 0.3,
               linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(values = pal_source2,
                     limits = setdiff(names(pal_source2), "Control"),
                     labels = lab_source2) +
  scale_shape_manual(values = shp_dataset, labels = lab_dataset) +
  labs(x = ax_lab(1), y = ax_lab(2),
       color = "Sample type", shape = "Dataset",
       subtitle = "All samples") +
  guides(color = guide_legend(ncol = 2, order = 1, override.aes = list(size = 2)),
         shape = guide_legend(ncol = 1, order = 2,
                              override.aes = list(size = 2, alpha = 1,
                                                  color = "gray30"))) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank())

# --- 7) Air-only subset for Panel B + Mantel --------------------------
air_lvls <- c("Antarctic Air", "Patagonian Air",
              "SO Air", "SO Air over land")
ps_air <- prune_samples(sample_data(ps_r)$source2 %in% air_lvls, ps_r)
ps_air <- prune_taxa(taxa_sums(ps_air) > 0, ps_air)
cat("Air-only subset: ", nsamples(ps_air), " samples\n", sep = "")

unif_w_air <- UniFrac(ps_air, weighted = TRUE, parallel = FALSE)
ord_air    <- ordinate(ps_air, method = "PCoA", distance = unif_w_air)
eig_air    <- ord_air$values$Relative_eig
ax_lab_air <- function(i) sprintf("PCoA%d (%.1f%%)", i, 100 * eig_air[i])

md_air <- as(sample_data(ps_air), "data.frame") |>
  mutate(sample_id = sample_names(ps_air))
scores_air <- ord_air$vectors[, 1:2] |>
  as.data.frame() |>
  tibble::rownames_to_column("sample_id") |>
  rename(Axis1 = Axis.1, Axis2 = Axis.2) |>
  left_join(md_air, by = "sample_id")

# Mantel: geographic (Haversine, km) vs weighted UniFrac, air only
coords <- as.matrix(md_air[, c("longitude", "latitude")])
geo_d  <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
geo_d  <- as.dist(geo_d)
mantel_air <- vegan::mantel(geo_d, unif_w_air,
                            method = "spearman",
                            permutations = 999)

# --- 7b) PERMANOVA: three tests ---------------------------------------
# (a) Overall: how much weighted UniFrac variance does source2 explain
#     across all samples? Drops `dataset` (perfectly nested in source2).
# (b) Air-only: within air samples, does region (Antarctic / Patagonian /
#     SO) structure communities? Categorical companion to the Mantel.
# (c) LRT-only factorial: medium (Air vs Soil) crossed with region
#     (Antarctic vs Patagonian): the clean test of the air-vs-soil
#     claim with region adjusted for.
permanova_overall <- adonis2(
  as.dist(unif_w) ~ source2,
  data         = md_r,
  permutations = 999
)

md_air_perm <- md_air |>
  mutate(region_air = case_when(
    source2 == "Antarctic Air"                   ~ "Antarctic",
    source2 == "Patagonian Air"                  ~ "Patagonian",
    source2 %in% c("SO Air", "SO Air over land") ~ "SO",
    TRUE                                         ~ NA_character_
  ))
permanova_air <- adonis2(
  as.dist(unif_w_air) ~ region_air,
  data         = md_air_perm,
  permutations = 999
)

ufm <- as.matrix(unif_w)
md_lrt <- md_r |>
  filter(dataset == "LRT",
         source2 %in% c("Antarctic Air", "Antarctic Soil",
                        "Patagonian Air", "Patagonian Soil")) |>
  mutate(
    medium = if_else(grepl("Air$",  source2), "Air", "Soil"),
    region = if_else(grepl("^Antarctic", source2), "Antarctic", "Patagonian")
  )
unif_lrt <- as.dist(ufm[md_lrt$sample_id, md_lrt$sample_id])
# Additive model gives clean marginal main effects. The interaction
# isn't zero (visible in the ordination), but for the opener-figure
# claim ("air vs soil regardless of region") the main effects are the
# headline. Interaction reported separately for completeness.
permanova_lrt <- adonis2(
  unif_lrt ~ medium + region,
  data         = md_lrt,
  permutations = 999,
  by           = "margin"
)
permanova_lrt_int <- adonis2(
  unif_lrt ~ medium * region,
  data         = md_lrt,
  permutations = 999,
  by           = "terms"
)

# --- 8) Panel B: air-only, color by latitude --------------------------
pB <- ggplot(scores_air,
             aes(Axis1, Axis2, color = latitude, shape = dataset)) +
  geom_point(size = 1.1, alpha = 0.85) +
  scale_color_viridis_c(option = pal_latitude_opt,
                        direction = pal_latitude_dir,
                        breaks = pretty_breaks(5)) +
  scale_shape_manual(values = shp_dataset, labels = lab_dataset) +
  labs(x = ax_lab_air(1), y = ax_lab_air(2),
       color = "Latitude (°)", shape = "Dataset",
       subtitle = "Air samples") +
  guides(color = guide_colorbar(barwidth = unit(2.5, "cm"),
                                barheight = unit(0.25, "cm"),
                                title.vjust = 0.8),
         shape = "none") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank())

# --- 9) compose & save ------------------------------------------------
# Collect all legends into one bottom row so panel legends don't collide.
fig <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag        = element_text(face = "bold", size = 10),
        plot.subtitle   = element_text(size = 8),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title    = element_text(size = 7),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

pdf_path <- file.path(fig_dir, paste0("global_comparison_", date_tag, ".pdf"))
png_path <- file.path(fig_dir, paste0("global_comparison_", date_tag, ".png"))
ggsave(pdf_path, fig, width = 183, height = 110, units = "mm")
ggsave(png_path, fig, width = 183, height = 110, units = "mm", dpi = 300)
cat("Saved:\n  ", pdf_path, "\n  ", png_path, "\n", sep = "")

# --- 10) write stats --------------------------------------------------
stats_path <- file.path(fig_dir,
                        paste0("global_comparison_stats_", date_tag, ".txt"))
sink(stats_path)
cat("Global comparison - opener figure stats\n")
cat("Date: ", date_tag, "\n", sep = "")
cat("Input: ", basename(ps_path), "\n", sep = "")
cat("Min-depth filter: ", min_depth, " reads (no rarefaction)\n\n", sep = "")
cat("=== PERMANOVA (a): all samples ~ source2 ===\n")
print(permanova_overall)
cat("\n=== PERMANOVA (b): air-only ~ region (Ant/Pat/SO) ===\n")
print(permanova_air)
cat("\n=== PERMANOVA (c): LRT-only ~ medium + region (additive, marginal) ===\n")
print(permanova_lrt)
cat("\n=== PERMANOVA (c'): LRT-only ~ medium * region (sequential, for interaction) ===\n")
print(permanova_lrt_int)
cat("\n=== Mantel: air-only weighted UniFrac vs Haversine, Spearman ===\n")
print(mantel_air)
cat("\nSample counts by source2 (post-min-depth, post-control):\n")
print(table(md_r$source2, md_r$dataset))
cat("\nLRT factorial cell sizes (medium x region):\n")
print(table(md_lrt$medium, md_lrt$region))
sink()
cat("Stats: ", stats_path, "\n", sep = "")
