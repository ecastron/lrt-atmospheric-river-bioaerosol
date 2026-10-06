# figure5_ar_event_01OCT26.R
# Figure 5. Bacterial community response to the 7-8 February 2022 atmospheric river at Risopatron (24 daily air
# samples): weighted UniFrac PCoA trajectory, alpha diversity (Chao1, Shannon, Faith's PD), phylogenetic structure
# (MPD, NRI, NTI), genus heatmap, PERMANOVA by phase, and the Stegen et al. (2013) assembly-process partition
# (betaNTI followed by abundance-based Raup-Crick on Bray-Curtis). betaNTI null permutations are drawn serially from
# set.seed(100) so results are reproducible. The AR core window (7-8 February, UTC) follows the ERA5 IVT analysis.
# Input  : 16S/ps_merged_withTree_noorg_01OCT26.RDS, 16S/palettes.R
# Outputs: 16S/figures/figure5_ar_event_01OCT26.{pdf,png}, 16S/figures/figure5_ar_event_stats_01OCT26.txt,
#          pairwise betaNTI and RC_Bray table; cached null results are written next to the outputs
# Run    : Rscript 16S/figure5_ar_event_01OCT26.R

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
  library(ggrepel)
  library(tibble)
  library(readr)
  library(parallel)
})

set.seed(100)
base       <- "./16S"
date_tag   <- "01OCT26"   # output tag
ps_tag     <- "noorg_01OCT26"   # organelle-free object
fig_dir    <- file.path(base, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

source(file.path(base, "palettes.R"))

ps_path <- file.path(base, paste0("ps_merged_withTree_", ps_tag, ".RDS"))
stopifnot(file.exists(ps_path))
ps <- readRDS(ps_path)

# Inline log10(x+1) transform.
log1p10 <- scales::trans_new(
  name      = "log1p10",
  transform = function(x) log10(x + 1),
  inverse   = function(x) 10^x - 1
)

# ============================================================
# 1) Subset to Risopatron Air (24 daily samples, RA1-RA24)
# ============================================================
md <- as(sample_data(ps), "data.frame")
keep <- md$dataset == "LRT" & md$region == "Risopatron" &
        md$sampler == "Coriolis_mu" & !md$is_control
ps_ra <- prune_samples(keep, ps)
ps_ra <- prune_taxa(taxa_sums(ps_ra) > 0, ps_ra)

sd_ra <- as(sample_data(ps_ra), "data.frame") |>
  mutate(date  = as.Date(date),
         phase = case_when(
           date <  ar_event_start ~ "Before",
           date <= ar_event_end   ~ "During",
           TRUE                   ~ "After"),
         phase = factor(phase, c("Before", "During", "After"))) |>
  arrange(date)

ps_ra <- prune_samples(rownames(sd_ra), ps_ra)
sample_data(ps_ra) <- sample_data(sd_ra)

stopifnot(
  nsamples(ps_ra) == 24,
  all(c("RA1", "RA16", "RA17", "RA24") %in% sample_names(ps_ra)),
  sum(sd_ra$phase == "During") == 2,
  sum(sd_ra$atm_river == "yes") == 8,   # RA13-RA20, Feb 5-11
  ape::is.rooted(phy_tree(ps_ra))
)
cat("Risopatron Air: ", nsamples(ps_ra), " samples, ",
    ntaxa(ps_ra), " ASVs (after pruning zeros)\n", sep = "")

# Reusable bits
sample_levels <- sd_ra$sample_id        # chronological order
phase_by_id   <- setNames(as.character(sd_ra$phase), sd_ra$sample_id)
date_by_id    <- setNames(sd_ra$date,                sd_ra$sample_id)
# 2-level grouping (Before vs During+After) used by βNTI, PERMANOVA-pairwise
# and per-genus Wilcoxon. Defined here so all downstream blocks can use it.
sd_ar2 <- sd_ra |> mutate(ar2 = if_else(phase == "Before",
                                        "Before", "During+After"))

# ============================================================
# 2) Panel A: Weighted UniFrac PCoA trajectory
# ============================================================
ps_rel <- transform_sample_counts(ps_ra, function(x) {
  s <- sum(x); if (s == 0) x else x / s
})
D_wu  <- phyloseq::UniFrac(ps_rel, weighted = TRUE, normalized = TRUE)
ord   <- ape::pcoa(D_wu)

ax_pct <- 100 * ord$values$Relative_eig[1:2]
coords <- ord$vectors[, 1:2] |>
  as.data.frame() |>
  rownames_to_column("sample_id") |>
  rename(Axis1 = Axis.1, Axis2 = Axis.2) |>
  left_join(sd_ra |> select(sample_id, date, phase, atm_river),
            by = "sample_id") |>
  arrange(date) |>
  mutate(sample_id = factor(sample_id, levels = sample_levels))

pA <- ggplot(coords, aes(Axis1, Axis2)) +
  geom_path(arrow = arrow(length = unit(2.5, "pt"), type = "closed"),
            color = "gray40", linewidth = 0.4) +
  geom_point(aes(fill = phase), shape = 21, color = "black",
             size = 2, stroke = 0.3) +
  scale_fill_manual(values = pal_ar_phase, name = "Phase") +
  labs(x = sprintf("Axis 1 (%.1f%%)", ax_pct[1]),
       y = sprintf("Axis 2 (%.1f%%)", ax_pct[2])) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank())

# ============================================================
# 3) Panel B: Alpha diversity time series (Chao1, Shannon)
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

# Shade only the ERA5-defined AR core (7-8 Feb), half a day either
# side so both sample days sit inside the band (one AR definition throughout).
ar_rect <- annotate("rect",
                    xmin = ar_event_start - 0.5, xmax = ar_event_end + 0.5,
                    ymin = -Inf, ymax = Inf,
                    fill = "#C44E52", alpha = 0.10)

pB <- ggplot(alpha_long, aes(date, value)) +
  ar_rect +
  geom_line(aes(group = metric), color = "gray50", linewidth = 0.4) +
  geom_point(aes(fill = phase), shape = 21, color = "black",
             size = 1.8, stroke = 0.25) +
  scale_fill_manual(values = pal_ar_phase, name = "Phase", guide = "none") +
  scale_x_date(date_breaks = "6 days", date_labels = "%b %d") +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             strip.position = "left") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(strip.placement    = "outside",
        strip.background.y = element_rect(fill = "gray92", color = NA),
        strip.text.y.left  = element_text(angle = 90, face = "bold"),
        panel.grid.minor   = element_blank())

# ============================================================
# 4) Panel C: Genus heatmap (top 20 x RA1-RA24)
# ============================================================
ps_g <- tax_glom(ps_ra, taxrank = "Genus", NArm = FALSE)
ps_g <- transform_sample_counts(ps_g, function(x) {
  s <- sum(x); if (s == 0) x else x / s
})

# Top 20 genera by total abundance; lump rest as "Other"
genus_tax  <- as.data.frame(tax_table(ps_g))
genus_tax$Genus_lbl <- ifelse(is.na(genus_tax$Genus) | genus_tax$Genus == "",
                              "Unassigned", genus_tax$Genus)
genus_sums <- taxa_sums(ps_g)
top_idx    <- order(genus_sums, decreasing = TRUE)[1:min(20, length(genus_sums))]
top_genera <- genus_tax$Genus_lbl[top_idx]

heat_long <- psmelt(ps_g) |>
  mutate(Genus_lbl = ifelse(is.na(Genus) | Genus == "",
                            "Unassigned", Genus),
         Taxon     = ifelse(Genus_lbl %in% top_genera,
                            Genus_lbl, "Other")) |>
  group_by(Sample, Taxon) |>
  summarize(Abundance = sum(Abundance), .groups = "drop") |>
  filter(Sample %in% sample_levels) |>
  mutate(Sample = factor(Sample, levels = sample_levels))

# Order taxa by total abundance (descending), pin "Other" to bottom
tax_order <- heat_long |>
  group_by(Taxon) |>
  summarize(tot = sum(Abundance), .groups = "drop") |>
  arrange(desc(tot)) |>
  pull(Taxon)
tax_order <- c(setdiff(tax_order, "Other"), intersect("Other", tax_order))
heat_long <- heat_long |>
  mutate(Taxon = factor(Taxon, levels = rev(tax_order)))

# AR window box: indices of atm_river=="yes" samples in the X order
yes_idx <- which(sd_ra$phase == "During")   # AR core samples (RA16, RA17)
ar_box  <- annotate("rect",
                    xmin = min(yes_idx) - 0.5,
                    xmax = max(yes_idx) + 0.5,
                    ymin = 0.5,
                    ymax = length(levels(heat_long$Taxon)) + 0.5,
                    color = "black", fill = NA,
                    linewidth = 0.7)

# Date labels for the X axis, matched to the chronological sample order.
# Keep positional X (24 columns, one per sample) and relabel breaks
# with dates so duplicate-date samples (RA14/RA15, RA23/RA24) still
# render as distinct tiles.
date_labels <- format(sd_ra$date, "%b %d")
names(date_labels) <- sd_ra$sample_id

pC_main <- ggplot(heat_long, aes(Sample, Taxon, fill = Abundance * 100)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#08306B", "#6BAED6", "#FFFFE5", "#FB6A4A", "#A50F15"),
    trans   = log1p10,
    breaks  = c(0.01, 0.1, 1, 10),
    labels  = function(x) paste0(x, "%"),
    limits  = c(0.01, 10),
    oob     = scales::squish,
    name    = "Relative\nabundance\n(log10+1)"
  ) +
  scale_x_discrete(labels = date_labels) +
  ar_box +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "right",
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.5, "mm"))

# --- Source attribution annotation (left-side panel) -----------------
# For each of the top genera, compute prevalence in 5 candidate source
# groups using the FULL merged phyloseq (not just Risopatron Air).
# Source groups mirror ar_patagonia_attribution_v3_01OCT26.R.
md_full <- as(sample_data(ps), "data.frame")
patagonia_regions <- c("PA_POR", "PN_SG", "Puerto Williams")
src_groups_attr <- list(
  "Patagonia air"     = rownames(md_full)[md_full$dataset == "LRT" &
                                          md_full$region %in% patagonia_regions &
                                          md_full$sampler == "Coriolis_mu" &
                                          !md_full$is_control],
  "Patagonia soil"    = rownames(md_full)[md_full$dataset == "LRT" &
                                          md_full$region %in% patagonia_regions &
                                          md_full$sampler %in% c("soil","rhizosphere") &
                                          !md_full$is_control],
  "Southern Ocean air" = rownames(md_full)[md_full$dataset == "Malard2022" &
                                          md_full$source2 == "SO Air"],
  "Southern Ocean air over land" = rownames(md_full)[md_full$dataset == "Malard2022" &
                                          md_full$source2 == "SO Air over land"],
  "Risopatrón soil" = rownames(md_full)[md_full$dataset == "LRT" &
                                          md_full$region == "Risopatron" &
                                          md_full$sampler %in% c("soil","rhizosphere") &
                                          !md_full$is_control]
)

tax_full <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
otu_full_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu_full_mat <- t(otu_full_mat)

genus_prev <- function(g, src_ids) {
  if (length(src_ids) == 0) return(NA_real_)
  asvs <- rownames(tax_full)[!is.na(tax_full$Genus) & tax_full$Genus == g]
  if (length(asvs) == 0) return(NA_real_)
  m <- otu_full_mat[src_ids, asvs, drop = FALSE]
  mean(rowSums(m) > 0)
}

attr_rows <- expand.grid(Taxon = top_genera,
                         Source = names(src_groups_attr),
                         stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(prev = genus_prev(Taxon, src_groups_attr[[Source]])) |>
  ungroup() |>
  mutate(Taxon  = factor(Taxon,  levels = levels(heat_long$Taxon)),
         Source = factor(Source, levels = names(src_groups_attr)))
# Keep the empty "Other" row so the two heatmaps share row positions

pC_attr <- ggplot(attr_rows, aes(Source, Taxon, fill = prev)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(prev) | prev == 0, "",
                               sprintf("%.0f", prev * 100))),
            size = 1.8, color = "gray20") +
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603",
                      limits = c(0, 1),
                      labels = scales::percent_format(),
                      na.value = "gray92",
                      name = "Source\nprevalence") +
  # italic only for genus names (incl. composites such as Methylobacterium-Methylorubrum); "Unassigned",
  # "Other" and clade names such as "OM60(NOR5) clade" stay upright
  scale_y_discrete(drop = FALSE, labels = function(x) as.expression(lapply(x, function(s)
    if (grepl("^[A-Z][a-z]+(-[A-Z][a-z]+)*$", s) && !s %in% c("Unassigned", "Other"))
      bquote(italic(.(s))) else bquote(plain(.(s)))))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 6.5),
        axis.text.y      = element_text(size = 6.5),
        panel.grid       = element_blank(),
        legend.position  = "right",
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.5, "mm"))

# Combine attribution + main heatmap with patchwork. wrap_elements()
# treats the pair as a single tagged unit so outer A/B/C/D tags stay
# correct.
pC <- patchwork::wrap_elements(
  pC_attr + pC_main + plot_layout(widths = c(0.3, 1.0))
)

# ============================================================
# 5) Panel D: Phylogenetic measures (Faith's PD + MPD + NRI)
# ============================================================
# Three complementary metrics:
#   * Faith's PD: total branch length (richness scale)
#   * MPD: mean pairwise phylogenetic distance among taxa
#                 (deep-tree breadth, richness-independent)
#   * NRI: -ses.mpd.z; standardized clustering vs taxa.labels null
#                 (positive = clustered, negative = overdispersed)
# NTI (-ses.mntd.z) is computed and reported in the stats file but is
# not shown in panel D: it is dominated by nearest-taxon distances
# among dense reagent-contaminant ASV clusters (in the decontam analysis,
# NTI loses significance after removing 408 reagent-contaminant ASVs,
# while MPD/NRI/Faith's PD remain robust).
otu_mat <- as(otu_table(ps_ra), "matrix")
if (taxa_are_rows(ps_ra)) otu_mat <- t(otu_mat)
otu_pa  <- (otu_mat > 0) + 0L
otu_pa  <- otu_pa[sample_levels, , drop = FALSE]   # chronological

tr <- phy_tree(ps_ra)
co <- cophenetic(tr)

cat("Computing Faith's PD ...\n")
pd_df <- picante::pd(otu_pa, tr, include.root = TRUE) |>
  rownames_to_column("sample_id") |>
  rename(Faith_PD = PD) |>
  select(sample_id, Faith_PD)

# ses.mpd / ses.mntd run serially after set.seed(100), so they are deterministic;
# cached for layout-only reruns (delete the file to recompute).
ses_cache <- file.path(fig_dir, paste0("figure5_sesmodels_", date_tag, ".rds"))
if (file.exists(ses_cache)) {
  cat("Loading cached ses.mpd / ses.mntd: ", ses_cache, "\n", sep = "")
  sc <- readRDS(ses_cache); ses_mpd <- sc$ses_mpd; ses_mntd <- sc$ses_mntd
} else {
  cat("Computing MPD + NRI via ses.mpd (999 perms) ...\n")
  set.seed(100)
  ses_mpd <- picante::ses.mpd(otu_pa, co,
                              null.model = "taxa.labels",
                              runs = 999,
                              abundance.weighted = FALSE)
  cat("Computing NTI via ses.mntd (999 perms; supplementary) ...\n")
  set.seed(100)
  ses_mntd <- picante::ses.mntd(otu_pa, co,
                                null.model = "taxa.labels",
                                runs = 999,
                                abundance.weighted = FALSE)
  saveRDS(list(ses_mpd = ses_mpd, ses_mntd = ses_mntd), ses_cache)
}
mpd_df <- data.frame(sample_id = rownames(ses_mpd),
                     MPD       = ses_mpd$mpd.obs)
nri_df <- data.frame(sample_id = rownames(ses_mpd),
                     NRI       = -1 * ses_mpd$mpd.obs.z)

nti_df <- data.frame(sample_id = rownames(ses_mntd),
                     NTI       = -1 * ses_mntd$mntd.obs.z)

phy_long <- pd_df |>
  inner_join(mpd_df, by = "sample_id") |>
  inner_join(nri_df, by = "sample_id") |>
  inner_join(sd_ra |> select(sample_id, date, phase),
             by = "sample_id") |>
  pivot_longer(c(Faith_PD, MPD, NRI),
               names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric,
                         levels = c("Faith_PD", "MPD", "NRI"),
                         labels = c("Faith's PD", "MPD",
                                    "NRI")))

pD <- ggplot(phy_long, aes(date, value)) +
  ar_rect +
  geom_hline(data = data.frame(metric = factor("NRI",
                                               levels = levels(phy_long$metric))),
             aes(yintercept = 0),
             linetype = 3, color = "gray40") +
  geom_line(aes(group = metric), color = "gray50", linewidth = 0.4) +
  geom_point(aes(fill = phase), shape = 21, color = "black",
             size = 1.8, stroke = 0.25) +
  scale_fill_manual(values = pal_ar_phase, name = "Phase", guide = "none") +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d") +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             strip.position = "left") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(strip.placement    = "outside",
        strip.background.y = element_rect(fill = "gray92", color = NA),
        strip.text.y.left  = element_text(angle = 90, face = "bold"),
        panel.grid.minor   = element_blank())

# ============================================================
# 5b) βNTI computation (used in Panel E and in stats output)
# ============================================================
# Stegen et al. 2013: βNTI = (βMNTD_obs - mean(βMNTD_null)) / sd(βMNTD_null).
# βMNTD = abundance-weighted comdistnt; null = tip-label shuffle on the
# cophenetic distance matrix, 999 perms. |βNTI|>2 = selection (positive
# = variable / heterogeneous; negative = homogeneous); |βNTI|<=2 =
# stochastic processes dominate.
cat("Computing βNTI (Stegen 2013; abundance-weighted; 999 perms) ...\n")
otu_full <- as(otu_table(ps_ra), "matrix")
if (taxa_are_rows(ps_ra)) otu_full <- t(otu_full)
otu_full <- otu_full[sample_levels, , drop = FALSE]

ncores_bnti <- max(1, parallel::detectCores() - 2)
NB_RUNS     <- 999

# The βNTI and RC_Bray null models take ~10 min; both are deterministic given
# seed 100, so their results are cached and reused for layout-only reruns.
# Delete the cache file to recompute.
null_cache <- file.path(fig_dir, paste0("figure5_nullmodels_", date_tag, ".rds"))
if (file.exists(null_cache)) {
  cat("Loading cached βNTI / RC_Bray null models: ", null_cache, "\n", sep = "")
  cache    <- readRDS(null_cache)
  bnti_mat <- cache$bnti_mat
} else {
  obs_bmntd <- as.matrix(picante::comdistnt(otu_full, co,
                                            abundance.weighted = TRUE,
                                            exclude.conspecifics = FALSE))
  # The tip-label permutations are drawn serially from one seed and only then
  # handed to mclapply. Drawing them inside the forked workers would give a
  # different null (and slightly different βNTI) for every run and core count.
  set.seed(100)
  perm_list <- replicate(NB_RUNS, sample(colnames(co)), simplify = FALSE)
  null_list <- parallel::mclapply(perm_list, function(rn) {
    rc <- co
    rownames(rc) <- rn; colnames(rc) <- rn
    as.matrix(picante::comdistnt(otu_full, rc,
                                 abundance.weighted = TRUE,
                                 exclude.conspecifics = FALSE))
  }, mc.cores = ncores_bnti, mc.preschedule = TRUE)
  null_arr <- simplify2array(null_list)
  bnti_mat <- (obs_bmntd - apply(null_arr, c(1, 2), mean)) /
              apply(null_arr, c(1, 2), sd)
}
bnti_mat[is.na(bnti_mat)] <- 0; diag(bnti_mat) <- NA   # NA: ignore self

# Per-sample summary: median βNTI to all Before samples (leave-one-out
# when the focal sample itself is in Before).
before_ids <- sd_ra$sample_id[sd_ra$phase == "Before"]
bnti_vs_before <- sapply(sample_levels, function(s) {
  others <- setdiff(before_ids, s)
  median(bnti_mat[s, others], na.rm = TRUE)
})
bnti_panel_df <- data.frame(
  sample_id        = sample_levels,
  median_bNTI_Bef  = unname(bnti_vs_before),
  date  = sd_ra$date[match(sample_levels, sd_ra$sample_id)],
  phase = sd_ra$phase[match(sample_levels, sd_ra$sample_id)],
  stringsAsFactors = FALSE)

# Pair-level long form (for stats summary section later)
ut       <- which(upper.tri(bnti_mat), arr.ind = TRUE)
bnti_df  <- data.frame(
  A     = rownames(bnti_mat)[ut[, 1]],
  B     = colnames(bnti_mat)[ut[, 2]],
  bNTI  = bnti_mat[upper.tri(bnti_mat)],
  A_ar2 = sd_ar2[rownames(bnti_mat)[ut[, 1]], "ar2"],
  B_ar2 = sd_ar2[colnames(bnti_mat)[ut[, 2]], "ar2"],
  stringsAsFactors = FALSE) |>
  mutate(pair_type = case_when(
    A_ar2 == "Before" & B_ar2 == "Before" ~ "Before-Before",
    A_ar2 != "Before" & B_ar2 != "Before" ~ "D+A-D+A",
    TRUE                                  ~ "Before-D+A"),
    pair_type = factor(pair_type,
                       c("Before-Before","Before-D+A","D+A-D+A")))

# ============================================================
# 5c) RC_Bray (Stegen et al. 2013, step 2)
# ============================================================
# Abundance-based Raup-Crick on Bray-Curtis (Stegen's raup_crick_abundance):
# regional pool = ASVs in the 24 samples; each null community keeps the
# observed richness (ASVs drawn without replacement, prob. ~ occupancy) and
# read total (one read per ASV, the rest ~ regional abundance).
# RC_Bray = 2 * (P(null < obs) + 0.5 P(null == obs)) - 1.
# Null communities are drawn per sample within each replicate, which gives
# each pair the same null distribution as Stegen's per-pair draws.
cat("Computing RC_Bray (999 randomizations) ...\n")
otu_int <- otu_full; storage.mode(otu_int) <- "integer"
bc_obs  <- as.matrix(vegan::vegdist(otu_int, method = "bray"))
occ     <- colSums(otu_int > 0); abund <- colSums(otu_int)
rich    <- rowSums(otu_int > 0); depth <- rowSums(otu_int)
n_asv   <- ncol(otu_int)
null_community <- function(S, N) {
  x  <- integer(n_asv)
  sp <- sample.int(n_asv, S, replace = FALSE, prob = occ)
  x[sp] <- 1L
  if (N > S) x[sp] <- x[sp] + as.integer(rmultinom(1, N - S, abund[sp]))
  x
}
if (exists("cache")) {
  rc_mat <- cache$rc_mat
} else {
  set.seed(100)
  n_lt <- n_eq <- matrix(0, nrow(otu_int), nrow(otu_int), dimnames = dimnames(bc_obs))
  for (r in seq_len(NB_RUNS)) {
    nul <- t(vapply(seq_len(nrow(otu_int)),
                    function(i) null_community(rich[i], depth[i]), integer(n_asv)))
    bc_null <- as.matrix(vegan::vegdist(nul, method = "bray"))
    n_lt <- n_lt + (bc_null < bc_obs); n_eq <- n_eq + (bc_null == bc_obs)
  }
  rc_mat <- 2 * ((n_lt + 0.5 * n_eq) / NB_RUNS) - 1
  saveRDS(list(bnti_mat = bnti_mat, rc_mat = rc_mat), null_cache)
}

process_levels <- names(pal_stegen_process)
bnti_df <- bnti_df |>
  mutate(RC_bray = rc_mat[cbind(A, B)],
         process = case_when(
           bNTI < -2       ~ "Homogeneous selection",
           bNTI >  2       ~ "Variable selection",
           RC_bray < -0.95 ~ "Homogenizing dispersal",
           RC_bray >  0.95 ~ "Dispersal limitation",
           TRUE            ~ "Drift (undominated)"),
         process = factor(process, process_levels))
stopifnot(!anyNA(bnti_df$RC_bray), nrow(bnti_df) == choose(24, 2))

proc_frac <- bnti_df |>
  count(pair_type, process, .drop = FALSE) |>
  group_by(pair_type) |>
  mutate(n_pairs = sum(n), pct = 100 * n / n_pairs) |>
  ungroup()
pair_lbl <- proc_frac |> distinct(pair_type, n_pairs) |>
  mutate(lbl = paste0(recode(as.character(pair_type),
                             "Before-Before" = "Within\nBefore",
                             "Before-D+A"    = "Before vs\nDuring+After",
                             "D+A-D+A"       = "Within\nDuring+After"),
                      "\n(n = ", n_pairs, ")"))
proc_frac <- proc_frac |>
  mutate(pair_type = factor(pair_type, rev(as.character(pair_lbl$pair_type)),
                            labels = rev(pair_lbl$lbl)))   # rev: first pair on top after coord_flip

# Panel E: Stegen process fractions per phase pair (Stegen estimates process
# influence as the fraction of pairwise comparisons, not as medians)
pE <- ggplot(proc_frac, aes(pair_type, pct, fill = process)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.2,
           position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(pct >= 8, sprintf("%.0f", pct), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.2,
            color = "black") +
  scale_fill_manual(values = pal_stegen_process, name = "Process",
                    drop = FALSE, guide = guide_legend(nrow = 2)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(x = NULL, y = "Pairwise comparisons") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(3, "mm"))

# ============================================================
# 6) Compose & save
# ============================================================
fig <- ((pA | pB) / pC / pD / pE) +
  plot_layout(heights = c(0.9, 2.1, 1.1, 0.75), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag        = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        legend.box      = "vertical",
        legend.text     = element_text(size = 7),
        legend.title    = element_text(size = 7.5))

pdf_path <- file.path(fig_dir, paste0("figure5_ar_event_", date_tag, ".pdf"))
png_path <- file.path(fig_dir, paste0("figure5_ar_event_", date_tag, ".png"))
ggsave(pdf_path, fig, width = 183, height = 260, units = "mm", device = cairo_pdf)
ggsave(png_path, fig, width = 183, height = 260, units = "mm", dpi = 300)
cat("Saved:\n  ", pdf_path, "\n  ", png_path, "\n", sep = "")

stopifnot(file.exists(pdf_path), file.size(png_path) > 200e3)

# ============================================================
# 7) Stats output
# ============================================================
stats_path <- file.path(fig_dir,
                        paste0("figure5_ar_event_stats_", date_tag, ".txt"))
sink(stats_path)
cat("Figure 5 stats - AR-event signal at Risopatron Air (script figure5_ar_event_01OCT26.R)\n")
cat("Date: ", date_tag, "  | input ps: ", ps_tag, "\n", sep = "")

cat("\n=== Sample counts by phase (strict UTC) ===\n")
print(table(sd_ra$phase))
cat("atm_river=='yes' (broad shading): ",
    sum(sd_ra$atm_river == "yes"), " samples\n", sep = "")

# --- 7a) PERMANOVA on weighted UniFrac ------------------------------
# vegan::adonis2 position-matches distance-row -> data-row when D's
# labels and data's rownames are not in the same order. sd_ra is
# chronological (after arrange(date)); D_wu inherits the phyloseq's
# alphabetical sample order. Align sd to D's labels before adonis2.
sd_perm <- sd_ra[attr(D_wu, "Labels"), ]
cat("\n=== PERMANOVA: wUniFrac ~ phase (9999 permutations) ===\n")
set.seed(100)
pm <- vegan::adonis2(D_wu ~ phase, data = sd_perm,
                     permutations = 9999, by = "margin")
print(pm)

cat("\n--- Pairwise PERMANOVA (Before vs During+After) ---\n")
sd_perm$ar2 <- if_else(sd_perm$phase == "Before",
                       "Before", "During+After")
set.seed(100)
pm2 <- vegan::adonis2(D_wu ~ ar2, data = sd_perm,
                      permutations = 9999, by = "margin")
print(pm2)

# --- 7b) Stepwise wUniFrac with cyclic-shift null --------------------
cat("\n=== Stepwise wUniFrac with cyclic-shift null (1999 perms) ===\n")
D_mat <- as.matrix(D_wu)[sample_levels, sample_levels]
n     <- length(sample_levels)
step_obs <- vapply(2:n,
                   function(i) D_mat[i, i - 1],
                   numeric(1))
set.seed(100)
nperm <- 1999
step_null <- replicate(nperm, {
  shift <- sample.int(n - 1, 1)
  idx   <- ((seq_len(n) - 1 + shift) %% n) + 1
  Ds    <- D_mat[idx, idx]
  vapply(2:n, function(i) Ds[i, i - 1], numeric(1))
})
# null: matrix (n-1) x nperm
mu  <- rowMeans(step_null)
sdv <- apply(step_null, 1, sd)
zsc <- (step_obs - mu) / sdv
pval <- rowMeans(step_null >= step_obs) + 1 / (nperm + 1)
step_tab <- tibble(
  step       = paste0(sample_levels[1:(n - 1)], "->",
                      sample_levels[2:n]),
  to_date    = sd_ra$date[2:n],
  to_phase   = sd_ra$phase[2:n],
  obs_dist   = round(step_obs, 4),
  null_mean  = round(mu, 4),
  z          = round(zsc, 3),
  one_sided_p = round(pval, 4)
)
print(step_tab, n = Inf)
cat("\nFlagged transitions:\n")
print(step_tab |>
        filter(grepl("RA15->RA16|RA16->RA17", step)))

# --- 7c) Wilcoxon Before vs During+After for alpha + phylo ----------
cat("\n=== Wilcoxon Before vs During+After (alpha + phylo) ===\n")
metrics <- alpha |>
  select(sample_id, Chao1, Shannon) |>
  inner_join(pd_df,  by = "sample_id") |>
  inner_join(mpd_df, by = "sample_id") |>
  inner_join(nri_df, by = "sample_id") |>
  inner_join(nti_df, by = "sample_id") |>
  inner_join(sd_ar2 |> select(sample_id, ar2), by = "sample_id")

for (m in c("Chao1", "Shannon", "Faith_PD", "MPD", "NRI", "NTI")) {
  v <- metrics[[m]]
  g <- factor(metrics$ar2, levels = c("Before", "During+After"))
  w  <- suppressWarnings(wilcox.test(v ~ g, exact = FALSE))
  n1 <- sum(g == "Before")
  n2 <- sum(g == "During+After")
  z  <- qnorm(1 - w$p.value / 2) * sign(diff(tapply(v, g, median)))
  r  <- abs(z) / sqrt(n1 + n2)
  cat(sprintf("%-12s  W=%.1f  p=%.4f  effect size r=%.2f  ",
              m, w$statistic, w$p.value, r))
  cat(sprintf("(median Before=%.3f, During+After=%.3f)\n",
              median(v[g == "Before"]),
              median(v[g == "During+After"])))
}
cat("\nNote (11MAY26): NTI is reported for transparency but loses\n")
cat("significance after removing 408 reagent-contaminant ASVs flagged\n")
cat("by decontam (see decontam_phylo_metrics_11MAY26.txt). NTI is\n")
cat("dominated by nearest-taxon distances among dense contaminant ASV\n")
cat("clusters (Pseudomonas, Stenotrophomonas). MPD and NRI rely on\n")
cat("the full pairwise tree and are robust to decontamination.\n")

# --- 7c-bis) βNTI summary (computed earlier in §5b; see Panel E) --------
cat("\n=== βNTI (Stegen et al. 2013, abundance-weighted, 999 perms) ===\n")
bnti_summary <- bnti_df |> group_by(pair_type) |>
  summarize(n        = n(),
            median   = round(median(bNTI), 2),
            q25      = round(quantile(bNTI, 0.25), 2),
            q75      = round(quantile(bNTI, 0.75), 2),
            hetero_pct = round(100 * mean(bNTI >   2), 1),
            homo_pct   = round(100 * mean(bNTI <  -2), 1),
            nonsel_pct = round(100 * mean(abs(bNTI) <= 2), 1),
            .groups = "drop")
cat("βNTI summary by phase-pair ",
    "(|βNTI|>2 = selection; |βNTI|≤2 = not selection-dominated):\n", sep = "")
print(bnti_summary, n = Inf)
w_bnti <- suppressWarnings(wilcox.test(
  bNTI ~ pair_type,
  data = bnti_df |> filter(pair_type %in% c("Before-Before","Before-D+A")),
  exact = FALSE))
cat(sprintf("Wilcoxon Before-Before vs Before-D+A: W = %.0f, p = %.4f\n",
            w_bnti$statistic, w_bnti$p.value))

cat("\n=== Stegen et al. 2013 process partition (βNTI + RC_Bray, 999 each) ===\n")
print(as.data.frame(proc_frac |> mutate(pct = round(pct, 1)) |>
  select(pair_type, process, n, n_pairs, pct)), row.names = FALSE)
cat("\nRC_Bray among |βNTI| <= 2 pairs (median, IQR):\n")
print(as.data.frame(bnti_df |> filter(abs(bNTI) <= 2) |> group_by(pair_type) |>
  summarize(n = n(), median = round(median(RC_bray), 3),
            q25 = round(quantile(RC_bray, 0.25), 3),
            q75 = round(quantile(RC_bray, 0.75), 3))), row.names = FALSE)

bnti_tsv <- file.path(fig_dir,
                      paste0("figure5_betaNTI_pairs_", date_tag, ".tsv"))
readr::write_tsv(bnti_df, bnti_tsv)
cat("βNTI pairs written to: ", bnti_tsv, "\n", sep = "")

# --- 7d) Per-genus Wilcoxon (top 20) --------------------------------
cat("\n=== Per-genus Wilcoxon: Before vs During+After (BH-adj) ===\n")
heat_with_phase <- heat_long |>
  rename(sample_id = Sample) |>
  mutate(sample_id = as.character(sample_id)) |>
  inner_join(sd_ar2 |> select(sample_id, ar2), by = "sample_id")

genus_tab <- heat_with_phase |>
  filter(Taxon != "Other") |>
  group_by(Taxon) |>
  summarize(
    mean_Before  = mean(Abundance[ar2 == "Before"]) * 100,
    mean_DurAft  = mean(Abundance[ar2 == "During+After"]) * 100,
    p_raw        = suppressWarnings(
      wilcox.test(Abundance ~ ar2, exact = FALSE)$p.value),
    .groups      = "drop") |>
  mutate(p_BH = p.adjust(p_raw, method = "BH")) |>
  arrange(p_BH)
print(genus_tab, n = Inf)

sink()
cat("Stats: ", stats_path, "\n", sep = "")
