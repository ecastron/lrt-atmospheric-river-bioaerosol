#!/usr/bin/env Rscript
# air_nife_hyddb_29SEP26.R
# Supplementary Figure 4. Number of [NiFe]-hydrogenase large-subunit ORF calls from the air assemblies per HydDB group
# (DIAMOND best hit against group-labeled HydDB references), summed over the single-cell and metagenome assemblies.
# Input  : shotgun_metagenomics/data/curation_evidence/nife_group_calls.tsv, 16S/palettes.R
# Output : shotgun_metagenomics/figures/air_nife_hyddb_29SEP26.{pdf,png}
# Run    : Rscript shotgun_metagenomics/scripts/air_nife_hyddb_29SEP26.R

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(stringr); library(forcats)})
LRT  <- "."
ROOT <- file.path(LRT, "shotgun_metagenomics")
source(file.path(LRT, "16S", "palettes.R"))

d <- read.table(file.path(ROOT, "data", "curation_evidence", "nife_group_calls.tsv"), sep = "\t", header = TRUE) |>
  mutate(group = str_replace(group, "\\[NiFe\\]_Group_", ""))
stopifnot(nrow(d) == 79)

gd <- d |> count(group) |>
  mutate(class = case_when(group %in% c("1h", "1l") ~ "Group 1l (atmospheric, high-affinity)",
                           str_starts(group, "4")   ~ "Group 4 (energy-converting, fermentative)",
                           str_starts(group, "3")   ~ "Group 3 (bidirectional, cofactor-coupled)",
                           TRUE                     ~ "Groups 1a-1g (uptake and other)"),
         class = factor(class, levels = names(pal_hyddb)),
         group = fct_reorder(group, n))
stopifnot(!any(gd$group == "1h"))   # legend states no group 1h ORFs

p <- ggplot(gd, aes(n, group, fill = class)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = n), hjust = -0.3, size = 2.6) +
  scale_fill_manual(values = pal_hyddb, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Air [NiFe]-hydrogenase ORFs (n = 79)", y = "HydDB group", fill = NULL) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(legend.position = "right", legend.text = element_text(size = 7),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(ROOT, "figures", "air_nife_hyddb_29SEP26.pdf"), p, width = 183, height = 80, units = "mm", device = cairo_pdf)
ggsave(file.path(ROOT, "figures", "air_nife_hyddb_29SEP26.png"), p, width = 183, height = 80, units = "mm", dpi = 300)
print(as.data.frame(gd))
