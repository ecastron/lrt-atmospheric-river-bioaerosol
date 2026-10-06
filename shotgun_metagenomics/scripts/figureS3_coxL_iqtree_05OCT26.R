# figureS3_coxL_iqtree_05OCT26.R
# Supplementary Figure 3. Maximum-likelihood tree of the 89 dereplicated CoxL-like genes from the air assemblies with
# 67 references (33 form I, motif AYXCSFR; 30 form II, motif AYXGAGR; 4 molybdenum hydroxylase outgroups); MAFFT v7.490,
# IQ-TREE 3.0.1 (LG+G4), 1,000 ultrafast bootstraps. The four outgroups are not monophyletic, so the tree is rooted on
# outgroup tip P80457, as in 07b_classify_v2.R. Tip classes come from candidates_classified_v2.tsv and ref_coxL.tsv.
# Inputs : shotgun_metagenomics/data/curation_evidence/coxL_reanalysis_01OCT26/{coxL_iqtree.contree,ref_coxL.tsv,
#          candidates_classified_v2.tsv}, 16S/palettes.R
# Output : shotgun_metagenomics/figures/figureS3_coxL_iqtree_05OCT26.{pdf,png}
# Run    : Rscript shotgun_metagenomics/scripts/figureS3_coxL_iqtree_05OCT26.R

suppressPackageStartupMessages({library(ape); library(ggtree); library(ggplot2); library(dplyr)})
root <- "."
d    <- file.path(root, "shotgun_metagenomics/data/curation_evidence/coxL_reanalysis_01OCT26")
source(file.path(root, "16S", "palettes.R"))

tr   <- read.tree(file.path(d, "coxL_iqtree.contree"))
ref  <- read.delim(file.path(d, "ref_coxL.tsv"))
cand <- read.delim(file.path(d, "candidates_classified_v2.tsv"))
stopifnot(Ntip(tr) == nrow(ref) + nrow(cand), all(tr$tip.label %in% c(ref$id, cand$cand)))

out_ids <- ref$id[ref$form == "outgroup"]
# outgroups are not monophyletic; root on the P80457 outgroup tip as in 07b_classify_v2.R
tr <- root(tr, grep("P80457", out_ids, value = TRUE), resolve.root = TRUE)
tr <- ladderize(tr)

cls <- bind_rows(
  ref |> transmute(label = id,
                   class = case_when(form == "I" ~ "Form I reference (AYXCSFR)",
                                     form == "II" ~ "Form II reference (AYXGAGR)",
                                     TRUE ~ "Mo-hydroxylase outgroup"),
                   tip = ifelse(set == "named_formI", sapply(strsplit(organism, " "), function(w) paste(w[1:2], collapse = " ")), NA)),
  cand |> transmute(label = cand,
                    class = case_when(verdict == "form I placement; motif not form I/not covered" ~ "Air gene, form I clade, motif not confirmed",
                                      verdict == "form I placement, weak support" ~ "Air gene, form I clade, weak support",
                                      verdict == "form II" ~ "Air gene, form II",
                                      TRUE ~ "Air gene, not CoxL (other Mo-hydroxylase)"),
                    tip = ifelse(grepl("^form I placement", verdict),   # label only form I placements (not form II)
                                 paste0(library, " (", ifelse(motif == "" | is.na(motif), "motif not covered", motif), ")"), NA)))
stopifnot(all(cls$class %in% names(pal_coxL)))
cls$ff <- ifelse(grepl("reference", cls$class), "italic", "plain")   # species names in italics

p <- ggtree(tr, linewidth = 0.25) %<+% cls +
  geom_nodepoint(aes(subset = !isTip & suppressWarnings(as.numeric(label)) >= 95), size = 0.6, color = "gray30") +
  geom_tippoint(aes(color = class), size = 1.1) +
  geom_tiplab(aes(label = tip, fontface = ff), size = 1.9, offset = 0.02, na.rm = TRUE) +
  scale_color_manual(values = pal_coxL, name = NULL) +
  geom_treescale(width = 0.5, fontsize = 2, linesize = 0.3) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5)) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
  xlim(0, max(node.depth.edgelength(tr)) * 1.35)

out <- file.path(root, "shotgun_metagenomics/figures/figureS3_coxL_iqtree_05OCT26")
ggsave(paste0(out, ".pdf"), p, width = 183, height = 230, units = "mm")
ggsave(paste0(out, ".png"), p, width = 183, height = 230, units = "mm", dpi = 300)
cat("tips by class:\n"); print(table(cls$class))
