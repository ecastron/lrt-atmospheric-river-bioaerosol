# 07b_classify_v2.R
# Classifies each CoxL candidate from the rooted IQ-TREE tree, its active-site motif and its gene context. The four
# outgroup molybdenum hydroxylases are not monophyletic, so the tree is rooted on outgroup tip P80457. For each
# candidate the script records:
#   in_coxL   inside the MRCA clade of all form I and form II references (CoxL clade)
#   fI_sup    highest UFBoot of any ancestor clade whose references (>= 2) are all form I
#   fII_sup   same for form II
#   motif     residues aligned to P19919 AYRCSFR (columns from the MAFFT alignment)
# A gene is called form I CO dehydrogenase when it carries the form I motif (AYXCSFR) and fI_sup >= 95.
# Inputs : coxL_iqtree.contree, tree_input.aln, ref_coxL.tsv, candidates_dereplicated.tsv, candidates_context.tsv,
#          ../../air_orf_lca_taxonomy.tsv
# Output : candidates_classified_v2.tsv
# Run    : Rscript 07b_classify_v2.R (paths are set relative to the repository root)

suppressPackageStartupMessages({library(ape); library(phangorn); library(Biostrings); library(dplyr); library(readr); library(tidyr)})
setwd("./shotgun_metagenomics/data/curation_evidence/coxL_reanalysis_01OCT26")
aln <- readAAStringSet("tree_input.aln"); names(aln) <- sub(" .*$", "", names(aln))
p19 <- grep("_P19919$", names(aln), value = TRUE); ps <- as.character(aln[[p19]])
m <- regexpr("A-*Y-*R-*C-*S-*F-*R", ps); stopifnot(m > 0)
cols <- m:(m + attr(m, "match.length") - 1)
p19cols <- which(strsplit(ps, "")[[1]] != "-")
motif_of <- function(nm) gsub("-", "", substring(as.character(aln[[nm]]), min(cols), max(cols)))
cov_of <- function(nm) { s <- strsplit(as.character(aln[[nm]]), "")[[1]][p19cols]; r <- range(which(s != "-")); paste0(r[1], "-", r[2]) }

tr <- read.tree("coxL_iqtree.contree")
out <- grep("_OUT_", tr$tip.label, value = TRUE)
r <- root(tr, grep("P80457", out, value = TRUE), resolve.root = TRUE)
ntip <- Ntip(r); sup <- suppressWarnings(as.numeric(r$node.label))
refs_of <- function(n) { d <- r$tip.label[Descendants(r, n, "tips")[[1]]]; d[!grepl("^CAND_", d)] }
coxL_node <- getMRCA(r, grep("_FI+_", r$tip.label, value = TRUE))
fI_node <- getMRCA(r, grep("_FI_", r$tip.label, value = TRUE))
coxL_tips <- r$tip.label[Descendants(r, coxL_node, "tips")[[1]]]
cat("CoxL clade (all form I + II refs): UFBoot", sup[coxL_node - ntip], "| contains outgroup:", any(grepl("_OUT_", coxL_tips)), "\n")
cat("Form I refs monophyletic: UFBoot", sup[fI_node - ntip], "| FII refs inside it:", sum(grepl("_FII_", refs_of(fI_node))), "\n")

cand <- grep("^CAND_", r$tip.label, value = TRUE)
res <- lapply(cand, function(cn) {
  anc <- Ancestors(r, match(cn, r$tip.label), "all")
  pure <- function(tag) { s <- vapply(anc, function(a) { rf <- refs_of(a); if (length(rf) >= 2 && all(grepl(tag, rf))) sup[a - ntip] else NA_real_ }, 0)
                          if (all(is.na(s))) NA_real_ else max(s, na.rm = TRUE) }
  data.frame(cand = cn, in_coxL = cn %in% coxL_tips, fI_sup = pure("_FI_"), fII_sup = pure("_FII_"))
}) |> bind_rows()
res$motif <- vapply(res$cand, motif_of, ""); res$p19919_cov <- vapply(res$cand, cov_of, "")
res$motif_class <- case_when(grepl("^AY.CSFR$", res$motif) ~ "form I (AYXCSFR)",
                             grepl("^AY.GAGR$", res$motif) ~ "form II (AYXGAGR)",
                             grepl("^.Y.CSFR$", res$motif) ~ "CSFR core, non-canonical",
                             nchar(res$motif) < 7 ~ "not covered", TRUE ~ "other")
d <- cophenetic(r); reftips <- r$tip.label[!grepl("^CAND_", r$tip.label)]
res$nearest_ref <- vapply(res$cand, function(cn) names(which.min(d[cn, reftips])), "")
res <- res |> mutate(placement = case_when(!in_coxL ~ "outside CoxL clade",
                                           !is.na(fI_sup) & fI_sup >= 95 ~ "form I (>=95)",
                                           !is.na(fI_sup) ~ "form I (<95)",
                                           !is.na(fII_sup) & fII_sup >= 95 ~ "form II (>=95)",
                                           !is.na(fII_sup) ~ "form II (<95)", TRUE ~ "CoxL, unresolved form"),
                     verdict = case_when(motif_class == "form I (AYXCSFR)" & placement == "form I (>=95)" ~ "form I, CO-oxidation capacity",
                                         placement == "form I (>=95)" ~ "form I placement; motif not form I/not covered",
                                         placement == "form I (<95)" ~ "form I placement, weak support",
                                         motif_class == "form II (AYXGAGR)" | placement == "form II (>=95)" ~ "form II",
                                         placement == "outside CoxL clade" ~ "not CoxL (other Mo hydroxylase)",
                                         TRUE ~ "CoxL, form unresolved"))

# join dereplication, context and taxonomy (via the representative ORF; members share a gene)
der <- read_tsv("candidates_dereplicated.tsv", show_col_types = FALSE) |>
  mutate(cand = paste0("CAND_", library, "_", sub("^.*cov_[0-9]+\\.", "", rep)))
stopifnot(all(res$cand %in% der$cand))
ctx <- read_tsv("candidates_context.tsv", show_col_types = FALSE) |> select(rep, contig_len, assessable, coxM, coxS, coxMSL_order, coxG, coxD, coxE, coxF)
tax <- read_tsv("../../air_orf_lca_taxonomy.tsv", show_col_types = FALSE) |> select(orf_id, supported_name, lineage)
taxm <- der |> select(cand, members) |> separate_rows(members, sep = ";") |> left_join(tax, by = c(members = "orf_id")) |>
  group_by(cand) |> summarize(taxonomy = paste(unique(na.omit(supported_name)), collapse = "|"),
                              lineage = paste(unique(na.omit(lineage)), collapse = "|"), .groups = "drop")
qc_pass <- c("PWA2", "RA3", "RA8", "RA14", "RA15", "RA16", "RA17", "RA18", "RA22", "RA23", "RA24")
tab <- res |> left_join(der |> select(cand, rep, library, n_members, len_aa, best_ref, best_pid), by = "cand") |>
  left_join(ctx, by = "rep") |> left_join(taxm, by = "cand") |>
  mutate(qc_pass = library %in% qc_pass,
         context = ifelse(assessable, paste0("contig ", contig_len, " bp; coxM ", coxM, ", coxS ", coxS, ", MSL order ", coxMSL_order,
                                              ", coxG ", coxG, ", coxD/E/F ", coxD | coxE | coxF), paste0("not assessable (contig ", contig_len, " bp)"))) |>
  arrange(factor(verdict, c("form I, CO-oxidation capacity", "form I placement; motif not form I/not covered", "form I placement, weak support",
                            "form II", "CoxL, form unresolved", "not CoxL (other Mo hydroxylase)")), library)
write_tsv(tab, "candidates_classified_v2.tsv")
print(table(tab$verdict)); print(table(tab$library, tab$verdict))
print(as.data.frame(tab |> filter(grepl("form I", verdict)) |> select(cand, len_aa, best_pid, p19919_cov, motif, fI_sup, nearest_ref, context, taxonomy)), right = FALSE)
