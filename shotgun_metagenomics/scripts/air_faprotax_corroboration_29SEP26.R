#!/usr/bin/env Rscript
# air_faprotax_corroboration_29SEP26.R
# For each function inferred from 16S taxonomy by FAPROTAX, tests whether the encoding marker gene is directly
# detected in the air metagenomes. Only libraries passing the uvrA single-copy check (0.6-1.6 copies per genome after
# recA normalization) are counted, and a function is detected when its marker genes recruit at least MIN_READS reads.
# Also tabulates functions FAPROTAX does not model (DMS/DMSP organosulfur cycling, methylated-amine oxidation).
# Inputs : shotgun_metagenomics/data/air_func_local/ (read recruitment), shotgun_metagenomics/data/faprotax/,
#          shotgun_metagenomics/data/curation_evidence/, shotgun_metagenomics/data/atmos_marker_curation.tsv
# Outputs: shotgun_metagenomics/tables/air_faprotax_corroboration_29SEP26.csv, air_shotgun_only_dms_29SEP26.csv,
#          air_shotgun_only_functions_29SEP26.csv, air_uvrA_gate_29SEP26.csv
# Run    : Rscript shotgun_metagenomics/scripts/air_faprotax_corroboration_29SEP26.R

suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(purrr);library(forcats)})

ROOT <- "./shotgun_metagenomics"
L    <- file.path(ROOT,"data","air_func_local")

# ---- shotgun air marker reads (presence) ----
rd <- function(dir){
  do.call(rbind,lapply(list.files(dir,"markercounts.tsv$",full.names=TRUE),function(f){
    ln<-readLines(f); samp<-strsplit(ln[1],"\t")[[1]][2]
    b<-read.table(f,sep="\t",comment.char="#"); names(b)<-c("marker","reads")
    b$short<-str_extract(samp,"PWA\\d+|RA\\d+"); b}))}
air <- rd(file.path(L,"reads_recruit"))
MIN_READS <- 5

# ---- uvrA single-copy gate (identical to air_category_enrichment_curated.R) ----
len  <- read.table(file.path(L,"marker_lengths.tsv"),sep="\t",header=TRUE)
rk   <- air|>left_join(len|>select(marker,aa=median_aa_len),by="marker")|>
  filter(!is.na(aa))|>mutate(rk=reads/aa)
recA <- rk|>filter(marker=="recA")|>group_by(short)|>summarize(recA_rk=sum(rk),.groups="drop")
gate <- rk|>filter(marker=="uvrA")|>group_by(short)|>summarize(u=sum(rk),.groups="drop")|>
  right_join(tibble(short=sort(unique(air$short))),by="short")|>
  left_join(recA,by="short")|>
  mutate(uvrA_cpg=u/recA_rk, pass=!is.na(uvrA_cpg) & recA_rk>0 & uvrA_cpg>=0.6 & uvrA_cpg<=1.6)
cat("=== uvrA gate ===\n"); print(as.data.frame(gate|>mutate(uvrA_cpg=round(uvrA_cpg,2))))
airids <- sort(gate$short[gate$pass])                   # QC-passing air libs
stopifnot(length(airids) == 11)

# ---- FAPROTAX inferred (16S), restricted to the matched air libs ----
fap <- read.table(file.path(ROOT,"data","faprotax","faprotax_output_norm.tsv"),
                  sep="\t",header=TRUE,check.names=FALSE); names(fap)[1]<-"Function"
fap_air <- fap|>select(Function,any_of(airids))|>
  pivot_longer(-Function,names_to="short",values_to="fap")

# ---- HydDB-resolved group-1h/1l (atmospheric high-affinity) [NiFe] presence per sample ----
nife <- read.table(file.path(ROOT,"data","curation_evidence","nife_group_calls.tsv"),sep="\t",header=TRUE)
nife <- nife|>mutate(grp=str_replace(group,"\\[NiFe\\]_Group_",""),short=str_extract(orf,"PWA\\d+|RA\\d+"))
g1l_samples <- nife|>filter(grp %in% c("1h","1l"))|>distinct(short)|>pull(short)

# ---- function group -> (FAPROTAX functions, validated shotgun markers) ----
GR <- list(
 list(grp="Nitrogen fixation", funcs=c("nitrogen_fixation"),
      markers=c("NifH","Nitrogenase_NifH_99")),
 list(grp="Nitrification / NH3 ox", funcs=c("nitrification","aerobic_ammonia_oxidation","aerobic_nitrite_oxidation"),
      markers=c("Ammonia_monooxygenase_amoA","amoA","AmoB","AmoC","Nitrite_oxidoreductase_NxrA")),
 list(grp="Denitrification / NO3 red", funcs=c("denitrification","nitrate_reduction","nitrate_respiration","nitrate_denitrification","nitrite_denitrification","nitrous_oxide_denitrification"),
      markers=c("Dissimilatory_nitrate_reductase_NarG","NarG","Periplasmic_nitrate_reductase_NapA","NapA","Copper_containing_nitrite_reductase_NirK","Cytochrome_cd1_nitrite_reductase_NirS_sequence","Nitric_oxide_reductase_NorB","Nitrous_oxide_reductase_NosZ")),
 list(grp="Sulfur oxidation", funcs=c("dark_sulfide_oxidation","dark_sulfur_oxidation","dark_sulfite_oxidation","dark_thiosulfate_oxidation","dark_oxidation_of_sulfur_compounds","respiration_of_sulfur_compounds"),
      markers=c("Thiosulfohydrolase_SoxB","Sulfide_quinone_oxidoreductase_Sqr","Flavocytochrome_c_sulfide_dehydrogenase_FCC","Dissimilatory_sulfite_reductase_DsrA")),
 list(grp="Methanol ox / methylotrophy", funcs=c("methanol_oxidation","methylotrophy"),
      markers=c("xoxF","mxaF","mauA")),
 list(grp="Methanotrophy", funcs=c("methanotrophy"),
      markers=c("Particulate_methane_monooxygenase_PmoA","PmoA","Soluble_methane_monooxygenase_MmoA","MmoA")),
 list(grp="Dark H2 oxidation (HydDB grp 1l)", funcs=c("dark_hydrogen_oxidation","knallgas_bacteria"),
      markers="__HYDDB_1l__")   # HydDB-resolved group-1l atmospheric high-affinity (generic NiFe dropped)
)

n_air <- length(airids)
tab <- map_dfr(GR, function(g){
  inf <- fap_air|>filter(Function %in% g$funcs)|>group_by(short)|>summarize(v=sum(fap),.groups="drop")
  inferred_n <- sum(inf$v>0)
  if(identical(g$markers,"__HYDDB_1l__")){
    detected_n <- sum(airids %in% g1l_samples); reads<-NA
    note<-paste0("HydDB group-1l ORFs in {",paste(sort(intersect(airids,g1l_samples)),collapse=","),"}; atmospheric high-affinity")
  } else if(length(g$markers)==0){
    det <- tibble(short=airids,r=0); detected_n<-NA; reads<-NA; note<-"no specific marker (generic NiFe dropped)"
  } else {
    det <- air|>filter(marker %in% g$markers)|>group_by(short)|>summarize(r=sum(reads),.groups="drop")
    det <- tibble(short=airids)|>left_join(det,by="short")|>mutate(r=coalesce(r,0))
    detected_n<-sum(det$r>=MIN_READS); reads<-sum(det$r); note<-""
  }
  tibble(group=g$grp, fap_inferred_n=inferred_n, shotgun_detected_n=detected_n,
         shotgun_total_reads=reads, n_air=n_air,
         corroborated=ifelse(is.na(detected_n),"no marker",ifelse(detected_n>0,"YES","no")),
         note=note)
})
cat("=== FAPROTAX inferred vs shotgun gene detected (>= MIN_READS reads; uvrA-passing air libraries) ===\n")
print(as.data.frame(tab))

# ---- shotgun-ONLY functions FAPROTAX does not model (organosulfur + methylated amines) ----
# Pull markers straight from the curation table (decision==keep) so the figure
# stays in sync with curation and silently drops nothing. Two real atmospheric
# functions with NO FAPROTAX category -> can only be seen by the shotgun.
cur <- read.table(file.path(ROOT,"data","atmos_marker_curation.tsv"),sep="\t",header=TRUE,quote="")
SO_LAB <- c(DMS_DMSP_organosulfur = "DMS/DMSP cycling (organosulfur)",
            Methylated_amines     = "Methylated-amine oxidation")
so_markers <- cur|>filter(decision=="keep", subcategory %in% names(SO_LAB))|>
  transmute(marker, group=unname(SO_LAB[subcategory]))

# per-marker detection (supplementary table, as before)
dms <- air|>filter(short %in% airids, marker %in% so_markers$marker)|>group_by(marker)|>
  summarize(detected_in=sum(reads>=MIN_READS),total_reads=sum(reads),.groups="drop")|>arrange(desc(total_reads))
cat("\n=== shotgun-ONLY: DMS/DMSP + methylated-amine genes (no FAPROTAX function exists) ===\n")
print(as.data.frame(dms))

# per-function presence across the 15 air libraries (union of the group's markers);
# fap_inferred_n = NA because FAPROTAX cannot infer these at all
so_tab <- so_markers|>inner_join(air|>filter(short %in% airids),by="marker",relationship="many-to-many")|>
  group_by(group,short)|>summarize(r=sum(reads),.groups="drop")|>
  group_by(group)|>summarize(shotgun_detected_n=sum(r>=MIN_READS), shotgun_total_reads=sum(r),.groups="drop")|>
  mutate(fap_inferred_n=NA_integer_, n_air=n_air, corroborated="shotgun-only",
         note="no FAPROTAX function exists")

# CO oxidation is not scored here: CoxL genes are classified separately by active-site motif and phylogeny
# (shotgun_metagenomics/data/curation_evidence/coxL_reanalysis_01OCT26/), and no air gene met the form I criteria.
cat("\n=== shotgun-ONLY functions promoted to the figure ===\n")
print(as.data.frame(so_tab|>select(group,shotgun_detected_n,note)))

write.csv(tab,   file.path(ROOT,"tables","air_faprotax_corroboration_29SEP26.csv"),row.names=FALSE)
write.csv(dms,   file.path(ROOT,"tables","air_shotgun_only_dms_29SEP26.csv"),row.names=FALSE)
write.csv(so_tab,file.path(ROOT,"tables","air_shotgun_only_functions_29SEP26.csv"),row.names=FALSE)

write.csv(gate, file.path(ROOT,"tables","air_uvrA_gate_29SEP26.csv"), row.names=FALSE)
cat("\nwrote tables/*_29SEP26.csv\n")
