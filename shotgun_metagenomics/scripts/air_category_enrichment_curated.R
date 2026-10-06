#!/usr/bin/env Rscript
# air_category_enrichment_curated.R
# Compares atmosphere-associated marker genes between the air metagenomes and Byers Peninsula bulk soil. Read counts
# per marker are normalized by median protein length and by recA (single copy) to give copies per genome; libraries
# are kept when the single-copy control uvrA falls between 0.6 and 1.6 copies per genome. Only atmosphere-specific,
# family-coherent markers are used (generic CoxL, [NiFe]/[FeFe]-hydrogenase, MDH and FdhA buckets are excluded; see
# data/atmos_marker_curation.tsv). Air DNA was MDA-amplified and soil DNA was not, so only large differences are
# interpreted.
# Inputs : shotgun_metagenomics/data/air_func_local/ (read recruitment counts, marker lengths),
#          shotgun_metagenomics/data/atmos_marker_curation.tsv
# Outputs: shotgun_metagenomics/tables/air_curated_category_persample.csv, air_curated_subcategory.csv,
#          shotgun_metagenomics/figures/air_curated_enrichment.{png,pdf}
# Run    : Rscript shotgun_metagenomics/scripts/air_category_enrichment_curated.R

suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(forcats)})

ROOT <- "./shotgun_metagenomics"
L    <- file.path(ROOT,"data","air_func_local")
cur  <- read.table(file.path(ROOT,"data","atmos_marker_curation.tsv"),sep="\t",header=TRUE,quote="")
keep <- cur |> filter(decision=="keep")
len  <- read.table(file.path(L,"marker_lengths.tsv"),sep="\t",header=TRUE)

rd <- function(dir,matrix){
  do.call(rbind,lapply(list.files(dir,"markercounts.tsv$",full.names=TRUE),function(f){
    ln<-readLines(f); samp<-strsplit(ln[1],"\t")[[1]][2]
    b<-read.table(f,sep="\t",comment.char="#"); names(b)<-c("marker","reads")
    b$sample<-samp; b$matrix<-matrix; b}))}
d <- rbind(rd(file.path(L,"reads_recruit"),"Air"), rd(file.path(L,"byers_recruit"),"Soil")) |>
  left_join(len|>select(marker,aa=median_aa_len),by="marker") |> filter(!is.na(aa)) |> mutate(rk=reads/aa)

recA <- d|>filter(marker=="recA")|>group_by(sample)|>summarize(recA_rk=sum(rk),.groups="drop")
ctrl <- d|>filter(marker=="uvrA")|>group_by(sample,matrix)|>summarize(u=sum(rk),.groups="drop")|>
  left_join(recA,by="sample")|>mutate(uvrA_cpg=u/recA_rk)
ok <- ctrl|>filter(recA_rk>0,uvrA_cpg>=0.6,uvrA_cpg<=1.6)|>pull(sample)

dk <- d|>filter(sample %in% ok, marker %in% keep$marker)|>left_join(recA,by="sample")|>
  mutate(cpg=rk/recA_rk)|>left_join(keep|>select(marker,subcategory),by="marker")

# ---- whole curated-Atmospheric category, air vs soil ----
tot <- dk|>group_by(matrix,sample)|>summarize(cpg=sum(cpg),.groups="drop")
ptot <- wilcox.test(cpg~matrix,data=tot)$p.value
cat("=== CURATED atmosphere-specific category: copies/genome (mean over samples) ===\n")
cat(sprintf("  Air  = %.2f   Soil = %.2f   fold = %.2f   Wilcoxon p = %.3g\n",
            mean(tot$cpg[tot$matrix=="Air"]), mean(tot$cpg[tot$matrix=="Soil"]),
            mean(tot$cpg[tot$matrix=="Air"])/mean(tot$cpg[tot$matrix=="Soil"]), ptot))

# ---- per subcategory ----
sub <- dk|>group_by(matrix,sample,subcategory)|>summarize(cpg=sum(cpg),.groups="drop")
cat("\n=== per subcategory (mean cpg, fold, Wilcoxon p) ===\n")
subsum <- sub|>group_by(subcategory,matrix)|>summarize(mean=mean(cpg),.groups="drop")|>
  pivot_wider(names_from=matrix,values_from=mean)|>mutate(fold=round(Air/Soil,2))
for(sc in subsum$subcategory){s<-sub|>filter(subcategory==sc)
  p<-tryCatch(wilcox.test(cpg~matrix,data=s)$p.value,error=function(e)NA)
  subsum$p[subsum$subcategory==sc]<-p}
print(subsum|>mutate(Air=round(Air,3),Soil=round(Soil,3),p=signif(p,2))|>as.data.frame())

write.csv(tot,file.path(ROOT,"tables","air_curated_category_persample.csv"),row.names=FALSE)
write.csv(subsum,file.path(ROOT,"tables","air_curated_subcategory.csv"),row.names=FALSE)

# ============================ FIGURE ============================
sub2 <- sub|>mutate(subcategory=fct_reorder(subcategory,cpg,.fun=median,.desc=TRUE))
labp <- subsum|>mutate(lab=sprintf("%s\nfold %.2f, p=%.2g",subcategory,fold,p))
sub2 <- sub2|>left_join(labp|>select(subcategory,lab,fold),by="subcategory")|>
  mutate(lab=fct_reorder(lab,fold,.desc=TRUE))
p <- ggplot(sub2,aes(matrix,cpg,fill=matrix))+
  geom_boxplot(outlier.size=.5,alpha=.85,width=.6)+geom_jitter(width=.12,size=.9,alpha=.5)+
  facet_wrap(~lab,scales="free_y",nrow=2)+
  scale_fill_manual(values=c("Air"="#D55E00","Soil"="#009E73"))+
  labs(x=NULL,y="copies per genome (recA-normalized)",fill=NULL,
       title="Recurated atmosphere-specific gene enrichment: Antarctic air vs Byers soil",
       subtitle=sprintf("Generic enzyme-family buckets (CoxL/NiFe/FeFe/MDH/FdhA) excluded. Curated category: Air %.1f vs Soil %.1f, fold %.1f, p=%.2g",
                        mean(tot$cpg[tot$matrix=="Air"]),mean(tot$cpg[tot$matrix=="Soil"]),
                        mean(tot$cpg[tot$matrix=="Air"])/mean(tot$cpg[tot$matrix=="Soil"]),ptot),
       caption="recA-normalized; uvrA~1 gated. air MDA vs soil non-MDA -> large effects. See data/atmos_marker_curation.tsv.")+
  theme_minimal(base_size=11,base_family="Helvetica")+
  theme(legend.position="top",strip.text=element_text(face="bold",size=8.5),
        plot.title=element_text(face="bold",size=13))
ggsave(file.path(ROOT,"figures","air_curated_enrichment.png"),p,width=11,height=6.5,dpi=300)
ggsave(file.path(ROOT,"figures","air_curated_enrichment.pdf"),p,width=11,height=6.5)
cat("\nwrote figures/air_curated_enrichment.png/.pdf + tables/air_curated_*.csv\n")
