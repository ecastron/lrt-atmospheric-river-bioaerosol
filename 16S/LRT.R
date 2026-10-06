# LRT.R
# DADA2 processing of the LRT amplicon libraries: 16S rRNA gene V4 (first block) and 18S rRNA gene V9 (second block).
# Each block filters and trims reads, learns error rates, infers ASVs with pooled sample inference, merges pairs,
# removes chimeras, assigns taxonomy and builds a phyloseq object. Intended to be run interactively, block by block.
# Inputs : demultiplexed FASTQ files (NCBI SRA BioProject PRJNA1540148) in the folders set by miseq_path;
#          metadata_LRT.csv and metadata_LRT18S.csv; SILVA 138.1 and DADA2 EUK SSU v1.9 reference files; trees
#          (lrt.tre, 18S/18S_aln.tre) inferred from the written alignments.
# Outputs: 16S/ps_16S_LRT_withSeqs.RDS, 16S/ps_16S_LRT.RDS, 18S/ps_18S_LRT.RDS
# Run    : interactively in R from the 16S/ folder (set miseq_path to the FASTQ folder).

# Packages used in both blocks
library(dada2); library(DECIPHER); library(phangorn); library(phyloseq); library(Biostrings); library(ape)

##### 16S rRNA gene V4: DADA2 inference #####
set.seed(100)

miseq_path <- "../raw_data/16S/" # folder with the demultiplexed 16S FASTQ files
list.files(miseq_path)

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_L001_R1_001.fastq.gz"))
fnRs <- sort(list.files(miseq_path, pattern="_L001_R2_001.fastq.gz"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]

library(dada2)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF)
plotErrors(errR)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))

seqtabNoC <- removeBimeraDenovo(seqtabAll)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtabAll), rowSums(seqtabNoC))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sampleNames
head(track)

fastaRef <- "silva_nr99_v138.1_train_set.fa"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
taxTabExtra <- addSpecies(taxTab, "silva_species_assignment_v138.1.fa", verbose=TRUE)
unname(head(taxTabExtra))

View(taxTabExtra)

seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)


phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

##### Generate phyloseq object ####

samdf <- read.csv("metadata_LRT.csv", header=TRUE, row.names = 1, stringsAsFactors = T)
rownames(seqtabNoC) %in% rownames(samdf)
  all(rownames(seqtabAll) %in% samdf$sample_name)
library(ape)
tree<- read.tree(file = "lrt.tre")
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTabExtra),
               phy_tree(tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
saveRDS(ps, file = "ps_16S_LRT_withSeqs.RDS")
rank_names(ps)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, file = "ps_16S_LRT.RDS")
readRDS(file = "ps_16S_LRT.RDS") -> ps

##### DADA2 18S rRNA Inference #####
library(dada2)
set.seed(100)

setwd("../18S/data/")
miseq_path <- "./" # 18S/data/ holds the demultiplexed 18S FASTQ files
list.files(miseq_path)

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_L001_R1_001.fastq.gz"))
fnRs <- sort(list.files(miseq_path, pattern="_L001_R2_001.fastq.gz"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]

library(dada2)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF)
plotErrors(errR)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))

seqtabNoC <- removeBimeraDenovo(seqtabAll)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtabAll), rowSums(seqtabNoC))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sampleNames
head(track)

fastaRef <- "../../DADA2_EUK_SSU_v1.9.fasta"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

View(taxTab)

seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

library(Biostrings)
writeXStringSet(alignment, filepath = "../../18S/18S_aln.fa", format = "fasta")
tree<- read.tree(file = "../../18S/18S_aln.tre")

##### Generate phyloseq object ####
setwd("..")
samdf <- read.csv("../metadata_LRT18S.csv", header=TRUE, row.names = 1, stringsAsFactors = T)
rownames(seqtabNoC) %in% rownames(samdf)
all(rownames(seqtabAll) %in% samdf$sample_name)

ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab),
               phy_tree(tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

rank_names(ps)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, file = "ps_18S_LRT.RDS")
readRDS(file = "ps_18S_LRT.RDS") -> ps
