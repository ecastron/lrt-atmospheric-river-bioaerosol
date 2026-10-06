#!/bin/bash
#SBATCH --job-name=coxL_search_01OCT26
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=02_search_%j.log

# 02_search_orfs.sh
# Searches all assembled air ORFs (single-cell and metagenome assemblies, 15 libraries) against the motif-labeled CoxL
# references plus outgroup, and against the accessory Cox proteins for gene-context calls (DIAMOND blastp).
# Set PROJ to the shotgun project directory and COXL_DIR to this folder (defaults: current directory).
# Run    : sbatch 02_search_orfs.sh

set -euo pipefail
cd "${COXL_DIR:-$PWD}"
ORF="${PROJ:-$PWD}/air_func/all_orfs.faa"
set +u; module purge; module load Diamond-2.1.23/diamond-2.1.23; set -u
diamond version
diamond makedb --in ref_coxL.faa -d ref_coxL --quiet
diamond makedb --in ref_cox_accessory.faa -d ref_cox_accessory --quiet
diamond blastp -d ref_coxL -q "$ORF" --sensitive -e 1e-10 --max-target-seqs 5 --threads "$SLURM_CPUS_PER_TASK" --quiet \
  -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen -o coxL_hits.tsv
diamond blastp -d ref_cox_accessory -q "$ORF" --sensitive -e 1e-10 --max-target-seqs 3 --threads "$SLURM_CPUS_PER_TASK" --quiet \
  -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen -o cox_accessory_hits.tsv
# extract candidate ORF sequences and all ORF headers (for contig order and length)
cut -f1 coxL_hits.tsv | sort -u > coxL_hit_ids.txt
awk 'NR==FNR{k[$1]=1;next} /^>/{id=substr($1,2); p=(id in k)} p' coxL_hit_ids.txt "$ORF" > coxL_hit_orfs.faa
grep "^>" "$ORF" | sed 's/^>//' | awk '{print $1"\t"$3"\t"$5"\t"$7}' > all_orf_coords.tsv
echo "coxL hit ORFs: $(grep -c '^>' coxL_hit_orfs.faa)"
echo done
