#!/bin/bash
#SBATCH --job-name=air_mplca
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/mplca_%j.out
#SBATCH --error=logs/mplca_%j.err

# 55_air_mplca.sh
# Taxonomic assignment of marker-hit ORFs by minimum-support lowest common ancestor (MetaPathways v3.5
# LCAComputation; bitscore >= 50, top 10% of hits, support >= 5) from DIAMOND blastp hits against NCBI nr.
# Set PROJ to the shotgun project directory (default: current directory).
# Run    : sbatch 55_air_mplca.sh

set -euo pipefail
PROJ="${PROJ:-$PWD}"
AF="$PROJ/air_func"

module load MetaPathways-LCA/3.5
echo "[$(date +%F_%T)] tree=$MPLCA_TREE  node=$(hostname)"

mp_lca_run.py \
  --hits "$AF/nr_marker_hits.tsv" \
  --out  "$AF/orf_lca_taxonomy.tsv" \
  --min-score 50 --top-percent 10 --min-support 5

echo "[$(date +%F_%T)] done -> $AF/orf_lca_taxonomy.tsv ($(wc -l < "$AF/orf_lca_taxonomy.tsv") ORFs)"
