#!/bin/bash
#SBATCH --job-name=prodigal_air
#SBATCH --partition=basic
#SBATCH --array=0-29%6
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/prodigal_%A_%a.out
#SBATCH --error=logs/prodigal_%A_%a.err

# 20_prodigal_air.sh
# Predicts ORFs with Prodigal (metagenome mode) on the 30 air assemblies after removing contigs shorter than 500 bp.
# Set PROJ to the shotgun project directory (default: current directory).
# Run    : sbatch 20_prodigal_air.sh

set -euo pipefail
shopt -s nullglob
PROJ="${PROJ:-$PWD}"
OUT="$PROJ/air_func/orfs"; mkdir -p "$OUT" logs

FILES=( "$PROJ"/assembly_air/sc/*/*.sc.contigs.fasta "$PROJ"/assembly_air/meta/*/*.meta.contigs.fasta )
(( ${#FILES[@]} == 30 )) || { echo "ERROR: ${#FILES[@]} contig files, expected 30" >&2; exit 1; }
F="${FILES[$SLURM_ARRAY_TASK_ID]}"
name="$(basename "$F" .contigs.fasta)"           # e.g. 006_6_RA8_GMCF3029a.meta

SCRATCH="/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_prod"
mkdir -p "$SCRATCH"; cleanup(){ rm -rf "$SCRATCH"; }; trap cleanup EXIT

module purge 2>/dev/null; module load BBtools/bbtools-39.49 2>/dev/null
reformat.sh in="$F" out="$SCRATCH/filt.fasta" minlength=500 ow=t 2>"$SCRATCH/reformat.log"
nkept=$(grep -c "^>" "$SCRATCH/filt.fasta" 2>/dev/null || echo 0)

module purge 2>/dev/null; module load Prodigal/prodigal-2.6.3 2>/dev/null
echo "[$(date +%F_%T)] prodigal $name ($nkept contigs >=500bp) on $(hostname)"
prodigal -p meta -q -i "$SCRATCH/filt.fasta" \
  -a "$OUT/${name}.faa" -d "$OUT/${name}.fna" -f gff -o "$OUT/${name}.gff"
# tag ORF headers with assembly id so they stay unique across samples in the pool
sed -i "s/^>/>${name}__/" "$OUT/${name}.faa" "$OUT/${name}.fna"
echo "[$(date +%F_%T)] $name: $(grep -c '^>' "$OUT/${name}.faa") ORFs"
