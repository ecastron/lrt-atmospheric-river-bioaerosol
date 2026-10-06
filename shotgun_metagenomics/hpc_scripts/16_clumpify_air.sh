#!/bin/bash
#SBATCH --job-name=clumpify_air
#SBATCH --partition=basic
#SBATCH --array=0-14%3
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --time=06:00:00
#SBATCH --output=logs/clumpify_%A_%a.out
#SBATCH --error=logs/clumpify_%A_%a.err

# 16_clumpify_air.sh
# Removes amplification duplicates from the trimmed MDA air libraries with clumpify (BBTools) before assembly (Roux
# et al. 2019). Reads are decompressed to node-local scratch first, and memory is set so clumpify runs in core.
# Set PROJECT to the shotgun project directory (default: current directory).
# Run    : sbatch 16_clumpify_air.sh

set -euo pipefail
shopt -s nullglob
PROJECT="${PROJECT:-$PWD}"
TRIM="$PROJECT/trimmed"
OUT="$PROJECT/air_dedup"
mkdir -p "$OUT" logs

R1S=( "$TRIM"/*PWA*_R1.trimmed.fastq.gz "$TRIM"/*RA[0-9]*_R1.trimmed.fastq.gz )
(( ${#R1S[@]} == 15 )) || { echo "ERROR: matched ${#R1S[@]} air R1, expected 15" >&2; exit 1; }
R1="${R1S[$SLURM_ARRAY_TASK_ID]}"
R2="${R1/_R1.trimmed.fastq.gz/_R2.trimmed.fastq.gz}"
S="$(basename "$R1" | sed -E 's/_R1\.trimmed\.fastq\.gz//')"
O1="$OUT/${S}_R1.dedup.fastq.gz"; O2="$OUT/${S}_R2.dedup.fastq.gz"

SCRATCH="/scratch/$USER/${SLURM_JOB_ID}_clump_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SCRATCH"; cleanup(){ rm -rf "$SCRATCH"; }; trap cleanup EXIT
echo "[$(date +%F_%T)] staging + clumpify dedupe $S on $(hostname)"
zcat "$R1" > "$SCRATCH/r1.fastq"; zcat "$R2" > "$SCRATCH/r2.fastq"

module purge 2>/dev/null; module load BBtools/bbtools-39.49 2>/dev/null
# dedupe=t removes duplicate read pairs; optical=f so ALL duplicates (PCR +
# optical) go, not just optically adjacent ones: the right choice for MDA.
clumpify.sh -Xmx72g in1="$SCRATCH/r1.fastq" in2="$SCRATCH/r2.fastq" \
  out1="$O1" out2="$O2" dedupe=t optical=f threads="$SLURM_CPUS_PER_TASK"

b=$(( $(wc -l < "$SCRATCH/r1.fastq") / 4 )); a=$(( $(zcat "$O1" | wc -l) / 4 ))
echo "[$(date +%F_%T)] $S reads: $b -> $a ($(awk -v a=$a -v b=$b 'BEGIN{printf "%.1f", 100*a/b}')% kept)"
