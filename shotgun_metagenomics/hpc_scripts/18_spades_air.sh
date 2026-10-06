#!/bin/bash
#SBATCH --job-name=spades_air
#SBATCH --partition=basic
#SBATCH --array=0-29
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/spades_%A_%a.out
#SBATCH --error=logs/spades_%A_%a.err

# 18_spades_air.sh
# Assembles each deduplicated air library twice: SPAdes --sc (single-cell mode, robust to MDA coverage bias) and
# --meta (metaSPAdes). Array of 30 tasks: 0-14 single-cell, 15-29 metagenome; sample = task index modulo 15.
# Set PROJECT to the shotgun project directory (default: current directory).
# Run    : sbatch 18_spades_air.sh

set -uo pipefail
shopt -s nullglob
PROJECT="${PROJECT:-$PWD}"
DEDUP="$PROJECT/air_dedup"
ASM="$PROJECT/assembly_air"
mkdir -p "$ASM" logs

R1S=( "$DEDUP"/*PWA*_R1.dedup.fastq.gz "$DEDUP"/*RA[0-9]*_R1.dedup.fastq.gz )
(( ${#R1S[@]} == 15 )) || { echo "ERROR: matched ${#R1S[@]} deduped R1, expected 15" >&2; exit 1; }

idx=$SLURM_ARRAY_TASK_ID
s=$(( idx % 15 ))
if (( idx < 15 )); then MODE="sc"; FLAG="--sc"; else MODE="meta"; FLAG="--meta"; fi
R1="${R1S[$s]}"; R2="${R1/_R1.dedup.fastq.gz/_R2.dedup.fastq.gz}"
S="$(basename "$R1" | sed -E 's/_R1\.dedup\.fastq\.gz//')"
OUTDIR="$ASM/${MODE}/${S}"
mkdir -p "$(dirname "$OUTDIR")"

module purge 2>/dev/null; module load SPAdes/SPAdes-4.2.0 2>/dev/null
echo "[$(date +%F_%T)] SPAdes $FLAG  sample=$S  on $(hostname)"

# node-local NVMe scratch for SPAdes' heavy temp IO; copy back only contigs.
SCRATCH="/scratch/$USER/${SLURM_JOB_ID}_${MODE}_${S}"
mkdir -p "$SCRATCH"; cleanup(){ rm -rf "$SCRATCH"; }; trap cleanup EXIT
MEM_GB=$(( SLURM_MEM_PER_NODE / 1024 ))

spades.py $FLAG -1 "$R1" -2 "$R2" \
  -t "$SLURM_CPUS_PER_TASK" -m "$MEM_GB" \
  --tmp-dir "$SCRATCH/tmp" -o "$SCRATCH/out"
rc=$?
mkdir -p "$OUTDIR"
cp "$SCRATCH/out/spades.log" "$OUTDIR/${S}.${MODE}.spades.log" 2>/dev/null || true   # keep log even on failure
if [[ $rc -ne 0 || ! -s "$SCRATCH/out/contigs.fasta" ]]; then
  echo "[$(date +%F_%T)] SPAdes FAILED ($MODE $S) rc=$rc: log kept at $OUTDIR/${S}.${MODE}.spades.log" >&2; exit 1
fi
cp "$SCRATCH/out/contigs.fasta" "$OUTDIR/${S}.${MODE}.contigs.fasta"
cp "$SCRATCH/out/scaffolds.fasta" "$OUTDIR/${S}.${MODE}.scaffolds.fasta" 2>/dev/null || true
cp "$SCRATCH/out/spades.log" "$OUTDIR/" 2>/dev/null || true
echo "[$(date +%F_%T)] done -> $OUTDIR/${S}.${MODE}.contigs.fasta"
