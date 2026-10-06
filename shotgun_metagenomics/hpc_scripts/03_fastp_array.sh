#!/bin/bash
#SBATCH --job-name=fastp_GMCF3029a
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --array=0-74%8
#SBATCH --output=logs/fastp_%A_%a.out
#SBATCH --error=logs/fastp_%A_%a.err

# 03_fastp_array.sh
# Adapter and quality trimming with fastp (one paired-end library per array task). Air libraries (PWA*, RA*) were
# MDA-amplified and get an additional 5' hard trim (--trim_front 12) to remove random-hexamer priming bias; other
# libraries get the standard trim. fastp detects paired-end adapters and enables poly-G trimming for two-color
# chemistry; quality filtering uses fastp defaults.
# Set PROJECT to the shotgun project directory (default: current directory) and FASTQ_ROOT to the raw read folder.
# Run    : sbatch 03_fastp_array.sh

set -euo pipefail

PROJECT="${PROJECT:-$PWD}"
FASTQ_ROOT="${FASTQ_ROOT:-$PROJECT/raw_reads}"
TRIMDIR="$PROJECT/trimmed"
REPDIR="$PROJECT/qc/fastp"
mkdir -p "$TRIMDIR" "$REPDIR" logs

# --- Resolve this task's sample ---------------------------------------------
mapfile -t R1_FILES < <(find "$FASTQ_ROOT" -name '*_R1_001.fastq.gz' | sort)
if [[ ${#R1_FILES[@]} -ne 75 ]]; then
  echo "ERROR: expected 75 R1 files, found ${#R1_FILES[@]}" >&2; exit 1
fi
R1="${R1_FILES[$SLURM_ARRAY_TASK_ID]}"
R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
SAMPLE="$(basename "$R1" | sed -E 's/_S[0-9]+_L[0-9]+_R1_001\.fastq\.gz//')"
[[ -s "$R1" ]] || { echo "ERROR: missing/empty R1: $R1" >&2; exit 1; }
[[ -s "$R2" ]] || { echo "ERROR: missing/empty R2: $R2" >&2; exit 1; }

# --- Pick the profile from the library-type token (3rd underscore field) -----
TYPE="$(echo "$SAMPLE" | cut -d_ -f3)"
EXTRA=()
if [[ "$TYPE" =~ ^(PWA|RA) ]]; then
  EXTRA=(--trim_front1 12 --trim_front2 12)   # MDA random-hexamer bias
  PROFILE="MDA-air (front-trim 12)"
else
  PROFILE="standard"
fi
echo "[$(date +%F_%T)] task $SLURM_ARRAY_TASK_ID  sample=$SAMPLE  type=$TYPE  profile=$PROFILE  node=$(hostname)"

# --- Node-local scratch (higher IO than NFS /home), cleaned on exit ----------
SCRATCH="/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$SCRATCH"
cleanup() { rm -rf "$SCRATCH"; }
trap cleanup EXIT

cp "$R1" "$R2" "$SCRATCH"/
L_R1="$SCRATCH/$(basename "$R1")"
L_R2="$SCRATCH/$(basename "$R2")"
O_R1="$SCRATCH/${SAMPLE}_R1.trimmed.fastq.gz"
O_R2="$SCRATCH/${SAMPLE}_R2.trimmed.fastq.gz"

# --- Run --------------------------------------------------------------------
module purge
module load Fastp/fastp

fastp \
  -i "$L_R1" -I "$L_R2" \
  -o "$O_R1" -O "$O_R2" \
  --detect_adapter_for_pe \
  --length_required 50 \
  --thread "$SLURM_CPUS_PER_TASK" \
  --json "$SCRATCH/${SAMPLE}.fastp.json" \
  --html "$SCRATCH/${SAMPLE}.fastp.html" \
  "${EXTRA[@]}"

# --- Persist trimmed reads + reports back to /home --------------------------
cp "$O_R1" "$O_R2" "$TRIMDIR"/
cp "$SCRATCH/${SAMPLE}.fastp.json" "$SCRATCH/${SAMPLE}.fastp.html" "$REPDIR"/

echo "[$(date +%F_%T)] done sample=$SAMPLE"
