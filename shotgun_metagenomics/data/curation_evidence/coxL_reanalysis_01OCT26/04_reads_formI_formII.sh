#!/bin/bash
#SBATCH --job-name=coxL_reads_01OCT26
#SBATCH --partition=basic
#SBATCH --array=0-29%6
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --output=logs/reads_%A_%a.out

# 04_reads_formI_formII.sh
# Read-level recruitment of form I versus form II CoxL, with the same read sets and DIAMOND settings as
# 26_air_reads_func.sh (air, deduplicated) and 28_byers_reads_func.sh (Byers soil, trimmed); database = motif-labeled
# CoxL references (form I AYXCSFR, form II AYXGAGR) plus molybdenum hydroxylase outgroup. Tasks 0-14 air, 15-29 soil.
# Set PROJ to the shotgun project directory and COXL_DIR to this folder (defaults: current directory).
# Run    : sbatch 04_reads_formI_formII.sh

set -euo pipefail
shopt -s nullglob
PROJ="${PROJ:-$PWD}"; W="${COXL_DIR:-$PWD}"; mkdir -p "$W/reads_out" "$W/logs"
DB="$W/ref_coxL.dmnd"
SCRATCH="/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_cox"; mkdir -p "$SCRATCH"; trap 'rm -rf "$SCRATCH"' EXIT
i=$SLURM_ARRAY_TASK_ID
if (( i < 15 )); then
  R1S=( "$PROJ"/air_dedup/*PWA*_R1.dedup.fastq.gz "$PROJ"/air_dedup/*RA[0-9]*_R1.dedup.fastq.gz )
  (( ${#R1S[@]} == 15 )) || { echo "air R1 count ${#R1S[@]}"; exit 1; }
  R1="${R1S[$i]}"; R2="${R1/_R1.dedup.fastq.gz/_R2.dedup.fastq.gz}"
  S="air_$(basename "$R1" | sed -E 's/_R1\.dedup\.fastq\.gz//')"
  zcat "$R1" "$R2" > "$SCRATCH/reads.fq"
else
  SOIL="$PROJ/byers_soil"; mapfile -t RUNS < <(cut -f1 "$SOIL/byers_bulksoil.tsv")
  RUN="${RUNS[$((i-15))]}"; S="byers_${RUN}"
  set +u; module purge; module load Fastp/fastp 2>/dev/null || module load fastp 2>/dev/null; set -u
  fastp -i "$SOIL/${RUN}_1.fastq.gz" -I "$SOIL/${RUN}_2.fastq.gz" -o "$SCRATCH/t1.fq.gz" -O "$SCRATCH/t2.fq.gz" \
    --length_required 100 -q 20 -w "$SLURM_CPUS_PER_TASK" -j /dev/null -h /dev/null 2>"$SCRATCH/fastp.log"
  zcat "$SCRATCH/t1.fq.gz" "$SCRATCH/t2.fq.gz" > "$SCRATCH/reads.fq"
fi
nreads=$(( $(wc -l < "$SCRATCH/reads.fq") / 4 ))
set +u; module purge; module load Diamond-2.1.23/diamond-2.1.23; set -u
diamond blastx -d "$DB" -q "$SCRATCH/reads.fq" --threads "$SLURM_CPUS_PER_TASK" \
  --id 60 --query-cover 70 --max-target-seqs 1 --max-hsps 1 --quiet \
  -f 6 qseqid sseqid pident sstart send -o "$SCRATCH/hits.tsv"
cp "$SCRATCH/hits.tsv" "$W/reads_out/${S}.hits.tsv"
python3 - "$SCRATCH/hits.tsv" "$W/reads_out/${S}.counts.tsv" "$S" "$nreads" <<'PY'
import sys, collections
h, out, s, n = sys.argv[1:5]; c = collections.Counter()
for ln in open(h):
    sid = ln.split('\t')[1]
    c['formI' if '_FI_' in sid else ('formII' if '_FII_' in sid else 'outgroup')] += 1
open(out, 'w').write(f"sample\tnreads\tformI\tformII\toutgroup\n{s}\t{n}\t{c['formI']}\t{c['formII']}\t{c['outgroup']}\n")
PY
echo "$S done: $nreads reads"
