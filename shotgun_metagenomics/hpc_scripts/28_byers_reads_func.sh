#!/bin/bash
#SBATCH --job-name=byers_func
#SBATCH --partition=basic
#SBATCH --array=0-14%5
#SBATCH --cpus-per-task=16
#SBATCH --mem=44G
#SBATCH --time=14:00:00
#SBATCH --output=logs/byers_%A_%a.out
#SBATCH --error=logs/byers_%A_%a.err

# 28_byers_reads_func.sh
# Soil comparator: the same read-level marker recruitment on the Byers Peninsula (Livingston Island, Antarctica)
# bulk-soil metagenomes (BioProject PRJNA746701). These libraries were not MDA-amplified; reads are trimmed with fastp
# (-q 20, minimum length 100) and searched with the same DIAMOND settings as the air reads.
# Set PROJ to the shotgun project directory (default: current directory).
# Run    : sbatch 28_byers_reads_func.sh

set -euo pipefail
PROJ="${PROJ:-$PWD}"
SOIL="$PROJ/byers_soil"
DB="$PROJ/atmos_markers/atmos_all_markers.dmnd"
OUT="$PROJ/air_func/byers_recruit"; mkdir -p "$OUT" logs

mapfile -t RUNS < <(cut -f1 "$SOIL/byers_bulksoil.tsv")
mapfile -t LIBS < <(cut -f2 "$SOIL/byers_bulksoil.tsv")
(( ${#RUNS[@]} == 15 )) || { echo "ERROR: ${#RUNS[@]} runs, expected 15" >&2; exit 1; }
RUN="${RUNS[$SLURM_ARRAY_TASK_ID]}"; LIB="${LIBS[$SLURM_ARRAY_TASK_ID]}"
R1="$SOIL/${RUN}_1.fastq.gz"; R2="$SOIL/${RUN}_2.fastq.gz"
[[ -s "$R1" && -s "$R2" ]] || { echo "ERROR: missing reads for $RUN" >&2; exit 1; }

SCRATCH="/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_by"
mkdir -p "$SCRATCH"; cleanup(){ rm -rf "$SCRATCH"; }; trap cleanup EXIT

set +u; module purge; module load Fastp/fastp 2>/dev/null || module load fastp 2>/dev/null; set -u
fastp -i "$R1" -I "$R2" -o "$SCRATCH/t1.fq.gz" -O "$SCRATCH/t2.fq.gz" \
  --length_required 100 -q 20 -w "$SLURM_CPUS_PER_TASK" -j /dev/null -h /dev/null 2>"$SCRATCH/fastp.log"
zcat "$SCRATCH/t1.fq.gz" "$SCRATCH/t2.fq.gz" > "$SCRATCH/reads.fq"
nreads=$(( $(wc -l < "$SCRATCH/reads.fq")/4 ))

set +u; module purge; module load Diamond-2.1.23/diamond-2.1.23; set -u
diamond blastx -d "$DB" -q "$SCRATCH/reads.fq" --threads "$SLURM_CPUS_PER_TASK" \
  --id 60 --query-cover 70 --max-target-seqs 1 --max-hsps 1 --quiet \
  -f 6 qseqid sseqid pident -o "$SCRATCH/hits.tsv"

python3 - "$SCRATCH/hits.tsv" "$OUT/byers_soil__${LIB}.markercounts.tsv" "$LIB" "$nreads" <<'PY'
import sys, re, collections
hits,out,samp,nreads = sys.argv[1:5]
def marker(s): return s.split("__")[0] if "__" in s else re.sub(r'_(sequences_)?\d+$','',s)
c=collections.Counter()
for line in open(hits):
    q,s,pid=line.rstrip("\n").split("\t"); c[marker(s)]+=1
with open(out,"w") as fh:
    fh.write(f"#sample\t{samp}\tmatrix\tbyers_soil\tnreads\t{nreads}\n#marker\treads\n")
    for m,n in sorted(c.items(),key=lambda x:-x[1]): fh.write(f"{m}\t{n}\n")
PY
echo "[$(date +%F_%T)] byers $LIB ($RUN) done, $nreads reads"
