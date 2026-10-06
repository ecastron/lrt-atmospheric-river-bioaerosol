#!/bin/bash
#SBATCH --job-name=air_readsfunc
#SBATCH --partition=basic
#SBATCH --array=0-14%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=logs/readsfunc_%A_%a.out
#SBATCH --error=logs/readsfunc_%A_%a.err

# 26_air_reads_func.sh
# Read-level marker recruitment for the air libraries: DIAMOND blastx of the deduplicated reads against the curated
# marker protein database (>= 60% identity over >= 70% of the read), counting reads per marker. Counts are later
# normalized by protein length and by recA to give copies per genome; MDA preserves within-genome gene ratios, and
# read-level counting avoids assembly fragmentation.
# Set PROJ to the shotgun project directory (default: current directory).
# Run    : sbatch 26_air_reads_func.sh

set -euo pipefail
shopt -s nullglob
PROJ="${PROJ:-$PWD}"
DEDUP="$PROJ/air_dedup"
DB="$PROJ/atmos_markers/atmos_all_markers.dmnd"
OUT="$PROJ/air_func/reads_recruit"; mkdir -p "$OUT" logs

R1S=( "$DEDUP"/*PWA*_R1.dedup.fastq.gz "$DEDUP"/*RA[0-9]*_R1.dedup.fastq.gz )
(( ${#R1S[@]} == 15 )) || { echo "ERROR: ${#R1S[@]} air R1, expected 15" >&2; exit 1; }
R1="${R1S[$SLURM_ARRAY_TASK_ID]}"; R2="${R1/_R1.dedup.fastq.gz/_R2.dedup.fastq.gz}"
S="$(basename "$R1" | sed -E 's/_R1\.dedup\.fastq\.gz//')"

SCRATCH="/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_rf"
mkdir -p "$SCRATCH"; cleanup(){ rm -rf "$SCRATCH"; }; trap cleanup EXIT
zcat "$R1" "$R2" > "$SCRATCH/reads.fq"
nreads=$(( $(wc -l < "$SCRATCH/reads.fq") / 4 ))

set +u; module purge; module load Diamond-2.1.23/diamond-2.1.23; set -u
# best marker hit per read; >=60% aa id over >=70% of the (translated) read.
diamond blastx -d "$DB" -q "$SCRATCH/reads.fq" --threads "$SLURM_CPUS_PER_TASK" \
  --id 60 --query-cover 70 --max-target-seqs 1 --max-hsps 1 --quiet \
  -f 6 qseqid sseqid pident -o "$SCRATCH/hits.tsv"

# collapse to read counts per MARKER (parse marker name from reference id)
python3 - "$SCRATCH/hits.tsv" "$OUT/${S}.markercounts.tsv" "$S" "$nreads" <<'PY'
import sys, re, collections
hits, out, samp, nreads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
def marker(s):
    return s.split("__")[0] if "__" in s else re.sub(r'_(sequences_)?\d+$','',s)
c = collections.Counter()
for line in open(hits):
    q,s,pid = line.rstrip("\n").split("\t"); c[marker(s)] += 1
with open(out,"w") as fh:
    fh.write(f"#sample\t{samp}\tnreads\t{nreads}\n#marker\treads\n")
    for m,n in sorted(c.items(), key=lambda x:-x[1]): fh.write(f"{m}\t{n}\n")
PY
echo "[$(date +%F_%T)] $S done ($nreads reads) -> $OUT/${S}.markercounts.tsv"
