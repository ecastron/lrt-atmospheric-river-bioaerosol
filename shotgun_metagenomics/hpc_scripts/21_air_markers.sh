#!/bin/bash
#SBATCH --job-name=air_markers
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/air_markers_%j.out
#SBATCH --error=logs/air_markers_%j.err

# 21_air_markers.sh
# Searches all predicted air ORFs against the curated marker protein database (Supplementary Data 2) with DIAMOND
# blastp and builds a non-redundant nucleotide catalogue of marker-hit ORFs (cd-hit-est, 95% identity).
# Set PROJ to the shotgun project directory (default: current directory).
# Run    : sbatch 21_air_markers.sh

set -euo pipefail
PROJ="${PROJ:-$PWD}"
ORF="$PROJ/air_func/orfs"
DB="$PROJ/atmos_markers/atmos_all_markers.dmnd"
OUT="$PROJ/air_func"; mkdir -p "$OUT" logs
cd "$OUT"

cat "$ORF"/*.faa > all_orfs.faa
cat "$ORF"/*.fna > all_orfs.fna
echo "[$(date +%F_%T)] total ORFs: $(grep -c '^>' all_orfs.faa)"

set +u; module purge 2>/dev/null; module load Diamond-2.1.23/diamond-2.1.23 2>/dev/null; set -u
# protein-protein; thresholds match the methods (>=40% id, >=70% query cover);
# one best marker per ORF.
diamond blastp -d "$DB" -q all_orfs.faa --threads "$SLURM_CPUS_PER_TASK" \
  --id 40 --query-cover 70 --max-target-seqs 1 --max-hsps 1 --quiet \
  -f 6 qseqid sseqid pident length qcovhsp bitscore -o marker_blastp.tsv
echo "[$(date +%F_%T)] marker-hit ORFs: $(wc -l < marker_blastp.tsv)"

# parse marker name from the reference id; write ORF->marker map + presence
python3 - <<'PY'
import re, csv, collections
def marker(s):
    if "__" in s: return s.split("__")[0]                 # supplement: gene__category__acc
    return re.sub(r'_(sequences_)?\d+$','', s)             # Guajardo: Name_..._N
orf2m={}
for line in open("marker_blastp.tsv"):
    q,s,pid,ln,qc,bs=line.rstrip("\n").split("\t")
    orf2m[q]=marker(s)
# sample id is the assembly tag before "__" in the ORF id (we tagged headers)
pres=collections.defaultdict(set)
with open("air_orf_marker_map.tsv","w") as fh:
    fh.write("orf\tsample_assembly\tmarker\n")
    for orf,m in orf2m.items():
        samp=orf.split("__")[0]            # e.g. 006_6_RA8_GMCF3029a.meta
        fh.write(f"{orf}\t{samp}\t{m}\n"); pres[samp].add(m)
print("samples with marker hits:", len(pres))
# keep the ORF id list for fna extraction
open("marker_orf_ids.txt","w").write("\n".join(orf2m)+"\n")
PY

# extract marker-ORF nucleotide seqs, then dereplicate to a gene catalogue
set +u; module purge 2>/dev/null; module load SAMTOOLS/samtools.1.20 2>/dev/null; set -u
# (use seqkit-free extraction via awk: pull records whose header id is in the list)
awk 'NR==FNR{keep[$1]=1; next} /^>/{id=substr($1,2); p=(id in keep)} p' marker_orf_ids.txt all_orfs.fna > marker_orfs.fna
echo "[$(date +%F_%T)] marker ORF nt seqs: $(grep -c '^>' marker_orfs.fna)"

set +u; module purge 2>/dev/null; module load CD-HIT/cd-hit 2>/dev/null; set -u
cd-hit-est -i marker_orfs.fna -o marker_gene_catalog.fna -c 0.95 -d 0 -T "$SLURM_CPUS_PER_TASK" -M 40000 >/dev/null
echo "[$(date +%F_%T)] dereplicated catalogue: $(grep -c '^>' marker_gene_catalog.fna) genes -> $OUT/marker_gene_catalog.fna"
