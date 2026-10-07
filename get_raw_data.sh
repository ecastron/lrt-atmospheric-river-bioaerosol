#!/usr/bin/env bash
# get_raw_data.sh
# Downloads the raw reads of BioProject PRJNA1540148 from the NCBI SRA and restores the original MiSeq/NovaSeq file
# names in the folders the pipeline expects, using data/sra_run_manifest.tsv (one row per SRA run: run accession,
# BioSample, sample, marker, destination folder, original R1/R2 file names).
# Requires: sra-tools (prefetch, fasterq-dump) and gzip (pigz is used when available).
# Run    : bash get_raw_data.sh [16S|18S|ITS|shotgun ...]   (from the repository root; no argument = all markers)
# Sizes  : amplicons about 8 GB compressed; the 15 shotgun libraries about 57 GB compressed.

set -euo pipefail

MANIFEST="data/sra_run_manifest.tsv"
THREADS="${THREADS:-4}"
TMP="${TMPDIR:-/tmp}/lrt_sra_$$"
[[ -s "$MANIFEST" ]] || { echo "ERROR: $MANIFEST not found; run from the repository root" >&2; exit 1; }
command -v fasterq-dump >/dev/null || { echo "ERROR: fasterq-dump (sra-tools) not found" >&2; exit 1; }
GZ=gzip; command -v pigz >/dev/null && GZ="pigz -p $THREADS"
WANT="${*:-16S 18S ITS shotgun}"
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

# columns: run biosample sample_name marker library_name destination fastq_R1 fastq_R2
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r run biosample sample marker library dest r1 r2; do
  [[ " $WANT " == *" $marker "* ]] || continue
  if [[ -s "$dest/$r1" && -s "$dest/$r2" ]]; then echo "skip $run ($library): already present"; continue; fi
  mkdir -p "$dest"
  echo "get  $run ($library) -> $dest/"
  prefetch --max-size u -O "$TMP" "$run" >/dev/null
  fasterq-dump --split-files --threads "$THREADS" -O "$TMP" "$TMP/$run" >/dev/null
  # paired runs come out as <run>_1.fastq and <run>_2.fastq; restore the original names
  [[ -s "$TMP/${run}_1.fastq" && -s "$TMP/${run}_2.fastq" ]] || { echo "ERROR: $run did not yield paired reads" >&2; exit 1; }
  $GZ -c "$TMP/${run}_1.fastq" > "$dest/$r1"
  $GZ -c "$TMP/${run}_2.fastq" > "$dest/$r2"
  rm -rf "$TMP/$run" "$TMP/${run}"_*.fastq
done
echo "done"
