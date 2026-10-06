# CoxL form I / form II classification

Convention (UniProt P19919, P19913): form I = AYXCSFR (bona fide CO dehydrogenase large subunit); form II = AYXGAGR (putative).

Pipeline, in order:

1. `01_build_references.py`: downloads 67 reference sequences (33 form I, including 11 named CO oxidizers; 30 form II; 4 Mo-hydroxylase outgroup).
2. `02_search_orfs.sh`: DIAMOND 2.1.23 blastp of all air ORFs against the references.
3. `03_candidates.py`: keeps ORFs whose best hit is CoxL with alignment length of at least 150 aa, then dereplicates (89 genes).
4. `05_tree.sh`: MAFFT 7.490 (`--auto`) and IQ-TREE 3.0.1 (LG+G4, 1,000 ultrafast bootstraps).
5. `06_context.py`: extracts the active-site motif region and genomic context for each candidate.
6. `07b_classify_v2.R`: combines motif and clade placement into per-gene calls.
7. `04_reads_formI_formII.sh` and `08_reads_norm.R`: read-level recruitment to form I and form II genes, normalized per recA.

A gene is called form I only if it carries the form I motif and sits in a form I clade with UFBoot of at least 95.

Main outputs: `candidates_classified_v2.tsv` (per-gene calls) and `reads_formI_formII_per_library.tsv` (read counts per library).

Shell scripts run on a Slurm cluster. Set `PROJ` to the shotgun project directory and `COXL_DIR` to this folder (both default to the current directory).
