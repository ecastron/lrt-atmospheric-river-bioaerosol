# Atmospheric river-driven delivery of Patagonian bioaerosol to the Antarctic Peninsula: analysis code

## Overview

This repository holds the analysis code for a study of airborne microbial transport from southern Patagonia to the
Antarctic Peninsula during the austral summer of 2021-2022. We sampled near-surface air, bulk soil and rhizosphere
along a 1,500 km transect (Puerto Natales, San Gregorio, Punta Arenas, Porvenir and Puerto Williams in Patagonia;
Risopatron Base and Yelcho Base in Antarctica) and sequenced 16S rRNA gene V4, ITS1 and 18S rRNA gene amplicons.
Fifteen air samples were also shotgun sequenced after multiple displacement amplification. The analyses test whether
Patagonian bioaerosol reaches the Antarctic Peninsula, with a focus on the atmospheric river that reached Risopatron
Base on 7-8 February 2022.

The code covers amplicon processing (DADA2), placement of the samples among Southern Ocean air from Malard et al.
(2022), contamination assessment (decontam), community and phylogenetic analyses, source tracking (FEAST), functional
inference (FAPROTAX), shotgun marker-gene recruitment and the figure scripts. The atmospheric transport simulations
(HYSPLIT and CHIMERE) were run separately; this repository only contains the script that assembles their output into
Figure 4.

Each script starts with a header that lists its purpose, inputs, outputs and the command to run it.

## Data sources

| Data | Source |
|---|---|
| Raw amplicon (16S, ITS, 18S) and shotgun reads | NCBI SRA, BioProject PRJNA1540148 |
| Southern Ocean air 16S reads (Malard et al. 2022) | NCBI SRA, BioProject PRJNA697829 |
| Bacterial and fungal phyloseq objects | Supplementary Data 1 of the paper |
| Curated atmospheric marker protein database | Supplementary Data 2 of the paper |
| Sample metadata and small intermediate tables (decontam results, FAPROTAX output, shotgun marker counts, CoxL candidates and tree) | This repository, at the paths the scripts read |
| Risopatron weather-station record | Not included |

Reference databases are downloaded separately: SILVA 138.1 DADA2 training and species files, DADA2 EUK SSU v1.9,
the UNITE general FASTA release (19.02.2025) and FAPROTAX 1.2.12.

## Expected directory layout

Scripts are run from the repository root unless their header says otherwise. The repository already contains the
metadata files and the intermediate tables. Add the remaining inputs as follows:

```
./
  DADA2_EUK_SSU_v1.9.fasta                  18S reference
  raw_data/16S/                             demultiplexed 16S FASTQ files (PRJNA1540148)
  18S/data/                                 demultiplexed 18S FASTQ files (PRJNA1540148)
  16S/
    silva_nr99_v138.1_train_set.fa, silva_species_assignment_v138.1.fa
    ACEDATA/Allfastq/                       Malard et al. 2022 FASTQ files (PRJNA697829, uncompressed)
    FAPROTAX_1.2.12/                        FAPROTAX release
    ps_merged_withTree_noorg_01OCT26.RDS    Supplementary Data 1: bacteria_16S_phyloseq.rds, renamed
  ITS/
    trimmed/                                ITS R1 FASTQ files after cutadapt (*_R1.fastq.gz)
    sh_general_release_all_19.02.2025/      UNITE reference
    ps_LRT_ITS_R1only_decontam_withTree_29MAY26.RDS   Supplementary Data 1: fungi_ITS_phyloseq.rds, renamed
  wind_rose/DatosAll_CR1000X_RPT.xlsx       weather-station record (not included)
  New_Hysplit/, New_Chimere/                HYSPLIT and CHIMERE figure panels (not included)
```

With the two Supplementary Data 1 objects in place, every bacterial and fungal figure script can be run without
repeating the upstream steps. `feast_stats_noorg_01OCT26.R` and `build_ps_noorg_01OCT26.R` also read
`16S/ps_merged_withTree_07MAY26.RDS`, the object produced by the upstream 16S chain (steps 1 to 5 below).
`shotgun_metagenomics/data/curation_evidence/coxL_reanalysis_01OCT26/06_context.py` reads `all_orf_coords.tsv`
(ORF coordinates from the Prodigal step, 207 MB), which is not included; rerun `20_prodigal_air.sh` to create it.

## Run order

### 16S rRNA gene

1. `16S/LRT.R`: DADA2 for the LRT 16S and 18S libraries. Run interactively from the `16S/` folder.
2. `16S/ACE_dada2.R`: DADA2 for the Malard et al. (2022) reads.
3. `16S/harmonize_metadata.R`: harmonizes LRT and Malard metadata.
4. `16S/merge_LRT_ACE.R`: merges both sequence tables, removes chimeras and assigns SILVA taxonomy.
5. `16S/build_tree.R`: MAFFT alignment and FastTree phylogeny of the merged ASVs.
6. `16S/build_ps_noorg_01OCT26.R`: removes chloroplast and mitochondrial ASVs. Its output is the input for all
   bacterial analyses.
7. `16S/decontam_controls_11MAY26.R`: decontam prevalence test against the buffer controls.
8. Analysis and figure scripts (see the table below).

Steps 2 to 5 tag their outputs `07MAY26`, the tag of the objects used in the manuscript, so each script finds the
output of the previous one.

The FEAST scripts in `16S/sensitivity_organelles_01OCT26/` read the organelle-free object from that folder and write
their output there. Copy `ps_merged_withTree_noorg_01OCT26.RDS` into the folder before running them, and copy
`feast_AR_event_01OCT26_proportions.tsv` to `16S/figures/` before running `figure6_feast_apportionment_01OCT26.R`.

### ITS (forward reads only)

1. Trim the conserved small-subunit primer region from R1 with cutadapt 5.2 into `ITS/trimmed/`.
2. `ITS/dada2_its_R1only.R`
3. `ITS/its_r1only_decontam_29MAY26.R`
4. `ITS/its_r1only_assigntax_29MAY26.R`
5. Build the ASV tree with VeryFastTree 4.0.5 (`ITS/asvs_r1only_decontam_29MAY26.nwk`).
6. `ITS/its_r1only_finalize_29MAY26.R`

### Shotgun metagenomics

The scripts in `shotgun_metagenomics/hpc_scripts/` run on a Slurm cluster in numeric order (fastp, clumpify, SPAdes,
Prodigal, marker search, read recruitment for air and for the Byers Peninsula soil comparator, MetaPathways LCA).
Set `PROJECT` or `PROJ` to the project directory (default: current directory) and `FASTQ_ROOT` to the raw read folder.
The CoxL classification pipeline is in `shotgun_metagenomics/data/curation_evidence/coxL_reanalysis_01OCT26/` and
has its own README. The R scripts in `shotgun_metagenomics/scripts/` run locally on the recruitment tables.

### Weather station

1. `wind_rose/fix_decimal_separator_02OCT26.py` (run from `wind_rose/`)
2. `wind_rose/wind_rose_phases_02OCT26.R`
3. `16S/figureS5_synoptic_02OCT26.py`

## Figures and tables

| Item | Script(s) |
|---|---|
| Figure 1 | `16S/figure_global_comparison_02OCT26.R` |
| Figure 2 | `16S/figure2_lrt_overview_01OCT26.R` |
| Figure 3 | `16S/figure3_export_faprotax_02OCT26.R`, then FAPROTAX 1.2.12 `collapse_table.py` (command given in that script), then `16S/figure3_faprotax_plot_byMedium_02OCT26.R` |
| Figure 4 | `16S/figure4_atmospheric_transport_17SEP26.py` (composite of HYSPLIT and CHIMERE panels from separate model runs) |
| Figure 5 | `16S/figure5_ar_event_01OCT26.R` (includes the betaNTI and RC_Bray null models) |
| Figure 6 | `16S/sensitivity_organelles_01OCT26/feast_AR_event_noorg_01OCT26.R`, `feast_AR_event_postprocess_noorg_01OCT26.R`, `feast_stats_noorg_01OCT26.R`, then `16S/figure6_feast_apportionment_01OCT26.R` |
| Supplementary Figure 1 | `16S/figureS_study_sites_17SEP26.py` |
| Supplementary Figure 2 | `shotgun_metagenomics/scripts/air_faprotax_corroboration_29SEP26.R`, `shotgun_metagenomics/scripts/figureS_shotgun_faprotax_05OCT26.R` |
| Supplementary Figure 3 | `shotgun_metagenomics/data/curation_evidence/coxL_reanalysis_01OCT26/` pipeline, `shotgun_metagenomics/scripts/figureS3_coxL_iqtree_05OCT26.R` |
| Supplementary Figure 4 | `shotgun_metagenomics/scripts/air_nife_hyddb_29SEP26.R` |
| Supplementary Figure 5 | `wind_rose/fix_decimal_separator_02OCT26.py`, `wind_rose/wind_rose_phases_02OCT26.R`, `16S/figureS5_synoptic_02OCT26.py` |
| Supplementary Figure 6 | HYSPLIT back-trajectories (separate model runs, not in this repository) |
| Supplementary Figure 7 | CHIMERE simulations (separate model runs, not in this repository) |
| Supplementary Figure 8 | `16S/decontam_controls_11MAY26.R`, `16S/decontam_betaNTI_01OCT26.R`, `16S/figureS8_decontam_checks_01OCT26.R` |
| Supplementary Figure 9 | `16S/figureS9_fungi_source_29SEP26.R` |
| Supplementary Table 1 | `16S/patagonia_transfer_sweep_02OCT26.R`, `16S/supp_table1_patagonia_shared_asvs_02OCT26.py` |

`16S/palettes.R` defines the colors used by all R figure scripts. The fungal AR-event panel is drawn by
`16S/figure5_ar_event_fungi_r1only_29MAY26.R`, and `16S/ar_patagonia_attribution_v3_01OCT26.R` and
`16S/decontam_AR_arrival_crosscheck_01OCT26.R` provide the source attribution of AR-arrival taxa and its
contamination check.

## Software

R 4.5.3 with: ape 5.8.1, Biostrings 2.78.0, dada2 1.38.0, DECIPHER 3.6.0, decontam 1.30.0, dplyr 1.2.1,
FEAST 0.1.0, forcats 1.0.1, geosphere 1.6.8, ggh4x 0.3.1, ggplot2 4.0.3, ggrepel 0.9.8, ggtree 4.0.5,
lubridate 1.9.5, multcompView 0.1.11, openair 3.1.0, patchwork 1.3.2, phangorn 2.12.1, phyloseq 1.54.2,
picante 1.8.2, purrr 1.2.2, readr 2.2.0, rstatix 0.7.3, scales 1.4.0, stringr 1.6.0, tibble 3.3.1, tidyr 1.3.2,
vegan 2.7.3 (parallel is part of base R).

Python 3.13 with: cartopy 0.25.0, matplotlib 3.10.9, numpy 2.4.4, openpyxl 3.1.5, pandas 3.0.2, Pillow 12.2.0,
pypdf 6.13.2.

Command-line tools: DADA2 1.38, cutadapt 5.2, MAFFT 7.490, FastTree 2.1.11, VeryFastTree 4.0.5, IQ-TREE 3.0.1,
DIAMOND 2.1.23, fastp 1.3.3, BBTools 39.91, SPAdes 4.2.0, Prodigal 2.6.3, MetaPathways 3.5, FAPROTAX 1.2.12,
FEAST, decontam 1.30, HYSPLIT 5.2.2 and CHIMERE v2023r3.

The cluster scripts load software with `module load`; adjust the module names to your system.

## Reproducibility

Every stochastic step uses `set.seed(100)`. In the betaNTI null models, permutations are drawn serially from that
seed before parallel execution, so results do not depend on the number of cores.

## License

MIT. See `LICENSE`.

## Citation

Castro-Nallar E, Guajardo-Leiva S, Poblete-Castro I, Bozkurt D, Opazo C, Molina-Montenegro M, Díez B, Huneeus N,
Mailler S, Lapere R, Galbán-Malagón C. Atmospheric river-driven delivery of Patagonian bioaerosol to the Antarctic
Peninsula. In preparation.
