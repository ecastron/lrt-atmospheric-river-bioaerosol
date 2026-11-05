Supplementary Data 2. Marker protein database used for read-level functional profiling of the air metagenomes.

marker_database_proteins.faa.gz
    Protein reference searched with DIAMOND blastx (>= 60% identity over >= 70% of the read).
    211,249 sequences: 88 biogeochemical marker genes curated for this study (headers
    <marker>_sequences_<n>) plus 2,305 reviewed, full-length bacterial UniProtKB/Swiss-Prot
    proteins for 18 additional markers (headers <marker>__<category>__<UniProt accession>).
marker_sources.tsv
    Each of the 106 markers with its source (curated set or Swiss-Prot supplement) and category.
swissprot_supplement_annotation.tsv
    One row per Swiss-Prot supplement protein: accession, marker, category, description.
marker_lengths.tsv
    Median protein length per marker (amino acids), used for length normalization.
marker_curation.tsv
    Markers kept for or dropped from the air versus soil enrichment test, with the evidence
    for each decision. Dropped generic families (CoxL, [NiFe]- and [FeFe]-hydrogenase, generic
    methanol dehydrogenase, formate dehydrogenase) were resolved separately on assembled ORFs (CoxL by
    form I active-site motif AYXCSFR plus form I clade placement; [NiFe]-hydrogenase by HydDB group).

recA and uvrA are single-copy controls: markers are expressed as copies per genome relative to
recA, and libraries were retained when uvrA fell between 0.6 and 1.6 copies per genome.
