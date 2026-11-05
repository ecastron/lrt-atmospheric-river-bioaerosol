#!/usr/bin/env python3
# make_Additional_file_3_mixs_05NOV26.py
# Builds Additional file 3, the per-sample metadata sheet formatted to the MIxS MIMARKS survey
# standard of the Genomic Standards Consortium (https://genomicsstandardsconsortium.github.io/mixs/).
#
# Inputs (from the NCBI submission package, ../data/ in this repository):
#   biosample_attributes_06OCT26.tsv       the BioSample attribute sheet submitted to NCBI
#   biosample_accessions_SAMN_06OCT26.tsv  BioSample accessions returned by NCBI
# Output:
#   Additional_file_3_MIxS_metadata.tsv    one row per sequenced sample, 24 MIxS fields
#
# The NCBI submission used the generic BioSample package (12 attributes), which does not carry the
# MIxS core fields. Those are derived here: env_broad_scale from the sampling region, env_local_scale
# and env_medium from the sample type, and the sequencing fields from the Methods. Sample identifiers
# are unchanged, so rows join to BioProject PRJNA1540148 on sample_name.
#
# ENVO terms were verified against the EBI Ontology Lookup Service (ols4) on 2026-11-05.
#
# Run: python3 make_Additional_file_3_mixs_05NOV26.py

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
SRC = os.path.join(DATA, 'biosample_attributes_06OCT26.tsv')
ACC = os.path.join(DATA, 'biosample_accessions_SAMN_06OCT26.tsv')
OUT = os.path.join(HERE, 'Additional_file_3_MIxS_metadata.tsv')

# env_broad_scale must be a biome term; the Antarctic sites are polar, the Patagonian sites
# terrestrial, and the four seawater samples marine.
POLAR = 'polar biome [ENVO:01000339]'
TERR = 'terrestrial biome [ENVO:00000446]'
MARINE = 'marine biome [ENVO:00000447]'

LOCAL = {'air': 'atmospheric boundary layer [ENVO:01000323]',
         'soil': 'coast [ENVO:01000687]',
         'rhizosphere': 'coast [ENVO:01000687]',
         'seawater': 'coast [ENVO:01000687]',
         'sediment': 'glacier [ENVO:00000133]',
         'control': 'not applicable'}

MEDIUM = {'air': 'air [ENVO:00002005]',
          'soil': 'soil [ENVO:00001998]',
          'rhizosphere': 'rhizosphere [ENVO:00005801]',
          'seawater': 'coastal sea water [ENVO:00002150]',
          'sediment': 'sediment [ENVO:00002007]',
          'control': 'not applicable'}

SAMP_SIZE = {'air': '54 m3 air', 'seawater': '1.5 L water', 'soil': '0.25 g',
             'rhizosphere': '0.25 g', 'sediment': '0.25 g', 'control': 'not applicable'}

PRIMERS = ('16S rRNA gene 515F/806R FWD:GTGYCAGCMGCCGCGGTAA;REV:GGACTACNVGGGTWTCTAAT'
           ' | ITS1 ITS1f/ITS2 FWD:CTTGGTCATTTAGAGGAAGTAA;REV:GCTGCGTTCTTCATCGATGC')

COLS = ['sample_name', 'biosample_accession', 'bioproject_accession', 'investigation_type',
        'project_name', 'env_broad_scale', 'env_local_scale', 'env_medium', 'geo_loc_name',
        'lat_lon', 'collection_date', 'samp_collect_device', 'samp_size', 'isolation_source',
        'organism', 'target_gene', 'target_subfragment', 'pcr_primers', 'seq_meth', 'lib_layout',
        'nucl_acid_ext', 'nucl_acid_amp', 'replicate', 'description']


def sample_kind(row):
    """Map a BioSample row onto the six categories used for the ENVO assignments."""
    if row['isolation_source'].startswith('negative control'):
        return 'control'
    organism = row['organism']
    if 'air metagenome' in organism:
        return 'air'
    if 'rhizosphere' in organism:
        return 'rhizosphere'
    if 'seawater' in organism:
        return 'seawater'
    if 'sediment' in organism:
        return 'sediment'
    return 'soil'


def environments(kind, geo_loc_name):
    if kind == 'control':
        return 'not applicable', LOCAL[kind], MEDIUM[kind]
    if kind == 'seawater':
        broad = MARINE
    elif geo_loc_name.startswith('Antarctica'):
        broad = POLAR
    else:
        broad = TERR
    return broad, LOCAL[kind], MEDIUM[kind]


def load_accessions(path):
    """sample_name -> BioSample accession, from the sheet NCBI returns after processing."""
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        acc_col = next((c for c in reader.fieldnames if 'accession' in c.lower()), None)
        name_col = next((c for c in reader.fieldnames
                         if c.lower() in ('sample_name', 'sample name', 'isolate', 'library_id')), None)
        if not (acc_col and name_col):
            return {}
        return {r[name_col]: r[acc_col] for r in reader}


def main():
    rows = list(csv.DictReader(open(SRC), delimiter='\t'))
    accessions = load_accessions(ACC)

    with open(OUT, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t')
        writer.writeheader()
        for row in rows:
            kind = sample_kind(row)
            broad, local, medium = environments(kind, row['geo_loc_name'])
            writer.writerow({
                'sample_name': row['sample_name'],
                'biosample_accession': accessions.get(row['sample_name'], ''),
                'bioproject_accession': row['bioproject_accession'],
                'investigation_type': 'mimarks-survey',
                'project_name': ('Atmospheric river-driven delivery of Patagonian bioaerosol '
                                 'to the Antarctic Peninsula'),
                'env_broad_scale': broad,
                'env_local_scale': local,
                'env_medium': medium,
                'geo_loc_name': row['geo_loc_name'],
                'lat_lon': row['lat_lon'],
                'collection_date': row['collection_date'],
                'samp_collect_device': row['samp_collect_device'] or 'not applicable',
                'samp_size': SAMP_SIZE[kind],
                'isolation_source': row['isolation_source'],
                'organism': row['organism'],
                'target_gene': '16S rRNA gene; ITS',
                'target_subfragment': 'V4; ITS1',
                'pcr_primers': PRIMERS,
                'seq_meth': 'Illumina MiSeq',
                'lib_layout': 'paired',
                'nucl_acid_ext': 'DNeasy PowerSoil Pro kit (Qiagen)',
                'nucl_acid_amp': ('single-step PCR with barcoded primers '
                                  '(Earth Microbiome Project protocol)'),
                'replicate': row['replicate'],
                'description': row['description'],
            })

    matched = sum(1 for r in rows if accessions.get(r['sample_name']))
    print(f'wrote {OUT}: {len(rows)} rows, {len(COLS)} columns')
    print(f'BioSample accessions matched: {matched}/{len(rows)}')


if __name__ == '__main__':
    main()
