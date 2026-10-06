# 01_build_references.py
# Builds CoxL reference sets labeled by active-site motif (form I = AYRCSFR, as in UniProt P19919 and P19913; form II
# = AYRGAGR), including named form I CoxL from organisms with demonstrated CO oxidation, an outgroup of molybdopterin
# hydroxylase large subunits, and CoxM/CoxS/CoxG/CoxD/CoxE/CoxF references for gene-context calls (UniProt REST API).
# Outputs: ref_coxL.faa, ref_coxL.tsv, ref_cox_accessory.faa
# Run    : python3 01_build_references.py (from this folder)

import urllib.request, urllib.parse, re, collections, json, csv, time
UA = {'User-Agent': 'LRT-coxL-reanalysis/1.0 (mailto:ecastron@utalca.cl)'}
F = 'accession,reviewed,protein_name,gene_names,organism_name,lineage,length,sequence'
def uq(q, size=None):
    url = 'https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=' + F + '&query=' + urllib.parse.quote(q)
    for attempt in range(3):
        try:
            rows = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300).read().decode().splitlines()
            break
        except Exception as e:
            time.sleep(5)
    h = rows[0].split('\t'); return [dict(zip(h, r.split('\t'))) for r in rows[1:]]
def motif(s):
    m = re.search(r'AY[A-Z]CSFR', s) or re.search(r'AY[A-Z]GAGR', s) or re.search(r'AY[A-Z]{4}R', s)
    return m.group(0) if m else 'NA'
def form(s):
    # form I = AYXCSFR (M. smegmatis, a verified atmospheric CO oxidizer, carries AYACSFR); form II = AYXGAGR
    if re.search(r'AY[A-Z]CSFR', s): return 'I'
    if re.search(r'AY[A-Z]GAGR', s): return 'II'
    return None
def phylum(lin):
    return next((x.strip().split(' (')[0] for x in lin.split(',') if '(phylum)' in x), '')
refs = []   # dicts: id, set, form, acc, org, phylum, motif, seq
# form convention: AY[X]CSFR form I (P19919, P19913, M. smegmatis AYACSFR); AY[X]GAGR form II
# (a) named organisms with demonstrated CO oxidation (form I expected)
named_q = {
 'Oligotropha carboxidovorans (P19919)': 'accession:P19919',
 'Hydrogenophaga pseudoflava (P19913)': 'accession:P19913',
 'Mycolicibacterium smegmatis': '(gene:coxL OR gene:cutL OR protein_name:"carbon monoxide dehydrogenase large") AND (organism_name:"smegmatis") AND (length:[700 TO 900])',
 'Thermomicrobium roseum': '(gene:coxL OR protein_name:"carbon monoxide dehydrogenase large") AND (organism_name:"Thermomicrobium roseum") AND (length:[700 TO 900])',
 'Thermogemmatispora': '(gene:coxL OR protein_name:"carbon monoxide dehydrogenase large") AND (organism_name:"Thermogemmatispora") AND (length:[700 TO 900])',
 'Ktedonobacter racemifer': '(gene:coxL OR protein_name:"carbon monoxide dehydrogenase large") AND (organism_name:"Ktedonobacter racemifer") AND (length:[700 TO 900])',
 'Paraburkholderia': '(gene:coxL OR protein_name:"carbon monoxide dehydrogenase large") AND (organism_name:"Paraburkholderia xenovorans") AND (length:[700 TO 900])',
 'Cupriavidus': '(gene:coxL OR protein_name:"carbon monoxide dehydrogenase large") AND (organism_name:"Cupriavidus") AND (length:[700 TO 900])',
}
for label, q in named_q.items():
    rows = uq(q); keep = [r for r in rows if form(r['Sequence']) == 'I']
    print(f'named {label}: {len(rows)} hits, {len(keep)} with form I motif', [motif(r["Sequence"]) for r in rows][:6])
    for r in keep[:2]:
        refs.append(dict(set='named_formI', form='I', acc=r['Entry'], org=r['Organism'], phylum=phylum(r.get('Taxonomic lineage', '')),
                         motif=motif(r['Sequence']), seq=r['Sequence']))
# (b) broad motif-labeled sets, one per genus per form, up to 30 per form spread across phyla
rows = uq('((gene:coxL) OR (gene:cutL)) AND (taxonomy_id:2) AND (length:[750 TO 850]) AND (fragment:false)')
print('broad UniProt CoxL hits:', len(rows), collections.Counter(form(r['Sequence']) for r in rows))
rows.sort(key=lambda r: (r['Reviewed'] != 'reviewed', r['Organism']))
named_acc = {r['acc'] for r in refs}
pick = collections.OrderedDict()
for r in rows:
    f = form(r['Sequence'])
    if not f or r['Entry'] in named_acc: continue
    g = r['Organism'].split()[0].strip('[]')
    if (g, f) not in pick:
        pick[(g, f)] = dict(set='broad', form=f, acc=r['Entry'], org=r['Organism'], phylum=phylum(r.get('Taxonomic lineage', '')),
                            motif=motif(r['Sequence']), seq=r['Sequence'])
for f in ('I', 'II'):
    L = [v for v in pick.values() if v['form'] == f]; byph = collections.defaultdict(list)
    for v in L: byph[v['phylum']].append(v)
    out = []
    while len(out) < 30 and any(byph.values()):
        for ph in list(byph):
            if byph[ph] and len(out) < 30: out.append(byph[ph].pop(0))
    refs += out
    print(f'broad form {f}: {len(L)} genera available, {len(out)} kept; phyla', dict(collections.Counter(v["phylum"] for v in out)))
# (c) outgroup: molybdopterin hydroxylase large subunits (xanthine dehydrogenase / aldehyde oxidoreductase)
og = uq('(accession:Q46509 OR accession:O54051 OR accession:Q46799 OR accession:P80457)')  # D. gigas MOP, R. capsulatus XdhB, E. coli XdhA, bovine XDH, human AOX? (whatever resolves)
for r in og:
    refs.append(dict(set='outgroup', form='outgroup', acc=r['Entry'], org=r['Organism'], phylum=phylum(r.get('Taxonomic lineage', '')),
                     motif=motif(r['Sequence']), seq=r['Sequence']))
print('outgroup:', [(r['Entry'], r['Organism'][:30], len(r['Sequence'])) for r in og])
# (d) accessory Cox proteins for gene-context calls (Bacteria, any form), capped
acc_sets = {}
for gene in ['coxM', 'coxS', 'coxG', 'coxD', 'coxE', 'coxF']:
    rr = uq(f'(gene:{gene}) AND (taxonomy_id:2) AND (fragment:false)')
    rr.sort(key=lambda r: (r['Reviewed'] != 'reviewed', r['Organism']))
    seen = set(); keep = []
    for r in rr:
        g = r['Organism'].split()[0]
        if g in seen: continue
        seen.add(g); keep.append(r)
        if len(keep) >= 60: break
    acc_sets[gene] = keep
    print(f'{gene}: {len(rr)} hits, {len(keep)} kept')
# write outputs
with open('ref_coxL.faa', 'w') as fh, open('ref_coxL.tsv', 'w', newline='') as th:
    w = csv.writer(th, delimiter='\t'); w.writerow(['id', 'set', 'form', 'accession', 'organism', 'phylum', 'motif', 'length'])
    for i, r in enumerate(refs):
        tag = re.sub(r'[^A-Za-z0-9]+', '_', r['org'])[:35]
        rid = f"REF{i:03d}_{'F' + r['form'] if r['form'] in ('I', 'II') else 'OUT'}_{tag}_{r['acc']}"
        r['id'] = rid
        fh.write(f'>{rid}\n{r["seq"]}\n'); w.writerow([rid, r['set'], r['form'], r['acc'], r['org'], r['phylum'], r['motif'], len(r['seq'])])
with open('ref_cox_accessory.faa', 'w') as fh:
    for gene, L in acc_sets.items():
        for r in L:
            fh.write(f">{gene}__{r['Entry']}__{re.sub(r'[^A-Za-z0-9]+', '_', r['Organism'])[:30]}\n{r['Sequence']}\n")
print('written', len(refs), 'CoxL refs')
