# 03_candidates.py
# Selects candidate CoxL ORFs (best hit is a CoxL reference, alignment >= 150 aa) and dereplicates them per library
# (single-cell versus metagenome copies and overlapping ORFs: >= 98% identity over >= 100 aa, DIAMOND all-versus-all).
# Outputs: cand_raw.faa, candidates_dereplicated.tsv, candidates_dereplicated.faa
# Run    : python3 03_candidates.py (from this folder)

import csv, re, collections, subprocess
seqs = {}; h = None
for ln in open('coxL_hit_orfs.faa'):
    if ln.startswith('>'): h = ln[1:].split()[0]; seqs[h] = ''
    else: seqs[h] += ln.strip().rstrip('*')
best = {}
for r in csv.reader(open('coxL_hits.tsv'), delimiter='\t'):
    q, s, pid, L, bs = r[0], r[1], float(r[2]), int(r[3]), float(r[9])
    if q not in best or bs > best[q]['bits']:
        best[q] = dict(ref=s, pid=pid, aln=L, qlen=int(r[10]), bits=bs, sstart=int(r[6]), send=int(r[7]))
lib = lambda q: re.search(r'(PWA\d+|RA\d+)', q).group(1)
cand = [q for q, b in best.items() if ('_FI_' in b['ref'] or '_FII_' in b['ref']) and b['aln'] >= 150]
print('candidate ORFs (pre-dereplication):', len(cand))
with open('cand_raw.faa', 'w') as fh:
    for q in cand: fh.write(f'>{q}\n{seqs[q]}\n')
if False: subprocess.run('diamond makedb --in cand_raw.faa -d cand_raw --quiet && diamond blastp -d cand_raw -q cand_raw.faa --quiet '
               '--max-target-seqs 500 -e 1e-20 -f 6 qseqid sseqid pident length -o cand_self.tsv', shell=True, check=True)
parent = {q: q for q in cand}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for r in csv.reader(open('cand_self.tsv'), delimiter='\t'):
    a, b, pid, L = r[0], r[1], float(r[2]), int(r[3])
    if a != b and lib(a) == lib(b) and pid >= 98 and L >= 100: parent[find(a)] = find(b)
groups = collections.defaultdict(list)
for q in cand: groups[find(q)].append(q)
rows = []
for g, mem in groups.items():
    rep = max(mem, key=lambda q: len(seqs[q]))     # longest member represents the gene
    b = best[rep]
    rows.append(dict(rep=rep, library=lib(rep), n_members=len(mem), members=';'.join(mem), len_aa=len(seqs[rep]),
                     best_ref=b['ref'], best_pid=round(b['pid'], 1), aln=b['aln'],
                     best_form='I' if '_FI_' in b['ref'] else 'II'))
rows.sort(key=lambda r: (r['library'], r['rep']))
with open('candidates_dereplicated.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
with open('candidates_dereplicated.faa', 'w') as fh:
    for r in rows: fh.write(f">CAND_{r['library']}_{r['rep'].split('.')[-1][:60]}\n{seqs[r['rep']]}\n")
print('dereplicated genes:', len(rows))
print('per library x best-hit form:', dict(collections.Counter((r['library'], r['best_form']) for r in rows)))
