# 06_context.py
# Gene context of each dereplicated CoxL candidate across all assembly copies. Neighbors are ORFs on the same contig
# within +-10 positions (Prodigal order); accessory genes are called by best DIAMOND hit (e-value <= 1e-10, identity
# >= 30%, alignment >= 50 aa) to UniProt CoxM/CoxS/CoxG/CoxD/CoxE/CoxF. coxMSL: coxM and coxS both within +-3 ORFs of
# coxL, ordered M-S-L along the coding strand.
# Output : candidates_context.tsv
# Run    : python3 06_context.py (from this folder)

import csv, re, collections
cand = list(csv.DictReader(open('candidates_dereplicated.tsv'), delimiter='\t'))
need = set(); contig_of = lambda o: o.rsplit('_', 1)[0]
for c in cand:
    for m in c['members'].split(';'): need.add(contig_of(m))
coords = collections.defaultdict(list)   # contig -> [(idx, orf, start, end, strand)]
with open('all_orf_coords.tsv') as fh:
    for ln in fh:
        o, s, e, st = ln.rstrip('\n').split('\t')[:4]
        ctg = contig_of(o)
        if ctg in need: coords[ctg].append((int(o.rsplit('_', 1)[1]), o, int(s), int(e), int(st)))
acc = {}
for r in csv.reader(open('cox_accessory_hits.tsv'), delimiter='\t'):
    q, s, pid, L, bs = r[0], r[1], float(r[2]), int(r[3]), float(r[9])
    if pid >= 30 and L >= 50 and (q not in acc or bs > acc[q][1]): acc[q] = (s.split('__')[0], bs, pid)
out = []
for c in cand:
    best_ctx = None
    for m in c['members'].split(';'):
        ctg = contig_of(m); L = int(re.search(r'length_(\d+)', ctg).group(1)); orfs = sorted(coords[ctg])
        idx = int(m.rsplit('_', 1)[1]); strand = next((x[4] for x in orfs if x[0] == idx), 0)
        near = {x[0]: acc.get(x[1], (None,))[0] for x in orfs if abs(x[0] - idx) <= 10 and x[0] != idx}
        genes = {g for g in near.values() if g}
        pos = lambda g: [i for i, v in near.items() if v == g and abs(i - idx) <= 3]
        mP, sP = pos('coxM'), pos('coxS')
        # on + strand, M-S-L means indices increase M < S < L; on - strand they decrease
        msl = bool(mP and sP and any((strand > 0 and a < b < idx) or (strand < 0 and a > b > idx) for a in mP for b in sP))
        n_orfs = len(orfs)
        ctx = dict(contig_len=L, n_orfs_on_contig=n_orfs, coxM=bool(mP), coxS=bool(sP), coxMSL_order=msl,
                   coxG='coxG' in genes, coxD='coxD' in genes, coxE='coxE' in genes, coxF='coxF' in genes,
                   assessable=n_orfs >= 3, member=m)
        score = (ctx['coxMSL_order'], ctx['coxM'] + ctx['coxS'], ctx['coxG'] + ctx['coxD'] + ctx['coxE'] + ctx['coxF'], L)
        if best_ctx is None or score > best_ctx[0]: best_ctx = (score, ctx)
    out.append({**{k: c[k] for k in ('rep', 'library', 'len_aa', 'best_form', 'best_pid')}, **best_ctx[1]})
with open('candidates_context.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter='\t'); w.writeheader(); w.writerows(out)
print('context assessable (>=3 ORFs on contig):', sum(o['assessable'] for o in out), 'of', len(out))
print('coxM adjacent:', sum(o['coxM'] for o in out), '| coxS adjacent:', sum(o['coxS'] for o in out), '| coxMSL order:', sum(o['coxMSL_order'] for o in out))
print('coxG nearby:', sum(o['coxG'] for o in out), '| coxD/E/F nearby:', sum(o['coxD'] or o['coxE'] or o['coxF'] for o in out))
