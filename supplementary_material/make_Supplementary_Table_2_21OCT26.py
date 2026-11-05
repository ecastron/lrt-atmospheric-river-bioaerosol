# make_Supplementary_Table_2_21OCT26.py
# Builds Supplementary_Table_2_permdisp_pairwise.xlsx from 16S/figures/permdisp_pairwise_21OCT26.tsv
# (written by 16S/permdisp_pairwise_21OCT26.R). Sheet "Summary" follows the Supplementary Table 1 layout:
# title, notes, blank row, header, body.
import csv, os
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
TSV = os.path.join(HERE, '..', '16S', 'figures', 'permdisp_pairwise_21OCT26.tsv')

lab_group = {'Antarctic Air': 'Antarctic Peninsula air', 'Antarctic Soil': 'Antarctic Peninsula soil',
             'Patagonian Air': 'Patagonian air', 'Patagonian Soil': 'Patagonian soil',
             'SO Air': 'Southern Ocean air, over sea', 'SO Air over land': 'Southern Ocean air, over land',
             'SO Precipitation': 'Southern Ocean precipitation',
             'Antarctic': 'Antarctic Peninsula', 'Patagonian': 'Patagonia', 'SO': 'Southern Ocean'}
lab_cmp = {'Fig1 all samples': 'Fig. 1, all samples', 'Fig1 air only': 'Fig. 1, air only',
           'Fig1 this study': 'Fig. 1, this study', 'Fig5 Risopatron air': 'Fig. 5, Risopatrón air'}
lab_factor = {'source2': 'sample type', 'region_air': 'region', 'medium': 'medium (air vs soil)',
              'region': 'region (Antarctic Peninsula vs Patagonia)', 'cell': 'region × medium',
              'phase': 'AR phase', 'ar2': 'Before vs During+After'}

def num(v, nd):
    return None if v in ('NA', '') else round(float(v), nd)

rows = []
with open(TSV) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        pair = r['pair']
        if pair != 'NA':
            pair = ' vs '.join(lab_group.get(g, g) for g in pair.split(' vs '))
        else:
            pair = 'all groups'
        rows.append((lab_cmp[r['comparison']], 'pairwise' if r['test'] == 'pairwise PERMANOVA' else r['test'], lab_factor[r['factor']], pair, r['n'],
                     num(r['F'], 2), num(r['R2'], 3), num(r['p'], 4), num(r['p_BH'], 4)))

wb = openpyxl.Workbook()
ws = wb.active
ws.title = 'Summary'
ws.append(['Supplementary Table 2. Dispersion tests and pairwise PERMANOVA for the PERMANOVAs reported in Figures 1 and 5.'])
ws.append(['PERMDISP: vegan::betadisper (distance to centroid) with vegan::permutest on weighted UniFrac distances; '
           'a nonsignificant result means group spreads do not differ, so a significant PERMANOVA reflects a difference in composition.'])
ws.append(['Pairwise (pairwise PERMANOVA): vegan::adonis2 on each pair of groups; p_BH is Benjamini-Hochberg adjusted within each family of pairs. '
           'All tests used 9,999 permutations with set.seed(100). Distances and sample sets are those of Figures 1 and 5.'])
ws.append(['With two AR-core (During) samples, pairs involving During have few distinct permutations (minimum p of about 0.007 for Before vs During) and little power.'])
ws.append([])
ws.append(['Analysis', 'Test', 'Factor', 'Groups compared', 'n', 'F', 'R²', 'p', 'p_BH'])
for r in rows:
    ws.append(list(r))
wb.save(os.path.join(HERE, 'Supplementary_Table_2_permdisp_pairwise.xlsx'))
print(len(rows), 'rows')
