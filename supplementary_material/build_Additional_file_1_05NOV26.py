# build_Additional_file_1_05NOV26.py
# Assembles Additional file 1 (supplementary figures and tables) as a PDF from:
#   Supplementary_Figure_Legends_21OCT26.docx        figure legends (edit there, then rerun)
#   Supplementary_Figure_<n>_*.pdf                   figures; files marked "superseded" are skipped
#   Supplementary_Table_1_patagonia_shared_ASVs.xlsx "Summary" sheet is typeset here; full table ships as Excel
#   Supplementary_Table_2_permdisp_pairwise.xlsx     "Summary" sheet (built by make_Supplementary_Table_2_21OCT26.py)
# Figures are labelled Fig. S1..Fig. S10 and tables Table S1..S2, matching the manuscript.
# The figure PDFs and table workbooks are the supplementary deliverables themselves and are not
# duplicated in this repository; place them beside this script to rebuild the PDF.
# Writes Additional_file_1_05NOV26.tex and compiles it with xelatex (Unicode legends, vector figures).
# Run from anywhere: python3 build_Additional_file_1_05NOV26.py

import glob, os, re, subprocess, zipfile
from lxml import etree
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = 'Additional_file_1_05NOV26'
TITLE = 'Atmospheric river-driven delivery of Patagonian bioaerosol to the Antarctic Peninsula'
AUTHORS = ('Eduardo Castro-Nallar, Sergio Guajardo-Leiva, Ignacio Poblete-Castro, Deniz Bozkurt, Charlie Opazo, '
           'Marco Molina-Montenegro, Beatriz Díez, Nicolás Huneeus, Sylvain Mailler, Rémy Lapere, '
           'Cristóbal Galbán-Malagón')
CORRESPONDING = 'ecastron@utalca.cl'


def tex(s):
    """Escape LaTeX special characters in plain text."""
    s = s.replace('\\', r'\textbackslash{}')
    for a, b in [('&', r'\&'), ('%', r'\%'), ('$', r'\$'), ('#', r'\#'), ('_', r'\_'),
                 ('{', r'\{'), ('}', r'\}'), ('~', r'\textasciitilde{}'), ('^', r'\textasciicircum{}')]:
        s = s.replace(a, b)
    s = s.replace('textbackslash\\{\\}', 'textbackslash{}')
    s = s.replace('--', '-{}-')  # keep literal double hyphens (e.g. --ultra-sensitive) from becoming en dashes
    # Arial lacks the superscript minus, so typeset superscript runs (m⁻³, 10⁴) with \textsuperscript
    sup = dict(zip('⁰¹²³⁴⁵⁶⁷⁸⁹⁻', '0123456789-'))
    s = re.sub('[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+', lambda m: r'\textsuperscript{' + ''.join(sup[c] for c in m.group()) + '}', s)
    # italic spans marked by docx_paragraphs() / italic_genera() become \textit
    return s.replace(IT_ON, r'\textit{').replace(IT_OFF, '}')


IT_ON, IT_OFF = '\x01', '\x02'


def docx_paragraphs(path):
    """Paragraph texts of a .docx, with italic runs wrapped in IT_ON/IT_OFF so genus names stay italic."""
    W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
    root = etree.fromstring(zipfile.ZipFile(path).read('word/document.xml'))
    out = []
    for p in root.iter(W + 'p'):
        txt = ''
        for r in p.iter(W + 'r'):
            t = ''.join(x.text or '' for x in r.findall(W + 't'))
            i = r.find(W + 'rPr/' + W + 'i')
            italic = i is not None and i.get(W + 'val') not in ('0', 'false')
            txt += IT_ON + t + IT_OFF if italic and t.strip() else t
        out.append(txt.replace(IT_OFF + IT_ON, ''))
    return out


def italic_genera(s):
    """Italicize genus names in a 'Genus (n), Genus (n)' list; clades, placeholders and (Family) entries stay roman."""
    def one(item):
        m = re.match(r'(Candidatus )?([A-Z][a-z]+(?:-[A-Z][a-z]+)*)( \(\d+\))$', item)
        if not m:
            return item
        if m.group(1):  # Candidatus italic, the provisional name upright (as in the main text)
            return IT_ON + 'Candidatus' + IT_OFF + ' ' + m.group(2) + m.group(3)
        return IT_ON + m.group(2) + IT_OFF + m.group(3)
    return ', '.join(one(x) for x in s.split(', '))


# Figure legends: one paragraph per figure, starting "Supplementary Figure N."
legends = {}
for par in docx_paragraphs(os.path.join(HERE, 'Supplementary_Figure_Legends_21OCT26.docx')):
    m = re.match(r'Supplementary Figure (\d+)\.\s+(.*)', par.strip(), re.S)
    if m:
        legends[int(m.group(1))] = m.group(2).strip()

# Current figure files, one per number
figs = {}
for f in glob.glob(os.path.join(HERE, 'Supplementary_Figure_*_*.pdf')):
    if 'superseded' in f:
        continue
    n = int(re.search(r'Supplementary_Figure_(\d+)_', f).group(1))
    assert n not in figs, f'two current files for Supplementary Figure {n}'
    figs[n] = os.path.basename(f)
assert sorted(figs) == sorted(legends), f'figures {sorted(figs)} vs legends {sorted(legends)}'

# Supplementary Table 1 summary sheet: title, three note rows, blank, header, group rows, total
ws = openpyxl.load_workbook(os.path.join(HERE, 'Supplementary_Table_1_patagonia_shared_ASVs.xlsx'))['Summary']
rows = [r for r in ws.iter_rows(values_only=True)]
t1_title = re.sub(r'^Supplementary Table 1\.\s*', '', rows[0][0])
t1_notes = [r[0] for r in rows[1:4]]
t1_header = rows[5]
t1_body = rows[6:]

# Supplementary Table 2 (PERMDISP + pairwise PERMANOVA): title, three note rows, blank, header, body
ws2 = openpyxl.load_workbook(os.path.join(HERE, 'Supplementary_Table_2_permdisp_pairwise.xlsx'))['Summary']
rows2 = [r for r in ws2.iter_rows(values_only=True)]
t2_title = re.sub(r'^Supplementary Table 2\.\s*', '', rows2[0][0])
t2_notes = [r[0] for r in rows2[1:4]]
t2_header = rows2[5]
t2_body = rows2[6:]


def cell(v):
    if v is None:
        return ''
    if isinstance(v, float):
        return f'{v:.2f}'
    return tex(str(v))


def fig_width(name):
    """Full text width, except small figures, which are enlarged at most 1.3x so their fonts stay in proportion."""
    out = subprocess.run(['pdfinfo', os.path.join(HERE, name)], capture_output=True, text=True, check=True).stdout
    w = float(re.search(r'Page size:\s+([\d.]+)', out).group(1))
    return rf'\textwidth' if w * 1.3 >= 470 else f'{w * 1.3:.0f}pt'  # 470 pt = A4 text width at 2.2 cm margins


L = [r'''\documentclass[11pt]{article}
\usepackage[a4paper,margin=2.2cm]{geometry}
\usepackage{fontspec}
\setmainfont{Arial}
\usepackage{graphicx,tabularx,booktabs,array,longtable}
\usepackage[hidelinks]{hyperref}
% Caption style without the caption package (not in the local TinyTeX): small, bold label, period separator
\renewcommand{\figurename}{Fig.}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\tablename}{Table}
\renewcommand{\thetable}{S\arabic{table}}
\makeatletter
\renewcommand{\@makecaption}[2]{\vskip 8pt\small\textbf{#1.} #2\par\vskip 4pt}
\makeatother
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt}
\begin{document}
\thispagestyle{empty}
\begin{center}
{\Large\bfseries Additional file 1\par}
\vspace{0.4em}
{\large\bfseries Supplementary figures and tables\par}
\vspace{1em}
{\large ''' + tex(TITLE) + r'''\par}
\vspace{1em}
''' + tex(AUTHORS) + r'''\par
\vspace{0.8em}
Corresponding author: \href{mailto:''' + CORRESPONDING + '}{' + CORRESPONDING + r'''}
\end{center}
\vspace{2em}
\textbf{Contents}

\begin{tabularx}{\textwidth}{@{}l>{\raggedright\arraybackslash}X@{}}  % second column wraps within the margin''']
for n in sorted(figs):
    first = re.split(r'(?<=\.)\s', legends[n], maxsplit=1)[0]
    L.append(rf'Fig. S{n} & \hyperref[fig:{n}]{{{tex(first.rstrip("."))}}} \\')
L.append(rf'Table S1 & \hyperref[tab:1]{{{tex(t1_title.rstrip("."))}}} \\')
L.append(rf'Table S2 & \hyperref[tab:2]{{{tex(t2_title.rstrip("."))}}} \\')
L.append(r'''\end{tabularx}

\vspace{1em}
{\raggedright Provided as separate files: Table S1 (full table, Excel; summary on page \pageref{tab:1}),
Additional file 2 (processed amplicon data) and Additional file 3 (MIxS sample metadata). The curated marker
protein database is archived on Zenodo (https://doi.org/10.5281/zenodo.23986366).\par}
''')
for n in sorted(figs):
    L.append(rf'''\clearpage
\begin{{figure}}[p]
\centering
\includegraphics[width={fig_width(figs[n])},height=0.6\textheight,keepaspectratio]{{{figs[n]}}}
\caption{{{tex(legends[n])}}}
\label{{fig:{n}}}
\end{{figure}}''')

ncol = len(t1_header)
L.append(r'''\clearpage
\begin{table}[p]
\caption{''' + tex(t1_title) + r'''}
\label{tab:1}
\small
\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{3.1cm}rr>{\raggedright\arraybackslash}p{1.9cm}>{\raggedright\arraybackslash}p{1.9cm}>{\raggedright\arraybackslash}X@{}}
\toprule''')
L.append(' & '.join(r'\textbf{' + cell(h) + '}' for h in t1_header) + r' \\ \midrule')
for r in t1_body:
    if r[0] == 'Total':
        L.append(r'\midrule')
    vals = list(r[:ncol])
    if isinstance(vals[-1], str):
        vals[-1] = italic_genera(vals[-1])  # last column lists genera
    L.append(' & '.join(cell(v) for v in vals) + r' \\')
L.append(r'''\bottomrule
\end{tabularx}
\vspace{0.5em}

\footnotesize''')
for note in t1_notes:
    L.append(tex(note) + '\n')
L.append(r'The full table (207 ASVs, with taxonomy and per-group detection counts) is provided as a separate Excel file.')
L.append(r'''\end{table}''')

def cell2(v, j):
    if v is None:
        return ''
    if j == 7 or j == 8:                      # p and p_BH
        return '< 0.001' if v < 0.001 else f'{v:.3f}'
    if j == 6:                                # R2
        return f'{v:.3f}'
    if j == 5:                                # F
        return f'{v:.2f}'
    if j == 4:                                # n per group: allow line breaks after each slash
        return tex(str(v)).replace('/', r'/\allowbreak{}')
    return tex(str(v))

# Table 2 runs to about two pages, so it is a longtable (header repeats on each page). The caption is typeset by hand
# to match the float captions; \refstepcounter makes it Supplementary Table 2 for \label/\hyperref.
L.append(r'''\clearpage
{\small\refstepcounter{table}\label{tab:2}\textbf{Table S2.} ''' + tex(t2_title) + r'''\par}
\vspace{4pt}
{\footnotesize''')
for note in t2_notes:
    L.append(tex(note) + '\n')
L.append(r'''}
\scriptsize
\setlength{\LTpre}{2pt}\setlength{\tabcolsep}{3pt}
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{1.9cm}>{\raggedright\arraybackslash}p{1.5cm}>{\raggedright\arraybackslash}p{1.9cm}>{\raggedright\arraybackslash}p{4.0cm}>{\raggedright\arraybackslash}p{1.7cm}rrrr@{}}
\toprule''')
head2 = ' & '.join(r'\textbf{' + tex(str(h)) + '}' for h in t2_header) + r' \\ \midrule'
L.append(head2 + r'''
\endfirsthead
\toprule
''' + head2 + r'''
\endhead
\bottomrule
\endfoot''')
for r in t2_body:
    L.append(' & '.join(cell2(v, j) for j, v in enumerate(r[:len(t2_header)])) + r' \\')
L.append(r'''\end{longtable}
\end{document}''')

with open(os.path.join(HERE, OUT + '.tex'), 'w') as fh:
    fh.write('\n'.join(L) + '\n')

# Two passes so the contents page references resolve
for _ in range(2):
    subprocess.run(['xelatex', '-interaction=nonstopmode', '-halt-on-error', OUT + '.tex'],
                   cwd=HERE, check=True, capture_output=True)
for ext in ('aux', 'log', 'out'):
    p = os.path.join(HERE, f'{OUT}.{ext}')
    if os.path.exists(p):
        os.remove(p)
print('wrote', OUT + '.pdf', '|', len(figs), 'figures + Supplementary Tables 1-2 summaries')
