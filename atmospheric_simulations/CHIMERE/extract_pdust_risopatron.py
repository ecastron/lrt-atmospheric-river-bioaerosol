# extract_pdust_risopatron.py
# Recovers the CHIMERE surface mineral dust (pDUST) series at Risopatron from the vector time-series plot of the
# simulation (ts_Risopatron_pDUST.pdf, Matplotlib 3.3.2) used for Supplementary Figure 7.
# Each plotted line is stored as a path of vertices; the axis ticks give the mapping from page coordinates to
# date (x) and concentration (y). Matplotlib's path simplification drops vertices that do not change the drawn line,
# so the series is an irregular subset of the hourly model output that reproduces the plotted curve, peaks included.
# Output: pdust_risopatron_from_plot.csv (domain, time_utc, pdust_ug_m3)
# Run   : python3 atmospheric_simulations/CHIMERE/extract_pdust_risopatron.py

import csv, os, re, zlib
from datetime import datetime, timedelta

HERE = os.path.dirname(os.path.abspath(__file__))
pdf = open(os.path.join(HERE, 'ts_Risopatron_pDUST.pdf'), 'rb').read()
stream = zlib.decompress(re.search(rb'stream\r?\n(.*?)endstream', pdf, re.S).group(1)).decode('latin1')

# x axis: tick marks (short vertical segments below the axes) followed by their "MM-DD" labels
xt = re.findall(r'([\d.]+) 55\.250166 m\n[\d.]+ 51\.750166 l\n\nB\n.*?\((\d\d-\d\d)\) Tj', stream, re.S)
xt = [(float(x), datetime.strptime('2022-' + lab, '%Y-%m-%d')) for x, lab in xt]
# y axis: tick marks (short horizontal segments left of the axes) followed by their numeric labels
yt = re.findall(r'34\.23125 ([\d.]+) m\n30\.73125 [\d.]+ l\n\nB\n.*?\(([\d.]+)\) Tj', stream, re.S)
yt = [(float(y), float(v)) for y, v in yt]
assert len(xt) >= 2 and len(yt) >= 2, 'axis ticks not found'

# Linear calibration from the first and last tick on each axis; check every intermediate tick fits
(x0, d0), (x1, d1) = xt[0], xt[-1]
pt_per_day = (x1 - x0) / ((d1 - d0).total_seconds() / 86400)
(y0, v0), (y1, v1) = yt[0], yt[-1]
val_per_pt = (v1 - v0) / (y1 - y0)
for x, d in xt:
    assert abs(x0 + (d - d0).total_seconds() / 86400 * pt_per_day - x) < 0.01, 'x ticks not linear'
for y, v in yt:
    assert abs(v0 + (y - y0) * val_per_pt - v) < 1e-3, 'y ticks not linear'

# Data lines: the three clipped paths drawn inside the axes, in legend order d01, d02, d03
clip = '34.23125 55.2501656841 223.2 166.32 re W n'
paths = []
for chunk in stream.split(clip)[1:]:
    verts = re.findall(r'([\d.]+) ([\d.]+) [ml]\n', chunk.split('\nS\n')[0])
    if len(verts) > 10:
        paths.append(verts)
labels = re.findall(r'\((d0\d)\) Tj', stream)
assert len(paths) == len(labels) == 3, (len(paths), labels)

rows = []
for dom, verts in zip(labels, paths):
    seen = set()
    for xs, ys in verts:
        t = d0 + timedelta(days=(float(xs) - x0) / pt_per_day)
        t = t.replace(second=0, microsecond=0) + timedelta(minutes=round(t.minute / 60) * 60 - t.minute)  # hourly model steps
        if t in seen:  # Matplotlib repeats the closing vertex
            continue
        seen.add(t)
        rows.append((dom, t.strftime('%Y-%m-%d %H:%M'), round(v0 + (float(ys) - y0) * val_per_pt, 4)))

with open(os.path.join(HERE, 'pdust_risopatron_from_plot.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['domain', 'time_utc', 'pdust_ug_m3'])
    w.writerows(rows)

for dom in labels:
    r = [x for x in rows if x[0] == dom]
    peak = max(r, key=lambda x: x[2])
    pre = max(x[2] for x in r if x[1] < '2022-02-06')
    print(f'{dom}: {len(r)} points, {r[0][1]} to {r[-1][1]}, peak {peak[2]:.2f} at {peak[1]}, max before 6 Feb {pre:.2f}')
