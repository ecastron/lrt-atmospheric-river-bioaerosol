#!/usr/bin/env python3
# figureS5_synoptic_02OCT26.py
# Supplementary Figure 5. Combines the Risopatron wind roses by phase (panel a, wind_rose/wind_rose_phases_02OCT26.R)
# with the ERA5 integrated vapor transport and 500 hPa geopotential maps (panel b). Output width 183 mm at 300 dpi.
# Inputs : 16S/figures/figureS5a_wind_roses_openair_02OCT26.png, 16S/figures/figureS_atmospheric_synoptic_17SEP26.png
# Output : 16S/figures/figureS5_synoptic_02OCT26.{png,pdf}
# Run    : python3 16S/figureS5_synoptic_02OCT26.py

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

FIG = Path(__file__).resolve().parent / "figures"
A_SRC = FIG / "figureS5a_wind_roses_openair_02OCT26.png"  # openair version
B_SRC = FIG / "figureS_atmospheric_synoptic_17SEP26.png"
OUT = FIG / "figureS5_synoptic_02OCT26"
B_TOP = 1000  # first row below the panel a area of the source composite (blank band at rows 995-1054 of the 300 dpi composite)

a = Image.open(A_SRC).convert("RGB")
# trim the blank band ggsave leaves under the roses
import numpy as np
a = a.resize((2161, round(a.height * 2161 / a.width)), Image.LANCZOS)
rows = np.where((np.asarray(a) < 240).any(axis=2).any(axis=1))[0]
a = a.crop((0, 0, a.width, min(a.height, rows.max() + 20)))
b_full = Image.open(B_SRC).convert("RGB")
b = b_full.crop((0, B_TOP, b_full.width, b_full.height))
b = b.resize((a.width, round(b.height * a.width / b.width)), Image.LANCZOS)

pad = 40
canvas = Image.new("RGB", (a.width, a.height + pad + b.height), "white")
canvas.paste(a, (0, 0))
canvas.paste(b, (0, a.height + pad))

# panel tag for a, matching the bold "B" baked into panel b
draw = ImageDraw.Draw(canvas)
try:
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 60)
except OSError:
    font = ImageFont.load_default()
draw.text((20, 20), "A", fill="black", font=font)

canvas.save(f"{OUT}.png", dpi=(300, 300))
import subprocess
subprocess.run(["sips", "-s", "format", "pdf", f"{OUT}.png", "--out", f"{OUT}.pdf"], check=True, capture_output=True)
print(f"wrote {OUT}.png / .pdf  size {canvas.size}")
