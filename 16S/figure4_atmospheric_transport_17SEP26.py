#!/usr/bin/env python3
# figure4_atmospheric_transport_17SEP26.py
# Figure 4. Assembles two vector panels produced by the atmospheric modeling: (a) HYSPLIT trajectories driven by
# Polar-WRF and (b) CHIMERE pISOPA1 surface concentration with site time series. Panels are placed without
# rasterization; only panel letters are added.
# Inputs : New_Hysplit/hysplit_ar_transport_dual_receptor.pdf, New_Chimere/figure_chimere1.pdf
# Output : 16S/figures/figure4_atmospheric_transport_17SEP26.{pdf,png}
# Run    : python3 16S/figure4_atmospheric_transport_17SEP26.py

import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject

ROOT = Path(__file__).resolve().parent.parent            # repository root
OUT = ROOT / "16S" / "figures" / "figure4_atmospheric_transport_17SEP26"
PANEL_A = ROOT / "New_Hysplit" / "hysplit_ar_transport_dual_receptor.pdf"
PANEL_B = ROOT / "New_Chimere" / "figure_chimere1.pdf"

WIDTH = 510.0        # pt, ~180 mm double-column width
GAP = 6.0            # vertical gap between panels
LAB_H = 13.0         # strip reserved above each panel for its letter
PAD = 4.0            # outer margin

# CHIMERE carries wide white margins; trim them so the two panels balance.
CROP_B = dict(left=0.055, right=0.015, top=0.0, bottom=0.05)   # fractions


def placed(path, crop=None):
    """Read a one-page PDF and return (page, width, height) after optional crop."""
    page = PdfReader(str(path)).pages[0]
    box = page.mediabox
    x0, y0 = float(box.left), float(box.bottom)
    w, h = float(box.width), float(box.height)
    if crop:
        x0 += w * crop["left"]
        y0 += h * crop["bottom"]
        w -= w * (crop["left"] + crop["right"])
        h -= h * (crop["top"] + crop["bottom"])
        page.mediabox = RectangleObject((x0, y0, x0 + w, y0 + h))
        page.cropbox = RectangleObject((x0, y0, x0 + w, y0 + h))
    return page, x0, y0, w, h


page_a, ax0, ay0, aw, ah = placed(PANEL_A)
page_b, bx0, by0, bw, bh = placed(PANEL_B, CROP_B)

scale_a = (WIDTH - 2 * PAD) / aw
scale_b = (WIDTH - 2 * PAD) / bw
h_a, h_b = ah * scale_a, bh * scale_b
height = PAD + h_b + LAB_H + GAP + h_a + LAB_H + PAD

canvas = PdfWriter().add_blank_page(width=WIDTH, height=height)

# panel a on top, panel b below (y grows upward)
y_a = PAD + h_b + LAB_H + GAP
y_b = PAD
canvas.merge_transformed_page(
    page_a, Transformation().scale(scale_a).translate(PAD - ax0 * scale_a, y_a - ay0 * scale_a))
canvas.merge_transformed_page(
    page_b, Transformation().scale(scale_b).translate(PAD - bx0 * scale_b, y_b - by0 * scale_b))

# panel letters as vector text on a transparent overlay
fig = plt.figure(figsize=(WIDTH / 72, height / 72))
fig.patch.set_alpha(0)
for label, y in [("a", y_a + h_a + 1.5), ("b", y_b + h_b + 1.5)]:
    fig.text(PAD / WIDTH, y / height, label, fontsize=11, fontweight="bold",
             ha="left", va="bottom")
overlay_path = OUT.parent / "_fig4_labels_tmp.pdf"
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(overlay_path, transparent=True)
plt.close(fig)
canvas.merge_page(PdfReader(str(overlay_path)).pages[0])

writer = PdfWriter()
writer.add_page(canvas)
with open(f"{OUT}.pdf", "wb") as fh:
    writer.write(fh)
overlay_path.unlink()

subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile",
                f"{OUT}.pdf", str(OUT)], check=True)
print(f"wrote {OUT}.pdf / .png  ({WIDTH:.0f} x {height:.0f} pt)")
