#!/usr/bin/env python3
# figureS_study_sites_17SEP26.py
# Supplementary Figure 1. Map of the 1,500 km Patagonia to Antarctic Peninsula transect with the five Patagonian
# sites and the two Antarctic bases. Colors are read from 16S/palettes.R.
# Output : 16S/figures/figureS_study_sites_17SEP26.{pdf,png}
# Run    : python3 16S/figureS_study_sites_17SEP26.py

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT = Path(__file__).resolve().parent.parent          # repository root
OUT = ROOT / "16S" / "figures" / "figureS_study_sites_17SEP26"

# --- colors from palettes.R (single source of truth) ------------------
def read_palette(name):
    """Parse a `name <- c("key" = "#hex", ...)` block out of palettes.R."""
    txt = (ROOT / "16S" / "palettes.R").read_text()
    block = re.search(rf"{name}\s*<-\s*c\((.*?)\)", txt, re.S).group(1)
    return dict(re.findall(r'"([^"]+)"\s*=\s*"(#[0-9A-Fa-f]{6})"', block))

OCEAN, LAND, COAST = "#eef3f7", "#e8e6df", "#7f8c96"
PAT_C, ANT_C = "#a6cee3", "#1f78b4"                   # pal_region Patagonian / Antarctic
HALO = [pe.withStroke(linewidth=2.5, foreground="white")]

# --- sites: lon, lat, label, label offset (dx, dy), alignment ---------
PATAGONIA = {
    "Puerto Natales":  (-72.51, -51.73, "Puerto Natales",  0.35,  0.30, "left"),
    "San Gregorio":    (-70.13, -52.63, "San Gregorio",    0.35,  0.28, "left"),
    "Punta Arenas":    (-70.92, -53.16, "Punta Arenas",   -0.45, -0.45, "right"),
    "Porvenir":        (-70.37, -53.29, "Porvenir",        0.35, -0.30, "left"),
    "Puerto Williams": (-67.62, -54.93, "Puerto Williams", 0.45,  0.10, "left"),
}
ANTARCTIC = {
    "Risopatron": (-59.68, -62.40, "Risopatrón", 0.45,  0.20, "left"),
    "Yelcho":     (-63.58, -64.87, "Yelcho",    -0.55, -0.10, "right"),
}

# --- figure -----------------------------------------------------------
fig = plt.figure(figsize=(5.6, 7.0), dpi=300)
proj = ccrs.LambertAzimuthalEqualArea(central_longitude=-67, central_latitude=-58)
axA = plt.axes(projection=proj)
axA.set_extent([-80, -55, -66.5, -50.0], crs=ccrs.PlateCarree())
axA.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor=OCEAN, zorder=0)
axA.add_feature(cfeature.LAND.with_scale("50m"), facecolor=LAND, zorder=1)
axA.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor=COAST, linewidth=0.5, zorder=2)
axA.add_feature(cfeature.BORDERS.with_scale("50m"), edgecolor=COAST, linewidth=0.35,
                linestyle=":", zorder=2)

gl = axA.gridlines(draw_labels=True, linewidth=0.3, color="#c3ccd4", alpha=0.8, zorder=1)
gl.top_labels = gl.right_labels = False
gl.xlabel_style = gl.ylabel_style = {"size": 7, "color": COAST}

# transect line, Puerto Williams -> Risopatrón, as a distance cue only
axA.plot([-67.62, -59.68], [-54.93, -62.40], transform=ccrs.Geodetic(),
         color=COAST, linewidth=1.0, linestyle=(0, (5, 3)), zorder=3)
axA.text(-62.6, -58.7, "≈1,500 km", transform=ccrs.PlateCarree(), fontsize=8,
         style="italic", color="#43505c", ha="left", va="center",
         path_effects=HALO, zorder=6)

for group, color in [(PATAGONIA, PAT_C), (ANTARCTIC, ANT_C)]:
    for lon, lat, label, dx, dy, ha in group.values():
        axA.plot(lon, lat, "o", color=color, markersize=6.5, markeredgecolor="#2c3e50",
                 markeredgewidth=0.7, transform=ccrs.PlateCarree(), zorder=5)
        axA.text(lon + dx, lat + dy, label, transform=ccrs.PlateCarree(), fontsize=8,
                 color="#1c2b38", ha=ha, va="center", path_effects=HALO, zorder=6)

axA.text(-78.5, -50.7, "PATAGONIA", transform=ccrs.PlateCarree(), fontsize=9,
         fontweight="bold", color="#43505c", ha="left", path_effects=HALO, zorder=6)
axA.text(-59.5, -64.2, "ANTARCTIC\nPENINSULA", transform=ccrs.PlateCarree(), fontsize=9,
         fontweight="bold", color="#43505c", ha="center", va="top",
         path_effects=HALO, zorder=6)

axA.legend(handles=[
    Line2D([0], [0], marker="o", color="w", markerfacecolor=PAT_C, markersize=7,
           markeredgecolor="#2c3e50", markeredgewidth=0.7, label="Patagonia (source region)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=ANT_C, markersize=7,
           markeredgecolor="#2c3e50", markeredgewidth=0.7, label="Antarctic Peninsula (sink)"),
], loc="lower left", fontsize=7.5, frameon=True, framealpha=0.95, edgecolor="#cfd8e0")
axA.spines["geo"].set_edgecolor(COAST)
axA.spines["geo"].set_linewidth(0.8)

OUT.parent.mkdir(parents=True, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}.{ext}", facecolor="white", bbox_inches="tight", pad_inches=0.08)
print(f"wrote {OUT}.pdf / .png")
