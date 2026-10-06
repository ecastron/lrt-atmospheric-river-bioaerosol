# palettes.R
# Color and shape palettes shared by all figure scripts; sourced at the start of each figure script.

# --- source2: 9 sample classes (LRT + Malard 2022) -------------------
# Air = blue family, Soil = green family, Precipitation = brown,
# SO ocean-side = orange family, Control = gray.
pal_source2 <- c(
  "Antarctic Air"    = "#1f78b4",   # mid blue
  "Patagonian Air"   = "#a6cee3",   # light blue
  "SO Air"           = "#ffcc80",   # light orange
  "SO Air over land" = "#ff8c00",   # dark orange
  "SO Precipitation" = "#8c510a",   # brown
  "Antarctic Soil"   = "#33a02c",   # mid green
  "Patagonian Soil"  = "#b2df8a",   # light green
  "Control"          = "#999999"    # gray
)

# --- medium: 3 sample compartments -----------------------------------
pal_medium <- c(
  "Air"           = "#1f78b4",
  "Soil"          = "#33a02c",
  "Precipitation" = "#8c510a"
)

# --- region: 3 geographic blocks -------------------------------------
pal_region <- c(
  "Antarctic"  = "#1f78b4",
  "Patagonian" = "#a6cee3",
  "SO"         = "#ff8c00"
)

# --- dataset: shapes (color is reserved for biology) -----------------
shp_dataset <- c(
  "LRT"        = 16,    # filled circle
  "Malard2022" = 17     # filled triangle
)

# --- latitude: continuous sequential --------------------------------
# Use as: scale_color_viridis_c(option = pal_latitude_opt,
#                               direction = pal_latitude_dir)
pal_latitude_opt <- "C"
pal_latitude_dir <- -1

# --- phylum palette ------------------
# Top 8 dominant phyla in the LRT data + an "Other" catch-all. Colors
# chosen for contrast (mustard/orange
# Proteobacteria, light-blue Actinobacteriota, navy Acidobacteriota,
# rust Verrucomicrobiota, etc.).
pal_phylum <- c(
  "Proteobacteria"    = "#E59B26",   # mustard/orange
  "Actinobacteriota"  = "#7BB6DA",   # light blue
  "Bacteroidota"      = "#2E8B57",   # green
  "Firmicutes"        = "#F4D03F",   # yellow
  "Acidobacteriota"   = "#2874A6",   # dark blue
  "Verrucomicrobiota" = "#D35400",   # rust orange
  "Cyanobacteria"     = "#E91E63",   # pink
  "Chloroflexi"       = "#6E7B7F",   # dark gray
  "Other"             = "#BDC3C7"    # light gray
)

# --- medium palette: 3-way for Figure 2 -----------------------------
# (re-uses pal_medium for Air/Soil and adds Rhizosphere; works alongside
# pal_medium above for figures that lump bulk-soil and rhizosphere)
pal_medium3 <- c(
  "Air"         = "#1f78b4",
  "Bulk soil"   = "#33a02c",
  "Rhizosphere" = "#b15928"
)

# --- site (region) palette, order, and display labels ----------------
# Used for site-level facets in LRT figures. Order is geographic
# (Patagonia north -> south, then Antarctic Peninsula north -> south).
# `site_labels` gives compact display strings: use as_labeller(site_labels)
# in facet calls so panel strips don't get truncated.
site_order <- c("PA_POR", "PN_SG", "Puerto Williams",
                "Risopatron", "Yelcho")
site_labels <- c(
  "PA_POR"          = "PA_POR",
  "PN_SG"           = "PN_SG",
  "Puerto Williams" = "PW",
  "Risopatron"      = "Risopatron",
  "Yelcho"          = "Yelcho"
)
pal_site <- c(
  "PA_POR"          = "#a6cee3",
  "PN_SG"           = "#1f78b4",
  "Puerto Williams" = "#6a3d9a",
  "Risopatron"      = "#33a02c",
  "Yelcho"          = "#b15928"
)

# --- AR phase: 3-level Before / During / After -----------------------
# Used in Figure 4 (Risopatron AR event). "During" is the loud color
# because it's the focal window; Before/After are cool/neutral.
pal_ar_phase <- c(
  "Before" = "#4C72B0",   # steel blue
  "During" = "#C44E52",   # brick red (AR window)
  "After"  = "#8C8C8C"    # slate gray
)
pal_atm_river <- c("yes" = "#C44E52", "no" = "#4C72B0")

# AR window dates (UTC). "Strict" matches the ERA5 IVT-defined window
# (Feb 7-8). "Affected" matches the metadata
# atm_river=="yes" flag (Feb 5-11) and is used for visual shading on
# time-series panels.
ar_event_start    <- as.Date("2022-02-07")
ar_event_end      <- as.Date("2022-02-08")
ar_affected_start <- as.Date("2022-02-05")
ar_affected_end   <- as.Date("2022-02-11")

# --- fungal source habitat (Supp Fig S1) -----------------------------
# Coarse origin class assigned from FungalTraits (Polme et al. 2020):
# Aquatic_habitat marine -> Marine; freshwater/aquatic -> Freshwater;
# land lifestyles (soil/litter/wood/dung saprotroph, plant pathogen,
# mycorrhizal, lichenized, endophyte ...) -> Terrestrial; else Unknown.
# Green = land, blue = sea, teal = freshwater, gray = unresolved.
pal_fungal_origin <- c(
  "Terrestrial" = "#2E8B57",   # sea green (land source = Patagonia)
  "Marine"      = "#3B6FB6",   # blue (Southern Ocean source)
  "Freshwater"  = "#5AB4C5",   # teal
  "Unknown"     = "#BDBDBD"    # gray
)

# --- fungal mean abundance ramp (Supp Fig S1 source-attribution) ------
# Color encodes log10 mean relative abundance; bubble size encodes
# prevalence. Use as: scale_fill_viridis_c(option = pal_fungal_abund_opt,
#                                           direction = pal_fungal_abund_dir,
#                                           trans = "log10")
# mako (blue-green), reversed so high abundance reads dark/heavy. Chosen
# distinct from the green prevalence ramp it replaces, so size and color
# are not mistaken for the same signal.
pal_fungal_abund_opt <- "mako"
pal_fungal_abund_dir <- -1

# Per-class strip highlight for Supp Fig S1: draw the eye to the
# long-range transport band; other bands stay neutral gray.
pal_s1_strip <- c(
  "Long-range transport candidate" = "#FCE3B8",   # warm amber (signal)
  "Ambiguous / cosmopolitan"        = "#ECECEC",   # neutral gray
  "Local Antarctic background"      = "#ECECEC"     # neutral gray
)

# --- FEAST source pools (Figure 6). Patagonia = red family, Southern Ocean reference = blue,
# local Antarctic soil = brown, Unknown = gray.
pal_feast_source <- c(
  "Patagonia_Air"         = "#C44E52",
  "Patagonia_Terrestrial" = "#8C2A2A",
  "Malard_SO_Air"         = "#4C72B0",
  "Risopatron_Local"      = "#937860",
  "Unknown"               = "#BDBDBD"
)

# --- CoxL tree tip classes (Supp Fig 3) ---
# Okabe-Ito colors. Classes follow coxL_reanalysis_01OCT26/07b_classify_v2.R.
pal_coxL <- c(
  # form I = AYXCSFR (UniProt P19919/P19913), form II = AYXGAGR
  "Form I reference (AYXCSFR)"                       = "#000000",
  "Form II reference (AYXGAGR)"                      = "#999999",
  "Mo-hydroxylase outgroup"                          = "#CCCCCC",
  "Air gene, form I clade, motif not confirmed"      = "#0072B2",
  "Air gene, form I clade, weak support"             = "#56B4E9",
  "Air gene, form II"                                = "#D55E00",
  "Air gene, not CoxL (other Mo-hydroxylase)"        = "#E69F00"
)

# --- HydDB [NiFe]-hydrogenase group classes (Supp Fig 4)
pal_hyddb <- c(
  "Group 1l (atmospheric, high-affinity)"     = "#D55E00",
  "Groups 1a-1g (uptake and other)"           = "#E69F00",
  "Group 3 (bidirectional, cofactor-coupled)" = "#56B4E9",
  "Group 4 (energy-converting, fermentative)" = "#999999"
)

# --- negative-control types (Supp Fig 8a). Hex values are the ggplot2
# default hues for three levels.
pal_ctrl_type <- c(
  "buffer_control"     = "#F8766D",
  "coriolis_control"   = "#00BA38",
  "deposition_control" = "#619CFF"
)

# --- betaNTI phase pairs (Supp Fig 8b). Before-Before and Before-D+A reuse the AR
# phase blue/red.
pal_bnti_pair <- c(
  "Before-Before" = "#4C72B0",
  "Before-D+A"    = "#C44E52",
  "D+A-D+A"       = "#55A868"
)

# --- Stegen et al. 2013 process partition (Figure 5e) -------
# Selection in blues, dispersal in green/orange, drift in neutral gray
# (Okabe-Ito hues, colorblind-safe).
pal_stegen_process <- c(
  "Homogeneous selection"  = "#0072B2",
  "Variable selection"     = "#56B4E9",
  "Homogenizing dispersal" = "#009E73",
  "Dispersal limitation"   = "#E69F00",
  "Drift (undominated)"    = "#BBBBBB"
)
