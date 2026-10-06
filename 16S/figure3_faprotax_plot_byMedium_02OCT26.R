# figure3_faprotax_plot_byMedium_02OCT26.R
# Figure 3. FAPROTAX functional potential by medium (air, bulk soil, rhizosphere): mean relative abundance per
# function, grouped by biogeochemical cycle, with pairwise Wilcoxon tests (Benjamini-Hochberg) and letters.
# Inputs : 16S/figure3_faprotax_02OCT26.tsv (FAPROTAX output), 16S/figure3_lrt_metadata_02OCT26.csv, 16S/palettes.R,
#          16S/faprotax_function_cycle_mapping.csv (optional)
# Outputs: 16S/figures/figure3_faprotax_bubble_byMedium_02OCT26.{pdf,png} and _stats_02OCT26.txt
# Run    : Rscript 16S/figure3_faprotax_plot_byMedium_02OCT26.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(rstatix)
  library(multcompView)
})

set.seed(100)
base      <- "./16S"
date_tag  <- "02OCT26"
fapro_tag <- "02OCT26"
fig_dir   <- file.path(base, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

source(file.path(base, "palettes.R"))

fapro_path <- file.path(base, paste0("figure3_faprotax_",     fapro_tag, ".tsv"))
md_path    <- file.path(base, paste0("figure3_lrt_metadata_", fapro_tag, ".csv"))
cyc_path   <- file.path(base, "faprotax_function_cycle_mapping.csv")
stopifnot(file.exists(fapro_path), file.exists(md_path))

# --- 1) read FAPROTAX output ----------------------------------------
fapro <- read_tsv(fapro_path, show_col_types = FALSE)
if ("OTUID" %in% names(fapro)) fapro$OTUID <- NULL

fapro_long <- fapro |>
  rename(`function` = group) |>
  pivot_longer(-`function`, names_to = "sample_id", values_to = "rel_ab")

# --- 2) join metadata ------------------------------------------------
md <- read_csv(md_path, show_col_types = FALSE)
fapro_long <- fapro_long |>
  inner_join(md |> select(sample_id, medium), by = "sample_id") |>
  mutate(medium = factor(medium,
                         levels = c("Air", "Bulk soil", "Rhizosphere")))

# --- 3) cycle classification (same as the main figure) -------------
classify_cycle <- function(fn) {
  case_when(
    str_detect(fn, "methan|methylot|methanol")     ~ "Methane cycle",
    str_detect(fn, "sulf|thiosulf|sulphur")        ~ "Sulfur cycle",
    str_detect(fn, "nitr|ammon|denitrif|anammox") ~ "Nitrogen cycle",
    str_detect(fn, "phototro|photosynth|photohet|cyanob|photoautotroph|oxygenic") ~ "Phototrophy",
    str_detect(fn, "fermen|chemoheterotro|hydrocarbon|cellulol|aromat|acetogen|chitinol|xylanol|fixation") ~ "Carbon cycle",
    str_detect(fn, "pathogen|parasite|chloroplast|associated|symbiont|predator") ~ "Host-associated",
    TRUE                                           ~ "Other redox"
  )
}
cyc_existing <- if (file.exists(cyc_path)) {
  read_csv(cyc_path, show_col_types = FALSE) |>
    select(`function` = Function, Cycle)
} else tibble(`function` = character(0), Cycle = character(0))

cycle_df <- tibble(`function` = unique(fapro_long$`function`)) |>
  left_join(cyc_existing, by = "function") |>
  mutate(Cycle = if_else(is.na(Cycle), classify_cycle(`function`), Cycle))

fapro_long <- fapro_long |>
  left_join(cycle_df, by = "function") |>
  filter(Cycle != "Host-associated")

# --- 4) function filter (same as main figure) -----------------------
per_medium_mean <- fapro_long |>
  group_by(`function`, medium) |>
  summarize(mean_pct = mean(rel_ab, na.rm = TRUE) * 100, .groups = "drop")

keep_functions <- per_medium_mean |>
  group_by(`function`) |>
  summarize(max_mean = max(mean_pct), .groups = "drop") |>
  filter(max_mean >= 0.5) |>
  pull(`function`)

focal_functions <- c(
  "methylotrophy", "methanol_oxidation", "methanotrophy",
  "methanogenesis", "hydrogenotrophic_methanogenesis",
  "sulfate_respiration", "sulfur_respiration",
  "dark_sulfite_oxidation", "sulfite_respiration",
  "thiosulfate_respiration", "respiration_of_sulfur_compounds",
  "dark_sulfide_oxidation", "dark_thiosulfate_oxidation",
  "dark_oxidation_of_sulfur_compounds"
)
keep_functions <- intersect(union(keep_functions, focal_functions),
                            unique(fapro_long$`function`))

fapro_keep <- fapro_long |> filter(`function` %in% keep_functions)
cat("Functions kept: ", length(keep_functions),
    " (of ", length(unique(fapro_long$`function`)), ")\n", sep = "")

# --- 5) aggregate to medium (drop region) ---------------------------
group_mean <- fapro_keep |>
  group_by(`function`, Cycle, medium) |>
  summarize(mean_pct = mean(rel_ab, na.rm = TRUE) * 100,
            n        = sum(!is.na(rel_ab)),
            .groups  = "drop") |>
  filter(n > 0) |>
  mutate(group_label = factor(medium,
                              levels = c("Air", "Bulk soil", "Rhizosphere")))

# --- 6) function ordering -------------------------------------------
cycle_order_top_down <- c("Methane cycle", "Sulfur cycle", "Phototrophy",
                          "Nitrogen cycle", "Carbon cycle", "Other redox")
cycle_order_top_down <- intersect(cycle_order_top_down,
                                  unique(group_mean$Cycle))
cycle_order_bottom_up <- rev(cycle_order_top_down)

func_order <- group_mean |>
  group_by(`function`, Cycle) |>
  summarize(top = max(mean_pct), .groups = "drop") |>
  mutate(Cycle = factor(Cycle, levels = cycle_order_bottom_up)) |>
  arrange(Cycle, top) |>
  pull(`function`)

group_mean <- group_mean |>
  mutate(`function` = factor(`function`, levels = func_order),
         Cycle      = factor(Cycle, levels = cycle_order_top_down))

# --- 6b) per-function pairwise Wilcoxon -> letters ------------------
# Test rel_ab ~ medium for each function; BH-adjusted; compact-letter
# display via multcompView. Same machinery as Figure 2.
make_letters <- function(sub) {
  sub <- sub |> mutate(medium = droplevels(factor(medium)))
  mediums <- levels(sub$medium)
  if (length(mediums) < 2 || all(sub$rel_ab == 0, na.rm = TRUE)) {
    return(tibble(medium = mediums, cld = rep("a", length(mediums))))
  }
  pw <- tryCatch(
    rstatix::pairwise_wilcox_test(sub, rel_ab ~ medium,
                                  p.adjust.method = "BH"),
    error = function(e) NULL
  )
  if (is.null(pw) || nrow(pw) == 0) {
    return(tibble(medium = mediums, cld = rep("a", length(mediums))))
  }
  pmat <- with(pw, setNames(p.adj, paste(group1, group2, sep = "-")))
  cld  <- multcompView::multcompLetters(pmat)$Letters
  tibble(medium = names(cld), cld = unname(cld))
}

letters_df <- fapro_keep |>
  group_by(`function`) |>
  group_modify(~ make_letters(.x)) |>
  ungroup() |>
  mutate(medium = factor(medium,
                         levels = c("Air", "Bulk soil", "Rhizosphere")))

group_mean <- group_mean |>
  left_join(letters_df, by = c("function", "medium"))

# --- 7) two-panel split ---------------------------------------------
focal_cycles <- intersect(c("Methane cycle", "Sulfur cycle", "Nitrogen cycle"),
                          cycle_order_top_down)
# "Other redox" reduces to ureolysis at LRT scale: drop it.
bulk_cycles  <- setdiff(cycle_order_top_down,
                        c(focal_cycles, "Other redox"))

plotdat <- group_mean |> filter(mean_pct > 0)
plotdat_focal <- plotdat |>
  filter(Cycle %in% focal_cycles) |>
  mutate(Cycle = factor(Cycle, levels = focal_cycles))
plotdat_bulk  <- plotdat |>
  filter(Cycle %in% bulk_cycles) |>
  mutate(Cycle = factor(Cycle, levels = bulk_cycles))

base_theme <- theme_bw(base_size = 8) +
  theme(
    axis.text.x        = element_text(angle = 30, hjust = 1),
    panel.grid.major   = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.spacing.y    = unit(0.15, "lines"),
    strip.placement    = "outside",
    strip.background.y = element_rect(fill = "gray92", color = NA),
    strip.text.y.left  = element_text(angle = 0, face = "bold", hjust = 0),
    legend.position    = "right",
    strip.text.y       = element_text(size = 7),
    legend.title       = element_text(size = 7),
    legend.text        = element_text(size = 7)
  )

build_panel <- function(dat, size_breaks, size_labels, max_size,
                        size_limits, size_name, drop_x = FALSE) {
  p <- ggplot(dat,
              aes(x = group_label, y = `function`,
                  size = mean_pct, color = medium)) +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = cld),
              position = position_nudge(y = 0.32),
              size = 2.1, color = "gray20",
              fontface = "bold",
              show.legend = FALSE) +
    scale_size_area(max_size = max_size,
                    breaks   = size_breaks,
                    labels   = size_labels,
                    limits   = size_limits,
                    oob      = scales::squish,
                    name     = size_name) +
    scale_color_manual(values = pal_medium3, name = "Medium") +
    facet_grid(Cycle ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    # Readable function names for display only (underscores to spaces).
    scale_y_discrete(labels = function(x) str_replace_all(x, "_", " ")) +
    labs(x = NULL, y = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3), order = 1)) +
    base_theme
  if (drop_x) {
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  }
  p
}

pA <- build_panel(
  dat         = plotdat_focal,
  size_breaks = c(0.1, 0.5, 1, 2, 5),
  size_labels = c("0.1%", "0.5%", "1%", "2%", "5%"),
  max_size    = 8,
  size_limits = c(0, 5),
  size_name   = "Mean relative\nabundance (a)",
  drop_x = TRUE
)

pB <- build_panel(
  dat         = plotdat_bulk,
  size_breaks = c(1, 5, 10, 20, 40),
  size_labels = c("1%", "5%", "10%", "20%", "40%"),
  max_size    = 8,
  size_limits = c(0, 40),
  size_name   = "Mean relative\nabundance (b)"
)

n_focal <- length(unique(plotdat_focal$`function`))
n_bulk  <- length(unique(plotdat_bulk$`function`))

fig <- (pA / pB) +
  plot_layout(heights = c(n_focal, n_bulk), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10),
        legend.position = "right")

# Three-column variant: narrower than the by-region figure.
w <- 120                                            # mm
h <- min(240, 5.2 * (n_focal + n_bulk) + 25)        # mm

pdf_path <- file.path(fig_dir, paste0("figure3_faprotax_bubble_byMedium_", date_tag, ".pdf"))
png_path <- file.path(fig_dir, paste0("figure3_faprotax_bubble_byMedium_", date_tag, ".png"))
ggsave(pdf_path, fig, width = w, height = h, units = "mm")
ggsave(png_path, fig, width = w, height = h, units = "mm", dpi = 300)
cat("Saved:\n  ", pdf_path, "\n  ", png_path, "\n", sep = "")

# --- 8) stats output -------------------------------------------------
stats_path <- file.path(fig_dir,
                        paste0("figure3_faprotax_bubble_byMedium_stats_",
                               date_tag, ".txt"))
sink(stats_path)
cat("Figure 3 stats - FAPROTAX functional potential, LRT only (by medium)\n")
cat("Date: ", date_tag, "  | FAPROTAX run: ", fapro_tag, "\n", sep = "")
cat("Input: ", basename(fapro_path), "\n\n", sep = "")

cat("Sample counts by medium:\n")
print(table(md$medium))

cat("\nFunctions retained (",
    length(keep_functions), "):\n", sep = "")
print(sort(keep_functions))

cat("\nMedium-level means (mean rel-ab %, by function x medium) with letters:\n")
print(group_mean |>
        arrange(Cycle, `function`, medium) |>
        select(Cycle, `function`, medium, n, mean_pct, cld),
      n = Inf)

cat("\nPer-function pairwise Wilcoxon (BH-adjusted) across medium:\n")
pw_all <- fapro_keep |>
  group_by(`function`) |>
  group_modify(~ {
    sub <- .x |> mutate(medium = droplevels(factor(medium)))
    if (length(unique(sub$medium)) < 2) return(tibble())
    rstatix::pairwise_wilcox_test(sub, rel_ab ~ medium,
                                  p.adjust.method = "BH")
  }) |>
  ungroup() |>
  left_join(cycle_df, by = "function") |>
  arrange(Cycle, `function`, group1, group2)
print(pw_all, n = Inf)
sink()
cat("Stats: ", stats_path, "\n", sep = "")
