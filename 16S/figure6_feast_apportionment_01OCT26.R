# figure6_feast_apportionment_01OCT26.R
# Figure 6. FEAST source apportionment at Risopatron across the atmospheric river: per-sample source proportions
# through time (a) and proportions by phase (Before versus During+After) with Benjamini-Hochberg-adjusted Wilcoxon
# p-values (b).
# Inputs : 16S/figures/feast_AR_event_01OCT26_proportions.tsv (aggregated FEAST output from
#          feast_AR_event_postprocess_noorg_01OCT26.R), 16S/ps_merged_withTree_noorg_01OCT26.RDS (sample dates)
# Outputs: 16S/figures/figure6_feast_apportionment_01OCT26.{pdf,png}
# Run    : Rscript 16S/figure6_feast_apportionment_01OCT26.R

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(patchwork); library(scales)
})

setwd("./16S")
date_tag <- "01OCT26"   # input tag: organelle-free FEAST run
out_tag  <- "01OCT26"
source("palettes.R")
fig_dir  <- "figures"
in_tsv   <- file.path(fig_dir,
                      paste0("feast_AR_event_", date_tag, "_proportions.tsv"))

prop <- read_tsv(in_tsv, show_col_types = FALSE)

## Attach date + phase from the phyloseq sample metadata
ps <- readRDS("ps_merged_withTree_noorg_01OCT26.RDS")
sd <- as(sample_data(ps), "data.frame")
prop$date <- as.Date(sd$date[match(prop$sample_id, rownames(sd))])
prop <- prop |>
  mutate(phase = case_when(
           date <= as.Date("2022-02-06")                                ~ "Before",
           date >= as.Date("2022-02-07") & date <= as.Date("2022-02-08") ~ "During",
           TRUE                                                         ~ "After"),
         phase = factor(phase, c("Before","During","After")),
         ar2   = factor(if_else(phase == "Before", "Before", "During+After"),
                        levels = c("Before","During+After"))) |>
  arrange(date)

pool_levels <- c("Patagonia_Air", "Patagonia_Terrestrial",
                 "Malard_SO_Air", "Risopatron_Local", "Unknown")
## Single-line labels for Panel A legend
pool_labels <- c("Patagonia_Air"         = "Patagonia air",
                 "Patagonia_Terrestrial" = "Patagonia terrestrial",
                 "Malard_SO_Air"         = "Southern Ocean air (Malard 2022)",
                 "Risopatron_Local"      = "Risopatrón local soil/rhizosphere",
                 "Unknown"               = "Unknown")
## Two-line labels for Panel B facet strips (so they're readable)
pool_labels_facet <- c(
  "Patagonia_Air"         = "Patagonia\nair",
  "Patagonia_Terrestrial" = "Patagonia\nterrestrial",
  "Malard_SO_Air"         = "Southern\nOcean air",
  "Risopatron_Local"      = "Risopatrón\nlocal soil",
  "Unknown"               = "Unknown")

pal_src <- pal_feast_source   # from palettes.R

long <- prop |>
  pivot_longer(cols = all_of(pool_levels),
               names_to = "source", values_to = "proportion") |>
  mutate(source = factor(source, levels = pool_levels))

## ---- Panel A: stacked area time series ------------------------------------
## AR dates (ar_event_*, ar_affected_*) come from palettes.R
## Shade only the ERA5-defined AR core (7-8 Feb), half a day either side so
## both sample days sit inside the band (one AR definition throughout).
ar_rect <- annotate("rect",
                    xmin = ar_event_start - 0.5, xmax = ar_event_end + 0.5,
                    ymin = -Inf, ymax = Inf,
                    fill = "#C44E52", alpha = 0.10)

pA <- ggplot(long, aes(date, proportion, fill = source)) +
  geom_area(color = "white", linewidth = 0.25) +
  ## AR overlays drawn on top of the areas so they stay visible
  annotate("rect", xmin = ar_event_start - 0.5, xmax = ar_event_end + 0.5,
           ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.15) +
  scale_fill_manual(values = pal_src, labels = pool_labels, name = "Source") +
  scale_x_date(date_breaks = "3 days", date_labels = "%b %d",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Source proportion") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.key.size  = unit(3.5, "mm"),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 7.5))

## ---- Panel B: per-source boxplots by phase, with Wilcoxon p ---------------
wilcox_tbl <- long |>
  group_by(source) |>
  summarize(
    medB = median(proportion[ar2 == "Before"],       na.rm = TRUE),
    medA = median(proportion[ar2 == "During+After"], na.rm = TRUE),
    p    = suppressWarnings(
      wilcox.test(proportion ~ ar2, exact = FALSE)$p.value),
    .groups = "drop") |>
  mutate(p_BH  = p.adjust(p, "BH"),
         label = paste0("p[BH] == ", format(signif(p_BH, 2),
                                            scientific = FALSE)))

p_text_df <- long |>
  group_by(source) |>
  summarize(y = max(proportion, na.rm = TRUE) * 1.12, .groups = "drop") |>
  left_join(wilcox_tbl |> select(source, label), by = "source")

pB <- ggplot(long, aes(ar2, proportion, fill = source)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85, width = 0.6) +
  geom_text(data = p_text_df,
            aes(x = 1.5, y = y, label = label),
            parse = TRUE, size = 2.4, inherit.aes = FALSE) +
  scale_fill_manual(values = pal_src, guide = "none") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0.05, 0.12))) +
  facet_wrap(~ source, ncol = 5, scales = "free_y",
             labeller = labeller(source = pool_labels_facet)) +
  labs(x = NULL, y = "Source proportion") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(size = 7, face = "bold",
                                        lineheight = 0.95),
        strip.background = element_rect(fill = "gray92", color = NA),
        axis.text.x      = element_text(angle = 30, hjust = 1))

## ---- Compose --------------------------------------------------------------
fig <- (pA / pB) +
  plot_layout(heights = c(1.2, 1.0)) +
  plot_annotation(
    tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))

pdf_out <- file.path(fig_dir,
                     paste0("figure6_feast_apportionment_", out_tag, ".pdf"))
png_out <- file.path(fig_dir,
                     paste0("figure6_feast_apportionment_", out_tag, ".png"))
print(wilcox_tbl)
ggsave(pdf_out, fig, width = 183, height = 140, units = "mm", device = cairo_pdf)
ggsave(png_out, fig, width = 183, height = 140, units = "mm", dpi = 300)
cat("Saved:\n  ", pdf_out, "\n  ", png_out, "\n", sep = "")
