#!/usr/bin/env Rscript
# figureS_shotgun_faprotax_05OCT26.R
# Supplementary Figure 2. Percentage of QC-passing air metagenomes in which each function is inferred by FAPROTAX and
# in which its marker genes are detected; functions that FAPROTAX does not model are shown separately. CO oxidation is
# not shown because no air CoxL gene met the form I criteria (coxL_reanalysis_01OCT26).
# Inputs : shotgun_metagenomics/tables/air_faprotax_corroboration_29SEP26.csv,
#          shotgun_metagenomics/tables/air_shotgun_only_functions_29SEP26.csv
# Output : 16S/figures/figureS_shotgun_faprotax_05OCT26.{pdf,png}
# Run    : Rscript shotgun_metagenomics/scripts/figureS_shotgun_faprotax_05OCT26.R

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})

ROOT <- "."
SG   <- file.path(ROOT, "shotgun_metagenomics")
OUT  <- file.path(ROOT, "16S", "figures", "figureS_shotgun_faprotax_05OCT26")

modeled <- read.csv(file.path(SG, "tables", "air_faprotax_corroboration_29SEP26.csv")) |>
  filter(!is.na(shotgun_detected_n)) |>
  transmute(group, panel = "FAPROTAX modeled",
            FAPROTAX_inferred = fap_inferred_n, shotgun_detected = shotgun_detected_n,
            n_air)
beyond <- read.csv(file.path(SG, "tables", "air_shotgun_only_functions_29SEP26.csv")) |>
  filter(group != "CO oxidation (form-I CoxL)") |>   # CoxL forms are reported separately (coxL_reanalysis_01OCT26)
  transmute(group, panel = "Shotgun only",
            FAPROTAX_inferred = NA_integer_, shotgun_detected = shotgun_detected_n,
            n_air)

relabel <- c("Methanol ox / methylotrophy"      = "Methanol oxidation / methylotrophy",
             "Denitrification / NO3 red"        = "Denitrification / nitrate reduction",
             "Nitrification / NH3 ox"           = "Nitrification / ammonia oxidation",
             "Dark H2 oxidation (HydDB grp 1l)" = "Dark H2 oxidation (HydDB group 1l)")
pl <- bind_rows(modeled, beyond) |>
  mutate(group = dplyr::recode(group, !!!relabel)) |>
  pivot_longer(c(FAPROTAX_inferred, shotgun_detected),
               names_to = "layer", values_to = "n_samples") |>
  filter(!is.na(n_samples)) |>
  mutate(frac = n_samples / n_air)

stopifnot(nrow(pl) > 0, all(pl$frac <= 1))

# order rows by their highest value so the two panels read as one ranking
ord <- pl |> group_by(group) |> summarize(m = max(frac), .groups = "drop") |>
  arrange(m) |> pull(group)
pl <- pl |> mutate(
  group = factor(group, levels = ord),
  panel = factor(panel, levels = c("FAPROTAX modeled",
                                   "Shotgun only")))

p <- ggplot(pl, aes(frac, group, color = layer)) +
  # shotgun-only rows carry a single point, so their connector draws nothing
  geom_line(aes(group = group), color = "gray78", linewidth = 0.6) +
  geom_point(size = 2.8) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_color_manual(values = c(FAPROTAX_inferred = "#999999",
                                shotgun_detected  = "#0072B2"),
                     labels = c("FAPROTAX inferred (16S rRNA)",
                                "shotgun gene detected")) +
  labs(x = "QC-passing air libraries with the function (% of 11)", y = NULL, color = NULL) +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(legend.position = "top",
        legend.margin = margin(b = -4),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.y = element_text(face = "bold", size = 8, angle = -90),
        axis.text.y = element_text(size = 8.5),
        plot.margin = margin(4, 6, 4, 4))

ggsave(paste0(OUT, ".pdf"), p, width = 7.1, height = 4.0)
ggsave(paste0(OUT, ".png"), p, width = 7.1, height = 4.0, dpi = 300)
cat("wrote", paste0(OUT, ".pdf / .png"), "\n")
