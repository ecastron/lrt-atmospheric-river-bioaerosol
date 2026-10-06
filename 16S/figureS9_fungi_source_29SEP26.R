# figureS9_fungi_source_29SEP26.R
# Supplementary Figure 9. Source attribution of the Risopatron airborne fungi: for each genus present in the
# atmospheric-river window, prevalence (bubble size) and mean relative abundance (color, log10) across five pools,
# separating genera associated with Patagonia that appear only during and after the event from genera prevalent in
# local soil.
# Input  : ITS/ps_LRT_ITS_R1only_decontam_withTree_29MAY26.RDS (ITS chain), 16S/palettes.R
# Outputs: 16S/figures/figureS9_fungi_source_29SEP26.{pdf,png,tsv}
# Run    : Rscript 16S/figureS9_fungi_source_29SEP26.R

suppressPackageStartupMessages({
  library(phyloseq); library(dplyr); library(tidyr); library(ggplot2)
  library(forcats); library(ggh4x); library(scales)
})
set.seed(100)
base <- "."
its_dir <- file.path(base, "ITS"); fig_dir <- file.path(base, "16S", "figures")
source(file.path(base, "16S", "palettes.R"))

ps <- readRDS(file.path(its_dir, "ps_LRT_ITS_R1only_decontam_withTree_29MAY26.RDS"))
# Drop samples with no reads left after yeast removal (RA2), so they
# do not count as absences in prevalence denominators (Before: 14, not 15)
ps <- prune_samples(sample_sums(ps) > 0, ps)
sd <- as(sample_data(ps), "data.frame"); sd$d <- as.Date(sd$date)

# --- define the five pools --------------------------------------------
pool <- rep(NA_character_, nsamples(ps)); names(pool) <- sample_names(ps)
air <- sd$sampler == "Coriolis_mu" & !sd$is_control
ter <- sd$category == "Terrestrial"
pat <- sd$region %in% c("PA_POR", "PN_SG", "Puerto Williams")
pool[air & pat]                       <- "Patagonia air"
pool[ter & pat]                       <- "Patagonia soil"
pool[ter & sd$region == "Risopatron"] <- "Risopatron local soil"
ra <- air & sd$region == "Risopatron"
pool[ra & sd$d <  as.Date("2022-02-07")] <- "Ris. air Before"
pool[ra & sd$d >= as.Date("2022-02-07")] <- "Ris. air Dur+Aft"

# --- per-genus prevalence + mean relative abundance per pool ----------
g  <- tax_glom(ps, "Genus", NArm = TRUE)          # genus-attributable only
gr <- transform_sample_counts(g, function(x){ t <- sum(x); if (t==0) x else x/t })
m  <- psmelt(gr) |> mutate(Genus = gsub("^g__", "", Genus), pool = pool[Sample]) |>
  filter(!is.na(pool))

stat <- m |> group_by(Genus, pool) |>
  summarize(prev = mean(Abundance > 0),
            mean_ab = mean(Abundance) * 100, .groups = "drop")

# --- focal set: genera present in >=2 of the 9 Dur+Aft air samples ----
focal <- stat |> filter(pool == "Ris. air Dur+Aft", prev >= 0.2) |> pull(Genus)
stat  <- stat |> filter(Genus %in% focal)

# --- classify each focal genus (transparent, vectorized rule) ---------
# transport  : present Dur+Aft, prevalent in Patagonia, rare in local soil
# local       : prevalent in local Antarctic soil
# ambiguous   : everything else
wide <- stat |> select(Genus, pool, prev) |>
  pivot_wider(names_from = pool, values_from = prev, values_fill = 0) |>
  mutate(pat_max = pmax(`Patagonia air`, `Patagonia soil`),
         local   = `Risopatron local soil`,
         durAft  = `Ris. air Dur+Aft`,
         class = factor(case_when(
           durAft >= 0.2 & pat_max >= 0.3 & local <= 0.15 ~ "Long-range transport candidate",
           local >= 0.3                                   ~ "Local Antarctic background",
           TRUE                                           ~ "Ambiguous / cosmopolitan"),
           levels = c("Long-range transport candidate",
                      "Ambiguous / cosmopolitan",
                      "Local Antarctic background")))
plotdf <- stat |> left_join(wide |> select(Genus, class), by = "Genus") |>
  mutate(pool = factor(pool, levels = c("Ris. air Before", "Ris. air Dur+Aft",
                                        "Patagonia air", "Patagonia soil",
                                        "Risopatron local soil")),
         pool_grp = factor(ifelse(grepl("^Ris. air", pool),
                                  "Antarctic air (sink)", "Sources / local"),
                           levels = c("Antarctic air (sink)", "Sources / local")))
# order genera within class by Dur+Aft prevalence
ord <- plotdf |> filter(pool == "Ris. air Dur+Aft") |>
  arrange(class, prev) |> pull(Genus)
plotdf$Genus <- factor(plotdf$Genus, levels = unique(ord))

write.table(wide, file.path(fig_dir, "figureS9_fungi_source_29SEP26.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# present cells only (no ghost dots for absent genera); label heavy cells
present <- plotdf |> filter(prev > 0)
labs    <- present |> filter(prev >= 0.30) |>
  mutate(txt_col = ifelse(mean_ab >= 0.1, "white", "gray15"))

# nicer display labels for the x axis
pool_lab <- c("Ris. air Before" = "Risopatrón air\nBefore",
              "Ris. air Dur+Aft" = "Risopatrón air\nDuring+After",
              "Patagonia air" = "Patagonia\nair",
              "Patagonia soil" = "Patagonia\nsoil",
              "Risopatron local soil" = "Risopatrón\nlocal soil")

# per-class strip styling: tint the transport band, neutral elsewhere
strip <- strip_themed(
  background_y = elem_list_rect(fill = unname(pal_s1_strip[levels(plotdf$class)])),
  text_y       = elem_list_text(face = "bold"))

# --- balloon heatmap (size = prevalence, color = mean abundance) ------
p <- ggplot(present, aes(pool, Genus)) +
  geom_point(aes(size = prev, fill = mean_ab),
             shape = 21, color = "gray30", stroke = 0.35) +
  geom_text(data = labs,
            aes(label = sprintf("%.2f", prev), color = txt_col),
            size = 1.9, fontface = "bold") +
  scale_color_identity() +
  scale_size_area(max_size = 8, name = "Prevalence", limits = c(0, 1),
                  breaks = c(0.25, 0.5, 0.75, 1)) +
  scale_fill_viridis_c(option = pal_fungal_abund_opt,
                       direction = pal_fungal_abund_dir, trans = "log10",
                       name = "Mean rel.\nabundance (%)",
                       labels = function(x) format(x, drop0trailing = TRUE, scientific = FALSE),
                       guide = guide_colorbar(barheight = 6)) +
  scale_x_discrete(labels = pool_lab) +
  scale_y_discrete(labels = function(x) sub("_gen_Incertae_sedis$", " (genus incertae sedis)", x)) +
  facet_grid2(class ~ pool_grp, scales = "free", space = "free",
              strip = strip,
              labeller = labeller(class = label_wrap_gen(16))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x   = element_text(size = 7, lineheight = 0.9),
        axis.text.y   = element_text(face = "italic", size = 7),
        strip.text.y  = element_text(angle = 0, size = 7),
        strip.text.x  = element_text(size = 7.5, face = "bold"),
        panel.grid.major = element_line(color = "gray92"),
        panel.grid.minor = element_blank(),
        legend.text   = element_text(size = 7),
        legend.title  = element_text(size = 7.5))

ggsave(file.path(fig_dir, "figureS9_fungi_source_29SEP26.pdf"), p, width = 183, height = 170, units = "mm", device = cairo_pdf)
ggsave(file.path(fig_dir, "figureS9_fungi_source_29SEP26.png"), p, width = 183, height = 170, units = "mm", dpi = 300)
cat("focal genera:", length(focal), "\n")
print(table(wide$class))
cat("Saved figureS9_fungi_source_29SEP26.{pdf,png,tsv}\n")
