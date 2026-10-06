# figureS7_pdust_risopatron.R
# Supplementary Figure 7: CHIMERE surface mineral dust (pDUST) at Risopatron, 1 January to 25 February 2022, in the
# three nested model domains (one panel each, own y scale). Input: atmospheric_simulations/CHIMERE/
# pdust_risopatron_from_plot.csv (from extract_pdust_risopatron.py). Output: 16S/figures/figureS7_pdust_risopatron.{pdf,png}
# Run   : Rscript atmospheric_simulations/CHIMERE/figureS7_pdust_risopatron.R  (from the repository root)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(tibble)})
LRT <- getwd()  # repository root
source(file.path(LRT, "16S", "palettes.R"))

d <- read_csv(file.path(LRT, "atmospheric_simulations", "CHIMERE", "pdust_risopatron_from_plot.csv"),
              show_col_types = FALSE) |>
  mutate(time_utc = as.POSIXct(time_utc, tz = "UTC"),
         domain = factor(domain, levels = c("d01", "d02", "d03")))
stopifnot(setequal(levels(d$domain), names(pal_chimere_domain)), !anyNA(d$domain))

# AR core shaded as in the other time-series panels (ERA5 IVT window, 7-8 February UTC)
ar <- tibble(xmin = as.POSIXct(ar_event_start, tz = "UTC"), xmax = as.POSIXct(ar_event_end + 1, tz = "UTC"))

p <- ggplot(d, aes(time_utc, pdust_ug_m3, color = domain)) +
  geom_rect(data = ar, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
            fill = pal_ar_phase[["During"]], alpha = 0.12) +
  geom_line(linewidth = 0.45) +
  # one panel per domain with its own y scale: the single-hour AR peaks overlap on a shared axis, and the
  # outer domain's peak (0.45) would be invisible next to the inner domains (6.5-8.1)
  facet_wrap(~domain, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(d01 = "d01 (outer domain)", d02 = "d02", d03 = "d03 (inner domain)"))) +
  scale_color_manual(values = pal_chimere_domain, guide = "none") +
  scale_x_datetime(date_breaks = "1 week", date_labels = "%d %b", expand = expansion(mult = 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  labs(x = NULL, y = expression("Mineral dust (pDUST, " * mu * "g m"^-3 * ")")) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(strip.text = element_text(size = 7, hjust = 0), panel.grid.minor = element_blank())

out <- file.path(LRT, "16S", "figures", "figureS7_pdust_risopatron")
ggsave(paste0(out, ".pdf"), p, width = 120, height = 100, units = "mm", device = cairo_pdf)
ggsave(paste0(out, ".png"), p, width = 120, height = 100, units = "mm", dpi = 300)
