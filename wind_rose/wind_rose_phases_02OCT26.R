# wind_rose_phases_02OCT26.R
# Temperature-colored wind roses at Risopatron for the Before (22 January to 6 February), During (7-8 February) and
# After (9-14 February) phases (openair::pollutionRose; 10 degree sectors, calm threshold 0.5 m/s). Records before
# 2022-01-22 19:40 (station setup) are excluded.
# Input  : wind_rose/DatosAll_CR1000X_RPT_decimalfix_02OCT26.csv (fix_decimal_separator_02OCT26.py)
# Output : 16S/figures/figureS5a_wind_roses_openair_02OCT26.{pdf,png}
# Run    : Rscript wind_rose/wind_rose_phases_02OCT26.R

library(dplyr)
library(readr)
library(lubridate)
library(openair)

base <- "."
dat <- read_csv(file.path(base, "wind_rose", "DatosAll_CR1000X_RPT_decimalfix_02OCT26.csv"), show_col_types = FALSE) |>
  transmute(date = force_tz(TIMESTAMP, "UTC"), ws = WS_ms_S_WVT, wd = WindDir_D1_WVT, temp = AirTC_Avg) |>
  filter(date >= ymd_hms("2022-01-22 19:40:00", tz = "UTC")) |>   # station setup records excluded
  filter(!is.na(ws), !is.na(wd), !is.na(temp), ws >= 0, wd >= 0, wd <= 360)

lv <- c("Before (22 Jan to 6 Feb, UTC)", "During (7 to 8 Feb, UTC)", "After (9 to 14 Feb, UTC)")
dat_periods <- dat |>
  mutate(period = case_when(
    date <  ymd_hms("2022-02-07 00:00:00", tz = "UTC") ~ lv[1],
    date <  ymd_hms("2022-02-09 00:00:00", tz = "UTC") ~ lv[2],
    date <  ymd_hms("2022-02-15 00:00:00", tz = "UTC") ~ lv[3],
    TRUE ~ NA_character_)) |>
  filter(!is.na(period)) |>
  mutate(period = factor(period, levels = lv), wd = ifelse(wd >= 360, 0, wd))
print(table(dat_periods$period))

calm_thr  <- 0.5
brks_temp <- unique(pretty(range(dat_periods$temp, na.rm = TRUE), n = 6))

draw <- function() {
  pollutionRose(
    mydata     = dat_periods,
    pollutant  = "temp",
    ws         = "ws",
    wd         = "wd",
    type       = "period",
    angle      = 10,
    breaks     = brks_temp,
    paddle     = TRUE,
    calm       = calm_thr,
    key.title  = "Air temperature (°C)",   # openair 3.x uses key.title in place of key.header/key.footer
    main       = ""
  )
}

out <- file.path(base, "16S", "figures", "figureS5a_wind_roses_openair_02OCT26")
pdf(paste0(out, ".pdf"), width = 12, height = 4.6); draw(); dev.off()
png(paste0(out, ".png"), width = 12, height = 4.6, units = "in", res = 300); draw(); dev.off()
cat("wrote", out, "\n")
