# harmonize_metadata.R
# Builds one sample metadata table for the LRT and ACE (Malard et al. 2022) samples with harmonized sample classes.
# Inputs : 16S/metadata_LRT3.csv (semicolon-separated), 16S/ACEDATA/metadata_final_all.csv (semicolon-separated)
# Output : 16S/metadata_merged_<DATE>.csv (DATE is today's date; merge_LRT_ACE.R reads metadata_merged_07MAY26.csv)
# Run    : Rscript 16S/harmonize_metadata.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
})

base     <- "./16S"
date_tag <- "07MAY26"   # tag of the outputs used in the manuscript; downstream scripts read *_07MAY26 files

lrt <- read_delim(file.path(base, "metadata_LRT3.csv"),
                  delim = ";", show_col_types = FALSE)
ace <- read_delim(file.path(base, "ACEDATA", "metadata_final_all.csv"),
                  delim = ";", show_col_types = FALSE)

# --- LRT --------------------------------------------------------------
# Source2 already matches the figure-legend categories (Antarctic Air, Antarctic Soil, Patagonian Air,
# Patagonian Soil, Control). Keep as-is.
lrt_h <- lrt |>
  transmute(
    sample_id  = sample_name,
    dataset    = "LRT",
    source2    = Source2,
    category   = Category2,
    latitude   = Latitude,
    longitude  = Longitude,
    date       = dmy(collection_date),
    region     = region,
    sampler    = if_else(Source == "air", "Coriolis_mu", Source),
    is_control = Source %in% c("buffer_control", "coriolis_control",
                               "deposition_control", "sea_water_control"),
    atm_river  = atm_river
  )

# --- ACE / Malard 2022 -------------------------------------------------
# Environment -> source2 mapping (figure-legend categories):
#   Marine               -> SO Air                  (over open ocean)
#   Marine Precipitation -> SO Precipitation
#   Terrestrial          -> SO Air over land        (sub-Antarctic islands)
#   NEG / NEG1 / POS     -> Control
ace_h <- ace |>
  transmute(
    sample_id  = as.character(Sample_ID),
    dataset    = "Malard2022",
    source2    = case_when(
      Environment == "Marine"               ~ "SO Air",
      Environment == "Marine Precipitation" ~ "SO Precipitation",
      Environment == "Terrestrial"          ~ "SO Air over land",
      Environment %in% c("NEG", "NEG1", "POS") ~ "Control",
      TRUE                                  ~ NA_character_
    ),
    category   = case_when(
      Environment == "Marine"               ~ "Air",
      Environment == "Marine Precipitation" ~ "Marine Precipitation",
      Environment == "Terrestrial"          ~ "Air",
      Environment %in% c("NEG", "NEG1", "POS") ~ "Control",
      TRUE                                  ~ NA_character_
    ),
    latitude   = suppressWarnings(as.numeric(Latitude)),
    longitude  = suppressWarnings(as.numeric(Longitude)),
    date       = dmy(Date),
    region     = NA_character_,
    sampler    = Sampler,
    is_control = Environment %in% c("NEG", "NEG1", "POS") |
                 Sampler %in% c("NEG", "POS"),
    atm_river  = NA_character_
  )

# Drop ACE samples whose Environment was blank in the source metadata
# (3 SKC samples: 351, 347a, 350a): unclassified at source.
n_drop <- sum(is.na(ace_h$source2))
if (n_drop > 0) {
  message("Dropping ", n_drop,
          " ACE rows with NA source2 (blank Environment in source metadata)")
  ace_h <- filter(ace_h, !is.na(source2))
}

merged <- bind_rows(lrt_h, ace_h)

# Sample IDs: LRT uses character codes (GUS1, RA1, ...), ACE uses
# integers. They should not collide: assert.
dups <- merged$sample_id[duplicated(merged$sample_id)]
if (length(dups) > 0) {
  warning("Duplicated sample_ids across datasets: ",
          paste(unique(dups), collapse = ", "),
          ": prefixing ACE ids with 'ACE_'")
  merged <- merged |>
    mutate(sample_id = if_else(dataset == "Malard2022",
                               paste0("ACE_", sample_id),
                               sample_id))
}

# Any unmapped Source2? Worth knowing.
if (any(is.na(merged$source2))) {
  warning(sum(is.na(merged$source2)),
          " rows have NA source2: check Environment values:\n",
          paste(unique(ace$Environment[is.na(ace_h$source2)]),
                collapse = ", "))
}

out_path <- file.path(base, paste0("metadata_merged_", date_tag, ".csv"))
write_csv(merged, out_path)

message("Wrote ", out_path)
message(sprintf("  %d LRT + %d Malard2022 = %d total samples",
                sum(merged$dataset == "LRT"),
                sum(merged$dataset == "Malard2022"),
                nrow(merged)))
print(merged |> count(dataset, source2))
