library(data.table)
library(ggplot2)
library(readxl)
library(emmeans)
library(patchwork)
# stats packages 
library(car)
library(DEoptim) # for getting the best starting values for fitting
# font etc and customised functions
source("scripts/utilities.R")

# ==============================================================================
# TREATMENT DATES
# ==============================================================================
# Planting and key phenological dates for each season
planting_dates <- list(`2024` = as.Date("2024-05-13"), 
                       `2025` = as.Date("2025-05-07"))
heating_starts <- list(`2024_HE` = as.Date("2024-06-07"), 
                       `2025_HTI` = as.Date("2025-06-13"),
                       `2024_HTI` = as.Date("2024-06-28"), 
                       `2025_HB` = as.Date("2025-07-12"))
heating_ends <- list(`2024_HE` = as.Date("2024-06-21"), 
                     `2025_HTI` = as.Date("2025-06-29"),
                     `2024_HTI` = as.Date("2024-07-15"), 
                     `2025_HB` = as.Date("2025-07-28"))
leaf_removal_dates <- list(`2024` = as.Date("2024-09-02"),
                           `2025` = as.Date("2025-08-18"))
Season_start_end <- list(`2024_start` = as.Date("2024-05-13"),
                         `2024_end` = as.Date("2024-09-09"),
                         `2025_start` = as.Date("2025-05-07"),
                         `2025_end` = as.Date("2025-08-28"))
emergence_date <- list(`2024` = as.Date("2024-06-07"),
                       `2025` = as.Date("2025-05-23"))
# Tuber initiation dates used for supplementary reference
# 2025 has a single shared date; 2024 differs by treatment due to early heat stress
tuber_initation_dates <- list(`2024_AC` = as.Date("2024-06-24"),
                              `2024_HE` = as.Date("2024-06-20"),
                              `2025` = as.Date("2025-06-10"))

# ==============================================================================
# HARVEST SCHEDULE
# ==============================================================================
Harvest_dates <- data.table(
  Harvest_dates = as.Date(c("2024-06-24", "2024-07-05", "2024-07-15", "2024-07-29",
                    "2024-08-19",
                    "2025-07-02", "2025-07-14", "2025-07-29", "2025-08-11",
                    "2024-09-10", "2025-08-26")))
Harvest_dates[, Season := year(Harvest_dates)]
Harvest_dates[Season == 2024, DAP := as.integer(Harvest_dates - planting_dates$`2024`)]
Harvest_dates[Season == 2025, DAP := as.integer(Harvest_dates - planting_dates$`2025`)]
Harvest_dates[, HarvestNr := 1:.N, by = .(Season)]

# ==============================================================================
# VISUAL STYLING
# ==============================================================================
# Treatment colour palette: AC = black (control), heat treatments = distinct hues
colors_temp <- c(
  "AC"  = "#000000",
  "HE"  = "#314856",
  "HTI" = "#8B1C1C",
  "HB"  = "#fd8700"
)
linetype_temp <- c(
  "AC"  = "solid",    # control: solid, most prominent
  "HE"  = "dashed",
  "HTI" = "solid",    # focal: solid, matches control prominence
  "HB"  = "dotdash"
)
point_shape <- c("AC" = 16, "HE" = 1, "HTI" = 17, "HB" = 2)  # AC filled circle, HTI filled triangle, others open
legend_order <-  c("AC", "HE", "HTI", "HB")
# Base plot parameters
lw = 0.8        # default line width
ps = 2          # default point size
xlims <- c(0, 120)
fontsize     <- 12
fontsize_ppt <- 16

# ==============================================================================
# DERIVED DATE STRUCTURES (for plot annotations)
# ==============================================================================
# Long-form date tables used to draw treatment shading rectangles on figures
heating_starts_dt <- as.data.table(heating_starts) |>
  melt.data.table(variable.factor = FALSE)
heating_ends_dt <- as.data.table(heating_ends) |>
  melt.data.table(variable.factor = FALSE)
planting_dt <- as.data.table(planting_dates) |>
  melt.data.table(variable.factor = FALSE,
                  variable.name = "Season",
                  value.name = "Planting")
planting_dt[, Season := as.integer(Season)]

# Combine starts and ends into a single long table; compute DAP from planting
treatments_DAP <- rbindlist(list(heating_starts_dt, heating_ends_dt),
                            use.names = TRUE, fill = TRUE)
treatments_DAP[, Season := year(value)]
treatments_DAP <- treatments_DAP[planting_dt, on = "Season"]
treatments_DAP[, DAP := as.integer(value - Planting)]
treatments_DAP[, Treatment := gsub("(.+)(H.+)", "\\2", variable)]

# DAP positions for start/end arrows annotating treatment periods on figures
arrow_map <- data.table::data.table(
  Season = rep(c(2024, 2025), each = 4),
  DAP    = c(25, 46, 39, 63, 37, 66, 53, 82),
  arrow_color = c("red", "red", "blue", "blue", "red", "red", "blue", "blue"))

treatments_DAP_arrows <- merge(treatments_DAP, arrow_map, by = c("Season", "DAP"))
som_sampling <- data.table(Harvest = 1:4,
                           Date = as.Date(c("2025-07-02", "2025-07-14",
                                            "2025-07-29", "2025-08-11")))

# Treatment shading rectangles: one row per Season x Treatment window
treatment_rect_dt <- copy(treatments_DAP_arrows)
treatment_rect_dt <- treatment_rect_dt[
  , .(DAP_min = min(DAP), DAP_max = max(DAP)),
  by = .(Season, Treatment)
]

# ==============================================================================
# CURVE FITTING CONSTANTS
# ==============================================================================
# Plot 8 in 2025 was excluded due to sensor malfunction
EXCLUDE_PLOT_2025 <- 8
# DAP at which canopy cover dropped below 0.05 (used to truncate predictions)
CANOPY_DEATH_DAP_2025 <- 98
CANOPY_DEATH_DAP_2024 <- 108 # From canopy cover curve: drop below 0.05


