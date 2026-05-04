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
# Experimental constants (same across seasons)
get_experiment_constants <- function() {
  list(
    total_plants = 12,                    # total per plot
    row_spacing = 0.75,                   # meters
    plant_spacing = 0.3,                  # meters
    plant_area = 0.75 * 0.3,             # m² per plant
    plot_area = 12 * 0.75 * 0.3          # m² per plot
  )
}

# Season-specific constants
get_season_constants <- function(season) {
  constants <- get_experiment_constants()
  
  if (season == "2024") {
    constants$plants_detailed <- 3
    constants$plants_remaining <- 9
    constants$planting_date <- as.Date("2024-05-13")
    constants$final_harvest_DAP <- 120
    constants$HE_start <- as.Date("2024-06-07")
    constants$HE_end <- as.Date("2024-06-21")
    constants$HTI_start <- as.Date("2024-06-28")
    constants$HTI_end <- as.Date("2024-07-15")
  } else if (season == "2025") {
    constants$plants_harvested_intermediate <- 6
    constants$plants_harvested_final <- 12
    constants$planting_date <- as.Date("2025-05-07")
    constants$final_harvest_DAP <- 111
    constants$intermediate_scale_factor <- 12 / 6  # 2
    constants$final_scale_factor <- 12 / 12        # 1
    constants$HTI_start <- as.Date("2025-06-13")
    constants$HTI_end <- as.Date("2025-06-29")
    constants$HB_start <- as.Date("2025-07-12")
    constants$HB_end <- as.Date("2025-07-28")
  }
  
  return(constants)
}
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
tuber_initation_dates <- list(`2024_AC` = as.Date("2024-06-24"),
                              `2024_HE` = as.Date("2024-06-20"),
                              `2025` = as.Date("2025-06-10"))
Harvest_dates <- data.table(
  Harvest_dates = as.Date(c("2024-06-24", "2024-07-05", "2024-07-15", "2024-07-29",
                    "2024-08-19",
                    "2025-07-02", "2025-07-14", "2025-07-29", "2025-08-11",
                    "2024-09-10", "2025-08-26")))
Harvest_dates[, Season := year(Harvest_dates)]
Harvest_dates[Season == 2024, DAP := as.integer(Harvest_dates - planting_dates$`2024`)]
Harvest_dates[Season == 2025, DAP := as.integer(Harvest_dates - planting_dates$`2025`)]
Harvest_dates[, HarvestNr := 1:.N, by = .(Season)]
# theme scheme
colors_temp <- c(
  "AC"  = "#000000",  # black (control reference)
  "HE"  = "#999999",  # mid grey (muted)
  "HTI" = "#D55E00",  # vermillion (focal heat treatment)
  "HB"  = "#BDBDBD"   # light grey (muted)
)
linetype_temp <- c(
  "AC"  = "solid",    # control: solid, most prominent
  "HE"  = "dashed",
  "HTI" = "solid",    # focal: solid, matches control prominence
  "HB"  = "dotdash"
)
point_shape <- c("AC" = 16, "HE" = 1, "HTI" = 17, "HB" = 2)  # AC filled circle, HTI filled triangle, others open
legend_order <-  c("AC", "HE", "HTI", "HB")
colors_temp_old <- c(
  "AC" = "#2166AC",   # Cool blue (Ambient Control)
  "HE" = "#F46D43",   # Warm orange (Heat Early - 2024)
  "HTI" = "#FDAE61",  # Warm yellow (Heat Tuber Initiation)
  "HB" = "#A50026"    # Hot red (Heat Bulk - 2025)
)
# Define the colors as a named character vector
colors_temp <- c(
  "AC"          = "#000000",
  "HE" = "#314856",
  "HTI"  = "#8B1C1C",
  "HB" = "#fd8700"
)
lw = 0.8
ps = 2
xlims <- c(0,120)
fontsize <-  12
fontsize_ppt <- 16
heating_starts_dt <- as.data.table(heating_starts) |> 
  melt.data.table(variable.factor = FALSE)
heating_ends_dt <- as.data.table(heating_ends) |> 
  melt.data.table(variable.factor = FALSE)
planting_dt <- as.data.table(planting_dates) |> 
  melt.data.table(variable.factor = FALSE,
                  variable.name = "Season", 
                  value.name = "Planting")
planting_dt[, Season := as.integer(Season)]
treatments_DAP <- rbindlist(list(heating_starts_dt, heating_ends_dt), 
                            use.names = TRUE, fill = TRUE)
treatments_DAP[, Season := year(value)]
treatments_DAP <- treatments_DAP[planting_dt, on = "Season"]
treatments_DAP[, DAP := as.integer(value - Planting)]
treatments_DAP[, Treatment := gsub("(.+)(H.+)", "\\2", variable)]

# Create arrow map for treatment annotations
arrow_map <- data.table::data.table(
  Season = rep(c(2024, 2025), each = 4),
  DAP    = c(25,46,39,63,37,66,53,82),
  arrow_color = c("red","red","blue","blue","red","red","blue","blue"))

treatments_DAP_arrows <- merge(treatments_DAP, arrow_map, by = c("Season","DAP"))
som_sampling <- data.table(Harvest = 1:4, 
                           Date = as.Date(c("2025-07-02", "2025-07-14",
                                            "2025-07-29", "2025-08-11")))
# Treatment windows in DAP, recoded to the same treatment names used in the canopy plot
treatment_rect_dt <- copy(treatments_DAP_arrows)

# One start/end DAP window per season and treatment
treatment_rect_dt <- treatment_rect_dt[
  , .(DAP_min = min(DAP), DAP_max = max(DAP)),
  by = .(Season, Treatment)
]

# constants for curve fitting  --------------------------------------------

# Define constants -------------------------------------------------------------
EXCLUDE_PLOT_2025 <- 8  # Plot to exclude in 2025 season
CANOPY_DEATH_DAP_2025 <- 98  # Day after planting when canopy died in 2025
CANOPY_DEATH_DAP_2024 <- 108# VALUE from Canopy cover curve drop below 0.05

# key time points  --------------------------------------------------------

# treatments_DAP
# leaf_removal_dates
# Season_start_end
# tar_load(daily_td)
# daily_td[Date == leaf_removal_dates$`2024`]
# daily_td[Date == leaf_removal_dates$`2025`]






