source("scripts/00_setup.R")

daily_temp_data <- fread("data/daily_weather.csv", skip = 1)
# rain  -------------------------------------------------------------------

p1 <- ggplot(daily_temp_data, aes(x = DAP)) +
  geom_col(aes(y = Rain), alpha = 0.7, fill = "black") +
  geom_line(aes(y = Tmean, color = "Tmean"), linewidth = lw - 0.3, show.legend = FALSE) +
  geom_line(aes(y = Tmax, color = "Tmax"), linewidth = lw - 0.3, linetype = "dashed") +
  geom_line(aes(y = Tmin, color = "Tmin"), linewidth = lw - 0.3, linetype = "dotted") +
  scale_y_continuous(
    name = "Temperature (°C)",expand = c(0,0),
    sec.axis = sec_axis(~ ., name = "Rainfall (mm)"), 
    limits = c(0, 40)
  ) +
  scale_x_continuous(expand = c(0,0),limits = c(1,120)) +
  geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colors_temp, guide = "none") +
  scale_color_manual(values = c("Tmean" = "black", "Tmax" = "black", "Tmin" = "black")) +
  labs(color = "Legend") +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  facet_grid(~Season) +
  theme(# strip.background = element_blank(),
        legend.box.margin = margin(b = -5, unit = "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(5, "mm"))

daily_temp_data[, cum_GDD := cumsum(Tmean - 2), by = Season]
# calcualte summarise stats for paper 
# Weather summary statistics for methods section
weather_summary <- daily_temp_data[, .(
  Mean_Temp = round(mean(Tmean, na.rm = TRUE), 1),
  Min_Temp = round(min(Tmin, na.rm = TRUE), 1),
  Max_Temp = round(max(Tmax, na.rm = TRUE), 1),
  cum_GDD  = max(cum_GDD),
  Total_Rain = round(sum(Rain, na.rm = TRUE), 1),
  Rain_Days = sum(Rain > 0, na.rm = TRUE),
  cum_Rad = round(sum(IRRAD, na.rm = TRUE), 1),
  Growing_Days = .N
), by = Season]
print(weather_summary)

daily_temp_data[, Cumulative_Radiation := cumsum(IRRAD), by = Season]

# Radiation plot
p2 <- ggplot(daily_temp_data, aes(x = DAP, y = IRRAD)) +
  geom_line(color = "grey50", alpha = 0.8, linewidth = lw) +
  geom_line(aes(y = Cumulative_Radiation/80), color = "black", linewidth = lw - 0.3) +
  scale_y_continuous(
    name = expression(atop(Daily~Radiation, (MJ~m^{-2}))),
    sec.axis = sec_axis(~ . * 80, name = expression(atop(Cumulative~Radiation, (MJ~m^{-2})))), 
    limits = c(0, 35), expand = c(0,0)
  ) + 
  scale_x_continuous(name = "Days After Planting", expand = c(0,0), limits = c(0,120)) +
  geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE) +
  scale_fill_manual(values = colors_temp, guide = "none") +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(panel.spacing.x = unit(5, "mm"),
        panel.grid = element_blank())+
  facet_grid(~Season)

weather_plot <- (p1 +
                   theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank()))/ (
                           p2 + 
                             theme(strip.text = element_blank(),
                                   strip.background = element_blank())) + 
  plot_layout(guides = "collect", axes = "collect", heights = c(0.95, 1)) & 
  theme(legend.position = "top")
save_plot(weather_plot, "fig1_weather")

# treatment period weather metrics ----------------------------------------
# Get DAP ranges for each treatment period

trt_ranges <- treatments_DAP_arrows[
  variable %in% c("2024_HE", "2024_HTI", "2025_HTI", "2025_HB"),
  .(DAP_start = DAP[1], DAP_end = DAP[2]), by = variable
][, Season := as.integer(substr(variable, 1, 4))]

# Non-equi join: match daily_temp_data rows falling within each period
air_periods <- daily_temp_data[trt_ranges, on = .(Season, DAP >= DAP_start, DAP <= DAP_end), nomatch = NULL]
# Mean temperature per treatment period
mean_temps <- air_periods[, .(Mean_Temp = round(mean(Tmean, na.rm = TRUE), 1),
                              sd_Temp = sd(Tmean), 
                              Rain = sum(Rain) |> round(),
                              Rain_mean = sum(Rain)/.N,
                              Days = .N), by = .(Treatment = variable)]
print(mean_temps)
rad_periods <- daily_temp_data[trt_ranges, on = .(Season, DAP >= DAP_start, DAP <= DAP_end), nomatch = NULL]
sum_rad_periods <- rad_periods[, .(rad_mean = mean(IRRAD) |> round(),
                                   rad_sum= sum(IRRAD) |> round()), by = .(Treatment = variable)]
weather_metric_periods <- mean_temps[sum_rad_periods, on = "Treatment"]

# soil and treatment verification -----------------------------------------
# Summary stats #!!! Plot 8 
soil_temperature_daily <- fread("data/soil_temperature_daily.csv")
soil_temperature_daily <- soil_temperature_daily[!(Season == 2025 & 
                                                     Treatment == "HB" & 
                                                     PlotID == 8)]
# calculate mean and sd time series by Treatment/Depth/Season
soil_temp_ts <- soil_temperature_daily[, .(
  mean_temp = mean(mean_temp, na.rm = TRUE),
  sd_temp = sd(mean_temp, na.rm = TRUE)
), by = .(Season, Depth, Treatment, Date, DAP)
][order(Season, Depth, Treatment, Date, DAP)]

# 2-day moving average for mean and upper/lower bounds
soil_temp_ts[, `:=`(
  mean_ma = frollmean(mean_temp, n = 2, align = "right"),
  upper_ma = frollmean(mean_temp + sd_temp, n = 2, align = "right"),
  lower_ma = frollmean(mean_temp - sd_temp, n = 2, align = "right")
), by = .(Season, Depth, Treatment)]
soil_temp_ts[, Depth := paste0(Depth, " cm")]
fsize <- 12
# Plot mean +/- SD  -----------------
panel_tag_soil <- data.table(Season = c(2024, 2025), Depth = c("20 cm", "20 cm"),
                             label = c("a", "b"))
soil_temp_sd_p <- soil_temp_ts[Depth != "40 cm"] |>
  ggplot(aes(DAP, mean_temp)) +
  geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = c(14, 22), color = "#e74c3c")+
  geom_ribbon(aes(ymin = mean_temp - sd_temp, 
                  ymax = mean_temp + sd_temp,
                  fill = Treatment),
              alpha = 0.5) +
  geom_line(aes(linetype = Treatment, color = Treatment), linewidth = lw - .3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 120)) +
  scale_y_continuous(limits = c(10, 35))+
  facet_grid( ~ Season) +
  scale_linetype_manual(
    values = linetype_temp,
    breaks = legend_order
  ) +
  scale_color_manual(
    values = colors_temp,
    breaks = legend_order
  ) +
  scale_fill_manual(
    values = colors_temp,
    breaks = legend_order
  ) +
  geom_text(data = panel_tag_soil,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.5,    vjust = 1.3,    size = 5,  inherit.aes = FALSE) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(legend.position = "top",
        panel.spacing.x = unit(5, "mm"),
        # axis.text.x = element_blank(),
        # axis.title.x =  element_blank(),
        panel.grid = element_blank(),
        legend.key.width = unit(15, "mm"),
      legend.title = element_text(margin = margin(r = 15))) +
  labs(y = "Soil temperature (°C)", x = "Days after planting")
soil_temp_sd_p
save_plot(soil_temp_sd_p, "fig2_soil_daily", width = 6, height = 3)

air_temp_p <- ggplot(data = daily_temp_data)+
  geom_line(aes(x= DAP, y = Tmean), inherit.aes = FALSE, 
          linewidth = lw - .3) +
  facet_grid( ~ Season) +
  geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colors_temp, guide = "none") +
  theme_bw(base_size = fontsize_ppt, base_family = "Times New Roman") +
  theme(legend.position = "top",
        panel.spacing.x = unit(5, "mm"),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.key.width = unit(15, "mm")) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 120)) +
  scale_y_continuous(limits = c(9, 35))+
  labs(y = "Air Temperature (°C)", x = "Days After Planting")

soil_air_p <- soil_temp_sd_p /air_temp_p + plot_layout(heights = c(2,1.5))
save_plot(soil_air_p, "fig2_soil_air_daily")
## the soil temperature is already averaged across blocks!!!!

## filtering  daily value ------------------------------------------------------
HE_24 <- soil_temperature_daily[
      Season == 2024 &
        Date >= heating_starts$`2024_HE` &
        Date <= heating_ends$`2024_HE` &
        Treatment %in% c("AC", "HE")
    ]
HTI_24 <- soil_temperature_daily[
      Season == 2024 &
        Date >= heating_starts$`2024_HTI` &
        Date <= heating_ends$`2024_HTI` &
        Treatment %in% c("AC", "HTI")
    ]
HTI_25 <- soil_temperature_daily[
      Season == 2025 &
        Date >= heating_starts$`2025_HTI` &
        Date <= heating_ends$`2025_HTI` &
        Treatment %in% c("AC", "HTI")
    ]
HB_25 <- soil_temperature_daily[
      Season == 2025 &
        Date >= heating_starts$`2025_HB` &
        Date <= heating_ends$`2025_HB` &
        Treatment %in% c("AC", "HB")
    ]

## calculating the mean, max and min, and difference between treatments for each period
calculate_temp_stats <- function(data, 
                                 grouping_vars = c("Season","Treatment", "Depth")) {
  data[, .(
    mean_temp = round(mean(mean_temp, na.rm = TRUE)),
# the max, min etc doesn't make sense anymore since it would calculate temporal max and mins
    max_temp = round(max(mean_temp, na.rm = TRUE)),
    min_temp = round(min(mean_temp, na.rm = TRUE)),
    n = .N,
    sd = round(sd(mean_temp, na.rm = TRUE), 1)
  ), by = grouping_vars]
}

stats_h1_24 <- calculate_temp_stats(HE_24)
stats_h2_24 <- calculate_temp_stats(HTI_24)
stats_h1_25 <- calculate_temp_stats(HTI_25)
stats_h2_25 <- calculate_temp_stats(HB_25)

# Format for manuscript: mean ± sd (all results)
all_stats <- rbindlist(list(
  stats_h1_24[, Period := "2024_HE"],
  stats_h2_24[, Period := "2024_HTI"],
  stats_h1_25[, Period := "2025_HTI"],
  stats_h2_25[, Period := "2025_HB"]
))
all_stats[, formatted := paste0(mean_temp, "±", sd)]
fig1_stats <- all_stats[, .(Season, Period, Treatment, Depth, formatted)] |>
  dcast(Season + Period + Treatment ~ Depth, value.var = "formatted") |>
  setnames(as.character(c(20)), c("formatted_20")) 
col_control <- fig1_stats[Treatment == "AC",.(Season, Period, Treatment, formatted_20)]
col_trts_soilT <- fig1_stats[Treatment != "AC",.(Season, Period, Treatment, formatted_20)]
TableS1 <- col_trts_soilT[weather_metric_periods, on = "Period == Treatment"
               ][col_control, on = "Period"]
irrigation <- fread("data/irrigation_data.csv")

HE_2024 <- irrigation[Season == 2024 & DAP >= 25 & DAP <= 39
                      ][,.(Irrigation = round(sum(mm))), 
                        by = .(Treatment)]
HTI_2024 <- irrigation[Season == 2024 & DAP >= 46 & DAP <= 63
                       ][,.(Irrigation = round(sum(mm))), 
                         by = .(Treatment)]
HTI_2025 <- irrigation[Season == 2025 & DAP >= 37 & DAP <= 53
                       ][,.(Irrigation = round(sum(mm))), 
                         by = .(Treatment)]
HB_2025 <- irrigation[Season == 2025 & DAP >= 66 & DAP <= 82
                      ][,.(Irrigation = round(sum(mm))), 
                        by = .(Treatment)]
# value selected from the calculation above 
TableS1[, Irrigation := c(5,4, 25, 90)]
TableS1[, Total_water := Irrigation + Rain]
TableS1[, Treatment := factor(Treatment, levels = legend_order)]
setorder(TableS1, Season, Treatment)

# rename the col and change order 
TableS1[, .(Season, Treatment, Days, `Soil temperature at 20 cm` = formatted_20, 
             `Ambient soil temperature` = i.formatted_20,
             `Mean air temperature` = sprintf("%.0f%s%.1f", Mean_Temp, "±", sd_Temp),
             Rain, Irrigation, Irradiance = rad_sum)] 

# canopy temperature as a verification of non stressed canopy -------------

canopy_temp <- fread("data/canopy_temperature.csv")
# Add date as character for labeling
canopy_temp[, Date_label := format(Date, "%b %d")]
canopy_temp_sum <- canopy_temp[, .(mean = mean(canopy_temp),
                                   sd = sd(canopy_temp)), 
                               by = .(Season, Date, DAP, Treatment)]
# visualise the canopy temperature time series 
canopy_temp_sum |>
  ggplot(aes(DAP, mean, color = Treatment)) +
  geom_point(size = ps, alpha = 0.6, position = position_dodge(width = 5)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 3,
                position = position_dodge(width = 5)) +
  facet_grid(~Season) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(legend.position = "top",
        strip.background = element_blank()) +
  labs(x = "Days After Planting",
       y = "Canopy Temperature (°C)",
       color = "Treatment")
# Timeseries visualization of canopy temperature with air temperature
canopy_temp_plot <- canopy_temp_sum |>
  ggplot(aes(DAP, mean, shape = Treatment)) +
  geom_point(size = ps, alpha = 0.6, position = position_dodge(width = 6)) +
  geom_rect(
    data = treatment_rect_dt,
    aes(xmin =  DAP_min, xmax = DAP_max, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 3, 
                position = position_dodge(width = 6), alpha = 0.6) +
  # geom_line(aes(DAP, air_temp_mean, linetype = "Air Temperature"), 
            # color = "black",
            # linewidth = lw, alpha = 0.7) +
  facet_wrap(~Season, ncol = 1) +
  scale_y_continuous(
    name = "Surface canopy temperature (°C)",
    limits = c(0, 35),
    expand = c(0, 0)
  ) +
  scale_shape_manual(values = point_shape, breaks = c("AC", "HE", "HTI", "HB")) +
  # scale_color_manual(values = colors_temp, breaks = c("AC", "HE", "HTI", "HB")) +
  scale_fill_manual(values = colors_temp, breaks = c("AC", "HE", "HTI", "HB")) +
  # scale_linetype_manual(name = "",values = c("Air Temperature" = "solid")) +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(panel.grid = element_blank()) +
  labs(x = "Days after planting")

print(canopy_temp_plot)
save_plot(canopy_temp_plot, "FigS9_canopy_temp_timeseries", height = 4, width = 5)

























































