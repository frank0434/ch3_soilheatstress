source("scripts/00_setup.R")

# ==============================================================================
# LOAD & PREPARE DATA
# ==============================================================================
li_td_sum <- fread("data/LI_cleaned.csv")

# Thermal days at leaf removal: used as te upper bound during curve fitting
# (te cannot exceed leaf removal since canopy is forcibly terminated at this point)
leaf_removal_td <- data.table(
  Season = c(2024, 2025),
  te_fixed = c(59, 58))
li_td_sum <- merge(li_td_sum, leaf_removal_td, by = "Season", all.x = TRUE)

li_td_sum[, ':='(
  Season = factor(Season), 
  Treatment = factor(Treatment, levels = legend_order), 
  Block = factor(Block),
  PlotID = factor(PlotID)
)]
setorder(li_td_sum, Season, Treatment, Block, PlotID)

li_plot_fits <- li_td_sum[
  LI > 0 & is.finite(cumtd_DAE),
  .(data = list(.SD)),
  by = .(Season, Treatment, Block, PlotID, te_fixed)
]
setorder(li_plot_fits, Season, Treatment, Block, PlotID)

# ==============================================================================
# PER-PLOT CURVE FITTING
# ==============================================================================
# DEoptim finds robust starting values; nlsLM fits the combined canopy equation
li_plot_fits[, starts := mapply(
  function(x, te) optimize_start_values_combined(x, y_col = "LI", x_col = "cumtd_DAE"),
  data,
  SIMPLIFY = FALSE)]
# fitting with nlsLM with bounds
li_plot_fits[, fit := mapply(function(data, start, te) 
  fit_combined(data, y_col = "LI", x_col = "cumtd_DAE", start_vals = start), 
  data, starts,SIMPLIFY = FALSE)]
# Calculate R² per plot fit (stored in li_plot_fits for optional diagnostics)
li_plot_fits[, r2 := sapply(fit, r2_nls)]
li_plot_params <- li_plot_fits[, .(sapply(fit, extract_combined_params, simplify = FALSE)),
   by = .(Season, Treatment, Block, PlotID)
  ][, unlist(V1, recursive = FALSE), by = .(Season, Treatment, Block, PlotID)]

prediction_grid <- generate_time_grid(li_td_sum[!is.na(cumtd_DAE)], end_time = 65, step = 0.1)
predictions_all <- merge(
  prediction_grid,
  li_plot_params[, .(Season, Treatment, Block, PlotID, tm1, t1, t2, te, vmax)],
  by = c("Season", "Treatment", "Block", "PlotID"),
  all.x = TRUE
)
predictions_all[, prediction := combined_equation(cumtd_DAE, tm1, t1, t2, te, vmax)]
predictions_all[!is.finite(prediction) | prediction < 0, prediction := 0]

li_plot_params[, delta_t := t2 - t1]  # duration of plateau phase (thermal days)
canopy_param_cols <- c("vmax", "t1", "tm1", "t2", "te", "delta_t")

# ==============================================================================
# PARAMETER-LEVEL STATISTICAL MODELS
# ==============================================================================
# Separate linear models per season; treatment effects on each canopy parameter
li_param_models_24_lm <- setNames(lapply(canopy_param_cols, function(param) {
  fit_treatment_lm(li_plot_params[Season == 2024], param)
}), canopy_param_cols)
li_param_models_25_lm <- setNames(lapply(canopy_param_cols, function(param) {
  fit_treatment_lm(li_plot_params[Season == 2025], param)
}), canopy_param_cols)
shapiro.test(residuals(li_param_models_24_lm$te$model))

shapiro.test(residuals(li_param_models_25_lm$te$model))

param_LI_tables_24 <- rbindlist(lapply(canopy_param_cols, function(param) {
  summarise_emmeans_table(
    emm = li_param_models_24_lm[[param]]$emmeans,
    season = 2024,
    model_name = "LI_deoptim_nls",
    parameter = param,
    output_component = "treatment_means"
  )
})) 

param_LI_pair_24 <- rbindlist(lapply(canopy_param_cols, function(param) {
  summarise_pairwise_table(
    emm = li_param_models_24_lm[[param]]$emmeans,
    season = 2024,
    model_name = "LI_deoptim_nls",
    parameter = param
  )
}))

param_LI_tables_25 <- rbindlist(lapply(canopy_param_cols, function(param) {
  summarise_emmeans_table(
    emm = li_param_models_25_lm[[param]]$emmeans,
    season = 2025,
    model_name = "LI_deoptim_nls",
    parameter = param,
    output_component = "treatment_means"
  )
})) 

param_LI_pair_25 <- rbindlist(lapply(canopy_param_cols, function(param) {
  summarise_pairwise_table(
    emm = li_param_models_25_lm[[param]]$emmeans,
    season = 2025,
    model_name = "LI_deoptim_nls",
    parameter = param
  )
}))
rbindlist(list(param_LI_pair_24[,.(contrast, p.value, Season, Parameter, Significant)], 
               param_LI_pair_25[,.(contrast, p.value, Season, Parameter, Significant)]), use.names = TRUE, fill = TRUE) |> 
  fwrite("results/LI_combined_model_parameters_pairwise_comparisons.csv", bom = TRUE)

# Test whether the HTI treatment effect on canopy parameters (t2, te) differs
# between seasons — i.e., does the timing of heat stress interact with season?
hti_ac_li <- li_plot_params[Treatment %in% c("HTI", "AC")]
hti_ac_li_mod <- lm(t2 ~ Block + Treatment * Season, data = hti_ac_li)
shapiro.test(residuals(hti_ac_li_mod))
anova(hti_ac_li_mod)
TukeyHSD(aov(t2 ~ Block + Treatment * Season, data = hti_ac_li), "Treatment", adjust = "tukey")
hti_ac_li_mod <- lm(te ~ Block + Treatment * Season, data = hti_ac_li)
shapiro.test(residuals(hti_ac_li_mod))
anova(hti_ac_li_mod)
TukeyHSD(aov(te ~ Block + Treatment * Season, data = hti_ac_li), "Treatment", adjust = "tukey")
hti_ac_li[, mean(te), by = .( Season)]
# supplementary table with the parameter estimates for each plot
li_plot_params_table <- li_plot_params[, .(Season, Treatment, Block, PlotID, tm1, t1, t2, te, vmax, delta_t)] |> 
  melt(id.vars = c("Season", "Treatment", "Block", "PlotID"))
li_params_summary <- li_plot_params_table[, .(mean = mean(value), se = sd(value)/sqrt(.N)),
   by = .(Season, Treatment, variable)
  ][, formatted := sprintf("%.2f ± %.2f", mean, se)] 
li_params_summary[variable == "vmax", formatted := sprintf("%.2f ± %.2f", mean, se)]
li_params_summary[variable != "vmax", formatted := sprintf("%.0f ± %.1f", mean, se)]
li_params_summary |> 
  dcast(Season + Treatment ~ variable, value.var = "formatted") |> 
  fwrite("results/LI_combined_model_parameters_summary.csv", bom = TRUE)

# ==============================================================================
# TREATMENT-LEVEL FITTING & PREDICTIONS
# ==============================================================================
# Re-fit the canopy curve at treatment level (pooled data) for plotting smooth
# mean curves; plot-level parameters above are used for statistical testing
li_trt_fits <- li_td_sum[
  LI > 0 & is.finite(cumtd_DAE),
  .(data = list(.SD)),
  by = .(Season, Treatment, te_fixed)]
setorder(li_trt_fits, Season, Treatment)

li_trt_fits[, starts := mapply(
  function(x) optimize_start_values_combined(x, y_col = "LI", x_col = "cumtd_DAE"),
  data,
  SIMPLIFY = FALSE)]
# fitting with nlsLM with bounds
li_trt_fits[, fit := mapply(function(data, start) 
  fit_combined(data, y_col = "LI", x_col = "cumtd_DAE", start_vals = start), 
  data, starts,SIMPLIFY = FALSE)]
# Calculate R² for treatment-level fits (stored in li_trt_fits for diagnostics)
li_trt_fits[, r2 := sapply(fit, r2_nls)]
li_trt_params <- li_trt_fits[, .(sapply(fit, extract_combined_params, simplify = FALSE)),
   by = .(Season, Treatment)][, unlist(V1, recursive = FALSE), by = .(Season, Treatment)]
li_trt_params[, curve_dat := mapply(function(tm1, t1, t2, te, vmax) {
  time <- seq(1, 60, 0.1)
  y <- combined_equation(time, tm1, t1, t2, te, vmax)
  data.table(cumtd_DAE = time, predictions = y)
}, tm1, t1, t2, te, vmax, SIMPLIFY = FALSE)]

LI_trt_curve_dt <- li_trt_params[, unlist(curve_dat, recursive = FALSE), by = .(Season, Treatment)]
# summarise raw data for plotting
li_summary <- li_td_sum[, .(mean_LI = mean(LI, na.rm = TRUE),
                      se_LI = sd(LI, na.rm = TRUE) / sqrt(.N)),
                  by = .(Season, Treatment, cumtd_DAE)][mean_LI > 0 & is.finite(cumtd_DAE)]

mean_LI_curve <- LI_trt_curve_dt[predictions >= 0
  ][, cumtd_DAE := round(cumtd_DAE) |> as.integer()
  ][, .(predictions = mean(predictions, na.rm = TRUE)), 
    by = .(Season, Treatment, cumtd_DAE)
  ][, ':=' (Season = as.integer(as.character(Season)))]
# Add corresponding DAP values back from daily_td to enable a dual x-axis
# (thermal days primary; DAP shown as secondary labels in parentheses)
daily_td <- fread("data/daily_td.csv")
x_for_fPAR <- daily_td[,.(Season, DAP, cumtd_DAE = round (cumtd_DAE))
  ][ !is.na(cumtd_DAE)][, .SD[1], by = .(Season, cumtd_DAE)]
mean_LI_curve <- merge(mean_LI_curve, x_for_fPAR, by = c("Season", "cumtd_DAE"), all.x = TRUE)
dap_axis <- unique(
  mean_LI_curve[
    cumtd_DAE %in% seq(0, 60, 10),
    .(Season, cumtd_DAE, DAP)
  ]
)
dap_axis <- rbindlist(list(dap_axis,
   daily_td[DAE == 0][,.(Season, DAP, cumtd_DAE = 0)]
), use.names = TRUE, fill = TRUE)
rect_td_min <- daily_td[treatment_rect_dt, on = c("Season", "DAP == DAP_min")
                         ][, .(Season, Treatment, DAP_min = DAP, DAP_max, xmin = cumtd_DAE)]
rect_td <- daily_td[rect_td_min, on = c("Season", "DAP == DAP_max")
                     ][, .(Season, Treatment, xmin, xmax = cumtd_DAE)]
p_fPAR <- mean_LI_curve |>
  ggplot(aes(cumtd_DAE, predictions, colour = Treatment)) +
  geom_rect(
    data = rect_td,
    aes(xmin =  xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colors_temp, 
                    breaks = legend_order, guide = "none") +
  geom_line(aes(linetype = Treatment), linewidth = lw ) +
  geom_point(data = li_summary, aes(x = cumtd_DAE, y = mean_LI), size = 1, alpha = 0.5) +
  geom_errorbar(data = li_summary, aes(x = cumtd_DAE, ymin = mean_LI - se_LI, ymax = mean_LI + se_LI), 
                width = 1, alpha = 0.5, inherit.aes = FALSE) +
  facet_grid(. ~ Season) +
  geom_text(
    data = dap_axis,
    aes(x = cumtd_DAE, y = 0, label = paste0("(", DAP, ")")),
    inherit.aes = FALSE,
    vjust = 3.2, size = floor(fontsize / ggplot2::.pt) - 1 ) +
  scale_color_manual(name = "Treatment", values = colors_temp, breaks = legend_order) +
  scale_x_continuous(name = "Thermal days\n(Days after planting)", limits = c(0, 65), 
                     breaks = seq(0, 65, 10)) +
  scale_linetype_manual(name = "Treatment", values = linetype_temp, breaks = legend_order) +
  scale_y_continuous(name = "Fraction of light intercepted", limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0,0)) +
    coord_cartesian(clip = "off") +
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(
    legend.position = "top",
    legend.box.margin = margin(-1, -1, -1, 0, "mm"),
    legend.key.size = unit(15, "mm"),
    legend.margin = margin(1, 0, 0, 0, "mm"),
    legend.background = element_blank(),
    legend.title = element_text(
      margin = margin(r = 6, unit = "mm")),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1, "mm"),
    panel.spacing.y = unit(-10, "mm"),
    plot.margin = margin(1, 1, 12, 1, "mm"),
    axis.title.x = element_text(vjust = - 3)
  )

save_plot(p_fPAR, "fig_6_light_interception", height = 4)

# ==============================================================================
# RADIATION USE EFFICIENCY (RUE)
# ==============================================================================
# Daily intercepted PAR (Rint) = predicted LI x daily radiation (IRRAD)
# Cumulative intercepted PAR (Rcum) used as the x-axis in RUE regressions
# RUE (slope) is estimated per plot as: DW ~ PAR, then tested with ANOVA
# RUE ---------------------------------------------------------------------
# total Light interception  -----------------------------------------------
daily_radiation_data <- fread("data/daily_weather.csv", skip = 1)
daily_radiation_data[, Date := as.Date(Date, format = "%d/%m/%Y")]
daily_td_Rad <- merge.data.table(daily_td[DAE >= 0], daily_radiation_data,
                                    by = c("Season","Date"), all.x = TRUE)
# 2024-08-24 no radiation data, fill with the last observation carried forward
daily_td_Rad[, IRRAD := nafill(IRRAD, type = "locf")]
data_td <- daily_td_Rad[, .(data_daily = list(.SD)) , by = .(Season)]
li_plot_fits[, Season := as.integer(as.character(Season))]

daily_li <- merge(li_plot_fits, data_td, by = "Season", all.x = TRUE)
daily_li[, prediction := mapply(function(fit, data) {
  data[, prediction := predict(fit, newdata = data)]
  data[prediction < 0, prediction := 0]
  return(data)
}, fit, data_daily, SIMPLIFY = FALSE)]
daily_li_rad <- daily_li[, unlist(prediction, recursive = FALSE),
   by = .(Season, Treatment, Block, PlotID)]
daily_li_rad[, Rint := prediction * IRRAD]
daily_li_rad[, Rcum := cumsum(Rint), 
  by = .(Season, Treatment, Block, PlotID)]

Rcum_24 <- daily_li_rad[Season == 2024, .SD[.N], by = .(Season, Treatment, Block, PlotID)]
model <- lm(Rcum ~ Block + Treatment, Rcum_24)
anova(model)
summary(model)
shapiro.test(residuals(model))
Rcum_25 <- daily_li_rad[Season == 2025, .SD[.N], by = .(Season, Treatment, Block, PlotID)]
model <- lm(Rcum ~ Block + Treatment, Rcum_25)
Anova(model)
shapiro.test(residuals(model))

# TOTAL light interception by the end of the season (Rcum) is 
# the cumulative sum of daily intercepted radiation (Rint) for each plot, 
# which is calculated as the product of the predicted light interception (LI) 
# and the daily radiation (IRRAD). This gives us an estimate of the total amount
#  of light energy intercepted by the crop canopy over the growing season, 
# which can be used to calculate radiation use efficiency (RUE) when combined 
# with biomass production data.


# Bring over the total weight of the crop  --------------------------------
dma_data <- fread("data/biomass_cleaned.csv")

dma_data[, ':='(
  Treatment = factor(Treatment, levels = legend_order),
  Block = factor(Block),
  PlotID = factor(PlotID)
)]
biomass_Rcum <- merge(dma_data, daily_li_rad, all.x = TRUE,
   by = c("Season", "Treatment", "Block", "PlotID","Date"))

# Labels for equation + R2 (per treatment, simple y~x for visual annotation)
biomass_Rcum[, `:=`(DM_g_m2 = total_dw_g_m2, PAR = Rcum / 2)]
ann <- biomass_Rcum[, {
  mm <- lm(total_dw_g_m2 ~ PAR, data = .SD)
  b0 <- coef(mm)[1]
  b1 <- coef(mm)[2]
  r2 <- summary(mm)$r.squared
  .(
    label = sprintf("DW = %.3f %+.2fPAR\nR² = %.3f", b0, b1, r2),
    x = 100,
    y = 1000,
    b0 = b0,
    b1 = b1,
    r2 = r2
  )
}, by = .(Season, Treatment, Block, PlotID)]
# test the slopes or the real RUE 
RUE_24 <- lm(b1 ~ Block +  Treatment, data = ann[Season == 2024])
anova(RUE_24)
shapiro.test(residuals(RUE_24))
RUE_25 <- lm(b1 ~ Block +  Treatment, data = ann[Season == 2025])
Anova(RUE_25)
shapiro.test(residuals(RUE_25))
# acoss season 
RUE_all <- lm(b1 ~ Block + Season * Treatment, data = ann)

Anova(RUE_all)
shapiro.test(residuals(RUE_all))
# output 
RUE_summary <- calculate_summary_stats(ann, "b1", c("Season", "Treatment"))[]
RUE_summary[, formatted := sprintf("%.2f ± %.2f", mean, se)] |> 
  fwrite("results/RUE_summary.csv", bom = TRUE)
# ==============================================================================
# RUE FIGURE (Supplementary)
# ==============================================================================
# fcase() label positions are treatment/season-specific to avoid overlapping
# annotation text on the plot
ann_trt <- biomass_Rcum[, {
  mm <- lm(total_dw_g_m2 ~ PAR, data = .SD)
  b0 <- coef(mm)[1]
  b1 <- coef(mm)[2]
  r2 <- summary(mm)$r.squared
  .(
    label = sprintf("DW = %.0f %+.2fPAR\nR² = %.2f", b0, b1, r2),
    x = 10,
    y = 1000,
    b0 = b0,
    b1 = b1,
    r2 = r2
  )
}, by = .(Season, Treatment)]

ann_trt[, y := fcase(Treatment == "HE", 1400, 
                 Treatment == "HTI" & Season == 2024, 1200,
                 Treatment == "HTI" & Season == 2025, 1400,
                 Treatment == "HB", 1200,
                 Treatment == "AC" & Season == 2024, 1600,
                 Treatment == "AC" & Season == 2025, 1600)]
figS_RUE <- ggplot(biomass_Rcum, aes(x = PAR, y = DM_g_m2, color = Treatment)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  # facet_wrap(~Treatment, scales = "free") +
  geom_text(
    data = unique(ann_trt),
    aes(x = x, y = y+ 100, label = paste(Treatment, label)),
    hjust = 0, vjust = 1, color = "black", size = 3,
    inherit.aes = FALSE
  ) +
  labs(
    x = expression(Accumulated~intercepted~PAR~(MJ~~m^{-2})),
    y = expression(Total~dry~weight~(g~~m^{-2}))) +
  facet_grid(~Season)+
  scale_color_manual(values = colors_temp, breaks = legend_order) +
  theme_bw() +
  theme(panel.grid = element_blank())

save_plot(figS_RUE, "figS_RUE", width = 7, height = 4)
