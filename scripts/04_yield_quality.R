source("scripts/00_setup.R")

# ==============================================================================
# 1. TUBER DRY MATTER CONTENT (DMC) - Final Harvest
# ==============================================================================

dmc <- fread("data/dmc_final_harvest.csv")
dmc[, tuber_dmc := tuber_dmc * 100]  # Convert to percentage
dmc[, `:=`(PlotID = as.factor(PlotID), Treatment = factor(Treatment))]

final_data_2024 <- dmc[Season == 2024]
final_data_2024 <- droplevels(final_data_2024)
tuber_model_2024 <- lm(tuber_dmc ~ Treatment,
                        data = final_data_2024)
anova(tuber_model_2024)

# 2025 DMC analysis ----
final_data_2025 <- dmc[Season == 2025 & PlotID != 8]
final_data_2025 <- droplevels(final_data_2025)

tuber_model_2025 <- lm(tuber_dmc ~ Treatment,
                        data = final_data_2025)
Anova(tuber_model_2025, type = "II")
# checked DMC for across seasons with lm - consistent with lme results
dmc_ac_hti <- droplevels(dmc[Treatment %in% c("AC", "HTI")])
dmc_model_ac_hti <- lm(tuber_dmc ~ Block + Treatment + Season + Treatment:Season, data = dmc_ac_hti)
Anova(dmc_model_ac_hti, type = "II")
shapiro.test(dmc_model_ac_hti$residuals)

# 2. TUBER SIZE DISTRIBUTION - Weibull Scale Parameter (2025 only) -------------
# Note: Weibull fitting only available for Season 2025 
# (2024 data lacks individual tuber measurements)
# plot 8 already taken out 
weibull_fit_results <- fread("data/weibull_fit_results_2025.csv")

weibull_params <- weibull_fit_results[!is.na(shape) & !is.na(scale)]
scale_model_2025 <- lm(scale ~ Treatment, data = weibull_params)
shapiro.test(scale_model_2025$residuals)  # Check normality of residuals
shape_model_2025 <- lm(shape ~ Treatment, data = weibull_params)
shapiro.test(shape_model_2025$residuals)  # Check normality of residuals
Anova(shape_model_2025)
# emmeans(scale_model_2025, pairwise ~ Treatment) is called inline where needed

# ==============================================================================
# SUPPLEMENTARY FIGURES: TSD and Deformation
# ==============================================================================
# Create publication-quality figures for tuber size distribution and deformation

# --- Prepare deformation data for visualization -----
tuber_q <- fread("data/tuber_quality_prepared.csv")
tuber_q[, total_w := sum(Weight.FW, na.rm = TRUE), by = .(PlotID)]
tuber_q[, pct_w := Weight.FW / total_w]
tuber_q[, pct_in_integer := round(pct_w * 100)]

# Filter deformed tubers (Secondary and Others), exclude Plot 8
deformation <- tuber_q[Classes %in% c("Secondary", "Others") & !PlotID %in% 8,
                       .(pct_deformed = pct_in_integer),
                       by = .(PlotID, Treatment, Block, Classes)]

deformation[, `:=`(
  Treatment = factor(Treatment),
  Classes = factor(Classes, levels = c("Secondary", "Others"))
)]

# Create significance letters for Secondary class (from previous analysis)
sig_letters <- data.table(
  Treatment = factor(c("AC", "HTI", "HB"), levels = c("AC", "HTI", "HB")),
  Classes = factor("Second\nGrowth", levels = c("Second\nGrowth", "Others")),
  label = c("a", "b", "a")  # HTI significantly different from AC and HB
)

deformation_stats <- deformation[Classes != "Others", .(
  mean = mean(pct_deformed),
  max = max(pct_deformed)
), by = Treatment]
sig_letters <- sig_letters[deformation_stats, on = "Treatment"]

# --- Load TSD visualization data ---
weibull_density_curves <- fread("data/weibull_density_curves.csv")
tuber_size_bounds_with_midpoints <- fread("data/tuber_size_bounds_with_midpoints.csv")

# --- Create Panel A: TSD histogram with Weibull density overlay ---
size_summary <- calculate_summary_stats(tuber_size_bounds_with_midpoints, "pct_in_integer", 
                                  c("Treatment", "cls", "Midpoint"))[]
size_summary[, Treatment := factor(Treatment, levels = c("AC", "HTI", "HB"))]

seasonal_curves <- calculate_summary_stats(weibull_density_curves, "density", 
                                     c("Treatment", "x_axis"))[]
seasonal_curves[, density := mean]
seasonal_curves[, Treatment := factor(Treatment, levels = c("AC", "HTI", "HB"))]

# Scale factor to align histogram and density curve axes
scale_factor <- max(size_summary$mean, na.rm = TRUE) / 
                max(seasonal_curves$density, na.rm = TRUE)

fig_TSD_w <- seasonal_curves |>
  ggplot() +
  geom_col(data = size_summary,
           aes(x = Midpoint, y = mean, fill = Treatment),
           position = position_dodge(width = 5), width = 4, alpha = 0.5) +
  geom_errorbar(data = size_summary,
                aes(x = Midpoint,
                    ymin = pmax(mean - se, 0),
                    ymax = mean + se,
                    color = Treatment),
                position = position_dodge(width = 5), width = 3) +
  geom_line(aes(x = x_axis, y = density * scale_factor, linetype = Treatment, color = Treatment),
            linewidth = lw) +
  scale_x_continuous(breaks = seq(0, 110, 20), limits = c(0, 110)) +
  scale_y_continuous(
    name = "Fresh yield (%)", expand = c(0, 0), limits = c(0, 65),
    sec.axis = sec_axis(~ . / scale_factor, name = "Density")) +
  labs(x = "Tuber size (mm)", tag = "A") +
  scale_color_manual(values = colors_temp, name = "Treatment") +
  scale_fill_manual(values = colors_temp, name = "Treatment") +
  scale_linetype_manual(values = linetype_temp, name = "Treatment") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.key.width = unit(15, "mm"),
    plot.tag = element_text(face = "bold", size = 14)
  )

# --- Create Panel B: Deformation boxplot with significance letters ---
dodge_width <- 0.75
deformation[, Classes := ifelse(Classes == "Secondary", "Second\nGrowth", "Others")]
fig_deformation <- deformation |>
  ggplot(aes(x = Classes, y = pct_deformed, fill = Treatment)) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = 16, 
               position = position_dodge(width = dodge_width)) +
  geom_text(data = sig_letters, 
            aes(x = Classes, y = max + 1.5, label = label, group = Treatment),
            size = 5, fontface = "bold", 
            position = position_dodge(width = dodge_width)) +
  scale_y_continuous(
    name = "Fresh yield (%)", 
    expand = c(0, 0), limits = c(0, 65)) +
  labs(x = "Tuber class", tag = "B") +
  scale_fill_manual(values = colors_temp, name = "") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.key.width = unit(15, "mm"),
    plot.tag = element_text(face = "bold", size = 14)
  )

# --- Combine and save figures ---
fig_combined <- (fig_TSD_w + plot_spacer()+(fig_deformation + theme(legend.position = "none"))) + 
  plot_layout(widths = c(2, -0.15,0.6))

print(fig_combined)
save_plot(fig_combined, "figS_TSD_deformation_combined", width = 7, height = 5)

# ==============================================================================
# 3. MARKETABLE YIELD AND SECONDARY GROWTH PERCENTAGES
# ==============================================================================
# Market: ≥50 mm tubers; Secondary: deformed/secondary growth; Undersized: <50 mm
# AC/HTI cross-season model tests whether the HTI effect is consistent
# across both seasons (Treatment:Season interaction)
tuber_w_pct <- fread("data/tuber_w_pct_combined.csv")

market_2024 <- droplevels(tuber_w_pct[Season == "2024"])
market_model_2024 <- lm(market_pct ~ Treatment, data = market_2024)
shapiro.test(market_model_2024$residuals)  # Check normality of residuals
results_market_2024 <- extract_emmeans_results(market_model_2024, "2024")

market_2025 <- droplevels(tuber_w_pct[Season == "2025" & PlotID != 8])
market_model_2025 <- lm(market_pct ~ Treatment, data = market_2025)
shapiro.test(market_model_2025$residuals)  # Check normality of residuals
results_market_2025 <- extract_emmeans_results(market_model_2025, "2025")

# across season analysis for marketable % with lm 
market_ac_hti <- droplevels(tuber_w_pct[Treatment %in% c("AC", "HTI")])
market_ac_hti_model <- lm(market_pct ~ Block + Treatment + Season + Treatment:Season, data = market_ac_hti)

Anova(market_ac_hti_model, type = "II")
shapiro.test(market_ac_hti_model$residuals)
# Secondary % analysis ----
sec_2024 <- droplevels(tuber_w_pct[Season == "2024"])
sec_model_2024 <- lm(Secondary_pct ~ Treatment, data = sec_2024)
shapiro.test(sec_model_2024$residuals)  # Check normality of residuals
results_sec_2024 <- extract_emmeans_results(sec_model_2024, "2024")
sec_2025 <- droplevels(tuber_w_pct[Season == "2025" & PlotID != 8])
sec_model_2025 <- lm(Secondary_pct ~ Treatment, data = sec_2025)
shapiro.test(sec_model_2025$residuals)  # Check normality of residuals
results_sec_2025 <- extract_emmeans_results(sec_model_2025, "2025")

# analysis HTI vs AC for secondary growth with season 
sec_ac_hti <- droplevels(tuber_w_pct[Treatment %in% c("AC", "HTI")])
sec_ac_hti_model <- lm(Secondary_pct ~ Block + Season + Treatment + Season:Treatment, data = sec_ac_hti)
Anova(sec_ac_hti_model, type = "II")
shapiro.test(sec_ac_hti_model$residuals)
# 4. AVERAGE TUBER WEIGHT (AVT) -----------------------------------------------
# AVT = Total tuber weight / Total tuber number per plot
# Also derived: tuber number per m² (density)
# AC/HTI cross-season model checks whether HTI effect on AVT is consistent
tuber_avt <- fread("data/tuber_avt_combined.csv")
tuber_avt[, avt := weight_of_tubers / number_of_tubers]
tuber_avt[, tuber_n_per_m2 := number_of_tubers / (number_of_plants * 0.75 * 0.3)]
tuber_avt[, `:=`(Season = factor(Season), Treatment = factor(Treatment), 
                 Block = as.factor(Block), PlotID = as.factor(PlotID))]

avt_2024 <- droplevels(tuber_avt[Season == "2024"])
avt_model_2024 <- lm(avt ~  Treatment, data = avt_2024)
anova(avt_model_2024)
shapiro.test(avt_model_2024$residuals)  # Check normality of residuals
results_avt_2024 <- extract_emmeans_results(avt_model_2024, "2024")

avt_2025 <- droplevels(tuber_avt[Season == "2025"])
avt_model_2025 <- lm(avt ~ Treatment, data = avt_2025)
summary(avt_model_2025)
Anova(avt_model_2025, type = "II")
shapiro.test(avt_model_2025$residuals)  # Check normality of residuals
results_avt_2025 <- extract_emmeans_results(avt_model_2025, "2025")

# checked AVT for across seasons with lm - consistent with lme results
avt_ac_hti <- droplevels(tuber_avt[Treatment %in% c("AC", "HTI")])
avt_model_ac_hti <- lm(avt ~ Block + Treatment + Season + Treatment:Season, data = avt_ac_hti)
summary(avt_model_ac_hti)
Anova(avt_model_ac_hti, type = "II")
shapiro.test(avt_model_ac_hti$residuals)

# tuber number 
number_2024 <- droplevels(tuber_avt[Season == "2024"])
number_model_2024 <- lm(tuber_n_per_m2 ~ Block + Treatment, data = number_2024)
summary(number_model_2024)
results_number_2024 <- extract_emmeans_results(number_model_2024, "2024")

number_2025 <- droplevels(tuber_avt[Season == "2025"])
number_model_2025 <- lm(tuber_n_per_m2 ~ Block + Treatment, data = number_2025)
results_number_2025 <- extract_emmeans_results(number_model_2025, "2025")
tuber_avt[, .(tuber_n_m2 = mean(tuber_n_per_m2),
              se = sd(tuber_n_per_m2)/.N^0.5), by = .(Season, Treatment)][order(Season, Treatment)]

# Tuber number - AC and HTI only with Season as random effect
number_ac_hti <- droplevels(tuber_avt[Treatment %in% c("AC", "HTI")])
number_model_ac_hti <- lm(tuber_n_per_m2 ~ Block + Season + Treatment + Season:Treatment, 
               data = number_ac_hti)
shapiro.test(number_model_ac_hti$residuals)
Anova(number_model_ac_hti, type = "II")
results_number_ac_hti <- extract_emmeans_results(number_model_ac_hti, "AC_HTI")

# 5. FRYING COLOUR QUALITY------------------------------------------------------
# Lower frying index = better colour quality (darker fry = higher index)
# AC/HTI cross-season model assesses whether HTI effect on fry quality persists
# across seasons
frying_colour_data <- fread("data/frying_colour_data.csv")

fc_2024 <- droplevels(frying_colour_data[Season == "2024" & !is.na(Fry.index)])
fc_model_2024 <- lm(Fry.index ~ Treatment, data = fc_2024)
anova(fc_model_2024)
fc_2025 <- droplevels(frying_colour_data[Season == "2025" & PlotID != 8 & !is.na(Fry.index)])
fc_model_2025 <- lm(Fry.index ~ Treatment, data = fc_2025)
Anova(fc_model_2025)

# testing frying colour with lm across seasons with AC and HTI only - consistent with lme results
fc_ac_hti <- droplevels(frying_colour_data[Treatment %in% c("AC", "HTI") & !is.na(Fry.index)])
fc_model_ac_hti <- lm(Fry.index ~ Block + Treatment + Season + Treatment:Season, data = fc_ac_hti)
summary(fc_model_ac_hti)
Anova(fc_model_ac_hti, type = "II")
emm_fc_ac_hti <- emmeans(fc_model_ac_hti, pairwise ~ Treatment | Season)

# ==============================================================================
# 6. FRESH TUBER YIELD (t/ha)
# ==============================================================================
# Fresh yield derived from AVT dataset; converted to g/m² using plot plant area
# AC/HTI cross-season model included to test consistency of yield response
# 0. FRESH TUBER YIELD (t/ha) -------------------------------------------------
tuber_avt[, area := number_of_plants * 0.75 * 0.3]
tuber_avt[, fresh_yield_g_m2 := weight_of_tubers / area ]  # Convert to kg/ha

yield_data <- tuber_avt[, .(Season, Block, Treatment, PlotID, fresh_yield_g_m2)]

yield_model_2024 <- lm(fresh_yield_g_m2 ~ Block + Treatment,
                        data    = droplevels(yield_data[Season == 2024]))
summary(yield_model_2024)
anova(yield_model_2024)
results_yield_2024 <- extract_emmeans_results(yield_model_2024, "2024")

yield_model_2025 <- lm(fresh_yield_g_m2 ~ Block + Treatment,
                        data    = droplevels(yield_data[Season == 2025 & PlotID != 8]))
summary(yield_model_2025)
Anova(yield_model_2025, type = "II")
results_yield_2025 <- extract_emmeans_results(yield_model_2025, "2025")
# checked raw means and se for fresh yield - consistent with biomass results
calculate_summary_stats(yield_data, "fresh_yield_g_m2", c("Season","Treatment"))[]

# tuber lm fresh yield model
tuber_fresh_yield_model <- lm(fresh_yield_g_m2 ~ Season + Block +  Treatment + Season:Treatment, 
  data = yield_data[Treatment %in% c("AC", "HTI")])
summary(tuber_fresh_yield_model)
Anova(tuber_fresh_yield_model, type = "II")
shapiro.test(tuber_fresh_yield_model$residuals)
calculate_summary_stats(yield_data[Treatment %in% c("AC", "HTI")], "fresh_yield_g_m2", c("Season"))[]

