# ==============================================================================
# GENERIC PLOTTING UTILITIES
# ==============================================================================

#' Save Publication-Quality Plot
#'
#' Consistent method for saving plots with standard dimensions and DPI
#'
save_plot <- function(plot, name, width = 7, height = 5, dpi = 300, format = "svg") {
  filename <- paste0("figures/", name, ".", format)
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}

# ==============================================================================
# MISSING VALUE HANDLING
# ==============================================================================

#' Replace NA with 0
#'
#' Simple utility to replace missing values with zero.
#'
#' @param col numeric vector
#'
#' @return numeric vector with NA replaced by 0
#'
#' @keywords internal
replace_NA <- function(col) {
  ifelse(is.na(col), 0, col)
}

#' Replace 0 with NA
#'
#' Replace zero values with missing value indicator.
#'
#' @param col numeric vector
#'
#' @return numeric vector with 0 replaced by NA
#'
#' @keywords internal
replace_0 <- function(col) {
  ifelse(col == 0, NA, col)
}
# ==============================================================================
# SUMMARY STATISTICS AND DATA TRANSFORMATION
# ==============================================================================
#' Calculate Summary Statistics
#' 
#' @param dt data.table - input data
#' @param value_col character - column name to calculate statistics for
#' @param group_vars character vector - grouping variables
#' @return data.table with mean, se, sd, n, and formatted columns
calculate_summary_stats <- function(dt, value_col, group_vars) {
  summary_dt <- dt[, .(
    mean = mean(get(value_col), na.rm = TRUE),
    se = sd(get(value_col), na.rm = TRUE) / sqrt(sum(!is.na(get(value_col)))),
    sd = sd(get(value_col), na.rm = TRUE),
    n = sum(!is.na(get(value_col)))
  ), by = group_vars]
  
  summary_dt[, formatted := sprintf("%.0f ± %.1f", mean, se)]
  return(summary_dt)
}
fit_treatment_lm <- function(data, response) {
  model <- stats::lm(
    formula = stats::as.formula(paste(response, "~ Block + Treatment")),
    data = data
  )
  emm <- emmeans::emmeans(model, specs = ~ Treatment)
  list(
    model = model,
    anova = anova(model),
    emmeans = emm,
    pairs_tukey = pairs(emm, adjust = "tukey")
  )
}


# fitting curves for light interception data 

#' Optimize Starting Values for Combined Canopy Equation using DEoptim
#'
#' Uses Differential Evolution to find stable starting values for the
#' `combined_equation()` canopy model before fitting with `nlsLM` or `nlme`.
#'
#' @param clean_data data.table with response and thermal time columns
#' @param global_vmax numeric optional upper-scale reference for `vmax`
#' @param y_col character response column name
#' @param x_col character thermal time column name
#' @param seed numeric optional seed for reproducibility
#'
#' @return list with optimized tm1, t1, t2, te, and vmax values
optimize_start_values_combined <- function(clean_data,
                                          global_vmax = NULL,
                                          y_col = "CC",
                                          x_col = "cumtd_DAE",
                                          tmax = 58, # tmax cannot exceed leaf removal cumtd_DAE, which is 58 in 2024 and 2025
                                          seed = NULL) {

  required_cols <- c(y_col, x_col)
  if (!all(required_cols %in% names(clean_data))) {
    stop(sprintf("Missing required columns: %s", paste(setdiff(required_cols, names(clean_data)), collapse = ", ")))
  }

  dat <- copy(clean_data)[!is.na(get(y_col)) & !is.na(get(x_col))]
  if (nrow(dat) < 5) {
    return(NULL)
  }

  setorderv(dat, x_col)

  if (is.null(seed)) {
    seed <- as.integer(abs(sum(dat[[y_col]], na.rm = TRUE) * 1000) %% 10000)
  }
  set.seed(seed)

  observed_vmax <- max(dat[[y_col]], na.rm = TRUE)
  if (!is.finite(observed_vmax)) {
    return(NULL)
  }

  if (is.null(global_vmax) || !is.finite(global_vmax) || global_vmax <= 0) {
    global_vmax <- observed_vmax
  }

  x_min <- min(dat[[x_col]], na.rm = TRUE)
  x_max <- max(dat[[x_col]], na.rm = TRUE)
  x_range <- x_max - x_min
  if (!is.finite(x_range) || x_range <= 0) {
    return(NULL)
  }

  min_gap <- max(5, x_range * 0.1)

  tm1_lower <- max(0, x_min)
  tm1_upper <- max(tm1_lower + min_gap, x_min + x_range * 0.35)

  t1_lower <- tm1_lower + min_gap
  t1_upper <- max(t1_lower + min_gap, x_min + x_range * 0.6)

  t2_lower <- t1_lower + min_gap
  t2_upper <- max(t2_lower + min_gap, x_min + x_range * 0.9)

  te_lower <- t2_lower + min_gap
  te_upper <- tmax

  vmax_lower <- max(0.5 * observed_vmax, 0.5)
  vmax_upper <- max(vmax_lower + 1e-3, 1)

  rss_combined <- function(params, data, y_col, x_col, vmax_lower, vmax_upper, min_gap) {
    tm1 <- params[1]
    t1 <- params[2]
    t2 <- params[3]
    te <- params[4]
    vmax <- params[5]

    if (!is.finite(tm1) || !is.finite(t1) || !is.finite(t2) || !is.finite(te) || !is.finite(vmax) ||
        t1 <= (tm1 + min_gap) || t2 <= (t1 + min_gap) || te <= (t2 + min_gap) ||
        vmax < vmax_lower || vmax > vmax_upper) {
      return(1e12)
    }

    pred <- tryCatch(
      combined_equation(data[[x_col]], tm1 = tm1, t1 = t1, t2 = t2, te = te, vmax = vmax),
      error = function(e) rep(NA_real_, nrow(data))
    )

    if (any(!is.finite(pred))) {
      return(1e12)
    }

    rss <- sum((data[[y_col]] - pred)^2, na.rm = TRUE)
    if (!is.finite(rss)) 1e12 else rss
  }

  de_result <- DEoptim(
    fn = rss_combined,
    lower = c(tm1_lower, t1_lower, t2_lower, te_lower, vmax_lower),
    upper = c(tm1_upper, t1_upper, t2_upper, te_upper, vmax_upper),
    control = DEoptim.control(itermax = 400, F = 0.8, CR = 0.9, trace = FALSE),
    data = dat,
    y_col = y_col,
    x_col = x_col,
    vmax_lower = vmax_lower,
    vmax_upper = vmax_upper,
    min_gap = min_gap
  )

  list(
    tm1 = unname(de_result$optim$bestmem[1]),
    t1 = unname(de_result$optim$bestmem[2]),
    t2 = unname(de_result$optim$bestmem[3]),
    te = unname(de_result$optim$bestmem[4]),
    vmax = unname(de_result$optim$bestmem[5])
  )
}

## Combined-equation fit for season 2024 
fit_combined <- function(data) {
  # data <- cc_td_season[Season == "2024" & Treatment == "H0"]
  vmax_obs <- max(data$CC, na.rm = TRUE)

  start_vals <- optimize_start_values_combined(
    clean_data = data,
    global_vmax = vmax_obs,
    y_col = "CC",
    x_col = "td_sum"
  )
  if (is.null(start_vals)) {
    return(NULL)
  }

  x_min <- min(data$td_sum, na.rm = TRUE)
  x_max <- max(data$td_sum, na.rm = TRUE)
  x_range <- x_max - x_min
  min_gap <- max(0.5, x_range * 0.03)

  lower <- c(
    tm1 = max(0, x_min),
    t1 = max(0, x_min) + min_gap,
    t2 = max(0, x_min) + (2 * min_gap),
    te = max(0, x_min) + (3 * min_gap),
    vmax = max(vmax_obs * 0.5, 0.05)
  )
  upper <- c(
    tm1 = max(lower[["tm1"]] + min_gap, x_min + x_range * 0.35),
    t1 = max(lower[["t1"]] + min_gap, x_min + x_range * 0.6),
    t2 = max(lower[["t2"]] + min_gap, x_min + x_range * 0.9),
    te = max(lower[["te"]] + min_gap, x_max + x_range * 0.25),
    vmax = max(lower[["vmax"]] + 1e-6, vmax_obs * 1.2)
  )

  tryCatch(
    minpack.lm::nlsLM(
      CC ~ combined_equation(td_sum, tm1, t1, t2, te, vmax),
      data = data,
      start = start_vals,
      lower = lower,
      upper = upper,
      control = minpack.lm::nls.lm.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
}
# BASELINE APPROACH (Standard starting values):
#   1. get_start_values() - Calculate starting parameters per plot (both eq.)
#   2. fit_logistic_curve() - Fit logistic curve to plot data
#   4. fit_gompertz_curve() - Fit Gompertz curve to plot data
# 
# OPTIMIZED APPROACH (DEoptim-based starting values):
#   5. optimize_start_values_logistic() - Find optimal starting values using DEoptim
#   6. optimize_start_values_gompertz() - Find optimal starting values using DEoptim
# ==============================================================================
# Fit the combined canopy model per plot to get robust starts for nlme.
fit_combined <- function(data, y_col = "CC", x_col = "cumtd_DAE", te_cap = 58, start_vals = NULL) {
  vmax_obs <- max(data[[y_col]], na.rm = TRUE)
  x_min <- min(data[[x_col]], na.rm = TRUE)
  x_max <- max(data[[x_col]], na.rm = TRUE)
  x_range <- x_max - x_min
  min_gap <- max(0.5, x_range * 0.03)

  tm1_lower <- max(0, x_min)
  tm1_upper <- max(tm1_lower + min_gap, x_min + x_range * 0.35)

  t1_lower <- tm1_lower + min_gap
  t1_upper <- max(t1_lower + min_gap, x_min + x_range * 0.6)

  t2_lower <- t1_lower + min_gap
  t2_upper <- max(t2_lower + min_gap, x_min + x_range * 0.9)

  te_lower <- t2_lower + min_gap
  te_upper <- te_cap

  vmax_lower <- max(0.5 * vmax_obs, 0.5)
  vmax_upper <- max(vmax_lower + 1e-3, 1)

  model_formula <- stats::as.formula(
    sprintf("%s ~ combined_equation(%s, tm1, t1, t2, te, vmax)", y_col, x_col)
  )

  tryCatch(
    minpack.lm::nlsLM(
      formula = model_formula,
      data = data,
      start = start_vals,
      lower = c(tm1_lower, t1_lower, t2_lower, te_lower, vmax_lower),
      upper = c(tm1_upper, t1_upper, t2_upper, te_upper, vmax_upper),
      control = minpack.lm::nls.lm.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
}
#' Combined Canopy Cover Equation with Plateau
#'
#' Calculates canopy cover dynamics using a logistic growth function with
#' a plateau phase followed by senescence. Suitable for potato canopy modeling.
#'
#' @param td_sum numeric vector - cumulative thermal dose (accumulated over season)
#' @param tm1 numeric - thermal dose at emergence threshold (start of canopy expansion)
#' @param t1 numeric - thermal dose at end of canopy expansion (onset of senescence)
#' @param t2 numeric - thermal dose at end of plateau (beginning of senescence phase)
#' @param te numeric - thermal dose at crop termination
#' @param vmax numeric - maximum canopy cover (0-1 or 0-100% depending on scale)
#'
#' @return numeric vector - canopy cover values following the combined equation
#'
#' @details
#' The equation has three phases:
#'   1. Exponential growth (tm1 to t1): logistic growth phase
#'   2. Plateau (t1 to t2): maximum canopy maintained
#'   3. Senescence decay (t2 to te): gradual surface area decline
#'
#' @examples
#' td_sequence <- seq(0, 800, by=50)
#' cc <- combined_equation(td_sequence, tm1=50, t1=200, t2=400, te=600, vmax=1.0)
#'
#' @keywords internal
combined_equation <- function(td_sum, tm1, t1, t2, te, vmax) {
  ifelse(
    td_sum < t1,
    # Growth phase
    vmax * (1 + (t1 - td_sum) / (t1 - tm1)) * (td_sum / t1) ^ (t1 / (t1 - tm1)),
    ifelse(
      td_sum <= t2,
      # Plateau phase
      vmax,
      # Senescence phase
      vmax * ((te - td_sum) / (te - t2)) * ((td_sum + t1 - t2) / t1) ^ (t1 / (te - t2))
    )
  )
}
# calculate r2 for each fit
r2_nls <- function(fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }

  y <- as.numeric(fit$m$lhs())
  yhat <- as.numeric(stats::fitted(fit))

  sst <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  sse <- sum((y - yhat)^2, na.rm = TRUE)

  if (!is.finite(sst) || sst == 0) {
    return(NA_real_)
  }

  1 - sse / sst
}
extract_combined_params <- function(fit) {
  if (is.null(fit)) {
    return(data.table(
      tm1 = NA_real_,
      t1 = NA_real_,
      t2 = NA_real_,
      te = NA_real_,
      vmax = NA_real_
    ))
  }
  
  est <- coef(fit)
  
  data.table(
    tm1 = unname(est[["tm1"]]),
    t1 = unname(est[["t1"]]),
    t2 = unname(est[["t2"]]),
    te = unname(est[["te"]]),
    vmax = unname(est[["vmax"]])
  )
}
# Build a daily prediction grid spanning the observed range for each plot.
generate_time_grid <- function(data, end_time = 60, step = 0.01) {
  time_grid <- data[, .(cumtd_DAE = seq(0, end_time, step)),
                    
                    by = .(Season, Treatment, Block, PlotID)]
  time_grid[, `:=`(
    Treatment = factor(Treatment, levels = levels(Treatment)),
    Block = factor(Block, levels = levels(Block)),
    PlotID = factor(PlotID, levels = levels(PlotID)))]
  
  setorder(time_grid, Season, Block, PlotID, cumtd_DAE)
  time_grid
}

significance_label <- function(p_value) {
  ifelse(
    is.na(p_value),
    NA_character_,
    ifelse(p_value < 0.001, "**", ifelse(p_value < 0.05, "*", "ns"))
  )
}
summarise_pairwise_table <- function(emm, season, model_name, parameter) {
  out <- as.data.table(summary(pairs(emm, adjust = "tukey"), infer = c(TRUE, TRUE)))
  out[, `:=`(
    Season = season,
    Model = model_name,
    Parameter = parameter,
    Significant = significance_label(p.value)
  )]
  out
}

summarise_emmeans_table <- function(emm, season, model_name, parameter, output_component) {
  out <- as.data.table(summary(emm, infer = c(TRUE, TRUE)))
  
  if ("emmean" %in% names(out)) {
    setnames(out, "emmean", "estimate")
  }
  
  out[, `:=`(
    Season = season,
    Model = model_name,
    Output_component = output_component,
    Parameter = parameter
  )]
  
  if ("p.value" %in% names(out)) {
    out[, Significant := significance_label(p.value)]
  } else {
    out[, Significant := NA_character_]
  }
  
  out
}
# Helper function to extract emmeans results
extract_emmeans_results <- function(model, season_label) {
  # Get treatment means with SE and CI
  emm <- emmeans(model, specs = ~ Treatment)
  means_table <- as.data.table(summary(emm, infer = c(TRUE, TRUE)))
  means_table[, Season := season_label]
  means_table[, Mean_SE := sprintf("%.0f ± %.1f", emmean, SE)]
  # Get pairwise comparisons
  pairs_table <- as.data.table(summary(pairs(emm, adjust = "tukey"), infer = c(TRUE, TRUE)))
  pairs_table[, Season := season_label]
  
  list(means = means_table, pairwise = pairs_table)
}
