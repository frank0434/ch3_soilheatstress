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
