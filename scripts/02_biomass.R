source("scripts/00_setup.R")
dma_long <- fread("data/drymatter_organs_long.csv")
#### leaf, stem, tuber, root, stolon in one graph ------------------------------
# ignore stolon - already in stem
dma_organs <- dma_long[Organ %in% c("leaf", "stem", "tuber", "dead_canopy", "leaf_stem")] 
# === 1.1 Summary statistics by organ, treatment, and DAP ===
dma_summary <- calculate_summary_stats(dma_organs[!is.na(value)], 
                                       "value", c("Season", "Treatment", "Organ", "DAP"))
print(dma_summary)

# break HE (2024) line between DAP 53 and 77 by inserting NA
dma_summary_plot <- data.table::rbindlist(list(
  dma_summary,
  data.table(Season = 2024, Treatment = "HE",
             DAP = 63, Organ = c("leaf", "stem","tuber","leaf_stem"),
             mean = NA_real_, se = NA_real_ )), use.names = TRUE, fill = TRUE)


treatments_DAP_arrows[, Organ := "Leaf"]
dma_summary_plot[, Organ := sub("^(.)", "\\U\\1", Organ, perl = TRUE)]
# order the organs 
dma_summary_plot[Organ == "Leaf_stem", Organ := "Leaf:Stem"]
dma_summary_plot[, Organ := factor(
  Organ,
  levels = c("Leaf", "Stem", "Leaf:Stem", "Tuber", "Dead_canopy")
)]
# output 
# season-DAP combinations to keep
dap_keep <- data.table(
  Season = c(2024L, 2024L, 2025L, 2025L),
  DAP    = c(42L,   63L,   56L,   83L)
)
dma_summary_plot[dap_keep, on = .(Season, DAP)]
# assemble the labels for annotate the fig
# all signficance in 2025 - leaf 56 and 83; stem 56; tuber 83
# Create dataframe for asterisk annotations on 2025 organ figure


significance_annotations_25 <- data.table(
  DAP = c(56, 83, 83),
  Organ = c("Stem","Leaf", "Tuber"),
  Season = 2025,
  x_pos = c(56, 83, 83),
  y_pos = c(150, 220, 1300), # adjust y positions as needed
  label = "*"
)
significance_annotations_25[, Organ := factor(
  Organ,
  levels = c("Leaf", "Stem",  "Tuber")
)]
# assemble the final organ figure with annotations an-----------------
panel_tags <- data.table(
  Organ = factor(
    c("Leaf", "Leaf", "Stem", "Stem", "Tuber", "Tuber"),
    levels = c("Leaf", "Stem", "Tuber")
  ),
  Season = c(2024, 2025, 2024, 2025, 2024, 2025),
  tag = c("a", "b", "c", "d", "e", "f")
)

# visualise summary statistics with enhanced styling ---------------
p_organs <- dma_summary_plot[Organ != "Leaf:Stem" & Organ != "Dead_canopy"] |>
  ggplot(aes(DAP, mean, linetype = Treatment, shape = Treatment)) +
  geom_errorbar(aes(x = DAP, y = mean, ymin = mean - se, ymax = mean + se), 
                width = ps, inherit.aes = FALSE) +
  geom_rect(
    data = treatments_DAP_arrows[, .(.N, xmin = min(DAP), xmax = max(DAP)), by = .(Season, Treatment)],
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Treatment),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colors_temp, guide = "none") +
  geom_line(aes(color = Treatment), linewidth = lw -.3) +
  geom_point(size = ps, alpha = 0.7 ) +
  scale_color_manual(values = colors_temp,
                     limits = legend_order) +
  scale_linetype_manual(values = linetype_temp,
                        limits = legend_order) +
  scale_shape_manual(values = point_shape,
                     limits = legend_order) +
  facet_grid(Organ ~ Season, scales = "free_y") +
  scale_x_continuous(name = "Days after planting",
                     expand = c(0, 0), limits = c(0, 120)) +
  scale_y_continuous(name = expression(Dry~weight~(g~m^{-2})))+
  theme_bw(base_size = fontsize, base_family = "Times New Roman") +
  theme(
    legend.position = "top",
    panel.spacing.y = unit(0, "mm"),
    legend.key.width = grid::unit(15, "mm"),
    legend.title = element_text(margin = margin(r = 15)),
    panel.grid = element_blank(),
    panel.spacing.x = unit(7 , "mm")
  )
p_organs_anno <- p_organs +
  geom_text(
    data = panel_tags,
    aes(x = -Inf, y = Inf, label = tag),
    hjust = -0.5,
    vjust = 1.3,
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_text(data = significance_annotations_25, 
            aes(x = x_pos, y = y_pos, label = label),
            size = 8, color = "black", inherit.aes = FALSE)+
  ggh4x::facetted_pos_scales(
    y = list(
      Organ == "Leaf" ~ scale_y_continuous(limits = c(0, 250)),
      Organ == "Stem" ~ scale_y_continuous(limits = c(0, 250)),
      # Organ == "Leaf:Stem" ~ scale_y_continuous(limits = c(0, 3.5)),
      Organ == "Tuber" ~ scale_y_continuous(limits = c(0, 1500))
    )
  ) + ylab(expression(Dry~weight~(g~m^{-2})))



save_plot(p_organs_anno, "fig3_all_organs", height = 6)
