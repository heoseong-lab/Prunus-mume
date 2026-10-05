# =============================================================================
# 06_figures.R  Manuscript figures (PNG 400 dpi and PDF)
#   Figure 1  between-year reproducibility of each descriptor
#   Figure 2  PCoA of the weighted distance, colored by cluster, shaped by flower type
#   Figure 3  dendrogram of the reference partition
#   Figure 4  leave-one-year-out transfer of the distances for the three weighting schemes
# =============================================================================
# Inside-panel legend placement: ggplot2 >= 3.5 uses legend.position.inside; older versions take numeric legend.position
GG35 <- packageVersion("ggplot2") >= "3.5.0"
legend_inside <- function(x, y, hx, hy) {
  if (GG35) theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification.inside = c(hx, hy))
  else theme(legend.position = c(x, y), legend.justification = c(hx, hy))
}
ci <- read_csv(file.path(OUT, "tables", "trait_bootstrap_intervals.csv"), show_col_types = FALSE)
rel <- main$reliability %>% filter(!fixed) %>% left_join(ci[, c("trait", "coef_low", "coef_high")], by = "trait") %>%
  mutate(label = LABEL[trait], panel = ifelse(type == "quantitative", "Quantitative\nvariables\n(Spearman ρ)", "Categorical\nvariables\n(Fleiss' κ)"),
         panel = factor(panel, levels = c("Categorical\nvariables\n(Fleiss' κ)", "Quantitative\nvariables\n(Spearman ρ)")))
rel$label <- factor(rel$label, levels = rel$label[order(rel$panel, rel$coefficient)])
f1 <- ggplot(rel, aes(x = coefficient, y = label, color = panel)) +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_errorbarh(aes(xmin = coef_low, xmax = coef_high), height = 0.25, color = "grey45") +
  geom_point(size = 3) + facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("#0072B2", "#D55E00"), guide = "none") +
  labs(x = "Between-year reproducibility", y = NULL) +
  theme(strip.text.y = element_text(size = 9, angle = 0))
save_plot(f1, "Figure_1_reproducibility", 7.3, 6.0)

pc_tab$cluster <- factor(pc_tab$cluster); pc_tab$flower_type <- factor(pc_tab$flower_type, levels = c("single", "double"))
lab <- sprintf("%s (n = %d)", levels(pc_tab$cluster), as.integer(table(pc_tab$cluster)))
f2 <- ggplot(pc_tab, aes(PCo1, PCo2, color = cluster, shape = flower_type)) +
  geom_point(size = 2.3, alpha = 0.9) +
  scale_color_manual(values = CL_COL[seq_len(K)], labels = lab, name = "Cluster") +
  scale_shape_manual(values = c(single = 16, double = 17), name = "Flower type") +
  labs(x = sprintf("PCoA 1 (%.1f%%)", pc_pct[1]), y = sprintf("PCoA 2 (%.1f%%)", pc_pct[2])) +
  coord_equal() +
  guides(
    color = guide_legend(
      order = 1, ncol = 1, position = "inside",
      theme = theme(legend.position.inside = c(0.02, 0.98), legend.justification.inside = c(0, 1))
    ),
    shape = guide_legend(
      order = 2, position = "inside",
      theme = theme(legend.position.inside = c(0.98, 0.98), legend.justification.inside = c(1, 1))
    )
  ) +
  theme(
    legend.spacing.y = unit(0, "pt"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "transparent", color = NA)
  )
if (!GG35) {
  # ggplot2 < 3.5 cannot place two legends in different corners: keep the cluster legend inside at the top-left and
  # add the flower-type legend as a separate grob anchored to the top-right corner of the panel
  f2_tl <- f2 + guides(shape = "none") + legend_inside(0.02, 0.98, 0, 1)
  f2_sh <- ggplot(pc_tab, aes(PCo1, PCo2, shape = flower_type)) + geom_point(size = 2.3) +
    scale_shape_manual(values = c(single = 16, double = 17), name = "Flower type") +
    guides(shape = guide_legend(title.position = "top", ncol = 1)) +
    theme(legend.position = "right", legend.direction = "vertical", legend.background = element_blank(), legend.key = element_rect(fill = "transparent", color = NA))
  gs <- ggplotGrob(f2_sh); leg <- gs$grobs[[which(gs$layout$name == "guide-box")]]
  g2 <- ggplotGrob(f2_tl); pl <- g2$layout[g2$layout$name == "panel", ]
  leg_vp <- grid::viewport(x = unit(1, "npc") - unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"), just = c("right", "top"),
                           width = sum(leg$widths), height = sum(leg$heights))
  f2 <- gtable::gtable_add_grob(g2, grid::grobTree(leg, vp = leg_vp), t = pl$t, l = pl$l, b = pl$b, r = pl$r, name = "flower-type-legend", clip = "off")
}
save_plot(f2, "Figure_2_PCoA", 7.3, 6.8)

hc <- base_run$hc; dend <- as.dendrogram(hc); ord <- order.dendrogram(dend)
seg <- function(dd, x0 = 0) {   # recursive extraction of dendrogram segments with cluster colors
  out <- list()
  walk <- function(node, xpos) {
    h <- attr(node, "height")
    if (is.leaf(node)) return(list(x = xpos, h = 0, cl = base_run$cl[attr(node, "label")]))
    kids <- lapply(seq_along(node), function(i) NULL)
    xs <- c(); cls <- c(); cur <- xpos
    for (i in seq_along(node)) {
      nleaf <- attr(node[[i]], "members"); r <- walk(node[[i]], cur); xs <- c(xs, r$x); cls <- c(cls, r$cl)
      out[[length(out) + 1]] <<- data.frame(x = r$x, xend = r$x, y = r$h, yend = h, cl = r$cl)
      cur <- cur + nleaf
    }
    cl <- if (length(unique(cls)) == 1) cls[1] else NA
    out[[length(out) + 1]] <<- data.frame(x = min(xs), xend = max(xs), y = h, yend = h, cl = cl)
    list(x = mean(range(xs)), h = h, cl = cl)
  }
  walk(dd, 1); bind_rows(out)
}
segs <- seg(dend); segs$cl <- factor(segs$cl, levels = seq_len(K))
f3 <- ggplot(segs) + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = cl), linewidth = 0.35) +
  scale_color_manual(values = CL_COL[seq_len(K)], na.value = "grey35", name = "Cluster", breaks = levels(segs$cl)) +
  labs(x = NULL, y = "Ward merge height") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank()) +
  guides(color = guide_legend(ncol = 1, position = "inside")) +
  legend_inside(0.98, 0.98, 1, 1) +
  theme(
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", color = NA)
  )
save_plot(f3, "Figure_3_dendrogram", 7.3, 4.8)

# Figure 4: between-year transfer of the distances (leave-one-year-out)
SCHEME_LAB <- c(equal = "Equal weights (15 measured)", retained = "Equal weights (retained)", weighted = "Reproducibility weights")
SCHEME_COL <- c(equal = "#E69F00", retained = "#56B4E9", weighted = "#0072B2")
lo <- temporal %>% mutate(scheme = factor(scheme, levels = names(SCHEME_LAB), labels = SCHEME_LAB), held_out_year = factor(held_out_year))
f4 <- ggplot(lo, aes(held_out_year, distance_spearman, fill = scheme)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", distance_spearman)), position = position_dodge(width = 0.78), vjust = -0.35, size = 2.9) +
  scale_fill_manual(values = unname(SCHEME_COL), name = NULL) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Held-out year (LOYO)", y = "Spearman \u03c1 (training vs. held-out distances)") +
  legend_inside(0.02, 0.98, 0, 1) +
  theme(
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(10, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.grid.major.x = element_blank()
  ) + guides(fill = guide_legend(nrow = 1, byrow = TRUE))
save_plot(f4, "Figure_4_between_year_transfer", 6.2, 4.2)
ch <- data.frame(iteration = seq_along(sc$coverage_history), coverage = 100 * sc$coverage_history)
f5 <- ggplot(ch, aes(iteration, coverage)) +
  geom_hline(yintercept = 100 * c(CORE_COV_MIN, CORE_COV_MAX), linetype = c("dashed", "solid"), color = "grey40") +
  geom_vline(xintercept = sc$N0 + 0.5, linetype = "dotted", color = "grey40") +
  geom_line(color = "#0072B2") + geom_point(color = "#0072B2", size = 1.8) +
  annotate("text", x = sc$N0 + 0.3, y = 20, label = sprintf("N0 = %d (CV > %.0f%%)", sc$N0, 100 * CORE_COV_MIN), hjust = 1, size = 3.2) +
  annotate("text", x = 1, y = 100 * CORE_COV_MAX - 5, label = sprintf("N1 = %d (CV \u2265 %.0f%%)", sc$N1, 100 * CORE_COV_MAX), hjust = 0, size = 3.2) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) + scale_x_continuous(breaks = seq(0, length(sc$history), 5)) +
  labs(x = "Accessions added in the covering phase", y = "Coverage of descriptor classes (%)")
save_plot(f5, "Figure_5_core_coverage", 6.0, 4.2)
message("Figures written")
