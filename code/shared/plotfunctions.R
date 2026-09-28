# Shared plotting and annotation helpers.
#
# Source this file from analysis scripts with:
#   source("shared/plotfunctions.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin",
  ggplot2::GeomViolin,
  draw_group = function(data, panel_scales, coord, draw_quantiles = NULL) {
    data <- transform(
      data,
      xminv = x - violinwidth * (x - xmin),
      xmaxv = x + violinwidth * (xmax - x)
    )
    group <- data[1, "group"]
    newdata <- plyr::arrange(
      transform(
        data,
        x = if (group %% 2 == 1) xminv else xmaxv
      ),
      if (group %% 2 == 1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

    if (length(draw_quantiles) > 0) {
      stop("draw_quantiles is not implemented for split violins.", call. = FALSE)
    }

    ggplot2::GeomPolygon$draw_panel(newdata, panel_scales, coord)
  }
)

geom_split_violin <- function(
  mapping = NULL,
  data = NULL,
  stat = "ydensity",
  position = "identity",
  ...,
  draw_quantiles = NULL,
  trim = TRUE,
  scale = "area",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      ...
    )
  )
}

get_marker <- function(marker_vecs, marker, cells) {
  marker_values <- marker_vecs[[marker]]
  if (is.null(marker_values)) {
    return(stats::setNames(rep(FALSE, length(cells)), cells))
  }
  out <- marker_values[cells]
  out[is.na(out)] <- FALSE
  out
}

sample_mean_violin <- function(
  obj,
  feature,
  group_col,
  sample_col = "sample",
  assay = "RNA",
  layer = "data",
  group_levels = NULL,
  fill_cols = NULL,
  dot_cols = fill_cols,
  test_method = "wilcox.test"
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("sample_mean_violin() requires the Seurat package.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("sample_mean_violin() requires the dplyr package.", call. = FALSE)
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("sample_mean_violin() requires the ggpubr package.", call. = FALSE)
  }

  df <- Seurat::FetchData(
    obj,
    vars = c(feature, group_col, sample_col),
    assay = assay,
    layer = layer
  )
  colnames(df) <- c("value", "group", "sample")
  df <- dplyr::filter(df, !is.na(group), !is.na(sample))
  if (!is.null(group_levels)) {
    df$group <- factor(df$group, levels = group_levels)
  }

  df_sample <- df |>
    dplyr::group_by(sample, group) |>
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = 0.9) +
    ggplot2::geom_point(
      data = df_sample,
      ggplot2::aes(x = group, y = mean_value, color = group),
      position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 42),
      size = 3,
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    ggplot2::stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.4,
      color = "black",
      linewidth = 0.8
    ) +
    ggpubr::stat_compare_means(
      method = test_method,
      label = "p.signif",
      label.x = 1.5,
      label.y = max(df$value, na.rm = TRUE) * 0.98,
      size = 6
    ) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    ggplot2::labs(title = feature, y = "Expression")

  if (!is.null(fill_cols)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_cols)
  }
  if (!is.null(dot_cols)) {
    p <- p + ggplot2::scale_color_manual(values = dot_cols)
  }

  p
}
