##########################################################
# CODEX immune marker overlap audit
##########################################################
#
# Purpose:
#   Count how often immune markers overlap before deciding whether the
#   priority order in immcelltype assignment is acceptable.
#
# This script does not modify Seurat metadata.
#
# Example usage after loading either obj_codex or seurat_list_codex:
#
#   source("CODEX/audits/check_immune_marker_overlaps.R")
#   audit <- audit_codex_immune_overlaps(obj = obj_codex)
#   audit$summary
#   audit$lineage_combo_counts
#   audit$marker_pair_counts
#
#   audit <- audit_codex_immune_overlaps(
#     seurat_list = seurat_list_codex,
#     output_dir = "../results/audits/codex_immune_overlap_audit"
#   )

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

get_marker_positive_table <- function(
  obj,
  markers = c("CD45", "CD4", "CD8", "CD14", "CD68", "FOXP3", "CXCR3"),
  assay = "RNA",
  layer = "counts",
  sample = NA_character_
) {
  DefaultAssay(obj) <- assay
  counts <- tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )

  cells <- colnames(obj)
  out <- data.frame(
    sample = sample,
    cell = cells,
    stringsAsFactors = FALSE
  )

  for (marker in markers) {
    out[[marker]] <- if (marker %in% rownames(counts)) {
      as.numeric(counts[marker, cells]) > 0
    } else {
      rep(FALSE, length(cells))
    }
  }

  out
}

make_combo_label <- function(df, markers) {
  apply(df[, markers, drop = FALSE], 1, function(row) {
    positive <- markers[as.logical(row)]
    if (length(positive) == 0) "none" else paste(positive, collapse = ";")
  })
}

get_current_codex_immune_assignment <- function(marker_table) {
  rule_hits <- lapply(seq_len(nrow(marker_table)), function(i) {
    hits <- character()
    if (isTRUE(marker_table$CD4[i]) && isTRUE(marker_table$FOXP3[i])) hits <- c(hits, "CD4+ FOXP3+")
    if (isTRUE(marker_table$CD4[i]) && isTRUE(marker_table$CXCR3[i])) hits <- c(hits, "CD4+ CXCR3+")
    if (isTRUE(marker_table$CD14[i])) hits <- c(hits, "CD14+")
    if (isTRUE(marker_table$CD68[i])) hits <- c(hits, "CD68+")
    if (isTRUE(marker_table$CD8[i])) hits <- c(hits, "CD8+")
    if (isTRUE(marker_table$CD4[i])) hits <- c(hits, "CD4+")
    if (isTRUE(marker_table$CD45[i])) hits <- c(hits, "CD45+")
    hits
  })

  assigned <- vapply(rule_hits, function(hits) {
    if (length(hits) == 0) "none" else hits[1]
  }, character(1))

  hit_labels <- vapply(rule_hits, function(hits) {
    if (length(hits) == 0) "none" else paste(hits, collapse = ";")
  }, character(1))

  suppressed_labels <- vapply(rule_hits, function(hits) {
    if (length(hits) <= 1) "none" else paste(hits[-1], collapse = ";")
  }, character(1))

  data.frame(
    assigned_by_current_priority = assigned,
    assignment_rule_hits = hit_labels,
    suppressed_by_current_priority = suppressed_labels,
    n_assignment_rule_hits = lengths(rule_hits),
    assignment_affecting_overlap = lengths(rule_hits) > 1,
    stringsAsFactors = FALSE
  )
}

plot_assignment_affecting_overlaps <- function(overlap_counts) {
  if (nrow(overlap_counts) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No assignment-affecting immune marker overlaps found.") +
        theme_void()
    )
  }

  overlap_counts$marker_combo <- factor(
    overlap_counts$marker_combo,
    levels = rev(unique(overlap_counts$marker_combo[order(overlap_counts$n_cells)]))
  )

  ggplot(
    overlap_counts,
    aes(x = marker_combo, y = n_cells, fill = assigned_by_current_priority)
  ) +
    geom_col(width = 0.75, color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Cells",
      fill = "Assigned to",
      title = "CODEX immune marker overlaps that affect assignment"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}

plot_suppressed_assignment_rules <- function(suppressed_counts) {
  if (nrow(suppressed_counts) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No suppressed immune assignment rules found.") +
        theme_void()
    )
  }

  suppressed_counts$marker_combo <- factor(
    suppressed_counts$marker_combo,
    levels = rev(unique(suppressed_counts$marker_combo[order(suppressed_counts$n_cells)]))
  )

  ggplot(
    suppressed_counts,
    aes(x = marker_combo, y = n_cells, fill = suppressed_assignment_rule)
  ) +
    geom_col(width = 0.75, color = "black", linewidth = 0.2) +
    coord_flip() +
    facet_wrap(~ assigned_by_current_priority, scales = "free_y") +
    labs(
      x = NULL,
      y = "Cells",
      fill = "Also matched, but not assigned",
      title = "Immune rules hidden by current priority order"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}

audit_codex_immune_overlaps <- function(
  obj = NULL,
  seurat_list = NULL,
  assay = "RNA",
  layer = "counts",
  lineage_markers = c("CD4", "CD8", "CD14", "CD68"),
  context_markers = c("CD45", "FOXP3", "CXCR3"),
  output_dir = NULL
) {
  markers <- unique(c(context_markers, lineage_markers))

  if (!is.null(seurat_list)) {
    marker_table <- do.call(rbind, lapply(names(seurat_list), function(sample_name) {
      get_marker_positive_table(
        obj = seurat_list[[sample_name]],
        markers = markers,
        assay = assay,
        layer = layer,
        sample = sample_name
      )
    }))
  } else if (!is.null(obj)) {
    marker_table <- get_marker_positive_table(
      obj = obj,
      markers = markers,
      assay = assay,
      layer = layer,
      sample = "merged"
    )
  } else {
    stop("Provide either obj or seurat_list.", call. = FALSE)
  }

  marker_table$lineage_n_positive <- rowSums(marker_table[, lineage_markers, drop = FALSE])
  marker_table$lineage_combo <- make_combo_label(marker_table, lineage_markers)
  marker_table$all_marker_combo <- make_combo_label(marker_table, markers)
  marker_table$has_lineage_overlap <- marker_table$lineage_n_positive > 1
  marker_table <- cbind(marker_table, get_current_codex_immune_assignment(marker_table))

  summary <- data.frame(
    n_cells = nrow(marker_table),
    n_any_lineage_positive = sum(marker_table$lineage_n_positive > 0),
    n_lineage_overlap = sum(marker_table$has_lineage_overlap),
    pct_lineage_overlap_all_cells = 100 * mean(marker_table$has_lineage_overlap),
    pct_lineage_overlap_lineage_positive = 100 * (
      sum(marker_table$has_lineage_overlap) /
        max(1, sum(marker_table$lineage_n_positive > 0))
    )
  )

  lineage_combo_counts <- as.data.frame(table(marker_table$lineage_combo), stringsAsFactors = FALSE)
  colnames(lineage_combo_counts) <- c("lineage_combo", "n_cells")
  lineage_combo_counts <- lineage_combo_counts[order(-lineage_combo_counts$n_cells), ]

  all_marker_combo_counts <- as.data.frame(table(marker_table$all_marker_combo), stringsAsFactors = FALSE)
  colnames(all_marker_combo_counts) <- c("all_marker_combo", "n_cells")
  all_marker_combo_counts <- all_marker_combo_counts[order(-all_marker_combo_counts$n_cells), ]

  marker_pair_counts <- do.call(rbind, combn(lineage_markers, 2, simplify = FALSE, FUN = function(pair) {
    data.frame(
      marker_1 = pair[1],
      marker_2 = pair[2],
      n_cells = sum(marker_table[[pair[1]]] & marker_table[[pair[2]]]),
      stringsAsFactors = FALSE
    )
  }))
  marker_pair_counts <- marker_pair_counts[order(-marker_pair_counts$n_cells), ]

  assignment_affecting <- marker_table[marker_table$assignment_affecting_overlap, ]
  if (nrow(assignment_affecting) > 0) {
    assignment_affecting_overlap_counts <- aggregate(
      cell ~ all_marker_combo + assigned_by_current_priority + assignment_rule_hits,
      data = assignment_affecting,
      FUN = length
    )
    colnames(assignment_affecting_overlap_counts) <- c(
      "marker_combo",
      "assigned_by_current_priority",
      "assignment_rule_hits",
      "n_cells"
    )
    assignment_affecting_overlap_counts <- assignment_affecting_overlap_counts[
      order(
        -assignment_affecting_overlap_counts$n_cells,
        assignment_affecting_overlap_counts$assigned_by_current_priority,
        assignment_affecting_overlap_counts$marker_combo
      ),
    ]
  } else {
    assignment_affecting_overlap_counts <- data.frame(
      marker_combo = character(),
      assigned_by_current_priority = character(),
      assignment_rule_hits = character(),
      n_cells = integer()
    )
  }

  if (nrow(assignment_affecting) > 0) {
    suppressed_long <- do.call(rbind, lapply(seq_len(nrow(assignment_affecting)), function(i) {
      suppressed <- strsplit(assignment_affecting$suppressed_by_current_priority[i], ";", fixed = TRUE)[[1]]
      suppressed <- suppressed[suppressed != "none"]
      if (length(suppressed) == 0) return(NULL)
      data.frame(
        marker_combo = assignment_affecting$all_marker_combo[i],
        assigned_by_current_priority = assignment_affecting$assigned_by_current_priority[i],
        suppressed_assignment_rule = suppressed,
        cell = assignment_affecting$cell[i],
        stringsAsFactors = FALSE
      )
    }))
  } else {
    suppressed_long <- NULL
  }

  if (!is.null(suppressed_long) && nrow(suppressed_long) > 0) {
    suppressed_assignment_counts <- aggregate(
      cell ~ marker_combo + assigned_by_current_priority + suppressed_assignment_rule,
      data = suppressed_long,
      FUN = length
    )
    colnames(suppressed_assignment_counts) <- c(
      "marker_combo",
      "assigned_by_current_priority",
      "suppressed_assignment_rule",
      "n_cells"
    )
    suppressed_assignment_counts <- suppressed_assignment_counts[
      order(
        suppressed_assignment_counts$suppressed_assignment_rule,
        -suppressed_assignment_counts$n_cells,
        suppressed_assignment_counts$marker_combo
      ),
    ]
  } else {
    suppressed_assignment_counts <- data.frame(
      marker_combo = character(),
      assigned_by_current_priority = character(),
      suppressed_assignment_rule = character(),
      n_cells = integer()
    )
  }

  by_sample <- aggregate(
    has_lineage_overlap ~ sample,
    data = marker_table,
    FUN = function(x) c(n_overlap = sum(x), pct_overlap = 100 * mean(x))
  )
  by_sample <- data.frame(
    sample = by_sample$sample,
    n_overlap = by_sample$has_lineage_overlap[, "n_overlap"],
    pct_overlap = by_sample$has_lineage_overlap[, "pct_overlap"],
    row.names = NULL
  )

  result <- list(
    summary = summary,
    by_sample = by_sample,
    lineage_combo_counts = lineage_combo_counts,
    all_marker_combo_counts = all_marker_combo_counts,
    marker_pair_counts = marker_pair_counts,
    assignment_affecting_overlap_counts = assignment_affecting_overlap_counts,
    suppressed_assignment_counts = suppressed_assignment_counts,
    per_cell_marker_table = marker_table
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
    write.csv(by_sample, file.path(output_dir, "by_sample.csv"), row.names = FALSE)
    write.csv(lineage_combo_counts, file.path(output_dir, "lineage_combo_counts.csv"), row.names = FALSE)
    write.csv(all_marker_combo_counts, file.path(output_dir, "all_marker_combo_counts.csv"), row.names = FALSE)
    write.csv(marker_pair_counts, file.path(output_dir, "marker_pair_counts.csv"), row.names = FALSE)
    write.csv(
      assignment_affecting_overlap_counts,
      file.path(output_dir, "assignment_affecting_overlap_counts.csv"),
      row.names = FALSE
    )
    write.csv(
      suppressed_assignment_counts,
      file.path(output_dir, "suppressed_assignment_counts.csv"),
      row.names = FALSE
    )
    write.csv(marker_table, file.path(output_dir, "per_cell_marker_table.csv"), row.names = FALSE)

    p <- plot_assignment_affecting_overlaps(assignment_affecting_overlap_counts)
    plot_height <- max(4, 0.28 * nrow(assignment_affecting_overlap_counts) + 2)
    ggsave(
      filename = file.path(output_dir, "assignment_affecting_overlaps_barplot.pdf"),
      plot = p,
      width = 10,
      height = plot_height
    )
    ggsave(
      filename = file.path(output_dir, "assignment_affecting_overlaps_barplot.png"),
      plot = p,
      width = 10,
      height = plot_height,
      dpi = 300
    )

    p_suppressed <- plot_suppressed_assignment_rules(suppressed_assignment_counts)
    suppressed_plot_height <- max(4, 0.28 * nrow(suppressed_assignment_counts) + 2)
    ggsave(
      filename = file.path(output_dir, "suppressed_assignment_rules_barplot.pdf"),
      plot = p_suppressed,
      width = 12,
      height = suppressed_plot_height
    )
    ggsave(
      filename = file.path(output_dir, "suppressed_assignment_rules_barplot.png"),
      plot = p_suppressed,
      width = 12,
      height = suppressed_plot_height,
      dpi = 300
    )

    message("CODEX immune overlap audit saved to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
    message("Files: ", paste(list.files(output_dir), collapse = ", "))
  }

  result
}

if (!exists("seurat_list_codex")) {
  load("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2/data/codex/seurat_list_codex.RData")
}

if (!exists("seurat_list_codex")) {
  stop("Could not find seurat_list_codex. Load it before running this script.", call. = FALSE)
}

audit <- audit_codex_immune_overlaps(
  seurat_list = seurat_list_codex,
  output_dir = "//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2/results/audits/codex_immune_overlap_audit"
)

print(audit$summary)
print(audit$assignment_affecting_overlap_counts)
print(audit$suppressed_assignment_counts)
