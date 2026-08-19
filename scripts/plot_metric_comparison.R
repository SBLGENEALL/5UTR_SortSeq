#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_metric_comparison.R METRIC_COMPARISON_DIR OUTPUT_DIR")
}

required_packages <- c("ggplot2", "scales", "viridisLite", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

supported_path <- file.path(input_dir, "metric_comparison_supported.csv")
summary_path <- file.path(input_dir, "metric_comparison_summary.csv")
overlap_path <- file.path(input_dir, "metric_comparison_topn_overlap.csv")
if (!file.exists(supported_path) || !file.exists(summary_path) || !file.exists(overlap_path)) {
  stop("Metric-comparison CSVs are missing under: ", input_dir)
}

comparison <- read.csv(
  supported_path, stringsAsFactors = FALSE, check.names = FALSE
)
summary_table <- read.csv(
  summary_path, stringsAsFactors = FALSE, check.names = FALSE
)
overlap_table <- read.csv(
  overlap_path, stringsAsFactors = FALSE, check.names = FALSE
)

numeric_columns <- c(
  "expected_bin_score", "high15_probability",
  "score_rank_filtered", "top15_rank_filtered",
  "score_rank_minus_top15_rank", "total_6bin_count",
  "high15_probability", paste0("bin", 1:6, "_probability")
)
for (column in numeric_columns) {
  if (column %in% colnames(comparison)) {
    comparison[[column]] <- suppressWarnings(as.numeric(comparison[[column]]))
  }
}
as_bool <- function(values) {
  tolower(as.character(values)) %in% c("true", "t", "1")
}
if ("is_reference_variant" %in% colnames(comparison)) {
  comparison$is_reference_variant <- as_bool(comparison$is_reference_variant)
} else {
  comparison$is_reference_variant <-
    tolower(trimws(comparison$variant_id)) %in% c("original", "orginal")
}

colors <- c(
  consensus = "#D55E00", score_only = "#0072B2", top15_only = "#009E73",
  neither = "#C7CED3", reference = "#B2182B", ink = "#24323D",
  muted = "#667580", grid = "#DDE5EA"
)
membership_levels <- c("neither", "consensus", "score_only", "top15_only")
comparison$top_list_membership <- factor(
  comparison$top_list_membership, levels = membership_levels
)

theme_comparison <- function() {
  ggplot2::theme_minimal(base_size = 14, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = colors[["ink"]], face = "bold", size = 20,
        margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        color = colors[["muted"]], size = 13, margin = ggplot2::margin(b = 10)
      ),
      plot.caption = ggplot2::element_text(
        color = colors[["muted"]], size = 11, hjust = 0,
        margin = ggplot2::margin(t = 9)
      ),
      axis.title = ggplot2::element_text(
        color = colors[["ink"]], face = "bold", size = 14
      ),
      axis.text = ggplot2::element_text(color = colors[["ink"]], size = 11),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = colors[["grid"]], linewidth = 0.35
      ),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = ggplot2::margin(18, 28, 16, 22)
    )
}

save_png <- function(filename, plot, width = 11.8, height = 6.6) {
  ragg::agg_png(
    filename = file.path(output_dir, filename), width = width, height = height,
    units = "in", res = 320, background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

summary_value <- function(metric, scope = NULL) {
  selected <- summary_table$metric == metric
  if (!is.null(scope)) {
    selected <- selected & summary_table$scope == scope
  }
  values <- suppressWarnings(as.numeric(summary_table$value[selected]))
  if (length(values) == 0) return(NA_real_)
  values[[1]]
}

rho <- summary_value(
  "spearman_expected_score_vs_high15_probability", "robust_read_support"
)
configured_top_n <- summary_value("configured_top_n")
actual_top_n <- summary_value("actual_top_n")
reference_score <- summary_value("reference_score")
reference_top15 <- summary_value("reference_high15_probability")
consensus_n <- summary_value("top_list_consensus_count")
score_only_n <- summary_value("top_list_score_only_count")
top15_only_n <- summary_value("top_list_top15_only_count")

plots <- list()
highlight <- comparison[
  comparison$top_list_membership != "neither", , drop = FALSE
]
reference <- comparison[comparison$is_reference_variant, , drop = FALSE]

# 21: direct metric relationship on the exact same read-supported UTR universe.
p_metric <- ggplot2::ggplot(
  comparison,
  ggplot2::aes(
    x = expected_bin_score,
    y = high15_probability
  )
) +
  ggplot2::geom_point(color = colors[["neither"]], alpha = 0.48, size = 1.5) +
  ggplot2::geom_point(
    data = highlight,
    ggplot2::aes(color = top_list_membership, shape = top_list_membership),
    alpha = 0.88, size = 2.4
  ) +
  ggplot2::scale_color_manual(values = colors[membership_levels], drop = FALSE) +
  ggplot2::scale_shape_manual(
    values = c(neither = 16, consensus = 16, score_only = 17, top15_only = 15),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Supporting six-bin score versus primary High15 probability",
    subtitle = sprintf(
      "Same read-supported UTRs; Spearman rho = %.3f",
      rho
    ),
    caption = paste0(
      "Total six-bin count > 200; unsorted >= 50; bin1+2 raw count >= 20. ",
      "High15 is primary; score summarizes the complete six-bin distribution."
    ),
    x = "Expected bin score (whole six-bin position)",
    y = "High15 probability: P(bin1 or bin2 | UTR, target gate)",
    color = sprintf("Top %.0f membership", actual_top_n),
    shape = sprintf("Top %.0f membership", actual_top_n)
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  theme_comparison()
if (is.finite(reference_score)) {
  p_metric <- p_metric + ggplot2::geom_vline(
    xintercept = reference_score, linetype = "dashed",
    color = colors[["reference"]], linewidth = 0.65
  )
}
if (is.finite(reference_top15)) {
  p_metric <- p_metric + ggplot2::geom_hline(
    yintercept = reference_top15, linetype = "dashed",
    color = colors[["reference"]], linewidth = 0.65
  )
}
if (nrow(reference) == 1) {
  p_metric <- p_metric + ggplot2::geom_point(
    data = reference, shape = 8, size = 5.2, stroke = 1.2,
    color = colors[["reference"]]
  )
}
save_png("21_expected_score_vs_high15_probability.png", p_metric)
plots[["21"]] <- p_metric

# 22: rank agreement makes candidate reordering visible without assuming a
# linear relation between the two differently scaled metrics.
p_rank <- ggplot2::ggplot(
  comparison,
  ggplot2::aes(x = score_rank_filtered, y = top15_rank_filtered)
) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", color = colors[["muted"]]
  ) +
  ggplot2::geom_point(color = colors[["neither"]], alpha = 0.48, size = 1.5) +
  ggplot2::geom_point(
    data = highlight,
    ggplot2::aes(color = top_list_membership, shape = top_list_membership),
    alpha = 0.88, size = 2.4
  ) +
  ggplot2::scale_x_reverse(labels = scales::comma) +
  ggplot2::scale_y_reverse(labels = scales::comma) +
  ggplot2::scale_color_manual(values = colors[membership_levels], drop = FALSE) +
  ggplot2::scale_shape_manual(
    values = c(neither = 16, consensus = 16, score_only = 17, top15_only = 15),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Rank agreement between the two phenotype summaries",
    subtitle = "The upper-right corner is best by both metrics; the diagonal is identical rank",
    x = "Expected-score rank (1 = best)",
    y = "High15-probability rank (1 = best)",
    color = sprintf("Top %.0f membership", actual_top_n),
    shape = sprintf("Top %.0f membership", actual_top_n)
  ) +
  theme_comparison()
save_png("22_score_rank_vs_top15_rank.png", p_rank)
plots[["22"]] <- p_rank

# 23: exact Top-N list overlap.
overlap_counts <- data.frame(
  category = factor(
    c("Score only", "Consensus", "Top15 only"),
    levels = c("Score only", "Consensus", "Top15 only")
  ),
  variants = c(score_only_n, consensus_n, top15_only_n),
  stringsAsFactors = FALSE
)
p_overlap <- ggplot2::ggplot(
  overlap_counts,
  ggplot2::aes(x = category, y = variants, fill = category)
) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(variants)),
    vjust = -0.25, size = 4, color = colors[["ink"]]
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Score only" = colors[["score_only"]],
      "Consensus" = colors[["consensus"]],
      "Top15 only" = colors[["top15_only"]]
    )
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.14))
  ) +
  ggplot2::labs(
    title = sprintf("Overlap of the two Top %.0f candidate lists", actual_top_n),
    subtitle = sprintf(
      "Consensus = %.0f of %.0f (%.1f%%)",
      consensus_n, actual_top_n, 100 * consensus_n / actual_top_n
    ),
    x = NULL, y = "UTRs"
  ) +
  theme_comparison()
save_png("23_top_candidate_overlap.png", p_overlap)
plots[["23"]] <- p_overlap

# 24: only discordant Top-N candidates are shown because these profiles decide
# which endpoint is selecting a biologically different phenotype.
discordant <- comparison[
  comparison$top_list_membership %in% c("score_only", "top15_only"),
  , drop = FALSE
]
if (nrow(discordant) > 0) {
  discordant$relevant_rank <- ifelse(
    discordant$top_list_membership == "score_only",
    discordant$score_rank_filtered,
    discordant$top15_rank_filtered
  )
  discordant <- discordant[
    order(discordant$top_list_membership, discordant$relevant_rank),
    , drop = FALSE
  ]
  keep_index <- unlist(lapply(
    c("score_only", "top15_only"),
    function(group) {
      index <- which(discordant$top_list_membership == group)
      head(index, 20)
    }
  ))
  discordant <- discordant[keep_index, , drop = FALSE]
  discordant$display_id <- paste0(
    discordant$variant_id, " [", discordant$top_list_membership, "]"
  )
  probability_columns <- paste0("bin", 1:6, "_probability")
  heat <- do.call(
    rbind,
    lapply(seq_len(nrow(discordant)), function(row_number) {
      data.frame(
        variant_id = discordant$display_id[[row_number]],
        bin = factor(paste0("bin", 1:6), levels = paste0("bin", 1:6)),
        probability = as.numeric(discordant[row_number, probability_columns]),
        stringsAsFactors = FALSE
      )
    })
  )
  heat$variant_id <- factor(heat$variant_id, levels = rev(discordant$display_id))
  p_discordant <- ggplot2::ggplot(
    heat,
    ggplot2::aes(x = bin, y = variant_id, fill = probability)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradientn(
      colors = viridisLite::viridis(256, option = "B"),
      labels = scales::percent
    ) +
    ggplot2::labs(
      title = "Six-bin profiles of candidates selected by only one metric",
      subtitle = "Up to 20 per group; score-only often favors broad/intermediate shifts, while top15-only favors the high tail",
      x = "FACS bin", y = "UTR", fill = "Estimated\ncell fraction"
    ) +
    theme_comparison() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7.5)
    )
  save_png(
    "24_discordant_candidate_bin_heatmap.png", p_discordant,
    10.5, max(7.2, nrow(discordant) * 0.22)
  )
  plots[["24"]] <- p_discordant
}

plot_statistics <- data.frame(
  metric = c(
    "read_supported_utr", "spearman_score_vs_high15_probability", "configured_top_n",
    "actual_top_n", "top_n_consensus", "top_n_score_only",
    "top_n_top15_only"
  ),
  value = c(
    nrow(comparison), rho, configured_top_n, actual_top_n,
    consensus_n, score_only_n, top15_only_n
  ),
  stringsAsFactors = FALSE
)
write.csv(
  plot_statistics,
  file.path(output_dir, "metric_comparison_plot_statistics.csv"),
  quote = FALSE, row.names = FALSE
)

plots <- plots[order(names(plots))]
grDevices::pdf(
  file.path(output_dir, "metric_comparison_figures.pdf"),
  width = 13.333, height = 7.5, onefile = TRUE,
  family = "Helvetica", paper = "special"
)
for (plot in plots) {
  print(plot)
}
grDevices::dev.off()

message("Metric-comparison figures written to: ", normalizePath(output_dir))
