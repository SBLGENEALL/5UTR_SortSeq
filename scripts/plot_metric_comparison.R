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
  "step6_high15_probability", "step7_expected_bin_score",
  "step8_high15_relative_enrichment",
  "step8_high15_log2_relative_enrichment",
  "step6_rank_filtered", "step7_rank_filtered", "step8_rank_filtered",
  "score_rank_filtered", "top15_rank_filtered",
  "score_rank_minus_top15_rank", "total_6bin_count",
  paste0("bin", 1:6, "_probability"),
  paste0("bin", 1:6, "_relative_enrichment"),
  paste0("bin", 1:6, "_log2_relative_enrichment")
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
configured_min_total <- summary_value("configured_min_total_count")
configured_min_unsorted <- summary_value("configured_min_unsorted_count")
configured_min_high <- summary_value("configured_min_high_bin_count")
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
      "Total six-bin count >= ", configured_min_total,
      "; bin1+2 raw count >= ", configured_min_high,
      ifelse(configured_min_unsorted > 0,
        paste0("; unsorted >= ", configured_min_unsorted, ". "),
        "; unsorted is QC only and is not used as a cutoff. "
      ),
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
        probability = as.numeric(unlist(
          discordant[row_number, probability_columns], use.names = FALSE
        )),
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

# 25: pairwise Step 6/7/8 rank correlations. Step 6 versus Step 8 must be 1
# because Step 8 High15 enrichment is Step 6 divided by the constant 0.15.
correlation_heat <- summary_table[
  summary_table$section == "correlation" &
    summary_table$scope == "robust_read_support" &
    summary_table$method == "spearman" &
    summary_table$comparison %in% c(
      "step6_high15_vs_step7_expected_score",
      "step6_high15_vs_step8_high15_relative",
      "step7_expected_score_vs_step8_high15_relative"
    ),
  , drop = FALSE
]
if (nrow(correlation_heat) > 0) {
  correlation_labels <- c(
    step6_high15_vs_step7_expected_score = "Step 6 vs Step 7",
    step6_high15_vs_step8_high15_relative = "Step 6 vs Step 8",
    step7_expected_score_vs_step8_high15_relative = "Step 7 vs Step 8"
  )
  correlation_heat$pair <- unname(correlation_labels[correlation_heat$comparison])
  correlation_heat$rho <- suppressWarnings(as.numeric(correlation_heat$value))
  correlation_heat$pair <- factor(
    correlation_heat$pair,
    levels = c("Step 6 vs Step 7", "Step 6 vs Step 8", "Step 7 vs Step 8")
  )
  p_correlation <- ggplot2::ggplot(
    correlation_heat,
    ggplot2::aes(x = pair, y = "Spearman", fill = rho)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("rho = %.3f", rho)),
      color = "white", fontface = "bold", size = 5
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1)
    ) +
    ggplot2::labs(
      title = "Rank correlation among Steps 6, 7, and 8",
      subtitle = paste0(
        "Step 6 vs 7 is the informative comparison; Step 6 vs 8 is an identity check"
      ),
      x = NULL, y = NULL, fill = "Spearman\nrho"
    ) +
    theme_comparison() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  save_png("25_step6_step7_step8_correlation.png", p_correlation, 10.5, 5.8)
  plots[["25"]] <- p_correlation
}

# 26: candidate-list overlap at Top20/50/100 for all pairwise comparisons.
if ("pair" %in% colnames(overlap_table)) {
  overlap_plot <- overlap_table[
    overlap_table$pair %in% c("step6_vs_step7", "step6_vs_step8", "step7_vs_step8"),
    , drop = FALSE
  ]
  if (nrow(overlap_plot) > 0) {
    overlap_plot$requested_top_n <- suppressWarnings(
      as.numeric(overlap_plot$requested_top_n)
    )
    overlap_plot$overlap_percent_of_each_list <- suppressWarnings(
      as.numeric(overlap_plot$overlap_percent_of_each_list)
    )
    overlap_plot$pair <- factor(
      overlap_plot$pair,
      levels = c("step6_vs_step7", "step6_vs_step8", "step7_vs_step8"),
      labels = c("Step 6 vs 7", "Step 6 vs 8", "Step 7 vs 8")
    )
    p_pair_overlap <- ggplot2::ggplot(
      overlap_plot,
      ggplot2::aes(
        x = factor(requested_top_n),
        y = overlap_percent_of_each_list,
        fill = pair
      )
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.72) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.0f%%", overlap_percent_of_each_list)),
        position = ggplot2::position_dodge(width = 0.78),
        vjust = -0.25, size = 3.6
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Step 6 vs 7" = "#D55E00",
          "Step 6 vs 8" = "#009E73",
          "Step 7 vs 8" = "#0072B2"
        )
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 108), breaks = seq(0, 100, 20),
        labels = function(x) paste0(x, "%")
      ) +
      ggplot2::labs(
        title = "Top candidate overlap across Steps 6, 7, and 8",
        subtitle = "Top20, Top50, and Top100 are calculated on the same filtered UTR set",
        x = "Candidate-list size", y = "Overlap within each list", fill = NULL
      ) +
      theme_comparison()
    save_png("26_step6_step7_step8_topn_overlap.png", p_pair_overlap, 11.2, 6.5)
    plots[["26"]] <- p_pair_overlap
  }
}

# 27: Step-8 profile heatmap. Bins follow their experimental numbering,
# bin1-to-bin6, so the highest-expression population appears on the left.
# Neutral relative enrichment is 0.
relative_columns <- paste0("bin", 1:6, "_log2_relative_enrichment")
if (all(relative_columns %in% colnames(comparison)) && nrow(comparison) > 0) {
  display_n <- min(50, nrow(comparison))
  profile_candidates <- comparison[
    order(comparison$step6_rank_filtered, comparison$step7_rank_filtered),
    , drop = FALSE
  ]
  profile_candidates <- head(profile_candidates, display_n)
  profile_candidates$display_id <- paste0(
    profile_candidates$variant_id,
    " [H rank ", profile_candidates$step6_rank_filtered, "]"
  )
  profile_heat <- do.call(
    rbind,
    lapply(seq_len(nrow(profile_candidates)), function(row_number) {
      data.frame(
        variant_id = profile_candidates$display_id[[row_number]],
        bin = factor(paste0("bin", 1:6), levels = paste0("bin", 1:6)),
        log2_relative_enrichment = as.numeric(unlist(
          profile_candidates[row_number, relative_columns], use.names = FALSE
        )),
        stringsAsFactors = FALSE
      )
    })
  )
  profile_heat$variant_id <- factor(
    profile_heat$variant_id, levels = rev(profile_candidates$display_id)
  )
  heat_limit <- max(abs(profile_heat$log2_relative_enrichment), na.rm = TRUE)
  heat_limit <- max(1, min(heat_limit, 4))
  p_relative_heat <- ggplot2::ggplot(
    profile_heat,
    ggplot2::aes(x = bin, y = variant_id, fill = log2_relative_enrichment)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-heat_limit, heat_limit), oob = scales::squish
    ) +
    ggplot2::labs(
      title = "Step-8 relative-enrichment profiles of Top High15 candidates",
      subtitle = "bin1 (highest fluorescence) is left; white means neutral (p/w = 1)",
      x = "FACS bin (bin1 = highest mCherry)", y = "UTR",
      fill = "log2(p / w)"
    ) +
    theme_comparison() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6.8)
    )
  save_png(
    "27_top_high15_relative_enrichment_heatmap.png", p_relative_heat,
    10.8, max(7.4, display_n * 0.20)
  )
  plots[["27"]] <- p_relative_heat
}

# 28: median Step-8 profiles. This summarizes whether the Step-6 and Step-7
# top lists both shift toward the high-expression bins without averaging raw
# sequencing depth into the phenotype profile.
if (all(relative_columns %in% colnames(comparison)) && nrow(comparison) > 0) {
  top_k <- min(as.integer(actual_top_n), nrow(comparison))
  group_rows <- list(
    "All eligible UTRs" = comparison,
    "Step 6 Top list" = comparison[
      comparison$step6_rank_filtered <= top_k, , drop = FALSE
    ],
    "Step 7 Top list" = comparison[
      comparison$step7_rank_filtered <= top_k, , drop = FALSE
    ]
  )
  if (nrow(reference) == 1) {
    group_rows[["original"]] <- reference
  }
  median_profiles <- do.call(
    rbind,
    lapply(names(group_rows), function(group_name) {
      table <- group_rows[[group_name]]
      data.frame(
        group = group_name,
        bin = factor(paste0("bin", 1:6), levels = paste0("bin", 1:6)),
        median_log2_relative_enrichment = vapply(
          relative_columns,
          function(column) stats::median(table[[column]], na.rm = TRUE),
          numeric(1)
        ),
        stringsAsFactors = FALSE
      )
    })
  )
  median_profiles$group <- factor(
    median_profiles$group,
    levels = c("All eligible UTRs", "Step 7 Top list", "Step 6 Top list", "original")
  )
  p_median_profile <- ggplot2::ggplot(
    median_profiles,
    ggplot2::aes(
      x = bin, y = median_log2_relative_enrichment,
      color = group, group = group
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = colors[["muted"]]) +
    ggplot2::geom_line(linewidth = 1.15) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_color_manual(
      values = c(
        "All eligible UTRs" = "#7F8C8D",
        "Step 7 Top list" = "#0072B2",
        "Step 6 Top list" = "#D55E00",
        original = "#B2182B"
      ), drop = FALSE
    ) +
    ggplot2::labs(
      title = "Median relative-enrichment profiles of candidate groups",
      subtitle = "A high-expression profile is elevated toward bin1 on the left",
      x = "FACS bin (bin1 = highest mCherry)",
      y = "Median log2(relative enrichment p/w)", color = NULL
    ) +
    theme_comparison()
  save_png("28_group_median_relative_enrichment_profiles.png", p_median_profile)
  plots[["28"]] <- p_median_profile
}

# 29: individual top-candidate profiles. Unlike the group median, this figure
# exposes broad, monotonic, or irregular shapes for each selected UTR.
if (all(relative_columns %in% colnames(comparison)) && nrow(comparison) > 0) {
  individual_n <- min(24, nrow(comparison))
  individual_candidates <- comparison[
    order(comparison$step6_rank_filtered, comparison$step7_rank_filtered),
    , drop = FALSE
  ]
  individual_candidates <- head(individual_candidates, individual_n)
  individual_candidates$display_id <- paste0(
    individual_candidates$variant_id,
    " | H rank ", individual_candidates$step6_rank_filtered
  )
  individual_profiles <- do.call(
    rbind,
    lapply(seq_len(nrow(individual_candidates)), function(row_number) {
      data.frame(
        variant_id = individual_candidates$display_id[[row_number]],
        bin = factor(paste0("bin", 1:6), levels = paste0("bin", 1:6)),
        bin_number = 1:6,
        log2_relative_enrichment = as.numeric(unlist(
          individual_candidates[row_number, relative_columns], use.names = FALSE
        )),
        stringsAsFactors = FALSE
      )
    })
  )
  individual_profiles$variant_id <- factor(
    individual_profiles$variant_id,
    levels = individual_candidates$display_id
  )
  p_individual_profiles <- ggplot2::ggplot(
    individual_profiles,
    ggplot2::aes(
      x = bin, y = log2_relative_enrichment,
      group = variant_id
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = colors[["muted"]],
      linewidth = 0.45
    ) +
    ggplot2::geom_line(color = "#D55E00", linewidth = 0.85) +
    ggplot2::geom_point(color = "#D55E00", size = 1.7) +
    ggplot2::facet_wrap(~variant_id, ncol = 4) +
    ggplot2::labs(
      title = "Individual Step-8 profiles of Top High15 candidates",
      subtitle = "Top 24 by Step 6; bin1 (highest mCherry) is shown at the left",
      x = "FACS bin (bin1 = highest mCherry)",
      y = "log2(relative enrichment p/w)"
    ) +
    theme_comparison() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 8, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
  save_png(
    "29_top_high15_individual_relative_enrichment_profiles.png",
    p_individual_profiles, 13.2, max(7.5, ceiling(individual_n / 4) * 1.65)
  )
  plots[["29"]] <- p_individual_profiles
}

plot_statistics <- data.frame(
  metric = c(
    "read_supported_utr", "spearman_score_vs_high15_probability", "configured_top_n",
    "actual_top_n", "top_n_consensus", "top_n_score_only",
    "top_n_top15_only", "step6_step8_rank_mismatch_count"
  ),
  value = c(
    nrow(comparison), rho, configured_top_n, actual_top_n,
    consensus_n, score_only_n, top15_only_n,
    summary_value("step6_step8_rank_mismatch_count")
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
