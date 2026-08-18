#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_sortseq.R PROJECT_RESULTS_DIR OUTPUT_DIR")
}

required_packages <- c("ggplot2", "scales", "patchwork", "viridisLite", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them before plotting."
  )
}

project_results <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sortseq_dir <- file.path(project_results, "sortseq")
rescue_dir <- file.path(project_results, "index_rescue")
libraryqc_dir <- file.path(project_results, "library_qc")

read_tsv <- function(path) {
  read.delim(
    path, header = TRUE, sep = "\t", quote = "", check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

full_path <- file.path(sortseq_dir, "utr_results_full.tsv")
sample_qc_path <- file.path(sortseq_dir, "sample_qc.tsv")
if (!file.exists(full_path) || !file.exists(sample_qc_path)) {
  stop("Sort-seq results are missing under: ", sortseq_dir)
}
results <- read_tsv(full_path)
sample_qc <- read_tsv(sample_qc_path)
results$variant_id <- as.character(results$variant_id)

required_result_columns <- c(
  "variant_id", "pass_coverage", "expected_bin_score", "high15_enrichment",
  "expression_tier", "estimated_rank", "unsorted_frequency",
  "reconstructed_gate_frequency"
)
missing_result_columns <- setdiff(required_result_columns, colnames(results))
if (length(missing_result_columns) > 0) {
  stop("Missing Sort-seq columns: ", paste(missing_result_columns, collapse = ", "))
}
for (column in c(
  "expected_bin_score", "high15_enrichment", "estimated_rank",
  "high15_probability", "unsorted_frequency", "reconstructed_gate_frequency",
  "unsorted_count", "total_6bin_count", "detected_in_n_bins",
  paste0("bin", 1:6, "_count"),
  paste0("bin", 1:6, "_probability")
)) {
  results[[column]] <- suppressWarnings(as.numeric(results[[column]]))
}
as_bool <- function(values) {
  tolower(as.character(values)) %in% c("true", "t", "1")
}
results$pass_coverage <- tolower(as.character(results$pass_coverage)) %in% c("true", "t", "1")
if ("high_candidate_flag" %in% colnames(results)) {
  results$high_candidate_flag <- as_bool(results$high_candidate_flag)
}
if ("is_reference_variant" %in% colnames(results)) {
  results$is_reference_variant <- as_bool(results$is_reference_variant)
} else {
  results$is_reference_variant <- tolower(trimws(results$variant_id)) %in% c("original", "orginal")
}
if (!("reference_display_label" %in% colnames(results))) {
  results$reference_display_label <- results$variant_id
}
passing <- results[results$pass_coverage & is.finite(results$expected_bin_score), , drop = FALSE]
if (nrow(passing) == 0) {
  stop("No UTR passed the coverage filter.")
}
reference <- passing[passing$is_reference_variant, , drop = FALSE]
if (nrow(reference) > 1) {
  stop("More than one reference variant was identified: ", paste(reference$variant_id, collapse = ", "))
}

colors <- c(
  orange = "#D55E00", gold = "#E69F00", sky = "#56B4E9",
  blue = "#0072B2", green = "#009E73", purple = "#7B61A8",
  ink = "#24323D", muted = "#667580", grid = "#DDE5EA"
)

theme_sortseq <- function() {
  ggplot2::theme_minimal(base_size = 14, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = colors[["ink"]], face = "bold", size = 20,
        margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        color = colors[["muted"]], size = 13, lineheight = 1.15,
        margin = ggplot2::margin(b = 10)
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
      panel.grid.major = ggplot2::element_line(color = colors[["grid"]], linewidth = 0.35),
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

plots <- list()

# 01: Undetermined rescue QC
rescue_summary_path <- file.path(rescue_dir, "rescue_summary.csv")
if (file.exists(rescue_summary_path)) {
  rescue <- read.csv(rescue_summary_path, stringsAsFactors = FALSE, check.names = FALSE)
  rescue_samples <- rescue[rescue$category == "sample", , drop = FALSE]
  rescue_samples$name <- factor(rescue_samples$name, levels = rescue_samples$name)
  p_rescue <- ggplot2::ggplot(
    rescue_samples,
    ggplot2::aes(x = name, y = read_pairs_or_reads, fill = name)
  ) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%s\n(%.1f%%)", scales::comma(read_pairs_or_reads), percent_of_undetermined)
      ),
      vjust = -0.2, size = 3.7, color = colors[["ink"]]
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.16))
    ) +
    ggplot2::scale_fill_manual(values = rep(viridisLite::viridis(7), length.out = nrow(rescue_samples))) +
    ggplot2::labs(
      title = "Dual-index rescue from Undetermined FASTQ",
      subtitle = "Only unique nearest i7+i5 pairs within the configured mismatch limit were reassigned",
      x = NULL, y = "Rescued read pairs / reads"
    ) +
    theme_sortseq()
  save_png("01_undetermined_rescue.png", p_rescue)
  plots[[length(plots) + 1]] <- p_rescue
}

# 02: sample mapping/depth QC
sample_qc$assigned_reference_reads <- suppressWarnings(as.numeric(sample_qc$assigned_reference_reads))
sample_qc$detected_reference_variants <- suppressWarnings(as.numeric(sample_qc$detected_reference_variants))
sample_qc$sample_id <- factor(sample_qc$sample_id, levels = sample_qc$sample_id)
p_sample <- ggplot2::ggplot(
  sample_qc,
  ggplot2::aes(x = sample_id, y = assigned_reference_reads, fill = sample_type)
) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(assigned_reference_reads)),
    vjust = -0.25, size = 3.5, color = colors[["ink"]]
  ) +
  ggplot2::scale_y_log10(labels = scales::comma, expand = ggplot2::expansion(mult = c(0.03, 0.16))) +
  ggplot2::scale_fill_manual(values = c(bin = colors[["blue"]], unsorted = colors[["green"]])) +
  ggplot2::labs(
    title = "Assigned 5'UTR reads by sample",
    subtitle = "Counts after dual-index rescue and reference-defined amplicon assignment",
    x = NULL, y = "Assigned reference reads (log scale)", fill = "Sample type"
  ) +
  theme_sortseq()
save_png("02_sample_assigned_reads.png", p_sample)
plots[[length(plots) + 1]] <- p_sample

# 03: score and high-bin enrichment
tier_order <- c(
  "top_1pct", "top_5pct", "top_10pct", "upper_25pct",
  "middle_50pct", "lower_25pct", "bottom_10pct"
)
passing$expression_tier <- factor(passing$expression_tier, levels = tier_order)
reference <- passing[passing$is_reference_variant, , drop = FALSE]
p_score <- ggplot2::ggplot(
  passing,
  ggplot2::aes(x = expected_bin_score, y = high15_enrichment, color = expression_tier)
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = colors[["muted"]]) +
  ggplot2::geom_point(alpha = 0.65, size = 1.8) +
  ggplot2::scale_color_manual(
    values = setNames(viridisLite::viridis(length(tier_order), option = "C"), tier_order),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "5'UTR fluorescence score and high-bin enrichment",
    subtitle = "bin1 is highest mCherry; bin6 is lowest within the mCherry+/GFP- gate",
    caption = "The score is an exploratory fluorescence-associated rank, not an exact translation-rate measurement.",
    x = "Expected bin score (6 = high, 1 = low)",
    y = "Top-15% enrichment (bin1 + bin2; 1 = pool average)",
    color = "Rank tier"
  ) +
  theme_sortseq()
if (nrow(reference) == 1) {
  p_score <- p_score +
    ggplot2::geom_point(
      data = reference,
      ggplot2::aes(x = expected_bin_score, y = high15_enrichment),
      inherit.aes = FALSE, shape = 8, size = 5.2, stroke = 1.2,
      color = "#B2182B"
    ) +
    ggplot2::geom_text(
      data = reference,
      ggplot2::aes(
        x = expected_bin_score, y = high15_enrichment,
        label = paste0(reference_display_label, " [reference]")
      ),
      inherit.aes = FALSE, nudge_x = 0.06, nudge_y = 0.08,
      hjust = 0, color = "#B2182B", fontface = "bold", size = 3.8
    )
}
save_png("03_score_vs_high15_enrichment.png", p_score)
plots[[length(plots) + 1]] <- p_score

# 04: top-UTR bin probability heatmap
top_n <- min(50, nrow(passing))
top <- passing[order(passing$estimated_rank), , drop = FALSE][seq_len(top_n), , drop = FALSE]
reference_added <- FALSE
if (nrow(reference) == 1 && !(reference$variant_id[[1]] %in% top$variant_id)) {
  top <- rbind(top, reference)
  reference_added <- TRUE
}
top$display_id <- ifelse(
  top$is_reference_variant,
  paste0(top$reference_display_label, " [reference]"),
  top$variant_id
)
probability_columns <- paste0("bin", 1:6, "_probability")
heat <- do.call(
  rbind,
  lapply(seq_len(nrow(top)), function(row_number) {
    data.frame(
      variant_id = top$display_id[[row_number]],
      bin = factor(paste0("bin", 1:6), levels = paste0("bin", 1:6)),
      probability = as.numeric(top[row_number, probability_columns]),
      rank = top$estimated_rank[[row_number]],
      is_reference = top$is_reference_variant[[row_number]],
      stringsAsFactors = FALSE
    )
  })
)
heat$variant_id <- factor(heat$variant_id, levels = rev(top$display_id))
p_heat <- ggplot2::ggplot(
  heat,
  ggplot2::aes(x = bin, y = variant_id, fill = probability)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.15) +
  ggplot2::scale_fill_gradientn(
    colors = viridisLite::viridis(256, option = "B"), labels = scales::percent
  ) +
  ggplot2::labs(
    title = if (reference_added) {
      sprintf("Top %d UTR fluorescence-bin distributions + reference", top_n)
    } else {
      sprintf("Top %d UTR fluorescence-bin distributions", top_n)
    },
    subtitle = "Cell-fraction and sequencing-depth corrected; bin1 = highest mCherry",
    x = "FACS bin", y = "UTR", fill = "Estimated\ncell fraction"
  ) +
  theme_sortseq() +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 7.5))
if (any(heat$is_reference)) {
  p_heat <- p_heat + ggplot2::geom_tile(
    data = heat[heat$is_reference, , drop = FALSE],
    fill = NA, color = "#B2182B", linewidth = 0.9
  )
}
save_png("04_top_utr_bin_heatmap.png", p_heat, 10.2, max(7.2, top_n * 0.19))
plots[[length(plots) + 1]] <- p_heat

# 05: tier counts
tier_counts <- as.data.frame(table(passing$expression_tier), stringsAsFactors = FALSE)
colnames(tier_counts) <- c("tier", "variants")
tier_counts <- tier_counts[tier_counts$variants > 0, , drop = FALSE]
p_tier <- ggplot2::ggplot(
  tier_counts,
  ggplot2::aes(x = tier, y = variants, fill = tier)
) +
  ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = scales::comma(variants)), vjust = -0.3, size = 4) +
  ggplot2::scale_fill_manual(values = setNames(viridisLite::viridis(length(tier_order), option = "C"), tier_order)) +
  ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(
    title = "Exploratory expression-rank tiers",
    subtitle = sprintf("%s UTRs passed coverage filters", scales::comma(nrow(passing))),
    x = NULL, y = "Number of UTRs"
  ) +
  theme_sortseq() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 28, hjust = 1))
save_png("05_expression_tiers.png", p_tier)
plots[[length(plots) + 1]] <- p_tier

# 06: whole unsorted versus reconstructed target-gate composition
composition <- passing[
  passing$unsorted_frequency > 0 & passing$reconstructed_gate_frequency > 0,
  , drop = FALSE
]
p_gate <- ggplot2::ggplot(
  composition,
  ggplot2::aes(x = unsorted_frequency, y = reconstructed_gate_frequency)
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = colors[["muted"]]) +
  ggplot2::geom_point(alpha = 0.48, size = 1.5, color = colors[["blue"]]) +
  ggplot2::scale_x_log10(labels = scales::label_scientific()) +
  ggplot2::scale_y_log10(labels = scales::label_scientific()) +
  ggplot2::labs(
    title = "Whole unsorted versus reconstructed mCherry+/GFP- gate",
    subtitle = "Deviation from the diagonal reflects gate representation as well as technical sampling",
    x = "UTR frequency in whole unsorted",
    y = "Reconstructed UTR frequency in target gate"
  ) +
  theme_sortseq()
reference_composition <- composition[composition$is_reference_variant, , drop = FALSE]
if (nrow(reference_composition) == 1) {
  p_gate <- p_gate +
    ggplot2::geom_point(
      data = reference_composition,
      ggplot2::aes(x = unsorted_frequency, y = reconstructed_gate_frequency),
      inherit.aes = FALSE, shape = 8, size = 5.2, stroke = 1.2,
      color = "#B2182B"
    ) +
    ggplot2::geom_text(
      data = reference_composition,
      ggplot2::aes(
        x = unsorted_frequency, y = reconstructed_gate_frequency,
        label = paste0(reference_display_label, " [reference]")
      ),
      inherit.aes = FALSE,
      hjust = -0.08, vjust = -0.5,
      color = "#B2182B", fontface = "bold", size = 3.8
    )
}
save_png("06_unsorted_vs_target_gate.png", p_gate)
plots[[length(plots) + 1]] <- p_gate

# 07: reference position in the score distribution
if (nrow(reference) == 1) {
  p_reference <- ggplot2::ggplot(
    passing,
    ggplot2::aes(x = expected_bin_score)
  ) +
    ggplot2::geom_histogram(bins = 45, fill = colors[["sky"]], color = "white") +
    ggplot2::geom_vline(
      xintercept = reference$expected_bin_score[[1]],
      color = "#B2182B", linewidth = 1.2, linetype = "dashed"
    ) +
    ggplot2::annotate(
      "label",
      x = reference$expected_bin_score[[1]], y = Inf,
      label = sprintf(
        "%s [reference]\nscore %.3f | rank %.0f | percentile %.2f",
        reference$reference_display_label[[1]], reference$expected_bin_score[[1]],
        reference$estimated_rank[[1]], reference$expression_percentile[[1]]
      ),
      vjust = 1.25, hjust = ifelse(reference$expected_bin_score[[1]] > median(passing$expected_bin_score), 1.05, -0.05),
      color = "#B2182B", fill = "white", fontface = "bold", size = 3.8
    ) +
    ggplot2::labs(
      title = "Reference position in the fluorescence-score distribution",
      subtitle = "The dashed red line marks the original/orginal control UTR",
      x = "Expected bin score", y = "Coverage-passing UTRs"
    ) +
    theme_sortseq()
  save_png("07_reference_score_position.png", p_reference)
  plots[[length(plots) + 1]] <- p_reference
}

# Strict-candidate figures. These are generated when either the v0.1.5 strict
# flags are present in the full result or a server-local
# high_confidence_candidates.tsv was created with the documented filter.
strict_candidate_path <- file.path(sortseq_dir, "high_confidence_candidates.tsv")
strict_candidate_ids <- character(0)
strict_candidates_recomputed <- FALSE
if (file.exists(strict_candidate_path)) {
  strict_candidate_file <- read_tsv(strict_candidate_path)
  if ("variant_id" %in% colnames(strict_candidate_file)) {
    strict_candidate_ids <- unique(as.character(strict_candidate_file$variant_id))
  }
} else if ("high_confidence_candidate_flag" %in% colnames(results)) {
  strict_candidate_ids <- as.character(
    results$variant_id[as_bool(results$high_confidence_candidate_flag)]
  )
}

# Backward-compatible fallback for v0.1.4 result tables and for an absent or
# empty manually generated candidate file. Recompute the documented strict
# filter directly from utr_results_full.tsv so the R-only plot step is enough.
if (length(strict_candidate_ids) == 0) {
  fallback_unsorted_cutoff <- max(
    100, round(0.10 * median(results$unsorted_count, na.rm = TRUE))
  )
  fallback_total_cutoff <- max(
    500, round(0.10 * median(results$total_6bin_count, na.rm = TRUE))
  )
  fallback_high_bin_count <- results$bin1_count + results$bin2_count
  fallback_max_probability <- apply(
    results[, paste0("bin", 1:6, "_probability"), drop = FALSE],
    1, max, na.rm = TRUE
  )
  fallback_strict_pass <-
    results$unsorted_count >= fallback_unsorted_cutoff &
    results$total_6bin_count >= fallback_total_cutoff &
    fallback_high_bin_count >= 50 &
    results$detected_in_n_bins >= 3
  fallback_jackpot <-
    fallback_max_probability >= 0.85 & results$detected_in_n_bins <= 2
  fallback_high_candidate <- if ("high_candidate_flag" %in% colnames(results)) {
    as_bool(results$high_candidate_flag)
  } else {
    rep(FALSE, nrow(results))
  }
  strict_candidate_ids <- as.character(
    results$variant_id[
      fallback_strict_pass & fallback_high_candidate & !fallback_jackpot
    ]
  )
  strict_candidates_recomputed <- TRUE
  message(
    "Read-supported exploratory candidates recomputed from utr_results_full.tsv: ",
    length(strict_candidate_ids)
  )
}

if (length(strict_candidate_ids) > 0) {
  strict_output_dir <- file.path(output_dir, "strict")
  dir.create(strict_output_dir, recursive = TRUE, showWarnings = FALSE)
  save_strict_png <- function(filename, plot, width = 11.8, height = 6.6) {
    ragg::agg_png(
      filename = file.path(strict_output_dir, filename), width = width, height = height,
      units = "in", res = 320, background = "white"
    )
    print(plot)
    grDevices::dev.off()
  }

  if (!strict_candidates_recomputed && "strict_coverage_pass" %in% colnames(results)) {
    results$strict_coverage_pass <- as_bool(results$strict_coverage_pass)
    strict_passing <- results[
      results$strict_coverage_pass & is.finite(results$expected_bin_score),
      , drop = FALSE
    ]
    strict_unsorted_cutoff <- suppressWarnings(max(results$strict_unsorted_cutoff, na.rm = TRUE))
    strict_total_cutoff <- suppressWarnings(max(results$strict_total_6bin_cutoff, na.rm = TRUE))
  } else {
    strict_unsorted_cutoff <- max(100, round(0.10 * median(results$unsorted_count, na.rm = TRUE)))
    strict_total_cutoff <- max(500, round(0.10 * median(results$total_6bin_count, na.rm = TRUE)))
    high_bin_raw_count <- results$bin1_count + results$bin2_count
    strict_mask <-
      results$unsorted_count >= strict_unsorted_cutoff &
      results$total_6bin_count >= strict_total_cutoff &
      high_bin_raw_count >= 50 &
      results$detected_in_n_bins >= 3
    strict_passing <- results[
      strict_mask & is.finite(results$expected_bin_score),
      , drop = FALSE
    ]
  }
  strong_mask <-
    results$unsorted_count >= 200 &
    results$total_6bin_count >= 1000 &
    (results$bin1_count + results$bin2_count) >= 100 &
    results$detected_in_n_bins >= 3
  strong_passing <- results[
    strong_mask & is.finite(results$expected_bin_score),
    , drop = FALSE
  ]
  strict_candidates <- results[
    results$variant_id %in% strict_candidate_ids & is.finite(results$expected_bin_score),
    , drop = FALSE
  ]
  strict_candidates <- strict_candidates[
    order(strict_candidates$expected_bin_score, decreasing = TRUE),
    , drop = FALSE
  ]
  strict_plots <- list()

  # 08: keep the full passing set as faint context, while the strict universe
  # and final high-confidence candidates are visually separated.
  p_strict_score <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = colors[["muted"]]) +
    ggplot2::geom_point(
      data = passing,
      ggplot2::aes(x = expected_bin_score, y = high15_enrichment),
      color = "#C7CED3", alpha = 0.38, size = 1.25
    ) +
    ggplot2::geom_point(
      data = strict_passing,
      ggplot2::aes(x = expected_bin_score, y = high15_enrichment),
      color = colors[["blue"]], alpha = 0.52, size = 1.55
    ) +
    ggplot2::geom_point(
      data = strong_passing,
      ggplot2::aes(x = expected_bin_score, y = high15_enrichment),
      color = colors[["purple"]], alpha = 0.72, size = 1.85
    ) +
    ggplot2::geom_point(
      data = strict_candidates,
      ggplot2::aes(x = expected_bin_score, y = high15_enrichment),
      color = colors[["orange"]], alpha = 0.90, size = 2.5
    ) +
    ggplot2::labs(
      title = "Read-supported exploratory 5'UTR candidates",
      subtitle = sprintf(
        "Gray = permissive; blue = supported (%s); purple = strong support (%s); orange = supported high candidates (%s)",
        scales::comma(nrow(strict_passing)), scales::comma(nrow(strong_passing)),
        scales::comma(nrow(strict_candidates))
      ),
      caption = sprintf(
        "Supported cutoffs: unsorted >= %s; total six bins >= %s; bin1+2 >= 50; detected bins >= 3. No strong-support high candidate was assumed.",
        scales::comma(strict_unsorted_cutoff), scales::comma(strict_total_cutoff)
      ),
      x = "Expected bin score (6 = high, 1 = low)",
      y = "Top-15% enrichment (bin1 + bin2; 1 = pool average)"
    ) +
    theme_sortseq()
  if (nrow(reference) == 1) {
    p_strict_score <- p_strict_score +
      ggplot2::geom_point(
        data = reference,
        ggplot2::aes(x = expected_bin_score, y = high15_enrichment),
        shape = 8, size = 5.2, stroke = 1.2, color = "#B2182B"
      ) +
      ggplot2::geom_text(
        data = reference,
        ggplot2::aes(
          x = expected_bin_score, y = high15_enrichment,
          label = paste0(reference_display_label, " [reference]")
        ),
        nudge_x = 0.06, nudge_y = 0.08, hjust = 0,
        color = "#B2182B", fontface = "bold", size = 3.8
      )
  }
  save_strict_png("08_strict_score_vs_high15_enrichment.png", p_strict_score)
  strict_plots[[length(strict_plots) + 1]] <- p_strict_score

  # Use the same UTR order for the probability and enrichment heatmaps.
  top_strict_n <- min(50, nrow(strict_candidates))
  top_strict <- strict_candidates[seq_len(top_strict_n), , drop = FALSE]
  strict_reference_added <- FALSE
  if (nrow(reference) == 1 && !(reference$variant_id[[1]] %in% top_strict$variant_id)) {
    top_strict <- rbind(top_strict, reference)
    strict_reference_added <- TRUE
  }
  top_strict$display_id <- ifelse(
    top_strict$is_reference_variant,
    paste0(top_strict$reference_display_label, " [reference]"),
    top_strict$variant_id
  )
  strict_heat <- do.call(
    rbind,
    lapply(seq_len(nrow(top_strict)), function(row_number) {
      data.frame(
        variant_id = top_strict$display_id[[row_number]],
        bin = factor(paste0("bin", 1:6), levels = paste0("bin", 1:6)),
        probability = as.numeric(top_strict[row_number, probability_columns]),
        is_reference = top_strict$is_reference_variant[[row_number]],
        stringsAsFactors = FALSE
      )
    })
  )
  strict_heat$variant_id <- factor(strict_heat$variant_id, levels = rev(top_strict$display_id))
  p_strict_probability <- ggplot2::ggplot(
    strict_heat,
    ggplot2::aes(x = bin, y = variant_id, fill = probability)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradientn(
      colors = viridisLite::viridis(256, option = "B"), labels = scales::percent
    ) +
    ggplot2::labs(
      title = sprintf("Top %d read-supported high-score UTR bin probabilities", top_strict_n),
      subtitle = if (strict_reference_added) {
        "Strict candidates plus reference; absolute estimated cell fraction"
      } else {
        "Strict candidates; absolute estimated cell fraction"
      },
      x = "FACS bin", y = "UTR", fill = "Estimated\ncell fraction"
    ) +
    theme_sortseq() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 7.5))
  if (any(strict_heat$is_reference)) {
    p_strict_probability <- p_strict_probability + ggplot2::geom_tile(
      data = strict_heat[strict_heat$is_reference, , drop = FALSE],
      fill = NA, color = "#B2182B", linewidth = 0.9
    )
  }
  save_strict_png(
    "09_strict_top_utr_bin_probability_heatmap.png",
    p_strict_probability, 10.2, max(7.2, nrow(top_strict) * 0.19)
  )
  strict_plots[[length(strict_plots) + 1]] <- p_strict_probability

  # Relative enrichment removes the visual advantage of wider bins. A value
  # of zero means p(bin|UTR) equals the sorter population fraction.
  sample_input_path <- file.path(sortseq_dir, "input", "samples.tsv")
  default_bin_fraction <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.20)
  bin_fraction <- default_bin_fraction
  if (file.exists(sample_input_path)) {
    sample_input <- read_tsv(sample_input_path)
    sample_bins <- sample_input[tolower(sample_input$sample_type) == "bin", , drop = FALSE]
    sample_bins <- sample_bins[order(as.numeric(sample_bins$bin_number)), , drop = FALSE]
    loaded_fraction <- suppressWarnings(as.numeric(sample_bins$population_fraction))
    if (length(loaded_fraction) == 6 && all(is.finite(loaded_fraction)) && all(loaded_fraction > 0)) {
      bin_fraction <- loaded_fraction / sum(loaded_fraction)
    }
  }
  strict_heat$bin_fraction <- bin_fraction[as.integer(strict_heat$bin)]
  strict_heat$log2_enrichment <- log2(
    pmax(strict_heat$probability, 1e-9) / strict_heat$bin_fraction
  )
  p_strict_enrichment <- ggplot2::ggplot(
    strict_heat,
    ggplot2::aes(x = bin, y = variant_id, fill = log2_enrichment)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-3, 3), oob = scales::squish
    ) +
    ggplot2::labs(
      title = sprintf("Top %d read-supported high-score UTR bin enrichment", top_strict_n),
      subtitle = "log2(probability / sorter-bin fraction); red = enriched, blue = depleted",
      x = "FACS bin", y = "UTR", fill = "log2\nbin enrichment"
    ) +
    theme_sortseq() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 7.5))
  if (any(strict_heat$is_reference)) {
    p_strict_enrichment <- p_strict_enrichment + ggplot2::geom_tile(
      data = strict_heat[strict_heat$is_reference, , drop = FALSE],
      fill = NA, color = "#111111", linewidth = 0.9
    )
  }
  save_strict_png(
    "10_strict_top_utr_bin_enrichment_heatmap.png",
    p_strict_enrichment, 10.2, max(7.2, nrow(top_strict) * 0.19)
  )
  strict_plots[[length(strict_plots) + 1]] <- p_strict_enrichment

  p_strict_gate <- ggplot2::ggplot() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = colors[["muted"]]) +
    ggplot2::geom_point(
      data = composition,
      ggplot2::aes(x = unsorted_frequency, y = reconstructed_gate_frequency),
      color = "#C7CED3", alpha = 0.35, size = 1.2
    ) +
    ggplot2::geom_point(
      data = strict_candidates,
      ggplot2::aes(x = unsorted_frequency, y = reconstructed_gate_frequency),
      color = colors[["orange"]], alpha = 0.88, size = 2.5
    ) +
    ggplot2::scale_x_log10(labels = scales::label_scientific()) +
    ggplot2::scale_y_log10(labels = scales::label_scientific()) +
    ggplot2::labs(
      title = "Read-supported exploratory candidates in the target gate",
      subtitle = "Gray = all coverage-passing UTRs; orange = supported high-score candidates",
      x = "UTR frequency in whole unsorted",
      y = "Reconstructed UTR frequency in target gate"
    ) +
    theme_sortseq()
  if (nrow(reference_composition) == 1) {
    p_strict_gate <- p_strict_gate + ggplot2::geom_point(
      data = reference_composition,
      ggplot2::aes(x = unsorted_frequency, y = reconstructed_gate_frequency),
      shape = 8, size = 5.2, stroke = 1.2, color = "#B2182B"
    )
  }
  save_strict_png("11_strict_unsorted_vs_target_gate.png", p_strict_gate)
  strict_plots[[length(strict_plots) + 1]] <- p_strict_gate

  if (nrow(reference) == 1 && nrow(strict_passing) > 0) {
    p_strict_reference <- ggplot2::ggplot(
      strict_passing,
      ggplot2::aes(x = expected_bin_score)
    ) +
      ggplot2::geom_histogram(bins = 45, fill = colors[["sky"]], color = "white") +
      ggplot2::geom_vline(
        xintercept = reference$expected_bin_score[[1]],
        color = "#B2182B", linewidth = 1.2, linetype = "dashed"
      ) +
      ggplot2::labs(
        title = "Reference position after strict coverage filtering",
        subtitle = sprintf(
          "%s strict-coverage UTRs; red line = original/orginal",
          scales::comma(nrow(strict_passing))
        ),
        x = "Expected bin score", y = "Strict-coverage UTRs"
      ) +
      theme_sortseq()
    save_strict_png("12_strict_reference_score_position.png", p_strict_reference)
    strict_plots[[length(strict_plots) + 1]] <- p_strict_reference
  }

  strict_statistics <- data.frame(
    metric = c(
      "strict_unsorted_cutoff", "strict_total_6bin_cutoff",
      "supported_coverage_passing_utr", "strong_coverage_passing_utr",
      "supported_high_candidate_count"
    ),
    value = c(
      strict_unsorted_cutoff, strict_total_cutoff,
      nrow(strict_passing), nrow(strong_passing), nrow(strict_candidates)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    strict_statistics,
    file.path(strict_output_dir, "strict_plot_statistics.csv"),
    quote = FALSE, row.names = FALSE
  )
  grDevices::pdf(
    file.path(strict_output_dir, "strict_candidate_figures.pdf"),
    width = 13.333, height = 7.5, onefile = TRUE,
    family = "Helvetica", paper = "special"
  )
  for (plot in strict_plots) {
    print(plot)
  }
  grDevices::dev.off()
  plots <- c(plots, strict_plots)
}

# Optional NGS_LibraryQC mapping summary used for tabular audit.
library_metrics_path <- file.path(libraryqc_dir, "combined", "all_sample_metrics.csv")
if (file.exists(library_metrics_path)) {
  library_metrics <- read.csv(library_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
  write.csv(
    library_metrics,
    file.path(output_dir, "libraryqc_sample_metrics_used.csv"),
    quote = FALSE, row.names = FALSE
  )
}

statistics <- data.frame(
  metric = c(
    "reference_utr_rows", "coverage_passing_utr", "coverage_passing_percent",
    "high_candidate_flag_count", "median_expected_bin_score",
    "median_high15_enrichment"
  ),
  value = c(
    nrow(results), nrow(passing), 100 * nrow(passing) / nrow(results),
    sum(tolower(as.character(results$high_candidate_flag)) %in% c("true", "t", "1")),
    median(passing$expected_bin_score, na.rm = TRUE),
    median(passing$high15_enrichment, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
if (nrow(reference) == 1) {
  reference_statistics <- data.frame(
    metric = c(
      "reference_variant_id", "reference_expected_bin_score",
      "reference_high15_probability", "reference_high15_enrichment",
      "reference_estimated_rank", "reference_expression_percentile",
      "variants_with_score_above_reference",
      "variants_with_score_and_high15_above_reference"
    ),
    value = c(
      reference$variant_id[[1]], reference$expected_bin_score[[1]],
      reference$high15_probability[[1]], reference$high15_enrichment[[1]],
      reference$estimated_rank[[1]], reference$expression_percentile[[1]],
      sum(passing$expected_bin_score > reference$expected_bin_score[[1]], na.rm = TRUE),
      sum(
        passing$expected_bin_score > reference$expected_bin_score[[1]] &
          passing$high15_probability > reference$high15_probability[[1]],
        na.rm = TRUE
      )
    ),
    stringsAsFactors = FALSE
  )
  statistics <- rbind(statistics, reference_statistics)
}
write.csv(statistics, file.path(output_dir, "plot_statistics.csv"), quote = FALSE, row.names = FALSE)

grDevices::pdf(
  file.path(output_dir, "sortseq_qc_figures.pdf"),
  width = 13.333, height = 7.5, onefile = TRUE, family = "Helvetica", paper = "special"
)
for (plot in plots) {
  print(plot)
}
grDevices::dev.off()

message("R figures written to: ", normalizePath(output_dir, mustWork = TRUE))
