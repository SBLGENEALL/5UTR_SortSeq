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
  "unsorted_frequency", "reconstructed_gate_frequency",
  paste0("bin", 1:6, "_probability")
)) {
  results[[column]] <- suppressWarnings(as.numeric(results[[column]]))
}
results$pass_coverage <- tolower(as.character(results$pass_coverage)) %in% c("true", "t", "1")
passing <- results[results$pass_coverage & is.finite(results$expected_bin_score), , drop = FALSE]
if (nrow(passing) == 0) {
  stop("No UTR passed the coverage filter.")
}
if ("is_reference_variant" %in% colnames(passing)) {
  passing$is_reference_variant <- tolower(as.character(passing$is_reference_variant)) %in% c("true", "t", "1")
} else {
  passing$is_reference_variant <- tolower(trimws(passing$variant_id)) %in% c("original", "orginal")
}
if (!("reference_display_label" %in% colnames(passing))) {
  passing$reference_display_label <- passing$variant_id
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
