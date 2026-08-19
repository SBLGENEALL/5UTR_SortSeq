#!/usr/bin/env bash
set -Eeuo pipefail

REPO_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-}"
OPTION="${2:-}"
CONFIG_FILE="${SORTSEQ_PROJECT_CONFIG:-${REPO_DIR}/config/project.env}"

usage() {
  cat >&2 <<'EOF'
Usage:
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh preflight
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh rescue
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh libraryqc
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh analyze
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh scoring-steps
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh bimodality-qc
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh compare-metrics
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh plot
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh status
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh resources
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh reanalyze-analysis
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh reanalyze
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh full [--replace]
EOF
}

if [[ ! "${MODE}" =~ ^(preflight|rescue|libraryqc|analyze|scoring-steps|bimodality-qc|compare-metrics|plot|status|resources|reanalyze-analysis|reanalyze|full)$ ]]; then
  usage
  exit 2
fi
if [[ -n "${OPTION}" && "${OPTION}" != "--replace" ]]; then
  usage
  exit 2
fi
if [[ "${MODE}" != "full" && -n "${OPTION}" ]]; then
  usage
  exit 2
fi
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "Project config not found: ${CONFIG_FILE}" >&2
  echo "Copy config/project.env.example to a private path and edit it." >&2
  exit 2
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

: "${RESULTS_DIR:?RESULTS_DIR is required}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
WORKERS="${WORKERS:-128}"
LIBRARYQC_BATCH_SIZE="${LIBRARYQC_BATCH_SIZE:-10000}"
MAX_TOTAL_INDEX_MISMATCHES="${MAX_TOTAL_INDEX_MISMATCHES:-2}"
MAX_PER_INDEX_MISMATCHES="${MAX_PER_INDEX_MISMATCHES:-2}"
MIN_INDEX_DISTANCE_MARGIN="${MIN_INDEX_DISTANCE_MARGIN:-1}"
I5_ORIENTATION="${I5_ORIENTATION:-auto}"
PROGRESS_INTERVAL_SECONDS="${PROGRESS_INTERVAL_SECONDS:-30}"
PROGRESS_CHECK_READS="${PROGRESS_CHECK_READS:-100000}"
RESCUE_COMPRESSION_BACKEND="${RESCUE_COMPRESSION_BACKEND:-auto}"
RESCUE_PIGZ_THREADS_PER_FILE="${RESCUE_PIGZ_THREADS_PER_FILE:-4}"
NUMA_INTERLEAVE="${NUMA_INTERLEAVE:-auto}"
NESTED_WORKER_THREADS="${NESTED_WORKER_THREADS:-1}"
OVERALL_GATE_FRACTION="${OVERALL_GATE_FRACTION:-0.90}"
MIN_UNSORTED_COUNT="${MIN_UNSORTED_COUNT:-50}"
MIN_TOTAL_BIN_COUNT="${MIN_TOTAL_BIN_COUNT:-100}"
TOP_HIT_MIN_UNSORTED_COUNT="${TOP_HIT_MIN_UNSORTED_COUNT:-50}"
TOP_HIT_MIN_TOTAL_BIN_COUNT="${TOP_HIT_MIN_TOTAL_BIN_COUNT:-200}"
TOP_HIT_MIN_HIGH_BIN_COUNT="${TOP_HIT_MIN_HIGH_BIN_COUNT:-20}"
TOP_HIT_MIN_ENRICHMENT="${TOP_HIT_MIN_ENRICHMENT:-1.0}"
HIGH15_BOOTSTRAP_REPLICATES="${HIGH15_BOOTSTRAP_REPLICATES:-1000}"
HIGH15_BOOTSTRAP_SEED="${HIGH15_BOOTSTRAP_SEED:-20260819}"
HIGH15_BOOTSTRAP_LOWER_QUANTILE="${HIGH15_BOOTSTRAP_LOWER_QUANTILE:-0.10}"
HIGH15_BOOTSTRAP_TOP_N="${HIGH15_BOOTSTRAP_TOP_N:-50}"
HIGH15_BOOTSTRAP_MIN_PROBABILITY="${HIGH15_BOOTSTRAP_MIN_PROBABILITY:-0.90}"
STRICT_MIN_UNSORTED_COUNT="${STRICT_MIN_UNSORTED_COUNT:-1000}"
STRICT_MIN_TOTAL_BIN_COUNT="${STRICT_MIN_TOTAL_BIN_COUNT:-5000}"
STRICT_RELATIVE_MEDIAN_FRACTION="${STRICT_RELATIVE_MEDIAN_FRACTION:-0.10}"
STRICT_MIN_HIGH_BIN_COUNT="${STRICT_MIN_HIGH_BIN_COUNT:-200}"
STRICT_MIN_DETECTED_BINS="${STRICT_MIN_DETECTED_BINS:-3}"
JACKPOT_MAX_BIN_PROBABILITY="${JACKPOT_MAX_BIN_PROBABILITY:-0.85}"
JACKPOT_MAX_DETECTED_BINS="${JACKPOT_MAX_DETECTED_BINS:-2}"
REFERENCE_VARIANT_ID="${REFERENCE_VARIANT_ID:-auto}"
BIMODALITY_MIN_TOTAL_COUNT="${BIMODALITY_MIN_TOTAL_COUNT:-200}"
BIMODALITY_MIN_EACH_TAIL_COUNT="${BIMODALITY_MIN_EACH_TAIL_COUNT:-20}"
BIMODALITY_MIN_HIGH_TAIL_PROBABILITY="${BIMODALITY_MIN_HIGH_TAIL_PROBABILITY:-0.20}"
BIMODALITY_MIN_LOW_TAIL_PROBABILITY="${BIMODALITY_MIN_LOW_TAIL_PROBABILITY:-0.20}"
BIMODALITY_MAX_MIDDLE_PROBABILITY="${BIMODALITY_MAX_MIDDLE_PROBABILITY:-0.30}"
BIMODALITY_MAX_VALLEY_RATIO="${BIMODALITY_MAX_VALLEY_RATIO:-0.75}"
METRIC_COMPARE_MIN_TOTAL_COUNT="${METRIC_COMPARE_MIN_TOTAL_COUNT:-201}"
METRIC_COMPARE_MIN_UNSORTED_COUNT="${METRIC_COMPARE_MIN_UNSORTED_COUNT:-50}"
METRIC_COMPARE_MIN_HIGH_BIN_COUNT="${METRIC_COMPARE_MIN_HIGH_BIN_COUNT:-20}"
METRIC_COMPARE_TOP_N="${METRIC_COMPARE_TOP_N:-50}"

# Multiprocessing stages own the parallelism. Prevent BLAS/OpenMP libraries
# loaded by individual workers from multiplying 128 workers by extra threads.
export OMP_NUM_THREADS="${NESTED_WORKER_THREADS}"
export OPENBLAS_NUM_THREADS="${NESTED_WORKER_THREADS}"
export MKL_NUM_THREADS="${NESTED_WORKER_THREADS}"
export NUMEXPR_NUM_THREADS="${NESTED_WORKER_THREADS}"

require_path() {
  local path="$1"
  local label="$2"
  if [[ ! -e "${path}" ]]; then
    echo "Required ${label} not found: ${path}" >&2
    exit 2
  fi
}

require_rescue_inputs() {
  : "${RAW_DIR:?RAW_DIR is required}"
  : "${SAMPLE_SHEET:?SAMPLE_SHEET is required}"
  require_path "${RAW_DIR}" "raw-data directory"
  require_path "${SAMPLE_SHEET}" "SampleSheet"
}

detect_logical_cpus() {
  getconf _NPROCESSORS_ONLN 2>/dev/null || echo "unknown"
}

detect_physical_cores() {
  if command -v lscpu >/dev/null 2>&1; then
    lscpu -p=CORE,SOCKET 2>/dev/null \
      | awk -F, '!/^#/ {print $1 "," $2}' \
      | sort -u \
      | wc -l \
      | tr -d ' '
  else
    echo "unknown"
  fi
}

detect_sockets() {
  if command -v lscpu >/dev/null 2>&1; then
    lscpu -p=SOCKET 2>/dev/null \
      | awk -F, '!/^#/ {print $1}' \
      | sort -u \
      | wc -l \
      | tr -d ' '
  else
    echo "unknown"
  fi
}

detect_numa_nodes() {
  if command -v lscpu >/dev/null 2>&1; then
    lscpu 2>/dev/null | awk -F: '/^NUMA node\(s\)/ {gsub(/ /, "", $2); print $2; found=1} END {if (!found) print "unknown"}'
  else
    echo "unknown"
  fi
}

detect_memory_gib() {
  awk '/^MemTotal:/ {printf "%.1f", $2 / 1024 / 1024}' /proc/meminfo 2>/dev/null || echo "unknown"
}

selected_compression_backend() {
  if [[ "${RESCUE_COMPRESSION_BACKEND}" == "auto" ]]; then
    if command -v pigz >/dev/null 2>&1; then
      echo "pigz"
    else
      echo "parallel_python"
    fi
  else
    echo "${RESCUE_COMPRESSION_BACKEND}"
  fi
}

numa_interleave_enabled() {
  case "${NUMA_INTERLEAVE}" in
    1|true|yes)
      command -v numactl >/dev/null 2>&1 && numactl --interleave=all true >/dev/null 2>&1
      ;;
    0|false|no)
      return 1
      ;;
    auto)
      local nodes
      nodes="$(detect_numa_nodes)"
      command -v numactl >/dev/null 2>&1 \
        && [[ "${nodes}" =~ ^[0-9]+$ ]] \
        && (( nodes > 1 )) \
        && numactl --interleave=all true >/dev/null 2>&1
      ;;
    *)
      echo "NUMA_INTERLEAVE must be auto, true/1, or false/0" >&2
      return 2
      ;;
  esac
}

run_compute() {
  if numa_interleave_enabled; then
    numactl --interleave=all "$@"
  else
    "$@"
  fi
}

show_resources() {
  local numa_policy="disabled_or_unavailable"
  if numa_interleave_enabled; then
    numa_policy="interleave_all"
  fi
  echo "Sort-seq compute profile:"
  echo "  logical_cpus: $(detect_logical_cpus)"
  echo "  physical_cores: $(detect_physical_cores)"
  echo "  sockets: $(detect_sockets)"
  echo "  numa_nodes: $(detect_numa_nodes)"
  echo "  memory_gib: $(detect_memory_gib)"
  echo "  libraryqc_workers: ${WORKERS}"
  echo "  libraryqc_batch_size: ${LIBRARYQC_BATCH_SIZE}"
  echo "  rescue_compression: $(selected_compression_backend)"
  echo "  pigz_threads_per_output: ${RESCUE_PIGZ_THREADS_PER_FILE}"
  echo "  numa_policy: ${numa_policy}"
}

write_resource_profile() {
  local target="${RESULTS_DIR}/resource_profile.tsv"
  local numa_policy="disabled_or_unavailable"
  if numa_interleave_enabled; then
    numa_policy="interleave_all"
  fi
  {
    printf 'field\tvalue\n'
    printf 'generated_at\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %z')"
    printf 'logical_cpus\t%s\n' "$(detect_logical_cpus)"
    printf 'physical_cores\t%s\n' "$(detect_physical_cores)"
    printf 'sockets\t%s\n' "$(detect_sockets)"
    printf 'numa_nodes\t%s\n' "$(detect_numa_nodes)"
    printf 'memory_gib\t%s\n' "$(detect_memory_gib)"
    printf 'libraryqc_workers\t%s\n' "${WORKERS}"
    printf 'libraryqc_batch_size\t%s\n' "${LIBRARYQC_BATCH_SIZE}"
    printf 'rescue_compression\t%s\n' "$(selected_compression_backend)"
    printf 'pigz_threads_per_output\t%s\n' "${RESCUE_PIGZ_THREADS_PER_FILE}"
    printf 'numa_policy\t%s\n' "${numa_policy}"
  } >"${target}"
}

mkdir -p "${RESULTS_DIR}"
LOCK_FILE="${RESULTS_DIR}/.sortseq_analysis.lock"
PIPELINE_STATUS_FILE="${RESULTS_DIR}/pipeline_status.tsv"
PIPELINE_HISTORY_FILE="${RESULTS_DIR}/pipeline_history.tsv"
CURRENT_STAGE=""
CURRENT_STAGE_LABEL=""
CURRENT_STAGE_STARTED_EPOCH=0
CURRENT_STAGE_STARTED_AT=""

format_seconds() {
  local total="$1"
  printf '%02d:%02d:%02d' "$((total / 3600))" "$(((total % 3600) / 60))" "$((total % 60))"
}

write_pipeline_status() {
  local state="$1"
  local elapsed="$2"
  local now
  local temporary="${PIPELINE_STATUS_FILE}.tmp"
  now="$(date '+%Y-%m-%d %H:%M:%S %z')"
  {
    printf 'field\tvalue\n'
    printf 'state\t%s\n' "${state}"
    printf 'stage\t%s\n' "${CURRENT_STAGE}"
    printf 'label\t%s\n' "${CURRENT_STAGE_LABEL}"
    printf 'pid\t%s\n' "$$"
    printf 'started_at\t%s\n' "${CURRENT_STAGE_STARTED_AT}"
    printf 'updated_at\t%s\n' "${now}"
    printf 'elapsed_seconds\t%s\n' "${elapsed}"
  } >"${temporary}"
  mv -- "${temporary}" "${PIPELINE_STATUS_FILE}"
  if [[ "${state}" == "completed" || "${state}" == "failed" ]]; then
    if [[ ! -e "${PIPELINE_HISTORY_FILE}" ]]; then
      printf 'state\tstage\tlabel\tstarted_at\tfinished_at\telapsed_seconds\n' \
        >"${PIPELINE_HISTORY_FILE}"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${state}" "${CURRENT_STAGE}" "${CURRENT_STAGE_LABEL}" \
      "${CURRENT_STAGE_STARTED_AT}" "${now}" "${elapsed}" \
      >>"${PIPELINE_HISTORY_FILE}"
  fi
}

stage_begin() {
  CURRENT_STAGE="$1"
  CURRENT_STAGE_LABEL="$2"
  CURRENT_STAGE_STARTED_EPOCH="$(date +%s)"
  CURRENT_STAGE_STARTED_AT="$(date '+%Y-%m-%d %H:%M:%S %z')"
  write_pipeline_status "running" 0
  echo "[${CURRENT_STAGE}] START ${CURRENT_STAGE_LABEL} at ${CURRENT_STAGE_STARTED_AT}"
}

stage_complete() {
  local elapsed
  elapsed="$(( $(date +%s) - CURRENT_STAGE_STARTED_EPOCH ))"
  write_pipeline_status "completed" "${elapsed}"
  echo "[${CURRENT_STAGE}] DONE ${CURRENT_STAGE_LABEL} in $(format_seconds "${elapsed}")"
  CURRENT_STAGE=""
  CURRENT_STAGE_LABEL=""
}

handle_error() {
  local exit_code="$?"
  trap - ERR
  if [[ -n "${CURRENT_STAGE}" ]]; then
    local elapsed="$(( $(date +%s) - CURRENT_STAGE_STARTED_EPOCH ))"
    write_pipeline_status "failed" "${elapsed}"
    echo "[${CURRENT_STAGE}] FAILED ${CURRENT_STAGE_LABEL} after $(format_seconds "${elapsed}")" >&2
  fi
  exit "${exit_code}"
}
trap handle_error ERR

show_status() {
  echo "Sort-seq project: ${RESULTS_DIR}"
  exec 8>"${LOCK_FILE}"
  if command -v flock >/dev/null 2>&1 && ! flock -n 8; then
    echo "analysis_process: RUNNING"
  else
    echo "analysis_process: IDLE"
  fi
  if [[ -s "${PIPELINE_STATUS_FILE}" ]]; then
    echo
    echo "Latest pipeline stage:"
    column -s $'\t' -t "${PIPELINE_STATUS_FILE}" 2>/dev/null || cat "${PIPELINE_STATUS_FILE}"
  fi
  if [[ -s "${RESULTS_DIR}/index_rescue/rescue_progress.json" ]]; then
    echo
    "${PYTHON_BIN}" - "${RESULTS_DIR}/index_rescue/rescue_progress.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    value = json.load(handle)

eta = value.get("eta_seconds")
eta_text = "unknown" if eta is None else f"{eta / 60:.1f} min"
print("Rescue progress:")
print(f"  status: {value.get('status')}")
print(f"  progress: {value.get('progress_percent', 0):.2f}%")
print(f"  chunk: {value.get('chunk')}/{value.get('chunks_total')}")
print(f"  processed: {value.get('processed_read_pairs_or_reads', 0):,}")
print(f"  rescued so far: {value.get('rescued_read_pairs_or_reads', 0):,} "
      f"({value.get('rescued_percent_so_far', 0):.2f}%)")
print(f"  speed: {value.get('read_pairs_or_reads_per_second', 0):,.0f} read pairs/s")
print(f"  ETA: {eta_text}")
print(f"  updated: {value.get('updated_at')}")
PY
  fi
  echo
  echo "Checkpoint files:"
  for item in \
    "1/5 preflight|${RESULTS_DIR}/index_rescue/rescue_inspection.json" \
    "2/5 rescue|${RESULTS_DIR}/index_rescue/rescue_manifest.json" \
    "3/5 libraryqc|${RESULTS_DIR}/library_qc/combined/variant_count_matrix.csv" \
    "4/5 analyze|${RESULTS_DIR}/sortseq/utr_results_full.tsv" \
    "5/5 plot|${RESULTS_DIR}/sortseq/figures/sortseq_qc_figures.pdf"; do
    local label="${item%%|*}"
    local path="${item#*|}"
    if [[ -s "${path}" ]]; then
      echo "  [DONE] ${label}"
    else
      echo "  [----] ${label}"
    fi
  done
}

if [[ "${MODE}" == "status" ]]; then
  show_status
  exit 0
fi

if [[ "${MODE}" == "resources" ]]; then
  show_resources
  exit 0
fi

exec 9>"${LOCK_FILE}"
if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
  echo "Another Sort-seq analysis is already running for ${RESULTS_DIR}." >&2
  exit 2
fi
write_resource_profile

archive_existing_outputs() {
  local archive_root="${PROJECT_DIR:-$(dirname "${RESULTS_DIR}")}/archive"
  local timestamp
  timestamp="$(date +%Y%m%d_%H%M%S)"
  mkdir -p "${archive_root}"
  for stage in index_rescue library_qc sortseq; do
    if [[ -e "${RESULTS_DIR}/${stage}" ]]; then
      mv -- "${RESULTS_DIR}/${stage}" "${archive_root}/${stage}_${timestamp}"
    fi
  done
}

archive_sortseq_output() {
  local archive_root="${PROJECT_DIR:-$(dirname "${RESULTS_DIR}")}/archive"
  local timestamp
  timestamp="$(date +%Y%m%d_%H%M%S)"
  if [[ -e "${RESULTS_DIR}/sortseq" ]]; then
    mkdir -p "${archive_root}"
    mv -- "${RESULTS_DIR}/sortseq" "${archive_root}/sortseq_before_reanalyze_${timestamp}"
    echo "Archived prior Sort-seq result: ${archive_root}/sortseq_before_reanalyze_${timestamp}"
  fi
}

run_preflight() {
  require_rescue_inputs
  stage_begin "1/5" "Preflight: FASTQ and dual-index inspection"
  run_compute "${PYTHON_BIN}" "${REPO_DIR}/scripts/rescue_dual_index.py" \
    --mode inspect \
    --input-dir "${RAW_DIR}" \
    --sample-sheet "${SAMPLE_SHEET}" \
    --outdir "${RESULTS_DIR}/index_rescue" \
    --max-total-mismatches "${MAX_TOTAL_INDEX_MISMATCHES}" \
    --max-per-index-mismatches "${MAX_PER_INDEX_MISMATCHES}" \
    --minimum-margin "${MIN_INDEX_DISTANCE_MARGIN}" \
    --i5-orientation "${I5_ORIENTATION}"
  stage_complete
}

run_rescue() {
  require_rescue_inputs
  if [[ -d "${RESULTS_DIR}/index_rescue/fastq" ]]; then
    echo "Rescued FASTQ directory already exists: ${RESULTS_DIR}/index_rescue/fastq" >&2
    echo "Use 'full --replace' to archive prior outputs and rerun safely." >&2
    exit 2
  fi
  stage_begin "2/5" "Rescue: uniquely assignable Undetermined reads"
  run_compute "${PYTHON_BIN}" "${REPO_DIR}/scripts/rescue_dual_index.py" \
    --mode rescue \
    --input-dir "${RAW_DIR}" \
    --sample-sheet "${SAMPLE_SHEET}" \
    --outdir "${RESULTS_DIR}/index_rescue" \
    --max-total-mismatches "${MAX_TOTAL_INDEX_MISMATCHES}" \
    --max-per-index-mismatches "${MAX_PER_INDEX_MISMATCHES}" \
    --minimum-margin "${MIN_INDEX_DISTANCE_MARGIN}" \
    --i5-orientation "${I5_ORIENTATION}" \
    --progress-interval-seconds "${PROGRESS_INTERVAL_SECONDS}" \
    --progress-check-reads "${PROGRESS_CHECK_READS}" \
    --compression-backend "${RESCUE_COMPRESSION_BACKEND}" \
    --pigz-threads-per-file "${RESCUE_PIGZ_THREADS_PER_FILE}"
  stage_complete
}

run_libraryqc() {
  : "${LIBRARYQC_PY:?LIBRARYQC_PY is required}"
  : "${LIBRARYQC_CONFIG:?LIBRARYQC_CONFIG is required}"
  require_path "${LIBRARYQC_PY}" "NGS_LibraryQC Python script"
  require_path "${LIBRARYQC_CONFIG}" "NGS_LibraryQC config"
  require_path "${RESULTS_DIR}/index_rescue/fastq" "rescued FASTQ directory"
  if [[ -e "${RESULTS_DIR}/library_qc" ]]; then
    echo "NGS_LibraryQC output already exists: ${RESULTS_DIR}/library_qc" >&2
    echo "Use 'full --replace' to archive prior outputs and rerun safely." >&2
    exit 2
  fi
  stage_begin "3/5" "Library QC: FASTQ QC and UTR counting"
  run_compute "${PYTHON_BIN}" "${LIBRARYQC_PY}" \
    --workers "${WORKERS}" \
    --batch-size "${LIBRARYQC_BATCH_SIZE}" \
    --config "${LIBRARYQC_CONFIG}" \
    --input-dir "${RESULTS_DIR}/index_rescue/fastq" \
    --outdir "${RESULTS_DIR}/library_qc"
  require_path "${RESULTS_DIR}/library_qc/combined/variant_count_matrix.csv" \
    "NGS_LibraryQC variant count matrix"
  stage_complete
}

run_analyze() {
  : "${SAMPLE_MAP:?SAMPLE_MAP is required}"
  require_path "${SAMPLE_MAP}" "sample map"
  local matrix="${RESULTS_DIR}/library_qc/combined/variant_count_matrix.csv"
  require_path "${matrix}" "NGS_LibraryQC variant count matrix"
  if [[ -e "${RESULTS_DIR}/sortseq/utr_results_full.tsv" ]]; then
    echo "Sort-seq result already exists: ${RESULTS_DIR}/sortseq" >&2
    echo "Use 'full --replace' to archive prior outputs and rerun safely." >&2
    exit 2
  fi
  stage_begin "4/5" "Analyze: High15 primary ranking + six-bin supporting score"
  mkdir -p "${RESULTS_DIR}/sortseq/input"
  "${PYTHON_BIN}" "${REPO_DIR}/scripts/prepare_libraryqc_sortseq.py" \
    --matrix "${matrix}" \
    --sample-map "${SAMPLE_MAP}" \
    --outdir "${RESULTS_DIR}/sortseq/input"
  "${PYTHON_BIN}" "${REPO_DIR}/scripts/analyze_sortseq.py" \
    --variants "${RESULTS_DIR}/sortseq/input/variants.tsv" \
    --samples "${RESULTS_DIR}/sortseq/input/samples.tsv" \
    --outdir "${RESULTS_DIR}/sortseq" \
    --overall-gate-fraction "${OVERALL_GATE_FRACTION}" \
    --min-unsorted-count "${MIN_UNSORTED_COUNT}" \
    --min-total-bin-count "${MIN_TOTAL_BIN_COUNT}" \
    --top-hit-min-unsorted-count "${TOP_HIT_MIN_UNSORTED_COUNT}" \
    --top-hit-min-total-bin-count "${TOP_HIT_MIN_TOTAL_BIN_COUNT}" \
    --top-hit-min-high-bin-count "${TOP_HIT_MIN_HIGH_BIN_COUNT}" \
    --top-hit-min-enrichment "${TOP_HIT_MIN_ENRICHMENT}" \
    --bootstrap-replicates "${HIGH15_BOOTSTRAP_REPLICATES}" \
    --bootstrap-seed "${HIGH15_BOOTSTRAP_SEED}" \
    --bootstrap-lower-quantile "${HIGH15_BOOTSTRAP_LOWER_QUANTILE}" \
    --bootstrap-top-n "${HIGH15_BOOTSTRAP_TOP_N}" \
    --bootstrap-min-probability-above-reference "${HIGH15_BOOTSTRAP_MIN_PROBABILITY}" \
    --strict-min-unsorted-count "${STRICT_MIN_UNSORTED_COUNT}" \
    --strict-min-total-bin-count "${STRICT_MIN_TOTAL_BIN_COUNT}" \
    --strict-relative-median-fraction "${STRICT_RELATIVE_MEDIAN_FRACTION}" \
    --strict-min-high-bin-count "${STRICT_MIN_HIGH_BIN_COUNT}" \
    --strict-min-detected-bins "${STRICT_MIN_DETECTED_BINS}" \
    --jackpot-max-bin-probability "${JACKPOT_MAX_BIN_PROBABILITY}" \
    --jackpot-max-detected-bins "${JACKPOT_MAX_DETECTED_BINS}" \
    --reference-variant-id "${REFERENCE_VARIANT_ID}"
  stage_complete
}

run_plot() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript is required for plotting." >&2
    exit 2
  fi
  require_path "${RESULTS_DIR}/sortseq/utr_results_full.tsv" "Sort-seq result table"
  stage_begin "5/5" "Plot: R QC figures and PDF"
  Rscript "${REPO_DIR}/scripts/plot_sortseq.R" \
    "${RESULTS_DIR}" "${RESULTS_DIR}/sortseq/figures"
  if [[ -e "${RESULTS_DIR}/sortseq/metric_comparison/metric_comparison_supported.csv" ]]; then
    Rscript "${REPO_DIR}/scripts/plot_metric_comparison.R" \
      "${RESULTS_DIR}/sortseq/metric_comparison" \
      "${RESULTS_DIR}/sortseq/figures/metric_comparison"
  else
    echo "Metric-comparison CSVs not found; run 'bash run_pipeline.sh compare-metrics' to add figures 21-24."
  fi
  stage_complete
}

run_scoring_steps() {
  : "${SAMPLE_MAP:?SAMPLE_MAP is required}"
  require_path "${SAMPLE_MAP}" "sample map"
  local matrix="${RESULTS_DIR}/library_qc/combined/variant_count_matrix.csv"
  require_path "${matrix}" "NGS_LibraryQC variant count matrix"
  echo "Exporting five scoring steps from: ${matrix}"
  "${PYTHON_BIN}" "${REPO_DIR}/scripts/export_scoring_steps.py" \
    --matrix "${matrix}" \
    --sample-map "${SAMPLE_MAP}" \
    --outdir "${RESULTS_DIR}/sortseq/scoring_steps" \
    --reference-variant-id "${REFERENCE_VARIANT_ID}"
}

run_bimodality_qc() {
  local result_table="${RESULTS_DIR}/sortseq/utr_results_full.tsv"
  require_path "${result_table}" "Sort-seq result table"
  echo "Quantifying high+low-tail polarization and read-count dependence"
  "${PYTHON_BIN}" "${REPO_DIR}/scripts/qc_bimodality.py" \
    --input "${result_table}" \
    --outdir "${RESULTS_DIR}/sortseq/bimodality_qc" \
    --min-total-count "${BIMODALITY_MIN_TOTAL_COUNT}" \
    --min-each-tail-count "${BIMODALITY_MIN_EACH_TAIL_COUNT}" \
    --min-high-tail-probability "${BIMODALITY_MIN_HIGH_TAIL_PROBABILITY}" \
    --min-low-tail-probability "${BIMODALITY_MIN_LOW_TAIL_PROBABILITY}" \
    --max-middle-probability "${BIMODALITY_MAX_MIDDLE_PROBABILITY}" \
    --max-valley-ratio "${BIMODALITY_MAX_VALLEY_RATIO}"
}

run_compare_metrics() {
  local result_table="${RESULTS_DIR}/sortseq/utr_results_full.tsv"
  require_path "${result_table}" "Sort-seq result table"
  echo "Comparing expected score and conditional High15 probability on a shared read-supported universe"
  "${PYTHON_BIN}" "${REPO_DIR}/scripts/compare_metrics.py" \
    --input "${result_table}" \
    --outdir "${RESULTS_DIR}/sortseq/metric_comparison" \
    --min-total-count "${METRIC_COMPARE_MIN_TOTAL_COUNT}" \
    --min-unsorted-count "${METRIC_COMPARE_MIN_UNSORTED_COUNT}" \
    --min-high-bin-count "${METRIC_COMPARE_MIN_HIGH_BIN_COUNT}" \
    --top-n "${METRIC_COMPARE_TOP_N}"
}

if [[ "${MODE}" == "full" && "${OPTION}" == "--replace" ]]; then
  archive_existing_outputs
fi

case "${MODE}" in
  preflight) run_preflight ;;
  rescue) run_rescue ;;
  libraryqc) run_libraryqc ;;
  analyze) run_analyze ;;
  scoring-steps) run_scoring_steps ;;
  bimodality-qc) run_bimodality_qc ;;
  compare-metrics) run_compare_metrics ;;
  plot) run_plot ;;
  status) show_status ;;
  resources) show_resources ;;
  reanalyze-analysis)
    archive_sortseq_output
    run_analyze
    echo "Python analysis completed; run 'bash run_pipeline.sh plot' in the R environment."
    ;;
  reanalyze)
    archive_sortseq_output
    run_analyze
    run_plot
    echo "Reference-aware analysis completed: ${RESULTS_DIR}/sortseq"
    ;;
  full)
    run_preflight
    run_rescue
    run_libraryqc
    run_analyze
    run_plot
    echo "Completed: ${RESULTS_DIR}"
    ;;
esac
