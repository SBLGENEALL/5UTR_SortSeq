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
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh plot
  SORTSEQ_PROJECT_CONFIG=/path/to/sortseq.env bash run_pipeline.sh full [--replace]
EOF
}

if [[ ! "${MODE}" =~ ^(preflight|rescue|libraryqc|analyze|plot|full)$ ]]; then
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
MAX_TOTAL_INDEX_MISMATCHES="${MAX_TOTAL_INDEX_MISMATCHES:-2}"
MAX_PER_INDEX_MISMATCHES="${MAX_PER_INDEX_MISMATCHES:-2}"
MIN_INDEX_DISTANCE_MARGIN="${MIN_INDEX_DISTANCE_MARGIN:-1}"
I5_ORIENTATION="${I5_ORIENTATION:-auto}"
OVERALL_GATE_FRACTION="${OVERALL_GATE_FRACTION:-0.90}"
MIN_UNSORTED_COUNT="${MIN_UNSORTED_COUNT:-50}"
MIN_TOTAL_BIN_COUNT="${MIN_TOTAL_BIN_COUNT:-100}"

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

mkdir -p "${RESULTS_DIR}"
LOCK_FILE="${RESULTS_DIR}/.sortseq_analysis.lock"
exec 9>"${LOCK_FILE}"
if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
  echo "Another Sort-seq analysis is already running for ${RESULTS_DIR}." >&2
  exit 2
fi

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

run_preflight() {
  require_rescue_inputs
  echo "[1/5] Inspecting FASTQ names, header dual indexes, and i5 orientation"
  "${PYTHON_BIN}" "${REPO_DIR}/scripts/rescue_dual_index.py" \
    --mode inspect \
    --input-dir "${RAW_DIR}" \
    --sample-sheet "${SAMPLE_SHEET}" \
    --outdir "${RESULTS_DIR}/index_rescue" \
    --max-total-mismatches "${MAX_TOTAL_INDEX_MISMATCHES}" \
    --max-per-index-mismatches "${MAX_PER_INDEX_MISMATCHES}" \
    --minimum-margin "${MIN_INDEX_DISTANCE_MARGIN}" \
    --i5-orientation "${I5_ORIENTATION}"
}

run_rescue() {
  require_rescue_inputs
  if [[ -d "${RESULTS_DIR}/index_rescue/fastq" ]]; then
    echo "Rescued FASTQ directory already exists: ${RESULTS_DIR}/index_rescue/fastq" >&2
    echo "Use 'full --replace' to archive prior outputs and rerun safely." >&2
    exit 2
  fi
  echo "[2/5] Rescuing uniquely assignable Undetermined reads"
  "${PYTHON_BIN}" "${REPO_DIR}/scripts/rescue_dual_index.py" \
    --mode rescue \
    --input-dir "${RAW_DIR}" \
    --sample-sheet "${SAMPLE_SHEET}" \
    --outdir "${RESULTS_DIR}/index_rescue" \
    --max-total-mismatches "${MAX_TOTAL_INDEX_MISMATCHES}" \
    --max-per-index-mismatches "${MAX_PER_INDEX_MISMATCHES}" \
    --minimum-margin "${MIN_INDEX_DISTANCE_MARGIN}" \
    --i5-orientation "${I5_ORIENTATION}"
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
  echo "[3/5] Running NGS_LibraryQC on assigned plus rescued FASTQs"
  "${PYTHON_BIN}" "${LIBRARYQC_PY}" \
    --workers "${WORKERS}" \
    --config "${LIBRARYQC_CONFIG}" \
    --input-dir "${RESULTS_DIR}/index_rescue/fastq" \
    --outdir "${RESULTS_DIR}/library_qc"
  require_path "${RESULTS_DIR}/library_qc/combined/variant_count_matrix.csv" \
    "NGS_LibraryQC variant count matrix"
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
  echo "[4/5] Correcting bin sizes and calculating UTR scores"
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
    --min-total-bin-count "${MIN_TOTAL_BIN_COUNT}"
}

run_plot() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript is required for plotting." >&2
    exit 2
  fi
  require_path "${RESULTS_DIR}/sortseq/utr_results_full.tsv" "Sort-seq result table"
  echo "[5/5] Creating R figures and a multi-page PDF"
  Rscript "${REPO_DIR}/scripts/plot_sortseq.R" \
    "${RESULTS_DIR}" "${RESULTS_DIR}/sortseq/figures"
}

if [[ "${MODE}" == "full" && "${OPTION}" == "--replace" ]]; then
  archive_existing_outputs
fi

case "${MODE}" in
  preflight) run_preflight ;;
  rescue) run_rescue ;;
  libraryqc) run_libraryqc ;;
  analyze) run_analyze ;;
  plot) run_plot ;;
  full)
    run_preflight
    run_rescue
    run_libraryqc
    run_analyze
    run_plot
    echo "Completed: ${RESULTS_DIR}"
    ;;
esac
