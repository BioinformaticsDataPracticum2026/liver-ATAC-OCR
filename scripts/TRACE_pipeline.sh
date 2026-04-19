#!/usr/bin/env bash
# =============================================================================
# TRACE_pipeline.sh — liver-ATAC-OCR Master Pipeline
# 03-713: Bioinformatics Data Integration Practicum, Spring 2026
#
# Runs Tasks 2–5 in sequence:
#   Step 1: HALPER cross-species liftover          (run_halper_mapping.sh)
#   Step 2: Promoter/Enhancer classification       (run_pe_classification.sh)
#   Step 3: rGREAT biological process enrichment   (run_rgreat.sh → rgreat_analysis.R)
#   Step 4: rGREAT plots                           (run_plots.sh  → rgreat_plots.R)
#   Step 5: Motif analysis (MEME-ChIP)             (run_motif_analysis.sh)
#
# Usage:
#   bash TRACE_pipeline.sh [OPTIONS]
#
# Options:
#   --root DIR          Project root directory (default: repo root, auto-detected)
#   --human FILE        Path to human ATAC-seq narrowPeak file (.gz)
#   --mouse FILE        Path to mouse ATAC-seq narrowPeak file (.gz)
#   --hal FILE          Path to HAL alignment file (.hal)
#   --tss FILE          Path to TSS annotation BED file
#   --genome FILE       Path to mm10 genome FASTA file
#   --jaspar FILE       Path to JASPAR motif database (.meme)
#   --source-species S  Source species for HALPER (default: Human)
#   --target-species T  Target species for HALPER (default: Mouse)
#   --skip-halper       Skip Step 1 (if liftover already done)
#   --skip-pe           Skip Step 2 (if PE classification already done)
#   --skip-great        Skip Steps 3–4 (if rGREAT already done)
#   --skip-motif        Skip Step 5
#   --conda-env ENV     Conda environment name for rGREAT (default: rgreat_env)
#   -h, --help          Show this help message
#
# Dependencies: bedtools, MEME-suite, R (rGREAT, ggplot2), conda/hal env
# Designed to run on Bridges-2 (PSC) or any Linux cluster with modules.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

HUMAN_PEAKS=""
MOUSE_PEAKS=""
HAL_FILE=""
TSS_FILE=""
MM10_GENOME=""
JASPAR_DB=""
SOURCE_SPECIES="Human"
TARGET_SPECIES="Mouse"
CONDA_ENV="rgreat_env"

SKIP_HALPER=false
SKIP_PE=false
SKIP_GREAT=false
SKIP_MOTIF=false

# -----------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $*" >&2; }
die()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; exit 1; }

step_banner() {
  echo ""
  echo "============================================================"
  echo "  $*"
  echo "============================================================"
}

# -----------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------
usage() {
  sed -n '/^# Usage:/,/^# Dependencies:/p' "$0" | sed 's/^# \?//'
  exit 0
}

# -----------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --root)            ROOT="$2";           shift 2 ;;
    --human)           HUMAN_PEAKS="$2";    shift 2 ;;
    --mouse)           MOUSE_PEAKS="$2";    shift 2 ;;
    --hal)             HAL_FILE="$2";       shift 2 ;;
    --tss)             TSS_FILE="$2";       shift 2 ;;
    --genome)          MM10_GENOME="$2";    shift 2 ;;
    --jaspar)          JASPAR_DB="$2";      shift 2 ;;
    --source-species)  SOURCE_SPECIES="$2"; shift 2 ;;
    --target-species)  TARGET_SPECIES="$2"; shift 2 ;;
    --conda-env)       CONDA_ENV="$2";      shift 2 ;;
    --skip-halper)     SKIP_HALPER=true;    shift ;;
    --skip-pe)         SKIP_PE=true;        shift ;;
    --skip-great)      SKIP_GREAT=true;     shift ;;
    --skip-motif)      SKIP_MOTIF=true;     shift ;;
    -h|--help)         usage ;;
    *) die "Unknown option: $1. Use -h for help." ;;
  esac
done

# -----------------------------------------------------------------------
# Resolve default file paths (from repo layout) if not passed explicitly
# -----------------------------------------------------------------------
DATA_DIR="$ROOT/data/raw"
HALPER_OUT="$ROOT/Mapping/outputs"
PE_INPUT="$ROOT/PE_classification/input"

HUMAN_PEAKS="${HUMAN_PEAKS:-$DATA_DIR/human_liver.narrowPeak.gz}"
MOUSE_PEAKS="${MOUSE_PEAKS:-$DATA_DIR/mouse_liver.narrowPeak.gz}"
HAL_FILE="${HAL_FILE:-$DATA_DIR/10plusway-master.hal}"
TSS_FILE="${TSS_FILE:-$PE_INPUT/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed}"
MM10_GENOME="${MM10_GENOME:-$DATA_DIR/mm10.fa}"
JASPAR_DB="${JASPAR_DB:-$DATA_DIR/motif_dbs/JASPAR2026_vertebrates_combined.meme}"

SCRIPTS_DIR="$ROOT/scripts"
RGREAT_SCRIPTS="$ROOT/rGREAT_Analysis/scripts"

LOG_DIR="$ROOT/logs/pipeline_$(date '+%Y%m%d_%H%M%S')"
mkdir -p "$LOG_DIR"
PIPELINE_LOG="$LOG_DIR/pipeline.log"

# Tee all output to log file
exec > >(tee -a "$PIPELINE_LOG") 2>&1

# -----------------------------------------------------------------------
# Environment / module helpers
# -----------------------------------------------------------------------
load_modules() {
  if command -v module >/dev/null 2>&1; then
    module purge
    module load bedtools   2>/dev/null || warn "Could not load bedtools module — assuming it is in PATH"
    module load MEME-suite 2>/dev/null || warn "Could not load MEME-suite module — assuming it is in PATH"
    module load anaconda3  2>/dev/null || warn "Could not load anaconda3 module — assuming conda is in PATH"
  else
    warn "module command not found — assuming all tools are already in PATH"
  fi
}

activate_conda() {
  local env="$1"
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate "$env" || warn "Could not activate conda env '$env' — proceeding anyway"
  fi
}

# -----------------------------------------------------------------------
# Pre-flight checks
# -----------------------------------------------------------------------
preflight_checks() {
  log "Running pre-flight checks..."

  local missing=0

  check_file() {
    local label="$1" path="$2" required="$3"
    if [[ ! -f "$path" ]]; then
      if [[ "$required" == "required" ]]; then
        warn "Missing required file [$label]: $path"
        missing=$((missing + 1))
      else
        warn "Optional file not found [$label]: $path"
      fi
    else
      log "  Found [$label]: $path"
    fi
  }

  $SKIP_HALPER || {
    check_file "human peaks"  "$HUMAN_PEAKS" required
    check_file "mouse peaks"  "$MOUSE_PEAKS" required
    check_file "HAL file"     "$HAL_FILE"    required
  }

  $SKIP_PE || {
    check_file "TSS annotation" "$TSS_FILE" required
  }

  $SKIP_MOTIF || {
    check_file "mm10 genome"  "$MM10_GENOME" required
    check_file "JASPAR DB"    "$JASPAR_DB"   required
  }

  if [[ $missing -gt 0 ]]; then
    die "$missing required input file(s) not found. Check paths above or pass them via flags."
  fi

  # Check tools
  for tool in bedtools meme-chip Rscript; do
    command -v "$tool" >/dev/null 2>&1 \
      && log "  Tool found: $tool" \
      || warn "  Tool not found in PATH: $tool (may be loaded as a module later)"
  done

  log "Pre-flight checks complete."
}

# -----------------------------------------------------------------------
# Step 1: HALPER Liftover
# -----------------------------------------------------------------------
run_halper() {
  step_banner "Step 1 / 5 — HALPER Cross-Species Liftover"
  log "Script: $SCRIPTS_DIR/run_halper_mapping.sh"

  # Override hardcoded paths inside the script via exports
  export HUMAN_PEAKS MOUSE_PEAKS HAL_FILE SOURCE_SPECIES TARGET_SPECIES
  export OUTPUT_DIR="$HALPER_OUT"
  export HALPER_SCRIPT="$HOME/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh"

  mkdir -p "$HALPER_OUT" "$LOG_DIR"

  activate_conda hal
  export PATH="$HOME/repos/hal/bin:$PATH"
  export PYTHONPATH="$HOME/repos/halLiftover-postprocessing:${PYTHONPATH:-}"

  bash "$SCRIPTS_DIR/run_halper_mapping.sh" \
    > "$LOG_DIR/01_halper.log" 2>&1 \
    && log "Step 1 complete. Output: $HALPER_OUT" \
    || die "Step 1 failed. See $LOG_DIR/01_halper.log"
}

# -----------------------------------------------------------------------
# Step 2: P/E Classification
# -----------------------------------------------------------------------
run_pe_classification() {
  step_banner "Step 2 / 5 — Promoter / Enhancer Classification"
  log "Script: $SCRIPTS_DIR/run_pe_classification.sh"

  bash "$SCRIPTS_DIR/run_pe_classification.sh" \
    > "$LOG_DIR/02_pe_classification.log" 2>&1 \
    && log "Step 2 complete. Summaries in $ROOT/PE_classification/output/" \
    || die "Step 2 failed. See $LOG_DIR/02_pe_classification.log"
}

# -----------------------------------------------------------------------
# Step 3: rGREAT Enrichment
# -----------------------------------------------------------------------
run_rgreat() {
  step_banner "Step 3 / 5 — rGREAT Biological Process Enrichment"
  log "Script: $RGREAT_SCRIPTS/rgreat_analysis.R"

  activate_conda "$CONDA_ENV"

  cd "$ROOT"
  Rscript "$RGREAT_SCRIPTS/rgreat_analysis.R" \
    > "$LOG_DIR/03_rgreat.log" 2>&1 \
    && log "Step 3 complete. CSVs in $ROOT/rGREAT_Analysis/outputs/" \
    || die "Step 3 failed. See $LOG_DIR/03_rgreat.log"
}

# -----------------------------------------------------------------------
# Step 4: rGREAT Plots
# -----------------------------------------------------------------------
run_rgreat_plots() {
  step_banner "Step 4 / 5 — rGREAT Plots"
  log "Script: $RGREAT_SCRIPTS/rgreat_plots.R"

  activate_conda "$CONDA_ENV"

  cd "$ROOT"
  Rscript "$RGREAT_SCRIPTS/rgreat_plots.R" \
    > "$LOG_DIR/04_rgreat_plots.log" 2>&1 \
    && log "Step 4 complete. Plots in $ROOT/rGREAT_Analysis/outputs/plots/" \
    || die "Step 4 failed. See $LOG_DIR/04_rgreat_plots.log"
}

# -----------------------------------------------------------------------
# Step 5: Motif Analysis (MEME-ChIP)
# -----------------------------------------------------------------------
run_motif_analysis() {
  step_banner "Step 5 / 5 — Motif Analysis (MEME-ChIP)"
  log "Script: $SCRIPTS_DIR/run_motif_analysis.sh"

  export MM10_GENOME JASPAR_DB

  bash "$SCRIPTS_DIR/run_motif_analysis.sh" \
    > "$LOG_DIR/05_motif_analysis.log" 2>&1 \
    && log "Step 5 complete. Results in $ROOT/Motif_analysis/output/" \
    || die "Step 5 failed. See $LOG_DIR/05_motif_analysis.log"
}

# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------
log "======================================================"
log "  liver-ATAC-OCR Pipeline — $(date)"
log "  Project root: $ROOT"
log "  Log directory: $LOG_DIR"
log "======================================================"

load_modules
preflight_checks

$SKIP_HALPER && log "Skipping Step 1 (--skip-halper)" || run_halper
$SKIP_PE     && log "Skipping Step 2 (--skip-pe)"     || run_pe_classification
$SKIP_GREAT  && log "Skipping Steps 3–4 (--skip-great)" || { run_rgreat; run_rgreat_plots; }
$SKIP_MOTIF  && log "Skipping Step 5 (--skip-motif)"  || run_motif_analysis

echo ""
log "======================================================"
log "  Pipeline complete!"
log "  Full log: $PIPELINE_LOG"
log "======================================================"