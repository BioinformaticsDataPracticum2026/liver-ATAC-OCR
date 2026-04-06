#!/usr/bin/bash
#SBATCH -J motif_analysis     # Job name
#SBATCH -p RM-shared          # Partition
#SBATCH -N 1                  # Number of nodes
#SBATCH --cpus-per-task=16     # Number of tasks (CPUs)
#SBATCH -t 12:00:00           # Walltime (hh:mm:ss)
#SBATCH --mem=30G             # Memory
#SBATCH -o /ocean/projects/bio230007p/sujathab/Motif_analysis/logs/task5_motif_%j.out
#SBATCH -e /ocean/projects/bio230007p/sujathab/Motif_analysis/logs/task5_motif_%j.err


set -euo pipefail

# MODULE LOADING
# Check available versions with: module spider bedtools && module spider MEME-suite

module purge
module load bedtools
module load MEME-suite/5.4.1

# PATHS — update MM10_GENOME and JASPAR_DB before running

ROOT="/ocean/projects/bio230007p/sujathab"
PE_OUT="$ROOT/PE_classification/output"
ROWCOUNT="$PE_OUT/rowcount"
UNIQUE="$PE_OUT/unique"

OUT="$ROOT/task5_motifs"
BEDS="$OUT/beds"
FASTAS="$OUT/fastas"

# Reference genome (mm10). All coordinates are in mouse space (HALPER output).
MM10_GENOME="/ocean/projects/bio230007p/sujathab/mm10.fa"

# CIS-BP motif database (MEME format).
CIS_BP_DB="/ocean/projects/bio230007p/ikaplow/CIS-BP_2.00/Mus_musculus.meme"

# Fixed window size for motif analysis (bp centered on peak midpoint)
WINDOW=100   # ±100bp = 200bp total; TF binding sites are typically 8–20bp

# SETUP

mkdir -p "$BEDS" "$FASTAS" "$OUT/logs" logs

echo "Performing Motif Analysis"
echo "Root:       $ROOT"
echo "Output:     $OUT"
echo "Genome:     $MM10_GENOME"
echo "CIS-BP:     $CIS_BP_DB"
echo "Window:     ±${WINDOW}bp"
echo ""

# Quick sanity checks
if [[ ! -f "$MM10_GENOME" ]]; then
  echo "ERROR: mm10 genome not found at $MM10_GENOME"
  echo "Update MM10_GENOME in this script and resubmit."
  exit 1
fi
if [[ ! -f "$CIS_BP_DB" ]]; then
  echo "ERROR: CIS-BP DB not found at $CIS_BP_DB"
  echo "Update CIS_BP_DB in this script and resubmit."
  exit 1
fi

# STEP 1: Generate missing mouse BED files
# Task 4 only classified human-projected OCRs. We need mouse-native
# enhancer/promoter calls and mouse-specific enhancers separately.

echo "[Step 1] Generating mouse-native and mouse-specific BED files..."

# Mouse-native enhancers: mouse peaks NOT near any TSS
bedtools intersect -v \
  -a "$ROWCOUNT/mouse_native.sorted.bed" \
  -b "$ROWCOUNT/mouse_promoter_pm2kb.merged.bed" \
  | cut -f1-3 | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/mouse_native_enhancer.bed"

# Mouse-native promoters: mouse peaks overlapping TSS ±2kb
bedtools intersect -u \
  -a "$ROWCOUNT/mouse_native.sorted.bed" \
  -b "$ROWCOUNT/mouse_promoter_pm2kb.merged.bed" \
  | cut -f1-3 | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/mouse_native_promoter.bed"

# Mouse-specific open = mouse peaks with no human ortholog (via HALPER)
bedtools intersect -v \
  -a "$ROWCOUNT/mouse_native.sorted.bed" \
  -b "$ROWCOUNT/halper_mouse.sorted.bed" \
  | cut -f1-3 | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/mouse_specific_open.bed"

# Split mouse-specific into enhancer and promoter
bedtools intersect -v \
  -a "$BEDS/mouse_specific_open.bed" \
  -b "$ROWCOUNT/mouse_promoter_pm2kb.merged.bed" \
  > "$BEDS/mouse_specific_enhancer.bed"

bedtools intersect -u \
  -a "$BEDS/mouse_specific_open.bed" \
  -b "$ROWCOUNT/mouse_promoter_pm2kb.merged.bed" \
  > "$BEDS/mouse_specific_promoter.bed"

# Human all-enhancers = shared + human-specific (all in mm10 coords via HALPER)
cat "$UNIQUE/shared_enhancer.bed" "$UNIQUE/human_specific_enhancer.bed" \
  | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/human_all_enhancer.bed"

# Human all-promoters = shared + human-specific
cat "$UNIQUE/shared_promoter.bed" "$UNIQUE/human_specific_promoter.bed" \
  | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/human_all_promoter.bed"

# Copy shared and species-specific enhancer beds for clarity
cp "$UNIQUE/shared_enhancer.bed"          "$BEDS/shared_enhancer.bed"
cp "$UNIQUE/human_specific_enhancer.bed"  "$BEDS/human_specific_enhancer.bed"

echo "[Step 1] Region counts:"
for f in "$BEDS"/*.bed; do
  printf "  %-35s %d regions\n" "$(basename "$f")" "$(wc -l < "$f")"
done
echo ""

# STEP 2: Resize peaks to fixed window centered on midpoint
# MEME-ChIP works best with uniform-length sequences. We take the midpoint of
# each peak and extend ±WINDOW bp. This also reduces noise from long peaks.

echo "[Step 2] Resizing peaks to $(( WINDOW * 2 ))bp windows..."

resize_bed() {
  local in_bed="$1"
  local out_bed="$2"
  awk -v w="$WINDOW" 'BEGIN{OFS="\t"} {
    mid = int(($2 + $3) / 2)
    s = mid - w
    if (s < 0) s = 0
    print $1, s, mid + w
  }' "$in_bed" | sort -k1,1 -k2,2n | uniq > "$out_bed"
}

declare -A PEAK_SETS
PEAK_SETS["human_enhancer"]="$BEDS/human_all_enhancer.bed"
PEAK_SETS["mouse_enhancer"]="$BEDS/mouse_native_enhancer.bed"
PEAK_SETS["human_promoter"]="$BEDS/human_all_promoter.bed"
PEAK_SETS["mouse_promoter"]="$BEDS/mouse_native_promoter.bed"
PEAK_SETS["shared_enhancer"]="$BEDS/shared_enhancer.bed"
PEAK_SETS["human_specific_enhancer"]="$BEDS/human_specific_enhancer.bed"
PEAK_SETS["mouse_specific_enhancer"]="$BEDS/mouse_specific_enhancer.bed"

for name in "${!PEAK_SETS[@]}"; do
  resize_bed "${PEAK_SETS[$name]}" "$BEDS/${name}.resized.bed"
  echo "  $name: $(wc -l < "$BEDS/${name}.resized.bed") windows"
done
echo ""

# STEP 3: Extract FASTA sequences from mm10 genome

echo "[Step 3] Extracting FASTA sequences from mm10..."

for name in "${!PEAK_SETS[@]}"; do
  bedtools getfasta \
    -fi "$MM10_GENOME" \
    -bed "$BEDS/${name}.resized.bed" \
    -fo "$FASTAS/${name}.fa"
  n_seqs=$(grep -c '^>' "$FASTAS/${name}.fa")
  echo "  $name: $n_seqs sequences"
done
echo ""

# STEP 4: Run MEME-ChIP on each set
# MEME-ChIP automatically subsamples large inputs for the MEME step (~600 seqs),
# then runs DREME (fast de novo) on the full set, plus TOMTOM against CIS-BP.

echo "[Step 4] Running MEME-ChIP (this will take a while)..."

run_meme_chip() {
  local name="$1"
  local fasta="$FASTAS/${name}.fa"
  local outdir="$OUT/meme_chip_${name}"
  local logfile="$OUT/logs/${name}.log"

  echo "  [$name] Starting..."
  meme-chip \
    -oc "$outdir" \
    -maxw 20 \
    -meme-nmotifs 10 \
    -meme-p "${SLURM_CPUS_PER_TASK:-8}" \
    -db "$CIS_BP_DB" \
    "$fasta" \
    > "$logfile" 2>&1
  echo "  [$name] Done. Results: $outdir"
}

# Run sequentially (MEME-ChIP uses -meme-p threads per job internally)
# To parallelise across sets, wrap each call in & and add wait at the end —
# but be mindful of total memory with 8 sets running concurrently.
for name in "${!PEAK_SETS[@]}"; do
  run_meme_chip "$name"
done

# STEP 5: Print summary of top motif hits

echo ""
echo "[Step 5] Top motif hits per set:"
for name in "${!PEAK_SETS[@]}"; do
  summary="$OUT/meme_chip_${name}/summary.tsv"
  if [[ -f "$summary" ]]; then
    echo "  === $name ==="
    head -4 "$summary"
    echo ""
  fi
done

echo "===== All done ====="
echo "Results in: $OUT"
echo "Open meme_chip_<name>/meme-chip.html in a browser for full results."