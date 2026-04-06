#!/usr/bin/env bash
#SBATCH --job-name=pe_classify
#SBATCH --output=slurm-pe_classify-%j.out
#SBATCH --error=slurm-pe_classify-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=bio230007p

set -euo pipefail

TOTAL_STEPS=7
CURRENT_STEP=0

print_progress() {
  local step="$1"
  local total="$2"
  local msg="$3"
  local width=30
  local filled=$(( step * width / total ))
  local empty=$(( width - filled ))

  printf "\n["
  printf "%${filled}s" "" | tr ' ' '#'
  printf "%${empty}s" "" | tr ' ' '-'
  printf "] %d/%d %s\n" "$step" "$total" "$msg"
}

next_step() {
  CURRENT_STEP=$((CURRENT_STEP + 1))
  print_progress "$CURRENT_STEP" "$TOTAL_STEPS" "$1"
}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IN="$ROOT/PE_classification/input"
OUT_ROOT="$ROOT/PE_classification/output"
ROW_OUT="$OUT_ROOT/rowcount"
UNIQ_OUT="$OUT_ROOT/unique"

HALPER="$IN/human_liver.HumanToMouse.HALPER.narrowPeak.gz"
MOUSE_PEAK="$IN/idr.conservative_peak.narrowPeak.gz"
TSS="$IN/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed"

mkdir -p "$ROW_OUT" "$UNIQ_OUT"

if ! command -v module >/dev/null 2>&1; then
  if [[ -f /etc/profile.d/modules.sh ]]; then
    source /etc/profile.d/modules.sh
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    source /usr/share/Modules/init/bash
  fi
fi

if command -v module >/dev/null 2>&1; then
  module load bedtools/2.31.1 2>/dev/null || module load bedtools 2>/dev/null || true
fi

if ! command -v bedtools >/dev/null 2>&1; then
  echo "ERROR: bedtools not found. On Bridges-2, run with: module load bedtools"
  exit 1
fi

for f in "$HALPER" "$MOUSE_PEAK" "$TSS"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing input file: $f"
    exit 1
  fi
done

dedup_by_coord() {
  local in_bed="$1"
  local out_bed="$2"
  cut -f1-3 "$in_bed" | sort -k1,1 -k2,2n -k3,3n | uniq > "$out_bed"
}

echo "Input directory: $IN"
echo "Output directory: $OUT_ROOT"

next_step "Step 1: sort HALPER output and mouse native peaks ..."
echo "Running: gzip + sort for HALPER and mouse native peaks"
gzip -dc "$HALPER" | sort -k1,1 -k2,2n > "$ROW_OUT/halper_mouse.sorted.bed"
gzip -dc "$MOUSE_PEAK" | sort -k1,1 -k2,2n > "$ROW_OUT/mouse_native.sorted.bed"
echo "Finished Step 1."

next_step "Step 2: create promoter regions from TSS (TSS ± 2 kb) ..."
echo "Running: awk + sort + bedtools merge"
awk 'BEGIN{OFS="\t"} {s=$2-2000; if(s<0)s=0; e=$3+2000; print $1,s,e}' "$TSS" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "$ROW_OUT/mouse_promoter_pm2kb.merged.bed"
echo "Finished Step 2."

next_step "Step 3: identify shared-open OCRs and human-open/mouse-closed OCRs ..."
echo "Running: bedtools intersect for shared_open"
bedtools intersect -sorted -u \
  -a "$ROW_OUT/halper_mouse.sorted.bed" \
  -b "$ROW_OUT/mouse_native.sorted.bed" \
  > "$ROW_OUT/shared_open.bed"

echo "Running: bedtools intersect for human_open_mouse_closed"
bedtools intersect -sorted -v \
  -a "$ROW_OUT/halper_mouse.sorted.bed" \
  -b "$ROW_OUT/mouse_native.sorted.bed" \
  > "$ROW_OUT/human_open_mouse_closed.bed"
echo "Finished Step 3."

next_step "Step 4: classify shared OCRs into promoter-like vs enhancer-like ..."
echo "Running: bedtools intersect for shared_promoter"
bedtools intersect -sorted -u \
  -a "$ROW_OUT/shared_open.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/shared_promoter.bed"

echo "Running: bedtools intersect for shared_enhancer"
bedtools intersect -sorted -v \
  -a "$ROW_OUT/shared_open.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/shared_enhancer.bed"
echo "Finished Step 4."

next_step "Step 5: classify human-specific OCRs into promoter-like vs enhancer-like ..."
echo "Running: bedtools intersect for human_specific_promoter"
bedtools intersect -sorted -u \
  -a "$ROW_OUT/human_open_mouse_closed.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/human_specific_promoter.bed"

echo "Running: bedtools intersect for human_specific_enhancer"
bedtools intersect -sorted -v \
  -a "$ROW_OUT/human_open_mouse_closed.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/human_specific_enhancer.bed"
echo "Finished Step 5."

next_step "Step 6: summarize row-based counts (reference) ..."
shared_total=$(wc -l < "$ROW_OUT/shared_open.bed")
shared_prom=$(wc -l < "$ROW_OUT/shared_promoter.bed")
shared_enh=$(wc -l < "$ROW_OUT/shared_enhancer.bed")

human_spec_total=$(wc -l < "$ROW_OUT/human_open_mouse_closed.bed")
human_spec_prom=$(wc -l < "$ROW_OUT/human_specific_promoter.bed")
human_spec_enh=$(wc -l < "$ROW_OUT/human_specific_enhancer.bed")

shared_prom_frac=$(awk -v p="$shared_prom" -v t="$shared_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
shared_enh_frac=$(awk -v e="$shared_enh" -v t="$shared_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

human_prom_frac=$(awk -v p="$human_spec_prom" -v t="$human_spec_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
human_enh_frac=$(awk -v e="$human_spec_enh" -v t="$human_spec_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

{
  echo -e "category\ttotal\tpromoter\tenhancer\tpromoter_fraction\tenhancer_fraction"
  echo -e "shared_open\t$shared_total\t$shared_prom\t$shared_enh\t$shared_prom_frac\t$shared_enh_frac"
  echo -e "human_open_mouse_closed\t$human_spec_total\t$human_spec_prom\t$human_spec_enh\t$human_prom_frac\t$human_enh_frac"
} > "$ROW_OUT/summary.tsv"

echo "Row-based summary:"
cat "$ROW_OUT/summary.tsv"
echo "Finished Step 6."

next_step "Step 7: deduplicate BED files by genomic coordinates and summarize unique OCR counts ..."
dedup_by_coord "$ROW_OUT/shared_open.bed" "$UNIQ_OUT/shared_open.bed"
dedup_by_coord "$ROW_OUT/human_open_mouse_closed.bed" "$UNIQ_OUT/human_open_mouse_closed.bed"
dedup_by_coord "$ROW_OUT/shared_promoter.bed" "$UNIQ_OUT/shared_promoter.bed"
dedup_by_coord "$ROW_OUT/shared_enhancer.bed" "$UNIQ_OUT/shared_enhancer.bed"
dedup_by_coord "$ROW_OUT/human_specific_promoter.bed" "$UNIQ_OUT/human_specific_promoter.bed"
dedup_by_coord "$ROW_OUT/human_specific_enhancer.bed" "$UNIQ_OUT/human_specific_enhancer.bed"

shared_total_u=$(wc -l < "$UNIQ_OUT/shared_open.bed")
shared_prom_u=$(wc -l < "$UNIQ_OUT/shared_promoter.bed")
shared_enh_u=$(wc -l < "$UNIQ_OUT/shared_enhancer.bed")

human_spec_total_u=$(wc -l < "$UNIQ_OUT/human_open_mouse_closed.bed")
human_spec_prom_u=$(wc -l < "$UNIQ_OUT/human_specific_promoter.bed")
human_spec_enh_u=$(wc -l < "$UNIQ_OUT/human_specific_enhancer.bed")

shared_prom_frac_u=$(awk -v p="$shared_prom_u" -v t="$shared_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
shared_enh_frac_u=$(awk -v e="$shared_enh_u" -v t="$shared_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

human_prom_frac_u=$(awk -v p="$human_spec_prom_u" -v t="$human_spec_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
human_enh_frac_u=$(awk -v e="$human_spec_enh_u" -v t="$human_spec_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

{
  echo -e "category\ttotal\tpromoter\tenhancer\tpromoter_fraction\tenhancer_fraction"
  echo -e "shared_open\t$shared_total_u\t$shared_prom_u\t$shared_enh_u\t$shared_prom_frac_u\t$shared_enh_frac_u"
  echo -e "human_open_mouse_closed\t$human_spec_total_u\t$human_spec_prom_u\t$human_spec_enh_u\t$human_prom_frac_u\t$human_enh_frac_u"
} > "$UNIQ_OUT/summary.tsv"

echo "Unique-based summary:"
cat "$UNIQ_OUT/summary.tsv"
echo "Finished Step 7."

echo
echo "[##############################] 7/7 Pipeline completed"
echo "Done."
echo "Row-based output folder: $ROW_OUT"
echo "Unique-based output folder: $UNIQ_OUT"