#!/bin/bash
#SBATCH -J halper_map       # Job name
#SBATCH -p RM-shared        # Partition
#SBATCH -N 1                # Number of nodes
#SBATCH -n 1                # Number of tasks (CPUs)
#SBATCH -t 02:00:00         # Walltime (hh:mm:ss)
#SBATCH --mem=16G           # Memory
#SBATCH -o logs/halper_%j.out
#SBATCH -e logs/halper_%j.err

# -------------------------------
# HALPER mapping sbatch script
# -------------------------------

# Load modules
module load anaconda3     

# Activate HAL conda environment
source activate hal

# Set file paths
HUMAN_PEAKS="/ocean/projects/bio230007p/ikaplow/HumanAtac/Liver/peak/rep1/SRR13439659_1.trim.merged.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz"
MOUSE_PEAKS="/ocean/projects/bio230007p/ikaplow/MouseAtac/Liver/peak/rep1/SRR8119852_1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz"
OUTPUT_DIR="/ocean/projects/bio230007p/nrajesh2/data/halper_output"

# Make output directory
mkdir -p $OUTPUT_DIR

# Copy or unzip peak files to working directory
echo "Unzipping human peaks..."
gunzip -c $HUMAN_PEAKS > $OUTPUT_DIR/human_liver.narrowPeak

echo "Unzipping mouse peaks..."
gunzip -c $MOUSE_PEAKS > $OUTPUT_DIR/mouse_liver.narrowPeak

# Run HALPER mapping
# Update paths to hal file, source, and target species
HAL_FILE="/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal"
SOURCE_SPECIES="Homo_sapiens"
TARGET_SPECIES="Mus_musculus"

# Run halper_map_peak_orthologs.sh
sbatch halper_map_peak_orthologs.sh \
    -b $OUTPUT_DIR/human_liver.narrowPeak \
    -o $OUTPUT_DIR \
    -s $SOURCE_SPECIES \
    -t $TARGET_SPECIES \
    -c $HAL_FILE

echo "HALPER mapping job submitted!"