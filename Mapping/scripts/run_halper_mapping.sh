#!/bin/bash
#SBATCH -J halper_map             # Job name
#SBATCH -p RM-shared              # Partition
#SBATCH -N 1                      # Number of nodes
#SBATCH -n 4                     # Number of tasks (CPUs)
#SBATCH -t 15:00:00               # Walltime (hh:mm:ss)
#SBATCH --mem=4G                 # Memory
#SBATCH -o /ocean/projects/bio230007p/nrajesh2/data/halper_output/halper_%j.out
#SBATCH -e /ocean/projects/bio230007p/nrajesh2/data/halper_output/halper_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nrajesh@andrew.cmu.edu

# -------------------------------
# HALPER mapping sbatch template
# -------------------------------

# Load modules
module load anaconda3     

# Activate HAL conda environment
conda activate hal

# Make sure halLiftover is in PATH
export PATH=$HOME/repos/hal/bin:$PATH

# Make sure Python can find orthologFind
export PYTHONPATH=$HOME/repos/halLiftover-postprocessing:$PYTHONPATH


# -------------------------------
# ===== USER PARAMETERS =====
# -------------------------------

# Input peak files (gzipped narrowPeak)
HUMAN_PEAKS="/ocean/projects/bio230007p/nrajesh2/data/human_liver.narrowPeak.gz"
MOUSE_PEAKS="/ocean/projects/bio230007p/nrajesh2/data/mouse_liver.narrowPeak.gz"

# Output directory (will be created if it doesn't exist)
OUTPUT_DIR="/ocean/projects/bio230007p/nrajesh2/data/halper_output"

# HAL file and species
HAL_FILE="/ocean/projects/bio230007p/nrajesh2/data/10plusway-master.hal"
SOURCE_SPECIES="Human"
TARGET_SPECIES="Mouse"

# Path to halper script
HALPER_SCRIPT="$HOME/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh"


# Make output directory
mkdir -p $OUTPUT_DIR
mkdir -p logs

# Unzip BED files
echo "Unzipping human peaks..."
gunzip -c $HUMAN_PEAKS > $OUTPUT_DIR/human_liver.narrowPeak

echo "Unzipping mouse peaks..."
gunzip -c $MOUSE_PEAKS > $OUTPUT_DIR/mouse_liver.narrowPeak

# Run HALPER mapping
echo "Running HALPER mapping..."
$HALPER_SCRIPT \
    -b $OUTPUT_DIR/human_liver.narrowPeak \
    -o $OUTPUT_DIR \
    -s $SOURCE_SPECIES \
    -t $TARGET_SPECIES \
    -c $HAL_FILE
echo "HALPER mapping finished!"
