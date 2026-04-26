#!/bin/bash
#SBATCH --job-name=rgreat_plots
#SBATCH --output=rgreat_plots_%j.out
#SBATCH --error=rgreat_plots_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000M
#SBATCH --account=bio230007p
#SBATCH --partition=RM-shared
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=surathas@andrew.cmu.edu

REPO_ROOT="/ocean/projects/bio230007p/ssriram6/liver-ATAC-OCR" 

cd "${REPO_ROOT}"
source activate rgreat_env
Rscript rGREAT_Analysis/scripts/rgreat_plots.R "${REPO_ROOT}"  