#!/bin/bash
#SBATCH --job-name=rgreat
#SBATCH --output=rgreat_%j.out
#SBATCH --error=rgreat_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000M
#SBATCH --account=bio230007p
#SBATCH --partition=RM-shared
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=surathas@andrew.cmu.edu

cd /ocean/projects/bio230007p/ssriram6/liver-ATAC-OCR

source activate rgreat_env

Rscript rgreat_analysis.R
























