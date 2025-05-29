#!/bin/bash
#SBATCH --job-name=tf_prior
#SBATCH --output=logs/tf_prior_%A_%a.out
#SBATCH --error=logs/tf_prior_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --array=0-$(ls data/meme_res_0.01/*_fimo.tsv.gz | wc -l)-1

# Create logs directory if it doesn't exist
mkdir -p logs

# Get the list of TF files
mapfile -t TF_FILES < <(ls data/meme_res_0.01/*_fimo.tsv.gz)

# Get the current TF file based on array index
CURRENT_TF=${TF_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract TF name from filename (remove path and _fimo.tsv.gz)
TF_NAME=$(basename "$CURRENT_TF" _fimo.tsv.gz)

# Load required modules (adjust these based on your HPC setup)
module load R/4.4.0
module load conda

# Activate conda environment
source activate tf_prior_activity

# Run the R script for this TF
echo "Processing $TF_NAME..."
Rscript src/prior_activity.R \
  --motifs "$CURRENT_TF" \
  --peaks data/atac_peaks/peaks.bed \
  --output "results/${TF_NAME}_prior" \
  --cores 8

# Check if the job was successful
if [ $? -eq 0 ]; then
    echo "Successfully processed $TF_NAME"
else
    echo "Error processing $TF_NAME"
    exit 1
fi 