#!/bin/bash

# Calculate number of TF files before SLURM directives
NUM_TFS=$(ls data/meme_res_0.01/*_fimo.tsv.gz | wc -l)
ARRAY_SIZE=$((NUM_TFS - 1))

#SBATCH --job-name=tf_prior_activity
#SBATCH --output=logs/tf_prior_%A_%a.log
#SBATCH --error=logs/tf_prior_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=0-${ARRAY_SIZE}

# Create logs directory if it doesn't exist
mkdir -p logs

# Get the list of TF files
mapfile -t TF_FILES < <(ls data/meme_res_0.01/*_fimo.tsv.gz)

# Get the current TF file based on array index
CURRENT_TF=${TF_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract TF name from filename (remove path and _fimo.tsv.gz)
TF_NAME=$(basename "$CURRENT_TF" _fimo.tsv.gz)

# Load conda environment
source ~/home_abl/miniconda3/etc/profile.d/conda.sh
conda activate tf_prior_activity

# Run the script for this TF
echo "Processing $TF_NAME..."
Rscript compute_tf_prior_activity.R --peaks data/sns_atac_phs003226/MCG001_ATAC.peaks.bed \
                                  --motifs "$CURRENT_TF" \
                                  --genome data/dna_sequence/chr_sizes.txt \
                                  --output "results/${TF_NAME}_prior"

# Check if the job was successful
if [ $? -eq 0 ]; then
    echo "Successfully processed $TF_NAME"
else
    echo "Error processing $TF_NAME"
    exit 1
fi 