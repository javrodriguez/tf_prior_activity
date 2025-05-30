#!/bin/bash
#SBATCH --job-name=tf_prior_activity
#SBATCH --output=tf_prior_activity_%j.log
#SBATCH --error=tf_prior_activity_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=general

# Load conda environment
source ~/home_abl/miniconda3/etc/profile.d/conda.sh
conda activate tf_prior_activity

# Run the script
Rscript compute_tf_prior_activity.R --peaks data/sns_atac_phs003226/MCG001_ATAC.peaks.bed \
                                  --motifs "data/meme_res_0.01/*_fimo.tsv.gz" \
                                  --genome data/dna_sequence/chr_sizes.txt \
                                  --output results 