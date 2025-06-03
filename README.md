# TF Prior Activity

This repository contains code for computing transcription factor (TF) prior activity based on motif analysis and ATAC-seq peaks.

## Environment Setup

### Prerequisites
- Conda (Miniconda or Anaconda) installed
- Git installed

### Installation Steps

1. Clone the repository:
```bash
git clone https://github.com/javrodriguez/tf_prior_activity.git
cd tf_prior_activity
```

2. Create and activate the conda environment:
```bash
# Create new environment
conda env create -f environment.yml

# Activate environment
conda activate tf_prior_activity
```

3. Install required R packages:
```bash
Rscript install_packages.R
```

## Usage

The main script `prior_activity.R` computes TF prior activity scores based on motif analysis and ATAC-seq peaks.

### Input Files Required
- ATAC-seq peaks in BED format
- Motif database in MEME format
- Genome sequence in FASTA format

### Running the Script
```bash
Rscript prior_activity.R --peaks <peaks.bed> --motifs <motifs.meme> --genome <genome.fa> --output <output_dir>
```

### Output
The script generates TF activity scores for each sample, stored in the specified output directory.

## Project Structure
```
tf_prior_activity/
├── src/
├── examples/                # Example data and usage
├── data/                    # Data directory for input files
├── environment.yml         # Conda environment specification
├── install_packages.R      # R package installation script
```

## Dependencies
- R >= 4.2.0
- Required R packages:
  - data.table
  - doParallel
  - foreach
  - optparse
  - GenomicRanges
  - rtracklayer
  - GenomeInfoDb
  - Biostrings
  - BSgenome.Hsapiens.UCSC.hg38

## License
This project is licensed under the MIT License - see the LICENSE file for details. 