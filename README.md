# Compute TF Prior Activity

This repository contains code to compute transcription factor (TF) prior activity based on motif analysis and ATAC-seq peaks.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/rodrij92/compute_tf_prior_activity.git
cd compute_tf_prior_activity
```

### 2. Create and Activate the Conda Environment

Install as many dependencies as possible via conda (this covers most packages):

```bash
conda config --add channels conda-forge
conda config --add channels bioconda

conda create -n tf_prior_activity -c bioconda -c conda-forge \
  r-base=4.3 \
  bioconductor-genomicranges \
  bioconductor-rtracklayer \
  r-doparallel \
  r-foreach \
  r-data.table \
  r-optparse \
  r-utils \
  r-essentials \
  -y

conda activate tf_prior_activity
```

### 3. Install Remaining R Packages

Some R packages (such as `BiocManager` and occasionally `R.utils`) may need to be installed from within R. You can do this in one R command:

```bash
R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('GenomicRanges', 'rtracklayer')); install.packages(c('R.utils', 'parallel', 'doParallel', 'foreach', 'data.table', 'optparse'), repos='https://cloud.r-project.org')"
```

**Note:**
- `parallel` is a base R package and should already be available.
- If you get errors about missing packages (e.g., `GenomicRanges`), try installing them from within R as shown above.
- If you still encounter issues, ensure you are running R from within the activated conda environment.

### 4. Test the Installation

You can test that all packages are available by running:

```bash
Rscript src/prior_activity.R --test
```

If you see errors about missing packages, repeat step 3 from within your conda environment.

## Usage

The script can be run in two modes:

### Test Mode
Processes a single TF (ASCL1) with downsampled data:
```bash
Rscript src/prior_activity.R --test
```

### Batch Mode
Processes all TF motif files in the data directory:
```bash
Rscript src/prior_activity.R --batch
```

## How the Code Works

The `prior_activity.R` script computes transcription factor (TF) prior activity based on motif analysis and ATAC-seq peaks. Here's a high-level overview of the workflow:

1. **Input Data**: The script reads motif files (from the `data/meme_res_0.01/` directory) and ATAC-seq peaks (from `data/sns_atac_phs003226/`). It also uses gene annotations (from `data/gene_annot/genes.csv`) to identify promoter regions.

2. **Motif Activity Calculation**: For each TF, the script:
   - Reads the motif file and, if in test mode, downsamples it to 100,000 motifs.
   - Reads the ATAC-seq peaks file.
   - Calculates motif activity by finding overlaps between motifs and peaks. If no peaks overlap with the TF's promoter, all motif scores are set to 0.

3. **Output Generation**: The script generates two types of output files for each TF:
   - A CSV file (`results/*_prior.csv`) containing motif activity scores for each chromosome, ensuring all chromosomes (chr1-22 and chrX) are included, even if no motifs are present.
   - A BigWig file (`results/*_prior.bw`) for visualization, also ensuring all chromosomes are represented.

4. **Parallel Processing**: The script uses parallel processing to speed up the motif activity calculation, especially in batch mode.

### Example Code

To run the script in test mode (processing ASCL1 only):

```bash
Rscript src/prior_activity.R --test
```

To run the script in batch mode (processing all TF motif files):

```bash
Rscript src/prior_activity.R --batch
```

### Key Functions

- `read_motifs()`: Reads and processes motif files, ensuring consistent chromosome information.
- `read_peaks()`: Reads and processes ATAC-seq peaks.
- `calculate_motif_activity()`: Computes motif activity by finding overlaps between motifs and peaks.
- `create_bigwig()`: Generates BigWig files for visualization, ensuring all chromosomes are included.

## Input Data Structure

The script expects the following directory structure:
```
compute_tf_prior_activity/
├── data/
│   ├── gene_annot/
│   │   └── genes.csv
│   ├── meme_res_0.01/
│   │   └── *_fimo.tsv.gz
│   ├── sns_atac_phs003226/
│   │   └── MCG001_ATAC.peaks.bed
│   └── dna_sequence/
│       └── chr_sizes.txt
├── src/
│   └── prior_activity.R
└── results/
    ├── *_prior.csv
    └── *_prior.bw
```

## Output

The script generates two types of output files for each TF:
1. CSV file with motif activity scores (`results/*_prior.csv`)
2. BigWig file for visualization (`results/*_prior.bw`)

## Dependencies

- R (>= 4.0)
- Bioconductor
- GenomicRanges (Bioconductor)
- rtracklayer (Bioconductor)
- parallel
- doParallel
- foreach
- data.table
- optparse
- BiocManager

You can install the required R packages with:

```r
if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')
BiocManager::install(c('GenomicRanges', 'rtracklayer'))
install.packages(c('parallel', 'doParallel', 'foreach', 'data.table', 'optparse'))
```

## License

MIT License 