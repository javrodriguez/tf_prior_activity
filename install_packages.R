#!/usr/bin/env Rscript

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install required CRAN packages
cran_packages <- c("data.table", "doParallel", "foreach", "optparse")
for (pkg in cran_packages) {
    if (!require(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

# Install required Bioconductor packages
bioc_packages <- c(
    "GenomicRanges",
    "rtracklayer",
    "GenomeInfoDb",
    "Biostrings",
    "BSgenome.Hsapiens.UCSC.hg38"
)

# Install Bioconductor packages
BiocManager::install(bioc_packages, update = FALSE)

# Verify installations
message("\nVerifying package installations...")
for (pkg in c(cran_packages, bioc_packages)) {
    if (require(pkg, quietly = TRUE)) {
        message(sprintf("✓ %s installed successfully", pkg))
    } else {
        message(sprintf("✗ Failed to install %s", pkg))
    }
} 