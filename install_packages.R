#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install required CRAN packages
cran_packages <- c("data.table", "doParallel", "foreach", "optparse")
for (pkg in cran_packages) {
    if (!require(pkg, quietly = TRUE)) {
        tryCatch({
            install.packages(pkg, repos = "https://cloud.r-project.org")
        }, error = function(e) {
            message(sprintf("Error installing %s: %s", pkg, e$message))
        })
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
tryCatch({
    BiocManager::install(bioc_packages, update = FALSE)
}, error = function(e) {
    message(sprintf("Error installing Bioconductor packages: %s", e$message))
})

# Verify installations
message("\nVerifying package installations...")
for (pkg in c(cran_packages, bioc_packages)) {
    if (require(pkg, quietly = TRUE)) {
        message(sprintf("✓ %s installed successfully", pkg))
    } else {
        message(sprintf("✗ Failed to install %s", pkg))
    }
} 