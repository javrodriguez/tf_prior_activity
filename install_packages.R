#!/usr/bin/env Rscript

# List of required packages
required_packages <- c(
  "GenomicRanges",
  "rtracklayer",
  "parallel",
  "doParallel",
  "foreach",
  "data.table",
  "optparse",
  "R.utils"  # Added for reading gzipped files
)

# Function to install packages if they're not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    message(sprintf("Installing package: %s", package))
    install.packages(package, repos = "https://cloud.r-project.org")
  } else {
    message(sprintf("Package already installed: %s", package))
  }
}

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install required packages
message("Installing required packages...")
for (package in required_packages) {
  if (package %in% c("GenomicRanges", "rtracklayer")) {
    # Install Bioconductor packages
    if (!require(package, character.only = TRUE)) {
      message(sprintf("Installing Bioconductor package: %s", package))
      BiocManager::install(package, update = FALSE)
    } else {
      message(sprintf("Bioconductor package already installed: %s", package))
    }
  } else {
    # Install CRAN packages
    install_if_missing(package)
  }
}

message("All required packages have been installed.") 