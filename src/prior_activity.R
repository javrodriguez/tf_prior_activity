#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(data.table)
  library(optparse)
})

# Get number of available cores
num_cores <- parallel::detectCores()

# Parse command line arguments
option_list <- list(
  make_option(c("-p", "--peaks"), type="character", default=NULL,
              help="Path to peaks file in BED format"),
  make_option(c("-m", "--motifs"), type="character", default=NULL,
              help="Path to motifs file in FIMO TSV format"),
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Path to genome file with chromosome sizes"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory for results"),
  make_option(c("--test"), action="store_true", default=FALSE,
              help="Run in test mode with downsampled motifs (default: FALSE)"),
  make_option(c("--batch"), action="store_true", default=FALSE,
              help="Process all motif files in batch mode (default: FALSE)"),
  make_option(c("--motifs-dir"), type="character", default=NULL,
              help="Directory containing motif files (required for batch mode)"),
  make_option(c("--gene-annot"), type="character", default=NULL,
              help="Path to gene annotations file (required for TF promoter analysis)"),
  make_option(c("--cores"), type="integer", default=num_cores,
              help=sprintf("Number of CPU cores to use (default: %d, all available cores)", num_cores))
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Debug: Print all arguments
message("Command line arguments:")
message("peaks: ", opt$peaks)
message("motifs: ", opt$motifs)
message("genome: ", opt$genome)
message("output: ", opt$output)
message("gene-annot: ", opt$`gene-annot`)
message("test: ", opt$test)
message("batch: ", opt$batch)
message("motifs_dir: ", opt$motifs_dir)
message("cores: ", opt$cores)

# Validate required arguments
if (is.null(opt$peaks)) {
  stop("Peaks file path is required. Use --peaks")
}
if (is.null(opt$motifs) && !opt$batch) {
  stop("Motifs file path is required. Use --motifs")
}
if (is.null(opt$genome)) {
  stop("Genome file path is required. Use --genome")
}
if (is.null(opt$output)) {
  stop("Output directory is required. Use --output")
}
if (opt$batch && is.null(opt$motifs_dir)) {
  stop("Motifs directory is required in batch mode. Use --motifs-dir")
}
if (is.null(opt$`gene-annot`)) {
  stop("Gene annotations file is required. Use --gene-annot")
}

# Create output directory if it doesn't exist
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

# Function to read chromosome sizes
read_chr_sizes <- function(file_path) {
  chr_sizes <- fread(file_path, header = FALSE, col.names = c("chr", "size"))
  # Add 'chr' prefix if not present
  chr_sizes$chr <- ifelse(grepl("^chr", chr_sizes$chr), chr_sizes$chr, paste0("chr", chr_sizes$chr))
  # Strip '.fa.gz' suffix
  chr_sizes$chr <- sub("\\.fa\\.gz$", "", chr_sizes$chr)
  
  # Filter to only include chr1-22 and chrX
  valid_chrs <- c(paste0("chr", 1:22), "chrX")
  chr_sizes <- chr_sizes[chr_sizes$chr %in% valid_chrs]
  
  # Sort chromosomes in the correct order
  chr_sizes$chr <- factor(chr_sizes$chr, levels = valid_chrs)
  chr_sizes <- chr_sizes[order(chr_sizes$chr), ]
  
  return(chr_sizes)
}

# Function to filter GRanges object to valid chromosomes
filter_to_valid_chrs <- function(gr) {
  valid_chrs <- c(paste0("chr", 1:22), "chrX")
  gr <- gr[seqnames(gr) %in% valid_chrs]
  seqlevels(gr) <- valid_chrs
  return(gr)
}

# Function to read and process motifs file
read_motifs <- function(file_path, test_mode = FALSE) {
  message(sprintf("Reading motifs file: %s", file_path))
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop(sprintf("Motifs file does not exist: %s", file_path))
  }
  
  # Check if file is empty
  if (file.size(file_path) == 0) {
    stop(sprintf("Motifs file is empty: %s", file_path))
  }
  
  # Try to read the file with error handling
  tryCatch({
    # Read the file with read.table, specifying tab as separator
    motifs <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    if (test_mode) {
      message("Test mode: Downsampling motifs to 100,000...")
      if (nrow(motifs) > 100000) {
        set.seed(42)  # For reproducibility
        motifs <- motifs[sample(nrow(motifs), 100000), ]
      }
    }
    
    message(sprintf("Using %d motifs", nrow(motifs)))
    
    # Add 'chr' prefix if not present
    motifs$sequence_name <- ifelse(grepl("^chr", motifs$sequence_name), motifs$sequence_name, paste0("chr", motifs$sequence_name))
    
    # Read chromosome sizes to filter motifs
    chr_sizes <- read_chr_sizes(opt$genome)
    valid_chrs <- chr_sizes$chr
    
    # Filter motifs to include only valid chromosomes
    motifs <- motifs[motifs$sequence_name %in% valid_chrs, ]
    
    # Convert to GRanges
    gr <- GRanges(
      seqnames = motifs$sequence_name,
      ranges = IRanges(start = motifs$start, end = motifs$stop),
      strand = motifs$strand,
      score = motifs$score
    )
    
    # Ensure consistent chromosome set
    gr <- filter_to_valid_chrs(gr)
    
    return(gr)
  }, error = function(e) {
    stop(sprintf("Failed to read motifs file %s. Error: %s", file_path, e$message))
  })
}

# Function to read and process peaks file
read_peaks <- function(file_path) {
  message("Reading peaks file...")
  peaks <- fread(file_path)
  message(sprintf("Using all %d peaks", nrow(peaks)))
  
  # Ensure strand values are valid
  peaks[, V4 := ifelse(V4 == ".", "*", V4)]
  peaks[, V4 := ifelse(!V4 %in% c("+", "-", "*"), "*", V4)]
  
  # Add 'chr' prefix if not present
  peaks$V1 <- ifelse(grepl("^chr", peaks$V1), peaks$V1, paste0("chr", peaks$V1))
  
  # Read chromosome sizes to filter peaks
  chr_sizes <- read_chr_sizes(opt$genome)
  valid_chrs <- chr_sizes$chr
  
  # Filter peaks to include only valid chromosomes
  peaks <- peaks[peaks$V1 %in% valid_chrs]
  
  gr <- GRanges(
    seqnames = peaks$V1,  # First column is chromosome
    ranges = IRanges(start = peaks$V2, end = peaks$V3),  # Second and third columns are start and end
    strand = peaks$V4  # Fourth column is strand
  )
  
  # Ensure consistent chromosome set
  gr <- filter_to_valid_chrs(gr)
  
  return(gr)
}

# Function to read gene annotations
read_gene_annotations <- function(file_path) {
  message(sprintf("Reading gene annotations from: %s", file_path))
  if (!file.exists(file_path)) {
    stop(sprintf("Gene annotations file does not exist: %s", file_path))
  }
  
  tryCatch({
    genes <- fread(file_path)
    message("Successfully read gene annotations file")
    message(sprintf("Number of genes: %d", nrow(genes)))
    message("Columns found: ", paste(names(genes), collapse=", "))
    
    # Check required columns
    required_cols <- c("chr", "start", "end", "strand", "gene_name")
    missing_cols <- setdiff(required_cols, names(genes))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns in gene annotations file: %s", 
                  paste(missing_cols, collapse=", ")))
    }
    
    # Calculate TSS based on strand
    genes[, tss := ifelse(strand == "+", start, end)]
    
    # Define promoter regions based on strand
    # For + strand: TSS-2000 to TSS+500
    # For - strand: TSS-500 to TSS+2000
    genes[, promoter_start := ifelse(strand == "+", tss - 2000, tss - 500)]
    genes[, promoter_end := ifelse(strand == "+", tss + 500, tss + 2000)]
    
    # Create GRanges object for promoter regions
    promoter_gr <- GRanges(
      seqnames = genes$chr,
      ranges = IRanges(start = genes$promoter_start, end = genes$promoter_end),
      gene_name = genes$gene_name,
      strand = genes$strand
    )
    
    # Ensure consistent chromosome set
    promoter_gr <- filter_to_valid_chrs(promoter_gr)
    
    return(promoter_gr)
  }, error = function(e) {
    stop(sprintf("Error reading gene annotations file: %s", e$message))
  })
}

# Function to calculate motif activity
calculate_motif_activity <- function(motifs_gr, peaks_gr, tf_name, cores = 1) {
  message("Calculating motif activity...")
  
  # Read gene annotations
  promoter_gr <- read_gene_annotations(opt$`gene-annot`)
  
  # Find the promoter for the transcription factor
  tf_promoter <- promoter_gr[promoter_gr$gene_name == tf_name]
  if (length(tf_promoter) == 0) {
    stop(sprintf("Could not find promoter for transcription factor %s", tf_name))
  }
  
  # Check if any peak overlaps with the TF's promoter
  promoter_overlaps <- findOverlaps(tf_promoter, peaks_gr)
  if (length(promoter_overlaps) == 0) {
    message(sprintf("No peaks overlap with %s promoter. Skipping this TF.", tf_name))
    quit(status = 0)  # Exit with success status
  }
  
  message(sprintf("Found %d peaks overlapping with %s promoter. Proceeding with motif activity calculation.", 
                 length(promoter_overlaps), tf_name))
  
  # Set up parallel processing
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Split motifs into chunks for parallel processing
  chunk_size <- ceiling(length(motifs_gr) / cores)
  motif_chunks <- split(motifs_gr, ceiling(seq_along(motifs_gr) / chunk_size))
  
  # Process chunks in parallel
  results <- foreach(chunk = motif_chunks, 
                    .combine = c, 
                    .packages = c("GenomicRanges", "IRanges")) %dopar% {
    # Find overlaps between motifs and peaks
    overlaps <- findOverlaps(chunk, peaks_gr)
    
    # Initialize activity scores for this chunk
    chunk_scores <- numeric(length(chunk))
    
    # For motifs with overlaps, sum their scores
    if (length(overlaps) > 0) {
      overlap_idx <- queryHits(overlaps)
      for (i in unique(overlap_idx)) {
        chunk_scores[i] <- chunk$score[i]
      }
    }
    return(chunk_scores)
  }
  
  # Clean up parallel processing
  stopCluster(cl)
  
  # Min-Max normalization directly on raw motif scores
  message("Applying Min-Max normalization directly to raw motif scores...")
  min_score <- min(results)
  max_score <- max(results)
  if (max_score == min_score) {
    results <- rep(0, length(results))
    message("All motif scores are identical; setting all normalized scores to 0.")
  } else {
    results <- (results - min_score) / (max_score - min_score)
  }
  message(sprintf("Score transformation statistics:"))
  message(sprintf("  Original score range: [%.2f, %.2f]", min_score, max_score))
  message(sprintf("  Transformed score range: [%.2f, %.2f]", min(results), max(results)))
  message(sprintf("  Number of non-zero scores: %d", sum(results > 0)))
  
  return(results)
}

# Function to create BigWig file
create_bigwig <- function(motifs_gr, activity_scores, output_file) {
  message("Creating BigWig file...")
  
  # Create a new GRanges object with activity scores
  activity_gr <- motifs_gr
  mcols(activity_gr)$score <- activity_scores
  
  # Sort by chromosome and start position
  activity_gr <- sort(activity_gr)
  
  # Process each chromosome separately to reduce memory usage
  message("Merging overlapping motifs by chromosome...")
  merged_chrs <- GRangesList()  # Use GRangesList instead of list
  for (chr in unique(seqnames(activity_gr))) {
    message(sprintf("Processing chromosome %s...", chr))
    chr_gr <- activity_gr[seqnames(activity_gr) == chr]
    
    # Merge overlapping ranges, keeping the maximum score, ignoring strand
    chr_merged <- reduce(chr_gr, with.revmap=TRUE, ignore.strand=TRUE)
    # For each merged range, take the maximum score from the original ranges
    chr_merged$score <- sapply(chr_merged$revmap, function(x) max(activity_scores[x]))
    chr_merged$revmap <- NULL  # Remove the revmap column
    
    merged_chrs[[chr]] <- chr_merged
    gc()  # Force garbage collection after each chromosome
  }
  
  # Combine all chromosomes
  activity_gr <- unlist(merged_chrs)
  
  # Read chromosome sizes and set seqlengths
  chr_sizes <- read_chr_sizes(opt$genome)
  
  # Ensure consistent chromosome set
  valid_chrs <- c(paste0("chr", 1:22), "chrX")
  activity_gr <- filter_to_valid_chrs(activity_gr)
  
  # Create empty ranges for chromosomes that don't have any motifs
  missing_chrs <- setdiff(valid_chrs, as.character(seqnames(activity_gr)))
  if (length(missing_chrs) > 0) {
    message(sprintf("Adding empty ranges for chromosomes: %s", paste(missing_chrs, collapse=", ")))
    # Create empty ranges for missing chromosomes
    empty_ranges <- GRanges(
      seqnames = missing_chrs,
      ranges = IRanges(start = 1, end = 1),
      score = 0
    )
    # Combine with existing ranges
    activity_gr <- c(activity_gr, empty_ranges)
  }

  # Set seqlevels to all valid chromosomes
  seqlevels(activity_gr) <- valid_chrs

  # Set chromosome lengths for all chromosomes
  seqlengths(activity_gr) <- chr_sizes$size[match(valid_chrs, chr_sizes$chr)]
  
  # Export as BigWig
  message(sprintf("Saving BigWig to %s...", output_file))
  export(activity_gr, output_file, format = "BigWig")
}

# Function to validate BigWig file
validate_bigwig <- function(bigwig_file) {
  message("Validating BigWig file...")
  bw <- import(bigwig_file)
  
  # Check if all expected chromosomes are present
  expected_chrs <- c(paste0("chr", 1:22), "chrX")
  missing_chrs <- setdiff(expected_chrs, unique(seqnames(bw)))
  if (length(missing_chrs) > 0) {
    warning(sprintf("Missing chromosomes in BigWig file: %s", paste(missing_chrs, collapse=", ")))
  }
  
  # Check if activity scores are non-negative
  neg_scores <- score(bw) < 0
  if (any(neg_scores)) {
    warning("Negative activity scores found in BigWig file.")
    neg_chrs <- unique(seqnames(bw)[neg_scores])
    message(sprintf("Negative scores found on chromosomes: %s", paste(neg_chrs, collapse=", ")))
    message(sprintf("Total number of negative scores: %d", sum(neg_scores)))
  }
  
  # Check if activity scores are within a reasonable range (e.g., < 1000)
  if (any(score(bw) > 1000)) {
    warning("Unusually high activity scores found in BigWig file.")
  }
  
  message("BigWig validation complete.")
}

# Main function
main <- function() {
  # Start timing
  start_time <- Sys.time()

  # Helper to process a single motif file
  process_one <- function(motifs_file) {
    tf_name <- gsub("_fimo.tsv.gz$", "", basename(motifs_file))
    output_csv <- sprintf("%s/%s_prior.csv", opt$output, tf_name)
    output_bw <- sprintf("%s/%s_prior.bw", opt$output, tf_name)

    motifs_gr <- read_motifs(motifs_file, test_mode = opt$test)
    peaks_gr <- read_peaks(opt$peaks)
    activity_scores <- calculate_motif_activity(motifs_gr, peaks_gr, tf_name, cores = opt$cores)

    # Create a data frame with all chromosomes, even if empty
    valid_chrs <- c(paste0("chr", 1:22), "chrX")
    all_chrs <- data.frame(
      chromosome = rep(valid_chrs, each = 1),
      start = 1,
      end = 1,
      activity_score = 0  # Initialize all scores to 0
    )
    
    # Update scores for chromosomes that have motifs
    if (length(activity_scores) > 0) {
      chr_scores <- tapply(activity_scores, seqnames(motifs_gr), max)
      all_chrs$activity_score[all_chrs$chromosome %in% names(chr_scores)] <- 
        chr_scores[match(all_chrs$chromosome[all_chrs$chromosome %in% names(chr_scores)], names(chr_scores))]
    }
    
    # Ensure all scores are numeric and replace NA with 0
    all_chrs$activity_score <- as.numeric(all_chrs$activity_score)
    all_chrs$activity_score[is.na(all_chrs$activity_score)] <- 0
    
    message(sprintf("Saving CSV results to %s...", output_csv))
    fwrite(all_chrs, output_csv)
    create_bigwig(motifs_gr, activity_scores, output_bw)
    
    # Validate the BigWig file
    validate_bigwig(output_bw)
  }

  if (opt$batch) {
    motif_files <- list.files(opt$motifs_dir, pattern = "_fimo.tsv.gz$", full.names = TRUE)
    message(sprintf("Batch mode: Found %d motif files.", length(motif_files)))
    for (motifs_file in motif_files) {
      message(sprintf("\nProcessing %s", motifs_file))
      process_one(motifs_file)
    }
  } else if (opt$test) {
    message("Test mode: Processing ASCL1 only...")
    if (!file.exists(opt$motifs)) {
      stop(sprintf("Test file %s not found", opt$motifs))
    }
    process_one(opt$motifs)
  } else {
    # Single file mode
    if (!file.exists(opt$motifs)) {
      stop(sprintf("Motif file %s not found", opt$motifs))
    }
    message(sprintf("Processing single motif file: %s", opt$motifs))
    process_one(opt$motifs)
  }

  # Print timing information
  end_time <- Sys.time()
  message(sprintf("\nTotal runtime: %.2f minutes", 
                 as.numeric(difftime(end_time, start_time, units = "mins"))))
}

# Run the main function
main() 