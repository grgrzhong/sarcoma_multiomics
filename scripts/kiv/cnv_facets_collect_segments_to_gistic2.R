#!/usr/bin/env Rscript

##############################################################################
## Author:  Zhong Guorui
## Created: 2025-06-12
## Updated: 2025-06-12
## Description: Collects multiple FACETS annotated TSV segment files from a directory
##              and converts them to a single custom TSV format for GISTIC2 input.
##              Output format: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean
## GISTIC2 is designed to identify significantly amplified and deleted regions across a cohort
## Note: GISTIC only works with regions that are covered by segments across all samples.
##############################################################################

# Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(optparse)
        library(fs)
        library(here)
        library(readr)
        library(tidyverse)
    })
)

# Parse command line arguments
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character", 
        default = NULL,
        help = "Input TSV directory ", 
        metavar = "DIR"
    ),
    make_option(
        c("-o", "--output"),
        type = "character", 
        default = NULL,
        help = "Output custom segment file path (columns: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean)", 
        metavar = "FILE"
    )
)

opt_parser <- OptionParser(
    option_list = option_list, 
    description = "Collect segment TSV files and convert to a custom format for GISTIC2 input."
)

opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input TSV directory must be specified with -i or --input", call. = FALSE)
}

if (is.null(opt$output)) {
    print_help(opt_parser)
    stop("Output custom segment file path must be specified with -o or --output", call. = FALSE)
}

# Function to process the input TSV and convert it to the custom format
collect_segments_to_gistic2 <- function(input_dir, output_file_path) {
    # message(paste("Loading input TSV file = ", input_file_path))

    if (!dir.exists(input_dir)) {
        stop(paste("Input directory not found:", input_dir), call. = FALSE)
    }

    # Create output directory if it doesn't exist
    output_dir <- dirname(output_file_path)
    
    if (!dir.exists(output_dir)) {
        message(paste("Output directory does not exist, creating:", output_dir))
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    tryCatch(
        {
            # Read input, treating "." as NA for relevant columns
            # https://github.com/dariober/cnv_facets

            # Get all TSV files
            tsv_files <- dir_ls(input_dir, recurse = TRUE, glob = "*.annotated.tsv")
            
            ## Get the vcf2tsv files
            # tsv_files <- tsv_files[grepl("/[^/]+\\.tsv$", tsv_files) & 
            #                     !grepl("/[^/]+\\.[^.]+\\.tsv$", tsv_files)]
            
            if (length(tsv_files) == 0) {
                stop(paste("No TSV files found in directory:", input_dir), call. = FALSE)
            }

            ## Load the data
            input_tsv <- purrr::map_dfr(
                tsv_files, 
                function(file) {
                    
                    sample_name <- path_file(file) |> str_remove(".annotated.tsv")
                    
                    message(paste("Processing sample = ", sample_name))

                    df <- readr::read_tsv(
                        here::here(file), 
                        show_col_types = FALSE
                        # col_types = cols(.default = "c")
                    )

                    ## Select relevant columns and rename them
                    required_columns <- c(
                        "Sample.ID",
                        "Chromosome", 
                        "Start_Position", 
                        "End_Position", 
                        "Num_Markers", 
                        "Seg.CN",
                        "SVTYPE"
                    )

                    df |>
                        dplyr::select(all_of(required_columns)) |>
                        filter(SVTYPE != "NEUTR") |>
                        distinct() |>
                        rename(
                            Sample = Sample.ID,
                            Chromosome = Chromosome,
                            Start = Start_Position,
                            End = End_Position,
                            Num_Probes = Num_Markers,
                            Segment_Mean = Seg.CN
                        ) |> 
                        mutate(
                            Chromosome = str_replace(Chromosome, "chr", ""),
                            Start = as.numeric(Start),
                            End = as.numeric(End),
                            Num_Probes = as.numeric(Num_Probes),
                            Segment_Mean = as.numeric(Segment_Mean)
                        ) |> 
                        select(-SVTYPE) |> 
                        ## Make sure the the segments to be start > end
                        mutate(
                            Start_fixed = if_else(Start > End, End, Start),
                            End_fixed = if_else(Start > End, Start, End)
                        ) |> 
                        select(
                            Sample, Chromosome, Start = Start_fixed, End = End_fixed, 
                            Num_Probes, Segment_Mean
                        )
                }
            )

            # Write the combined data to output file
            message(paste("Writing output to:", output_file_path))
            
            readr::write_tsv(input_tsv, output_file_path)

        },
        error = function(e) {
            stop(paste("Error processing TSV files:", e$message), call. = FALSE)
        }
    )
}

# If script is being run directly (not sourced)
if (!interactive()) {
    collect_segments_to_gistic2(opt$input, opt$output)
}

## Test
input_dir <- "/mnt/f/projects/250224_sarcoma_multiomics/data/WES/cnv_facets"
output_file <- "/mnt/f/projects/250224_sarcoma_multiomics/data/collect/dfsp_wes_cnvfacets_cohort_segments_gistic2.tsv"
collect_segments_to_gistic2(input_dir, output_file)
