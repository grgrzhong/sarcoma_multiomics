#!/usr/bin/env Rscript

##############################################################################
## Author:  Zhong Guorui
## Created: 2025-06-12
## Updated: 2025-06-12
## Description: Converts a segment TSV file (e.g., from FACETS VCF-to-TSV conversion)
##              to a custom TSV format with Chromosome, Start, End, nMajor, nMinor.
##############################################################################

# Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(optparse)
        library(tidyverse)
    })
)

# Function to process the input TSV and convert it to the custom format
convert_to_custom_format <- function(input_file_path, output_file_path) {
    # message(paste("Loading input TSV file = ", input_file_path))

    if (!file.exists(input_file_path)) {
        stop(paste("Input file not found:", input_file_path), call. = FALSE)
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
            required_input_cols <- c("CHROM", "POS", "END", "TCN_EM", "LCN_EM")

            input_tsv <- readr::read_tsv(
                input_file_path, 
                show_col_types = FALSE
            )

        }, 
        error = function(e) {
            stop(paste("Error reading input TSV file:", input_file_path, "-", e$message), call. = FALSE)
        }
    )

    missing_cols <- required_input_cols[!required_input_cols %in% names(input_tsv)]
    
    if (length(missing_cols) > 0) {
        stop(paste("Input TSV file is missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
    }

    output_tsv <- input_tsv |> 
                dplyr::select(all_of(required_input_cols)) |> 
                ## The not estimateable Lesser (minor) copy number segments were filtered out
                dplyr::filter(!(LCN_EM %in% c(".", NA, "NA"))) |>
                dplyr::rename(
                    Chromosome = CHROM,
                    Start = POS,
                    End = END
                ) |> 
                mutate(
                    Chromosome=str_replace(Chromosome, "chr", ""),
                    Start = as.integer(Start),
                    End = as.integer(End),
                    TCN_EM = as.integer(TCN_EM),
                    LCN_EM = as.integer(LCN_EM)
                ) |> 
                mutate(
                    nMajor = as.integer(round(TCN_EM - LCN_EM)),
                    nMinor = as.integer(round(LCN_EM))
                ) |>
                dplyr::select(
                    Chromosome, Start, End, nMajor, nMinor
                )

    if (nrow(output_tsv) == 0) {
        warning("No valid data rows remained after processing. Output file will be empty or contain only headers.")
    }
    
    tryCatch(
        {
            readr::write_tsv(
                output_tsv, 
                file = output_file_path, 
                na = "NA"
            )
            # message(paste("Saving output TSV file = ", output_file_path))
        },
        error = function(e) {
            stop(paste("Error writing output TSV file:", output_file_path, "-", e$message), call. = FALSE)
        }
    )
}

# Parse command line arguments
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character", 
        default = NULL,
        help = "Input segment TSV file path (e.g., from FACETS VCF-to-TSV)", 
        metavar = "FILE"
    ),
    make_option(
        c("-o", "--output"),
        type = "character", 
        default = NULL,
        help = "Output custom TSV file path (columns: Chromosome, Start, End, nMajor, nMinor)", 
        metavar = "FILE"
    )
)

opt_parser <- OptionParser(
    option_list = option_list, 
    description = "Converts a segment TSV file to a custom format with nMajor/nMinor."
)

opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input TSV file path must be specified with -i or --input", call. = FALSE)
}

if (is.null(opt$output)) {
    print_help(opt_parser)
    stop("Output TSV file path must be specified with -o or --output", call. = FALSE)
}

# If script is being run directly (not sourced)
if (!interactive()) {
    convert_to_custom_format(opt$input, opt$output)
}

## Testing
# input_file_path <- "/mnt/f/projects/250224_DFSP_Multiomics/outputs/facets_test/DFSP-001-T/DFSP-001-T.tsv"
# output_file_path <- "/mnt/f/projects/250224_DFSP_Multiomics/outputs/facets_test/DFSP-001-T/DFSP-001-T.pcgr.tsv"
# convert_to_custom_format(input_file_path, output_file_path)
