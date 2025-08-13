#!/bin/bash

#############################################################################
# This script runs GISTIC2 for copy number analysis on WES data.
# To find out the recurrence of copy number alterations in a set of samples.
###############################################################################

## Setup conda env
# conda create -n gistic2 -y
# conda install hcc::gistic2 
source "$(dirname "${BASH_SOURCE[0]}")/conf/config.sh"

# Activate the conda environment
conda activate gistic2

export GISTIC2_REFERENCE="/mnt/f/Reference/GISTIC2/hg38.UCSC.add_miR.160920.refgene.mat"

# export GISTIC2_DIR="/mnt/f/projects/250224_sarcoma_multiomics/data/wes/GISTIC2/somatic_matched_removed_neutral"
# export GISTIC2_DIR="/mnt/f/projects/250224_sarcoma_multiomics/data/wes/GISTIC2/somatic_unmatched"
export GISTIC2_DIR="/mnt/f/projects/250224_sarcoma_multiomics/data/temp/GISTIC2"

## create output directory if it does not exist
mkdir -p "$GISTIC2_DIR"

## Test file
# segment_file="/mnt/f/projects/sarcoma_multiomics/data/wes/GISTIC2/FS-DFSP/FS-DFSP.tsv"
segment_file="/mnt/f/projects/250224_sarcoma_multiomics/data/temp/GISTIC2/FGFR_mskimpact_segments/FGFR_mskimpact_segments.tsv"

# Function to run GISTIC2 on one file
run_gistic2() {

    local segment_file=$1
    
    # Extract the base name of the segment file
    dir_name=$(basename "$segment_file" .tsv)
    
    # Create output directory for GISTIC2 results
    output_dir="${GISTIC2_DIR}/${dir_name}"
    mkdir -p "$output_dir"
    # echo "$output_dir"

    echo "$(date +"%F") $(date +"%T") - (${dir_name}) Running GISTIC2 analysis ..."

    # Run GISTIC2
    gistic2 \
        -b "${output_dir}" \
        -seg "${segment_file}" \
        -refgene "${GISTIC2_REFERENCE}" \
        -qvt 0.5 \
        -ta 0.3 \
        -td 0.3 \
        -conf 0.99 \
        >& "${output_dir}/gistic2.log"
        # -ta 0.1 \
        # -td 0.1 \
        # -armpeel 1 \
        # -js 8 \
        # -genegistic 1 \
        # -conf 0.99 \
        # -brlen 0.8 \
}

# Export the function for GNU Parallel
export -f run_gistic2

# Find all ready-to-use segment files and run in parallel
find "$GISTIC2_DIR" -name "*.tsv" | 
    parallel \
    --jobs 6 \
    run_gistic2 {}
