#!/bin/bash

#############################################################################
# This script runs GISTIC2 for copy number analysis on WES data.
# To find out the recurrence of copy number alterations in a set of samples.
###############################################################################

## Setup conda env
# conda create -n gistic2 -y
# conda install hcc::gistic2 

# Activate the conda environment
conda activate gistic2

export GISTIC2_DIR="/mnt/f/projects/sarcoma_multiomics/data/wes/GISTIC2"
export GISTIC2_REFERENCE="/mnt/f/Reference/GISTIC2/hg38.UCSC.add_miR.160920.refgene.mat"

## create output directory if it does not exist
mkdir -p "$GISTIC2_DIR"

# segment_file=/mnt/f/projects/sarcoma_multiomics/data/wes/GISTIC2/gistic2_group/U-DFSP/U-DFSP.seg

# Function to run GISTIC2 on one file
run_gistic2() {

    local segment_file=$1
    
    # Extract the base name of the segment file
    file_name=$(basename "$segment_file" .tsv)
    
    # Create output directory for GISTIC2 results
    output_dir="${GISTIC2_DIR}/${file_name}"
    mkdir -p "$output_dir"

    echo "$(date +"%F") $(date +"%T") - (${file_name}) Running GISTIC2 analysis ..."

    # Run GISTIC2
    gistic2 \
        -b "${output_dir}" \
        -seg "${segment_file}" \
        -refgene "${GISTIC2_REFERENCE}" \
        -qvt 0.1 \
        -ta 0.3 \
        -td 0.3 \
        -conf 0.99 \
        >& "${output_dir}/gistic2.log"
}

# Export the function for GNU Parallel
export -f run_gistic2

# Find all ready-to-use segment files and run in parallel
find "$GISTIC2_DIR" -name "*.tsv" | 
    parallel \
    --jobs 3 \
    run_gistic2 {}
