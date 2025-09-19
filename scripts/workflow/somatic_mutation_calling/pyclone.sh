#!/bin/bash

## "==========================================================================="
## Credited from PyClone (https://github.com/Roth-Lab/pyclone)
## Authors: Zhong Guorui
## Date: 2025-09-05
## "==========================================================================="

## Activate the conda environment
conda activate pyclone

export PYCLONE_DATA_DIR=/mnt/f/projects/250224_DFSP_Multiomics/data/WES/PyClone_2

## Run PyClone for each case
case_ids=$(find "${PYCLONE_DATA_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)    

## Run PyClone for a given case
for case_id in $case_ids; do

    case_dir="${PYCLONE_DATA_DIR}/${case_id}"

    echo "$(date +"%F") $(date +"%T") - ${case_id} - Running PyClone ... "    
    # Change to the case directory
    cd "${case_dir}" || {
        echo "Error: Cannot change to directory ${case_dir}"
        continue
    }

    tsv_files=$(ls *.tsv)
    # echo "Files: $tsv_files"

    PyClone run_analysis_pipeline \
        --in_files ${tsv_files} \
        --working_dir . \
        >& pyclone.log

done
