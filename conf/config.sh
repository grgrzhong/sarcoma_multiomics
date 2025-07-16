#!/bin/bash

#############################################################################
# Somatic Mutation Workflow - Configuration Script
# This script sets up the environment and variables for the somatic mutation analysis
# workflow.
# It configures directories, container paths, and reference files.
#############################################################################

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate apptainer

# Exit on any error
set -e

# Set project directory relative to this script
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PROJECT_DIR

# Module directory containing scripts for processing data
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"

# Input and output directories
export INPUT_DIR="${1:-${PROJECT_DIR}/data/Raw}"
export OUTPUT_DIR="${2:-${PROJECT_DIR}/data}"

mkdir -p "$OUTPUT_DIR" 

export FASTQ_TRIM_DIR="${OUTPUT_DIR}/Fastq-trimmed"
export FASTQC_TRIM_DIR="${OUTPUT_DIR}/FastQC-trimmed"
export BAM_DIR="${OUTPUT_DIR}/BAM"
export MUTECT2_DIR="${OUTPUT_DIR}/Mutect2"
export CNV_DIR="${OUTPUT_DIR}/CNV"
export CNV_FACETS_DIR="${CNV_DIR}/cnv_facets"
export PCGR_DIR="${OUTPUT_DIR}/PCGR"

# Reference and annotation directories
export REFERENCE_DIR="${5:-/mnt/f/Reference}"
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
export DBSNP="${REFERENCE_DIR}/Population_database/dbSNP.vcf.gz"
export BAIT_INTERVALS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-probes-hg38.interval_list"
export TARGET_INTERVALS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-targets-hg38.interval_list"

# Print out the environment information
echo "========================================================================"
echo "Clinical RNA Fusion Analysis Workflow - Configuration"
echo "========================================================================"
echo "Project Directory:               $PROJECT_DIR"
echo "Input Directory:                 $INPUT_DIR"
echo "Output Directory:                $OUTPUT_DIR"
echo "Reference Directory:             $REFERENCE_DIR"
echo "Reference File:                  $REFERENCE"
echo "Interval File:                   $INTERVAL"
echo "DBSNP File:                      $DBSNP"
echo "Bait Intervals:                  $BAIT_INTERVALS"
echo "Target Intervals:                $TARGET_INTERVALS"
echo "Fastq Trimmed Directory:         $FASTQ_TRIM_DIR"
echo "FastQC Trimmed Directory:        $FASTQC_TRIM_DIR"
echo "BAM Directory:                   $BAM_DIR"
echo "Mutect2 Directory:               $MUTECT2_DIR"
echo "CNV Directory:                   $CNV_DIR"
echo "CNV Facets Directory:            $CNV_FACETS_DIR"
echo "PCGR Directory:                  $PCGR_DIR"
echo "Module Directory:                $MODULE_DIR"

# Setup containers if not already done
echo "======================================================================="
echo "Clinical RNA Fusion Analysis Workflow - Container Setup"
echo "======================================================================="
export CONTAINER_DIR="${6:-${PROJECT_DIR}}/containers"

if [ ! -d "$CONTAINER_DIR" ]; then
    echo "$(date +"%F") $(date +"%T") Container directory not found. Setting up containers..."
    mkdir -p "$CONTAINER_DIR"
    bash "${PROJECT_DIR}/conf/containers.sh"
    echo "$(date +"%F") $(date +"%T") ✓ Container setup completed"
else
    echo "$(date +"%F") $(date +"%T") Container directory already exists: $CONTAINER_DIR"
    echo "$(date +"%F") $(date +"%T") Skipping container setup step."
fi
