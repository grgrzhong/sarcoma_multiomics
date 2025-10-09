#!/bin/bash


conda activate renv
export CONTAINER_DIR="/mnt/f/projects/250224_sarcoma_multiomics/containers"
export PURECN="/home/zhonggr/miniforge3/envs/renv/lib/R/library/PureCN/extdata"

Rscript $PURECN/PureCN.R --help

export PROJECT_DIR="/mnt/f/projects/250224_sarcoma_multiomics"
export REFERENCE_DIR="/mnt/f/Reference"
export FASTA="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
export PURECN_REF_DIR="${REFERENCE_DIR}/PureCN"
export PURECN_DIR="${PROJECT_DIR}/data/WES/PureCN"
export BAM_DIR="${PROJECT_DIR}/data/WES/BAM"

export PON_DIR="${REFERENCE_DIR}/PON-Mutect"
export PON="${PON_DIR}/pon.vcf.gz"

export MUTECT2_DIR="${PROJECT_DIR}/data/WES/Mutect2"

## Generate reference files for PureCN
singularity exec \
    --bind ${REFERENCE_DIR} \
    ${CONTAINER_DIR}/purecn.sif \
    Rscript /usr/local/lib/R/site-library/PureCN/extdata/IntervalFile.R \
        --in-file ${INTERVAL} \
        --fasta ${FASTA} \
        --off-target \
        --genome hg38 \
        --force \
        --out-file $PURECN_REF_DIR/baits_hg38_intervals.txt \
        --export $PURECN_REF_DIR/baits_optimized_hg38.bed 

## For each sample, tumor and normal, calculate GC-normalized coverages:
export sample_id=DFSP-001-T
mkdir -p ${PURECN_DIR}/${sample_id}

calculate_coverage() {
    local sample_id=$1
    
    ## Create output directory if it doesn't exist
    mkdir -p "${PURECN_DIR}/${sample_id}"

    ## Calculate GC-normalized coverages
    singularity exec \
        --bind ${PROJECT_DIR} \
        --bind ${REFERENCE_DIR} \
        ${CONTAINER_DIR}/purecn.sif \
        Rscript $PURECN/Coverage.R \
            --bam "${BAM_DIR}/${sample_id}/${sample_id}_recalibrated.bam" \
            --out-dir "${PURECN_DIR}/${sample_id}" \
            --intervals "${PURECN_REF_DIR}/baits_hg38_intervals.txt" \
            --cores 4 \
            >& "${PURECN_DIR}/${sample_id}/coverage.log"
}

export -f calculate_coverage

sample_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)
PARALLEL_JOBS=4

echo "${sample_ids}" | parallel \
    --jobs "$PARALLEL_JOBS" \
    calculate_coverage {}

## build a normal database for coverage normalization
ls -a $OUT/normals/*_coverage.txt.gz | cat > normal_coverages.list

singularity exec \
    --bind ${PROJECT_DIR} \
    --bind ${REFERENCE_DIR} \
    ${CONTAINER_DIR}/purecn.sif \
    Rscript $PURECN/NormalDB.R \
        --out-dir "${PURECN_REF_DIR}" \
        --coverage-files normal_coverages.list \
        --normal-panel "${PON}" \
        --genome hg38

## Run the main PureCN analysis for each tumor sample
purecn() {
    
    local tumor_id=$1

    singularity exec \
        --bind ${PROJECT_DIR} \
        --bind ${REFERENCE_DIR} \
        Rscript $PURECN/PureCN.R \
            --out "${PURECN_DIR}/${tumor_id}" \
            --tumor "${PURECN_DIR}/${tumor_id}/${tumor_id}_coverage_loess.txt.gz" \
            --sampleid "${tumor_id}" \
            --vcf "${MUTECT2_DIR}/${tumor_id}/${tumor_id}.mutect2.vcf.gz" \
            --normaldb ${PURECN_REF_DIR}/normalDB_hg38.rds \
            --intervals ${PURECN_REF_DIR}/baits_hg38_intervals.txt \
            --genome hg38

}

tumour_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep "T" | sort)

export PARALLEL_JOBS=20

echo "${tumour_ids}" | parallel \
    --jobs "$PARALLEL_JOBS" \
    purecn {}