#!/bin/bash

export PROJECT_DIR="/mnt/f/projects/sarcoma_multiomics"
export REFERENCE_DIR="/mnt/m/Reference"

## ===========================================================================
## General configurations
## ===========================================================================
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export BAM_DIR="${PROJECT_DIR}/data/WES/BAM"
export MUTECT2_DIR="${PROJECT_DIR}/data/WES/Mutect2"
export CNA_DIR="${PROJECT_DIR}/data/WES/CNV/cnv_facets"

export PCGR_REFERENCE_DIR="${REFERENCE_DIR}/PCGR_reference/20250314"
export VEP_CACHE_DIR="${REFERENCE_DIR}/VEP_cache"
export PON_DIR="${REFERENCE_DIR}/PON-Mutect"
export PON="${PON_DIR}/pon.vcf.gz"
export PON_PCGR="${PON_DIR}/pon_pcgr.vcf.gz"

## Create PCGR-compatible PoN if it doesn't exist
if [ ! -f "${PON_PCGR}" ]; then
    echo "Creating PCGR-compatible PoN VCF..."
    { 
        bcftools view -h "${PON}" | head -n -1
        echo '##INFO=<ID=PANEL_OF_NORMALS,Number=0,Type=Flag,Description="Overlap with germline call among panel of normals">'
        bcftools view -h "${PON}" | tail -n 1
        bcftools view -H "${PON}" | awk 'BEGIN{OFS="\t"} {
            if ($8 == ".") {
                $8 = "PANEL_OF_NORMALS"
            } else {
                $8 = $8 ";PANEL_OF_NORMALS"
            }
            print $0
        }'
    } | bgzip > "${PON_PCGR}" && tabix -p vcf "${PON_PCGR}"
    echo "PCGR-compatible PoN created: ${PON_PCGR}"
fi

# bcftools view -h "${PON_PCGR}" | grep -i "^##"
# bcftools view -h "${PON_PCGR}" | head -20
# bcftools view -h "${PON_PCGR}" | grep -i "^#CHROM"
# bcftools view -h "${PON}" | grep -i "^##"

## Output directory for PCGR results
# PCGR_DIR="${PROJECT_DIR}/data/wes/pcgr"
PCGR_DIR="${PROJECT_DIR}/data/pcgr_test"
mkdir -p "${PCGR_DIR}"

# tumour_id="DFSP-001-T"
# tumour_id="DFSP-028-T"
tumour_id="DFSP-035-T-P1"
# tumour_id="DFSP-035-T"

## Function to run PCGR annotation
pcgr_annotation() {
    local tumour_id=$1
    local MUTECT2_DIR=$2
    local CNA_DIR=$3
    local rna_dir=$4
    local BAM_DIR=$5
    local PCGR_REFERENCE_DIR=$6
    local VEP_CACHE_DIR=$7
    local PCGR_DIR=$8
    local CONTAINER_DIR=${9}
    local PON_DIR=${10}
    local PON_PCGR=${11}

    ## Find the matched normal sample
    patient_id=$(echo "${tumour_id}" | cut -d'-' -f1,2)
    normal_id="${patient_id}-N"

    echo "$(date +"%F") $(date +"%T")" "Processing sample: ${tumour_id} for patient: ${patient_id}"

    ## Check if normal sample exists and set is_paired flag
    if [ ! -d "${BAM_DIR}/${normal_id}" ]; then
        is_paired="false"
    else
        is_paired="true"
    fi

    echo "${is_paired}"
    ## Check if input VCF exists
    input_vcf="${MUTECT2_DIR}/${tumour_id}/${tumour_id}.final.vcf.gz"

    if [ ! -f "${input_vcf}" ]; then
        echo "Input VCF file not found: ${input_vcf}"
        return 1
    fi
    
    input_cna="${CNA_DIR}/${tumour_id}/${tumour_id}.pcgr.tsv"
    
    if [ ! -f "${input_cna}" ]; then
        echo "Input CNA file not found: ${input_cna}"
        return 1
    fi

    ## Create output directory for this sample
    output_dir="${PCGR_DIR}/${tumour_id}"
    mkdir -p "${output_dir}"

    ## Define reformatted VCF path
    reformatted_vcf="${output_dir}/${tumour_id}.reformatted.vcf.gz"

    ## ========================================================================
    ## Tumour-Normal PCGR Annotation
    ## ========================================================================
    if [ "${is_paired}" = "true" ]; then

        echo "$(date +"%F") $(date +"%T")" "Found paired normal sample ${normal_id} for ${tumour_id}. Running tumour-normal analysis ..."

        ## Reformat VCF 
        singularity exec \
            --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "/tmp:/tmp" \
            "${CONTAINER_DIR}/pysam.sif" \
            python "${MODULE_DIR}/pcgr_reformat_vcf_tumour_normal.py" \
                --input "${input_vcf}" \
                --output "${reformatted_vcf}"
        
        ## Index it
        rm -f "${reformatted_vcf}.tbi"
        singularity exec \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix -p vcf "${reformatted_vcf}"
        
        ## Run PCGR annotation
        singularity exec \
            --bind "${PCGR_REFERENCE_DIR}:${PCGR_REFERENCE_DIR}" \
            --bind "${VEP_CACHE_DIR}:${VEP_CACHE_DIR}" \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            "${CONTAINER_DIR}/pcgr.sif" \
            pcgr \
            --input_vcf "${reformatted_vcf}" \
            --input_cna "${input_cna}" \
            --n_copy_gain 4 \
            --cna_overlap_pct 50 \
            --vep_dir "${VEP_CACHE_DIR}" \
            --refdata_dir "${PCGR_REFERENCE_DIR}" \
            --output_dir "${output_dir}" \
            --genome_assembly grch38 \
            --sample_id "${tumour_id}" \
            --assay WES \
            --effective_target_size_mb 34 \
            --tumor_dp_tag TDP \
            --tumor_af_tag TAF \
            --control_dp_tag NDP \
            --control_af_tag NAF \
            --tumor_dp_min 20 \
            --tumor_af_min 0.05 \
            --control_dp_min 10 \
            --control_af_max 0.01 \
            --estimate_tmb \
            --tmb_dp_min 20 \
            --tmb_af_min 0.05 \
            --estimate_msi \
            --estimate_signatures \
            --vcf2maf \
            --ignore_noncoding \
            --force_overwrite \
            >& "${output_dir}/pcgr.log"

    ## ========================================================================
    ## Tumour-Only PCGR Annotation
    ## ========================================================================
    else
        
        echo "$(date +"%F") $(date +"%T")" "Found No-paired normal sample for ${tumour_id}. Running tumour-only analysis ..."
        
        ## Reformat VCF 
        singularity exec \
            --bind "${MUTECT2_DIR}:${MUTECT2_DIR}" \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${MODULE_DIR}:${MODULE_DIR}" \
            --bind "/tmp:/tmp" \
            "${CONTAINER_DIR}/pysam.sif" \
            python "${MODULE_DIR}/pcgr_reformat_vcf_tumour_only.py" \
                --input "${input_vcf}" \
                --output "${reformatted_vcf}"
        
        ## Index it
        rm -f "${reformatted_vcf}.tbi"
        singularity exec \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${output_dir}:${output_dir}" \
            "${CONTAINER_DIR}/tabix.sif" \
            tabix -p vcf "${reformatted_vcf}"
        
        ## Run PCGR annotation
        singularity exec \
            --bind "${PCGR_REFERENCE_DIR}:${PCGR_REFERENCE_DIR}" \
            --bind "${VEP_CACHE_DIR}:${VEP_CACHE_DIR}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${PCGR_DIR}:${PCGR_DIR}" \
            --bind "${PON_DIR}:${PON_DIR}" \
            --bind "/tmp:/tmp" \
            "${CONTAINER_DIR}/pcgr.sif" \
            pcgr \
                --input_vcf "${reformatted_vcf}" \
                --input_cna "${input_cna}" \
                --pon_vcf "${PON_PCGR}" \
                --VEP_CACHE_DIR "${VEP_CACHE_DIR}" \
                --refdata_dir "${PCGR_REFERENCE_DIR}" \
                --output_dir "$output_dir" \
                --genome_assembly grch38 \
                --sample_id "${tumour_id}" \
                --assay WES \
                --effective_target_size_mb 34 \
                --tumor_only \
                --tumor_dp_tag TDP \
                --tumor_af_tag TAF \
                --tumor_dp_min 20 \
                --tumor_af_min 0.05 \
                --estimate_tmb \
                --tmb_dp_min 20 \
                --tmb_af_min 0.05 \
                --estimate_msi \
                --estimate_signatures \
                --vcf2maf \
                --ignore_noncoding \
                --force_overwrite \
                >& "${output_dir}/pcgr.log"

    fi
}

export -f pcgr_annotation

## Run PCGR annotation in parallel
## Run PCGR annotation in parallel
cat ${sample_list} | parallel \
    --jobs "${jobs}" \
    --progress \
    -k \
    pcgr_annotation {} \
    "${MUTECT2_DIR}" \
    "${CNA_DIR}" \
    "${rna_dir}" \
    "${BAM_DIR}" \
    "${PCGR_REFERENCE_DIR}" \
    "${VEP_CACHE_DIR}" \
    "${PCGR_DIR}" \
    "${CONTAINER_DIR}" \
    "${PON_DIR}" \
    "${PON_PCGR}" 

echo "$(date +"%F") $(date +"%T")" "PCGR annotation completed for all samples."