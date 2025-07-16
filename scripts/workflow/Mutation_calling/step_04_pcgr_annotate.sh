#!/bin/bash

###############################################################################
## PCGR annotation reference (https://sigven.github.io/pcgr/index.html)
## Authors: Zhong Guorui
## Date: 2025-06-10
## Description: This script was adapted from "M:\Scripts\DNA analysis\PCGR_annotate.sh"
## Key features: 
##      1. Added prallelization to run PCGR annotation for multiple samples; 
##      2. Use singularity container to keep the environment consistent;
###############################################################################

## ===========================================================================
## Enable and modify this if runing at local environment
## ===========================================================================
export project_dir="/mnt/f/projects/sarcoma_multiomics"
export ref_dir="/mnt/m/Reference"

## ===========================================================================
## General configurations
## ===========================================================================
export module_dir="${project_dir}/scripts/modules"
export container_dir="${project_dir}/containers"
export bam_dir="${project_dir}/data/WES/BAM"
export vcf_dir="${project_dir}/data/WES/Mutect2"
export cna_dir="${project_dir}/data/WES/CNV/cnv_facets"

export ref_data_dir="${ref_dir}/PCGR_reference/20250314"
export vep_dir="${ref_dir}/VEP_cache"
export panel_of_normals_dir="/mnt/m/WES/DFSP/PON-Mutect"
export panel_of_normals_original="${panel_of_normals_dir}/pon.vcf.gz"
export panel_of_normals="${panel_of_normals_dir}/pon_pcgr.vcf.gz"

## Create PCGR-compatible PoN if it doesn't exist
if [ ! -f "${panel_of_normals}" ]; then
    echo "Creating PCGR-compatible PoN VCF..."
    { 
        bcftools view -h "${panel_of_normals_original}" | head -n -1
        echo '##INFO=<ID=PANEL_OF_NORMALS,Number=0,Type=Flag,Description="Overlap with germline call among panel of normals">'
        bcftools view -h "${panel_of_normals_original}" | tail -n 1
        bcftools view -H "${panel_of_normals_original}" | awk 'BEGIN{OFS="\t"} {
            if ($8 == ".") {
                $8 = "PANEL_OF_NORMALS"
            } else {
                $8 = $8 ";PANEL_OF_NORMALS"
            }
            print $0
        }'
    } | bgzip > "${panel_of_normals}" && tabix -p vcf "${panel_of_normals}"
    echo "PCGR-compatible PoN created: ${panel_of_normals}"
fi

# bcftools view -h "${panel_of_normals}" | grep -i "^##"
# bcftools view -h "${panel_of_normals}" | head -20
# bcftools view -h "${panel_of_normals}" | grep -i "^#CHROM"
# bcftools view -h "${panel_of_normals_original}" | grep -i "^##"

## Output directory for PCGR results
# work_dir="${project_dir}/data/wes/pcgr"
work_dir="${project_dir}/data/pcgr_test"
mkdir -p "${work_dir}"

## tumour sample list
# sample_list="${project_dir}/data/wes/sample_info/tumour_all_samples.txt"
sample_list="${project_dir}/data/wes/sample_info/selected_samples.txt"
# cat ${sample_list}

## Parallel jobs
num_sample=$(wc -l < "${sample_list}")

if [ "${num_sample}" -ge 15 ]; then    
    jobs=15
else
    jobs=$num_sample
fi

## Print the configuration
echo "===================================================================="
echo "project_dir:             ${project_dir}"
echo "work_dir:                ${work_dir}"
echo "vcf_dir:                 ${vcf_dir}"
echo "cna_dir:                 ${cna_dir}"
echo "bam_dir:                 ${bam_dir}"
echo "container_dir:           ${container_dir}"
echo "ref_data_dir:            ${ref_data_dir}"
echo "vep_dir:                 ${vep_dir}"
echo "panel_of_normals:        ${panel_of_normals}"
echo "Number of samples:       ${num_sample}"
echo "Number of parallel jobs: ${jobs}"
echo "===================================================================="

# tumour_id="DFSP-001-T"
# tumour_id="DFSP-028-T"
# tumour_id="DFSP-031-T"

## Function to run PCGR annotation
pcgr_annotation() {
    local tumour_id=$1
    local vcf_dir=$2
    local cna_dir=$3
    local rna_dir=$4
    local bam_dir=$5
    local ref_data_dir=$6
    local vep_dir=$7
    local work_dir=$8
    local container_dir=${9}
    local panel_of_normals_dir=${10}
    local panel_of_normals=${11}

    ## Find the matched normal sample
    patient_id=$(echo "${tumour_id}" | cut -d'-' -f1,2)
    normal_id="${patient_id}-N"

    echo "$(date +"%F") $(date +"%T")" "Processing sample: ${tumour_id} for patient: ${patient_id}"

    ## Check if normal sample exists and set is_paired flag
    if [ ! -d "${bam_dir}/${normal_id}" ]; then
        is_paired="false"
    else
        is_paired="true"
    fi

    echo "${is_paired}"
    ## Check if input VCF exists
    input_vcf="${vcf_dir}/${tumour_id}/${tumour_id}.final.vcf.gz"

    if [ ! -f "${input_vcf}" ]; then
        echo "Input VCF not found: ${input_vcf}"
        return 1
    fi
    
    input_cna="${cna_dir}/${tumour_id}/${tumour_id}.pcgr.tsv"
    if [ ! -f "${input_cna}" ]; then
        echo "Input CNA file not found: ${input_cna}"
        return 1
    fi

    ## Create output directory for this sample
    output_dir="${work_dir}/${tumour_id}"
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
            --bind "${vcf_dir}:${vcf_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${module_dir}:${module_dir}" \
            --bind "/tmp:/tmp" \
            "${container_dir}/pysam.sif" \
            python "${module_dir}/pcgr_reformat_vcf_tumour_normal.py" \
                --input "${input_vcf}" \
                --output "${reformatted_vcf}"
        
        ## Index it
        rm -f "${reformatted_vcf}.tbi"
        singularity exec \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            "${container_dir}/tabix.sif" \
            tabix -p vcf "${reformatted_vcf}"
        
        ## Run PCGR annotation
        singularity exec \
            --bind "${ref_data_dir}:${ref_data_dir}" \
            --bind "${vep_dir}:${vep_dir}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${work_dir}:${work_dir}" \
            "${container_dir}/pcgr.sif" \
            pcgr \
            --input_vcf "${reformatted_vcf}" \
            --input_cna "${input_cna}" \
            --vep_dir "${vep_dir}" \
            --refdata_dir "${ref_data_dir}" \
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
            --bind "${vcf_dir}:${vcf_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${module_dir}:${module_dir}" \
            --bind "/tmp:/tmp" \
            "${container_dir}/pysam.sif" \
            python "${module_dir}/pcgr_reformat_vcf_tumour_only.py" \
                --input "${input_vcf}" \
                --output "${reformatted_vcf}"
        
        ## Index it
        rm -f "${reformatted_vcf}.tbi"
        singularity exec \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            "${container_dir}/tabix.sif" \
            tabix -p vcf "${reformatted_vcf}"
        
        ## Run PCGR annotation
        singularity exec \
            --bind "${ref_data_dir}:${ref_data_dir}" \
            --bind "${vep_dir}:${vep_dir}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${panel_of_normals_dir}:${panel_of_normals_dir}" \
            --bind "/tmp:/tmp" \
            "${container_dir}/pcgr.sif" \
            pcgr \
                --input_vcf "${reformatted_vcf}" \
                --input_cna "${input_cna}" \
                --pon_vcf "${panel_of_normals}" \
                --vep_dir "${vep_dir}" \
                --refdata_dir "${ref_data_dir}" \
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
    "${vcf_dir}" \
    "${cna_dir}" \
    "${rna_dir}" \
    "${bam_dir}" \
    "${ref_data_dir}" \
    "${vep_dir}" \
    "${work_dir}" \
    "${container_dir}" \
    "${panel_of_normals_dir}" \
    "${panel_of_normals}" 

echo "$(date +"%F") $(date +"%T")" "PCGR annotation completed for all samples."