#!/bin/bash
#SBATCH --job-name=cnv_facets
#SBATCH --partition=amd
#SBATCH --time=4:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_Multiomics/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_Multiomics/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

## "==========================================================================="
## cnv_facet analysis (https://github.com/dariober/cnv_facets)
## Authors: Zhong Guorui
## Date: 2025-06-10
## Description: This script runs CNV FACETS analysis for paired samples if available.
## Features:
##      1. Uses GNU Parallel for efficient processing of multiple samples
##      2. Reformat the output segment data
## "==========================================================================="

## General configuration
export project_dir="/mnt/f/projects/250224_sarcoma_multiomics"
export module_dir="${project_dir}/modules"
export container_dir="${project_dir}/containers"
export bam_dir="${project_dir}/data/wes/bam"

export ref_dir="/mnt/m/Reference"
export dbsnp="${ref_dir}/dbSNP.vcf.gz"
export dbsnp_index="${ref_dir}/dbSNP.vcf.gz.tbi"
export annotation="${ref_dir}/Gencode/annotation_protein_coding.bed"

## Define the working directory for CNV FACETS output
work_dir="${project_dir}/data/wes/cnv_facets"
mkdir -p "${work_dir}"

## sample list
# sample_list="${project_dir}/data/wes/sample_info/tumour_all_samples.txt"
# sample_list="${project_dir}/data/wes/sample_info/tumour_only_samples.txt"
sample_list="${project_dir}/data/wes/sample_info/selected_samples.txt"

# sample_list=${tumour_normal_samples}
cat "${sample_list}"

## Parallel jobs, modify this according to the available resources
num_sample=$(wc -l < "${sample_list}")
echo "Number of samples = ${num_sample}"

if [ "${num_sample}" -ge 15 ]; then    
    jobs=15
else
    jobs=$num_sample
fi

echo ====================================================================
echo "work_dir:                ${work_dir}"
echo "ref_dir:                 ${ref_dir}"
echo "dbsnp:                   ${dbsnp}"
echo "dbsnp_index:             ${dbsnp_index}"
echo "bam_dir:                 ${bam_dir}"
echo "sample_list:             ${sample_list}"
echo "jobs:                    ${jobs}"
echo ====================================================================

## Function to run CNV FACETS
cnv_facets() {
    local tumour_id="$1"
    local ref_dir="$2"
    local bam_dir="$3"
    local work_dir="$4"
    local dbsnp="$5"
    local dbsnp_index="$6"

    ## Extract patient_id from tumour_id (everything before the second hyphen)
    patient_id=$(echo "$tumour_id" | cut -d'-' -f1,2)
    
    ## Use matched normal or defined normal sample
    normal_id=${patient_id}-N
    # normal_id=DFSP-336-N
    
    echo "$(date +"%F")" "$(date +"%T")" "Processing tumour = ${tumour_id}, normal = ${normal_id}"

    ## Temporary output directory (WSL can not create makeinfo)
    out_dir="/tmp/cnv_facets/${tumour_id}"
    rm -rf "${out_dir}"
    mkdir -p "${out_dir}"

    ## Check if presence of paired normal samples
    if [ -d "${bam_dir}/${normal_id}" ]; then
        
        ## ===================================================================
        ## Step1: Run the CNV FACETS analysis
        ## ===================================================================
        echo "$(date +"%F")" "$(date +"%T")" "Step1: Running CNV FACETS for tumour = ${tumour_id}, normal = ${normal_id}"
        ## Check the bam files for tumour and normal samples
        tumour_bam="${bam_dir}/${tumour_id}/${tumour_id}_recalibrated.bam"
        normal_bam="${bam_dir}/${normal_id}/${normal_id}_recalibrated.bam"
        
        if [ ! -f "${tumour_bam}" ]; then
            echo "Error: Tumour BAM file not found: ${tumour_bam}"
            return 1
        fi
        
        if [ ! -f "${normal_bam}" ]; then
            echo "Error: Normal BAM file not found: ${normal_bam}"
            return 1
        fi

        prefix="${out_dir}/${tumour_id}"

        ## Run FACETS
        ## Note: when runing the singularity, should not use the mounted path 
        singularity exec \
            --bind "${project_dir}:${project_dir}" \
            --bind "${ref_dir}:${ref_dir}" \
            --bind "${bam_dir}:${bam_dir}" \
            --bind /tmp:/tmp \
            "${container_dir}/cnv_facets-0.16.1.sif" \
            cnv_facets.R \
                --snp-tumour "${tumour_bam}" \
                --snp-normal "${normal_bam}" \
                --snp-vcf "${dbsnp}" \
                --snp-nprocs 4 \
                --out "${prefix}"
                # >& "${prefix}.facets.log"

        ## ===================================================================
        ## Step2: Convert the VCF output to TSV format
        ## ===================================================================
        echo "$(date +"%F")" "$(date +"%T")" "Step2: Converting VCF output to TSV for tumour = ${tumour_id}"
        ## Move the outputs to the work directory
        mv "${out_dir}" "${work_dir}/"
        
        ## Convert the output vcf to TSV format 
        ## (https://github.com/sigven/vcf2tsvpy)
        singularity exec \
            --bind "${project_dir}:${project_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind /tmp:/tmp \
            "${container_dir}/vcf2tsvpy-0.6.1.sif" \
            vcf2tsvpy \
                --input_vcf "${work_dir}/${tumour_id}/${tumour_id}.vcf.gz" \
                --out_tsv "${work_dir}/${tumour_id}/${tumour_id}.tmp.tsv"
        
        ## Remove the first header line starting with # from the TSV output
        output_tsv="${work_dir}/${tumour_id}/${tumour_id}.tsv"
        grep -v "^#" "${work_dir}/${tumour_id}/${tumour_id}.tmp.tsv" > "${output_tsv}"
        
        rm -rf "${work_dir}/${tumour_id}/${tumour_id}.tmp.tsv"

        ## ===================================================================
        ## Step3: Prepare the Copy number segments for gpcr input
        ## ===================================================================
        ## Note: coordinates must be one-based; 
        ## https://sigven.github.io/pcgr/articles/input.html
        echo "$(date +"%F")" "$(date +"%T")" "Step3: Preparing Copy number segments for PCGR input for tumour = ${tumour_id}"
        
        pcgr_tsv="${work_dir}/${tumour_id}/${tumour_id}.pcgr.tsv"
        singularity exec \
            --bind "${project_dir}:${project_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${module_dir}:${module_dir}" \
            --bind /tmp:/tmp \
            "${container_dir}/r-4.4.2.sif" \
            Rscript ${module_dir}/cnv_facets_export_segments_to_pcgr.R \
                --input "${output_tsv}" \
                --output "${pcgr_tsv}"

        ## ===================================================================
        ## Step4: Annotate the segments with gene information
        ## ===================================================================
        echo "$(date +"%F")" "$(date +"%T")" "Step4: Annotating segments with gene information for tumour = ${tumour_id}"

        ## Create the segment file
        vcf_file="${work_dir}/${tumour_id}/${tumour_id}.vcf.gz"
        seg_file="${work_dir}/${tumour_id}/${tumour_id}.seg"

        singularity exec \
            --bind "${project_dir}:${project_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${module_dir}:${module_dir}" \
            --bind /tmp:/tmp \
            "${container_dir}/r-4.4.2.sif" \
            Rscript ${module_dir}/cnv_facets_export_vcf_to_segment.R \
                --input "${vcf_file}" \
                --output "${seg_file}"

        ## Annotate the segments with gene information
        seg_header_file="${work_dir}/${tumour_id}/${tumour_id}.seg_header.txt"
        head -n 1 "${seg_file}" > "${seg_header_file}"

        ## Create the header for the annotation file
        annotation_header_file="${work_dir}/${tumour_id}/${tumour_id}.annotation_header.txt"
        echo -e "Chr\tStart\tEnd\tGene" > "${annotation_header_file}"
        ## Combine the headers
        header_file="${work_dir}/${tumour_id}/${tumour_id}.header.txt"
        paste "${seg_header_file}" "${annotation_header_file}" > "${header_file}"

        ## Remove the header and convert to BED format
        bed_file="${work_dir}/${tumour_id}/${tumour_id}.bed"
        tail -n +2 "${seg_file}" > "${bed_file}"
        
        ## Run bedtools intersect using the converted BED file
        intersect_tsv="${work_dir}/${tumour_id}/${tumour_id}.intersect.tsv"
        
        singularity exec \
            --bind "${project_dir}:${project_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${ref_dir}:${ref_dir}" \
            --bind "${module_dir}:${module_dir}" \
            --bind /tmp:/tmp \
            "${container_dir}/bedtools-2.31.1.sif" \
            bedtools intersect \
                -wa \
                -wb \
                -a "${bed_file}" \
                -b "${annotation}" \
                > "${intersect_tsv}"
        
        ## Create the final annotated file
        final_file="${work_dir}/${tumour_id}/${tumour_id}.annotated.tsv"
        cat "${header_file}" "${intersect_tsv}" > "${final_file}"
        
        ## Clean up intermediate files
        rm "${seg_header_file}" "${annotation_header_file}" "${header_file}" "${bed_file}" "${intersect_tsv}"
    fi
}

export -f cnv_facets

## Run FACETS for each tumour sample
cat ${sample_list} | parallel \
    --jobs "$jobs" \
    --progress \
    -k \
    cnv_facets {} "${ref_dir}" "${bam_dir}" "${work_dir}" "${dbsnp}" "${dbsnp_index}"

echo "$(date +"%F")" "$(date +"%T")" "Finished cnv facet analysis for all samples."