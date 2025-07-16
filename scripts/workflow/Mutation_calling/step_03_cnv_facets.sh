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
export project_dir="/mnt/f/projects/sarcoma_multiomics"
export module_dir="${project_dir}/scripts/modules"
export container_dir="${project_dir}/containers"
export BAM_DIR="${project_dir}/data/wes/BAM"

export ref_dir="/mnt/f/Reference"
export dbsnp="${ref_dir}/dbSNP.vcf.gz"
export dbsnp_index="${ref_dir}/dbSNP.vcf.gz.tbi"
export annotation="${ref_dir}/Gencode/annotation_protein_coding.bed"

## Define the working directory for CNV FACETS output
export CNV_FACETS_DIR="${PROJECT_DIR}/data/wes/CNV/cnv_facets"
mkdir -p "${CNV_FACETS_DIR}"
export CONTAINER_DIR="/mnt/f/projects/sarcoma_multiomics/containers"

echo "CNV FACETS output directory: ${CNV_FACETS_DIR}"

## Find tumour samples (those containing 'T' in their name)


## Function to run CNV FACETS
cnv_facets() {
    
    local tumour_id="$1"

    ## Extract patient_id from tumour_id (everything before the second hyphen)
    patient_id=$(echo "$tumour_id" | cut -d'-' -f1,2)
    
    ## Use matched normal or defined normal sample
    normal_id=${patient_id}-N
    
    # normal_id=DFSP-336-N
    
    ## Temporary output directory (WSL can not create makeinfo)
    # out_dir="/tmp/cnv_facets/${tumour_id}"
    # rm -rf "${out_dir}"
    # mkdir -p "${out_dir}"

    ## Check if presence of paired normal samples
    # if [ -d "${BAM_DIR}/${normal_id}" ]; then
        
    #     echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Running CNV FACETS ..."
        
    #     ## Check the bam files for tumour and normal samples
    #     tumour_bam="${BAM_DIR}/${tumour_id}/${tumour_id}_recalibrated.bam"
    #     normal_bam="${BAM_DIR}/${normal_id}/${normal_id}_recalibrated.bam"
        
    #     if [ ! -f "${tumour_bam}" ]; then
    #         echo "Error: Tumour BAM file not found: ${tumour_bam}"
    #         return 1
    #     fi
        
    #     if [ ! -f "${normal_bam}" ]; then
    #         echo "Error: Normal BAM file not found: ${normal_bam}"
    #         return 1
    #     fi

    #     prefix="${out_dir}/${tumour_id}"

        # ## Run FACETS
        # ## Note: when runing the singularity, should not use the mounted path 
        # singularity exec \
        #     --bind "${project_dir}:${project_dir}" \
        #     --bind "${ref_dir}:${ref_dir}" \
        #     --bind "${BAM_DIR}:${BAM_DIR}" \
        #     --bind /tmp:/tmp \
        #     "${container_dir}/cnv_facets-0.16.1.sif" \
        #     cnv_facets.R \
        #         --snp-tumour "${tumour_bam}" \
        #         --snp-normal "${normal_bam}" \
        #         --snp-vcf "${dbsnp}" \
        #         --snp-nprocs 4 \
        #         --out "${prefix}"
        #         # >& "${prefix}.facets.log"
        
        # ## Convert the VCF output to TSV format
        # echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Converting VCF output to TSV format ..."

        # ## Move the outputs to the work directory
        # mv "${out_dir}" "${CNV_FACETS_DIR}/"
        
        # ## Convert the output vcf to TSV format 
        # ## (https://github.com/sigven/vcf2tsvpy)
        # singularity exec \
        #     --bind "${project_dir}:${project_dir}" \
        #     --bind "${CNV_FACETS_DIR}:${CNV_FACETS_DIR}" \
        #     --bind /tmp:/tmp \
        #     "${container_dir}/vcf2tsvpy-0.6.1.sif" \
        #     vcf2tsvpy \
        #         --input_vcf "${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.vcf.gz" \
        #         --out_tsv "${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.tmp.tsv"
        
        # ## Remove the first header line starting with # from the TSV output
        output_tsv="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.tsv"
        # grep -v "^#" "${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.tmp.tsv" > "${output_tsv}"
        
        # rm -rf "${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.tmp.tsv"

        # ## Prepare the Copy number segments for pgcr input
        # ## Note: coordinates must be one-based; 
        # ## https://sigven.github.io/pcgr/articles/input.html
        # echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Preparing Copy number segments for PCGR input ..."
        
        # pcgr_tsv="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.pcgr.tsv"
        # singularity exec \
        #     --bind "${project_dir}:${project_dir}" \
        #     --bind "${CNV_FACETS_DIR}:${CNV_FACETS_DIR}" \
        #     --bind "${module_dir}:${module_dir}" \
        #     --bind /tmp:/tmp \
        #     "${container_dir}/r-4.4.2.sif" \
        #     Rscript ${module_dir}/cnv_facets_export_segments_to_pcgr.R \
        #         --input "${output_tsv}" \
        #         --output "${pcgr_tsv}"

        ## Prepare the segments for GISTIC2 input
        echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Preparing segments for GISTIC2 input ..."
        
        gistic2_tsv="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.gistic2.tsv"

        singularity exec \
            --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
            --bind "${CNV_FACETS_DIR}:${CNV_FACETS_DIR}" \
            --bind /tmp:/tmp \
            "${CONTAINER_DIR}/r.sif" \
            Rscript ${MODULE_DIR}/cnv_facets_export_segments_to_gistic2.R \
                --input "${output_tsv}" \
                --output "${gistic2_tsv}"

        # ## Annotate the segments with gene information
        # echo "$(date +"%F") $(date +"%T") - (${tumour_id}) Annotating segments with gene information ..."

        # ## Create the segment file
        # vcf_file="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.vcf.gz"
        # seg_file="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.seg"

        # singularity exec \
        #     --bind "${project_dir}:${project_dir}" \
        #     --bind "${CNV_FACETS_DIR}:${CNV_FACETS_DIR}" \
        #     --bind "${module_dir}:${module_dir}" \
        #     --bind /tmp:/tmp \
        #     "${container_dir}/r-4.4.2.sif" \
        #     Rscript ${module_dir}/cnv_facets_export_vcf_to_segment.R \
        #         --input "${vcf_file}" \
        #         --output "${seg_file}"

        # ## Annotate the segments with gene information
        # seg_header_file="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.seg_header.txt"
        # head -n 1 "${seg_file}" > "${seg_header_file}"

        # ## Create the header for the annotation file
        # annotation_header_file="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.annotation_header.txt"
        # echo -e "Chr\tStart\tEnd\tGene" > "${annotation_header_file}"
        # ## Combine the headers
        # header_file="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.header.txt"
        # paste "${seg_header_file}" "${annotation_header_file}" > "${header_file}"

        # ## Remove the header and convert to BED format
        # bed_file="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.bed"
        # tail -n +2 "${seg_file}" > "${bed_file}"
        
        # ## Run bedtools intersect using the converted BED file
        # intersect_tsv="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.intersect.tsv"
        
        # singularity exec \
        #     --bind "${project_dir}:${project_dir}" \
        #     --bind "${CNV_FACETS_DIR}:${CNV_FACETS_DIR}" \
        #     --bind "${ref_dir}:${ref_dir}" \
        #     --bind "${module_dir}:${module_dir}" \
        #     --bind /tmp:/tmp \
        #     "${container_dir}/bedtools-2.31.1.sif" \
        #     bedtools intersect \
        #         -wa \
        #         -wb \
        #         -a "${bed_file}" \
        #         -b "${annotation}" \
        #         > "${intersect_tsv}"
        
        # ## Create the final annotated file
        # final_file="${CNV_FACETS_DIR}/${tumour_id}/${tumour_id}.annotated.tsv"
        # cat "${header_file}" "${intersect_tsv}" > "${final_file}"
        
        # ## Clean up intermediate files
        # rm "${seg_header_file}" "${annotation_header_file}" 
        # rm "${header_file}" "${bed_file}" "${intersect_tsv}"

    # fi
}

export -f cnv_facets

## Run FACETS for each tumour sample
tumour_ids=$(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | grep "T" | sort)

echo "${tumour_ids}" | parallel \
    --jobs 16 \
    cnv_facets {}
