#!/bin/bash
#SBATCH --job-name=mutect2_call
#SBATCH --partition=amd
#SBATCH --time=96:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_Multiomics/slurm/%x_%j.out
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_Multiomics/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

## ===========================================================================
## Mutect2 call
## Authors: Zhong Guorui
## Date: 2025-06-10
## Description: This script was adapted from "M:\Scripts\DNA analysis\PCGR_annotate.sh"
## Added features:
##      1. Prallelization using GNU parallel;
##      2. Use singularity container to reproduce the same environment
## ===========================================================================

## ===========================================================================
## Enable and modify this if runing at local environment
## ===========================================================================
export project_dir="/mnt/f/projects/250224_DFSP_Multiomics"
export bam_dir="${project_dir}/data/wes/bam"
# export bam_dir="/mnt/m/WES/DFSP/BAM"
export ref_dir="/mnt/m/Reference"
export reference="${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export germline="${ref_dir}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
export annotation_file="${ref_dir}/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/"
export interval="${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
export pon="${ref_dir}/WES/DFSP/PON-Mutect/pon.vcf.gz"

## Define work directory
export work_dir="${project_dir}/data/wes/mutect2"
mkdir -p ${work_dir}

## Create header files needed for repeatmasker and blacklist annotation
echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" >${work_dir}/vcf.rm.header
echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" >${work_dir}/vcf.map.header

## tumour sample list
tumour_all_samples="${project_dir}/data/wes/sample_info/tumour_all_samples.txt"

sample_list=${tumour_all_samples}
cat "${sample_list}"

## Parallel jobs, modify this according to the available resources
num_sample=$(echo "${sample_list}" | wc -l)
if [ "${num_sample}" -ge 15 ]; then    
    jobs=15
else
    jobs=$num_sample
fi

echo "===================================================================="
echo "Project directory:    ${project_dir}"
echo "Reference directory:  ${ref_dir}"
echo "Bam directory:        ${bam_dir}"
echo "Work directory:       ${work_dir}"
echo "Reference:            ${reference}"
echo "Germline:             ${germline}"
echo "Annotation file:      ${annotation_file}"
echo "Interval:             ${interval}"
echo "Panel of Normal:      ${pon}"
echo "===================================================================="

## Function to run mutect2 and filtering
mutect2_call() {
    local tumour_id=$1
    local ref_dir=$2
    local bam_dir=$3
    local work_dir=$4
    local reference=$5
    local germline=$6
    local annotation_file=$7
    local interval=$8
    local pon=${9}

    # Extract patient_id from tumour_id (everything before the second hyphen)
    patient_id=$(echo "${tumour_id}" | cut -d'-' -f1,2)
    normal_id=${patient_id}-N

    echo "$(date +"%F") $(date +"%T")" "Processing sample: ${tumour_id} for patient: ${patient_id}"

    # Create output directory
    output_dir=${work_dir}/${tumour_id}
    mkdir -p "${output_dir}"

    ## ========================================================================
    ## Step1. Get Pileup Summaries
    ## ========================================================================
    echo "$(date +"%F")" "$(date +"%T")" "Step1. Getting Pileup Summaries of tumour sample ..."
    gatk GetPileupSummaries \
        -I "${bam_dir}/${tumour_id}/${tumour_id}_recalibrated.bam" \
        -V "${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz" \
        -L "${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz" \
        -O "${output_dir}/${tumour_id}.getpileupsummaries.table"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: GetPileupSummaries failed for ${tumour_id}" >&2
        return 1
    fi

    # Get Pileup Summaries if presence of paired normal samples
    if [ -d "${bam_dir}/${normal_id}" ]; then
        echo "$(date +"%F")" "$(date +"%T")" "Getting Pileup Summaries of normal sample ..."
        gatk GetPileupSummaries \
            -I "${bam_dir}/${normal_id}/${normal_id}_recalibrated.bam" \
            -V "${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz" \
            -L "${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz" \
            -O "${output_dir}/${normal_id}.getpileupsummaries.table"

        if [[ $? -ne 0 ]]; then
            echo "ERROR: GetPileupSummaries failed for ${normal_id}" >&2
            return 1
        fi
    fi

    ## ====================================================================
    ## Step2. Calculate contamination
    ## ====================================================================
    ## Calculate contamination based on tumour samples only
    if [ -d "${bam_dir}/${normal_id}" ]; then
        echo "$(date +"%F")" "$(date +"%T")" "Step2. Calculating contamination of paired samples ..."
        gatk CalculateContamination \
            -I "${output_dir}/${tumour_id}.getpileupsummaries.table" \
            -matched "${output_dir}/${normal_id}.getpileupsummaries.table" \
            -O "${output_dir}/${tumour_id}.contamination.table" \
            -segments "${output_dir}/${tumour_id}.segments.table"

        if [[ $? -ne 0 ]]; then
            echo "ERROR: CalculateContamination failed for ${tumour_id}_vs_${normal_id}" >&2
            return 1
        fi

    else
        # Calculate contamination based on tumour samples only
        echo "$(date +"%F")" "$(date +"%T")" "Calculating contamination of tumour-only samples ..."
        gatk CalculateContamination \
            -I "${output_dir}/${tumour_id}.getpileupsummaries.table" \
            -O "${output_dir}/${tumour_id}.contamination.table" \
            -segments "${output_dir}/${tumour_id}.segments.table"

        if [[ $? -ne 0 ]]; then
            echo "ERROR: CalculateContamination failed for ${tumour_id}" >&2
            return 1
        fi
    fi

    ## ====================================================================
    ## Step3. Mutect2 call variants
    ## ====================================================================
    # Call variant on paired tumour samples
    if [ -d "${bam_dir}/${normal_id}" ]; then
        echo "$(date +"%F")" "$(date +"%T")" "Step3. Calling somatic variants on paired samples ..."
        gatk --java-options -Xmx4g Mutect2 \
            -I "$bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam" \
            -I "$bam_dir/${normal_id}/${normal_id}_recalibrated.bam" \
            -normal "$normal_id" \
            -R "${reference}" \
            -L "${interval}" \
            --germline-resource "${germline}" \
            --panel-of-normals "${pon}" \
            --f1r2-tar-gz "${output_dir}/${tumour_id}.f1r2.tar.gz" \
            --native-pair-hmm-threads 8 \
            --callable-depth 20 \
            -O "${output_dir}/${tumour_id}.mutect2.vcf.gz" \
            -bamout "${output_dir}/${tumour_id}.realigned.bam" \
            >& "${output_dir}/${tumour_id}.mutect2.log"

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Mutect2 failed for ${tumour_id}_vs_${normal_id}" >&2
            return 1
        fi

    else
        ## Mutect2 call variant on unpaired tumour samples
        echo "$(date +"%F")" "$(date +"%T")" "Step3. Calling somatic variants on unpaired samples ..."
        gatk --java-options -Xmx8g Mutect2 \
            -I "${bam_dir}/${tumour_id}/${tumour_id}_recalibrated.bam" \
            -R "${reference}" \
            -L "${interval}" \
            --germline-resource "${germline}" \
            --panel-of-normals "${pon}" \
            --f1r2-tar-gz "${output_dir}/${tumour_id}.f1r2.tar.gz" \
            --callable-depth 20 \
            -O "${output_dir}/${tumour_id}.mutect2.vcf.gz" \
            -bamout "${output_dir}/${tumour_id}.realigned.bam" \
            >& "${output_dir}/${tumour_id}.mutect2.log"

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Mutect2 failed for ${tumour_id}" >&2
            return 1
        fi
    fi

    ## ====================================================================
    ## Step4. Learn Read Orientation Model bias for artifacts
    ## ====================================================================
    echo "$(date +"%F")" "$(date +"%T")" "Step4. Learning Read Orientation Model ..."
    gatk LearnReadOrientationModel \
        -I "${output_dir}/${tumour_id}.f1r2.tar.gz" \
        -O "${output_dir}/${tumour_id}.readorientationmodel.tar.gz"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: LearnReadOrientationModel failed for ${tumour_id}" >&2
        return 1
    fi

    ## ====================================================================
    ## Step5. Learn Read Orientation Model bias for artifacts
    ## ====================================================================
    echo "$(date +"%F")" "$(date +"%T")" "Step5. Filtering Mutect calls ..."
    gatk --java-options -Xmx4g FilterMutectCalls \
        --variant "${output_dir}/${tumour_id}.mutect2.vcf.gz" \
        --stats "${output_dir}/${tumour_id}.mutect2.vcf.gz.stats" \
        --reference "${reference}" \
        --ob-priors "${output_dir}/${tumour_id}.readorientationmodel.tar.gz" \
        --contamination-table "${output_dir}/${tumour_id}.contamination.table" \
        --tumor-segmentation "${output_dir}/${tumour_id}.segments.table" \
        --min-allele-fraction 0.01 \
        --unique-alt-read-count 1 \
        --output "${output_dir}/${tumour_id}.filtermutectcalls.vcf.gz" \
        >& "${output_dir}/${tumour_id}.filtermutectcalls.log"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: FilterMutectCalls failed for ${tumour_id}" >&2
        return 1
    fi

    ## ====================================================================
    ## Step6. Normalize the variants
    ## ====================================================================
    ## Standardize the vcf, split multiallelic variants, and left align indels
    echo "$(date +"%F")" "$(date +"%T")" "Step6. Normalizing Mutect calls ..."
    bcftools norm \
        "${output_dir}/${tumour_id}.filtermutectcalls.vcf.gz" \
        -m-both -f "${reference}" \
        -Oz \
        -o "${output_dir}/${tumour_id}.normalized.vcf.gz"

    ## ====================================================================
    ## Step7. Filter out low-quality or failed variant calls
    ## ====================================================================
    ## Keep variants that have passed all filters (low-quality or failed variant calls)
    echo "$(date +"%F")" "$(date +"%T")" "Step7. Filtering PASS variants ..."
    bcftools view \
        -f PASS "${output_dir}/${tumour_id}.normalized.vcf.gz" \
        -o "${output_dir}/${tumour_id}.passed.vcf.gz"

    ## ====================================================================
    ## Step8. Filter out blacklist and repeatmasker regions
    ## ====================================================================
    ## Annotate repeatmasker and blacklist regions
    echo "$(date +"%F")" "$(date +"%T")" "Step8. Annotating repeatmasker regions ..."
    bcftools annotate \
        "${output_dir}/${tumour_id}.passed.vcf.gz" \
        --header-lines "${work_dir}/vcf.rm.header" \
        --annotations "${ref_dir}/RepeatMasker.bed.gz" \
        --columns CHROM,FROM,TO,RepeatMasker \
        --output "${output_dir}/${tumour_id}.repeatmasker.vcf.gz"

    echo "$(date +"%F")" "$(date +"%T")" "Step8 Annotating blacklist regions ..."
    bcftools annotate \
        "${output_dir}/${tumour_id}.repeatmasker.vcf.gz" \
        --header-lines "${work_dir}/vcf.map.header" \
        --annotations "${ref_dir}/blacklist.bed.gz" \
        --columns CHROM,FROM,TO,EncodeDacMapability \
        --output-type z \
        --output "${output_dir}/${tumour_id}.blacklist.vcf.gz"

    ## Filter out variants in RepeatMasker or Mapability
    echo "$(date +"%F")" "$(date +"%T")" "Step8. Filtering RepeatMasker and blacklist regions ..."
    bcftools filter \
        "${output_dir}/${tumour_id}.blacklist.vcf.gz" \
        -e 'INFO/RepeatMasker != "." || INFO/EncodeDacMapability != "."' \
        -Oz \
        -o "${output_dir}/${tumour_id}.filtered.vcf.gz"

    ## Remove the header lines from the final VCF file
    bcftools annotate \
        --remove INFO/RepeatMasker,INFO/EncodeDacMapability \
        --output-type z \
        --output "${output_dir}/${tumour_id}.final.vcf.gz" \
        "${output_dir}/${tumour_id}.filtered.vcf.gz"

    # Index the final VCF file
    tabix -p vcf "${output_dir}/${tumour_id}.final.vcf.gz"

    # Clean up intermediate files
    echo "$(date +"%F")" "$(date +"%T")" "Cleaning up intermediate files ..."
    rm -f "${output_dir}/${tumour_id}.normalized.vcf.gz"
    rm -f "${output_dir}/${tumour_id}.passed.vcf.gz"
    rm -f "${output_dir}/${tumour_id}.repeatmasker.vcf.gz"
    rm -f "${output_dir}/${tumour_id}.blacklist.vcf.gz"
    rm -f "${output_dir}/${tumour_id}.filtered.vcf.gz"

    ## ====================================================================
    ## Step9. Annotate variants
    ## ====================================================================
    ## Annotate variants by Funcotator
    echo "$(date +"%F")" "$(date +"%T")" "Step9. Annotating variants with Funcotator ..."
    gatk Funcotator \
        -R "${reference}" \
        -V "${output_dir}/${tumour_id}.final.vcf.gz" \
        -O "${output_dir}/${tumour_id}.funcotator.maf.gz" \
        -L "${interval}" \
        --output-file-format MAF \
        --data-sources-path "${annotation_file}" \
        --ref-version hg38 \
        --remove-filtered-variants true \
        >& "${output_dir}/${tumour_id}.funcotator.log"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Funcotator failed for ${tumour_id}" >&2
        return 1
    fi

    # Extract relevant columns from the Funcotator output
    less -S "${output_dir}/${tumour_id}.funcotator.maf.gz" |
        grep -v "#" > "${output_dir}/${tumour_id}.funcotator.tsv"

    # Annotate variant by Annovar
    echo "$(date +"%F")" "$(date +"%T")" "Step9. Annotating variants with Annovar ..."
    perl "${ref_dir}/annovar/table_annovar.pl" \
        "${output_dir}/${tumour_id}.final.vcf.gz" \
        "${ref_dir}/annovar/humandb/" \
        -buildver hg38 \
        -out "${output_dir}/${tumour_id}" \
        -remove \
        -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
        -operation gx,r,f,f,f,f,f \
        -nastring . \
        -polish \
        -xreffile "${ref_dir}/annovar/example/gene_fullxref.txt" \
        --otherinfo \
        --vcfinput \
        >& "${output_dir}/${tumour_id}.annovar.log"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Annovar failed for ${tumour_id}" >&2
        return 1
    fi

    # Extract relevant columns from the Annovar output
    less -S "${output_dir}/${tumour_id}.hg38_multianno.txt" |
        awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' > \
            "${output_dir}/${tumour_id}.annovar.txt"

    echo "$(date +"%F")" "$(date +"%T")" "Completed mutect2 variant calling  sample = $tumour_id"
}

## Export function to make it available to GNU parallel
export -f mutect2_call

## Run the processing in parallel
cat ${sample_list} | parallel \
    --jobs "${jobs}" \
    --progress \
    -k \
    mutect2_call {} "${ref_dir}" "${bam_dir}" "${work_dir}" "${reference}" "${germline}" "${annotation_file}" "${interval}" "${pon}"

rm ${work_dir}/vcf.rm.header
rm ${work_dir}/vcf.map.header

echo "$(date +"%F")" "$(date +"%T")" "Finished mutect2 call for all samples."
