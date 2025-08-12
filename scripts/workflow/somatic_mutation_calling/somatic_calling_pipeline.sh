#!/bin/bash
#SBATCH --job-name=somatic_mutation_pipeline
#SBATCH --partition=amd
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/somatic_mutation_pipeline/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/somatic_mutation_pipeline/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Load configuration
source "$(dirname "${BASH_SOURCE[0]}")/conf/config.sh"

# Add better error handling
set -e  # Exit on any error
set -u  # Exit on undefined variables
set -o pipefail  # Exit if any command in pipeline fails

# Log the start time to calculate total time taken
start_time=$(date +"%F %T")

echo "======================================================================="
echo "Somatic mutation calling Workflow - Pipeline"
echo "======================================================================="

# Step 1: Preprocessing
echo "$(date +"%F") $(date +"%T") Step 1: Preprocessing and QC ..."
bash "${PROJECT_DIR}/scripts/workflow/somatic_mutation_calling/step_01_preprocess.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: Preprocessing steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 1: Preprocessing and QC (✓) "

# Step 2: BWA alignment
echo "$(date +"%F") $(date +"%T") Step 2: BWA-mem alignment ..."
bash "${PROJECT_DIR}/scripts/workflow/somatic_mutation_calling/step_02_bwa_alignment.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: BWA alignment steps failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 2: BWA alignment (✓) "


# log the time of completion
end_time=$(date +"%F %T")
echo "$(date +"%F") $(date +"%T") Pipeline completed at: $end_time"

# Calculate total time taken
start_time_sec=$(date -d "$start_time" +%s)
end_time_sec=$(date -d "$end_time" +%s)
total_time=$((end_time_sec - start_time_sec))
hours=$((total_time / 3600))
minutes=$(((total_time % 3600) / 60))
seconds=$((total_time % 60))

# Print the pipeline summary
echo "======================================================================="
echo "Somatic mutation calling Workflow - Running Summary"
echo "======================================================================="
echo "Start Time:               $start_time"
echo "End Time:                 $end_time"
echo "Total Time:               ${hours}h ${minutes}m ${seconds}s"
echo "Output directory:         $OUTPUT_DIR"