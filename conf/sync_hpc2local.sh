#!/bin/bash

#############################################################################
## Sync remote server folders to local using rsync
## Author: Zhong Guorui
## Date: 2025-04-15
#############################################################################

# Remote server details
remote_user="zhonggr"
remote_host="hpcio1" # hpc2021-io1.hku.hk 
# remote_base_dir="/lustre1/g/path_my/250224_sarcoma_multiomics"

# Define folders to sync (absolute paths)
# Format: "local_absolute_path:remote_absolute_path"
folders_to_sync=(
    "/mnt/m/Reference:/lustre1/g/path_my/Reference"
    # "/mnt/f/projects/250224_sarcoma_multiomics/data/benchmark:/lustre1/g/path_my/250224_sarcoma_multiomics/data/benchmark"
    # "/mnt/f/projects/250224_sarcoma_multiomics/container:/lustre1/g/path_my/250224_sarcoma_multiomics/containers"
    # Add more folders as needed in the same format
)

# Sync each folder
for folder_pair in "${folders_to_sync[@]}"; do
    # Split the pair into local and remote paths
    local_path="${folder_pair%%:*}"
    remote_path="${folder_pair#*:}"
    
    echo "Syncing ${remote_user}@${remote_host}:${remote_path} to ${local_path}"
    
    # Create local directory structure
    mkdir -p "${local_path}"
    
    # Sync the folder contents from remote to local (trailing slash on source copies contents)
    rsync -av --update --progress \
        "${remote_user}@${remote_host}:${remote_path}/" \
        "${local_path}/"
        
    echo "-----------------------------------"
done

# Remove --dry-run flag when you've confirmed it works correctly