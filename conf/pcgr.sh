
singularity pull --force --dir "${CONTAINER_DIR}" pcgr.sif oras://ghcr.io/sigven/pcgr:2.2.1.singularity
apptainer pull --force --dir "${CONTAINER_DIR}" oras://ghcr.io/sigven/pcgr:2.2.3.singularity