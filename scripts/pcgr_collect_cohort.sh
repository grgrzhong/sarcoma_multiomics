#!/bin/bash


pcgr_vcf=/mnt/f/projects/250224_DFSP_Multiomics/outputs/pcgr_test/DFSP-001-T/DFSP-001-T.pcgr.grch38.vcf.gz

bcftools view -h ${pcgr_vcf} | tail -n +30

