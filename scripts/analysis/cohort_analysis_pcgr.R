#!/usr/bin/env Rscript
##############################################################################
## Description: This script processes PCGR data and performs cohort analysis
## Author:  Zhong Guorui
## Created: 2025-06-15
##############################################################################

## "==========================================================================="
## General configurations ----
## "==========================================================================="
## Load required libraries and functions
library(maftools)
source(here::here("scripts/lib/study_lib.R"))

## Output directories
out_data_dir <- "data/wes/collect"
out_plot_dir <- "figures/wes"
capture_size <- 34 

## "==========================================================================="
## Aggregate PCGR data --------------
## "==========================================================================="

## Merge all PCGR outputs
pcgr_dir <- "data/wes/pcgr"
sheet_name <- "SOMATIC_SNV_INDEL"
pcgr_data <- collectPCGRData(dir = pcgr_dir, sheet_name = sheet_name)

total_variants <- nrow(pcgr_data)
message(paste("Total number of variants = ", total_variants))

## Save the collected PCGR data
saveData(pcgr_data, dir = out_data_dir, filename = "dfsp_cohort_pcgr_merged")

##"==========================================================================="
## Filter PCGR data --------------
##"==========================================================================="
pcgr_data <- loadData(dir = out_data_dir, filename = "dfsp_cohort_pcgr_merged")

## Explore the threshold values for filtering
# pcgr_data$oncogenicity |> unique()
# pcgr_data$actionability |> unique()
# pcgr_data |> 
#     filter(sample_id == "DFSP-001-T") |> 
#     filter(gnom_a_de_af <= 0.001) |>
#     # filter(oncogenicity  %in% c("Oncogenic", "Likely_Oncogenic")) |> 
#     filter(actionability %in% c("Potential significance"))
#     select(
#         sample_id, 
#         symbol,
#         dp_tumor, 
#         vaf_tumor, 
#         dp_control, 
#         vaf_control, 
#         gnom_a_de_af
#     ) |> 
#     filter(
#         dp_tumor >= 20, 
#         vaf_tumor >= 0.05, 
#         dp_control >= 10, 
#         vaf_control <= 0.01, 
#         gnom_a_de_af <= 0.001
#     ) 

## Total variants = 7960, some samples don't have variants
## Might be too stingent, so we will relax the thresholds
# final_variants <- filterPCGRData(
#     pcgr_data,
#     min_dp_tumor = 20,
#     min_vaf_tumor = 0.05,
#     min_dp_control = 10,
#     max_vaf_control = 0.01,
#     max_population_frequency = 0.001
# )

## Total variants = 14136, all samples have variants
final_variants <- filterPCGRData(
    pcgr_data,
    min_dp_tumor = 20,
    min_vaf_tumor = 0.03,
    min_dp_control = 10,
    max_vaf_control = 0.01,
    max_population_frequency = 0.001
)

## Include all indels, actionable SNVs, and oncogenic SNVs, 
## and predictedly pathogenic SNVs
## total_variants = 197717, all samples
# final_variants <- filterPCGRData(
#     pcgr_data,
#     min_dp_tumor = NULL,
#     min_vaf_tumor = NULL,
#     min_dp_control = NULL,
#     max_vaf_control = NULL,
#     max_population_frequency = 0.001
# )

## Convert the PCGR data to a MAF format
maf_tbl <- convertPCGRToMaftools(pcgr_tbl = final_variants) 
clinical_info <- loadDFSPClinicalInfo()
maf_obj <- read.maf(maf = maf_tbl, clinicalData = clinical_info)

## "==========================================================================="
## Oncoplots ----
## "==========================================================================="
dfsp_colors <- loadDFSPColorConfigs()

oncoplot(
    maf = maf_obj,
    top = 20,
    clinicalFeatures = c("FST.Group"),
    annotationColor = dfsp_colors,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE
)

## plot titv summary
dfsp_titv = titv(maf = maf_obj, plot = FALSE, useSyn = TRUE)
plotTiTv(res = dfsp_titv)

lollipopPlot(
    maf = maf_obj,
    gene = "BRCA2",
    AACol = "protein_change",
    showMutationRate = TRUE,
    repel = TRUE
)

## "==========================================================================="
## Compare mutation load against TCGA cohort ----
## "==========================================================================="
tcgaCompare(
    maf = maf_obj,
    cohortName = "FS-DFSP",
    logscale = TRUE, 
    capture_size = capture_size
)

plotVaf(
    maf = maf_obj,
    vafCol = "vaf_tumor"
)

## "==========================================================================="
## Tumor mutational burden  ----
## "==========================================================================="
tmb_data <- tmb(maf_obj, captureSize = 34, logScale = TRUE)
tmb_tbl <- tmb_data |> 
    left_join(clinical_info, by = c("Tumor_Sample_Barcode" = "Sample.ID")) |> 
    as_tibble()

group_columns <- c("FST.Group", "Main.Met", "Metastasis")

for (group_column in group_columns) {

    ## Compare FST groups for TMB
    tmb_res <- compareDFSPGroups(
        data = tmb_tbl,
        value_col = "total_perMB",
        group_col = group_column,
        stat_method = "wilcox_test",
        p_adjust_method = "BH",
        is_paired = FALSE
    )

    savePlot(
        plot = tmb_res$plot,
        width = 3.5,
        height = 3,
        filename = paste0("tmb_pcgr_", group_column, "_groups_comparison"),
        only_png = TRUE,
        dir = out_plot_dir
    )
}

## Within FST groups, compare metastasis vs non-metastasis
tmb_tbl_fs_dfsp <- tmb_tbl |> filter(FST.Group %in% c("FS-DFSP"))
group_columns <- c("Main.Met", "Metastasis")

for (group_column in group_columns) {

    ## Compare FST groups for TMB
    tmb_res <- compareDFSPGroups(
        data = tmb_tbl_fs_dfsp,
        value_col = "total_perMB",
        group_col = group_column,
        stat_method = "wilcox_test",
        p_adjust_method = "BH",
        is_paired = FALSE
    )

    savePlot(
        plot = tmb_res$plot,
        width = 3.5,
        height = 3,
        filename = paste0("tmb_pcgr_", group_column, "_groups_comparison_within_fs_dfsp"),
        only_png = TRUE,
        dir = out_plot_dir
    )
}

## "========================================================================="
## Find significant variants FST groups ----
## "========================================================================="
group_column <- "FST.Group"
maf_obj@clinical.data[[group_column]] |> unique()
groups <- unique(maf_obj@clinical.data[[group_column]])
group_pairs <- combn(groups, 2, simplify = FALSE)

clinicalEnrichment(
    maf = maf_obj,
    clinicalFeature = group_column, 
    minMut = 5,
    pathways = FALSE
)

group_column <- "FST.Group"
sample_groups <- loadDFSPSampleGroups()


## Subset MAF for each group
maf_list <- map(
    sample_groups, 
    function(group) {
        clinical_query <- paste0(
            group_column, 
            " %in% c('", paste(group, collapse = "', '"), "')"
        )
        
        subsetMaf(
            maf = maf_obj, 
            clinQuery = clinical_query
        )
    }
)

## Find significant variants in each group
enrich_res <- map(
    maf_list,
    \(maf) {
        clinicalEnrichment(
            maf = maf,
            clinicalFeature = group_column, 
            minMut = 5,
            pathways = FALSE
        )
    }
)

plotEnrichmentResults(
    enrich_res = enrich_res[["U-DFSP vs Pre-FST"]]
)

enrichment_res$groupwise_comparision |> 
    filter(p_value < 0.05) |> 
    arrange(p_value) |> 
    as_tibble()

enrichment_res$groupwise_comparision %>%
    mutate(log10p = -log10(p_value)) %>%
    ggplot(aes(x = log(OR), y = log10p, color = Group1)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text_repel(
        data = . %>% filter(p_value < 0.001),
        aes(label = Hugo_Symbol), size = 3
    ) +
    facet_wrap(~Group1) +
    labs(x = "log(Odds Ratio)", y = "-log10(p-value)") +
    theme_minimal()

top_genes <- enrichment_res$groupwise_comparision %>%
  group_by(Group1) %>%
  slice_min(p_value, n = 5) %>%
  pull(Hugo_Symbol) %>%
  unique()

# Create a binary matrix of mutations
plotmafSummary(
    maf_obj, 
    showBarcodes = FALSE
)

oncoplot(
    maf = maf_obj,
    clinicalFeatures = clinical_features,
    annotationColor = annotation_colors,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE
)
