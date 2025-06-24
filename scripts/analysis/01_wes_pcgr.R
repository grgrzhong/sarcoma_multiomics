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
out_data_dir <- "data/collect"
out_plot_dir <- "figures/wes"
capture_size <- 34 

## Plot parameters
text_size <- 8
title_size <- 9

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
saveData(
    pcgr_data, 
    dir = out_data_dir, 
    filename = "dfsp_wes_pcgr_cohort_merged_tbl"
)

##"==========================================================================="
## Filter PCGR data --------------
##"==========================================================================="
pcgr_data <- loadData(
    dir = out_data_dir, 
    filename = "dfsp_wes_pcgr_cohort_merged_tbl"
)

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
# pcgr_data |> 
#     filter(dp_tumor >= min_dp_tumor & vaf_tumor >= min_vaf_tumor) |>
#     filter(dp_control >= min_dp_control & vaf_control <= max_vaf_control) |> 
#     filter(is.na(gnom_a_de_af) | gnom_a_de_af <= max_population_frequency)


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

table(final_variants$variant_class)
## Convert the PCGR data to a MAF format

maf_tbl <- convertPCGRToMaftools(pcgr_tbl = final_variants) 
clinical_info <- loadDFSPClinicalInfo()
maf_obj <- read.maf(maf = maf_tbl, clinicalData = clinical_info)
saveData(
    maf_obj, 
    dir = out_data_dir, 
    filename = "dfsp_wes_pcgr_cohort_maf_obj"
)

## "========================================================================="
## Mutated genes and variants summary  ----
## "========================================================================="
maf_obj <- loadData(
    dir = out_data_dir, 
    filename = "dfsp_wes_pcgr_cohort_maf_obj"
)

## Get the summary of mutated genes
gene_summary <- getGeneSummary(maf_obj) |> as_tibble()

## Mutated genes in the oncokb gene list
oncokb_genes <- loadOncoKBGeneList()

oncokb_genes_summary <- gene_summary |> 
    filter(Hugo_Symbol %in% oncokb_genes)
message("Total number of genes in OncoKB gene list: ", nrow(oncokb_genes_summary))

## Top genes affected
top_genes <- gene_summary |> head(20)
message(
    "Top 20 genes affected by variants:",
    paste(top_genes$Hugo_Symbol, collapse = ", ")
)

## Greater than 5% samples have mutations, n=161
n_gens_5pct <- gene_summary |> 
    filter(MutatedSamples > 8) |> 
    pull(Hugo_Symbol) |> 
    unique() |> 
    sort()

message(
    paste(
        "Genes with mutations in more than 5% of samples (n = ",
        length(n_gens_5pct), "):\n"
        # paste(n_gens_5pct, collapse = ", ")
    )
)

## Barplots of mutated genes
mafbarplot(
    maf_obj,
    n = 20,
    genes = NULL,
    fontSize = 0.6,
    includeCN = FALSE,
    legendfontSize = 1,
    borderCol = "#34495e",
    showPct = TRUE
)

## Plot the summary of variants
file <- here("figures/wes/dfsp_wes_pcgr_cohort_maf_summary.pdf")
pdf(file = file, width = 8, height = 6.5)
plotmafSummary(
    maf = maf_obj,
    rmOutlier = TRUE,
    addStat = "median",
    titvRaw = FALSE,
    fs = 0.8,
    # textSize = 0.3,
    # titleSize = c(0.6, 0.5),
    top = 20,
    dashboard = TRUE
)
message("Saving plot to: ", file)
dev.off()

## Total somatic variants
variants_summary <- maf_obj@summary  |> as_tibble()
variants_total <- variants_summary |> 
    filter(ID == "total") |> 
    pull(summary)

message(paste("Total number of variants = ", variants_total))

## Median of mutations for each tumour
variants_mean <- variants_summary |> 
    filter(ID == "total") |> 
    pull(Mean)

message(paste("Mean number of mutations per tumour = ", variants_mean))

variants_median <- variants_summary |> 
    filter(ID == "total") |> 
    pull(Median)

message(paste("Median number of mutations per tumour = ", variants_median))

## "==========================================================================="
## Tumor mutational burden  ----
## "==========================================================================="
file <- here("figures/wes/dfsp_wes_pcgr_cohort_tmb.png")
png(
    file = file, width = 6, height = 5, 
    type = "cairo", units = "in", res = 300
)
tmb_data <- tmb(
    maf_obj, 
    captureSize = 34, 
    logScale = TRUE, 
    plotType = "classic", 
    verbose = TRUE
)

message("Saving TMB plot to: ", file)
dev.off()

## Median TMB
message(
    paste(
        "Median TMB = ", 
        round(median(tmb_data$total_perMB, na.rm = TRUE), 2)
    )
)

## "==========================================================================="
## Compare TMB against TCGA cohort ----
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

## "========================================================================="
## Compare TMB across groups ----
## "========================================================================="

## Group comparsions
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
## Signaling pathways affected ----
## "========================================================================="
oncoplot(
    maf = maf_obj,
    pathways = "sigpw",
    titleText = "Signaling pathways affected",
    gene_mar = 8,
    fontSize = 0.8,
    topPathways = 10,
    collapsePathway = TRUE,
    showPct = FALSE
)

## "========================================================================="
## Compare mutated genes across cohort ----
## "========================================================================="
## DFSP-WGS (https://doi.org/10.1111/bjd.20976)
# peng_2022 <- read_xlsx(
#     here("data/public/bjd20976-sup-0013-tables2.xlsx"),
#     sheet = "Somatic genes",
#     skip = 1
# ) |> 
# sl
pull("Mutated genes")

peng_2022 <- c(
    ""
)

## "========================================================================="
## Find significant variants enriched in FST ----
## "========================================================================="
maf_obj <- loadData(
    dir = out_data_dir, 
    filename = "dfsp_wes_pcgr_cohort_maf_obj"
)

## Load the sample groups to compare
sample_groups <- loadDFSPSampleGroups()
message(
    "Sample groups to compare:\n", 
    paste("-", names(sample_groups), collapse = "\n")
)

## Loop through each group and perform clinical enrichment
enrich_res_sig <- list()

for (group_column in names(sample_groups)) {
    
    for (group_comparsion in names(sample_groups[[group_column]])) {

        message(
            paste0(
                "Processing group comparsion:", 
                group_column, ":", group_comparsion
            )
        )

        group <- sample_groups[[group_column]][[group_comparsion]]

        clinical_query <- paste0(
            group_column, 
            " %in% c('", paste(group, collapse = "', '"), "')"
        )

        maf_subset <- subsetMaf(
            maf = maf_obj, 
            clinQuery = clinical_query
        )

        enrich_res <- clinicalEnrichment(
            maf = maf_subset,
            clinicalFeature = group_column, 
            minMut = 5,
            pathways = FALSE
        )

        ## Save the enrichment results with sample sizes
        comparsion <- enrich_res$cf_sizes |> 
            mutate(name = paste0(cf, "(", N, ")")) |>
            pull(name)
        
        comparsion <- paste(comparsion, collapse = " vs ")

        sig_res <- enrich_res$pairwise_comparision |>
            filter(fdr < 0.05) |> 
            arrange(fdr) |> 
            as_tibble()

        sheet_name <- paste0(group_column, " | ", comparsion)

        enrich_res_sig[[sheet_name]] <- sig_res

    }
}

write_xlsx(
    enrich_res_sig, 
    here("results/dfsp_wes_pcgr_cohort_enrichment_results.xlsx"),
)

## Visualize the enrichment results
plot_data <- tibble(
    comparison = names(enrich_res_sig),
    n_sig_genes = map_int(enrich_res_sig, ~ nrow(.x))
)

plot <- plot_data |> 
    ggplot(aes(x = n_sig_genes, y = reorder(comparison, n_sig_genes))) +
    geom_col(fill = "steelblue") +
    geom_text(
        aes(label = n_sig_genes), 
        hjust = -0.2, 
        color = "black", 
        size = 8,
        size.unit = "pt"
    ) +
    labs(
        x = "Number of significant genes (FDR < 0.05)",
        y = NULL,
        title = "Number of significant genes in each comparison"
    ) +
    theme_minimal()+
    plot_theme()

savePlot(
    plot = plot,
    width = 6,
    height = 4,
    filename = "dfsp_wes_pcgr_cohort_enrichment_significant_genes_comparisons",
    only_png = TRUE,
    dir = out_plot_dir
)

## "==========================================================================="
## Oncoplots ----
## "==========================================================================="
file <- here("figures/wes/dfsp_wes_pcgr_cohort_oncoplot.png") 
png(file = file, width = 12, height = 8)
oncoplot(
    maf = maf_obj,
    clinicalFeatures = c("FST.Group", "Main.Met", "Metastasis"),
    # genes = all_sig_genes,
    top = 20,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE
)
dev.off()

getGeneSummary(maf_obj) |> 
    arrange(desc(Mutation_Count)) |> 
    slice_head(n = 30) |> 
    select(Hugo_Symbol, Mutation_Count)

plotmafSummary(
    maf = maf_obj,
    rmOutlier = TRUE,
    addStat = "median",
    dashboard = TRUE,
    titvRaw = FALSE
)

## Collect all the significant genes
all_sig_genes <- enrich_res_sig |> 
    bind_rows(.id = "comparison") |> 
    pull(Hugo_Symbol) |> 
    unique()

## Order the sig_genes by mutation frequency
ranked_sig_genes <- getGeneSummary(maf_obj) |> 
    filter(Hugo_Symbol %in% all_sig_genes) |>
    slice_head(n = 30) |> 
    pull(Hugo_Symbol)

oncoplot(
    maf = maf_obj,
    clinicalFeatures = c(
        "FST.Group"
        # "Main.Met"
        # "Metastasis"
    ),
    bgCol = "gray90",
    borderCol = "white",
    writeMatrix = TRUE,
    genes = ranked_sig_genes,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE
)

maf_obj@gene.summary
maf_obj@summary
maf_obj@variant.type.summary

maf_obj <- loadData(
    dir = out_data_dir,
    filename = "dfsp_wes_pcgr_cohort_maf_obj"
)

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
## Oncoplot of commonly mutated genes in DFSP ----
## "==========================================================================="
dfsp_common_mut_genes <- c(
    # Primary driver
    "COL1A1", "PDGFB", "PDGFRB", 
    
    # Secondary alterations
    "CDKN2A", "TP53", "RB1", "NF1", 
    "NRAS", "KRAS", "PTEN", "TERT",
    "MYC", "EGFR", "BRAF", "PIK3CA", 
    "AKT1", "DDIT4L"
)



