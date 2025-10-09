#!/usr/bin/env Rscript
##############################################################################
## Description: Find genomic drivers for FST and Metastasis
## Author:  Zhong Guorui
## Created: 2025-06-15
##############################################################################
## Load required functions and configurations
source(here::here("conf/study_conf.R"))

## "==========================================================================="
## Find genomic drivers contribute to FST/Metastatsis --------------
## "==========================================================================="
## Get the SNV/Indel table
pcgr_snv_indel_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_snv_indels_tbl"
)

## DFSP-087-T no Indels/SNV found
snv_indel_mat <- GetVariantMat(
    data = pcgr_snv_indel_tbl,
    variant_type = "symbol",
    variant_class = "variant_class",
    is_somatic_matched = FALSE # total sample = 161
)

snv_indel_freq <- GetVariantFreq(
    mat = snv_indel_mat,
    variant_type = "symbol"
)

message(
    "Top mutated genes:\n",
    paste(snv_indel_freq$symbol[1:30], collapse = ", ")
)

## Get the CNV tables
pcgr_cnv_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_cnv_tbl"
)

## DFSP-169-T, no CNV found
cn_cytoband_mat <- GetVariantMat(
    data = pcgr_cnv_tbl,
    variant_type = "cytoband",
    variant_class = "cn_status",
    is_somatic_matched = TRUE
)

cn_cytoband_freq <- GetVariantFreq(
    mat = cn_cytoband_mat,
    variant_type = "cytoband"
)

cn_gene_mat <- GetVariantMat(
    data = pcgr_cnv_tbl,
    variant_type = "symbol",
    variant_class = "cn_status",
    is_somatic_matched = TRUE
)

cn_gene_freq <- GetVariantFreq(
    mat = cn_gene_mat,
    variant_type = "symbol"
)

message(
    "Top CNV genes:\n",
    paste(cn_gene_freq$symbol[1:30], collapse = ", ")
)

## Get the CNV arm level data
ascets_arm_raw <- LoadData(
    dir = "data/processed/",
    filename = "wes_cnv_facets_ascets_arm_raw"
)

hg38_cytoband <- read_tsv(
    file = here(ascet_dir, "cytoband_coordinates_hg38.txt"),
    show_col_types = FALSE
)

## CNV arm matrix
ascets_arm_tbl <- map_dfr(
    ascets_arm_raw,
    \(sample) {
        df <- sample[["calls"]]
        df |> 
            pivot_longer(
                cols = -sample,
                names_to = "arm",
                values_to = "variant_class"
            ) |> 
            dplyr::rename(sample_id = sample) |> 
            dplyr::filter(
                !(variant_class %in% c(
                        ## insufficient data or coverage)
                        "NC",
                        ## low coverage, not enough data to make a reliable call
                        "LOWCOV",
                        ## no significant copy-number alteration
                        "NEUTRAL"
                    )
                )
            ) |> 
            ## Combine the chromosome and arm information
            left_join(
                hg38_cytoband |> select(arm, chrom),
                by = "arm"
            ) |> 
            mutate(
                arm = paste0("chr", chrom, ":", arm)
            ) |> 
            select(-chrom)
    }
)

## DFSP-042-T, no arm found
cn_arm_mat <- GetVariantMat(
    data = ascets_arm_tbl,
    variant_type = "arm",
    variant_class = "variant_class",
    is_somatic_matched = TRUE
)
cn_arm_mat[1:5, 1:5]

cn_arm_freq <- GetVariantFreq(
    mat = cn_arm_mat,
    variant_type = "arm"
)

## Find the genomic driver event for FST and metastasis
mat_list <- list(
    mutated_gene = snv_indel_mat,
    cnv_cytoband = cn_cytoband_mat,
    cnv_gene = cn_gene_mat,
    cnv_arm = cn_arm_mat
)

## Loop through different data type
stat_res_list <- list()

for (type in names(mat_list)) {

    message(paste0("\n *Processing ", type))

    if (type == "mutated_gene") {

        clinical_info <- LoadClinicalInfo(is_somatic_matched = FALSE)
        cohort_size <- nrow(clinical_info)

    } else {

        clinical_info <- LoadClinicalInfo(is_somatic_matched = TRUE)
        cohort_size <- nrow(clinical_info)
    }

    if (grepl("gene", type)) {
        
        variant_type <- "symbol"

    } else if (grepl("arm", type)) {

        variant_type <- "arm"

    } else {
        
        variant_type <- "cytoband"

    }

    ## Loop through different comparsions
    for (gc in names(group_comparisons)) {
        
        message(paste0(" - ", gc))
        
        group1 <- group_comparisons[[gc]][["group1"]]
        group2 <- group_comparisons[[gc]][["group2"]]

        if (gc == "Pre-FST_vs_Post-FST") {

            ## Paired  samples
            group1_samples <- clinical_info |>
                filter(Number.of.samples.per.patient >= 2) |>
                filter(FST.Group %in% group1) |>
                pull(Sample.ID)

            group2_samples <- clinical_info |>
                filter(Number.of.samples.per.patient >= 2) |>
                filter(FST.Group %in% group2) |>
                pull(Sample.ID)

        } else if (gc == "Primary_vs_Metastasis") {
            ## Get the case that have metastatsis
            case_ids <- clinical_info |>
                filter(Number.of.samples.per.patient >= 2) |>
                filter(grepl("Metastasis", Specimen.Nature)) |>
                pull(Case.ID)

            group1_samples <- clinical_info |>
                filter(Case.ID %in% case_ids) |>
                filter(Number.of.samples.per.patient >= 2) |>
                filter(Specimen.Nature %in% group1) |>
                pull(Sample.ID)

            group2_samples <- clinical_info |>
                filter(Case.ID %in% case_ids) |>
                filter(Number.of.samples.per.patient >= 2) |>
                filter(Specimen.Nature %in% group2) |>
                pull(Sample.ID)

        } else {
            group1_samples <- clinical_info |>
                filter(FST.Group %in% group1) |>
                pull(Sample.ID)

            group2_samples <- clinical_info |>
                filter(FST.Group %in% group2) |>
                pull(Sample.ID)
        }

        ## Mutated gene level analysis
        stat_res_list[["diff"]][[type]][[gc]] <- GetGroupStatDiff(
            mat = mat_list[[type]],
            group1_samples = group1_samples,
            group2_samples = group2_samples,
            variant_type = variant_type,
            group_comparison = gc,
            alternative = "two.sided",
            p_adjust_method = "BH"
        )
    }

    ## Save the shared genes/cytobands
    stat_res_list[["share"]][[type]] <- GetVariantFreq(
        mat = mat_list[[type]],
        variant_type = variant_type
    )
}

## No significant genes found if use the p_adj_thres
p_val_thres <- 0.05
p_adj_thres <- 0.1
min_altered_samples <- 3

stat_res_sig <- list()
for (type in names(stat_res_list$diff)) {

    stat_res_sig[[type]] <- map(
        stat_res_list$diff[[type]],
        ~ .x |>
            filter(fisher_p_val < p_val_thres) |>
            filter(fisher_p_adj < p_adj_thres) |>
            filter(
                group1_altered >= min_altered_samples |
                    group2_altered >= min_altered_samples
            )
    )

    output_dir <- "outputs/group_comparison"
    dir_create(here(output_dir))

    # ## Save signficant results
    # write_xlsx(
    #     stat_res_sig[[type]],
    #     path = here(
    #         output_dir,
    #         paste0("wes_group_comparison_stat_res_sig_", type, ".xlsx")
    #     )
        
    # )

    # ## Save all results
    # write_xlsx(
    #     stat_res_list$diff[[type]],
    #     path = here(
    #         output_dir,
    #         paste0("wes_group_comparison_stat_res_all_", type, ".xlsx")
    #     )
    # )
}

SaveData(
    stat_res_sig,
    dir = "data/processed",
    filename = "wes_pcgr_snv_indels_cnv_stat_res_sig"
)

# ## Convert the PCGR data to a MAF format
# maf_tbl <- ConvertPCGRToMaftools(pcgr_tbl = snv_indel_filtered)

# ## Create a MAF object
# maf_obj <- read.maf(
#     maf = maf_tbl, 
#     clinicalData = clinical_info,
#     # cnTable = cnv_data,
#     verbose = TRUE
# )
# maf_obj@variant.classification.summary
# # oncoplot(maf = maf_obj, top = 20)

# SaveData(
#     maf_obj, 
#     dir = "data/processed", 
#     filename = "wes_pcgr_DFSP_cohort_merged_snv_indels_obj"
# )

## "========================================================================="
## Figure1 ----
## "========================================================================="
## 
snv_indel_freq

## Load cBioportal sarcoma data
sarcoma_data <- LoadCBioPortalSarcomaData(
    mutate_gene_freq_thres = 0.01,
    cna_gene_freq_thres = 0.1,
    structural_variants_freq_thres = 0.01
)

## Other public DFSP studies
DFSP_study_data <- LoadPublicDFSPStudyData()

## Our study
snv_indel_freq |> 
    filter(freq > 0.05)


## "==========================================================================="
## CNV cytobands contribute to FST/Metastatsis ----
## "==========================================================================="
## Default to be somatic matched sampels n=147
clinical_info <- LoadClinicalInfo()
cohort_size <- nrow(clinical_info)

## Do the analysis for both somatic matched (147) and all samples (161)
gistic_dirs <- dir_ls(here("data/WES/GISTIC2"))
dir_names <- path_file(gistic_dirs)
gistic_dirs <- setNames(gistic_dirs, dir_names)
gistic_dir <- gistic_dirs[["somatic_matched"]]

## Find the recurrently altered cytobands in entire cohort
cn_data <- GetGisticCNData(
    gistic_dir = gistic_dir,
    group = "all_tumors",
    gistic_qval_thres = 0.5
)

cn_cytoband_tbl <- cn_data$cn_cytoband_tbl
cn_gene_tbl <- cn_data$cn_gene_tbl

cn_cytoband_freq <- cn_cytoband_tbl |> 
    pivot_longer(
        cols = -Unique_Name,
        names_to = "Sample.ID",
        values_to = "cn_status"
    ) |> 
    group_by(Unique_Name) |> 
    summarise(
        n_altered = sum(cn_status != ""),
        .groups = "drop"
    ) |> 
    arrange(desc(n_altered)) |> 
    mutate(freq = n_altered / cohort_size) |> 
    filter(freq > 0.3)

cn_gene_freq <- cn_gene_tbl |> 
    pivot_longer(
        cols = -symbol,
        names_to = "sample_id",
        values_to = "variant_class"
    ) |> 
    group_by(symbol) |> 
    summarise(
        n_altered = sum(variant_class != ""),
        .groups = "drop"
    ) |> 
    arrange(desc(n_altered)) |> 
    mutate(freq = n_altered / cohort_size) |> 
    filter(freq > 0.3)

## Find the cytoband driver event for FST and metastasis
cytoband_stat_res <- list()
cnv_gene_stat_res <- list()
for (gc in names(group_comparisons)) {

    message("Processing group comparison: ", gc)

    group1 <- group_comparisons[[gc]][["group1"]]
    group2 <- group_comparisons[[gc]][["group2"]]

    if (gc == "Pre-FST_vs_Post-FST") {

        ## Paired  samples
        group1_samples <- clinical_info |> 
            filter(Number.of.samples.per.patient >= 2) |>
            filter(FST.Group %in% group1) |> 
            pull(Sample.ID)
        
        group2_samples <- clinical_info |> 
            filter(Number.of.samples.per.patient >= 2) |>
            filter(FST.Group %in% group2) |> 
            pull(Sample.ID)

    } else if (gc == "Primary_vs_Metastasis") {
        
        ## Get the case that have metastatsis
        case_ids <- clinical_info |> 
            filter(Number.of.samples.per.patient >= 2) |> 
            filter(grepl("Metastasis", Specimen.Nature)) |> 
            pull(Case.ID)

        group1_samples <- clinical_info |> 
            filter(Case.ID %in% case_ids) |> 
            filter(Number.of.samples.per.patient >= 2) |>
            filter(Specimen.Nature %in% group1) |> 
            pull(Sample.ID)
        
        group2_samples <- clinical_info |> 
            filter(Case.ID %in% case_ids) |> 
            filter(Number.of.samples.per.patient >= 2) |>
            filter(Specimen.Nature %in% group2) |> 
            pull(Sample.ID)

    } else {

        group1_samples <- clinical_info |> 
            filter(FST.Group %in% group1) |> 
            pull(Sample.ID)
        
        group2_samples <- clinical_info |> 
            filter(FST.Group %in% group2) |> 
            pull(Sample.ID)
    }

    message(" - Performing cytoband analysis ... ")
    cytoband_stat_res[[gc]] <- GetGisticCytobandGroupStatRes(
        cn_cytoband_tbl = cn_cytoband_tbl,
        group1_samples = group1_samples,
        group2_samples = group2_samples,
        group_comparison = gc,
        alternative = "two.sided",
        p_adjust_method = "BH"
    )

    message(paste0(" - Performing CNV gene analysis ... \n"))
    cnv_gene_stat_res[[gc]] <- GetGisticCNVGeneGroupStatRes(
        cn_gene_tbl = cn_gene_tbl,
        group1_samples = group1_samples,
        group2_samples = group2_samples,
        group_comparison = gc,
        alternative = "two.sided",
        p_adjust_method = "BH"
    )
}

sig_cytoband <- map(
    cytoband_stat_res,
    ~ .x |> 
        filter(fisher_p_value < p_val_thres) |> 
        # filter(fisher_p_adj < p_adj_thres) |> 
        # filter(fisher_p_value < p_value_thres) |> 
        filter(
                group1_altered >= min_altered_samples |
                group2_altered >= min_altered_samples
        )
)

sig_cnv_gene <- map(
    cnv_gene_stat_res,
    ~ .x |> 
        filter(fisher_p_value < p_val_thres) |> 
        # filter(fisher_p_adj < p_adj_thres) |> 
        # filter(fisher_p_value < p_value_thres) |> 
        filter(
                group1_altered >= min_altered_samples |
                group2_altered >= min_altered_samples
        )
)

SaveData(
    cytoband_stat_res,
    dir = "data/processed",
    filename = "wes_cnvfacets_gistic2_group_comparisons_stat_res_all"
)

write_xlsx(
    cytoband_stat_res,
    path = here(
        "outputs",
        "wes_cnvfacets_gistic2_group_comparisons_stat_res_all.xlsx"
    )
)

## "========================================================================="
## Mutated genes and variants summary  ----
## "========================================================================="
maf_obj <- LoadData(
    dir = "data/processed", 
    filename = "wes_pcgr_DFSP_cohort_merged_snv_indels_obj"
)

## Get the summary of mutated genes
gene_summary <- getGeneSummary(maf_obj) |> as_tibble()
all_mut_genes <- gene_summary |> 
    pull(Hugo_Symbol)

message(
    "Top mutated genes:\n",
    paste(all_mut_genes[1:30], collapse = ", ")
)

## Mutated genes in the oncokb gene list
oncokb_genes <- LoadOncoKBGeneList()

oncokb_genes_summary <- gene_summary |> 
    filter(Hugo_Symbol %in% oncokb_genes)

message(
    "Total number of genes in OncoKB gene list: ", nrow(oncokb_genes_summary)
)

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
file <- here("figures/WES/pcgr/WES_DFSP_PCGR_cohort_maf_summary.png")
png(
    file = file, width = 8, height = 6.5,
    type = "cairo", units = "in", res = 300
)
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
## Estimate DFSP TMB  ----
## "==========================================================================="
file <- here("figures/WES/pcgr/WES_PCGR_DFSP_cohort_tmb.png")
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

tmb_data <- as_tibble(tmb_data)

## Save the TMB data
SaveData(
    tmb_data, 
    dir = "data/processed", 
    filename = "wes_pcgr_tmb_data"
)

# tmb_data <- LoadData(
#     dir = "data/processed", 
#     filename = "wes_pcgr_tmb_data"
# )


## "==========================================================================="
## Compare TMB (Epic meth groups) ----
## "==========================================================================="
file <- here("data/clinical/phenoData_DFSP_allTumour_CNV_Meth.csv")
meth_data <- read_csv(file = file, show_col_types = FALSE) |> 
    select(Sample.ID, Meth.Subtype.Main.2Class) |> 
    filter(!is.na(Meth.Subtype.Main.2Class))

meth_data <- meth_data |> 
    left_join(tmb_data, by = c("Sample.ID" = "Tumor_Sample_Barcode"))

plot <- BoxPlotCompareGroups(
    data = meth_data,
    x = "Meth.Subtype.Main.2Class",
    y = "total_perMB",
    y_lab = "TMB/mb",
    y_lim = NULL,
    stat_method = "wilcox_test",
    alternative = "two.sided"
)

SavePlot(
    plot = plot,
    width = 3.5,
    height = 2.5,
    dir = plot_dir,
    filename = "WES_PCGR_DFSP_cohort_tmb_EPIC_Meth.Subtype.Main.2Class",
)

## "==========================================================================="
## Compare TMB against TCGA cohort ----
## "==========================================================================="
file <- here(plot_dir, "WES_DFSP_PCGR_cohort_tmb_tcga_comparison.png")
png(
    file = file, width = 6, height = 4, 
    type = "cairo", units = "in", res = 300
)
tcgaCompare(
    maf = maf_obj,
    cohortName = "FS-DFSP",
    logscale = TRUE, 
    capture_size = capture_size
)
message("Saving TMB vs TCGA plot to: ", file)
dev.off()

plotVaf(
    maf = maf_obj,
    vafCol = "vaf_tumor"
)

## "========================================================================="
## Compare TMB across groups ----
## "========================================================================="
clinical_info <- LoadClinicalInfo()

## Group comparsions
tmb_tbl <- tmb_data |> 
    left_join(clinical_info, by = c("Tumor_Sample_Barcode" = "Sample.ID"))

sample_ids <- clinical_info$Sample.ID

group_columns <- c("FST.Group", "Main.Met", "Metastasis")

table(clinical_info$Main)
table(clinical_info$Main.Met)
table(clinical_info$Metastasis)
clinical_info |> 
    select(Case.ID, Sample.ID, Main, Main.Met, Metastasis) |>
    view()

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

    SavePlot(
        plot = tmb_res$plot,
        width = 3.5,
        height = 3,
        filename = paste0("tmb_pcgr_", group_column, "_groups_comparison"),
        only_png = TRUE,
        dir = plot_dir
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

    SavePlot(
        plot = tmb_res$plot,
        width = 3.5,
        height = 3,
        filename = paste0("tmb_pcgr_", group_column, "_groups_comparison_within_fs_dfsp"),
        only_png = TRUE,
        dir = plot_dir
    )
}

## "========================================================================="
## Known oncogenic singaling pathways ----
## "========================================================================="
file <- here("figures/WES/dfsp_wes_pcgr_cohort_oncogenic_pathways.png")
png(
    file = file, width = 8, height = 6, 
    type = "cairo", units = "in", res = 300
)
oncoplot(
    maf = maf_obj,
    pathways = "sigpw",
    titleText = "Enrichment of known oncogenic signaling pathways",
    gene_mar = 8,
    topBarData = NULL,
    fontSize = 0.8,
    topPathways = 10,
    collapsePathway = TRUE,
    showPct = TRUE
)
message("Saving oncogenic pathways plot to: ", file)
dev.off()

## "========================================================================="
## Exclusive or co-occurring set of genes ----
## "========================================================================="
file <- here("figures/WES/dfsp_wes_pcgr_cohort_somatic_interactions.png")
png(
    file = file, width = 8, height = 6, 
    type = "cairo", units = "in", res = 300
)
par(oma = c(0, 1, 0.5, 0))
somaticInteractions(
    maf = maf_obj,
    top = 30,
    pvalue = c(0.05, 0.1),
    nShiftSymbols = 2
)
message("Saving somatic interactions plot to: ", file)
dev.off()

## "========================================================================="
## Mutated genes among FST groups ----
## "========================================================================="
maf_obj <- LoadData(
    dir = data_dir, 
    filename = "dfsp_wes_pcgr_cohort_maf_obj"
)

## Load the sample groups to compare
sample_groups <- LoadDFSPSampleGroups()
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
            as_tibble() |> 
            mutate(group = group_column)

        enrich_res_sig[[comparsion]] <- sig_res

    }
}

write_xlsx(
    enrich_res_sig, 
    here("results/dfsp_wes_pcgr_cohort_enrichment_results.xlsx"),
)

## Visualize the number of significant mutated genes in each comparison
plot_data <- tibble(
    comparison = names(enrich_res_sig),
    group = map_chr(enrich_res_sig, ~ .x$group[1]),
    n_sig_genes = map_int(enrich_res_sig, ~ nrow(.x))
) |> 
    mutate(comparison = paste0(group, " | ", comparison)) |> 
    arrange(desc(n_sig_genes))

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

SavePlot(
    plot = plot,
    width = 6,
    height = 4,
    filename = "dfsp_wes_pcgr_cohort_enrichment_significant_genes_comparisons",
    only_png = TRUE,
    dir = plot_dir
)

## Extract gene lists from each comparison
gene_lists <- map(enrich_res_sig, ~ .x$Hugo_Symbol)
comparison_names <- names(gene_lists)

## Simple gene frequency table
gene_frequency <- enrich_res_sig |>
    bind_rows(.id = "comparison") |>
    count(Hugo_Symbol, sort = TRUE, name = "n_comparisons") |>
    left_join(
        enrich_res_sig |>
            bind_rows(.id = "comparison") |>
            group_by(Hugo_Symbol) |>
            summarise(
                comparisons = paste(comparison, collapse = " | "),
                groups = paste(unique(group), collapse = ", "),
                .groups = "drop"
            ),
        by = "Hugo_Symbol"
    ) |> 
    filter(n_comparisons > 1)

message("Genes appearing in multiple comparisons:")
print(gene_frequency)

write_xlsx(
    gene_frequency,
    here("results/dfsp_wes_pcgr_cohort_gene_overlap_summary.xlsx")
)

## Show specific overlapping genes
message("\nMost frequently overlapping genes:")
for (i in 1:min(10, nrow(gene_frequency))) {
    gene <- gene_frequency$Hugo_Symbol[i]
    n_comp <- gene_frequency$n_comparisons[i]
    comps <- gene_frequency$comparisons[i]
    
    message(paste0(i, ". ", gene, " (", n_comp, " comparisons): ", comps))
}

## "========================================================================="
## Oncoplots of significant genes in FST or metastasis----
## "========================================================================="
tmb_tbl <- tmb_data |>
    as_tibble() |> 
    dplyr::rename(tsb = Tumor_Sample_Barcode, tmb = total_perMB) |>
    dplyr::select(tsb, tmb) |> 
    as.data.frame()

## Genes enriched in FST groups or metastasis
file <- here("figures/WES/dfsp_wes_pcgr_cohort_oncoplot_sig_genes.png")
png(
    file = file, width = 12, height = 8,
    type = "cairo", units = "in", res = 300
)
oncoplot(
    maf = maf_obj,
    topBarData = tmb_tbl,
    clinicalFeatures = c("FST.Group", "Metastasis"),
    genes = gene_frequency$Hugo_Symbol,
    # top = 10,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE
)
message("Saving oncoplot of significant genes to: ", file)
dev.off()

## "========================================================================="
## Oncoplot public cohort mutated genes ----
## "========================================================================="
## Cbioportal sarcoma data
sarcoma_data <- LoadCBioPortalSarcomaData(
    mutate_gene_freq_thres = 0.01,
    cna_gene_freq_thres = 0.1,
    structural_variants_freq_thres = 0.01
)

## Other public DFSP studies
DFSP_study_data <- LoadPublicDFSPStudyData()

## Our study
snv_indel_freq |> 
    filter(freq > 0.05)

## Oncoplot the cytobands
cn_data
cn_freq

## Overlap genes between our cohort and smith_2025
all_mut_genes <- gene_summary |> 
    filter(MutatedSamples > 0) |> 
    pull(Hugo_Symbol)

recurrent_smith_2025 <- intersect(
    dfsp_public_cohorts$smith_2025$snv_indel, 
    all_mut_genes
)

if (length(recurrent_smith_2025) == 0) {
    message("No recurrent mutated genes found in public cohorts.")
} else {
    message(
        "Recurrent mutated genes found in current DFSP public cohorts (smith_2025): ",
        paste(recurrent_smith_2025, collapse = ", ")
    )
}

## Overlap genes between our cohort and peng_2022
recurrent_peng_2022 <- intersect(
    dfsp_public_cohorts$peng_2022$snv_indel, 
    all_mut_genes
)
if (length(recurrent_peng_2022) == 0) {
    message("No recurrent mutated genes found in public cohorts.")
} else {
    message(
        "Recurrent mutated genes found in current DFSP public cohorts (peng_2022): ",
        paste(recurrent_peng_2022, collapse = ", ")
    )
}

## Overlap genes between our cohort and cbioportal_sarcoma_genes
recurrent_cbioportal_sarcoma <- intersect(
    dfsp_public_cohorts$cbioportal_sarcoma$snv_indel, 
    all_mut_genes 
)
if (length(recurrent_cbioportal_sarcoma) == 0) {
    message("No recurrent mutated genes found in public cohorts.")
} else {
    message(
        "Recurrent mutated genes found in current DFSP public cohorts (cbioportal_sarcoma): ",
        paste(recurrent_cbioportal_sarcoma, collapse = ", ")
    )
}

## Plot the combined genes
show_genes <- c(
    gene_frequency$Hugo_Symbol, ## Overlap genes across FST groups
    dfsp_public_cohorts$smith_2025$snv_indel, ## Smith et al. 2025
    dfsp_public_cohorts$peng_2022$snv_indel, ## Peng et al. 2022
    dfsp_public_cohorts$cbioportal_sarcoma$snv_indel[1:20] ## Top 20 from cBioPortal
) |> 
    unique()

file <- here(
        paste0("figures/WES/dfsp_wes_pcgr_cohort_oncoplot_public_top.pdf")
    )
    pdf(file = file, width = 8, height = 18)
    # par(oma = c(0, 1, 0.5, 0))
    oncoplot(
        maf = maf_obj,
        topBarData = tmb_tbl,
        clinicalFeatures = c("FST.Group", "Metastasis"),
        genes = show_genes,
        sortByAnnotation = TRUE,
        titleText = "Mutated Genes in FST groups or in Public Cohorts",
        showTumorSampleBarcodes = FALSE,
        removeNonMutated = FALSE
    )
    message("Saving oncoplot to: ", file)
    dev.off()

for (i in names(dfsp_public_cohorts)) {
    
    cohort <- dfsp_public_cohorts[[i]]

    within_our_cohort_genes <- cohort$snv_indel[
        cohort$snv_indel %in% all_mut_genes] 
    
    message(
        paste0(
            "Mutated genes in ", cohort$title, 
            " that are also in our cohort: ", 
            paste(within_our_cohort_genes, collapse = ", ")
        )
    )

    file <- here(
        paste0("figures/WES/dfsp_wes_pcgr_cohort_oncoplot_public_", 
        i, ".png")
    )
    png(
        file = file, width = 8, height = 6,
        type = "cairo", units = "in", res = 300
    )
    par(oma = c(0, 1, 0.5, 0))
    oncoplot(
        maf = maf_obj,
        topBarData = tmb_tbl,
        clinicalFeatures = c("FST.Group", "Metastasis"),
        genes = cohort$snv_indel,
        # top = 20,
        sortByAnnotation = TRUE,
        titleText = cohort$title,
        showTumorSampleBarcodes = FALSE,
        removeNonMutated = FALSE
    )
    message("Saving oncoplot to: ", file)
    dev.off()
}

## "========================================================================="
## Compare freq across cohort for overlap genes ----
## "========================================================================="
our_cohort_freq <- gene_summary |> 
    mutate(mutation_frequency = MutatedSamples / 161) |> 
    select(Hugo_Symbol, mutation_frequency) |> 
    arrange(desc(mutation_frequency)) |> 
    mutate(cohort = "Our Cohort")

peng_2022_freq <- read_xlsx(
    here("data/public/bjd20976-sup-0013-tables2.xlsx"),
    sheet = "Somatic genes",
    skip = 1
) |> 
    janitor::clean_names() |> 
    select(mutated_genes, mutation_frequency) |> 
    filter(
        mutated_genes %in% recurrent_peng_2022
    ) |> 
    rename(
        Hugo_Symbol = mutated_genes
    ) |> 
    mutate(cohort = "Peng et al. 2022")

cbioportal_sarcoma_genes <- "data/public/cBioPortal_sarcoma_mutated_genes.txt"
cbioportal_sarcoma_data <- read.table(
    here(cbioportal_sarcoma_genes), 
    header = FALSE, 
    skip = 1,
    # sep = "\t", 
    stringsAsFactors = FALSE
)
colnames <- c(
    "Gene", "# Mut", "#", "Profiled Samples", "Freq", 
    "Is Cancer Gene (source: OncoKB)"
)

## Get the > 1% mutated genes
cbioportal_sarcoma_freq <- as_tibble(cbioportal_sarcoma_data) |> 
    setNames(colnames) |> 
    select(Gene, Freq) |>
    filter(Gene %in% recurrent_cbioportal_sarcoma) |>
    mutate(Freq = str_replace(Freq, "<", "")) |>
    mutate(Freq = str_replace(Freq, "%", "")) |>
    mutate(Freq = as.numeric(Freq)) |>
    arrange(desc(Freq)) |> 
    rename(
        Hugo_Symbol = Gene,
        mutation_frequency = Freq
    ) |>
    mutate(
        mutation_frequency = mutation_frequency / 100,
        cohort = "cBioPortal Sarcoma"
    )

plot <- cbioportal_sarcoma_freq |> 
    bind_rows(
        our_cohort_freq |> filter(
            Hugo_Symbol %in% cbioportal_sarcoma_freq$Hugo_Symbol
        )
    ) |> 
    mutate(
        Hugo_Symbol = factor(
            Hugo_Symbol, 
            levels = cbioportal_sarcoma_freq$Hugo_Symbol
        )
    ) |>
    # mutate(Hugo_Symbol = fct_reorder(
    #     Hugo_Symbol, 
    #     mutation_frequency, 
    #     .fun = max, 
    #     .desc = TRUE
    # )) |>
    ggplot(aes(
        x = Hugo_Symbol, 
        y = mutation_frequency, 
        fill = cohort
    )) +
    geom_col(position = "dodge") +
    labs(
        x = NULL, 
        y = "Mutation Frequency", 
        title = "Mutation Frequency of Commonly Mutated Genes"
    ) +
    plot_theme()+
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

filename <- "dfsp_wes_pcgr_cohort_mutation_frequency_comparison_cbioportal"
SavePlot(
    plot = plot,
    width = 8,
    height = 4,
    filename = filename,
    only_png = TRUE,
    dir = plot_dir
)

plot <- peng_2022_freq |> 
    bind_rows(
        our_cohort_freq |> filter(
            Hugo_Symbol %in% peng_2022_freq$Hugo_Symbol
        )
    ) |> 
    mutate(
        Hugo_Symbol = factor(
            Hugo_Symbol, 
            levels = peng_2022_freq$Hugo_Symbol
        )
    ) |>
    # mutate(Hugo_Symbol = fct_reorder(
    #     Hugo_Symbol, 
    #     mutation_frequency, 
    #     .fun = max, 
    #     .desc = TRUE
    # )) |>
    ggplot(aes(
        x = Hugo_Symbol, 
        y = mutation_frequency, 
        fill = cohort
    )) +
    geom_col(position = "dodge") +
    labs(
        x = NULL, 
        y = "Mutation Frequency", 
        title = "Mutation Frequency of Commonly Mutated Genes"
    ) +
    plot_theme()+
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

filename <- "dfsp_wes_pcgr_cohort_mutation_frequency_comparison_peng_2022"
SavePlot(
    plot = plot,
    width = 5,
    height = 4,
    filename = filename,
    only_png = TRUE,
    dir = plot_dir
)


## "==========================================================================="
## Oncoplots top mutated genes ----
## "==========================================================================="
tmb_tbl <- tmb_data |>
    as_tibble() |> 
    dplyr::rename(tsb = Tumor_Sample_Barcode, tmb = total_perMB) |>
    dplyr::select(tsb, tmb) |> 
    as.data.frame()

file <- here("figures/WES/dfsp_wes_pcgr_cohort_oncoplot.png") 
png(
    file = file, width = 12, height = 8,
    type = "cairo", units = "in", res = 300
)
oncoplot(
    maf = maf_obj,
    topBarData = tmb_tbl,
    clinicalFeatures = c("FST.Group", "Main.Met", "Metastasis"),
    # genes = all_sig_genes,
    top = 20,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE
)

message("Saving oncoplot to: ", file)
dev.off()

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

maf_obj <- LoadData(
    dir = data_dir,
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

## "==========================================================================="
## Mutated genes GSEA analysis ----
## "==========================================================================="
## Get the genes that are mutated in > 5% samples
mut_gene_list <- getGeneSummary(maf_obj) |> 
    as_tibble() |> 
    arrange(desc(MutatedSamples)) |> 
    mutate(mut_pct = MutatedSamples / 161) |> 
    filter(mut_pct > 0.05) |> 
    pull(Hugo_Symbol)

mut_gene_tbl <- bitr(
    geneID=mut_gene_list,
    fromType="SYMBOL",
    toType=c("ENTREZID", "ENSEMBL"),
    OrgDb="org.Hs.eg.db"
) |> 
    as_tibble() |> 
    filter(if_all(everything(), ~!is.na(.))) |> 
    mutate(ENTREZID = as.numeric(ENTREZID)) |> 
    dplyr::select(SYMBOL, ENTREZID) |>
    distinct()

## signature gene set
fgsea_sets <- LoadMsigdbGeneSet() |> cleanGeneSetName()
fgsea_sets <- split(fgsea_sets$gene_symbol, fgsea_sets$gs_name)

## GSEA multilevel
gsea_res <- fgseaMultilevel(
    fgsea_sets, 
    stats = deframe(mut_gene_tbl)
) |>
    as_tibble() |>
    arrange(desc(NES))

test <- gsea_res |> 
    # filter(padj < 0.05) |>
    filter(pval < 0.05)
view(test)

## "==========================================================================="
## Complextheatmap oncoplot ----
## "==========================================================================="
laml.maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml.clin = system.file("extdata", "tcga_laml_annot.tsv", package = "maftools")
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

matrix <- oncoplot(
    maf = laml,
    top = 20,
    writeMatrix = TRUE
)

# test <- mutCountMatrix(maf = laml)
# test[1:3, 1:3]

# class(test)

# test <- test |> as_tibble(rownames = "Hugo_Symbol")

# colnames(test)
# view(test)

matMut <- read.table("onco_matrix.txt", header = T, check.names = F, sep = "\t")
matMut[1:3, 1:3]

# Standardize mutation type names
matMut[matMut == "In-frame"] = "In_frame"
matMut[matMut == "Missense_Mutation"] = "Missense"
matMut[matMut == "Nonsense_Mutation"] = "Truncating"
matMut[matMut == "Frame_Shift_Del"] = "Truncating"
matMut[matMut == "Frame_Shift_Ins"] = "Truncating"
matMut[matMut == "In_Frame_Del"] = "In_frame"
matMut[matMut == "In_Frame_Ins"] = "In_frame"
matMut[matMut == "Splice_Site"] = "Truncating"

matMuttmp = matMut
matMuttmp$gene = row.names(matMuttmp)
mat_long <- reshape2::melt(matMuttmp, id.vars = "gene", value.name = "Variant_Classification")
levels(factor(mat_long$Variant_Classification))

pdata <- getClinicalData(laml)
pdata <- subset(pdata, pdata$Tumor_Sample_Barcode %in% colnames(matMut))
pdata = as.data.frame(pdata)
pdata$days_to_last_followup = ifelse(pdata$days_to_last_followup == "-Inf", 0, pdata$days_to_last_followup)
# 画图并去除无突变的样本和基因
pdata$days_to_last_followup = as.numeric(pdata$days_to_last_followup)
pdata$FAB_classification = factor(pdata$FAB_classification)
pdata$Overall_Survival_Status = factor(pdata$Overall_Survival_Status)
str(pdata)
## 'data.frame':	164 obs. of  4 variables:
##  $ Tumor_Sample_Barcode   : chr  "TCGA-AB-2802" "TCGA-AB-2804" "TCGA-AB-2805" "TCGA-AB-2806" ...
##  $ FAB_classification     : Factor w/ 8 levels "M0","M1","M2",..: 5 4 1 2 2 3 4 3 3 5 ...
##  $ days_to_last_followup  : num  365 2557 577 945 181 ...
##  $ Overall_Survival_Status: Factor w/ 2 levels "0","1": 2 1 2 2 2 1 2 2 2 2 ...

# Update matMut to match pdata samples
matMut <- matMut[, pdata$Tumor_Sample_Barcode]

# Check what mutation types are actually present in the data
print("Unique mutation types in matMut:")
print(unique(as.vector(as.matrix(matMut[rowSums(matMut != "0") > 0, ]))))

# Remove rows and columns with no mutations
matMut_filtered <- matMut[rowSums(matMut != "0") > 0, colSums(matMut != "0") > 0]
print(paste("Matrix dimensions after filtering:", nrow(matMut_filtered), "x", ncol(matMut_filtered)))

# Filter pdata to match the filtered matrix samples
pdata_filtered <- subset(pdata, pdata$Tumor_Sample_Barcode %in% colnames(matMut_filtered))

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "white", col = NA))
  },
  In_frame = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_frame"], col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),  
              gp = gpar(fill = col["Missense"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  },
  Truncating = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Truncating"], col = NA))
  }
)

heatmap_legend_param <- list(title = "Alternations", at = c("In_frame", "Missense",
    "Truncating", "Multi_Hit"), labels = c("In_frame", "Missense", "Truncating",
    "Multi_Hit"))

# 指定颜色, 调整颜色代码即可
col <- c(In_frame = "purple", Missense = "orange", Multi_Hit = "black", Truncating = "blue")
# 定义注释信息 自定义颜色 连续性变量设置颜色（外）
col_OS = colorRamp2(c(0, 973), c("white", "red"))

ha <- HeatmapAnnotation(
    OS = pdata_filtered$days_to_last_followup, 
    Status = pdata_filtered$Overall_Survival_Status,
    FAB_classification = pdata_filtered$FAB_classification, 
    col = list(OS = col_OS), 
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 7)
)

column_title <- "This is Oncoplot "

# Use the filtered matrix for oncoPrint
oncoPrint(matMut_filtered, alter_fun = alter_fun, col = col, 
          alter_fun_is_vectorized = FALSE,
          top_annotation = ha,
          column_title = column_title)

oncoPrint(matMut_filtered,
    bottom_annotation = ha, # 注释信息在底部
    #   top_annotation=top_annotation,
    # right_annotation=NULL,
    alter_fun = alter_fun,
    col = col,
    column_title = column_title,
    heatmap_legend_param = heatmap_legend_param,
    row_names_side = "left",
    pct_side = "right",
    # column_order=sample_order,
    #       column_split=3
    alter_fun_is_vectorized = FALSE
)

oncoplot_anno <- oncoPrint(matMut_filtered,
    bottom_annotation = ha, # 注释信息在底部
    #   top_annotation=top_annotation,
    # right_annotation=NULL,
    alter_fun = alter_fun,
    col = col,
    column_title = "",
    heatmap_legend_param = heatmap_legend_param,
    row_names_side = "left",
    pct_side = "right",
    # column_order=sample_order,
    #       column_split=3
    alter_fun_is_vectorized = FALSE
)
draw(oncoplot_anno, annotation_legend_side = "left", )