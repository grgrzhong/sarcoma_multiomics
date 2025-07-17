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
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSVA)
library(msigdbr)
library(maftools)

source(here::here("scripts/lib/study_lib.R"))

## Output directories
data_dir <- "data/WES/collect"
plot_dir <- "figures/wes"

## Capture size
capture_size <- 34 

## Clinical information
clinical_info <- LoadDFSPClinicalInfo()

# clinical_info |> 
#     select(
#         Case.ID, Sample.ID, Main, Main.Met, Metastasis
#     ) |> 
#     view()

## "==========================================================================="
## Collect PCGR CNV data --------------
## "==========================================================================="
## Merge all PCGR outputs
pcgr_dir <- "data/wes/PCGR"
pcgr_cna_data <- collectPCGRCNAData(
    dir = pcgr_dir, 
    sig_gain_thres = 4
)

## Paired samples
paired_samples <- loadDFSPWESSamples()[["paired"]]

## Filter the CNA data
pcgr_cna_data <- pcgr_cna_data |> 
    ## Exclude the undifined variant class
    filter(variant_class %in% c("gain", "homdel")) |> 
    filter(sample_id %in% paired_samples)

# Most frequently altered genes
gene_freq <- pcgr_cna_data %>%
    count(symbol, variant_class) %>%
    arrange(desc(n)) |> 
    filter(n >= 20)

head(gene_freq, 20)
unique(gene_freq$EVENT_TYPE)

# top amplified genes
gene_freq %>%
    filter(EVENT_TYPE == "gain") %>%
    top_n(20, n) %>%
    ggplot(aes(x = reorder(SYMBOL, n), y = n)) +
    geom_col() +
    coord_flip() +
    labs(title = "Top 20 Amplified Genes", x = "Gene", y = "Count")

## "==========================================================================="
## Collect PCGR SNV Indel data --------------
## "==========================================================================="
## Merge all PCGR outputs
pcgr_dir <- "data/WES/PCGR"
sheet_name <- "SOMATIC_SNV_INDEL"
pcgr_data <- collectPCGRSNVINDELData(dir = pcgr_dir, sheet_name = sheet_name)

total_variants <- nrow(pcgr_data)
message(paste("Total number of variants = ", total_variants))

## Save the collected PCGR data
SaveData(
    pcgr_data, 
    dir = data_dir, 
    filename = "dfsp_wes_pcgr_cohort_merged_tbl"
)

##"==========================================================================="
## Filter PCGR data --------------
##"==========================================================================="
pcgr_data <- LoadData(
    dir = data_dir, 
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

## Convert the PCGR data to a MAF format
maf_tbl <- convertPCGRToMaftools(pcgr_tbl = final_variants) 

## CNV facets data
facet_dir <- "/mnt/f/projects/250224_sarcoma_multiomics/data/wes/cnv_facets"
cnv_data <- collectCNVFacets(facet_dir = facet_dir)

## Create a MAF object
maf_obj <- read.maf(
    maf = maf_tbl, 
    clinicalData = clinical_info,
    # cnTable = cnv_data,
    verbose = TRUE
)
maf_obj@variant.classification.summary
# oncoplot(maf = maf_obj, top = 20)

SaveData(
    maf_obj, 
    dir = data_dir, 
    filename = "dfsp_wes_pcgr_cohort_maf_obj"
)

## "========================================================================="
## Mutated genes and variants summary  ----
## "========================================================================="
maf_obj <- LoadData(
    dir = data_dir, 
    filename = "dfsp_wes_pcgr_cohort_maf_obj"
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
file <- here("figures/wes/dfsp_wes_pcgr_cohort_maf_summary.png")
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
file <- here("figures/wes/dfsp_wes_pcgr_cohort_tmb_tcga_comparison.png")
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
file <- here("figures/wes/dfsp_wes_pcgr_cohort_oncogenic_pathways.png")
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
file <- here("figures/wes/dfsp_wes_pcgr_cohort_somatic_interactions.png")
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
file <- here("figures/wes/dfsp_wes_pcgr_cohort_oncoplot_sig_genes.png")
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
## Public cohorts of sarcoma from cBioPortal
cbioportal_sarcoma_genes <- "data/public/cBioPortal_sarcoma_mutated_genes.txt"
cbioportal_sarcoma_genes <- read.table(
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
cbioportal_sarcoma_genes <- as_tibble(cbioportal_sarcoma_genes) |> 
    setNames(colnames) |> 
    select(Gene, Freq) |>
    mutate(Freq = str_replace(Freq, "<", "")) |>
    mutate(Freq = str_replace(Freq, "%", "")) |>
    mutate(Freq = as.numeric(Freq)) |>
    arrange(desc(Freq)) |> 
    filter(Freq > 1)

## DFSP-peng_2022 (https://doi.org/10.1111/bjd.20976)
## smith_2025 (10.1016/j.modpat.2025.100737)
dfsp_public_cohorts <- list(
    cbioportal_sarcoma = list(
        title = "Top mutated genes (cBioPortal Sarcoma, Freq > 1%)",
        snv_indel = cbioportal_sarcoma_genes$Gene
    ),
    smith_2025 = list(
        title = "Top mutated genes (Smith et al. 2025 Reported, n=53)",
        snv_indel = c(
            "FANCC", "TERT", "CHEK2", 
            "KMT2D", "BRIP1", "ERCC2", "KMT2A", 
            "MGA", "MTOR", "PTEN", "RAD50",
            "RIT1", "TP53"
        )
    ),
    peng_2022 = list(
        title = "Top mutated genes (Peng et al. 2022 Reported, n=59)",
        snv_indel = c(
            "MUC6", "MUC4", "KMT2C", 
            "HERC2", "HLA-A", "PRSS3", "PDE4DIP", 
            "PRSS1", "RBMXL1", "BRCA1", "CTBP2",
            "HLA-B", "HLA-DQA2", "LILRA6", "OR8U1", 
            "PABPC1", "PCSK1", "ZXDB"
        )
    )
)

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
        paste0("figures/wes/dfsp_wes_pcgr_cohort_oncoplot_public_top.pdf")
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
        paste0("figures/wes/dfsp_wes_pcgr_cohort_oncoplot_public_", 
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

file <- here("figures/wes/dfsp_wes_pcgr_cohort_oncoplot.png") 
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
fgsea_sets <- loadMsigdbGeneSet() |> cleanGeneSetName()
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
