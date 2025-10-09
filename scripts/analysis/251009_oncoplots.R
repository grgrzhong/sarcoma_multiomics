#!/usr/bin/env Rscript
##############################################################################
## Description: Find genomic drivers for FST and Metastasis
## Author:  Zhong Guorui
## Created: 2025-06-15
##############################################################################
## Load required functions and configurations
source(here::here("conf/study_conf.R"))

library(RSQLite)
library(biomaRt)
library(clusterProfiler)
library(msigdbr)
library(ReactomePA)
library(KEGGREST)

## "==========================================================================="
## Load data -----
## "==========================================================================="
stat_res_list <- LoadData(
    dir = "data/processed",
    filename = "wes_snv_indels_cnv_stat_res"
)

stat_res_sig <- LoadData(
    dir = "data/processed",
    filename = "wes_snv_indels_cnv_stat_res_sig"
)

gistic_data <- LoadData(
    dir = "data/processed",
    filename = "wes_cnv_facets_gistic2_somatic_matched_all_tumors_qval0.5"
)

cnv_facet_gene_tbl <- LoadData(
    dir = "data/processed",
    filename = "cnv_facet_gene_tbl"
)

cosmic_gene_list <- LoadCOSMICGeneList()
oncokb_gene_list <- LoadOncoKBGeneList()
clinical_info <- LoadClinicalInfo()

## "==========================================================================="
## Figure1 -----
## "==========================================================================="
## SNV/Indels from PCGR annotation
pcgr_snv_indel_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_snv_indels_tbl"
)

pcgr_snv_indel_mat <- GetVariantMat(
    data = pcgr_snv_indel_tbl,
    variant_type = "symbol",
    variant_class = "variant_class",
    is_somatic_matched = FALSE
)

pcgr_snv_indel_freq <- GetVariantFreq(
    mat = pcgr_snv_indel_mat,
    variant_type = "symbol"
)

## Number of TERT mutations = 32
## TERT promoter mutations = 2 /32
n_tert_mut <- pcgr_snv_indel_tbl |> 
    filter(symbol == "TERT") |> 
    select(sample_id, symbol, consequence, alteration, protein_change) |>
    # filter(grepl("upstream", consequence)) |> 
    pull(sample_id) |> 
    unique() |> 
    length()

message(
    sprintf(
        "TERT mutated genes in our cohort: %s", 
        n_tert_mut
    )
)

## All mutated genes in our cohort
snv_indel_genes <- unique(pcgr_snv_indel_tbl$symbol)

## Top mutated genes in our cohort
snv_indel_genes <- pcgr_snv_indel_freq |> 
    filter(freq >=0.05) |>  
    pull(symbol)

## Gene mutated in high frequency in other DFSP papers
public_DFSP_study_data <- LoadPublicDFSPStudyData()
public_DFSP_snv_indel <- c(
    public_DFSP_study_data$smith_2025$mutated_genes,
    public_DFSP_study_data$peng_2022$mutated_genes
) |> unique()

public_DFSP_snv_indel <- public_DFSP_snv_indel[
    public_DFSP_snv_indel %in% snv_indel_genes
]

## Genes to show in oncoplot
snv_indel_genes <- unique(
    c(public_DFSP_snv_indel, snv_indel_genes)
)

# ## Group the genes by biological process
## One gene usually mapped to multiple biological processes
# ego <- enrichGO(
#     gene = snv_indel_genes,
#     OrgDb = org.Hs.eg.db,
#     keyType = "SYMBOL",
#     ont = "BP",
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     readable = TRUE
# )
# sim <- pairwise_termsim(ego)

## Manually group the genes according to their main biological process
pathway_gene_groups <- list(
    "Genome Maintenance" = c(
        "ATM", "ATR", "BRCA2", "BRIP1", "FANCA", "POLQ"
    ),
    "Oncogenes/TSG" = c(
        "TP53", "RB1", "NF1", "APC", "TERT", "NOTCH2", "ERBB4", "STAT5B"
    ),
    "Cell Cycle" = c(
        "STAG1", "STAG2", "SMC1A", "KNL1"
    ),
    "Transcription" = c(
        "EP300", "MED12", "NCOR1", "CREBBP", "SPEN", "ZFHX3"
    ),
    "Epigenetic" = c(
        "KMT2D", "KMT2C", "ATRX", "ARID1A", "KDM6A", "PBRM1", "SETD2", 
        "KAT6B", "TET1", "TET2", "TRIM24", "SMARCA1", "CHD5"
    ),
    "Cell Adhesion/Migration" = c(
        "FAT1", "FAT4", "ROBO2", "PLXNB3"
    ),
    "Signal Transduction" = c(
        "LRP1B", "PTPN13", "PTPRT", "PTPRB", "MAP3K13", "ARHGEF12", "WNK2"
    ),
    "Other/Unknown" = c(
        "CSMD3", "ZMYM3", "TTN", "CLTC", "DROSHA"
    )
)

## Covner to a tibble
pathway_gene_groups_tbl <- enframe(pathway_gene_groups, 
    name = "pathway", value = "genes"
) |> 
    unnest(cols = c(genes)) |> 
    left_join(
        pcgr_snv_indel_freq,
        by = c("genes" = "symbol")
    ) |> 
    mutate(
        pathway = factor(
            pathway, 
            levels = c(
                "Genome Maintenance", 
                "Oncogenes/TSG", 
                "Cell Cycle", 
                "Epigenetic", 
                "Cell Adhesion/Migration", 
                "Transcription", 
                "Signal Transduction", 
                "Other/Unknown"
            )
        )
    ) |> 
    group_by(pathway) |> 
    arrange(pathway, desc(freq), .by_group = TRUE)

snv_indel_genes_sort <- pathway_gene_groups_tbl$genes

snv_indel_mat <- pcgr_snv_indel_mat[snv_indel_genes_sort, ]

## CNV cytobands
gistic_cnv_cytoband_tbl <- stat_res_list$freq$gistic_cnv_cytoband |> 
    left_join(
        gistic_data$cn_cytoband_data,
        by =c("cytoband" = "Unique_Name")
    ) |> 
    mutate(
        cytoband_alteration = paste0(Cytoband, "_", Variant_Classification)
    ) |> 
    group_by(cytoband_alteration) |> 
    slice_max(freq, n = 1, with_ties = FALSE) |> 
    ungroup() |> 
    arrange(desc(freq)) |> 
    slice_head(n =30)

cnv_cytobands <- gistic_cnv_cytoband_tbl$cytoband
cnv_cytoband_mat <- gistic_data$cn_cytoband_mat[cnv_cytobands, ]
cnv_cytoband_mat[cnv_cytoband_mat == "Amp"] <- "AMP"
cnv_cytoband_mat[cnv_cytoband_mat == "Del"] <- "DEL"

## CNV genes
public_cnv_genes <- c(
    ## present in the recurrently amp/del cytoband in our cohort or reported 
    ## in the literature 
    "SPOP",
    "SPHK1",
    "ITGB4",
    "COL1A1",
    "FOXK2",
    "HLF",
    "RPS6KB1",
    "TBX2",
    "PPM1D",
    "SOX9",
    "CRKL",
    "BCR",
    "PDGFB",
    "MN1",

    ## 
    "TERT",
    "AKT1",
    "NFKBIA",
    "BTBD7",
    "CDKN2A",
    "CDKN2B"
)

public_cnv_genes %in% gistic_data$cn_gene_data$symbol
public_cnv_genes %in% stat_res_list$freq$cnv_gene_bedtools$symbol

high_freq_cnv_genes <- stat_res_list$freq$cnv_gene_bedtools |> 
    filter(alteration_type %in% c("AMP", "HOMDEL")) |> 
    filter(freq >=0.3) |> 
    pull(symbol)

## High level AMP/DEL only && well-known oncogenes and tumor suppressors
onco_cnv_genes <- stat_res_list$freq$cnv_gene_bedtools |> 
    filter(alteration_type %in% c("AMP", "HOMDEL")) |>
    ## Well-known oncogenes and tumor suppressors
    mutate(
        is_cosmic = if_else(
            symbol %in% cosmic_gene_list, TRUE, FALSE
        ),
        is_oncokb = if_else(
            symbol %in% oncokb_gene_list, TRUE, FALSE
        )
    ) |> 
    # filter(is_cosmic | is_oncokb) |>
    filter(is_cosmic) |>
    arrange(desc(freq)) |> 
    filter(freq >=0.2) |> 
    pull(symbol)

cnv_genes <- unique(
    c(public_cnv_genes, onco_cnv_genes)
)

cnv_gene_mat <- GetVariantMat(
    data = cnv_facet_gene_tbl |> 
        mutate(cytoband_symbol = paste(cytoband, symbol, sep = "_")),
    variant_type = "cytoband_symbol",
    variant_class = "cn_status"
)

cnv_gene_mat <- cnv_gene_mat[
    grepl(paste(cnv_genes, collapse = "|"), rownames(cnv_gene_mat)),
]

## Oncoplot
# plot_date <- Sys.Date() |> format("%Y%m%d")
plot_date <- "251009"

## Bottom annotations
sample_annotations <- c(
    "FST.Group", "Metastasis", "Specimen.Nature",
    "Age.cat", "Sex", "Site.Category.Broad",
    "Histology.subtype", "Molecular.subtype"
)

snv_indel_row_anno <- pathway_gene_groups_tbl |> 
        dplyr::select(genes, pathway) |> 
        column_to_rownames("genes")

## Barplots: TMB
tmb_data <- LoadData(
    dir = "data/processed", 
    filename = "wes_pcgr_tmb_data"
)

tmb_data <- tibble(sample_id = clinical_info$Sample.ID) |> 
    left_join(
        tmb_data,
        by = c("sample_id" = "Tumor_Sample_Barcode")
    ) |>
    dplyr::select(sample_id, total_perMB) |> 
    mutate(
        total_perMB = if_else(is.na(total_perMB), 0, total_perMB)
    ) |>
    column_to_rownames("sample_id")

msi_data <- clinical_info |> 
    select(Sample.ID, MSIsensor.Score) |> 
    column_to_rownames(var = "Sample.ID")

fga_data <- LoadData(
    dir = "data/processed",
    filename = "wes_cnv_facet_fga"
) |>
    dplyr::select(sample_id, fga) |>
    column_to_rownames(var = "sample_id")

hrd_data <- clinical_info |> 
    select(Sample.ID, HRDscore) |> 
    column_to_rownames(var = "Sample.ID")

VerticalCombinedOncoplot(
    snv_indel_mat = snv_indel_mat,
    cnv_gene_mat = cnv_gene_mat,
    cnv_cytoband_mat = cnv_cytoband_mat,
    snv_indel_row_sort = FALSE,
    snv_indel_row_anno = snv_indel_row_anno,
    snv_indel_top_anno = list(
        "TMB" = tmb_data,
        "MSI" = msi_data,
        "HRD" = hrd_data,
        "FGA" = fga_data
    ),
    pct_side = "right",
    main_type = "cnv_gene",
    sample_annotation = sample_annotations,
    column_title = NULL,
    text_size = 7,
    title_size = 9,
    width = 20,
    height = 22,
    dir = "figures/oncoplot",
    filename = paste0(plot_date, "_figure1_oncoplot", ".pdf")
)

## "==========================================================================="
## Figure2 -----
## "==========================================================================="