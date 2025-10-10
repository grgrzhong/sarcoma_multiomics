#!/usr/bin/env Rscript
##############################################################################
## Description: Find genomic drivers for FST and Metastasis
## Author:  Zhong Guorui
## Created: 2025-06-15
##############################################################################
## Load required functions and configurations
source(here::here("conf/study_conf.R"))

# library(RSQLite)
# library(biomaRt)
# library(clusterProfiler)
# library(msigdbr)
# library(ReactomePA)
# library(KEGGREST)

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

cnv_facet_gene_tbl |> 
    select(cytoband, cn_status) |> 
    distinct()

cosmic_gene_list <- LoadCOSMICGeneList()
oncokb_gene_list <- LoadOncoKBGeneList()
clinical_info <- LoadClinicalInfo()

# plot_date <- Sys.Date() |> format("%Y%m%d")
plot_date <- "251009"
plot_dir <- paste0("figures/", plot_date, "_oncoplot")

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
        filter(symbol %in% cnv_genes) |>
        mutate(cytoband_symbol = paste(cytoband, symbol, sep = "_")),
    variant_type = "cytoband_symbol",
    variant_class = "cn_status"
)

# cnv_gene_mat <- cnv_gene_mat[
#     grepl(paste(cnv_genes, collapse = "|"), rownames(cnv_gene_mat)),
# ]

## Oncoplot
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
    dir = plot_dir,
    filename = paste0(plot_date, "_figure1", ".pdf")
)

## "==========================================================================="
## Figure2 -----
## "==========================================================================="
### Plot parameters -----
title_size <- 7
text_size <- 6
column_title <- NULL
row_names_side <- "left"
pct_side <- "right"

sample_ids <- clinical_info |> 
    filter(FST.Group %in% c("Pre-FST", "Post-FST")) |>
    filter(Number.of.samples.per.patient >= 2) |> 
    # filter(Sample.ID %in% colnames(cnv_cytoband_mat)) |>
    pull(Sample.ID)

### "========================================================================"
### CNV gene matrix ----
### "========================================================================"
show_cnv_genes <- c(
    "KITLG",
    "SPRY4",
    "JS-2",
    "TERT",
    "CLPTM1L",
    "PRDM1",
    "AIM1",
    "FOXO3",
    "HACE1",
    
    "CDKN2A",
    "CDKN2B",
    "MTAP",
    "NOTCH2",
    "PAX7"
)

cnv_gene_mat <- GetVariantMat(
    data = cnv_facet_gene_tbl |> 
        filter(symbol %in% show_cnv_genes) |>
        mutate(cytoband_symbol = paste(cytoband, symbol, sep = "_")),
    variant_type = "cytoband_symbol",
    variant_class = "cn_status"
)

cnv_gene_mat <- cnv_gene_mat[, sample_ids, drop = FALSE]

### "========================================================================"
### CNV cytoband matrix ----
### "========================================================================"
## For same cytoband, we showly one alteration type with higher frequency
gistic_cn_cytoband_mat <- gistic_data$cn_cytoband_mat
gistic_cn_cytoband_mat[gistic_cn_cytoband_mat == "Amp"] <- "AMP"
gistic_cn_cytoband_mat[gistic_cn_cytoband_mat == "Del"] <- "DEL"

gc <- "Pre-FST_vs_Post-FST"

share_cytobands <- stat_res_list$share$gistic_cnv_cytoband[[gc]] |> 
    filter(group1_freq > 0.4 & group2_freq > 0.4) |> 
    mutate(
        cytoband_short = str_remove(cytoband, "^.*:")
    ) |> 
    mutate(freq = group1_freq + group2_freq) |>
    group_by(cytoband_short) |>
    slice_max(order_by = freq, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(group2_freq)) |> 
    pull(cytoband)

diff_cytobands <- stat_res_sig$gistic_cnv_cytoband[[gc]] |> 
    arrange(desc(group2_freq)) |> 
    pull(cytoband)

show_cytobands <- c(share_cytobands, diff_cytobands) |> unique()

cnv_cytoband_mat <- gistic_cn_cytoband_mat[
    show_cytobands, sample_ids, drop = FALSE
]

### "========================================================================"
### Top annotation ----
### "========================================================================"
top_annotation_list <- list(
    "TMB" = tmb_data[sample_ids, , drop = FALSE],
    "FGA" = fga_data[sample_ids, , drop = FALSE]
)

anno_list <- list()
            
for (anno_name in names(top_annotation_list)) {
    
    anno_data <- top_annotation_list[[anno_name]]

    numeric_cols <- sapply(anno_data, is.numeric)
    
    if (any(numeric_cols)) {
        
        plot_data <- anno_data[[which(numeric_cols)[1]]]

    } else {
        warning(paste("No numeric columns found in", anno_name))
        next
    }
    

    ## tick breaks
    data_range <- range(plot_data, na.rm = TRUE)
    data_max <- data_range[2]
    data_min <- data_range[1]

    # Create intelligent tick breaks
    if (data_max <= 1) {
        # For very small values, use finer steps
        tick_breaks <- pretty(c(data_min, data_max), n = 4)

    } else if (data_max <= 5) {
        # For small values, use integers
        tick_breaks <- pretty(c(data_min, data_max), n = 3)
    } else if (data_max <= 50) {
        # For medium values, use steps of 5 or 10
        tick_breaks <- pretty(c(data_min, data_max), n = 4)
    } else {
        # For large values, use larger steps
        tick_breaks <- pretty(c(data_min, data_max), n = 3)
    }
    
    # Ensure we include 0 and max
    tick_breaks <- unique(c(0, tick_breaks[tick_breaks > 0]))
    tick_breaks <- tick_breaks[tick_breaks <= data_max]
    
    # Format labels appropriately
    tick_labels <- ifelse(
        tick_breaks == floor(tick_breaks),
        as.character(as.integer(tick_breaks)),
        sprintf("%.1f", tick_breaks)
    )
    
    anno_list[[anno_name]] <- anno_barplot(
        plot_data,
        height = unit(2, "cm"),
        axis_param = list(
            at = tick_breaks,
            labels = tick_labels,
            gp = gpar(fontsize = text_size)
        ),
        gp = gpar(fill = "gray50")
    )

}

## Combine multiple annotations
if (length(anno_list) > 0) {
    # Use do.call instead of !!!
    top_annotation <- do.call(
        HeatmapAnnotation, 
        c(
            anno_list,
            list(
                annotation_name_gp = gpar(fontsize = text_size),
                annotation_name_side = "left",
                gap = unit(2, "mm")
            )
        )

    )
}

### "========================================================================"
### Bottom annotation ----
### "========================================================================"
sample_annotations <- c(
    "FST.Group", "Metastasis", "Specimen.Nature"
    # "Histology.subtype", "Molecular.subtype"
)

sample_annotation_df <- clinical_info |> 
    select(Sample.ID, all_of(sample_annotations)) |>
    filter(Sample.ID %in% sample_ids) |>
    column_to_rownames(var = "Sample.ID")

sample_annotation_colors <- list()

for (annotation in sample_annotations) {
    sample_annotation_colors[[annotation]] <- study_colors[[annotation]]
}

bottom_annotation <- HeatmapAnnotation(
    df = sample_annotation_df,
    col = sample_annotation_colors,
    annotation_height = unit(c(4, 4), "mm"),
    annotation_name_gp = gpar(fontsize = text_size),
    annotation_legend_param = list(
        labels_gp = gpar(fontsize = text_size), # Legend text
        title_gp = gpar(fontsize = title_size) # Legend title
    ),
    show_annotation_name = TRUE
)

### "========================================================================"
### Alteration colors ----
### "========================================================================"
alteration_colors <- study_colors$Alteration

alter_fun <- list(
    background = alter_graphic("rect", fill = "#CCCCCC"),
    AMP = alter_graphic(
        "rect", width = 0.9, height = 0.9,
        fill = study_colors$Alteration[["AMP"]]
    ),
    GAIN = alter_graphic(
        "rect", width = 0.9, height = 0.9,
        fill = study_colors$Alteration[["GAIN"]]
    ),
    HOMDEL = alter_graphic(
        "rect", width = 0.9, height = 0.9,
        fill = study_colors$Alteration[["HOMDEL"]]
    ),
    DEL = alter_graphic(
        "rect", width = 0.9, height = 0.9,
        fill = study_colors$Alteration[["DEL"]]
    )
)

alteration_legend <- Legend(
    labels = names(alter_fun)[-1],

    ## Custom graphics function to create rectangles with backgrounds and borders that match the heatmap
    graphics = list(
        # AMP - full size with background and border
        function(x, y, w, h) {
            # Background (grey)
            grid.rect(
                x, y, w, h,
                gp = gpar(
                    fill = "#CCCCCC", col = "white", lwd = 0.5
                )
            )
            # Foreground rectangle
            grid.rect(
                x, y, w * 0.9, h * 0.9,
                gp = gpar(
                    fill = study_colors$Alteration[["AMP"]],
                    col = "white", lwd = 0.5
                )
            )
        },
        # GAIN - full size with background and border
        function(x, y, w, h) {
            grid.rect(
                x, y, w, h,
                gp = gpar(
                    fill = "#CCCCCC", col = "white", lwd = 0.5
                )
            )
            grid.rect(
                x, y, w * 0.9, h * 0.9,
                gp = gpar(
                    fill = study_colors$Alteration[["GAIN"]],
                    col = "white", lwd = 0.5
                )
            )
        },
        # HOMDEL - full size with background and border
        function(x, y, w, h) {
            grid.rect(
                x, y, w, h,
                gp = gpar(fill = "#CCCCCC", col = "white", lwd = 0.5)
            )
            grid.rect(
                x, y, w * 0.9, h * 0.9,
                gp = gpar(
                    fill = study_colors$Alteration[["HOMDEL"]],
                    col = "white", lwd = 0.5
                )
            )
        },
        # DEL - full size with background and border
        function(x, y, w, h) {
            grid.rect(
                x, y, w, h,
                gp = gpar(fill = "#CCCCCC", col = "white", lwd = 0.5)
            )
            grid.rect(
                x, y, w * 0.9, h * 0.9,
                gp = gpar(
                    fill = study_colors$Alteration[["DEL"]],
                    col = "white", lwd = 0.5
                )
            )
        }
    ),
    title = "Alteration",
    title_gp = gpar(fontsize = title_size), # fontface = "bold"
    labels_gp = gpar(fontsize = text_size),
    grid_height = unit(4, "mm"),
    grid_width = unit(4, "mm"),
    direction = "vertical"
)
### "========================================================================"
### CNV gene oncoprint ----
### "========================================================================"
# cnv_gene_row_order <- rownames(cnv_gene_mat)

ht_cnv_gene <- oncoPrint(
    cnv_gene_mat,
    alter_fun = alter_fun,
    col = alteration_colors,
    column_title = column_title,
    bottom_annotation = NULL,
    top_annotation = top_annotation,
    # left_annotation = left_annotation,
    show_row_names = TRUE,
    row_names_side = row_names_side,
    # row_order  = cnv_gene_row_order,
    row_names_gp = gpar(fontsize = text_size),
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = text_size),
    column_title_gp = gpar(fontsize = title_size),
    show_heatmap_legend = FALSE,
    show_pct = TRUE,
    pct_gp = gpar(fontsize = text_size),
    pct_digits = 0,
    pct_side = pct_side,
)

### "========================================================================"
### CNV cyoband oncoprint ----
### "========================================================================"
row_order <- rownames(cnv_cytoband_mat)
ht_cnv_cytoband <- oncoPrint(
    cnv_cytoband_mat,
    alter_fun = alter_fun,
    col = alteration_colors,
    top_annotation = NULL,
    bottom_annotation = bottom_annotation,
    show_row_names = TRUE,
    row_names_side = "left",
    row_order = row_order,
    row_names_gp = gpar(fontsize = text_size),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = text_size),
    column_title_gp = gpar(fontsize = title_size),
    column_order = sample_ids,
    show_heatmap_legend = FALSE,
    show_pct = TRUE,
    pct_gp = gpar(fontsize = text_size),
    pct_digits = 0,
    pct_side = "right"
)

### "========================================================================"
### Combine and save oncoprint ----
### "========================================================================"
ht_list <- ht_cnv_gene %v% ht_cnv_cytoband

filename <- paste0(plot_date, "_figure2")

width <- 12
height <- 8

img_type <- c(".pdf", ".png")
for (img in img_type) {

    file <- here(plot_dir, paste0(filename, img))

    if (img == ".pdf") {

        # CairoPDF(file, width = width, height = height)
        pdf(file, width = width, height = height)

    } else if (img == ".png") {

        png(
            file, width = width, height = height, res = 600, units = "in",
            fonts = "Arial"
        )
    }

    draw(
        ht_list, 
        main_heatmap = 2,
        # ht_gap = unit(0.2, "cm"),

        ## Heatmap legend
        # show_heatmap_legend = TRUE,
        # heatmap_legend_side = "bottom",

        ## Annotation legend
        # annotation_legend_side = "bottom",
    
        merge_legends = TRUE,
        legend_gap = unit(1, "cm"),
        padding = unit(c(1, 1, 1, 1), "cm"),

        ## Add manual heatmap legend
        heatmap_legend_list = list(alteration_legend)
    )
    
    dev.off()

    message(
        paste0("Saved oncoplot: ", file)
    )
}

## "==========================================================================="
## Figure3 -----
## "==========================================================================="

## "==========================================================================="
## Oncoplot SNV/indels -----
## "==========================================================================="