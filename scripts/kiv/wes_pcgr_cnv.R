#!/usr/bin/env Rscript
##############################################################################
## Description: This script processes PCGR data and performs cohort analysis
## Author:  Zhong Guorui
## Created: 2025-06-15
##############################################################################
## Load required libraries and functions
suppressPackageStartupMessages(
    suppressWarnings({
        library(maftools)
        library(circlize)
        library(BSgenome.Hsapiens.UCSC.hg38)
        library(GenomicRanges)
        library(ggsci)
        library(ComplexHeatmap)
        library(reshape2)
    })
)
source(here::here("scripts/lib/study_lib.R"))
source(here::here("conf/study_conf.R"))

## "==========================================================================="
## Preprocess PCGR CNV data --------------
## "==========================================================================="
## Merge all PCGR outputs
pcgr_dir <- "data/WES/PCGR"
pcgr_cna_raw <- CollectPCGRCNVData(dir = pcgr_dir)

unique(pcgr_cna_raw$event_type)
unique(pcgr_cna_raw$variant_class)
length(unique(pcgr_cna_raw$sample_id))
pcgr_cna_raw |> 
    filter(variant_class != "undefined") |> 
    pull(sample_id) |> 
    unique() |> 
    length()

## Save the raw PCGR CNV data
SaveData(
    pcgr_cna_raw, 
    dir = "data/processed", 
    filename = "wes_pcgr_DFSP_cohort_cnv_raw_data"
)

## Load the PCGR annotated CNV data
pcgr_cna_raw <- LoadData(
    dir = "data/processed", 
    filename = "wes_pcgr_DFSP_cohort_cnv_raw_data"
)

## Load the CNV facets data
cnv_facet_raw <- LoadData(
    dir = "data/processed",
    filename = "wes_cnv_facets_DFSP_cohort_raw"
)

## 5 samples do not have facet estimated purity
## DFSP-042-T, DFSP-060-T, DFSP-331-T, DFSP-341-T, DFSP-354-T
cnv_facet_raw |> 
    filter(is.na(purity) | is.na(ploidy)) |> 
    select(sample, purity, ploidy) |> 
    distinct()

## Define the CNV events
## 10.1016/j.jtho.2021.02.023
pcgr_cna_tbl <- pcgr_cna_raw |> 
    separate_wider_delim(
            var_id, 
            delim = ":", 
            names = c(
                "CHROM", "position_range", "cn_major_parsed", "cn_minor_parsed"
            ),
            cols_remove = FALSE
    ) |> 
    select(-cn_major_parsed, -cn_minor_parsed) |>
    separate_wider_delim(
        position_range, 
        names = c("POS", "END"), 
        delim = "-", 
        cols_remove = FALSE
    ) |> 
    left_join(
        cnv_facet_raw |> 
            select(sample, CHROM, POS, END, SVTYPE, ploidy, purity),
        by = c("sample_id" = "sample", "CHROM", "POS", "END")
    ) |> 
    select(-CHROM, -POS, -END) |>
    mutate(
        ploidy = as.numeric(ploidy),
        purity = as.numeric(purity)
    ) |> 
    dplyr::rename(
        facet_svtype = SVTYPE,
        facet_ploidy = ploidy,
        facet_purity = purity
    ) |> 
    relocate(
        matches("facet"), .after = variant_class
        # facet_svtype, facet_ploidy, facet_purity, 
    ) |> 
    ## Define the CNV events
    mutate(
        cn_status = case_when(
            ## "Deep deletion"
            !is.na(facet_ploidy) & total_cn == 0 ~ "HOMDEL",
            
            ## Deletion
            !is.na(facet_ploidy) & (total_cn -facet_ploidy) <= -1 ~ "DEL",
            
            ## Neutral
            !is.na(facet_ploidy) & 
                ((total_cn -facet_ploidy) > -1) &
                ((total_cn -facet_ploidy) < 1) ~ "NEUTRAL",
            
            ## Amplification
            !is.na(facet_ploidy) & 
                ((total_cn - facet_ploidy) >= 1) &
                ((total_cn - facet_ploidy) < 3) ~ "GAIN",

            ## Deep amplification
            !is.na(facet_ploidy) & (total_cn - facet_ploidy) >=3 ~ "AMP",
        ),
        .after = variant_class
    ) |> 
    ## Keep only changed CNV events
    dplyr::filter(cn_status != "NEUTRAL") |> 
    separate_wider_delim(
        position_range,
        names = c("start", "end"),
        delim = "-"
    ) |> 
    mutate(
        chrom = str_extract(var_id, "chr[0-9XY]+"), .after = "sample_id",
        start = as.integer(start),
        end = as.integer(end)
    ) |> 
    mutate(
        cytoband_alteration = paste0(cytoband, "_", cn_status),
        .after = cn_status
    ) |> 
    mutate(
        gene_alteration = paste0(symbol, "_", cn_status),
        .after = symbol
    ) |> 
    arrange(desc(facet_ploidy))

## One sample do not have non-neutral CNA: "DFSP-169-T"
length(unique(pcgr_cna_tbl$sample_id))
setdiff(pcgr_cna_raw$sample_id, pcgr_cna_tbl$sample_id)

## Save the filtered Non-Neutral PCGR CNV data
SaveData(
    pcgr_cna_tbl,
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_cnv_filtered"
)

## "========================================================================="
## Different and shared cytobands---- 
## "========================================================================="
pcgr_cna_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_cnv_filtered"
)

## Analyze the somatic matched samples
clinical_info <- LoadClinicalInfo() |>
    filter(Somatic.Status %in% c("Matched"))

cohort_total_samples <- clinical_info |> nrow()

cna_data <- pcgr_cna_tbl |> 
    filter(sample_id %in% clinical_info$Sample.ID) |> 
    left_join(clinical_info, by = c("sample_id" = "Sample.ID")) |> 
    relocate(FST.Group, .after = sample_id)

## Get the significantly different cytobands across groups
all_diff_cytobands <- map_dfr(
    group_comparisons,
    function(comp) {
        GetPCGRCNVGroupDiffCytoband(
            data = cna_data,
            var = "FST.Group",
            cn_column = "cn_status",
            event_type = NULL,
            somatic_matched = TRUE,
            group1 = comp$group1,
            group2 = comp$group2,
            p_adjust_method = "fdr",
            min_altered_samples = 3,
            group_comparison = comp$name
        )
    }
)

p_val_thres <- 0.05

## At least 15% samples are mutated
freq_thres <- 0.15

sig_diff_cytobands <- all_diff_cytobands |> 
    filter(fisher_p < p_val_thres) |>
    filter(group1_freq > freq_thres | group2_freq > freq_thres) |>
    mutate(cytoband_alteration = paste0(cytoband, "_", cn_status)) |> 
    group_by(cytoband_alteration, cytoband, cn_status) |> 
    summarise(
        n_comparisons_present = n(),
        comparisons_present = paste(group_comparison, collapse = ", "),
        mean_genes_affected = mean(n_genes_affected, na.rm = TRUE),
        mean_segment_size = mean(mean_segment_size, na.rm = TRUE),
        .groups = "drop"
    ) |> 
    arrange(desc(n_comparisons_present))

## Get the shared cytobands across groups
shared_cytobands <- GetPCGRCNVGroupShareCytoband(
    data = cna_data,
    var = "FST.Group",
    cn_column = "cn_status",
    somatic_matched = TRUE,
    min_freq_per_group = 0.3
)

## Save the cytoband stat data
SaveData(
    list(
        all_diff_cytobands = all_diff_cytobands,
        sig_diff_cytobands = sig_diff_cytobands,
        shared_cytobands = shared_cytobands
    ),
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_somatic_matched_cytoband_stats"
)

## Save the different and shared regions across groups
write_xlsx(
    list(
        all_diff_cytobands = all_diff_cytobands,
        sig_diff_cytobands = sig_diff_cytobands,
        shared_cytobands = shared_cytobands$shared_cytobands
    ),
    here(
        "outputs", "wes_pcgr_DFSP_cohort_somatic_matched_diff_shared_cytoband.xlsx"
    )
)

## "==========================================================================="
## Chromosomal plot cytobands --------------
## "==========================================================================="
## Load the cytoband stat res
cytoband_stat_res <- LoadPCGRCytobandStatData()

all_diff_cytobands <- cytoband_stat_res$all_diff_cytobands 

all_diff_cytobands <- all_diff_cytobands |> 
    filter(fisher_p < p_value_thres) |> 
    filter(group1_freq > freq_thres | group2_freq > freq_thres) |> 
    mutate(cytoband_alteration = paste0(cytoband, "_", cn_status))

peak_label_list <- list(
    
    "U-DFSP" = all_diff_cytobands |> 
        filter(
            group_comparison %in% c("U-DFSP_vs_Pre-FST")
        ) |> 
        pull(cytoband_alteration) |> 
        unique(),
        
    "Pre-FST" = all_diff_cytobands |> 
        filter(
            group_comparison %in% c("U-DFSP_vs_Pre-FST")
        ) |> 
        pull(cytoband_alteration) |> 
        unique(),

    "Post-FST" = all_diff_cytobands |> 
        filter(
            group_comparison %in% c(
                "U-DFSP_vs_Post-FST",
                "U-DFSP_vs_FS-DFSP",
                "Pre-FST_vs_Post-FST",
                "Pre-FST_vs_FS-DFSP",
                "Post-FST_vs_FS-DFSP",
                "U-DFSP+Pre-FST_vs_Post-FST+FS-DFSP"
            )
        ) |> 
        pull(cytoband_alteration) |> 
        unique(),
    
    "FS-DFSP" = all_diff_cytobands |> 
        filter(
            group_comparison %in% c(
                "U-DFSP_vs_Post-FST",
                "U-DFSP_vs_FS-DFSP",
                "Pre-FST_vs_Post-FST",
                "Pre-FST_vs_FS-DFSP",
                "Post-FST_vs_FS-DFSP",
                "U-DFSP+Pre-FST_vs_Post-FST+FS-DFSP"
            )
        ) |> 
        pull(cytoband_alteration) |> 
        unique()
)

# peak_label <- c("chr7:q36.1_DEL", "chr6:p21.33_HOMDEL")

gistic_dir <- "data/WES/GISTIC2/somatic_matched"

plots <- list()
for (group in FST.Group) {
    
    ## Load the gistic score
    gistic_score <- LoadGisticScore(
        gistic_dir = gistic_dir,
        group = group
    )

    ## the sample info in entire cohort
    n_samples <- LoadClinicalInfo() |> 
        filter(Somatic.Status == "Matched") |> 
        filter(FST.Group == group) |> 
        nrow()

    peak_label <- peak_label_list[[group]]
    
    ## Get the peak coordinates
    peaks_data <- GetPCGRGisticOverlap2(
        peak_label = peak_label,
        pcgr_cna_data = cna_data, 
        gistic_score_data = gistic_score,
        min_overlap = 0.1
    )

    peaks_data <- AddPeaksPCGRCountData(
        peaks_data,
        pcgr_cna_data = cna_data,
        var = "FST.Group",
        group = group,
        is_somatic_matched = TRUE
    )

    ## Plot the enriched regions
    plots[[group]] <- GisticChromPlot(
        data = gistic_score,
        peaks_data = peaks_data,
        y_lim = c(-2, 3),
        x_expand = c(0.02, 0),
        title = paste0(group, " (", n_samples, ")"),
        show_horizontal_grid = TRUE,
        show_vertical_grid = TRUE,
        grid_line_type = "dotted",
        grid_line_color = "grey80",
        grid_line_alpha = 0.5,
        grid_line_width = 0.2,
        show_chromosome_ideogram = FALSE,
        ideogram_label_size = 4,
        rel_heights = c(4, 0.3)
    )
}

SavePlot(
    plot = wrap_plots(plots, ncol = 2) +
        plot_layout(guides = "collect"),
    width = 12,
    height = 5,
    dir = "figures/gistic2",
    filename = "wes_facet_pcgr_DFSP_cohort_gistic_plots_FST.Group"
)

## "========================================================================"
## Oncoplot cytobands ----
## "========================================================================"
## Cytoband freq data
cytoband_freq <- cna_data |> 
    select(sample_id, cytoband) |>
    distinct() |>
    group_by(cytoband) |> 
    summarise(
        n_samples_altered = n_distinct(sample_id),
        .groups = "drop"
    ) |> 
    mutate(
        cohort_freq = n_samples_altered / cohort_total_samples
    ) |> 
    arrange(desc(cohort_freq))

## Prepare the oncoplot data
cytoband_oncoplot_data <- cna_data  |> 
    select(sample_id, cytoband, cn_status) |>
    distinct() |>
    group_by(sample_id, cytoband) |>
    summarise(
        alterations = paste(unique(cn_status), collapse = ","),
        .groups = "drop"
    ) |> 
    arrange(desc(alterations)) |> 
    mutate(
        final_alteration = case_when(
            ## Count number of different alteration types (split by comma)
            str_count(alterations, ",") >= 1 ~ "MULTI",  # Multiple different alterations
            ## Single alterations
            alterations == "HOMDEL" ~ "HOMDEL",
            alterations == "DEL" ~ "DEL", 
            alterations == "GAIN" ~ "GAIN",
            alterations == "AMP" ~ "AMP",
            TRUE ~ alterations
        )
    ) |> 
    distinct()

## Select the cytobands to plot
selected_cytobands <- all_diff_cytobands |> 
    pull(cytoband) |> 
    unique()

## Prepare the oncoplot data
onco_data <- expand_grid(
    sample_id = clinical_info$Sample.ID,
    cytoband = selected_cytobands
) |> 
    left_join(
        cytoband_oncoplot_data |> select(-alterations), 
        by = c("sample_id", "cytoband")
    )

onco_mat <- onco_data |> 
    mutate(
        final_alteration = if_else(
            is.na(final_alteration), "", final_alteration
        )
    ) |> 
    ## Handle duplicates by keeping distinct combinations
    distinct(sample_id, cytoband, .keep_all = TRUE) |>
    pivot_wider(
        names_from = sample_id,
        values_from = final_alteration,
        values_fill = ""
    ) |>
    column_to_rownames("cytoband") |>
    as.matrix()

dim(onco_mat)

## Plot the oncoplot
plot_para <- list(
    entire_cohort_FST_group  = list(
        sample_sorted_by = "FST.Group",
        sample_sorted_level = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
    ),
    paired_samples = list(
        sample_sorted_by = NULL,
        sample_sorted_level = NULL
    ),
    EPIC_meth_group = list(
        sample_sorted_by = "Meth.Subtype.Main.2Class",
        sample_sorted_level = c("Meth1", "Meth2")
    ),
    HRD_group = list(
        sample_sorted_by = "HRD.cat",
        sample_sorted_level = c("High", "Low")
    ),
    MSI_group = list(
        sample_sorted_by = "MSI.cat",
        sample_sorted_level = c("High", "Low")
    )
)

sample_info <- LoadClinicalInfo() |> 
    filter(Somatic.Status == "Matched")

for (i in names(plot_para)) {

    filename <- paste0(
        "wes_cnvfacets_pcgr_oncoplot_somatic_matched_sort_by_", 
        plot_para[[i]]$sample_sorted_by
    )

    ## all samples
    all_sample_ids <- colnames(onco_mat)

    if (i == "paired_samples") {

        sample_ids <- sample_info |> 
            filter(grepl("Paired", Histology.Nature)) |> 
            pull(Sample.ID)

        sample_ids <- all_sample_ids[all_sample_ids %in% sample_ids]
        
        n_samples <- length(sample_ids)

    } else if (i == "EPIC_meth_group") {

        sample_ids <- sample_info |> 
            filter(!is.na(Meth.Subtype.Main.2Class)) |> 
            pull(Sample.ID)

        sample_ids <- all_sample_ids[all_sample_ids %in% sample_ids]

        n_samples <- length(sample_ids)
    }

    column_title <- paste0(
            "Somatic matched samples, N = ", n_samples,
            "\n", "(Sorted by ", plot_para[[i]]$sample_sorted_by, ")"
    )

    ## Subset the matrix
    mat <- onco_mat[, sample_ids, drop = FALSE]

    sample_annotation <- c(
        "FST.Group", "Metastasis", "Specimen.Nature",
        "Meth.Subtype.Main.2Class", "HRD.cat", "MSI.cat"
    )
    
    ## Plot the oncoplot
    PCGRComplexOncoplot(
        mat = mat,
        sample_annotation = sample_annotation,
        sample_sorted_by = plot_para[[i]]$sample_sorted_by,
        sample_sorted_level = plot_para[[i]]$sample_sorted_level,
        column_title = column_title,
        width = 20,
        height = 10,
        dir = "figures/oncoplot",
        filename = filename
    )
}


## "========================================================================="
## Different and shared genes -------
## "========================================================================="
pcgr_cna_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_cnv_filtered"
)

## Analyze the somatic matched samples
clinical_info <- LoadClinicalInfo() |>
    filter(Somatic.Status %in% c("Matched"))

cohort_total_samples <- clinical_info |> nrow()

cna_data <- pcgr_cna_tbl |> 
    filter(sample_id %in% clinical_info$Sample.ID) |> 
    left_join(clinical_info, by = c("sample_id" = "Sample.ID")) |> 
    relocate(FST.Group, .after = sample_id) |> 
    mutate(
        gene_alteration = paste0(symbol, "_", cn_status),
        .after = symbol
    )

## Get the significantly different genes across groups
all_diff_genes <- map_dfr(
    group_comparisons,
    function(comp) {
        GetPCGRCNVGroupDiffGene(
            data = cna_data,
            var = "FST.Group",
            cn_column = "cn_status",
            event_type = NULL,
            somatic_matched = TRUE,
            group1 = comp$group1,
            group2 = comp$group2,
            p_adjust_method = "fdr",
            min_altered_samples = 3,
            group_comparison = comp$name
        )
    }
)

sig_diff_genes <- all_diff_genes |> 
    filter(fisher_p < p_val_thres) |>
    filter(group1_freq > freq_thres | group2_freq > freq_thres) |>
    mutate(gene_alteration = paste0(symbol, "_", cn_status)) |>
    group_by(gene_alteration, symbol, cn_status) |>
    summarise(
        n_comparisons_present = n(),
        comparisons_present = paste(group_comparison, collapse = ", "),
        .groups = "drop"
    ) |> 
    arrange(desc(n_comparisons_present))

## Get the shared genes across groups
shared_genes <- GetPCGRCNVGroupShareGene(
    data = cna_data,
    var = "FST.Group",
    cn_column = "cn_status",
    somatic_matched = TRUE,
    min_freq_per_group = 0.3
)

## Save the gene stat data
SaveData(
    list(
        all_diff_genes = all_diff_genes,
        sig_diff_genes = sig_diff_genes,
        shared_genes = shared_genes
    ),
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_somatic_matched_gene_stats"
)

## Save the different and shared regions across groups
write_xlsx(
    list(
        all_diff_genes = all_diff_genes,
        sig_diff_genes = sig_diff_genes,
        shared_genes = shared_genes$shared_genes
    ),
    here(
        "outputs", "wes_pcgr_DFSP_cohort_somatic_matched_diff_shared_gene.xlsx"
    )
)

## "========================================================================"
## Oncoplot CNV genes puclic cohort ----
## "========================================================================"
## Calculate the frequency for all CNV genes
gene_freq <- cna_data |>
    select(sample_id, symbol) |>
    distinct() |> 
    group_by(symbol) |> 
    summarise(
        n_samples_altered = n_distinct(sample_id),
        .groups = "drop"
    ) |> 
    mutate(cohort_freq = n_samples_altered / cohort_total_samples) |>
    arrange(desc(cohort_freq))

## Prepare the oncoplot data
gene_oncoplot_data <- cna_data  |> 
    select(sample_id, symbol, cn_status) |>
    distinct() |>
    group_by(sample_id, symbol) |>
    summarise(
        alterations = paste(unique(cn_status), collapse = ","),
        .groups = "drop"
    ) |> 
    arrange(desc(alterations)) |> 
    mutate(
        final_alteration = case_when(
            ## Count number of different alteration types (split by comma)
            str_count(alterations, ",") >= 1 ~ "MULTI",  # Multiple different alterations
            ## Single alterations
            alterations == "HOMDEL" ~ "HOMDEL",
            alterations == "DEL" ~ "DEL", 
            alterations == "GAIN" ~ "GAIN",
            alterations == "AMP" ~ "AMP",
            TRUE ~ alterations
        )
    ) |> 
    distinct()

## cBioportal sarcoma CNV genes
sarcoma_cna_file <- here("data/public/cBioPortal_sarcoma_CNA_Genes.txt")
col_names <- c(
    "symbol", 
    "gistic_q_value",
    "cytoband",
    "cn_status",
    "profiled_samples",
    "#",
    "freq",
    "is_oncoKB"
)
sarcoma_cna_data <- read_tsv(sarcoma_cna_file, show_col_types = FALSE)
colnames(sarcoma_cna_data) <- col_names

sarcoma_cna_data <- sarcoma_cna_data |> 
    mutate(freq = parse_number(freq)) |> 
    arrange(desc(freq)) |>
    filter(freq >= 10) |>
    distinct(symbol, .keep_all = TRUE)

## Define public cohorts
dfsp_public_cohorts <- list(
    cbioportal_sarcoma = list(
        title = "Top mutated genes (cBioPortal Sarcoma, Freq > 1%)",
        cnv = sarcoma_cna_data$symbol
    ),
    ## Modern pathology
    smith_2025 = list(
        title = "Top mutated genes (Smith et al. 2025 Reported, n=53)",
        cnv = c(
            "PRKAR1A", "H3F3B", "MSI2", 
            "SRSF2", "BRIP1", "DDX5", "CLTC", 
            "GNA13", "CRKL", "MCL1", "MYH9",
            "SEPT9", "SOX10",
            "ATP1A1", "CANT1", "CLTCL1", "DAXX", "ERG", "EZR",
            "FGFR1OP",
            "LIFR",
            "MAF",
            "NOTCH2",
            "PDGFRB",
            "PRRX1",
            "RAC1",
            "RICTOR",
            "RNF43"
        )
    ),
    ## BJD
    peng_2022 = list(
        title = "Top mutated genes (Peng et al. 2022 Reported, n=59)",
        cnv = c(
            "TERT", "CDKN2A", "CDKN2B", 
            "AKT1", "NFKBIA", "BTBD7", "SPHK1", 
            "ITGB4", "COL1A1", "PDGFB", "PDGFD"
        )
    )
)

## Genes to plot
selected_genes <- c(
    dfsp_public_cohorts$smith_2025$cnv,
    dfsp_public_cohorts$peng_2022$cnv
) |> 
    unique()

## Construct the oncoplot data
onco_data <- expand_grid(
    sample_id = clinical_info$Sample.ID,
    symbol = selected_genes
) |> 
    left_join(
        gene_oncoplot_data |> select(-alterations), 
        by = c("sample_id", "symbol")
    ) |> 
    arrange(final_alteration)

onco_mat <- onco_data |> 
    mutate(
        final_alteration = if_else(
            is.na(final_alteration), "", final_alteration
        )
    ) |> 
    ## Handle duplicates by keeping distinct combinations
    distinct(sample_id, symbol, .keep_all = TRUE) |>
    pivot_wider(
        names_from = sample_id,
        values_from = final_alteration,
        values_fill = ""
    ) |>
    column_to_rownames("symbol") |>
    as.matrix()

## Plot the oncoplot
sample_sorted_by <- "FST.Group"  # Define which variable to sort by
sample_sorted_level <- c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
column_title <- paste0(
    "Somatic matched samples, N = ", cohort_total_samples, 
    "\n", "Public cohort (Peng et al. 2022; Smith et al. 2025)"
)

filename <- paste0(
    "wes_cnvfacets_pcgr_oncoplot_public_cohort_CNV_genes_somatic_matched_sort_by_", 
    sample_sorted_by
)

PCGRComplexOncoplot(
    mat = onco_mat,
    sample_annotation = c("FST.Group", "Metastasis", "Specimen.Nature"),
    sample_sorted_by = sample_sorted_by,
    sample_sorted_level = sample_sorted_level,
    column_title = column_title,
    width = 20,
    height = 10,
    dir = "figures/oncoplot",
    filename = filename
)
