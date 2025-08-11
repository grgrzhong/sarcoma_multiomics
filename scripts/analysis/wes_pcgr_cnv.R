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
suppressPackageStartupMessages(
    suppressWarnings({
        library(maftools)
        library(circlize)
        library(BSgenome.Hsapiens.UCSC.hg38)
        library(ggsci)
        library(ComplexHeatmap)
        library(reshape2)
    })
)
source(here::here("scripts/lib/study_lib.R"))
source(here::here("conf/study_conf.R"))

## Output directories
data_dir <- "data/processed"
plot_dir <- "figures/wes/pcgr"

## Capture size
capture_size <- 34 

## Clinical information
clinical_info <- LoadClinicalInfo()

## "==========================================================================="
## Gistic2 scrore plots --------------
## "==========================================================================="
gistic_dir <- "data/wes/GISTIC2/somatic_matched/U-DFSP"

        
## "==========================================================================="
## Collect PCGR CNV data --------------
## "==========================================================================="
## Merge all PCGR outputs
pcgr_dir <- "data/wes/PCGR"
pcgr_cna_raw <- CollectPCGRCNAData(dir = pcgr_dir)

unique(pcgr_cna_raw$event_type)
unique(pcgr_cna_raw$variant_class)

## Save the raw PCGR CNV data
SaveData(
    pcgr_cna_raw, 
    dir = "data/processed", 
    filename = "wes_pcgr_DFSP_cohort_cnv_raw_data"
)

## "==========================================================================="
## Preprocess PCGR CNV data --------------
## "==========================================================================="
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

## 5 samples do not have estimated purity
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
    arrange(desc(facet_ploidy))

## Save the filtered Non-Neutral PCGR CNV data
SaveData(
    pcgr_cna_tbl,
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_cnv_filtered"
)

## "========================================================================="
## Differentially mutated and shared regions across groups
## "========================================================================="
pcgr_cna_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_cnv_filtered"
)

length(unique(pcgr_cna_tbl$sample_id))

clinical_info <- LoadClinicalInfo() |> 
    filter(Somatic.Status %in% c("Matched"))

total_samples <- clinical_info |> nrow()

data <- pcgr_cna_tbl  |> 
    left_join(clinical_info, by = c("sample_id" = "Sample.ID")) |> 
    filter(sample_id %in% clinical_info$Sample.ID) |> 
    relocate(FST.Group, .after = sample_id)

## Get the significantly different cytobands across groups
all_diff_cytoands <- map_dfr(
    group_comparisons,
    function(comp) {
        GetPCGRCNVGroupDiffCytoband(
            data = data,
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

sig_diff_cytobands <- all_diff_cytoands |> 
    filter(fisher_p < 0.05) |>
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
    data,
    var = "FST.Group",
    cn_column = "cn_status",
    event_type = NULL,
    somatic_matched = TRUE,
    min_freq_per_group = 0.3,
    min_groups_present = 1
)

show_cytobands <- shared_cytobands$shared_cytobands |> 
    select(cytoband_alteration, cytoband, cn_status) |> 
    bind_rows(
        sig_diff_cytobands |> 
        select(cytoband_alteration, cytoband, cn_status)
    ) |> 
    distinct()

## Prepare the oncoplot data
cytoband_oncoplot_data <- data  |> 
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

onco_data <- expand_grid(
    sample_id = clinical_info$Sample.ID,
    cytoband = show_cytobands$cytoband
) |> 
    left_join(
        cytoband_oncoplot_data |> select(-alterations), 
        by = c("sample_id", "cytoband")
    )

## Calculate the frequency of show cytobands
onco_freq <- cytoband_oncoplot_data |> 
    select(-alterations) |> 
    distinct() |> 
    group_by(cytoband) |> 
    summarise(
        altered_samples = n(),
        freq = altered_samples / total_samples
    ) |> 
    arrange(desc(altered_samples))

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

## Plot the oncoplot
sample_sorted_level <- c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
sample_sorted_by <- "FST.Group"
column_title <- paste0(
    "Somatic matched samples, N = ", total_samples
)
filename <- paste0(
    "wes_cnvfacets_pcgr_oncoplot_somatic_matched_sort_by_", 
    sample_sorted_by
)

GenerateCytobandOncoplot(
    mat = mat,
    sample_annotation = c("FST.Group", "Metastasis"),
    sample_sorted_by = sample_sorted_by,
    sample_sorted_level = sample_sorted_level,
    column_title = column_title,
    width = 18,
    height = 12,
    dir = "figures/oncoplot",
    filename = filename
)


## "========================================================================="
## Differentially mutated and shared genes across groups
## "========================================================================="
data <- pcgr_cna_tbl  |> 
    left_join(
        clinical_info |> select(Sample.ID, FST.Group),
        by = c("sample_id" = "Sample.ID")
    ) |> 
    relocate(FST.Group, .after = sample_id)

gene_stat <- GetPCGRCNVGroupStatRes(
    data = data,
    var = "FST.Group",
    cn_column = "cn_status",
    event_type = NULL,
    somatic_matched = TRUE,
    group1 = "U-DFSP",
    group2 = "FS-DFSP",
    p_adjust_method = "fdr",
    group_comparison = "U-DFSP_vs_FS-DFSP"
)

gene_stat |> 
    filter(significantly_different)

length(unique(gene_stat$symbol))

## Perform all pairwise comparisons

## Get all differentially mutated genes across groups
all_diff_genes <- map_dfr(
    comparisons,
    function(comp) {
        GetPCGRCNVGroupStatRes(
            data = data,
            var = "FST.Group",
            cn_column = "cn_status",
            event_type = NULL,
            somatic_matched = TRUE,
            group1 = comp$group1,
            group2 = comp$group2,
            p_adjust_method = "fdr",
            group_comparison = comp$name
        )
    }
)

## Filter for significant results
sig_diff_genes <- all_diff_genes |>
    filter(
        fisher_q < 0.5,           # FDR < 10%
        abs(freq_diff) >= 0.2,    # At least 20% frequency difference
        total_altered >= 3         # At least 3 total alterations
    ) |>
    arrange(fisher_q, desc(abs(freq_diff)))

share_genes <- GetPCGRCNVGroupShareGene(
    data,
    var = "FST.Group",
    cn_column = "cn_status",
    event_type = NULL,
    somatic_matched = TRUE,
    min_freq_per_group = 0.3,
    min_groups_present = 2
)

## Get a list of genes from the shared and differentiall mutated
cnv_oncoplot_genes <- unique(
    c(sig_diff_genes$symbol, share_genes$symbol)
)

## Remove the genes without gene name annotation
cnv_oncoplot_genes <- cnv_oncoplot_genes[!grepl("ENSG", cnv_oncoplot_genes)]

## Public sarcoma genes
file <- here("data/public/cBioPortal_sarcoma_CNA_Genes.txt")
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
public_sarcoma_cna_data <- read_tsv(file, show_col_types = FALSE)
colnames(public_sarcoma_cna_data) <- col_names
unique(public_sarcoma_cna_data$cn_status)
public_sarcoma_cna_genes  <- public_sarcoma_cna_data |> 
    mutate(freq = str_replace(freq, "%", "")) |>
    mutate(freq = as.numeric(freq)) |> 
    arrange(desc(freq)) |> 
    # filter(freq >= 10) |> 
    pull(symbol)

cnv_oncoplot_genes[cnv_oncoplot_genes %in% public_sarcoma_cna_genes]

total_cohort_size <- clinical_info |> 
    filter(Somatic.Status == "Matched") |> 
    nrow()

## genes freq data in entire cohort
data |>
    filter(symbol %in% cnv_oncoplot_genes)


oncoplot_data <- data |>
    filter(symbol %in% cnv_oncoplot_genes) |>
    select(sample_id, symbol, cn_status, FST.Group) |>
    ## Simplify cn_status for plotting
    mutate(
        alteration = case_when(
            cn_status == "Deep deletion" ~ "HOMDEL",
            cn_status == "Deletion" ~ "DEL",
            cn_status %in% c("Amplification", "Deep amplification") ~ "AMP",
            TRUE ~ cn_status
        )
    ) |>
    ## Handle multiple alterations per gene per sample
    group_by(sample_id, symbol, FST.Group) |>
    summarise(
        alterations = paste(unique(alteration), collapse = ";"),
        .groups = "drop"
    ) |>
    mutate(
        final_alteration = case_when(
            str_detect(alterations, "AMP") & str_detect(alterations, "DEL") ~ "Multi",
            str_detect(alterations, "AMP") ~ "AMP",
            str_detect(alterations, "DEL") ~ "DEL",
            TRUE ~ alterations
        )
    ) |>
    select(sample_id, symbol, final_alteration)

all_samples <- clinical_info |>
    filter(Somatic.Status == "Matched") |>
    select(Sample.ID, FST.Group) |>
    rename(sample_id = Sample.ID)

complete_combinations <- expand_grid(
    sample_id = all_samples$sample_id,
    symbol = cnv_oncoplot_genes
) |>
    left_join(oncoplot_data, by = c("sample_id", "symbol")) |>
    mutate(
        final_alteration = if_else(is.na(final_alteration), "", final_alteration)
    )

onco_matrix <- complete_combinations |>
    pivot_wider(
        names_from = sample_id,
        values_from = final_alteration,
        values_fill = ""
    ) |>
    column_to_rownames("symbol") |>
    as.matrix()

dim(onco_matrix)
mat <- onco_matrix[1:20, ]

sample_annotations <- all_samples |>
    arrange(FST.Group) |>
    column_to_rownames("sample_id")

## Define colors
alter_colors <- c(
    "AMP" = "#E31A1C",      # Red for amplifications
    "DEL" = "#1F78B4"      # Blue for deletions
    # "" = "white"            # White for no alteration
)

alter_fun <- list(
    "AMP" = function(x, y, w, h) {
        grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = alter_colors[["AMP"]], col = NA))
    },
    "DEL" = function(x, y, w, h) {
        grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = alter_colors[["DEL"]], col = NA))
    },
    "background" = function(x, y, w, h) {
        grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey90"))
    }
)

col_anno <- HeatmapAnnotation(
    FST_Group = sample_annotations$FST.Group,
    col = list(FST_Group = study_colors$FST.Group),
    annotation_name_side = "left",
    annotation_legend_param = list(
        FST_Group = list(title = "FST Group", title_position = "topleft")
    ),
    show_annotation_name = TRUE
)

oncoPrint(
    mat,
    alter_fun = alter_fun,
    col = alter_colors,
    
    ## Order is already set in matrix
    row_order = NULL,
    column_order = NULL,
    
    ## Annotations
    top_annotation = col_anno,
    # right_annotation = right_anno,
    
    ## Gene names
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    
    ## Remove sample names for cleaner look
    show_column_names = FALSE,
    
    ## Title
    column_title = paste0("PCGR CNV Oncoplot - DFSP Cohort (n=", total_cohort_size, " samples, ", nrow(onco_matrix), " genes)"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    ## Show percentage
    show_pct = TRUE,
    pct_side = "right"
)

## "==========================================================================="
## Find significantly altered cytobands ----
## "==========================================================================="
## Binom.test for significant cytobands
background_rate <- 0.05 # Assume 5% background alteration rate

cytoband_stats <- pcgr_cna_tbl |>
    group_by(cytoband, variant_class) |>
    summarise(
        n_samples = n_distinct(sample_id),
        n_genes = n_distinct(symbol),
        total_samples = total_cohort_size,
        frequency = n_samples / total_cohort_size,
        mean_segment_size = mean(segment_length_mb, na.rm = TRUE),
        median_segment_size = median(segment_length_mb, na.rm = TRUE),
        .groups = "drop"
    ) |>
    ## Perform binomial test for significance
    mutate(
        background_rate = case_when(
            variant_class == "amp" ~ 0.05, # 5% background for gains
            variant_class == "del" ~ 0.03, # 3% background for losses
            TRUE ~ 0.04 # Default background rate
        ),
        p_value = pmap_dbl(
            list(n_samples, total_samples, background_rate),
            function(n, total, bg) {
                binom.test(n, total, p = bg, alternative = "greater")$p.value
            }
        ),
        q_value = p.adjust(p_value, method = "BH")
    ) |>
    arrange(desc(frequency)) |> 
    ## Filter for significant cytobands
    filter(
        q_value < 0.05,           # FDR < 5%
        frequency >= 0.10,        # At least 10% frequency
        n_samples >= 5            # At least 5 samples
    ) |> 
    arrange(q_value, desc(frequency))

write_csv(
    cytoband_stats, 
    file = here(data_dir, "wes_pcgr_cnv_cytoband_stats_all_tumors.csv"), 
)

pcgr_cna_tbl |> 
    select(sample_id, cytoband, var_id)

## "==========================================================================="
## Create GISTIC-style chromosome plot with peaks ----
## "==========================================================================="

## Load additional libraries for plotting
library(ggplot2)
library(ggsci)
library(dplyr)
library(ggrepel)
library(BSgenome.Hsapiens.UCSC.hg38)

## Get chromosome information
chrom_info <- tibble(
    chromName = paste0("chr", 1:22),
    chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:22]
)
chrom_info$chromNum <- 1:22

## Calculate cumulative chromosome lengths
chrom_info$chromlengthCumsum <- cumsum(as.numeric(chrom_info$chromlength))

## Calculate chromosome start positions from 0
chrom_info$chromStartPosFrom0 <- c(0, chrom_info$chromlengthCumsum[-nrow(chrom_info)])

## Calculate chromosome middle positions from 0
tmp_middle <- diff(c(0, chrom_info$chromlengthCumsum)) / 2
chrom_info$chromMiddlePosFrom0 <- chrom_info$chromStartPosFrom0 + tmp_middle

## Extract genomic coordinates from var_id for high-resolution plotting
cytoband_genome_data <- pcgr_cna_tbl %>%
    # Split var_id by ":" to get all components: chr12:25225628-25245334:0:0
    separate(var_id, into = c("chr", "pos_range", "cn_major_parsed", "cn_minor_parsed"), 
             sep = ":", remove = FALSE) %>%
    # Split position range to get start and end
    separate(pos_range, into = c("start", "end"), sep = "-") %>%
    mutate(
        chromosome = as.numeric(gsub("chr", "", chr)),
        start = as.numeric(start),
        end = as.numeric(end),
        midpoint = (start + end) / 2
    ) %>%
    # Remove rows with parsing issues
    filter(
        !is.na(chromosome) & 
        !is.na(start) & 
        !is.na(end) & 
        chromosome <= 22
    ) %>%
    # Add cumulative positions
    left_join(chrom_info, by = c("chromosome" = "chromNum")) %>%
    mutate(
        genome_start = start + chromStartPosFrom0,
        genome_end = end + chromStartPosFrom0,
        genome_midpoint = midpoint + chromStartPosFrom0
    )

## Create high-resolution bins for plotting (100kb bins for detail)
bin_size <- 1e5  # 100KB bins for high resolution

genomic_bins <- cytoband_genome_data %>%
    mutate(
        # Create fine-scale bins
        bin_start = floor(genome_start / bin_size) * bin_size,
        bin_end = bin_start + bin_size,
        bin_midpoint = bin_start + bin_size/2
    ) %>%
    group_by(bin_start, bin_end, bin_midpoint, chromosome, variant_class) %>%
    summarise(
        n_samples_bin = n_distinct(sample_id),
        .groups = "drop"
    ) %>%
    mutate(
        proportion = n_samples_bin / total_cohort_size,
        # Make deletions negative for plotting
        proportion_signed = case_when(
            variant_class == "amp" ~ proportion,
            variant_class == "del" ~ -proportion,
            TRUE ~ 0
        )
    ) %>%
    filter(proportion > 0)  # Only keep bins with alterations

## Prepare significant peaks for labeling
significant_peaks <- cytoband_stats %>%
    # Get top peaks for labeling
    slice_max(order_by = frequency, n = 20) %>%
    mutate(
        # Extract representative positions for labels
        chromosome = as.numeric(gsub("chr(\\d+):.*", "\\1", cytoband)),
        cytoband_label = gsub("chr\\d+:", "", cytoband),
        # Create signed frequency for positioning
        freq_signed = case_when(
            variant_class == "amp" ~ frequency,
            variant_class == "del" ~ -frequency,
            TRUE ~ 0
        )
    ) %>%
    # Add genome positions
    left_join(chrom_info, by = c("chromosome" = "chromNum")) %>%
    mutate(
        # Estimate positions (use chromosome midpoint as approximation)
        genome_pos = chromMiddlePosFrom0
    )

## Create the GISTIC-style chromosome plot
gistic_plot <- ggplot() +
    # Add the alteration frequency data
    geom_col(data = genomic_bins,
             aes(x = bin_midpoint, y = proportion_signed, fill = variant_class),
             width = bin_size * 0.8, alpha = 0.8) +
    
    # Add chromosome boundary lines
    geom_vline(data = chrom_info, 
               aes(xintercept = chromlengthCumsum), 
               color = "gray30", size = 0.5, alpha = 0.8) +
    
    # Add zero line
    geom_hline(yintercept = 0, color = "black", size = 0.3) +
    
    # Add labels for significant peaks
    geom_text_repel(data = significant_peaks,
                    aes(x = genome_pos, y = freq_signed * 1.1, 
                        label = paste0(cytoband_label, " (", n_samples, ")")),
                    size = 3, 
                    max.overlaps = 15,
                    force = 2,
                    nudge_y = 0.05,
                    segment.size = 0.3,
                    segment.alpha = 0.6) +
    
    # Color scheme
    scale_fill_manual(
        values = c("amp" = "#E31A1C", "del" = "#1F78B4"),  # Red for amp, blue for del
        labels = c("amp" = "Amplifications", "del" = "Deletions"),
        name = ""
    ) +
    
    # X-axis: chromosome positions
    scale_x_continuous(
        breaks = chrom_info$chromMiddlePosFrom0,
        labels = 1:22,
        expand = c(0.01, 0.01),
        name = "Chromosomes"
    ) +
    
    # Y-axis: proportions
    scale_y_continuous(
        name = "Proportion of Copy Number Gains/Losses",
        expand = expansion(mult = c(0.05, 0.15)),  # Extra space for labels
        labels = function(x) abs(x),  # Show absolute values
        breaks = seq(-0.8, 0.8, 0.2)
    ) +
    
    # Add chromosome number labels at bottom
    geom_text(data = chrom_info,
              aes(x = chromMiddlePosFrom0, y = -0.02, label = chromNum),
              size = 3, color = "black", vjust = 1) +
    
    # Theme and styling
    theme_minimal() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "gray90", size = 0.3),
        axis.text.x = element_blank(),  # Use our custom labels
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "top",
        legend.justification = "center",
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(20, 20, 20, 20)
    ) +
    
    # Add title
    labs(
        title = paste0("Copy Number Alterations Across Genome (n=", total_cohort_size, " samples)"),
        subtitle = "Significant peaks labeled with cytoband and sample count"
    )

## Save the plot
ggsave(
    filename = file.path(plot_dir, "gistic_style_chromosome_plot.png"),
    plot = gistic_plot,
    width = 18, height = 8, dpi = 300, bg = "white"
)

## Create a summary table of top alterations
top_alterations <- cytoband_stats %>%
    arrange(desc(frequency)) %>%
    slice_head(n = 10) %>%
    select(cytoband, variant_class, frequency, n_samples, q_value) %>%
    mutate(
        frequency = round(frequency, 3),
        q_value = formatC(q_value, format = "e", digits = 2)
    )

print("Top 10 most frequent alterations:")
print(top_alterations)

## Display the plot
print(gistic_plot)
