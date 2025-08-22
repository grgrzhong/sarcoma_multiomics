#!/usr/bin/env Rscript
##############################################################################
## Description: R libraries for data analysis and visualization
## Author:  Zhong Guorui
## Created: 2025-02-17
##############################################################################

## Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(qs)
        library(fs)
        library(here)
        library(readxl)
        library(writexl)
        library(kableExtra)
        library(gridExtra)
        library(ggpubr)
        library(rstatix)
        library(ggrepel)
        library(janitor)
        library(Cairo)
        library(patchwork)
        library(cowplot)
        library(maftools)
        library(circlize)
        library(ComplexHeatmap)
        library(reshape2)
        library(sigminer)
        library(tidyverse)
    })
)

options(scipen = 999) # Disable scientific notation globally


SavePlot <- function(
    plot, width = 8, height = 6, only_png = FALSE, dir, filename) {
    ### Save png or pdf plots using ggplot2

    fs::dir_create(here(dir))

    if (only_png) {
        ### Save only png
        img_type <- ".png"
        ggsave(
            filename = here(dir, paste0(filename, img_type)),
            plot = plot,
            width = width,
            height = height,
            units = "in",
            dpi = 300,
            device = ifelse(img_type == ".png", png, cairo_pdf)
        )
        
        message(
            "Saved plot: ", here(dir, paste0(filename, img_type))
        )

    } else {
        ### Save both png and pdf
        for (img_type in c(".png", ".pdf")) {
            ggsave(
                filename = here(dir, paste0(filename, img_type)),
                plot = plot,
                width = width,
                height = height,
                units = "in",
                dpi = 300,
                device = ifelse(img_type == ".png", png, cairo_pdf)
            )
            
            message(
                "Saved plot: ", here(dir, paste0(filename, img_type))
            )
        }
    }
}

plot_theme <- function(line_width = 0.3, base_size = 8) {
    ## General theme for publication figures
    # theme_classic(base_size = 10, base_family = "Arial") +
    # theme_classic(
    #     base_size = base_size,
    #     base_family = "Arial",
    #     base_line_size = line_width
    # ) +
    theme(
            plot.title = element_text(
                face = "plain", # "bold",
                # hjust = 0.5, 
                size = base_size
            ),
            # strip.background = element_blank(),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            # legend.position = "right",
            legend.title = element_text(
                # hjust = 0,
                size = base_size - 1
                # face = "bold"
            ),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = base_size - 1),
            # legend.box = "vertical",
            # legend.box = "horizontal",
            # legend.box.just = "left", # Align legend box to the left,
            # text = element_text(size = 11),
            # axis.text = element_text(color = "black"),
            # panel.grid = element_blank(),
            axis.ticks = element_line(
                # color = "black",
                linewidth = line_width
            ),
            # axis.line = element_line(color = "black", linewidth = line_width),
            # axis.line = element_blank(),
            # panel.border = element_rect(
            #     linewidth = line_width, color = "black", fill = NA
            # ),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")
        )
}

figure_theme <- function(line_width = 0.3, base_size = 8) {
    ## General theme for publication figures
    # theme_classic(base_size = 10, base_family = "Arial") +
    theme_classic(
        base_size = base_size,
        base_family = "Arial",
        base_line_size = line_width
    ) +
        theme(
            plot.title = element_text(
                face = "plain", # "bold",
                # hjust = 0.5, 
                size = base_size
            ),
            strip.background = element_blank(),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "right",
            legend.title = element_text(
                hjust = 0,
                size = base_size - 1
                # face = "bold"
            ),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = base_size - 1),
            legend.box = "vertical",
            # legend.box = "horizontal",
            legend.box.just = "left", # Align legend box to the left,
            # text = element_text(size = 11),
            axis.text = element_text(color = "black"),
            # panel.grid = element_blank(),
            axis.ticks = element_line(
                linewidth = line_width, color = "black"
            ),
            # axis.line = element_line(color = "black", linewidth = line_width),
            # axis.line = element_blank(),
            # panel.border = element_rect(
            #     linewidth = line_width, color = "black", fill = NA
            # ),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")
        )
}

figure_theme2 <- function() {
    ## General theme for publication figures
    # theme_classic(base_size = 10, base_family = "Arial") +
    line_width <- 0.3
    base_size <- 8
    theme_bw(
        base_size = base_size,
        base_family = "Arial",
        base_line_size = line_width
    ) +
        theme(
            plot.title = element_text(face = "plain", hjust = 0.5, size = 8),
            # strip.background = element_blank(),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "right",
            legend.title = element_text(
                hjust = 0,
                size = 8
                # face = "bold"
            ),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = base_size - 1),
            legend.box = "vertical",
            # legend.box = "horizontal",
            legend.box.just = "left", # Align legend box to the left,
            # text = element_text(size = 11),
            axis.text = element_text(color = "black"),
            # panel.grid = element_blank(),
            axis.ticks = element_line(
                linewidth = line_width, # color = "black"
            ),
            # axis.line = element_line(color = "black", linewidth = line_width),
            # axis.line = element_blank(),
            # panel.border = element_rect(
            #     linewidth = line_width, color = "black", fill = NA
            # ),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")
        )
}

plot_theme <- function() {
    # theme_classic(base_size = 10, base_family = "Arial") +
    # theme_classic() +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
        # strip.background = element_blank(),
        plot.background = element_rect(
            fill = "white", # Set plot background to white
            color = NA # No border
        ),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top",
        legend.title = element_text(
            hjust = 0,
            size = 8
            # face = "bold"
        ),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8),
        legend.box = "vertical",
        # legend.box = "horizontal",
        legend.box.just = "left", # Align legend box to the left,
        text = element_text(size = 11, family = "Arial"),
        axis.text = element_text(color = "black"),
        # axis.ticks = element_line(linewidth = 0.3, color = "black"),    # Set the thickness of axis ticks
        # axis.line = element_line(color = "black", linewidth = 0.5),
        # axis.line = element_blank(),
        # axis.line = element_line(linewidth = 0.5, color = "black"),     # Set the thickness of axis lines
        # panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), units = "cm")
    )
}

testVariantFilter <- function(
    filter_params,
    maf) {
    # Apply filtering with parameters
    maf_filtered <- maf |>
        ## Filter by sequencing quality
        filter(
            DP >= filter_params$min_depth | is.na(DP),
            (is.na(AF) | as.numeric(AF) >= filter_params$min_vaf),
            (is.na(VAD) | as.numeric(VAD) >= filter_params$min_vad),
        ) |>
        ## Filter by variant impact
        filter(
            !Variant_Classification %in% filter_params$exclude_classifications,
            !ExonicFunc.refGene %in% filter_params$exclude_exonic_func,
            (Func.refGene %in% filter_params$include_func | is.na(Func.refGene))
        ) |>
        ## Filter population variants
        mutate(
            gnomAD_exome_ALL_num = suppressWarnings(
                as.numeric(gnomAD_exome_ALL)
            )
        ) |>
        filter(
            is.na(gnomAD_exome_ALL_num) |
                gnomAD_exome_ALL_num <= filter_params$max_population_freq
        ) |>
        select(-gnomAD_exome_ALL_num)

    # Count variants per sample after filtering
    post_filter_count <- maf_filtered |>
        count(Tumor_Sample_Barcode, name = "n_variants_after")

    # Merge before and after counts
    filter_res <- post_filter_count |>
        mutate(
            n_variants_after = replace_na(n_variants_after, 0),
            # percent_retained = round(n_variants_after / n_variants_before * 100, 1)
        )

    list(
        maf_filtered = maf_filtered,
        filter_res = filter_res
    )
}

plotVariantFilter <- function(
    is_save = FALSE,
    fig_dir = "figures/variant_filtering",
    fig_name = "filter_params") {
    # Create comparison visualizations
    ## 1. Bar plot of variant counts before/after by sample
    p <- filter_res |>
        pivot_longer(
            cols = c(n_variants_before, n_variants_after),
            names_to = "filter_status",
            values_to = "variant_count"
        ) |>
        mutate(
            filter_status = ifelse(
                filter_status == "n_variants_before",
                "Before filtering", "After filtering"
            ),
            filter_status = factor(filter_status, levels = c("Before filtering", "After filtering"))
        ) |>
        ggplot(aes(x = variant_count, fill = filter_status)) +
        geom_histogram(
            binwidth = 200, alpha = 0.7, position = "identity", boundary = 0
        ) +
        facet_wrap(~filter_status, ncol = 1, scales = "free_y") +
        scale_fill_manual(
            values = c("Before filtering" = "#1f77b4", "After filtering" = "#ff7f0e")
        ) +
        # Add mean lines with labels
        geom_vline(
            data = . |> group_by(filter_status) |>
                summarize(mean_count = mean(variant_count)),
            aes(xintercept = mean_count),
            linetype = "dashed", linewidth = 0.5, color = "red"
        ) +
        # Add mean value labels
        geom_text(
            data = . |> group_by(filter_status) |>
                summarize(mean_count = mean(variant_count)),
            aes(x = mean_count, y = Inf, label = sprintf("Mean: %d", round(mean_count))),
            hjust = -0.1, vjust = 1.5, color = "black", size = 7, size.unit = "pt"
        ) +
        # Add filter parameters as annotation only on "After filtering" panel
        geom_text(
            data = data.frame(
                filter_status = "Before filtering",
                x = 20000, # Adjust this position as needed
                y = I(0.3) # Use relative position on y-axis
            ),
            aes(x = x, y = y, label = sprintf(
                "Filter params:\nread depth >= %d\nvaf >= %.2f\nvad >= %d\npop_freq <= %.4f",
                filter_params$min_depth,
                filter_params$min_vaf,
                filter_params$min_vad,
                filter_params$max_population_freq
            )),
            hjust = 0, # Left-aligned text
            vjust = 0, # Top-aligned text
            size = 6,
            size.unit = "pt",
            color = "black",
            fontface = "plain"
        ) +
        # Colors and theme
        labs(
            title = "Distribution of variant counts per sample",
            x = "Number of variants",
            y = "Number of samples",
            fill = "Filter status"
        ) +
        figure_theme() +
        theme(
            # strip.text = element_text(size = 12, face = "bold"),
            # axis.title = element_text(size = 11),
            legend.position = "none"
        )

    if (is_save) {
        SavePlot(
            plot = p,
            width = 4,
            height = 3,
            dir = fig_dir,
            filename = fig_name
        )
    }
}

mafOncoPlot <- function(
    maf,
    genes = NULL,
    top_n_genes = 30,
    clinicalFeatures = NULL,
    annotationColor = NULL,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = TRUE,
    removeNonMutated = FALSE,
    titleText = "OncoPlot",
    fontSize = 0.7,
    width = 10,
    height = 8,
    fig_dir = "figures/oncoplot",
    fig_name = "oncoplot") {
    # Create directory if it doesn't exist
    fs::dir_create(here(fig_dir))

    # If no genes specified, get top mutated ones
    if (is.null(genes)) {
        gene_summary <- getGeneSummary(maf)
        n_genes <- min(top_n_genes, nrow(gene_summary))

        # Check if we have enough genes
        if (n_genes < 2) {
            warning("Not enough genes for oncoplot (minimum 2 required)")
            return(NULL)
        }

        genes <- gene_summary[1:n_genes, "Hugo_Symbol"]
        genes <- genes$Hugo_Symbol
    }

    # Save PDF version
    pdf_file <- here(fig_dir, paste0(fig_name, ".pdf"))
    pdf(pdf_file, width = width, height = height)

    p <- oncoplot(
        maf = maf,
        genes = genes,
        clinicalFeatures = clinicalFeatures,
        annotationColor = annotationColor,
        sortByAnnotation = sortByAnnotation,
        showTumorSampleBarcodes = showTumorSampleBarcodes,
        removeNonMutated = removeNonMutated,
        fontSize = fontSize,
        titleText = titleText
    )

    dev.off()
    message(paste("Saved plots to:", pdf_file))

    # save PNG version
    png_file <- here(fig_dir, paste0(fig_name, ".png"))
    png(png_file, width = width, height = height, res = 300, units = "in")

    oncoplot(
        maf = maf,
        genes = genes,
        clinicalFeatures = clinicalFeatures,
        annotationColor = annotationColor,
        sortByAnnotation = sortByAnnotation,
        showTumorSampleBarcodes = showTumorSampleBarcodes,
        removeNonMutated = removeNonMutated,
        fontSize = fontSize,
        titleText = titleText
    )
    dev.off()

    message(paste("Saved plots to:", png_file))

    print(p)
}

LoadSampleInfo <- function() {
    read_xlsx(
        here(
            "data/clinical/DFSP_multiomics_sample_list_updated_20250509.xlsx"
        )
    )
}

LoadOncoKBGeneList <- function() {
    file_path <- here("data/public/cancerGeneList.tsv")

    if (!file_exists(file_path)) {
        stop("OncoKB gene list file not found: ", file_path)
    }

    data <- read_tsv(file_path, show_col_types = FALSE)

    data |> pull(`Hugo Symbol`)
}

addCancerHotspot <- function(
    maf,                           # maf data in dataframe format
    hotspot = NULL,                # path to the hotspot file
    qvalue = NULL,                 # qvalue threshold
    median_allele_freq_rank = NULL, # median Allele Frequency Rank threshold 
    log10_pvalue = NULL           # log10 pvalue threshold

) {
    if (is.null(hotspot)) {

        snv_hotspots <- read_xlsx(
            here("data/clinical/hotspots_v2.xlsx"),
            sheet = "SNV-hotspots"
        )
        

        indel_hotspots <- read_xlsx(
            here("data/clinical/hotspots_v2.xlsx"),
            sheet = "INDEL-hotspots"
        ) 
        
    } else {
        
        snv_hotspots <- read_xlsx(
            here(hotspot),
            sheet = "SNV-hotspots"
        )
        

        indel_hotspots <- read_xlsx(
            here(hotspot),
            sheet = "INDEL-hotspots"
        )
    }

    snv_hotspots <- snv_hotspots |>
        select(
            Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, qvalue, Median_Allele_Freq_Rank
        ) |>
        mutate(
            ref_aa = str_extract(Reference_Amino_Acid, "^[A-Z*]"),
            pos_aa = as.character(Amino_Acid_Position),
            var_aa = str_extract(Variant_Amino_Acid, "^[A-Z*]"),
            aaChange = paste0("p.", ref_aa, pos_aa, var_aa)
        ) |>
        distinct() |>
        mutate(
            snv_hotspot = paste(
                Hugo_Symbol, aaChange,
                sep = "_"
            )
        )
    
    indel_hotspots <- indel_hotspots |>
        select(
            Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, qvalue, Median_Allele_Freq_Rank
        ) |>
        mutate(
            aaChange = str_extract(Variant_Amino_Acid, "^[^:]+")
        ) |>
        mutate(
            aaChange = paste0("p.", aaChange)
        ) |>
        distinct() |>
        mutate(
            indel_hotspot = paste0(
                Hugo_Symbol, "_", aaChange
            )
        )
    
    ## Match the variant info from the aaChange column in maf and
    ## Variant_Amino_Acid column in the hotspot data
    hotspots <- list(
        snv_hotspots = snv_hotspots,
        indel_hotspots = indel_hotspots
    )

    message(
        sprintf(
            "Total SNV hotspots: %d", nrow(hotspots$snv_hotspots)
        ), "\n",
        sprintf(
            "Total INDEL hotspots: %d", nrow(hotspots$indel_hotspots)
        )
    )

    ## Apply filtering only if filter parameters are provided
    if (
        !is.null(qvalue) || !is.null(median_allele_freq_rank) ||
        !is.null(log10_pvalue)
    ) {
        ## build filter conditions
        filter_expr <- list()
    
        if (!is.null(qvalue)) {
            
            filter_expr <- append(
                filter_expr, rlang::expr(qvalue < !!qvalue)
            )
        }
        
        if (!is.null(median_allele_freq_rank)) {
            filter_expr <- append(
                filter_expr, 
                rlang::expr(Median_Allele_Freq_Rank > !!median_allele_freq_rank)
            )
        }
        
        if (!is.null(log10_pvalue)) {
            
            filter_expr <- append(
                filter_expr, rlang::expr(log10_pvalue > !!log10_pvalue)
            )
        }

        ## Apply filtering to select the high confidence hotspot variants
        hotspots <- map(
            hotspots,
            ~ .x |> filter(!!!filter_expr)
        )

        message("Applied filters:")
        if (!is.null(qvalue)) message("  - qvalue < ", qvalue)
        if (!is.null(median_allele_freq_rank)) message("  - Median_Allele_Freq_Rank > ", median_allele_freq_rank)
        if (!is.null(log10_pvalue)) message("  - log10_pvalue > ", log10_pvalue)

        message(
            sprintf(
                "After filtering SNV hotspots: %d", nrow(hotspots$snv_hotspots)
            ), "\n",
            sprintf(
                "After filtering INDEL hotspots: %d", nrow(hotspots$indel_hotspots)
            )
        )

    } else {

        message("No filters applied - using all hotspots")
    }

    all_hotspots <- c(
        hotspots$snv_hotspots |> pull(snv_hotspot),
        hotspots$indel_hotspots |> pull(indel_hotspot)
    )

    maf_hotspot <- maf |>
        # select(
        #     Tumor_Sample_Barcode, Variant_Classification, Hugo_Symbol, aaChange
        # ) |>
        # filter(
        #     # Variant_Classification == "Missense_Mutation"
        #     Variant_Classification == "Nonsense_Mutation"
        # ) |>
        mutate(match_change = paste0(Hugo_Symbol, "_", aaChange)) |>
        mutate(
            is_hotspot = if_else(
                match_change %in% all_hotspots,
                TRUE,
                FALSE
            )
        ) |>
        arrange(desc(is_hotspot)) |>
        select(-match_change) |>
        relocate(
            is_hotspot,
            .after = aaChange
        )
    
    is_hotspot_data <- maf_hotspot |>
        filter(is_hotspot) |>
        select(
            Tumor_Sample_Barcode, Variant_Classification, Hugo_Symbol, aaChange,
            is_hotspot
        )
    
    message(
        sprintf(
            "Number of hotspot variants in the maf: %d", nrow(is_hotspot_data)
        )
    )

    return(maf_hotspot)
}

LoadCancerHotspot <- function(
    hotspot = NULL,                # path to the hotspot file
    qvalue = NULL,                 # qvalue threshold
    median_allele_freq_rank = NULL, # median Allele Frequency Rank threshold 
    log10_pvalue = NULL           # log10 pvalue threshold

) {
    if (is.null(hotspot)) {

        snv_hotspots <- read_xlsx(
            here("data/public/hotspots_v2.xlsx"),
            sheet = "SNV-hotspots"
        )
        

        indel_hotspots <- read_xlsx(
            here("data/public/hotspots_v2.xlsx"),
            sheet = "INDEL-hotspots"
        ) 
        
    } else {
        
        snv_hotspots <- read_xlsx(
            here(hotspot),
            sheet = "SNV-hotspots"
        )
        

        indel_hotspots <- read_xlsx(
            here(hotspot),
            sheet = "INDEL-hotspots"
        )
    }

    snv_hotspots <- snv_hotspots |>
        select(
            Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, qvalue, Median_Allele_Freq_Rank
        ) |>
        mutate(
            ref_aa = str_extract(Reference_Amino_Acid, "^[A-Z*]"),
            pos_aa = as.character(Amino_Acid_Position),
            var_aa = str_extract(Variant_Amino_Acid, "^[A-Z*]"),
            aaChange = paste0("p.", ref_aa, pos_aa, var_aa)
        ) |>
        distinct() |>
        mutate(
            snv_hotspot = paste(
                Hugo_Symbol, aaChange,
                sep = "_"
            )
        )
    
    indel_hotspots <- indel_hotspots |>
        select(
            Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, qvalue, Median_Allele_Freq_Rank
        ) |>
        mutate(
            aaChange = str_extract(Variant_Amino_Acid, "^[^:]+")
        ) |>
        mutate(
            aaChange = paste0("p.", aaChange)
        ) |>
        distinct() |>
        mutate(
            indel_hotspot = paste0(
                Hugo_Symbol, "_", aaChange
            )
        )
    
    ## Match the variant info from the aaChange column in maf and
    ## Variant_Amino_Acid column in the hotspot data
    hotspots <- list(
        snv_hotspots = snv_hotspots,
        indel_hotspots = indel_hotspots
    )

    message(
        sprintf(
            "Total SNV hotspots: %d", nrow(hotspots$snv_hotspots)
        ), "\n",
        sprintf(
            "Total INDEL hotspots: %d", nrow(hotspots$indel_hotspots)
        )
    )

    ## Apply filtering only if filter parameters are provided
    if (
        !is.null(qvalue) || !is.null(median_allele_freq_rank) ||
        !is.null(log10_pvalue)
    ) {
        ## build filter conditions
        filter_expr <- list()
    
        if (!is.null(qvalue)) {
            
            filter_expr <- append(
                filter_expr, rlang::expr(qvalue < !!qvalue)
            )
        }
        
        if (!is.null(median_allele_freq_rank)) {
            filter_expr <- append(
                filter_expr, 
                rlang::expr(Median_Allele_Freq_Rank > !!median_allele_freq_rank)
            )
        }
        
        if (!is.null(log10_pvalue)) {
            
            filter_expr <- append(
                filter_expr, rlang::expr(log10_pvalue > !!log10_pvalue)
            )
        }

        ## Apply filtering to select the high confidence hotspot variants
        hotspots <- map(
            hotspots,
            ~ .x |> filter(!!!filter_expr)
        )

        message("Applied filters:")
        if (!is.null(qvalue)) message("  - qvalue < ", qvalue)
        if (!is.null(median_allele_freq_rank)) message("  - Median_Allele_Freq_Rank > ", median_allele_freq_rank)
        if (!is.null(log10_pvalue)) message("  - log10_pvalue > ", log10_pvalue)

        message(
            sprintf(
                "After filtering SNV hotspots: %d", nrow(hotspots$snv_hotspots)
            ), "\n",
            sprintf(
                "After filtering INDEL hotspots: %d", nrow(hotspots$indel_hotspots)
            )
        )

    } else {

        message("No filters applied - using all hotspots")
    }

    hotspots
}


CollectAnnovarData <- function(
    dir,
    is_save = FALSE,
    save_dir = "data/wes/annotation/merged"
) {
    
    ## The annovar directory should contain the annovar output files
    dir <- here(dir)

    input_files <- dir_ls(dir, recurse = TRUE, glob = "*annovar.txt")
    
    ## Collect all annovar output files
    message(
        "Collecting ", length(input_files),
        " annovar output files from: ", dir
    )
    
    maf_data <- annovarToMaf(input_files) |>
        as_tibble() |>
        mutate(
            ## If the tumour sample has a matched normal sample, Otherinfo13 = normal allell data, Otherinfo14 = tumour allell data
            ## if the tumour sample does not have a matched normal sample, Otherinfo13 = tumour allell data, Otherinfo14 = NA
            is_paired_normal = if_else(
                is.na(Otherinfo14) | Otherinfo14 == "NA" | Otherinfo14 == "",
                FALSE,
                TRUE
            )
        )
    
    # maf_data |> 
    #     select(
    #         Tumor_Sample_Barcode, Hugo_Symbol, aaChange, 
    #         Otherinfo13, Otherinfo14, is_paired_normal,
    #         AD, AF, DP
    #     )

    ## Extract the allele data for paired samples
    maf_data_paired <- maf_data |> filter(is_paired_normal)
    maf_data_unpaired <- maf_data |> filter(!is_paired_normal)

    n_paired_samples <- length(unique(maf_data_paired$Tumor_Sample_Barcode))
    n_unpaired_samples <- length(unique(maf_data_unpaired$Tumor_Sample_Barcode))

    if (n_paired_samples > 0) {

        message(
            sprintf(
                "Extracting allele data for tumor-normal sample: %d", 
                n_paired_samples
            )
        )
        
        maf_data_paired <- maf_data_paired |> 
            ## Extract the tumour allele data
            mutate(
                tumor_AD = sapply(
                    str_split(Otherinfo14, ":"),
                    function(x) x[2]
                )
            ) |> 
            separate(
                tumor_AD, 
                into = c("tumor_RAD", "tumor_VAD"), 
                sep = ",", 
                remove = TRUE,

            ) |> 
            mutate(
                tumor_RAD = as.numeric(tumor_RAD),
                tumor_VAD = as.numeric(tumor_VAD)
            ) |> 
            mutate(
                tumor_AF = as.numeric(
                        sapply(
                        str_split(Otherinfo14, ":"), 
                        function(x) x[3]
                    )
                ),
                tumor_DP = as.numeric(
                        sapply(
                        str_split(Otherinfo14, ":"), 
                        function(x) x[4]
                    )
                )
            ) |> 
            ## Extract the normal allele data
            mutate(
                normal_AD = sapply(
                    str_split(Otherinfo13, ":"),
                    function(x) x[2]
                )
            ) |> 
            separate(
                normal_AD, 
                into = c("normal_RAD", "normal_VAD"), 
                sep = ",", 
                remove = TRUE,

            ) |> 
            mutate(
                normal_RAD = as.numeric(normal_RAD),
                normal_VAD = as.numeric(normal_VAD)
            ) |> 
            mutate(
                normal_AF = as.numeric(
                        sapply(
                        str_split(Otherinfo13, ":"), 
                        function(x) x[3]
                    )
                ),
                normal_DP = as.numeric(
                        sapply(
                        str_split(Otherinfo13, ":"), 
                        function(x) x[4]
                    )
                )
            )

    } else {

        maf_data_paired <- maf_data_paired |> 
            mutate(
                tumor_RAD = NA_real_,
                tumor_VAD = NA_real_,
                tumor_AF = NA_real_,
                tumor_DP = NA_real_,
                normal_RAD = NA_real_,
                normal_VAD = NA_real_,
                normal_AF = NA_real_,
                normal_DP = NA_real_
            )
    }

    if (n_unpaired_samples > 0) {

        message(
            sprintf(
                "Extracting allele data for tumor-only sample: %d", 
                n_unpaired_samples
            )
        )

        maf_data_unpaired <- maf_data_unpaired |> 
            mutate(
                tumor_AD = sapply(
                    str_split(Otherinfo13, ":"),
                    function(x) x[2]
                )
            ) |> 
            separate(
                tumor_AD, 
                into = c("tumor_RAD", "tumor_VAD"), 
                sep = ",", 
                remove = TRUE,

            ) |> 
            mutate(
                tumor_RAD = as.numeric(tumor_RAD),
                tumor_VAD = as.numeric(tumor_VAD)
            ) |>
            mutate(
                tumor_AF = as.numeric(
                        sapply(
                        str_split(Otherinfo13, ":"), 
                        function(x) x[3]
                    )
                ),
                tumor_DP = as.numeric(
                        sapply(
                        str_split(Otherinfo13, ":"), 
                        function(x) x[4]
                    )
                )
            ) |> 
            mutate(
                normal_RAD = NA_real_,
                normal_VAD = NA_real_,
                normal_AF = NA_real_,
                normal_DP = NA_real_
            )
    } else {
        
        maf_data_unpaired <- maf_data_unpaired |> 
            mutate(
                tumor_RAD = NA_real_,
                tumor_VAD = NA_real_,
                tumor_AF = NA_real_,
                tumor_DP = NA_real_,
                normal_RAD = NA_real_,
                normal_VAD = NA_real_,
                normal_AF = NA_real_,
                normal_DP = NA_real_
            )
    }

    ## Combine the paired and unpaired samples
    all_columns <- colnames(maf_data_paired)
    
    maf_data_unpaired <- maf_data_unpaired |> 
        select(all_of(all_columns))
    
    stopifnot(
        all.equal(
            colnames(maf_data_paired),
            colnames(maf_data_unpaired)
        )
    )

    maf_tbl <- bind_rows(
        maf_data_paired,
        maf_data_unpaired
    ) |>
        mutate(
            gnomAD_exome_ALL = as.numeric(
                replace(gnomAD_exome_ALL, gnomAD_exome_ALL == ".", NA)
            )
        ) |> 
        ## reange the columns
        relocate(
            is_paired_normal, 
            Otherinfo13, Otherinfo14,
            AD, AF, DP, 
            tumor_RAD, tumor_VAD, tumor_AF, tumor_DP,
            normal_RAD, normal_VAD, normal_AF, normal_DP,
            .after = last_col()
        )

    # ## Save the data
    # if (is_save && !is.null(save_dir)) {
    #     # Save the merged maf file
    #     file_name <- "annovar_maf_merged.qs"
        
    #     message(
    #         "Saving merged annovar data: ", here(save_dir, file_name)
    #     )

    #     fs::dir_create(here(save_dir))
        
    #     qsave(maf_tbl, here(save_dir, file_name))
    # }
    
    # Return a tibble
    maf_tbl

}

FilterAnnovarData <- function(data) {

    annovar_tbl <- data
    ## VAF >= 5% and VAD >= 4 and DP >=20 in the tumor samples
    filter_data1 <- annovar_tbl |>
        filter(
            is.na(tumor_AF) | tumor_AF >= 0.05,
            is.na(tumor_DP) | tumor_DP >= 20,
            is.na(tumor_VAD) | tumor_VAD >= 4,
            is.na(gnomAD_exome_ALL) | gnomAD_exome_ALL < 0.001
        )

    # filter_data1 |> select(matches("tumor|normal|AF|VAD|DP"))

    ## VAF <=1% and VAD <=1 in the paired normal samples
    filter_data2 <- annovar_tbl |>
        filter(!is_paired_normal) |>
        filter(
            is.na(normal_AF) | normal_AF <= 0.01,
            is.na(normal_VAD) | normal_VAD <= 1
        )

    # filter_data2 |> select(matches("tumor|normal|AF|VAD|DP"))

    ## Specific for tert promotor mutations
    tert_gene <- filter_data1 |>
        filter(Hugo_Symbol == "TERT") |>
        filter(Func.refGene %in% c("exonic", "splicing", "upstream"))

    ## Combine the filtered data
    first_inclusion <- bind_rows(
        filter_data1,
        filter_data2,
        tert_gene
    ) |>
        distinct()

    ## All indel/splice site mutations
    unique(first_inclusion$Variant_Classification)
    unique(first_inclusion$ExonicFunc.refGene)
    unique(first_inclusion$Variant_Type)
    unique(first_inclusion$Func.refGene)

    first_inclusion_splice <- first_inclusion |>
        filter(Func.refGene %in% c("splicing"))

    first_inclusion_indels <- first_inclusion |>
        ## "SNP" "DEL" "INS" "DNP" "ONP"
        filter(Variant_Type %in% c("DEL", "INS")) |>
        ## "Missense_Mutation" "Silent" "In_Frame_Del" "Frame_Shift_Del" NA
        ## "Frame_Shift_Ins" "Nonsense_Mutation" "In_Frame_Ins" "Nonstop_Mutation"
        ## "Translation_Start_Site" "Unknown" "Inframe_INDEL"
        filter(
            Variant_Classification %in% c(
                "In_Frame_Del",
                "In_Frame_Ins",
                "Inframe_INDEL",
                "Frame_Shift_Del",
                "Frame_Shift_Ins",
                "Nonsense_Mutation",
                "Nonstop_Mutation",
                "Unknown",
                NA,
                "NA",
                "Translation_Start_Site"
            )
        )

    message(
        "Number of indels/splice site mutations: ",
        nrow(first_inclusion_indels) + nrow(first_inclusion_splice)
    )

    ##############################################################################
    ## SNV filtering  ----------------
    ##############################################################################
    ## Exclude the synonymous SNVs
    first_inclusion_snv <- first_inclusion |>
        ## >=2 base substitution also include
        filter(Variant_Type %in% c("SNP", "DNP", "ONP")) |>
        filter(!(ExonicFunc.refGene %in% c("synonymous SNV")))

    ## Deleteriousness functional impact or pathogenicity predictions
    ## Check if variant meets at least 3 deleterious functional impact or pathogenic predictions
    ## First options
    second_inclusion_snv1 <- first_inclusion_snv |>
        rowwise() |>
        mutate(
            deleterious_count = sum(
                # CADD >= 20
                (!is.na(CADD_phred) && CADD_phred >= 20),
                # VEST3 criteria
                (!is.na(VEST3_score) && VEST3_score >= 0.7) || (!is.na(VEST3_rankscore) && VEST3_rankscore >= 0.9),
                # DANN >= 0.9
                (!is.na(DANN_score) && DANN_score >= 0.9),
                # SIFT prediction = "D" (Deleterious)
                (!is.na(SIFT_pred) && SIFT_pred == "D"),
                # Polyphen2 prediction = "P"/"D" (Possibly/Probably damaging)
                (!is.na(Polyphen2_HVAR_pred) && Polyphen2_HVAR_pred %in% c("P", "D")) ||
                    (!is.na(Polyphen2_HDIV_pred) && Polyphen2_HDIV_pred %in% c("P", "D")),
                # MutationTaster prediction = "A"/"D" (Disease causing automatic/Disease causing)
                (!is.na(MutationTaster_pred) && MutationTaster_pred %in% c("A", "D")),
                # NA in at least one prediction algorithm qualifies
                is.na(CADD_phred) || is.na(VEST3_score) || is.na(DANN_score) ||
                    is.na(SIFT_pred) || is.na(Polyphen2_HVAR_pred) || is.na(Polyphen2_HDIV_pred) ||
                    is.na(MutationTaster_pred)
            )
        ) |>
        ungroup() |>
        filter(deleterious_count >= 3)

    message(
        "Number of SNVs with at least 3 deleterious functional impact or pathogenic predictions: ",
        nrow(second_inclusion_snv1)
    )

    ## Second options
    ## CLNSIG, cancer hotspot, oncoKB, COSMIC
    ## Filter by CLNSIG, cancer hotspot, oncoKB, COSMIC
    second_inclusion_snv_clnsig <- first_inclusion_snv |>
        filter(
            CLNSIG %in% c(
                "Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"
            )
        )

    ## COSMIC
    # unique(first_inclusion_snv$cosmic70)
    second_inclusion_snv_cosmic <- first_inclusion_snv |>
        filter(!(cosmic70 %in% c(".")))

    ## Cancer hotspot
    snv_hotspots <- LoadCancerHotspot()[["snv_hotspots"]]

    dfsp_snv_hotsopts <- first_inclusion_snv |>
        mutate(match_change = paste0(Hugo_Symbol, "_", aaChange)) |>
        mutate(
            is_hotspot = if_else(
                match_change %in% snv_hotspots$snv_hotspot,
                TRUE,
                FALSE
            )
        ) |>
        filter(is_hotspot)

    second_inclusion_snv_hotspot <- dfsp_snv_hotsopts |>
        select(-match_change, -is_hotspot)

    message(
        "Number of SNVs in cancer hotspots: ",
        nrow(second_inclusion_snv_hotspot)
    )

    ## OncoKB
    oncokb_gene_list <- "data/public/cancerGeneList.tsv"

    oncokb_gene_list <- read_tsv(
        here(oncokb_gene_list),
        show_col_types = FALSE
    ) |>
        pull(`Hugo Symbol`)

    second_inclusion_snv_oncokb <- first_inclusion_snv |>
        mutate(
            is_oncokb = if_else(
                Hugo_Symbol %in% oncokb_gene_list,
                TRUE,
                FALSE
            )
        ) |>
        filter(is_oncokb) |>
        select(-is_oncokb)

    message(
        "Number of variants with oncoKB ",
        nrow(second_inclusion_snv_oncokb)
    )

    snv_oncokb_uniq <- second_inclusion_snv_oncokb |>
        select(Hugo_Symbol, aaChange) |>
        distinct()

    message(
        "Number of unique OncoKB SNVs: ",
        nrow(snv_oncokb_uniq)
    )

    # Combine all included variants
    all_included_variants <- bind_rows(
        first_inclusion_splice,
        first_inclusion_indels,
        second_inclusion_snv1,
        second_inclusion_snv_clnsig,
        second_inclusion_snv_cosmic,
        second_inclusion_snv_hotspot,
        second_inclusion_snv_oncokb
    ) |>
        distinct()

    message(
        "Total number of filtered variants: ",
        nrow(annovar_tbl) - nrow(all_included_variants)
    )

    message(
        "Total number of included variants: ",
        nrow(all_included_variants)
    )
    
    all_included_variants
}


filterMergedMaf <- function(
    maf_tbl,
    # Filtering parameters by selected maf columns
    filter_params = list(
        # Minimum read depth
        DP = list(op = ">=", value = 8),
        # Minimum variant allele depth
        VAD = list(op = ">=", value = 4),
        # Minimum variant allele frequency
        AF = list(op = ">=", value = 0.05),
        # Maximum population frequency
        gnomAD_exome_ALL = list(op = "<=", value = 0.01),
        
        # Exclude these classifications
        Variant_Classification = list(op = "out", value = c("Silent")),
        
        # Exclude these functional annotations
        Func.refGene = list(op = "out", value = c("synonymous_SNV"))
    ),
    print_summary = TRUE
) {

    # Summary variants before filtering
    summary_before <- maf_tbl |>
        group_by(Tumor_Sample_Barcode) |>
        summarise(n_variants = n()) |>
        summarise(
            across(
                n_variants,
                list(
                    total = ~ sum(.x),
                    mean = ~ mean(.x, na.rm = TRUE),
                    median = ~ median(.x, na.rm = TRUE)
                ),
                .names = "{.fn}_variants_before"
            )
        )

    # # Glimplese the maf data
    # colnames(maf_tbl)
    # unique(maf_tbl$Variant_Classification)
    # unique(maf_tbl$ExonicFunc.refGene)
    # unique(maf_tbl$Func.refGene)

    # ## Functional prediction
    # unique(maf_tbl$SIFT_pred)
    # unique(maf_tbl$Polyphen2_HDIV_pred)
    # unique(maf_tbl$MutationTaster_pred)
    # unique(maf_tbl$MutationAssessor_pred)
    # unique(maf_tbl$FATHMM_pred)
    # unique(maf_tbl$PROVEAN_pred)
    # unique(maf_tbl$MetaSVM_pred)
    # unique(maf_tbl$`M-CAP_pred`)

    maf_filter <- maf_tbl

    # Dynamically apply filters based on the provided filter_params
    for (param in names(filter_params)) {
        
        filter_rule <- filter_params[[param]]
        
        op <- filter_rule$op
        
        value <- filter_rule$value
        
        # Apply numeric filters
        if (is.numeric(value)) {

            if (op == ">=") {

                maf_filter <- maf_filter |> 
                    dplyr::filter(
                        is.na(.data[[param]]) | .data[[param]] >= value
                    )

            } else if (op == "<=") {

                maf_filter <- maf_filter |> 
                    dplyr::filter(
                        is.na(.data[[param]]) | .data[[param]] <= value
                    )
            }
        }
        
        # Apply categorical filters
        if (is.character(value) || is.factor(value)) {
            
            if (op == "in") {
                
                maf_filter <- maf_filter |> 
                    dplyr::filter(.data[[param]] %in% value)

            } else if (op == "out") {

                maf_filter <- maf_filter |> 
                    dplyr::filter(!(.data[[param]] %in% value))
            }
        }
    }
    
    # maf_filtered <- maf_tbl |> 
    #     ## Filter out variants with low read depth, variant allele depth,
    #     dplyr::filter(
    #         is.na(DP) | DP >= min_rd,
    #         is.na(VAD) | as.numeric(VAD) >= min_vad,
    #         is.na(AF) | as.numeric(AF) >= min_vaf
    #     ) |> 
    #     ## Filter out variants with high population frequency
    #     dplyr::filter(
    #         is.na(gnomAD_exome_ALL) | gnomAD_exome_ALL <= max_pop_freq
    #     ) |> 
    #     arrange(gnomAD_exome_ALL) |> 
    #     ## Filter out silent (synonymous) variants
    #     dplyr::filter(
    #         (Variant_Classification %in% c("Silent")),
    #         (Func.refGene %in% c("synonymous_SNV")),
    #         (ExonicFunc.refGene %in% c("synonymous SNV"))

    #     ) |> 
    #     ## Inlcude Filter out non-exonic variants
    #     dplyr::filter(
    #         Func.refGene %in% c(
    #             "exonic", "splicing", "UTR5"
    #         )
    #     )

    summary_after <- maf_filter |>
        group_by(Tumor_Sample_Barcode) |>
        summarise(n_variants = n()) |>
        summarise(
            across(
                n_variants,
                list(
                    total = ~ sum(.x),
                    mean = ~ mean(.x, na.rm = TRUE),
                    median = ~ median(.x, na.rm = TRUE)
                ),
                .names = "{.fn}_variants_after"
            )
        )

    # Print summary if requested
    if (print_summary) {
        message("Variant filtering summary:")
        message("---------------------------------------------------")
        message(
            sprintf(
                "Total variants before filtering: %d",
                summary_before$total_variants_before
            )
        )
        message(
            sprintf(
                "Total variants after filtering:  %d",
                summary_after$total_variants_after
            )
        )
        message(
            sprintf(
                "\nMean variants before:            %.2f",
                summary_before$mean_variants_before
            )
        )
        message(
            sprintf(
                "Mean variants after:             %.2f",
                summary_after$mean_variants_after
            )
        )
        message(
            sprintf(
                "\nMedian variants before:          %.2f",
                summary_before$median_variants_before
            )
        )
        message(
            sprintf(
                "Median variants after:           %.2f",
                summary_after$median_variants_after
            )
        )
    }
    
    # Return the filtered maf
    return(maf_filter)
}

loadDFSPGroupInfo <- function() {
    
    sample_groups <- list(
        `U-DFSP` = c("Classic", "Myxoid", "Pigmented"),
        `Pre_FST` = c(
            "Pretransformed classic",
            "Pretransformed myxoid",
            "Paired classic",
            "Paired myxoid"
        ),
        `Post_FST` = c(
            "Posttransformed FST",
            "Paired FST",
            "Paired Pleomorphic"
        ),
        `FS_DFSP` = c("Unpaired FST")
    )

    sample_info <- LoadSampleInfo() |>
        filter(Specimen.Class == "Tumour") |>
        select(
            Sample.ID, Diagnosis, Specimen.Class, Specimen.Nature, Histology.Nature,
            Somatic.Status, purity, ploidy
        ) |>
        mutate(
            sample_group = case_when(
                Histology.Nature %in% sample_groups$`U-DFSP` ~ "U-DFSP",
                Histology.Nature %in% sample_groups$`Pre_FST` ~ "Pre-FST",
                Histology.Nature %in% sample_groups$`Post_FST` ~ "Post-FST",
                Histology.Nature %in% sample_groups$`FS_DFSP` ~ "FS-DFSP",
                TRUE ~ "Other"
            )
        ) |>
        rename(Tumor_Sample_Barcode = Sample.ID)
    
    sample_info
}

addDFSPGroupInfo <- function(maf_tbl) {

    ## Load the sample group info
    clinical_info <- LoadClinicalInfo()

    group_info <- clinical_info |> 
        dplyr::select(Sample.ID, FST.Group) |> 
        rename(
            Tumor_Sample_Barcode = Sample.ID
        )

    group_levels <- levels(group_info$FST.Group)

    ## Add the sample group to the maf data
    maf_tbl <- maf_tbl |>
        left_join(
            group_info,
            by = "Tumor_Sample_Barcode"
        )
    # maf_obj@clinical.data <- maf_obj@clinical.data |> 
    #     select(Tumor_Sample_Barcode) |>
    #     left_join(group_info, by = "Tumor_Sample_Barcode")

    # maf_obj
    maf_tbl
}

CollectPCGRCNVData <- function(
    dir = "data/wes/PCGR", 
    sig_gain_thres = NULL
) {

    ## Get all the PCGR files in the directory
    pcgr_files <- dir_ls(
        dir,
        glob = "*.cna_gene_ann.tsv.gz",
        recurse = TRUE
    )

    cna_list <- list()

    ## Loop through the files and read the data
    for (file in pcgr_files) {

        # file <- pcgr_files[[1]]
        message(" - Processing: ", file)
        
        if (!file.exists(file)) {
            warning("File does not exist: ", file)
            next
        }

        cna_data <- read_tsv(file, show_col_types = FALSE) |>
            clean_names() |> 
            mutate(
                total_cn = cn_major + cn_minor,
                .after = cn_minor
            )
        
        cna_data |> 
            select(sample_id, var_id, cn_major, cn_minor, variant_class) |>
            distinct()

        # if (!is.null(sig_gain_thres)) {
        #     cna_data <- cna_data |>
        #         filter(total_cnv >=4 | total_cnv == 0)
        # }   

        cna_list[[file]] <- cna_data

    }
    
    cna_tbl <- list_rbind(cna_list)

    ## Return
    cna_tbl
}

CollectPCGRSNVINDELData <- function(
    dir = "data/wes/PCGR", 
    sheet_name="SOMATIC_SNV_INDEL"
) {
    
    ## The sheet names in the PCGR files are:
    # SAMPLE_ASSAY
    # SOMATIC_SNV_INDEL
    # SOMATIC_SNV_INDEL_BIOMARKER
    # SOMATIC_CNA
    # TMB
    # MSI
    # MUTATIONAL_SIGNATURE

    ## Get all the PCGR files in the directory
    pcgr_files <- dir_ls(
        dir, 
        glob = "*.pcgr.grch38.xlsx", 
        recurse = TRUE
    )
    
    ## Get the sample names from the file names
    # pcgr_samples <- pcgr_files |> 
    #     path_file() |> 
    #     str_remove(".pcgr.grch38.xlsx")
    
    # message("Collecting PCGR data from ", length(pcgr_files), " files")

    # # Load all data from PCGR files
    # pcgr_data <- tibble(
    #     sample = pcgr_samples,
    #     file = pcgr_files
    # ) |>
    #     mutate(
    #         data = map2(
    #             file,
    #             sample,
    #             \(f, s) {
    #                 message(" - Processing sample = ", s)
    #                 read_xlsx(f, sheet = sheet_name)
    #             },
    #             .progress = TRUE
    #         )
    #     ) |> 
    #     unnest(data) |> 
    #     select(-sample, -file) |> 
    #     janitor::clean_names()
    
    data_list <- list()
    ## Loop through the files and read the data
    for (file in pcgr_files) {

        message(" - Processing: ", file)
        
        if (!file.exists(file)) {
            warning("File does not exist: ", file)
            next
        }

        data_list[[file]] <- read_xlsx(file, sheet = sheet_name) |>
            janitor::clean_names()

    }
    
    data <- list_rbind(data_list)

    ## Return
    data
}

FilterPCGRSNVIndels <- function(
    pcgr_data,
    min_dp_tumor = 20,
    min_vaf_tumor = 0.05,
    min_dp_control = 10,
    max_vaf_control = 0.01,
    max_population_frequency = 0.001
) {
    ## Check the columns in the PCGR data
    # all_columns <- colnames(pcgr_data)
    # print(all_columns)
    # anno_info <- c(
    #     "variant_class", "consequence", "loss_of_function",
    #     "exonic_status", "mutation_hotspot", "splice_effect",
    #     "actionability", "clinvar_classification"
    # )

    # str(pcgr_data)
    # unique(pcgr_data$variant_class)
    # unique(pcgr_data$consequence)
    # unique(pcgr_data$loss_of_function)
    # unique(pcgr_data$exonic_status)
    # unique(pcgr_data$mutation_hotspot)
    # unique(pcgr_data$splice_effect)
    # unique(pcgr_data$actionability)
    # unique(pcgr_data$clinvar_classification)
    # unique(pcgr_data$actionability)
    # unique(pcgr_data$effect)

    ## Sequencing depth and variant allele frequency filters
    # str(pcgr_data)
    message(
        "Number of samples in PCGR data = ", length(unique(pcgr_data$sample_id))
    )

    message(
        "Number of variants in PCGR data before filtering: ", nrow(pcgr_data)
    )
    
    oncogenic_variants <- pcgr_data |>
        filter(
            oncogenicity %in% c("Oncogenic", "Likely_Oncogenic")
        )
    
    message(
        "Number of oncogenic variants in PCGR data: ", nrow(oncogenic_variants)
    )
    
    actionable_variants <- pcgr_data |>
        filter(
            actionability %in% c("Potential significance")
        )
    message(
        "Number of actionable variants in PCGR data: ", nrow(actionable_variants)
    )
    
    tert_variants <- pcgr_data |>
        filter(symbol == "TERT")

    tp53_variants <- pcgr_data |>
        filter(symbol == "TP53")
    
    colnames(pcgr_data)
    
    unique(pcgr_data$consequence)
    unique(pcgr_data$loss_of_function)
    unique(pcgr_data$coding_status)
    unique(pcgr_data$exonic_status)

    ## OncoKB gene list
    oncokb_genes <- LoadOncoKBGeneList()
    oncokb_genes <- pcgr_data |> 
        filter(symbol %in% oncokb_genes) |>
        pull(symbol) |> 
        unique()

    message(
        "Number of OncoKB genes in PCGR data: ", length(oncokb_genes)
    )

    filtered_data <- pcgr_data

    ## Population frequency filtering
    if (!is.null(max_population_frequency)) {
        
        filtered_data <- filtered_data |>
            filter(
                is.na(gnom_a_de_af) | gnom_a_de_af < max_population_frequency
            )
    }

    ## Tumor filtering
    if (!is.null(min_dp_tumor)) {
        
        filtered_data <- filtered_data |>
            filter(dp_tumor >= min_dp_tumor)
    }
    
    if (!is.null(min_vaf_tumor)) {
        
        filtered_data <- filtered_data |>
            filter(vaf_tumor >= min_vaf_tumor)
    }

    ## Control filtering
    if (!is.null(min_dp_control)) {
        
        filtered_data <- filtered_data |>
            filter(dp_control >= min_dp_control)
    }

    if (!is.null(max_vaf_control)) {
        
        filtered_data <- filtered_data |>
            filter(
                is.na(vaf_control) | vaf_control <= max_vaf_control
            )
    }

    # filtered_data <- filtered_data |>
    #     filter(dp_tumor >= min_dp_tumor & vaf_tumor >= min_vaf_tumor) |>
    #     filter(dp_control >= min_dp_control & vaf_control <= max_vaf_control)

    message(
        "\nNumber of variants after Tumour/Control DP and VAF filtering = ", 
        nrow(filtered_data)
    )

    ## "====================================================================="
    ## Indel variants
    ## "====================================================================="
    filtered_indels <- filtered_data |>
        filter(exonic_status %in% c("exonic", "splicing")) |>
        filter(variant_class %in% c("insertion", "deletion"))

    message("Number of InDels variants = ", nrow(filtered_indels))

    ## Indel actionability variants
    indel_oncogenicity <- filtered_indels |>
        filter(oncogenicity %in% c("Oncogenic", "Likely_Oncogenic"))

    message(
        "   - Number of InDels with oncogenicity = ", nrow(indel_oncogenicity)
    )

    ## Indel actionability variants
    indel_actionability <- filtered_indels |>
        filter(actionability %in% c("Potential significance"))

    message(
        "   - Number of InDels with potential actionability = ", nrow(indel_actionability)
    )

    ## "====================================================================="
    ## SNV variants
    ## "====================================================================="
    filtered_snvs <- filtered_data |>
        filter(!(variant_class %in% c("insertion", "deletion"))) |> 
        filter(!consequence %in% c("synonymous_variant"))

    message("\nNumber of Non-indel variants (SNV, DNV, TNV, etc) = ", nrow(filtered_snvs))

    ## Actinability snv variants
    unique(filtered_data$actionability)
    snvs_actionability <- filtered_snvs |>
        filter(actionability %in% c("Potential significance"))

    message(
        "   - Number of SNVs with potential actionability = ", nrow(snvs_actionability)
    )

    ## Oncogenicity snv variants
    unique(filtered_snvs$oncogenicity)
    snvs_oncogenicity <- filtered_snvs |>
        filter(oncogenicity %in% c("Oncogenic", "Likely_Oncogenic"))

    message(
        "   - Number of SNVs with oncogenicity = ", nrow(snvs_oncogenicity)
    )

    ## ClinVar classification snv variants
    unique(filtered_snvs$clinvar_classification)
    snvs_clinvar <- filtered_snvs |>
        filter(clinvar_classification %in% c("Pathogenic", "Likely_pathogenic"))

    message(
        "   - Number of SNVs with ClinVar classification = ", nrow(snvs_clinvar)
    )

    ## Loss of function snv variants
    unique(filtered_snvs$loss_of_function)
    snvs_loss_of_function <- filtered_snvs |>
        filter(loss_of_function)

    message(
        "   - Number of SNVs with loss of function = ", nrow(snvs_loss_of_function)
    )

    ## Splice effect snv variants
    unique(filtered_snvs$splice_effect)
    snvs_splice_disrupting <- filtered_snvs |>
        filter(grepl("disrupting", splice_effect, ignore.case = TRUE))

    message(
        "   - Number of SNVs with splice disrupting effect = ", nrow(snvs_splice_disrupting)
    )

    ## Insilico predictions variant effect on protein function, damaging or pathogenic
    unique(filtered_snvs$effect_predictions)[20]
    predictions <- unique(filtered_snvs$effect_predictions)[20]
    ## Filter for variants predicted as damaging (D)
    snvs_damaging <- filtered_snvs |>
        filter(grepl("D", effect_predictions))

    message(
        "   - Number of SNVs with damaging effect predictions = ", nrow(snvs_damaging)
    )

    ## Combine all variants
    final_variants <- bind_rows(
        oncogenic_variants,
        actionable_variants,
        tert_variants,
        tp53_variants,
        filtered_indels,
        indel_actionability,
        indel_oncogenicity,
        snvs_actionability,
        snvs_oncogenicity,
        snvs_clinvar,
        snvs_loss_of_function,
        snvs_splice_disrupting,
        snvs_damaging
    ) |>
        distinct()
    
    message(
        paste0(
            "\n------------------------------------------------------",
            "\nAfter filtering and combining:\n"
        )
    )

    message(
        "   Number of included unique variants = ", nrow(final_variants)
    )

    sample_variants <- final_variants |> 
        group_by(sample_id) |> 
        summarise(
            n_variants = n(),
            .groups = "drop"
        ) |> 
        arrange(n_variants)
    

    message(
        "   Number of samples with variants: ", nrow(sample_variants)
    )

    summary_variants <- sample_variants |> 
        summarise(
            min_variants = min(n_variants, na.rm = TRUE),
            max_variants = max(n_variants, na.rm = TRUE),
            mean_variants = mean(n_variants, na.rm = TRUE),
            median_variants = median(n_variants, na.rm = TRUE)
        )
    # message(
    #     sprintf(
    #         "   Summary of variants per sample:\n%s",
    #         paste(capture.output(print(summary_variants)), collapse = "\n")
    #     )
    # )

    ## Return
    final_variants
}

ConvertPCGRToMaftools <- function(pcgr_tbl) {

    maf_tbl <- pcgr_tbl |>
        # Extract chromosome and position from genomic_change
        mutate(
            # Parse genomic change: e.g., "1:g.222633099TCCAG...>T"
            chr_pos = str_extract(genomic_change, "^(\\d+|X|Y):g\\.(\\d+)"),
            chromosome = str_extract(chr_pos, "^(\\d+|X|Y)"),
            start_pos = as.numeric(str_extract(chr_pos, "(\\d+)$")),
            
            # Extract reference and alternate alleles
            reference_allele = str_extract(genomic_change, "([ATCG]+)>", group = 1),
            alternate_allele = str_extract(genomic_change, ">([ATCG]+)$", group = 1)
        ) |> 
        ## Maftools requires the following columns
        mutate(Hugo_Symbol = symbol) |> 
        mutate(Chromosome = paste0("chr", chromosome)) |> 
        mutate(Start_Position = start_pos) |> 
        mutate(
            End_Position = case_when(
                variant_class == "deletion" ~ start_pos + nchar(reference_allele) - 1,
                variant_class == "insertion" ~ start_pos,
                variant_class == "SNV" ~ start_pos,
                variant_class == "substitution" ~ start_pos + nchar(reference_allele) - 1,
                TRUE ~ start_pos
            )
        ) |> 
        mutate(
            Reference_Allele = reference_allele,
            Tumor_Seq_Allele2 = alternate_allele
        ) |> 
        # mutate(
        #     Reference_Allele = case_when(
        #         variant_class == "insertion" ~ "-",
        #         is.na(reference_allele) ~ "-",
        #         TRUE ~ reference_allele
        #     )
            
        # ) |> 
        # mutate(
        #     Tumor_Seq_Allele2 = case_when(
        #         variant_class == "deletion" ~ "-",
        #         is.na(alternate_allele) ~ "-",
        #         TRUE ~ alternate_allele
        #     )
        # ) |> 
        mutate(
            Variant_Classification = case_when(
                str_detect(consequence, "stop_gained") ~ "Nonsense_Mutation",
                str_detect(consequence, "missense") ~ "Missense_Mutation",
                str_detect(consequence, "synonymous") ~ "Silent",
                str_detect(consequence, "frameshift") ~ ifelse(variant_class == "insertion", 
                                                            "Frame_Shift_Ins", "Frame_Shift_Del"),
                str_detect(consequence, "inframe_insertion") ~ "In_Frame_Ins",
                str_detect(consequence, "inframe_deletion") ~ "In_Frame_Del",
                str_detect(consequence, "splice") ~ "Splice_Site",
                str_detect(consequence, "start_lost") ~ "Translation_Start_Site",
                str_detect(consequence, "stop_lost") ~ "Nonstop_Mutation",
                str_detect(consequence, "coding_sequence_variant") ~ "Missense_Mutation",
                # Handle multiple consequences (take the first one)
                str_detect(consequence, "protein_altering") ~ "Missense_Mutation",
                TRUE ~ "Unknown"
            )
        ) |> 
        mutate(
            Variant_Type = case_when(
                variant_class == "SNV" ~ "SNP",
                variant_class == "insertion" ~ "INS", 
                variant_class == "deletion" ~ "DEL",
                variant_class == "substitution" & nchar(reference_allele) == 2 ~ "DNP",
                variant_class == "substitution" & nchar(reference_allele) > 2 ~ "ONP",
                TRUE ~ NA_character_
            ),
        ) |> 
        mutate(Tumor_Sample_Barcode = sample_id) |> 
        dplyr::select(
            -chr_pos, -chromosome, -start_pos,
            -reference_allele, -alternate_allele
        )

    maf_required_cols <- c(
        "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", 
        "Variant_Type"
    )
    
    ## Reorder the columns to have the required ones first
    maf_tbl |> 
        relocate(
            all_of(maf_required_cols), 
            everything()
        )

    maf_tbl

}

LoadClinicalInfo <- function(is_somatic_matched = TRUE) {

    wes_file <- "data/clinical/DFSP_WES_clinical.csv"

    clinical_info <- read_csv(
        here(wes_file), 
        show_col_types = FALSE
    ) |> 
        mutate(
            Tumor_Sample_Barcode = Sample.ID
        )
    
    ## Combine the Epic Meth group data
    meth_file <- here("data/clinical/phenoData_DFSP_allTumour_CNV_Meth.csv")
    meth_data <- read_csv(
        file = meth_file, 
        show_col_types = FALSE,
        name_repair = "unique_quiet"
    ) |> 
        select(Sample.ID, Meth.Subtype.Main.2Class) |> 
        filter(!is.na(Meth.Subtype.Main.2Class))

    clinical_info <- clinical_info |> 
        left_join(meth_data, by = "Sample.ID")

    # unique(clinical_info$FST.Group)
    # clinical_info <- clinical_info |> 
    #     mutate(
    #         FST.Group =factor(
    #             FST.Group, 
    #             levels = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
    #             # ordered = TRUE
    #         )
    #     )

    if (is_somatic_matched) {
        clinical_info <- clinical_info |> 
            filter(Somatic.Status == "Matched")
    }

    clinical_info
}

SaveData <- function(obj, dir, filename) {
    
    ## Save data to qs file

    fs::dir_create(here(dir))

    qsave(obj, here(dir, paste0(filename, ".qs")))
}

LoadData <- function(dir, filename) {
    
    ## Load qs data
    qread(here(dir, paste0(filename, ".qs")))
    
}

compareDFSPGroups <- function(
    data,
    group_col = "FST.Group",
    value_col = "total_perMB",
    stat_method = "wilcox_test",
    is_paired = FALSE,
    p_adjust_method = "BH",
    p_value_threshold = 0.05

) {

    ## Check if the group_col and value_col are in the data
    if (!(group_col %in% colnames(data))) {
        stop(sprintf("Column '%s' not found in data", group_col))
    }
    
    if (!(value_col %in% colnames(data))) {
        stop(sprintf("Column '%s' not found in data", value_col))
    }

    ## Get the summary statistics for each group
    summary_stats <- data |>
        group_by(!!sym(group_col)) |>
        get_summary_stats(!!sym(value_col))
    
    # print(summary_stats)

    # ## TMB is the count data, so we will use a non-parametric test
    # normality_test <- tmb_tbl |> 
    #     group_by(FST.Group) |>
    #     shapiro_test(total_perMB) |> 
    #     mutate(
    #         p = format(p, scientific = FALSE),
    #         is_normal = case_when(
    #             p < 0.05 ~ "No",
    #             TRUE ~ "Yes"
    #         )
    #     ) 
    # print(normality_test)

    ## statistical test
    if (stat_method == "wilcox_test") {

        stat_test <- data |> 
            mutate(value = !!sym(value_col), group = !!sym(group_col)) |>
            wilcox_test(
                value ~ group,
                p.adjust.method = p_adjust_method,
                paired = is_paired
            ) |> 
            add_significance(p.col = "p") |> 
            add_xy_position(
                x = group_col, 
                step.increase = 0.12,
                dodge = 0.8
            )

    } else {
        
        stat_test <- data |> 
            mutate(value = !!sym(value_col), group = !!sym(group_col)) |>
            t_test(
                value ~ group,
                p.adjust.method = p_adjust_method,
                paired = is_paired
            ) |> 
            add_significance(p.col = "p") |> 
            add_xy_position(
                x = group_col, 
                step.increase = 0.12,
                dodge = 0.8
            ) 
    }
    
    ## Boxplot groups
    plot <- data |> 
        ggplot(aes(x = !!sym(group_col), y = !!sym(value_col))) +
        geom_boxplot(outliers = FALSE, width = 0.5) +
        geom_jitter(aes(color = !!sym(group_col)), width = 0.2, size = 0.5) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
        stat_pvalue_manual(
            stat_test |> filter(p < p_value_threshold),
            label = "p = {scales::pvalue(p)}",
            y.position = "y.position",
            size = 0.1,
            label.size = 3,
            bracket.size = 0.2,
            tip.length = 0.01,
            bracket.nudge.y = 0.1,
            hide.ns = FALSE
        )+
        labs(
            title = "",
            x = "",
            y = "TMB/mb"
        ) +
        figure_theme2()

    ## Return
    list(
        summary_stats = summary_stats,
        stat_test = stat_test,
        plot = plot
    )
}

loadDFSPColorConfigs <- function() {
    list(
        FST.Group = c(
            "U-DFSP"    = "#3498db", # Blue - for untransformed DFSP
            "Pre-FST"   = "#2ecc71", # Green - for pre-transformation samples
            "Post-FST"  = "#e74c3c", # Red - for post-transformation samples
            "FS-DFSP"   = "#9b59b6" # Purple - for unpaired FST
        ),
        Specimen.Nature = c(
            "Primary"       = "#ff9800", # Orange
            "Recurrence"    = "#009688", # Teal
            "Metastasis"    = "#795548", # Brown
            "Residual"      = "#607d8b" # Blue-gray
        )
    )
}


LoadDFSPSampleGroups <- function() {

    file <- here("data/clinical/DFSP_WES_sample_groups.xlsx")
    
    groups <- excel_sheets(file)

    sample_groups <- map(
        groups,
        \(sheet) {
            read_excel(file, sheet = sheet, col_names = TRUE) |>
                pull(Sample.ID)
        }
    ) |> 
        set_names(groups)

    # return
    sample_groups
}

LoadSampleGroupInfo <- function(is_somatic_matched = FALSE) {

    ## Load the clinical info
    clinical_info <- LoadClinicalInfo()

    if (is_somatic_matched) {
        
        ## Filter out the samples that are not somatic matched
        clinical_info <- clinical_info |> 
            filter(Somatic.Status == "Matched")
    }

    ## Generate the WES sample groups
    ## Patient levele analysi:
    ## Main: Choosen the sample with the highest grade: presence of FST
    ## Main.Met

    ## "-----------------------------------------------------------------"
    ## All WES samples, EPIC=150, WES=161
    ## "-----------------------------------------------------------------"
    all_tumors <- clinical_info |> pull(Sample.ID)

    num_all_tumors <- length(all_tumors)

    ## "-----------------------------------------------------------------"
    ## Overall FST group
    ## "-----------------------------------------------------------------"
    FST <- clinical_info |> filter(FST == "Yes") |> pull(Sample.ID)

    `Non-FST` <- clinical_info |> filter(FST == "No") |> pull(Sample.ID)

    num_FST <- length(FST)
    num_Non_FST <- length(`Non-FST`)

    ## "-----------------------------------------------------------------"
    ## Overall Metastasis group
    ## "-----------------------------------------------------------------"
    ## The column "Metastasis", it means whether the patient has metastasis, 
    ## but the same patient can have samples coming from 
    ## primary site or metastatic site
    ## Case level metastasis
    Meta <- clinical_info |> filter(Metastasis == "Yes") |> pull(Sample.ID)

    `Non-Meta` <- clinical_info |> filter(Metastasis == "No") |> pull(Sample.ID)

    num_Meta <- length(Meta)
    num_Non_Meta <- length(`Non-Meta`)

    ## "-----------------------------------------------------------------"
    ## FST.subgroup
    ## "-----------------------------------------------------------------"

    ## U-DFSP, EPIC=48, WES=52
    `U-DFSP` = clinical_info |> 
        filter(FST.Group %in% c("U-DFSP")) |> 
        pull(Sample.ID)
    
    num_U_DFSP <- length(`U-DFSP`)

    ## Pre-FST, EPIC=40, WES=45
    `Pre-FST` = clinical_info |> 
        filter(FST.Group %in% c("Pre-FST")) |> 
        pull(Sample.ID)

    num_Pre_FST <- length(`Pre-FST`)

    ## Post-FST, EPIC=39, WES=41
    `Post-FST` = clinical_info |> 
        filter(FST.Group %in% c("Post-FST")) |> 
        pull(Sample.ID)
    
    num_Post_FST <- length(`Post-FST`)

    ## FS-DFSP group, EPIC=23, WES=23
    `FS-DFSP` = clinical_info |> 
        filter(FST.Group %in% c("FS-DFSP")) |> 
        pull(Sample.ID)
    
    num_FS_DFSP <- length(`FS-DFSP`)

    ## "-----------------------------------------------------------------"
    ## FS-DFSP subgroup
    ## "-----------------------------------------------------------------"
    ## FS-DFSP Meta, EPIC=6, WES=11
    `FS-DFSP_Meta` = clinical_info |> 
        filter(FST.Group %in% c("FS-DFSP")) |> 
        filter(Metastasis == "Yes") |> 
        pull(Sample.ID)
    
    num_FS_DFSP_Meta <- length(`FS-DFSP_Meta`)

    ## FS-DFSP Primary Rrecurrence, EPIC=17, WES=17
    `FS-DFSP_Pri_Rec` = clinical_info |> 
        filter(FST.Group %in% c("FS-DFSP")) |>
        filter(Specimen.Nature %in% c("Primary", "Recurrence")) |> 
        pull(Sample.ID)
    
    num_FS_DFSP_Pri_Rec <- length(`FS-DFSP_Pri_Rec`)

    ## "-----------------------------------------------------------------"
    ## Comparsion groups
    ## "-----------------------------------------------------------------"
    `Non-Meta_vs_Meta` <- c(`Non-Meta`, Meta)
    `Non-FST_vs_FST` <- c(`Non-FST`, FST)
    `U-DFSP_vs_Pre-FST` <- c(`U-DFSP`, `Pre-FST`)
    `U-DFSP_vs_Post-FST` <- c(`U-DFSP`, `Post-FST`)
    `U-DFSP_vs_FS-DFSP` <- c(`U-DFSP`, `FS-DFSP`)
    `Pre-FST_vs_Post-FST` <- c(`Pre-FST`, `Post-FST`)
    `Pre-FST_vs_FS-DFSP` <- c(`Pre-FST`, `FS-DFSP`)
    `Post-FST_vs_FS-DFSP` <- c(`Post-FST`, `FS-DFSP`)
    `U-DFSP+Pre-FST_vs_Post-FST+FS-DFSP` <- c(
        `U-DFSP`, `Pre-FST`, `Post-FST`, `FS-DFSP`
    )
    `FS-DFSP_Pri_Rec_vs_FS-DFSP_Meta` <- c(`FS-DFSP_Pri_Rec`, `FS-DFSP_Meta`)

    # return
    sample_group_list <- list(
        ## single groups
        groups = list(
            "all_tumors" = all_tumors,
            "FST" = FST,
            "Non-FST" = `Non-FST`,
            "Meta" = Meta,
            "Non-Meta" = `Non-Meta`,
            "U-DFSP" = `U-DFSP`,
            "Pre-FST" = `Pre-FST`,
            "Post-FST" = `Post-FST`,
            "FS-DFSP" = `FS-DFSP`,
            "FS-DFSP_Pri_Rec" = `FS-DFSP_Pri_Rec`,
            "FS-DFSP_Meta" = `FS-DFSP_Meta`
        )
    )

    ## Add the comparison groups to the list
    sample_group_list[["comparisons"]][[
        paste0("U-DFSP(", num_U_DFSP, ")_vs_Pre-FST(", num_Pre_FST, ")")
    ]] <- `U-DFSP_vs_Pre-FST`
    
    sample_group_list[["comparisons"]][[
        paste0("U-DFSP(", num_U_DFSP, ")_vs_Post-FST(", num_Post_FST, ")")
    ]] <- `U-DFSP_vs_Post-FST`
    
    sample_group_list[["comparisons"]][[
        paste0("U-DFSP(", num_U_DFSP, ")_vs_FS-DFSP(", num_FS_DFSP, ")")
    ]] <- `U-DFSP_vs_FS-DFSP`

    sample_group_list[["comparisons"]][[
        paste0("Pre-FST(", num_Pre_FST, ")_vs_Post-FST(", num_Post_FST, ")")
    ]] <- `Pre-FST_vs_Post-FST`
    
    sample_group_list[["comparisons"]][[
        paste0("Pre-FST(", num_Pre_FST, ")_vs_FS-DFSP(", num_FS_DFSP, ")")
    ]] <- `Pre-FST_vs_FS-DFSP`

    sample_group_list[["comparisons"]][[
        paste0("Post-FST(", num_Post_FST, ")_vs_FS-DFSP(", num_FS_DFSP, ")")
    ]] <- `Post-FST_vs_FS-DFSP`

    sample_group_list[["comparisons"]][[
        paste0(
            "FS-DFSP_Pri_Rec(", num_FS_DFSP_Pri_Rec, 
            ")_vs_FS-DFSP_Meta(", num_FS_DFSP_Meta, ")"
        )
    ]] <- `FS-DFSP_Pri_Rec_vs_FS-DFSP_Meta`

    sample_group_list[["comparisons"]][[
        paste0("Non-Meta(", num_Non_Meta, ")_vs_Meta(", num_Meta, ")")
    ]] <- `Non-Meta_vs_Meta`

    sample_group_list[["comparisons"]][[
        paste0("Non-FST(", num_Non_FST, ")_vs_FST(", num_FST, ")")
    ]] <- `Non-FST_vs_FST`

    ## return
    sample_group_list
}

loadMsigdbGeneSet <- function() {
    ## Get all the potential gene sets
    msigdb_gs <- bind_rows(
        ### Hallmark gene set
        msigdbr(
            species = "Homo sapiens", 
            collection = "H"
        ),

        ### Curated gene set, canonical pathways
        msigdbr(
            species = "Homo sapiens", 
            collection = "C2", 
            subcollection = "CP:KEGG_LEGACY"
        ),
        msigdbr(
            species = "Homo sapiens", 
            collection = "C2", 
            subcollection = "CP:REACTOME"
        ),

        ### Go gene sets
        msigdbr(
            species = "Homo sapiens", 
            # subcollection = "GO:BP",
            collection = "C5"
        )
    )
    
    msigdb_gs <- msigdb_gs |>
        dplyr::select(gs_name, gene_symbol, ncbi_gene, ensembl_gene) |> 
        mutate(db_name = str_extract(gs_name, "^[^_]+")) |>
        mutate(gs_name = str_extract(gs_name, "(?<=_).+"))
    
    msigdb_gs
}

cleanGeneSetName <- function(msigdb_gs) {
        
        # msigdb_gs |> 
        #     filter(grepl("mhc_class", gs_name, ignore.case = TRUE))

        msigdb_gs |> 
            mutate(
                gs_name = stringr::str_to_sentence(gs_name),
                gs_name = gsub("Atp", "ATP", gs_name),
                gs_name = gsub("atp", "ATP", gs_name),
                gs_name = gsub("Pd_1", "PD-1", gs_name),
                gs_name = gsub("orc_complex", "ORC_complex", gs_name),
                gs_name = gsub("Sirt1", "SIRT1", gs_name),
                gs_name = gsub("rrna_expression", "rRNA_expression", gs_name),
                gs_name = gsub("mrna", "mRNA", gs_name),
                gs_name = gsub("cap_binding_complex", "cap-binding_complex", gs_name),
                gs_name = gsub("and_eifs", "and_eIFs", gs_name),
                gs_name = gsub("Nonhomologous", "Non-homologous", gs_name),
                gs_name = gsub("nhej", "NHEJ", gs_name),
                gs_name = gsub("tca", "TCA", gs_name),
                gs_name = gsub("Mtorc1", "mTORC1", gs_name),
                gs_name = gsub("ar_androgen_receptor", "AR_androgen_receptor", gs_name),
                gs_name = gsub("klk2_and_klk3", "KLK2_and_KLK3", gs_name),
                gs_name = gsub("ap_site_formation", "AP_site_formation", gs_name),
                gs_name = gsub("Activated_pkn1", "Activated_PKN1", gs_name),
                gs_name = gsub("Hdms_demethylate", "HDMs_demethylate", gs_name),
                gs_name = gsub("Er_to_golgi", "ER_to_golgi", gs_name),
                gs_name = gsub("Gpcr", "GPCR", gs_name),
                gs_name = gsub("mhc_protein", "MHC_protein ", gs_name),
                gs_name = gsub("Mhc_class", "MHC_class", gs_name),
                gs_name = gsub("mhc_class", "MHC_class", gs_name),
                gs_name = gsub("Mhc_class_ii", "MHC_class_II", gs_name),
                gs_name = gsub("Mhc_protein", "MHC_protein", gs_name),
                gs_name = gsub("Apc_cdc20", "APC(Cdc20)", gs_name),
                gs_name = gsub("nek2a", "Nek2A", gs_name),
                gs_name = gsub("Il6", "IL6", gs_name),
                gs_name = gsub("jak", "JAK", gs_name),
                gs_name = gsub("jak_stat3", "JAK/STAT3", gs_name),
                gs_name = gsub("stat3", "STAT3", gs_name),
                gs_name = gsub("perk", "PERK", gs_name),
                gs_name = gsub("trna", "tRNA", gs_name),
                gs_name = gsub("Tnfa", "TNFα", gs_name),
                gs_name = gsub("Tnfa", "TNFα", gs_name),
                gs_name = gsub("nfkb", "NFkB", gs_name),
                gs_name = gsub("E2f", "E2F", gs_name),
                gs_name = gsub("e1", "E1", gs_name),
                gs_name = gsub("e2", "E2", gs_name),
                gs_name = gsub("Il2_stat5", "IL2_STAT5", gs_name),
                gs_name = gsub("G2m", "G2/M", gs_name),
                gs_name = gsub("eif2ak1_hri", "EIF2AK1 (HRI)", gs_name),
                gs_name = gsub("sars_cov_1", "SARS-CoV-1", gs_name),
                gs_name = gsub("sahf", "SAHF", gs_name),
                gs_name = gsub("h3_k27", "H3K27", gs_name),
                gs_name = gsub("Dna", "DNA", gs_name),
                gs_name = gsub("dna", "DNA", gs_name),
                gs_name = gsub("G1_s", "G1/S", gs_name),
                gs_name = gsub("Kras", "KRAS", gs_name),
                gs_name = gsub("cd4", "CD4", gs_name),
                gs_name = gsub("t_cell", "T_cell", gs_name),
                gs_name = gsub("gtpase", "GTPase", gs_name),
                gs_name = gsub("Fcgr", "FcγR", gs_name),
                gs_name = gsub("fcgr", "FcγR", gs_name),
                gs_name = gsub("Tp53", "TP53", gs_name),
                gs_name = gsub("Chk1_chk2_cds1", "Chk1/Chk2(Cds1)", gs_name),
                gs_name = gsub("cyclin_b_cdk1", "Cyclin_B:Cdk1", gs_name),
                gs_name = gsub("Nadph_", "NADPH_", gs_name),
                gs_name = gsub("Nadp_", "NADP_", gs_name)
            )
}

CollectCNVFacets <- function(facet_dir) {

    facet_files <- dir_ls(
        facet_dir, 
        glob = "*.annotated.tsv", 
        recurse = TRUE
    )

    message("Collecting Facets CNV data from ", length(facet_files), " files")

    cnv_facet_data <- tibble(
        sample = path_file(facet_files) |> str_remove(".annotated.tsv"),
        file = facet_files
    ) |>
        mutate(
            data = map2(
                sample,
                file,
                \(s, f) {
                    message("Processing sample = ", s)
                    read_tsv(f, show_col_types = FALSE)
                },
                .progress = TRUE
            )
        ) |> 
        unnest(data) |> 
        dplyr::select(-sample, -file) |> 
        janitor::clean_names()

    cnv_facet_data |>
        dplyr::select(
            "gene", "sample_id", "svtype", "chromosome", "start", "start_position", 
            "end_position", "end"
        ) |>
        dplyr::rename(Gene = gene, Sample_name = sample_id, CN = svtype) |> 
        distinct()

}

getDFSPSampleGroups <- function(sample_ids) {

    ## Load clinical info
    clinical_info <- LoadClinicalInfo()

    ## Get the sample groups for the given sample IDs
    sample_idx <- match(sample_ids, clinical_info$Sample.ID)

    ## Check if the sample_idx is valid
    if (any(is.na(sample_idx))) {
        stop("Some sample IDs are not found in the clinical info.")
    }

    ## "-----------------------------------------------------------------"
    ## Metastaisis vs Non-Metastasis
    ## "-----------------------------------------------------------------"
    

}

LoadDFSPWESSamples <- function() {

    sample_info_dir <- "data/clinical"

    paired_tumor_samples <- read_tsv(
        here(sample_info_dir, "tumour_paired_samples.txt"),
        col_names = FALSE,
        show_col_types = FALSE
    )

    unpaired_tumor_samples <- read_tsv(
        here(sample_info_dir, "tumour_unpaired_samples.txt"),
        col_names = FALSE,
        show_col_types = FALSE
    )

    list(
        paired = paired_tumor_samples$X1,
        unpaired = unpaired_tumor_samples$X1
    )
}

## DEPRECATED: Old comparison functions replaced by compare_gistic_transformation()
## Use compare_gistic_transformation() for comprehensive GISTIC comparison analysis

## Updated function to compare GISTIC results between groups
compare_gistic_groups <- function(untransformed_gistic, transformed_gistic,
                                  output_prefix = "fst_transformation",
                                  untrans_total_samples = NULL,
                                  trans_total_samples = NULL) {

    library(dplyr)
    
    ## Extract amplification and deletion data
    gistic_summary_untrans <- getCytobandSummary(untransformed_gistic)
    gistic_summary_trans <- getCytobandSummary(transformed_gistic)
    
    ## Separate amplifications and deletions
    untrans_amp <- gistic_summary_untrans |> 
        filter(Variant_Classification == "Amp")
    untrans_del <- gistic_summary_untrans |> 
        filter(Variant_Classification == "Del")
    
    trans_amp <- gistic_summary_trans |> 
        filter(Variant_Classification == "Amp")
    trans_del <- gistic_summary_trans |> 
        filter(Variant_Classification == "Del")

    ## Get actual sample counts if not provided
    if (is.null(untrans_total_samples)) {
        untrans_total_samples <- length(getSampleSummary(untransformed_gistic)$Tumor_Sample_Barcode)
    }
    if (is.null(trans_total_samples)) {
        trans_total_samples <- length(getSampleSummary(transformed_gistic)$Tumor_Sample_Barcode)
    }

    ## Compare amplifications
    amp_comparison <- compare_cna_by_type(untrans_amp, trans_amp, "Amplification",
                                         untrans_total_samples, trans_total_samples)

    ## Compare deletions
    del_comparison <- compare_cna_by_type(untrans_del, trans_del, "Deletion",
                                         untrans_total_samples, trans_total_samples)

    ## Combine results
    transformation_cnas <- rbind(amp_comparison, del_comparison)

    ## Create results directory if it doesn't exist
    dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

    ## Save results
    write.table(transformation_cnas,
        here("results", paste0(output_prefix, "_cna_comparison.txt")),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    return(transformation_cnas)
}

## Updated helper function with correct sample totals
compare_cna_frequencies_updated <- function(untrans_data, trans_data, alteration_type,
                                           untrans_total_samples, trans_total_samples) {
    ## Get all unique cytobands
    all_cytobands <- unique(c(untrans_data$Cytoband, trans_data$Cytoband))

    comparison_results <- data.frame()

    for (cytoband in all_cytobands) {
        ## Get sample counts for this cytoband
        untrans_count <- ifelse(cytoband %in% untrans_data$Cytoband,
                               untrans_data$nSamples[untrans_data$Cytoband == cytoband], 0
        )
        trans_count <- ifelse(cytoband %in% trans_data$Cytoband,
                             trans_data$nSamples[trans_data$Cytoband == cytoband], 0
        )

        ## Calculate actual frequencies
        untrans_freq <- untrans_count / untrans_total_samples
        trans_freq <- trans_count / trans_total_samples
        freq_diff <- trans_freq - untrans_freq

        ## Get q-values
        untrans_qval <- ifelse(cytoband %in% untrans_data$Cytoband,
                              untrans_data$qvalues[untrans_data$Cytoband == cytoband], 1
        )
        trans_qval <- ifelse(cytoband %in% trans_data$Cytoband,
                            trans_data$qvalues[trans_data$Cytoband == cytoband], 1
        )

        ## Get gene counts
        untrans_genes <- ifelse(cytoband %in% untrans_data$Cytoband,
                               untrans_data$nGenes[untrans_data$Cytoband == cytoband], 0
        )
        trans_genes <- ifelse(cytoband %in% trans_data$Cytoband,
                             trans_data$nGenes[trans_data$Cytoband == cytoband], 0
        )

        ## Get wide peak limits
        untrans_peak <- ifelse(cytoband %in% untrans_data$Cytoband,
                              as.character(untrans_data$Wide_Peak_Limits[untrans_data$Cytoband == cytoband]), ""
        )
        trans_peak <- ifelse(cytoband %in% trans_data$Cytoband,
                            as.character(trans_data$Wide_Peak_Limits[trans_data$Cytoband == cytoband]), ""

        )

        comparison_results <- rbind(comparison_results, data.frame(
            cytoband = cytoband,
            alteration_type = alteration_type,
            untransformed_count = untrans_count,
            transformed_count = trans_count,
            untransformed_total = untrans_total_samples,
            transformed_total = trans_total_samples,
            untransformed_freq = round(untrans_freq, 3),
            transformed_freq = round(trans_freq, 3),
            frequency_difference = round(freq_diff, 3),
            untransformed_qval = untrans_qval,
            transformed_qval = trans_qval,
            untransformed_genes = untrans_genes,
            transformed_genes = trans_genes,
            untransformed_peak = untrans_peak,
            transformed_peak = trans_peak,
            transformation_enriched = freq_diff > 0.2, # >20% difference
            significant_in_transformed = trans_qval < 0.25,
            significant_in_untransformed = untrans_qval < 0.25
        ))
    }

    return(comparison_results)
}

## Function to get gene-level differences (updated)
compare_gistic_genes <- function(untransformed_gistic, transformed_gistic) {
    ## Initialize results
    amp_specific <- character(0)
    del_specific <- character(0)
    amp_enriched <- data.frame()
    del_enriched <- data.frame()

    ## Try to extract gene data from GISTIC objects using getGeneSummary()
    tryCatch({
        ## Get gene-level data using getGeneSummary function
        untrans_gene_summary <- getGeneSummary(untransformed_gistic)
        trans_gene_summary <- getGeneSummary(transformed_gistic)
        
        ## Separate amplifications and deletions
        untrans_amp_genes <- untrans_gene_summary |> filter(Variant_Classification == "Amp")
        untrans_del_genes <- untrans_gene_summary |> filter(Variant_Classification == "Del")
        trans_amp_genes <- trans_gene_summary |> filter(Variant_Classification == "Amp")
        trans_del_genes <- trans_gene_summary |> filter(Variant_Classification == "Del")

        ## Compare amplified genes
        if (nrow(untrans_amp_genes) > 0 && nrow(trans_amp_genes) > 0) {
            amp_specific <- setdiff(trans_amp_genes$Hugo_Symbol, untrans_amp_genes$Hugo_Symbol)
            amp_enriched <- find_enriched_genes(untrans_amp_genes, trans_amp_genes)
        } else if (nrow(trans_amp_genes) > 0) {
            amp_specific <- trans_amp_genes$Hugo_Symbol
        }

        ## Compare deleted genes
        if (nrow(untrans_del_genes) > 0 && nrow(trans_del_genes) > 0) {
            del_specific <- setdiff(trans_del_genes$Hugo_Symbol, untrans_del_genes$Hugo_Symbol)
            del_enriched <- find_enriched_genes(untrans_del_genes, trans_del_genes)
        } else if (nrow(trans_del_genes) > 0) {
            del_specific <- trans_del_genes$Hugo_Symbol
        }
        
    }, error = function(e) {
        message("Could not extract gene-level data from GISTIC objects: ", e$message)
        message("Skipping gene-level comparison")
    })

    ## Find transformation-specific genes
    transformation_specific_genes <- list(
        amp_specific = amp_specific,
        del_specific = del_specific,
        amp_enriched = amp_enriched,
        del_enriched = del_enriched
    )

    return(transformation_specific_genes)
}

## Updated helper function to find enriched genes
find_enriched_genes <- function(untrans_genes, trans_genes) {
    
    if (nrow(untrans_genes) == 0 || nrow(trans_genes) == 0) {
        return(data.frame())
    }
    
    common_genes <- intersect(untrans_genes$Hugo_Symbol, trans_genes$Hugo_Symbol)

    if (length(common_genes) == 0) {
        return(data.frame())
    }

    enriched <- data.frame()

    for (gene in common_genes) {
        # Check if freq column exists, if not try other frequency-related columns
        if ("freq" %in% colnames(untrans_genes)) {
            untrans_freq <- untrans_genes$freq[untrans_genes$Hugo_Symbol == gene]
            trans_freq <- trans_genes$freq[trans_genes$Hugo_Symbol == gene]
        } else if ("Frequency" %in% colnames(untrans_genes)) {
            untrans_freq <- untrans_genes$Frequency[untrans_genes$Hugo_Symbol == gene]
            trans_freq <- trans_genes$Frequency[trans_genes$Hugo_Symbol == gene]
        } else {
            # Calculate frequency based on sample counts if available
            message("Frequency column not found, skipping gene-level comparison")
            return(data.frame())
        }

        if (length(trans_freq) > 0 && length(untrans_freq) > 0 && 
            trans_freq > untrans_freq + 0.2) { # >20% increase
            enriched <- rbind(enriched, data.frame(
                gene = gene,
                untransformed_freq = round(untrans_freq, 3),
                transformed_freq = round(trans_freq, 3),
                frequency_difference = round(trans_freq - untrans_freq, 3),
                fold_enrichment = ifelse(untrans_freq > 0, round(trans_freq / untrans_freq, 2), Inf),
                stringsAsFactors = FALSE
            ))
        }
    }

    return(enriched)
}

## Function to perform comprehensive GISTIC comparison between transformed and untransformed groups
compare_gistic_transformation <- function(untransformed_gistic, transformed_gistic,
                                        output_prefix = "fst_transformation",
                                        freq_threshold = 0.2, qval_threshold = 0.25) {
    
    library(dplyr)
    library(ggplot2)
    
    message("Starting comprehensive GISTIC transformation analysis...")
    
    ## Create results directory if it doesn't exist
    dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)
    dir.create(here("figures/wes/gistic2/comparative"), recursive = TRUE, showWarnings = FALSE)
    
    ## Get sample counts for proper frequency calculation
    untrans_total_samples <- length(getSampleSummary(untransformed_gistic)$Tumor_Sample_Barcode)
    trans_total_samples <- length(getSampleSummary(transformed_gistic)$Tumor_Sample_Barcode)
    
    message(paste("Untransformed samples:", untrans_total_samples))
    message(paste("Transformed samples:", trans_total_samples))
    
    ## Extract cytoband-level data
    gistic_summary_untrans <- getCytobandSummary(untransformed_gistic)
    gistic_summary_trans <- getCytobandSummary(transformed_gistic)
    
    ## Separate amplifications and deletions
    untrans_amp <- gistic_summary_untrans |> filter(Variant_Classification == "Amp")
    untrans_del <- gistic_summary_untrans |> filter(Variant_Classification == "Del")
    trans_amp <- gistic_summary_trans |> filter(Variant_Classification == "Amp")
    trans_del <- gistic_summary_trans |> filter(Variant_Classification == "Del")

    ## Compare amplifications and deletions
    amp_comparison <- compare_cna_by_type(untrans_amp, trans_amp, "Amplification",
                                         untrans_total_samples, trans_total_samples,
                                         freq_threshold, qval_threshold)
    
    del_comparison <- compare_cna_by_type(untrans_del, trans_del, "Deletion",
                                         untrans_total_samples, trans_total_samples,
                                         freq_threshold, qval_threshold)
    
    ## Combine cytoband results
    cytoband_results <- rbind(amp_comparison, del_comparison)
    
    ## Gene-level analysis
    gene_results <- compare_genes_between_groups(untransformed_gistic, transformed_gistic,
                                                freq_threshold)
    
    ## Statistical testing (Fisher's exact test for each cytoband)
    cytoband_results <- add_statistical_testing(cytoband_results)
    
    ## Save results
    write.table(cytoband_results,
               here("results", paste0(output_prefix, "_cytoband_comparison.txt")),
               sep = "\t", row.names = FALSE, quote = FALSE)
    
    saveRDS(gene_results, here("results", paste0(output_prefix, "_gene_comparison.rds")))
    
    ## Create summary tables
    transformation_enriched <- cytoband_results |>
        filter(transformation_enriched == TRUE, significant_in_transformed == TRUE) |>
        arrange(desc(frequency_difference))
    
    transformation_specific <- cytoband_results |>
        filter(untransformed_count == 0, transformed_count > 0, significant_in_transformed == TRUE) |>
        arrange(desc(transformed_freq))
    
    ## Save summary tables
    write.table(transformation_enriched,
               here("results", paste0(output_prefix, "_enriched_cytobands.txt")),
               sep = "\t", row.names = FALSE, quote = FALSE)
    
    write.table(transformation_specific,
               here("results", paste0(output_prefix, "_specific_cytobands.txt")),
               sep = "\t", row.names = FALSE, quote = FALSE)
    
    ## Create visualization
    create_transformation_plots(cytoband_results, output_prefix)
    
    ## Print summary
    message("=== FST Transformation Analysis Summary ===")
    message(paste("Total cytobands analyzed:", nrow(cytoband_results)))
    message(paste("Transformation-enriched cytobands:", nrow(transformation_enriched)))
    message(paste("Transformation-specific cytobands:", nrow(transformation_specific)))
    message(paste("Transformation-specific amplified genes:", length(gene_results$amp_specific)))
    message(paste("Transformation-specific deleted genes:", length(gene_results$del_specific)))
    message(paste("Enriched amplified genes:", nrow(gene_results$amp_enriched)))
    message(paste("Enriched deleted genes:", nrow(gene_results$del_enriched)))
    
    ## Return comprehensive results
    results <- list(
        cytoband_comparison = cytoband_results,
        transformation_enriched = transformation_enriched,
        transformation_specific = transformation_specific,
        gene_comparison = gene_results,
        summary = list(
            untransformed_samples = untrans_total_samples,
            transformed_samples = trans_total_samples,
            enriched_cytobands = nrow(transformation_enriched),
            specific_cytobands = nrow(transformation_specific)
        )
    )
    
    return(results)
}

## Helper function to compare CNAs by type
compare_cna_by_type <- function(untrans_data, trans_data, alteration_type,
                               untrans_total_samples, trans_total_samples,
                               freq_threshold, qval_threshold) {
    
    ## Get all unique cytobands
    all_cytobands <- unique(c(untrans_data$Cytoband, trans_data$Cytoband))
    
    comparison_results <- data.frame()
    
    for (cytoband in all_cytobands) {
        ## Get sample counts for this cytoband
        untrans_count <- ifelse(cytoband %in% untrans_data$Cytoband,
                               untrans_data$nSamples[untrans_data$Cytoband == cytoband], 0
        )
        trans_count <- ifelse(cytoband %in% trans_data$Cytoband,
                             trans_data$nSamples[trans_data$Cytoband == cytoband], 0
        )
        
        ## Calculate actual frequencies
        untrans_freq <- untrans_count / untrans_total_samples
        trans_freq <- trans_count / trans_total_samples
        freq_diff <- trans_freq - untrans_freq
        
        ## Get q-values
        untrans_qval <- ifelse(cytoband %in% untrans_data$Cytoband,
                              untrans_data$qvalues[untrans_data$Cytoband == cytoband], 1
        )
        trans_qval <- ifelse(cytoband %in% trans_data$Cytoband,
                            trans_data$qvalues[trans_data$Cytoband == cytoband], 1
        )
        
        ## Get additional information
        untrans_genes <- ifelse(cytoband %in% untrans_data$Cytoband,
                               untrans_data$nGenes[untrans_data$Cytoband == cytoband], 0
        )
        trans_genes <- ifelse(cytoband %in% trans_data$Cytoband,
                             trans_data$nGenes[trans_data$Cytoband == cytoband], 0
        )
        
        untrans_peak <- ifelse(cytoband %in% untrans_data$Cytoband,
                              as.character(untrans_data$Wide_Peak_Limits[untrans_data$Cytoband == cytoband]), "")
        trans_peak <- ifelse(cytoband %in% trans_data$Cytoband,
                            as.character(trans_data$Wide_Peak_Limits[trans_data$Cytoband == cytoband]), ""
        )

        comparison_results <- rbind(comparison_results, data.frame(
            cytoband = cytoband,
            alteration_type = alteration_type,
            untransformed_count = untrans_count,
            transformed_count = trans_count,
            untransformed_total = untrans_total_samples,
            transformed_total = trans_total_samples,
            untransformed_freq = round(untrans_freq, 3),
            transformed_freq = round(trans_freq, 3),
            frequency_difference = round(freq_diff, 3),
            fold_change = ifelse(untrans_freq > 0, round(trans_freq / untrans_freq, 2), Inf),
            untransformed_qval = untrans_qval,
            transformed_qval = trans_qval,
            untransformed_genes = untrans_genes,
            transformed_genes = trans_genes,
            untransformed_peak = untrans_peak,
            transformed_peak = trans_peak,
            transformation_enriched = freq_diff > freq_threshold,
            significant_in_transformed = trans_qval < qval_threshold,
            significant_in_untransformed = untrans_qval < qval_threshold,
            transformation_specific = untrans_count == 0 & trans_count > 0,
            stringsAsFactors = FALSE
        ))
    }
    
    return(comparison_results)
}

## Helper function to compare genes between groups
compare_genes_between_groups <- function(untransformed_gistic, transformed_gistic, freq_threshold) {
    
    ## Initialize results
    amp_specific <- character(0)
    del_specific <- character(0)
    amp_enriched <- data.frame()
    del_enriched <- data.frame()
    
    ## Try to extract gene data from GISTIC objects using getGeneSummary()
    tryCatch({
        ## Get gene-level data using getGeneSummary function
        untrans_gene_summary <- getGeneSummary(untransformed_gistic)
        trans_gene_summary <- getGeneSummary(transformed_gistic)
        

        ## Separate amplifications and deletions
        untrans_amp_genes <- untrans_gene_summary |> filter(Del == 0)
        untrans_del_genes <- untrans_gene_summary |> filter(Amp == 0)
        trans_amp_genes <- trans_gene_summary |> filter(Del == 0)
        trans_del_genes <- trans_gene_summary |> filter(Amp == 0)
        
        ## Compare amplified genes
        if (nrow(untrans_amp_genes) > 0 && nrow(trans_amp_genes) > 0) {
            amp_specific <- setdiff(trans_amp_genes$Hugo_Symbol, untrans_amp_genes$Hugo_Symbol)
            amp_enriched <- find_enriched_genes_helper(untrans_amp_genes, trans_amp_genes, freq_threshold)
        } else if (nrow(trans_amp_genes) > 0) {
            # If no untransformed amplified genes, all transformed genes are specific
            amp_specific <- trans_amp_genes$Hugo_Symbol
        }
        
        ## Compare deleted genes
        if (nrow(untrans_del_genes) > 0 && nrow(trans_del_genes) > 0) {
            del_specific <- setdiff(trans_del_genes$Hugo_Symbol, untrans_del_genes$Hugo_Symbol)
            del_enriched <- find_enriched_genes_helper(untrans_del_genes, trans_del_genes, freq_threshold)
        } else if (nrow(trans_del_genes) > 0) {
            # If no untransformed deleted genes, all transformed genes are specific
            del_specific <- trans_del_genes$Hugo_Symbol
        }
        
    }, error = function(e) {
        message("Could not extract gene-level data from GISTIC objects: ", e$message)
        message("Skipping gene-level comparison")
    })
    
    return(list(
        amp_specific = amp_specific,
        del_specific = del_specific,
        amp_enriched = amp_enriched,
        del_enriched = del_enriched
    ))
}

## Helper function to find enriched genes
find_enriched_genes_helper <- function(untrans_genes, trans_genes, freq_threshold) {
    
    if (is.null(untrans_genes) || is.null(trans_genes) || 
        nrow(untrans_genes) == 0 || nrow(trans_genes) == 0) {
        return(data.frame())
    }
    
    common_genes <- intersect(untrans_genes$Hugo_Symbol, trans_genes$Hugo_Symbol)
    if (length(common_genes) == 0) {
        return(data.frame())
    }
    
    enriched <- data.frame()
    
    for (gene in common_genes) {
        untrans_idx <- which(untrans_genes$Hugo_Symbol == gene)
        trans_idx <- which(trans_genes$Hugo_Symbol == gene)
        
        if (length(untrans_idx) > 0 && length(trans_idx) > 0) {
            # GISTIC gene data typically has frequency information
            # Check for different possible column names
            if ("freq" %in% colnames(untrans_genes)) {
                untrans_freq <- untrans_genes$freq[untrans_idx[1]]
                trans_freq <- trans_genes$freq[trans_idx[1]]
            } else if ("frequency" %in% colnames(untrans_genes)) {
                untrans_freq <- untrans_genes$frequency[untrans_idx[1]]
                trans_freq <- trans_genes$frequency[trans_idx[1]]
            } else if ("Frequency" %in% colnames(untrans_genes)) {
                untrans_freq <- untrans_genes$Frequency[untrans_idx[1]]
                trans_freq <- trans_genes$Frequency[trans_idx[1]]
            } else {
                # If no frequency column, skip this comparison
                message("No frequency column found in gene data, skipping enrichment analysis")
                return(data.frame())
            }
            
            # Calculate frequency difference
            freq_diff <- trans_freq - untrans_freq
            
            if (!is.na(freq_diff) && freq_diff > freq_threshold) {
                enriched <- rbind(enriched, data.frame(
                    gene = gene,
                    untransformed_freq = round(untrans_freq, 3),
                    transformed_freq = round(trans_freq, 3),
                    frequency_difference = round(freq_diff, 3),
                    fold_enrichment = ifelse(untrans_freq > 0, round(trans_freq / untrans_freq, 2), Inf),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }
    
    return(enriched)
}

## Helper function to add Fisher's exact test
add_statistical_testing <- function(cytoband_results) {
    
    cytoband_results$fisher_pvalue <- NA
    cytoband_results$fisher_odds_ratio <- NA
    
    for (i in 1:nrow(cytoband_results)) {
        row <- cytoband_results[i, ]
        
        ## Create contingency table
        contingency_table <- matrix(c(
            row$transformed_count, row$transformed_total - row$transformed_count,
            row$untransformed_count, row$untransformed_total - row$untransformed_count
        ), nrow = 2, byrow = TRUE)
        
        ## Perform Fisher's exact test
        if (sum(contingency_table) > 0) {
            fisher_result <- fisher.test(contingency_table)
            cytoband_results$fisher_pvalue[i] <- fisher_result$p.value
            cytoband_results$fisher_odds_ratio[i] <- fisher_result$estimate
        }
    }
    
    ## Add FDR correction
    cytoband_results$fisher_qvalue <- p.adjust(cytoband_results$fisher_pvalue, method = "BH")
    
    return(cytoband_results)
}

## Helper function to create visualization plots
create_transformation_plots <- function(cytoband_results, output_prefix) {
    
    ## Filter for significant results
    significant_results <- cytoband_results |>
        filter(transformation_enriched == TRUE, significant_in_transformed == TRUE) |>
        arrange(desc(frequency_difference)) |>
        head(20)  # Top 20 for visualization
    
    if (nrow(significant_results) > 0) {
        ## Create bar plot
        p <- ggplot(significant_results, aes(x = reorder(cytoband, frequency_difference), 
                                           y = frequency_difference)) +
            geom_col(aes(fill = alteration_type)) +
            coord_flip() +
            scale_fill_manual(values = c("Amplification" = "#D95F02", "Deletion" = "#1B9E77")) +
            labs(title = "FST Transformation-Enriched CNAs",
                 subtitle = paste("Top", nrow(significant_results), "cytobands with >20% frequency difference"),
                 x = "Cytogenetic Band",
                 y = "Frequency Difference (Transformed - Untransformed)",
                 fill = "Alteration Type") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 8))
        
        ## Save plot
        for (img in c("png", "pdf")) {
            file <- here("figures/wes/gistic2/comparative", 
                        paste0(output_prefix, "_enriched_cnas.", img))
            
            if (img == "png") {
                CairoPNG(file = file, width = 10, height = 8, res = 300, units = "in")
            } else {
                CairoPDF(file = file, width = 10, height = 8)
            }
            
            print(p)
            dev.off()
        }
        
        message(paste("Saved plot:", file))
    } else {
        message("No significant transformation-enriched CNAs found for plotting")
    }
}

GetStatResGistic2 <- function(
    gistic_dir,
    group1,
    group2,
    freq_thres = 0.2,
    qval_thres = 0.25
) {

    ## Load GISTIC data in maftools objects
    g1_obj <- readGistic(
        gisticAllLesionsFile = here(gistic_dir, group1, "all_lesions.conf_99.txt"),
        gisticAmpGenesFile = here(gistic_dir, group1, "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, group1, "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, group1, "scores.gistic"),
        verbose = FALSE
    )

    g2_obj <- readGistic(
        gisticAllLesionsFile = here(gistic_dir, group2, "all_lesions.conf_99.txt"),
        gisticAmpGenesFile = here(gistic_dir, group2, "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, group2, "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, group2, "scores.gistic"),
        verbose = FALSE
    )

    # head(g1_obj@data)
    # head(g1_obj@gis.scores)

    ## Get sample counts for proper frequency calculation
    g1_total_samples <- length(getSampleSummary(g1_obj)$Tumor_Sample_Barcode)
    g2_total_samples <- length(getSampleSummary(g2_obj)$Tumor_Sample_Barcode)

    message(paste("Group 1 (", group1, ") samples:", g1_total_samples))
    message(paste("Group 2 (", group2, ") samples:", g2_total_samples))

    ## "==================================================================="
    ## Cytoband-level comparison ----
    ## "==================================================================="
    ## Extract cytoband-level data
    g1_cytoband <- getCytobandSummary(g1_obj)
    g2_cytoband <- getCytobandSummary(g2_obj)

    ## Separate amplifications and deletions
    g1_cytoband_amp <- g1_cytoband |> filter(Variant_Classification == "Amp")
    g1_cytoband_del <- g1_cytoband |> filter(Variant_Classification == "Del")
    
    g2_cytoband_amp <- g2_cytoband |> filter(Variant_Classification == "Amp")
    g2_cytoband_del <- g2_cytoband |> filter(Variant_Classification == "Del")

    ## Cytoband comparison
    cytoband_comparison_amp <- GetGistic2GroupCytobandData(
        cytoband1 = g1_cytoband_amp,
        cytoband2 = g2_cytoband_amp,
        alteration_type = "Amp",
        g1_total_samples = g1_total_samples,
        g2_total_samples = g2_total_samples,
        freq_thres = freq_thres,
        qval_thres = qval_thres
    )
    
    cytoband_comparison_del <- GetGistic2GroupCytobandData(
        cytoband1 = g1_cytoband_del,
        cytoband2 = g2_cytoband_del,
        alteration_type = "Del",
        g1_total_samples = g1_total_samples,
        g2_total_samples = g2_total_samples,
        freq_thres = freq_thres,
        qval_thres = qval_thres
    )
    
    cytoband_comparison_results <- rbind(
        cytoband_comparison_amp, 
        cytoband_comparison_del
    ) |> 
        as_tibble()

    ## Statistical testing (Fisher's exact test for each cytoband)
    cytoband_stat_res <- Gistic2CytobandStat(
        cytoband_comparison_results, 
        adj_method = "BH"
    )

    ## "==================================================================="
    ## Gene level comparison ----
    ## "==================================================================="
    ## Get gene-level data
    g1_gene <- getGeneSummary(g1_obj)
    g2_gene <- getGeneSummary(g2_obj)

    ## Separate amplifications and deletions
    g1_gene_amp <- g1_gene |> filter(Amp !=0)
    g1_gene_del <- g1_gene |> filter(Del !=0)

    g2_gene_amp <- g2_gene |> filter(Amp !=0)
    g2_gene_del <- g2_gene |> filter(Del !=0)
    
    ## Compare gene-level amplifications and deletions
    gene_comparison_amp <- Gistic2CompareGenes(
        gene1 = g1_gene_amp,
        gene2 = g2_gene_amp,
        alteration_type = "Amp",
        g1_total_samples = g1_total_samples,
        g2_total_samples = g2_total_samples,
        freq_thres = freq_thres
    )
    
    gene_comparison_del <- Gistic2CompareGenes(
        gene1 = g1_gene_del,
        gene2 = g2_gene_del,
        alteration_type = "Del",
        g1_total_samples = g1_total_samples,
        g2_total_samples = g2_total_samples,
        freq_thres = freq_thres
    )
    
    gene_comparison_results <- rbind(
        gene_comparison_amp, 
        gene_comparison_del
    ) |> 
        as_tibble()
    
    ## Statistical testing for gene-level results
    gene_stat_res <- Gistic2GeneStat(
        gene_comparison_results, 
        adj_method = "BH"
    )
    
    cytoband_sig <- cytoband_stat_res |> 
        dplyr::filter(significantly_different)
    
    cytoband_sig_g2 <- cytoband_sig |> 
        filter(enriched_in_group2)
    
    # cytoband_sig_g2_bands <- cytoband_sig_g2$cytoband
    
    # g2_obj@data |> 
    #     as_tibble() |> 
    #     mutate(bands = str_extract(Cytoband, "(?<=:).+")) |>
    #     filter(bands %in% cytoband_sig_g2_bands) |>

    # gene_sig <- gene_stat_res |> 
    #     dplyr::filter(significantly_different)

    ## Return comprehensive results
    return(
        list(
            cytoband = cytoband_stat_res,
            cytoband_sig = cytoband_sig,
            gene = gene_stat_res,
            gene_sig = gene_sig
        )
    )

}

GetGistic2CytobandStatRes <- function(
    gistic_dir,
    group1,
    group2,
    freq_thres = 0.2,
    qval_thres = 0.25
) {

    ## Load GISTIC data in maftools objects
    g1_obj <- LoadGisticData(
        gistic_dir = gistic_dir,
        group = group1
    )

    g2_obj <- LoadGisticData(
        gistic_dir = gistic_dir,
        group = group2
    )

    # head(g1_obj@data)
    # head(g1_obj@gis.scores)
    # sample_groups <- LoadSampleGroupInfo(is_somatic_matched = is_somatic_matched)
    
    # g1_total_sampel = length(sample_groups[[group1]])

    ## Get sample counts for proper frequency calculation
    g1_total_samples <- length(getSampleSummary(g1_obj)$Tumor_Sample_Barcode)
    g2_total_samples <- length(getSampleSummary(g2_obj)$Tumor_Sample_Barcode)

    message(paste(" - Group 1 (", group1, ") samples:", g1_total_samples))
    message(paste(" - Group 2 (", group2, ") samples:", g2_total_samples))

    ## "==================================================================="
    ## Cytoband-level comparison ----
    ## "==================================================================="
    ## Extract cytoband-level data
    g1_cytoband <- getCytobandSummary(g1_obj)
    g2_cytoband <- getCytobandSummary(g2_obj)

    ## Separate amplifications and deletions
    g1_cytoband_amp <- g1_cytoband |> filter(Variant_Classification == "Amp")
    g1_cytoband_del <- g1_cytoband |> filter(Variant_Classification == "Del")
    
    g2_cytoband_amp <- g2_cytoband |> filter(Variant_Classification == "Amp")
    g2_cytoband_del <- g2_cytoband |> filter(Variant_Classification == "Del")

    ## Cytoband comparison
    cytoband_comparison_amp <- GetGistic2GroupCytobandData(
        cytoband1 = g1_cytoband_amp,
        cytoband2 = g2_cytoband_amp,
        alteration_type = "Amp",
        g1_total_samples = g1_total_samples,
        g2_total_samples = g2_total_samples,
        freq_thres = freq_thres,
        qval_thres = qval_thres
    )
    
    cytoband_comparison_del <- GetGistic2GroupCytobandData(
        cytoband1 = g1_cytoband_del,
        cytoband2 = g2_cytoband_del,
        alteration_type = "Del",
        g1_total_samples = g1_total_samples,
        g2_total_samples = g2_total_samples,
        freq_thres = freq_thres,
        qval_thres = qval_thres
    )
    
    cytoband_comparison_results <- rbind(
        cytoband_comparison_amp, 
        cytoband_comparison_del
    ) |> 
        as_tibble()

    ## Statistical testing (Fisher's exact test for each cytoband)
    cytoband_stat_res <- Gistic2CytobandStat(
        cytoband_comparison_results, 
        adj_method = "BH"
    )
    
    # cytoband_sig <- cytoband_stat_res |> 
    #     dplyr::filter(significantly_different)
    
    # cytoband_sig_g2 <- cytoband_sig |> 
    #     filter(enriched_in_group2)

    ## return
    cytoband_stat_res

}

GetGistic2GroupCytobandData <- function(
    cytoband1, 
    cytoband2, 
    alteration_type,
    g1_total_samples, 
    g2_total_samples, 
    freq_thres, 
    qval_thres
) {
    
    ## Get all unique cytobands
    all_cytobands <- unique(
        c(cytoband1$Cytoband, cytoband2$Cytoband)
    )

    comparison_results <- data.frame()

    ## Loop through each cytoband
    for (cytoband in all_cytobands) {
        ## Get sample counts for this cytoband - using safer indexing
        count1 <- 0
        count2 <- 0
        
        ## Safely extract counts for group1
        if (cytoband %in% cytoband1$Cytoband) {
            idx1 <- which(cytoband1$Cytoband == cytoband)[1]  # Take first match
            count1 <- cytoband1$nSamples[idx1]
        }
        
        ## Safely extract counts for group2
        if (cytoband %in% cytoband2$Cytoband) {
            idx2 <- which(cytoband2$Cytoband == cytoband)[1]  # Take first match
            count2 <- cytoband2$nSamples[idx2]
        }

        ## Calculate actual frequencies
        freq1 <- count1 / g1_total_samples
        freq2 <- count2 / g2_total_samples
        freq_diff <- freq2 - freq1

        ## Get q-values - using safer indexing
        qval1 <- 1
        qval2 <- 1
        
        ## Safely extract q-values for group1
        if (cytoband %in% cytoband1$Cytoband) {
            idx1 <- which(cytoband1$Cytoband == cytoband)[1]  # Take first match
            qval1 <- cytoband1$qvalues[idx1]
        }
        
        ## Safely extract q-values for group2
        if (cytoband %in% cytoband2$Cytoband) {
            idx2 <- which(cytoband2$Cytoband == cytoband)[1]  # Take first match
            qval2 <- cytoband2$qvalues[idx2]
        }

        ## Get gene counts - using safer indexing
        genes1 <- 0
        genes2 <- 0
        
        ## Safely extract gene counts for group1
        if (cytoband %in% cytoband1$Cytoband) {
            idx1 <- which(cytoband1$Cytoband == cytoband)[1]  # Take first match
            genes1 <- cytoband1$nGenes[idx1]
        }
        
        ## Safely extract gene counts for group2
        if (cytoband %in% cytoband2$Cytoband) {
            idx2 <- which(cytoband2$Cytoband == cytoband)[1]  # Take first match
            genes2 <- cytoband2$nGenes[idx2]
        }

        ## Get wide peak limits - using safer indexing
        peak1 <- ""
        peak2 <- ""
        
        ## Safely extract peak limits for group1
        if (cytoband %in% cytoband1$Cytoband) {
            idx1 <- which(cytoband1$Cytoband == cytoband)[1]  # Take first match
            peak1 <- as.character(cytoband1$Wide_Peak_Limits[idx1])
        }
        
        ## Safely extract peak limits for group2
        if (cytoband %in% cytoband2$Cytoband) {
            idx2 <- which(cytoband2$Cytoband == cytoband)[1]  # Take first match
            peak2 <- as.character(cytoband2$Wide_Peak_Limits[idx2])
        }

        comparison_results <- rbind(
            comparison_results,
            data.frame(
                cytoband = cytoband,
                alteration_type = alteration_type,
                group1_count = count1,
                group2_count = count2,
                group1_total = g1_total_samples,
                group2_total = g2_total_samples,
                group1_freq = round(freq1, 3),
                group2_freq = round(freq2, 3),
                freq_difference = round(freq_diff, 3),
                # fold_change = ifelse(freq1 > 0, round(freq2 / freq1, 2), Inf),
                group1_gistic2_qval = qval1,
                group2_gistic2_qval = qval2,
                group1_genes = genes1,
                group2_genes = genes2,
                group1_peak = peak1,
                group2_peak = peak2,
                # enriched = freq_diff > freq_thres,
                sig_in_group1 = qval1 < qval_thres,
                sig_in_group2 = qval2 < qval_thres,
                specific_in_group2 = count1 == 0 & count2 > 0,
                stringsAsFactors = FALSE
            )
        )        
    }
    ## return
    comparison_results
}

Gistic2CompareGenes <- function(
    gene1, gene2, alteration_type,
    g1_total_samples, g2_total_samples, freq_thres) {
    ## Get all unique genes
    all_genes <- unique(
        c(gene1$Hugo_Symbol, gene2$Hugo_Symbol)
    )

    comparison_results <- data.frame()

    ## Loop through each gene
    for (gene in all_genes) {
        
        ## Get sample counts for this gene - using AlteredSamples column
        count1 <- 0
        count2 <- 0
        
        ## Safely extract counts for group1
        if (gene %in% gene1$Hugo_Symbol) {
            idx1 <- which(gene1$Hugo_Symbol == gene)[1]  # Take first match
            count1 <- gene1$AlteredSamples[idx1]  # Use AlteredSamples instead of Amp/Del
        }
        
        ## Safely extract counts for group2
        if (gene %in% gene2$Hugo_Symbol) {
            idx2 <- which(gene2$Hugo_Symbol == gene)[1]  # Take first match
            count2 <- gene2$AlteredSamples[idx2]  # Use AlteredSamples instead of Amp/Del
        }
    

        ## Calculate actual frequencies
        freq1 <- count1 / g1_total_samples
        freq2 <- count2 / g2_total_samples
        freq_diff <- freq2 - freq1

        comparison_results <- rbind(
            comparison_results,
            data.frame(
                gene = gene,
                alteration_type = alteration_type,
                group1_count = count1,
                group2_count = count2,
                group1_total = g1_total_samples,
                group2_total = g2_total_samples,
                group1_freq = round(freq1, 3),
                group2_freq = round(freq2, 3),
                freq_difference = round(freq_diff, 3),
                fold_change = ifelse(freq1 > 0, round(freq2 / freq1, 2), Inf),
                enriched = freq_diff > freq_thres,
                specific_in_group2 = count1 == 0 & count2 > 0,
                stringsAsFactors = FALSE
            )
        )
    }
    ## return
    comparison_results
}

Gistic2GeneStat <- function(
    gene_comparison_results, 
    adj_method = "BH"
) {

    gene_comparison_results$fisher_pvalue <- NA
    gene_comparison_results$fisher_odds_ratio <- NA

    ## Loop through each gene comparison result
    for (i in 1:nrow(gene_comparison_results)) {
        
        row <- gene_comparison_results[i, ]
        
        ## Check for problematic values in the row
        if (!is.finite(row$group1_count) || row$group1_count < 0) {
            message("Warning: Gene ", row$gene, " has problematic group1_count: ", row$group1_count)
        }
        if (!is.finite(row$group2_count) || row$group2_count < 0) {
            message("Warning: Gene ", row$gene, " has problematic group2_count: ", row$group2_count)
        }
        if (!is.finite(row$group1_total) || row$group1_total < 1) {
            message("Warning: Gene ", row$gene, " has problematic group1_total: ", row$group1_total)
        }
        if (!is.finite(row$group2_total) || row$group2_total < 1) {
            message("Warning: Gene ", row$gene, " has problematic group2_total: ", row$group2_total)
        }
        
        ## Ensure all values are valid non-negative integers
        group1_count <- max(0, as.integer(row$group1_count))
        group2_count <- max(0, as.integer(row$group2_count))
        group1_total <- max(1, as.integer(row$group1_total))
        group2_total <- max(1, as.integer(row$group2_total))
        
        ## Create contingency table
        contingency_table <- matrix(
            c(
                group2_count, group2_total - group2_count,
                group1_count, group1_total - group1_count
            ),
            nrow = 2, 
            byrow = TRUE,
            dimnames = list(
                c("Group2", "Group1"),
                c("Altered", "Not_Altered")
            )
        )
        
        ## Check contingency table for negative values
        if (any(contingency_table < 0)) {
            message("Warning: Gene ", row$gene, " has negative values in contingency table:")
            print(contingency_table)
            ## Fix negative values
            contingency_table <- pmax(contingency_table, 0)
        }
        
        ## Perform Fisher's exact test only if contingency table is valid
        if (sum(contingency_table) > 0 && all(is.finite(contingency_table))) {
            tryCatch({
                fisher_result <- fisher.test(contingency_table)
                gene_comparison_results$fisher_pvalue[i] <- fisher_result$p.value
                gene_comparison_results$fisher_odds_ratio[i] <- fisher_result$estimate
            }, error = function(e) {
                message("Error in Fisher's test for gene ", row$gene, ": ", e$message)
                message("Contingency table:")
                print(contingency_table)
                gene_comparison_results$fisher_pvalue[i] <<- NA
                gene_comparison_results$fisher_odds_ratio[i] <<- NA
            })
        }
    }

    ## Add FDR correction
    gene_comparison_results$fisher_qvalue <- p.adjust(
        gene_comparison_results$fisher_pvalue,
        method = adj_method
    )
    
    ## Add interpretation columns for clarity
    gene_comparison_results <- gene_comparison_results |>
        mutate(
            significantly_different = fisher_qvalue < 0.05,
            enriched_in_group2 = fisher_qvalue < 0.05 & fisher_odds_ratio > 1 & freq_difference > 0,
            depleted_in_group2 = fisher_qvalue < 0.05 & fisher_odds_ratio < 1 & freq_difference < 0,
            effect_direction = case_when(
                !significantly_different ~ "No significant difference",
                enriched_in_group2 ~ "Enriched in Group2",
                depleted_in_group2 ~ "Depleted in Group2 (Enriched in Group1)",
                TRUE ~ "Unclear"
            )
        )

    gene_comparison_results
}

Gistic2CytobandStat <- function(
    cytoband_comparison_results, 
    adj_method = "BH"
) {

    cytoband_comparison_results$fisher_pvalue <- NA
    cytoband_comparison_results$fisher_odds_ratio <- NA

    ## Loop through each cytoband comparison result
    for (i in 1:nrow(cytoband_comparison_results)) {
        
        row <- cytoband_comparison_results[i, ]
        
        ## Create contingency table with group2 in first row for easier interpretation
        contingency_table <- matrix(
            c(
                row$group2_count, row$group2_total - row$group2_count,
                row$group1_count, row$group1_total - row$group1_count
            ),
            nrow = 2, 
            byrow = TRUE,
            dimnames = list(
                c("Group2", "Group1"),
                c("Altered", "Not_Altered")
            )
        )
        
        ## Perform Fisher's exact test
        if (sum(contingency_table) > 0) {
            fisher_result <- fisher.test(contingency_table)
            cytoband_comparison_results$fisher_pvalue[i] <- fisher_result$p.value
            cytoband_comparison_results$fisher_odds_ratio[i] <- fisher_result$estimate
        }
    }

    ## Add FDR correction
    cytoband_comparison_results$fisher_qvalue <- p.adjust(
        cytoband_comparison_results$fisher_pvalue,
        method = adj_method
    )
    
    ## Add interpretation columns for clarity
    cytoband_comparison_results <- cytoband_comparison_results |>
        mutate(
            significantly_different = fisher_qvalue < 0.05,
            enriched_in_group2 = fisher_qvalue < 0.05 & fisher_odds_ratio > 1 & freq_difference > 0,
            depleted_in_group2 = fisher_qvalue < 0.05 & fisher_odds_ratio < 1 & freq_difference < 0,
            effect_direction = case_when(
                !significantly_different ~ "No significant difference",
                enriched_in_group2 ~ "Enriched in Group2",
                depleted_in_group2 ~ "Depleted in Group2 (Enriched in Group1)",
                TRUE ~ "Unclear"
            )
        )

    cytoband_comparison_results
}

BoxPlotCompareGroups <- function(
    data, # a tibble
    x,
    y,
    y_lab,
    y_lim = NULL,
    title = NULL,
    stat_method = "wilcox_test", # t-test, wilcox_test
    alternative = "two.sided", # "greater", "less"
    step_increase = 0.03,
    p_vjust = 0
) {
    
    ## Do the statistical test
    data <- data |> 
        mutate(value = !!sym(y), group = !!sym(x))

    if (stat_method == "t_test") {

        stat_test <- data |>
            t_test(
                value ~ group,
                alternative = alternative,
                detailed = TRUE
            ) |>
            add_significance() |>
            add_xy_position(x = "group", dodge = 0.9)

    } else if (stat_method == "wilcox_test") {

        stat_test <- data |>
            wilcox_test(
                value ~ group,
                alternative = alternative,
                detailed = TRUE
            ) |>
            add_significance() |>
            add_xy_position(x = "group", dodge = 0.9)

    } else {

        stop("Unknown test")
    }

    stat_test <- stat_test |> 
        dplyr::mutate(
            p_label = ifelse(
                p < 0.05,
                paste0("p = ", scales::pvalue(p), "", p.signif),
                paste0("p = ", scales::pvalue(p))
            )
        )

    print(stat_test)

    ## Plot the boxplot
    ggboxplot(
        data,
        x = x, y = y,
        fill = x,
        size = 0.2,
        width = 0.5,
        alpha = 0.2,
        add = "jitter",
        add.params = list(size = 1, alpha = 0.8, shape = 21)
    ) +
        stat_pvalue_manual(
            stat_test,
            # label = "p.signif",
            # label = "p ={p}, n={n1 + n2}",
            # label = "p = {scales::pvalue(p)}, n = {n1 + n2}",
            label = "p_label",
            size = 0.1,
            label.size = 2,
            bracket.size = 0.3,
            tip.length = 0,
            # xmin = "xmin",
            # xmax = "xmax",
            # y.position = "y.position",
            step.increase = step_increase,
            bracket.nudge.y = 0.1,
            vjust = -0.2,
            hide.ns = FALSE
        ) +
        scale_y_continuous(
            limits = y_lim,
            expand = expansion(mult = c(0.15))
            # expand = c(0, 0.2)
        ) +
        scale_fill_manual(
            # name = NULL,
            # labels = c("Control", "PV"),
            values = c("#1f77b4", "#ff7f0e")
        ) +
        labs(
            x = NULL,
            y = y_lab,
            title = title
        ) +
        figure_theme2() +
        theme(
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title.y = element_text(vjust = 0, hjust = 0.5),
            # legend.position = "bottom",
            # legend.box.spacing = unit(0, "pt"),
            plot.margin = unit(c(0.3, 0.1, 0.1, 0.1), units = "cm")
        )

}

PlotGisticGroupComparsion <- function(
    gistic_dir = "data/wes/gistic2/somatic_matched",
    groups = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"),
    is_somatic_matched = FALSE,
    fdrCutOff = 0.25,
    markBands = NULL,
    color = c("#c82118", "#2c5496"),
    ref.build = "hg38",
    cytobandOffset = 0.01,
    txtSize = 0.6,
    cytobandTxtSize = 0.3,
    maf = NULL,
    mutGenes = NULL,
    y_lims = NULL,
    mutGenesTxtSize = 0.6,
    width = 8,
    height = 8,
    plot_dir = "figures/gistic2",
    filename = "wes_cnvfacets_gistic2_FST.subgroup"
) {
    
    ## Load group info
    sample_groups <- LoadSampleGroupInfo(
        is_somatic_matched = is_somatic_matched
    )

    dir_create(plot_dir)

    ## plot the chromosomal plots
    for (img in c("png", "pdf")) {
        
        file <- here(plot_dir, paste0(filename, ".", img))

        if (img == "png") {
            CairoPNG(
                file = file, width = width, height = height, res = 300,
                fonts = "Arial", units = "in"
            )
        } else {
            CairoPDF(
                file = file, width = width, height = height, family  = "Arial"
            )
        }

        # Set up 4 rows, 1 column layout
        nrow = length(groups)

        par(mfrow = c(nrow, 1), mar = c(3, 4, 3, 2))

        # Read GISTIC objects for both groups
        gistic_obj_list <- list()
        for (group in groups) {

            gistic_obj_list[[group]] <- readGistic(
                gisticAllLesionsFile = here(gistic_dir, group, "all_lesions.conf_99.txt"),
                gisticAmpGenesFile = here(gistic_dir, group, "amp_genes.conf_99.txt"),
                gisticDelGenesFile = here(gistic_dir, group, "del_genes.conf_99.txt"),
                gisticScoresFile = here(gistic_dir, group, "scores.gistic"),
                verbose = FALSE
            )

        }

        ## Plot the chromosomal plots for each group
        for (group in names(gistic_obj_list)) {
            
            # n_sample <- length(sample_groups[[group]])
            n_sample <- gistic_obj_list[[group]]@summary |> 
                filter(ID =="Samples") |> 
                pull(summary)

            # Plot the chromosomal plot
            gisticChromPlot(
                gistic = gistic_obj_list[[group]],
                fdrCutOff = fdrCutOff,
                markBands = markBands,
                color = color,
                ref.build = ref.build,
                cytobandOffset = cytobandOffset,
                txtSize = txtSize,
                cytobandTxtSize = cytobandTxtSize,
                maf = maf,
                mutGenes = mutGenes,
                mutGenesTxtSize = mutGenesTxtSize,
                y_lims = y_lims
            )
            title(main = paste0(group, " (n=", n_sample, ")"), family = "Arial")
        }
        # Reset par settings
        par(mfrow = c(1, 1))
        
        message(paste0("Saving combined plot: ", file))
        dev.off()
    }

}

GenerateCytobandOncoplot <- function(
    mat,
    sample_annotation = c("FST.Group", "Metastasis"),
    sample_sorted_by = "FST.Group",
    sample_sorted_level = c(
        "U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"
    ),
    width = 18,
    height = 12,
    dir = "figures/oncoplot",
    filename = "oncoplot"
) {
    ## Convert the matrix to a tibble
    data <- as_tibble(mat, rownames = "cytoband") |> 
        mutate(across(everything(), as.character)) |> 
        pivot_longer(
            cols = -cytoband, 
            names_to = "Sample.ID", 
            values_to = "alteration"
        ) |> 
        distinct() |> 
        mutate(
            alteration = if_else(
                alteration == "", NA_character_, alteration
            )
        ) |> 
        filter(!is.na(alteration))
    

    ## Get matched clinical info
    clinical_info <- LoadClinicalInfo()
    
    sample_idx <- match(
        colnames(mat), 
        clinical_info$Sample.ID
    )
    
    ## Sort samples
    if (!is.null(sample_sorted_by) && !is.null(sample_sorted_level)) {

        # data |> 
        #     left_join(
        #         clinical_info  |> select(Sample.ID, all_of(sample_sorted_by)), 
        #         by = "Sample.ID"
        #     ) |>
        #     group_by(FST.Group)
        #     arrange()

        sample_order <- clinical_info[sample_idx, ] |> 
            # group_by(across(all_of(sample_sorted_by))) |>
            mutate(
                sample_order = factor(
                    !!sym(sample_sorted_by), 
                    levels = sample_sorted_level
                )
            ) |>
            group_by(sample_order) |>
            arrange(!!sym(sample_sorted_by), Sample.ID, .by_group = TRUE) |> 
            select(Sample.ID, all_of(sample_sorted_by), sample_order) |> 
            pull(Sample.ID)

    } else {

        ## Sort by patient
        sample_order <- colnames(mat)
    }

    sort_mat <- mat[, sample_order]

    ## Alteration colors
    data <- as_tibble(sort_mat, rownames = "cytoband") |> 
        mutate(across(everything(), as.character)) |> 
        pivot_longer(
            cols = -cytoband, 
            names_to = "Sample.ID", 
            values_to = "alteration"
        ) |> 
        distinct() |> 
        mutate(
            alteration = if_else(
                alteration == "", NA_character_, alteration
            )
        ) |> 
        filter(!is.na(alteration))
    
    ## Sample annotations
    sort_sample_idx <- match(
        sample_order, 
        clinical_info$Sample.ID
    )

    sample_annotation_df <- clinical_info[sort_sample_idx, ] |> 
        select(Sample.ID, all_of(sample_annotation)) |> 
        column_to_rownames("Sample.ID")

    ## Define colors for the sample annotations
    sample_annotation_colors <- list()

    for (annotation in sample_annotation) {
        sample_annotation_colors[[annotation]] <- study_colors[[annotation]]
    }

    ## Create the matrix annotation colors
    alterations <- data |> pull(alteration) |> unique() |> sort()

    ## Get the matched alteration colors
    alteration_colors <- c()
    
    for (i in alterations) {
        alteration_colors <- c(alteration_colors, study_colors$cna_class[i])
    }

    alter_fun <- list(
        background = alter_graphic("rect", fill = "#CCCCCC"),
        homdel = alter_graphic("rect", fill = alteration_colors[["homdel"]]),
        hetdel = alter_graphic("rect", fill = alteration_colors[["hetdel"]]),
        gain = alter_graphic("rect", fill = alteration_colors[["gain"]])
    )

    ## Plot the oncoprint
    oncoplot <- oncoPrint(
        sort_mat,
        alter_fun = alter_fun,
        col = alteration_colors,
        
        bottom_annotation = HeatmapAnnotation(
            df = sample_annotation_df,
            col = sample_annotation_colors,
            annotation_height = unit(c(5, 5), "mm"),
            annotation_name_gp = gpar(fontsize = 8)
        ),

        show_pct = TRUE,
        pct_digits = 0,
        pct_side = "right",

        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 6),
        # row_title = "Cytobands (Most Frequent → Least Frequent)",

        column_title = "Significant altered cytobands",
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        
        column_order = sample_order,
        # Add frequency annotation on the right
        # right_annotation = rowAnnotation(
        #     Frequency = anno_barplot(
        #         cytoband_freq$alteration_frequency[
        #             match(rownames(oncoplot_matrix_sorted), cytoband_freq$matched_cytoband)
        #         ],
        #         bar_width = 1,
        #         gp = gpar(fill = "#666666"),
        #         axis_param = list(side = "bottom", labels_rot = 0)
        #     ),
        #     annotation_width = unit(2, "cm")
        # ),
        
        heatmap_legend_param = list(
            title = "CNV Alteration",
            labels = c(
                "Homozygous Deletion", 
                "Heterozygous Deletion", 
                "Gain", 
                "No Alteration"
            )
        )
    )

    ## Save the oncoplot
    file <- here(dir, paste0(filename, ".pdf"))
    pdf(file, width = width, height = height)

    draw(
        oncoplot,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legends = TRUE,
        padding = unit(c(1, 1, 1, 1), "cm")
    )

    dev.off()

    message(
        paste0("Saved oncoplot: ", file)
    )

}

LoadGisticData <- function(
    gistic_dir,
    group,
    gistic_qval_thres = 0.5
) {
    
    ## Load GISTIC data in maftools objects
    gisticAllLesionsFile <- here(
        gistic_dir, 
        group, 
        paste0("qval_", gistic_qval_thres),
        "all_lesions.conf_99.txt"
    )

    gisticAmpGenesFile <- here(
        gistic_dir,
        group,
        paste0("qval_", gistic_qval_thres),
        "amp_genes.conf_99.txt"
    )

    gisticDelGenesFile <- here(
        gistic_dir,
        group,
        paste0("qval_", gistic_qval_thres),
        "del_genes.conf_99.txt"
    )
    
    gisticScoresFile <- here(
        gistic_dir,
        group,
        paste0("qval_", gistic_qval_thres),
        "scores.gistic"
    )

    gistic_obj <- readGistic(
        gisticAllLesionsFile = gisticAllLesionsFile,
        gisticAmpGenesFile = gisticAmpGenesFile,
        gisticDelGenesFile = gisticDelGenesFile,
        gisticScoresFile = gisticScoresFile,
        verbose = FALSE
    )

    ## Get sample counts for proper frequency calculation
    # total_samples <- length(getSampleSummary(gistic_obj)$Tumor_Sample_Barcode)

    # message(paste(" - Group (", group, ") samples:", total_samples))

    ## Return
    gistic_obj
}




GetPCGRCNVGroupDiffGene <- function(
    data,
    var = "FST.Group",
    cn_column = "cn_status",
    event_type = NULL,
    somatic_matched = TRUE,
    group1 = "U-DFSP",
    group2 = "FS-DFSP",
    p_adjust_method = "BH", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
    min_altered_samples = 3,
    group_comparison = "U-DFSP_vs_FS-DFSP"
) {

    ## Filter by event type
    if (!is.null(event_type)) {
        data <- data |> 
            filter(event_type %in% {{event_type}})
    }

    ## We should use the total cohort size as the denominator for frequency 
    ## calculations as we are answer the question "What proportion of samples 
    ## have this alteration in my cohort?"
    
    if (somatic_matched) {

        clinical_info <- LoadClinicalInfo() |> 
            filter(Somatic.Status == "Matched")

    } else {

        clinical_info <- LoadClinicalInfo()

    }
    
    ## Total samples of each group in cohort
    n_sample <- nrow(clinical_info)
    
    n_group1 <- clinical_info |> 
        filter(!!sym(var) %in% c(group1)) |> 
        nrow()
    
    n_group2 <- clinical_info |> 
        filter(!!sym(var) %in% c(group2)) |> 
        nrow()

    ## Calculate gene alteration frequencies
    freq_data <- data |>
        filter(!!sym(var) %in% c(group1, group2)) |> 
        select(sample_id, symbol, !!sym(cn_column), !!sym(var)) |>
        distinct() |> 
        group_by(symbol, !!sym(cn_column)) |> 
        summarise(
            group1_altered = sum(!!sym(var) %in% group1),
            group1_total = n_group1,
            group2_altered = sum(!!sym(var) %in% group2),
            group2_total = n_group2,
            total_altered = group1_altered + group2_altered,
            .groups = "drop"
        ) |>
        arrange(desc(total_altered)) |>
        ## Add back gene-specific summary statistics
        left_join(
            data |> 
                group_by(symbol, !!sym(cn_column)) |>
                summarise(
                    mean_segment_size = mean(segment_length_mb, na.rm = TRUE),
                    n_cytobands_affected = n_distinct(cytoband),
                    .groups = "drop"
                ),
            by = c("symbol", cn_column)
        ) |> 
        ## Calculate frequencies
        mutate(
            group1_freq = group1_altered / group1_total,
            group2_freq = group2_altered / group2_total,
            freq_diff = group2_freq - group1_freq
        ) |> 
        arrange(desc(freq_diff))

    ## Filter genes with sufficient alterations
    freq_data |> 
        dplyr::filter(
            total_altered >= min_altered_samples,
            (group1_altered >= min_altered_samples | group2_altered >= min_altered_samples)
        ) |> 
        ## Perform Fisher's exact test
        rowwise() |>
        mutate(
            fisher_p = {
                contingency_table <- matrix(
                    c(
                        group1_altered, group1_total - group1_altered,
                        group2_altered, group2_total - group2_altered
                    ),
                    nrow = 2,
                    byrow = TRUE
                )
                colnames(contingency_table) <- c("Altered", "Not_Altered")
                rownames(contingency_table) <- c("Group1", "Group2")
                fisher.test(contingency_table)$p.value
            },
            fisher_odds_ratio = {
                contingency_table <- matrix(
                    c(
                        group1_altered, group1_total - group1_altered,
                        group2_altered, group2_total - group2_altered
                    ),
                    nrow = 2,
                    byrow = TRUE
                )
                colnames(contingency_table) <- c("Altered", "Not_Altered")
                rownames(contingency_table) <- c("Group1", "Group2")
                fisher.test(contingency_table)$estimate
            }
        ) |> 
        ungroup() |> 
        ## Adjust for multiple testing
        mutate(
            fisher_q = p.adjust(fisher_p, method = p_adjust_method),
            sig_diff = fisher_q < 0.05,
            enriched_in_group2 = fisher_q < 0.05 & 
                fisher_odds_ratio > 1 & freq_diff > 0,
            depleted_in_group2 = fisher_q < 0.05 & 
                fisher_odds_ratio < 1 & freq_diff < 0
        ) |> 
        arrange(fisher_q, desc(abs(freq_diff))) |> 
        mutate(
            group_comparison = {{group_comparison}},
            .before = "symbol"
        )
}

GetPCGRCNVGroupShareGene <- function(
    data,
    var = "FST.Group",
    cn_column = "cn_status",
    event_type = NULL,
    somatic_matched = TRUE,
    min_freq_per_group = 0.3,
    min_groups_present = 2
) {

    ## Filter by event type
    if (!is.null(event_type)) {
        data <- data |> 
            filter(event_type %in% {{event_type}})
    }

    ## We should use the total cohort size as the denominator for frequency calculations
    ## as we are answer the question "What proportion of samples have this alteration in my cohort?"
    clinical_info <- LoadClinicalInfo()
    
    if (somatic_matched) {

        sample_ids <- clinical_info |> 
            filter(Somatic.Status == "Matched") |> 
            pull(Sample.ID)
        
        total_sample <- length(sample_ids)

    } else {

        sample_ids <- clinical_info |> 
            pull(Sample.ID)

        total_sample <- length(sample_ids)

    }
    
    data <- data |> filter(sample_id %in% sample_ids)

    ## Calculate gene frequencies
    gene_group_freq <- data |>
        distinct(sample_id, symbol, !!sym(cn_column), !!sym(var)) |>
        group_by(symbol, !!sym(var), !!sym(cn_column)) |>
        summarise(
            n_samples_altered = n_distinct(sample_id),
            .groups = "drop"
        ) |>
        ## add total samples per group
        left_join(
            data |>
            group_by(!!sym(var)) |>
            summarise(
                group_total_samples = n_distinct(sample_id),
                .groups = "drop"
            ),
            by = var
        ) |> 
        mutate(
            group_freq = n_samples_altered / group_total_samples,
            cohort_freq = n_samples_altered / total_sample,
            gene_alteration = paste0(symbol, "_", !!sym(cn_column))
        ) |> 
        filter(group_freq >= min_freq_per_group)

    ## identify genes present in multiple groups
    shared_genes <- gene_group_freq |> 
        group_by(gene_alteration, symbol, !!sym(cn_column)) |>
        summarise(
            n_groups_present = n(),
            groups_present_freq = paste(
                paste0(!!sym(var), "(", round(group_freq, 2), ")"), collapse = ", "
            ),
            .groups = "drop"
        ) |> 
        filter(n_groups_present >= min_groups_present) |> 
        arrange(desc(n_groups_present)) 
    
    ## Return
    list(
        shared_genes = shared_genes,
        gene_group_freq = gene_group_freq
    )
}

GetPCGRCNVGroupDiffCytoband <- function(
    data,
    var = "FST.Group",
    cn_column = "cn_status",
    event_type = NULL,
    somatic_matched = TRUE,
    group1 = "U-DFSP",
    group2 = "Pre-FST",
    p_adjust_method = "BH", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")
    min_altered_samples = 3,
    group_comparison = "U-DFSP_vs_Pre-FST"
) {

    ## Filter by event type
    if (!is.null(event_type)) {
        data <- data |> 
            filter(event_type %in% {{event_type}})
    }

    ## We should use the total cohort size as the denominator for frequency 
    ## calculations as we are answer the question "What proportion of samples 
    ## have this alteration in my cohort?"
    
    if (somatic_matched) {

        clinical_info <- LoadClinicalInfo() |> 
            filter(Somatic.Status == "Matched")

    } else {

        clinical_info <- LoadClinicalInfo()

    }
    
    ## Total samples of each group in cohort
    n_sample <- nrow(clinical_info)
    
    n_group1 <- clinical_info |> 
        filter(!!sym(var) %in% c(group1)) |> 
        nrow()
    
    n_group2 <- clinical_info |> 
        filter(!!sym(var) %in% c(group2)) |> 
        nrow()

    ## Calculate alteration frequencies
    freq_data <- data |>
        filter(!!sym(var) %in% c(group1, group2)) |> 
        select(sample_id, cytoband, !!sym(cn_column), !!sym(var)) |>
        distinct() |> 
        group_by(cytoband, !!sym(cn_column)) |> 
        summarise (
            group1_altered = sum(!!sym(var) %in% group1),
            group2_altered = sum(!!sym(var) %in% group2),
            group1_total = n_group1,
            group2_total = n_group2,
            .groups = "drop"
        ) |>
        mutate(total_altered = group1_altered + group2_altered) |> 
        arrange(desc(total_altered)) |>
        ## Calculate frequencies
        mutate(
            group1_freq = group1_altered / group1_total,
            group2_freq = group2_altered / group2_total,
            freq_diff = group2_freq - group1_freq
        ) |> 
        arrange(desc(freq_diff))

    ## Filter genes with sufficient alterations
    freq_data |> 
        dplyr::filter(
            total_altered >= min_altered_samples,
            (group1_altered >= min_altered_samples | group2_altered >= min_altered_samples)
        ) |> 
        ## Perform Fisher's exact test
        rowwise() |>
        mutate(
            fisher_p = {
                contingency_table <- matrix(
                    c(
                        group1_altered, group1_total - group1_altered,
                        group2_altered, group2_total - group2_altered
                    ),
                    nrow = 2,
                    byrow = TRUE
                )
                colnames(contingency_table) <- c("Altered", "Not_Altered")
                rownames(contingency_table) <- c("Group1", "Group2")
                fisher.test(contingency_table)$p.value
            },
            fisher_odds_ratio = {
                contingency_table <- matrix(
                    c(
                        group1_altered, group1_total - group1_altered,
                        group2_altered, group2_total - group2_altered
                    ),
                    nrow = 2,
                    byrow = TRUE
                )
                colnames(contingency_table) <- c("Altered", "Not_Altered")
                rownames(contingency_table) <- c("Group1", "Group2")
                fisher.test(contingency_table)$estimate
            }
        ) |> 
        ungroup() |> 
        ## Adjust for multiple testing
        mutate(
            fisher_q = p.adjust(fisher_p, method = p_adjust_method),
            sig_diff = fisher_q < 0.05,
            enriched_in_group2 = fisher_q < 0.05 & 
                fisher_odds_ratio > 1 & freq_diff > 0,
            depleted_in_group2 = fisher_q < 0.05 & 
                fisher_odds_ratio < 1 & freq_diff < 0
        ) |> 
        arrange(fisher_q, desc(abs(freq_diff))) |> 
        mutate(
            group_comparison = {{group_comparison}},
            .before = "cytoband"
        )
}

GetPCGRCNVGroupShareCytoband <- function(
    data,
    var = "FST.Group",
    cn_column = "cn_status",
    event_type = NULL,
    somatic_matched = TRUE,
    min_freq_per_group = 0.3,
    min_groups_present = 1
) {

    ## Filter by event type
    if (!is.null(event_type)) {
        data <- data |> 
            filter(event_type %in% {{event_type}})
    }

    ## We should use the total cohort size as the denominator for frequency calculations
    ## as we are answer the question "What proportion of samples have this alteration in my cohort?"
    clinical_info <- LoadClinicalInfo()
    
    if (somatic_matched) {

        sample_ids <- clinical_info |> 
            filter(Somatic.Status == "Matched") |> 
            pull(Sample.ID)
        
        total_sample <- length(sample_ids)

    } else {

        sample_ids <- clinical_info |> 
            pull(Sample.ID)

        total_sample <- length(sample_ids)

    }
    
    data <- data |> filter(sample_id %in% sample_ids)

    ## Calculate gene frequencies
    cytoband_group_freq <- data |>
        distinct(sample_id, cytoband, !!sym(cn_column), !!sym(var)) |>
        group_by(cytoband, !!sym(var), !!sym(cn_column)) |>
        summarise(
            n_samples_altered = n_distinct(sample_id),
            .groups = "drop"
        ) |>
        ## add total samples per group
        left_join(
            data |>
            group_by(!!sym(var)) |>
            summarise(
                group_total_samples = n_distinct(sample_id),
                .groups = "drop"
            ),
            by = var
        ) |> 
        mutate(
            group_freq = n_samples_altered / group_total_samples,
            cohort_freq = n_samples_altered / total_sample,
            cytoband_alteration = paste0(cytoband, "_", !!sym(cn_column))
        ) |> 
        filter(group_freq >= min_freq_per_group)

    ## identify genes present in multiple groups
    shared_cytobands <- cytoband_group_freq |> 
        group_by(cytoband_alteration, cytoband, !!sym(cn_column)) |>
        summarise(
            n_groups_present = n(),
            groups_present_freq = paste(
                paste0(!!sym(var), "(", round(group_freq, 2), ")"), collapse = ", "
            ),
            # freq = ,
            .groups = "drop"
        ) |> 
        filter(n_groups_present >= min_groups_present) |> 
        arrange(desc(n_groups_present)) 
    
    ## Return
    list(
        shared_cytobands = shared_cytobands,
        cytoband_group_freq = cytoband_group_freq
    )
}

PCGRComplexOncoplot <- function(
    mat,
    sample_annotation = c("FST.Group", "Metastasis"),
    sample_sorted_by = "FST.Group",
    sample_sorted_level = c(
        "U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"
    ),
    column_title = NULL,
    width = 18,
    height = 12,
    dir = "figures/oncoplot",
    filename = "oncoplot"
) {
    ## Convert the matrix to a tibble
    data <- as_tibble(mat, rownames = "cytoband") |> 
        mutate(across(everything(), as.character)) |> 
        pivot_longer(
            cols = -cytoband, 
            names_to = "Sample.ID", 
            values_to = "alteration"
        ) |> 
        distinct() |> 
        mutate(
            alteration = if_else(
                alteration == "", NA_character_, alteration
            )
        ) |> 
        filter(!is.na(alteration))    

    ## Get matched clinical info
    clinical_info <- LoadClinicalInfo()
    
    sample_idx <- match(
        colnames(mat), 
        clinical_info$Sample.ID
    )
    
    ## Sort samples
    if (!is.null(sample_sorted_by) && !is.null(sample_sorted_level)) {

        # data |> 
        #     left_join(
        #         clinical_info  |> select(Sample.ID, all_of(sample_sorted_by)), 
        #         by = "Sample.ID"
        #     ) |>
        #     group_by(FST.Group)
        #     arrange()

        sample_order <- clinical_info[sample_idx, ] |> 
            # group_by(across(all_of(sample_sorted_by))) |>
            mutate(
                sample_order = factor(
                    !!sym(sample_sorted_by), 
                    levels = sample_sorted_level
                )
            ) |>
            group_by(sample_order) |>
            arrange(!!sym(sample_sorted_by), Sample.ID, .by_group = TRUE) |> 
            select(Sample.ID, all_of(sample_sorted_by), sample_order) |> 
            pull(Sample.ID)

    } else {

        ## Sort by patient
        sample_order <- colnames(mat)
    }

    sort_mat <- mat[, sample_order]

    ## Alteration colors
    data <- as_tibble(sort_mat, rownames = "cytoband") |> 
        mutate(across(everything(), as.character)) |> 
        pivot_longer(
            cols = -cytoband, 
            names_to = "Sample.ID", 
            values_to = "alteration"
        ) |> 
        distinct() |> 
        mutate(
            alteration = if_else(
                alteration == "", NA_character_, alteration
            )
        ) |> 
        filter(!is.na(alteration))
    
    ## Sample annotations
    sort_sample_idx <- match(
        sample_order, 
        clinical_info$Sample.ID
    )

    sample_annotation_df <- clinical_info[sort_sample_idx, ] |> 
        select(Sample.ID, all_of(sample_annotation)) |> 
        column_to_rownames("Sample.ID")

    ## Define colors for the sample annotations
    sample_annotation_colors <- list()

    for (annotation in sample_annotation) {
        sample_annotation_colors[[annotation]] <- study_colors[[annotation]]
    }

    ## Create the matrix annotation colors
    alterations <- data |> pull(alteration) |> unique() |> sort()

    ## Get the matched alteration colors
    alteration_colors <- c()
    
    for (i in alterations) {
        alteration_colors <- c(alteration_colors, study_colors$cna_class[i])
    }

    alter_fun <- list(
        background = alter_graphic("rect", fill = "#CCCCCC"),
        HOMDEL = alter_graphic("rect", fill = alteration_colors[["HOMDEL"]]),
        DEL = alter_graphic("rect", fill = alteration_colors[["DEL"]]),
        GAIN = alter_graphic("rect", fill = alteration_colors[["GAIN"]]),
        AMP = alter_graphic("rect", fill = alteration_colors[["AMP"]]),
        MULTI = alter_graphic("rect", height = 0.33, fill = alteration_colors[["MULTI"]])
    )

    ## Plot the oncoprint
    oncoplot <- oncoPrint(
        sort_mat,
        alter_fun = alter_fun,
        col = alteration_colors,
        
        bottom_annotation = HeatmapAnnotation(
            df = sample_annotation_df,
            col = sample_annotation_colors,
            annotation_height = unit(c(5, 5), "mm"),
            annotation_name_gp = gpar(fontsize = 8)
        ),

        show_pct = TRUE,
        pct_digits = 0,
        pct_side = "right",

        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 6),
        # row_title = "Cytobands (Most Frequent → Least Frequent)",

        column_title = column_title,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        
        column_order = sample_order,
        # Add frequency annotation on the right
        # right_annotation = rowAnnotation(
        #     Frequency = anno_barplot(
        #         cytoband_freq$alteration_frequency[
        #             match(rownames(oncoplot_matrix_sorted), cytoband_freq$matched_cytoband)
        #         ],
        #         bar_width = 1,
        #         gp = gpar(fill = "#666666"),
        #         axis_param = list(side = "bottom", labels_rot = 0)
        #     ),
        #     annotation_width = unit(2, "cm")
        # ),
        
        heatmap_legend_param = list(
            title = "CNV Alteration"
            # labels = c(
            #     "Homozygous Deletion", 
            #     "Heterozygous Deletion", 
            #     "Gain", 
            #     "No Alteration"
            # )
        )
    )

    ## Save the oncoplot
    dir_create(here(dir))
    img_type <- c(".pdf", ".png")
    
    for (img in img_type) {
        
        file <- here(dir, paste0(filename, img))
        
        if (img == ".pdf") {

            # CairoPDF(file, width = width, height = height)
            pdf(file, width = width, height = height)
    
        } else if (img == ".png") {
    
            # CairoPNG(
            #     file, width = width, height = height, res = 300, units = "in",
            #     fonts = "Arial"
            # )
            png(
                file, width = width, height = height, res = 600, units = "in",
                fonts = "Arial"
            )
        }
        
        draw(
            oncoplot,
            heatmap_legend_side = "right",
            annotation_legend_side = "right",
            merge_legends = TRUE,
            padding = unit(c(1, 1, 1, 1), "cm")
        )

        dev.off()

        message(
            paste0("Saved oncoplot: ", file)
        )
    }
    # file <- here(dir, paste0(filename, ".pdf"))
    # pdf(file, width = width, height = height)

    # draw(
    #     oncoplot,
    #     heatmap_legend_side = "right",
    #     annotation_legend_side = "right",
    #     merge_legends = TRUE,
    #     padding = unit(c(1, 1, 1, 1), "cm")
    # )

    # dev.off()

    # message(
    #     paste0("Saved oncoplot: ", file)
    # )

}

GisticComplexOncoplot <- function(
    mat,
    sample_annotation = c("FST.Group", "Metastasis"),
    sample_sorted_by = "FST.Group",
    sample_sorted_level = c(
        "U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"
    ),
    column_title = NULL,
    width = 18,
    height = 12,
    dir = "figures/oncoplot",
    filename = "oncoplot"
) {
    
    ## Get matched clinical info
    clinical_info <- LoadClinicalInfo(is_somatic_matched = FALSE)
    
    sample_idx <- match(
        colnames(mat), 
        clinical_info$Sample.ID
    )
    
    ## Sort samples
    if (!is.null(sample_sorted_by) && !is.null(sample_sorted_level)) {

        # data |> 
        #     left_join(
        #         clinical_info  |> select(Sample.ID, all_of(sample_sorted_by)), 
        #         by = "Sample.ID"
        #     ) |>
        #     group_by(FST.Group)
        #     arrange()

        sample_order <- clinical_info[sample_idx, ] |> 
            # group_by(across(all_of(sample_sorted_by))) |>
            mutate(
                sample_order = factor(
                    !!sym(sample_sorted_by), 
                    levels = sample_sorted_level
                )
            ) |>
            group_by(sample_order) |>
            arrange(!!sym(sample_sorted_by), Sample.ID, .by_group = TRUE) |> 
            select(Sample.ID, all_of(sample_sorted_by), sample_order) |> 
            pull(Sample.ID)

    } else {

        ## Sort by patient
        sample_order <- colnames(mat)
    }

    sort_mat <- mat[, sample_order]
    
    ## Sample annotations
    sort_sample_idx <- match(
        sample_order, 
        clinical_info$Sample.ID
    )

    sample_annotation_df <- clinical_info[sort_sample_idx, ] |> 
        select(Sample.ID, all_of(sample_annotation)) |> 
        column_to_rownames("Sample.ID")

    ## Define colors for the sample annotations
    sample_annotation_colors <- list()

    for (annotation in sample_annotation) {
        sample_annotation_colors[[annotation]] <- study_colors[[annotation]]
    }

    ## Create the matrix annotation colors
    alterations <- c("Amp", "Del")
    # sort_mat |> as_tibble(rownames = "cytoband") |> 
    #     pivot_longer(
    #         cols = -cytoband, 
    #         names_to = "Sample.ID", 
    #         values_to = "alteration"
    #     ) |> 
    #     pull(alteration) |> 
    #     unique()

    ## Get the matched alteration colors
    alteration_colors <- c()
    
    for (i in alterations) {
        alteration_colors <- c(alteration_colors, study_colors$gistic_cn_class[i])
    }

    alter_fun <- list(
        background = alter_graphic("rect", fill = "#CCCCCC"),
        Amp = alter_graphic("rect", fill = alteration_colors[["Amp"]]),
        Del = alter_graphic("rect", fill = alteration_colors[["Del"]])
    )

    ## Plot the oncoprint
    oncoplot <- oncoPrint(
        sort_mat,
        alter_fun = alter_fun,
        col = alteration_colors,
        
        bottom_annotation = HeatmapAnnotation(
            df = sample_annotation_df,
            col = sample_annotation_colors,
            annotation_height = unit(c(5, 5), "mm"),
            annotation_name_gp = gpar(fontsize = 8)
        ),

        show_pct = TRUE,
        pct_digits = 0,
        pct_side = "right",

        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 6),
        # row_title = "Cytobands (Most Frequent → Least Frequent)",

        column_title = column_title,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        
        column_order = sample_order,
        # Add frequency annotation on the right
        # right_annotation = rowAnnotation(
        #     Frequency = anno_barplot(
        #         cytoband_freq$alteration_frequency[
        #             match(rownames(oncoplot_matrix_sorted), cytoband_freq$matched_cytoband)
        #         ],
        #         bar_width = 1,
        #         gp = gpar(fill = "#666666"),
        #         axis_param = list(side = "bottom", labels_rot = 0)
        #     ),
        #     annotation_width = unit(2, "cm")
        # ),
        
        heatmap_legend_param = list(
            title = "CNV Alteration"
            # labels = c(
            #     "Homozygous Deletion", 
            #     "Heterozygous Deletion", 
            #     "Gain", 
            #     "No Alteration"
            # )
        )
    )

    ## Save the oncoplot
    dir_create(here(dir))
    img_type <- c(".pdf", ".png")
    
    for (img in img_type) {
        
        file <- here(dir, paste0(filename, img))
        
        if (img == ".pdf") {

            # CairoPDF(file, width = width, height = height)
            pdf(file, width = width, height = height)
    
        } else if (img == ".png") {
    
            # CairoPNG(
            #     file, width = width, height = height, res = 300, units = "in",
            #     fonts = "Arial"
            # )
            png(
                file, width = width, height = height, res = 600, units = "in",
                fonts = "Arial"
            )
        }
        
        draw(
            oncoplot,
            heatmap_legend_side = "right",
            annotation_legend_side = "right",
            merge_legends = TRUE,
            padding = unit(c(1, 1, 1, 1), "cm")
        )

        dev.off()

        message(
            paste0("Saved oncoplot: ", file)
        )
    }
    # file <- here(dir, paste0(filename, ".pdf"))
    # pdf(file, width = width, height = height)

    # draw(
    #     oncoplot,
    #     heatmap_legend_side = "right",
    #     annotation_legend_side = "right",
    #     merge_legends = TRUE,
    #     padding = unit(c(1, 1, 1, 1), "cm")
    # )

    # dev.off()

    # message(
    #     paste0("Saved oncoplot: ", file)
    # )

}

GisticComplexGroupOncoplot <- function(
    mat,
    is_somatic_matched = TRUE,
    split_by_group = "FST.Group",
    sort_group_level = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"),
    main_group = "FS-DFSP",
    sample_annotation = c("FST.Group", "Metastasis"),
    column_title = NULL,
    width = 18,
    height = 12,
    dir = "figures/oncoplot",
    filename = "oncoplot"
) {
    
    ## Get matched clinical info
    clinical_info <- LoadClinicalInfo(is_somatic_matched = is_somatic_matched)

    ## Prepare the plot data for each heatmap
    plot_data <- list()

    for (i in sort_group_level) {

        sample_ids <- clinical_info |> 
            filter(!!sym(split_by_group) == i) |> 
            pull(Sample.ID)

        ## Matrix
        cur_mat <- mat[, sample_ids, drop = FALSE]
        plot_data[[i]][["mat"]] <- cur_mat

        ## Sample annotation colors
        sample_idx <- match(
            colnames(cur_mat), 
            clinical_info$Sample.ID
        )

        cur_sample_annotation <- clinical_info[sample_idx, ] |> 
            select(Sample.ID, all_of(sample_annotation)) |> 
            column_to_rownames("Sample.ID")
        
        plot_data[[i]][["sample_annotation"]] <- cur_sample_annotation
        
        ## Calcualte the freq data
        # cn_mat[]
    }
    
    ## "---------------------------------------------------------------------"
    ## Sample annotation colors
    ## "---------------------------------------------------------------------"
    sample_annotation_colors <- list()

    for (annotation in sample_annotation) {
        sample_annotation_colors[[annotation]] <- study_colors[[annotation]]
    }

    ## "---------------------------------------------------------------------"
    ## Create the matrix annotation colors
    ## "---------------------------------------------------------------------"
    alterations <- c("Amp", "Del")
    alteration_colors <- c()
    for (i in alterations) {
        alteration_colors <- c(alteration_colors, study_colors$gistic_cn_class[i])
    }

    alter_fun <- list(
        background = alter_graphic("rect", fill = "#CCCCCC"),
        Amp = alter_graphic("rect", fill = alteration_colors[["Amp"]]),
        Del = alter_graphic("rect", fill = alteration_colors[["Del"]])
    )

    ## "---------------------------------------------------------------------"
    ## Generate the heatmap list
    ## "---------------------------------------------------------------------"
    ht_list <- list()
    n_group <- length(sort_group_level)
    main_ht_idx <- n_group

    title_size <- 7
    text_size <- 6

    for (i in seq_along(sort_group_level)) {

        if (i == n_group) {
            
            show_annotation_name <- TRUE
            show_row_names <- TRUE
            show_heatmap_legend <- TRUE

        } else {

            show_annotation_name <- FALSE
            show_row_names <- FALSE
            show_heatmap_legend <- FALSE
        }
        
        bottom_annotation <- HeatmapAnnotation(
                df = plot_data[[i]][["sample_annotation"]],
                col = sample_annotation_colors,
                annotation_height = unit(c(4, 4), "mm"),
                annotation_name_gp = gpar(fontsize = text_size),
                annotation_legend_param = list(
                    labels_gp = gpar(fontsize = text_size),       # Legend text
                    title_gp = gpar(fontsize = title_size)         # Legend title
                ),
                show_annotation_name = show_annotation_name
        )

        ht_list[[i]] <- oncoPrint(
            plot_data[[i]][["mat"]],
            alter_fun = alter_fun,
            col = alteration_colors,
            bottom_annotation = bottom_annotation,
            show_row_names = show_row_names,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = text_size),
            show_column_names = TRUE,
            column_names_gp = gpar(fontsize = text_size),
            column_title_gp = gpar(fontsize = title_size),
            show_heatmap_legend = FALSE,
            heatmap_legend_param = list(
                # labels_gp = gpar(fontsize = text_size), 
                title_gp = gpar(fontsize = title_size),
                title = "CNV Alteration"
            ),
            show_pct = TRUE,
            pct_digits = 0,
            pct_side = "left",
            pct_gp = gpar(fontsize = text_size)
        )
    }

    ht_combined <- Reduce(`+`, ht_list)
    
    heatmap_legend <- Legend(
            labels = names(alteration_colors),
            legend_gp = gpar(fill = alteration_colors),
            title = "CNV Alteration",
            title_gp = gpar(fontsize = title_size),
            labels_gp = gpar(fontsize = text_size),
            grid_height = unit(4, "mm"),
            grid_width = unit(4, "mm")
    )

    ## Save the oncoplot
    dir_create(here(dir))
    img_type <- c(".pdf", ".png")
    
    for (img in img_type) {
        
        file <- here(dir, paste0(filename, img))
        
        if (img == ".pdf") {

            # CairoPDF(file, width = width, height = height)
            pdf(file, width = width, height = height)
    
        } else if (img == ".png") {
    
            # CairoPNG(
            #     file, width = width, height = height, res = 300, units = "in",
            #     fonts = "Arial"
            # )
            png(
                file, width = width, height = height, res = 600, units = "in",
                fonts = "Arial"
            )
        }
        
        draw(
            ht_combined, 
            main_heatmap = main_ht_idx,
            ht_gap = unit(0.2, "cm"),

            ## Heatmap legend
            show_heatmap_legend = TRUE,
            heatmap_legend_side = "bottom",

            ## Annotation legend
            annotation_legend_side = "bottom",
        
            merge_legends = TRUE,
            legend_gap = unit(1, "cm"),
            padding = unit(c(1, 1, 1, 1), "cm"),

            ## Add manual heatmap legend
            heatmap_legend_list = list(heatmap_legend)
        )
        
        dev.off()

        message(
            paste0("Saved oncoplot: ", file)
        )
    }

}

GisticChromPlot <- function(
    data,
    peaks_data = NULL,
    y_lim = NULL,
    x_expand = c(0.02, 0),
    title = NULL,
    x_lab = NULL,
    y_lab = "G-Score",
    is_log_score = FALSE,
    log_base = 10,
    show_horizontal_grid = TRUE,
    show_vertical_grid = TRUE,
    grid_line_type = "dotted",
    grid_line_color = "grey80",
    grid_line_alpha = 1,
    grid_line_width = 0.5,
    peak_label_size = 1.5,
    peak_label_color = "black",
    show_chromosome_ideogram = FALSE,
    ideogram_label_size = 4,
    rel_heights = c(4, 0.3)
) {
    
    ## "--------------------------------------------------------------------"
    ## Get the chromosome information
    ## "--------------------------------------------------------------------"
    chrom_info <- data.frame(
        chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), 
        chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    )

    ## Keep only chromsome 1-22
    chrom_info <- chrom_info[1:22, ] |> as_tibble()

    chrom_info$chromName <- str_remove(chrom_info$chromName, "chr")

    ## make all chromosome coordinate are start from 0
    chrom_info$chromlengthCumsum <- cumsum(as.numeric(chrom_info$chromlength))
    chrom_info$chormStartPosFrom0 <- c(0, chrom_info$chromlengthCumsum[-nrow(chrom_info)])

    ## Calculate the middle coordinate for each chromosome
    tmp_middle <- diff(c(0, chrom_info$chromlengthCumsum)) / 2
    chrom_info$chromMidelePosFrom0 <- chrom_info$chormStartPosFrom0 + tmp_middle

    ## Use the start of coordinate as the Amp or Del
    chromID <- data$Chromosome
    data$StartPos <- data$Start + chrom_info$chormStartPosFrom0[chromID]
    data$EndPos <- data$End + chrom_info$chormStartPosFrom0[chromID]

    pseudocount <- 0.01

    ## Apply log transformation
    if (is_log_score) {
        
        data <- data |> 
            mutate(
                G.score = sign(G.score) * 
                    log10(abs(G.score) + pseudocount), base = log_base
            )
    }

    ## Convert the to Del to be negative
    data[data$Type == "Del", "G.score"] <- data[data$Type == "Del", "G.score"] * -1

    ## "--------------------------------------------------------------------"
    ## Process peaks data if provided
    ## "--------------------------------------------------------------------"
    if (!is.null(peaks_data) && nrow(peaks_data) > 0) {
        
        # Add positional information to peaks using the same coordinate system
        peaks_data$StartPos <- peaks_data$Start + chrom_info$chormStartPosFrom0[peaks_data$Chromosome]
        peaks_data$EndPos <- peaks_data$End + chrom_info$chormStartPosFrom0[peaks_data$Chromosome]
        peaks_data$MidPos <- (peaks_data$StartPos + peaks_data$EndPos) / 2
        
        # Find the G-score at peak positions by matching to nearest data point
        peaks_data$G_score_at_peak <- sapply(
            peaks_data$MidPos, function(pos) {
                closest_idx <- which.min(abs(data$StartPos - pos))
                data$G.score[closest_idx]
            }
        )

        # Position labels above/below peaks with some offset
        peaks_data <- peaks_data |>
            mutate(
                label_y = case_when(
                    Type == "Amp" ~ pmax(G_score_at_peak + 0.15, 0.1),  # Above positive peaks
                    Type == "Del" ~ pmin(G_score_at_peak - 0.15, -0.1), # Below negative peaks
                    TRUE ~ G_score_at_peak + 0.15
                )
            )
    }

    ## Determine y-axis limits if not provided
    if (is.null(y_lim)) {
        y_range <- range(data$G.score, na.rm = TRUE)
        y_padding <- diff(y_range) * 0.1
        y_lim <- c(y_range[1] - y_padding, y_range[2] + y_padding)

        ## Use the integer for y lim
        y_lim <- round(y_lim)     
    }

    ## Determine y-axis limits if not provided
    if (is.null(y_lim)) {
        y_range <- range(data$G.score, na.rm = TRUE)
        y_padding <- diff(y_range) * 0.1
        y_lim <- c(y_range[1] - y_padding, y_range[2] + y_padding)

        # Expand y_lim to accommodate peak labels if peaks_data exists
        if (!is.null(peaks_data) && nrow(peaks_data) > 0) {
            label_range <- range(peaks_data$label_y, na.rm = TRUE)
            y_lim[1] <- min(y_lim[1], label_range[1] - 0.1)
            y_lim[2] <- max(y_lim[2], label_range[2] + 0.1)
        }

        # Round appropriately based on log transformation
        if (is_log_score) {
            y_lim <- round(y_lim, 1)  # More precision for log scale
        } else {
            y_lim <- round(y_lim)     # Integer for original scale
        }

        # ## Use the integer for y lim
        y_lim <- round(y_lim)
        
    }

    y_breaks <- pretty(y_lim, n = y_lim[2] - y_lim[1])
    x_lim <- c(0, max(chrom_info$chromlengthCumsum))
    
    ## Create the base plot
    plot <- ggplot(data, aes(StartPos, G.score)) +
        geom_area(
            aes(group = Type, fill = factor(Type, levels = c("Del", "Amp"))),
            alpha = 0.8
        ) +
        scale_fill_lancet(
            guide = guide_legend(reverse = TRUE), name = "Type"
        )

    ## Add the horizontal grid lines
    if (show_horizontal_grid) {
        
        # Remove the horizontal line at y=0
        h_y_breaks <- y_breaks[y_breaks != 0]

        plot <- plot +
            geom_hline(
                yintercept = h_y_breaks, 
                linetype = grid_line_type, 
                color = grid_line_color, 
                alpha = grid_line_alpha,
                linewidth = grid_line_width
            )
    }

    ## Add vertical chromosome separator lines
    if (show_vertical_grid) {

        # Create vertical lines at chromosome boundaries
        ## excluding the first position and last
        chromosome_boundaries <- c(
            0,
            # chrom_info$chromlengthCumsum[-length(chrom_info$chromlengthCumsum)]
            chrom_info$chromlengthCumsum
        )
        
        plot <- plot +
            geom_vline(
                xintercept = chromosome_boundaries,
                linetype = grid_line_type, 
                color = grid_line_color, 
                alpha = grid_line_alpha,
                linewidth = grid_line_width
            )  
    }
    
    ## "--------------------------------------------------------------------"
    ## Add peak labels if peaks_data is provided
    ## "--------------------------------------------------------------------"
    if (!is.null(peaks_data) && nrow(peaks_data) > 0) {
        plot <- plot +
            # Add peak points
            # geom_point(
            #     data = peaks_data,
            #     aes(x = MidPos, y = G_score_at_peak),
            #     color = peak_label_color,
            #     size = 1.5,
            #     inherit.aes = FALSE
            # ) +
            # # Add connecting lines from peak to label
            # geom_segment(
            #     data = peaks_data,
            #     aes(
            #         x = MidPos, xend = MidPos,
            #         y = G_score_at_peak, yend = label_y
            #     ),
            #     color = peak_label_color,
            #     alpha = 0.6,
            #     size = 0.3,
            #     inherit.aes = FALSE
            # ) +
            # Add peak labels using ggrepel for better positioning
            ggrepel::geom_text_repel(
                data = peaks_data,
                aes(
                    x = MidPos, 
                    y = label_y, 
                    label = peak_label
                ),
                size = peak_label_size,
                color = peak_label_color,
                fontface = "plain",
                family = "Arial",
                # Repel parameters to avoid overlap
                min.segment.length = 0, # always show segments
                segment.color = peak_label_color,
                segment.alpha = 0.6,    # make segments more visible
                segment.size = 0.3,     # thicker segments
                segment.linetype = "solid", # solid line segments
                
                # Force labels further away to create longer segments
                box.padding = 0.3,   # creates more space around text boxes, forcing longer segments
                point.padding = 0.3, # creates more space around the original points
                force = 5,           # stronger repelling force pushes labels further away
                force_pull = 3,      # Added parameter that pulls labels toward empty areas
                max.overlaps = Inf,

                # Direction control to make segmetns longer
                direction = "both",   # # Allow both x and y movement
                nudge_y = ifelse(peaks_data$Type == "Del", -0.3, 0.3), # Pushes labels further up (for Amp) or down (for Del) from their original positions
                
                # Keep labels within plot area
                xlim = c(x_lim[1], x_lim[2]),
                ylim = c(y_lim[1], y_lim[2]),
                inherit.aes = FALSE
            )
    }

    ## Main plot styling
    plot <- plot +
        scale_x_continuous(
            expand = x_expand, 
            limits = x_lim, 
            name = NULL, 
            labels = NULL,
            breaks = NULL
        ) +
        scale_y_continuous(
            limits = y_lim,
            breaks = y_breaks,
            labels = scales::number_format(accuracy = 0.1),
            expand = c(0,0),  # Small padding
            oob = scales::oob_keep  # Keep out-of-bounds values instead of removing them
        ) +
        labs(y = y_lab, title = title) +
        theme_classic() +
        plot_theme() +
        theme(
            legend.position = "right",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank()
        )

    ## "--------------------------------------------------------------------"
    ## Add chromosome ideogram if needed
    ## "--------------------------------------------------------------------"
    if (show_chromosome_ideogram) {
        
        chrom_info$fill_color <- rep(
            c("black", "white"), 
            length.out = nrow(chrom_info)
        )
        chrom_info$label_color <- ifelse(
            chrom_info$fill_color == "black", "white", "black"
        )
        
        ideogram_data <- chrom_info |>
            mutate(
                xmin = chormStartPosFrom0,
                xmax = chromlengthCumsum,
                ymin = 0,
                ymax = 1
            )
        
        ideogram_plot <- ggplot() +
            geom_rect(
                data = ideogram_data,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = ideogram_data$fill_color,
                color = "black",
                size = 0.2
            ) +
            geom_text(
                data = chrom_info,
                aes(
                    x = chromMidelePosFrom0, 
                    y = 0.5, 
                    label = chromName, 
                    color = label_color
                ),
                size.unit = "pt",
                size = ideogram_label_size,
                inherit.aes = FALSE,
                show.legend = FALSE
            ) +
            scale_color_identity() +
            scale_x_continuous(
                expand = x_expand, 
                limits = x_lim,
                name = x_lab,
                labels = NULL,
                breaks = NULL
            ) +
            scale_y_continuous(
                expand = c(0, 0),
                limits = c(0, 1)
            ) +
            theme_void() +
            plot_theme() +
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                plot.margin = margin(t = 0, r = 10, b = 10, l = 10)
            )

        plot <- wrap_plots(
            plot + 
                theme(
                    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm"),
                    panel.spacing = unit(0, "lines")
                ), 
            ideogram_plot,
            ncol = 1,
            heights =  rel_heights
        )

    }

    ## Return
    plot
}

LoadGisticScore <- function(
    gistic_dir = "data/wes/GISTIC2/somatic_matched",
    group = "U-DFSP"
) {
    
    file <- here(gistic_dir, group, "scores.gistic")
    
    data <- read.table(
        file = file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    ) |> as_tibble()

    return(data)
}

GetPCGRGisticOverlap <- function(
    peak_label,
    pcgr_cna_data, 
    gistic_score_data, 
    min_overlap = 0.5
) {
    
    ## peak labels pcgr cytoband data
    
    ## Get the coordinate for each peak labels
    peaks_data <- list()
    
    for (i in peak_label) {
        
        message(paste(" - Processing peak:", i))

        data <- pcgr_cna_data |> 
            filter(cytoband_alteration %in% i) |> 
            select(
                cytoband_alteration, chrom, start, end, var_id, "cytoband", "cn_status"
            ) |>
            distinct()
    
        peak_type <- case_when(
            str_detect(i, "DEL|HOMDEL") ~ "Del",
            str_detect(i, "AMP|GAIN") ~ "Amp",
            TRUE ~ "Other"
        )
    
        ## Create GRanges objects
        pcgr_ranges <- GRanges(
            seqnames = data$chrom,
            ranges = IRanges(
                start = data$start,
                end = data$end
            ),
            mcols = data[, c("var_id", "cytoband", "cn_status")]
        )
    
        gistic_ranges <- GRanges(
            seqnames = paste0("chr", gistic_score$Chromosome),
            ranges = IRanges(
                start = gistic_score$Start,
                end = gistic_score$End
            ),
            mcols = gistic_score[, c("Type", "G.score", "X.log10.q.value.", "frequency")]
        )

        ## Find overlaps
        overlaps <- findOverlaps(pcgr_ranges, gistic_ranges)
        
        if (length(overlaps) == 0) {
            warning("No overlaps found between PCGR data and GISTIC scores")
            return(tibble())
        }

        ## Calculate overlap fractions
        overlap_widths <- width(
            pintersect(
                pcgr_ranges[queryHits(overlaps)],
                gistic_ranges[subjectHits(overlaps)]
            )
        )

        pcgr_widths <- width(pcgr_ranges[queryHits(overlaps)])
        
        overlap_fractions <- overlap_widths / pcgr_widths
        
        ## Filter by minimum overlap
        valid_overlaps <- overlaps[overlap_fractions >= min_overlap]
        
        if (length(valid_overlaps) == 0) {
            warning("No overlaps found meeting the minimum overlap threshold")
            return(tibble())
        }

        ## Create overlapped data and summarize by cytoband
        peaks_data[[i]] <- data[queryHits(valid_overlaps), ] |>
            bind_cols(
                gistic_score_data[
                    subjectHits(valid_overlaps), 
                    # c("Type", "G.score", "X.log10.q.value.", "frequency")
                ]
            ) |>
            mutate(
                overlap_fraction = overlap_fractions[overlap_fractions >= min_overlap]
            ) |>
            filter(Type %in% peak_type) |> 
            group_by(cytoband, cn_status, cytoband_alteration, Chromosome, Type) |> 
            summarise(
                ## Use mean positions for peak coordinates
                Start = as.integer(mean(Start)),
                End = as.integer(mean(End)),
                mean_gistic_score = mean(G.score, na.rm = TRUE),
                max_gistic_score = max(G.score, na.rm = TRUE),
                mean_overlap = mean(overlap_fraction),
                .groups = "drop"
            )
    }

    list_rbind(peaks_data)
}

GetPCGRGisticOverlap2 <- function(
    peak_label,
    pcgr_cna_data, 
    gistic_score_data, 
    min_overlap = 0.5,
    overlap_method = "any"  # "any", "pcgr_based", "gistic_based", "reciprocal"
) {
    
    ## Get the coordinate for each peak labels
    peaks_data <- list()
    successful_peaks <- character()
    failed_peaks <- character()
    
    for (i in peak_label) {
        
        message(paste(" - Processing peak:", i))
        
        data <- pcgr_cna_data |> 
            filter(cytoband_alteration %in% i) |> 
            select(
                cytoband_alteration, chrom, start, end, var_id, "cytoband", "cn_status"
            ) |>
            distinct()
        
        if (nrow(data) == 0) {
            warning("No PCGR data found for peak: ", i)
            failed_peaks <- c(failed_peaks, i)
            next
        }
    
        peak_type <- case_when(
            str_detect(i, "DEL|HOMDEL") ~ "Del",
            str_detect(i, "AMP|GAIN") ~ "Amp",
            TRUE ~ "Other"
        )
    
        ## Create GRanges objects
        pcgr_ranges <- GRanges(
            seqnames = data$chrom,
            ranges = IRanges(
                start = data$start,
                end = data$end
            ),
            mcols = data[, c("var_id", "cytoband", "cn_status")]
        )
    
        gistic_ranges <- GRanges(
            seqnames = paste0("chr", gistic_score_data$Chromosome),
            ranges = IRanges(
                start = gistic_score_data$Start,
                end = gistic_score_data$End
            ),
            mcols = gistic_score_data[, c("Type", "G.score", "X.log10.q.value.", "frequency")]
        )

        ## Find overlaps - use different strategies based on overlap_method
        if (overlap_method == "any") {
            overlaps <- findOverlaps(pcgr_ranges, gistic_ranges, minoverlap = 1L)
        } else {
            overlaps <- findOverlaps(pcgr_ranges, gistic_ranges)
        }
        
        if (length(overlaps) == 0) {
            warning("No overlaps found between PCGR data and GISTIC scores for peak: ", i)
            failed_peaks <- c(failed_peaks, i)
            next
        }

        ## Calculate overlap fractions only if needed
        if (overlap_method != "any") {
            overlap_widths <- width(
                pintersect(
                    pcgr_ranges[queryHits(overlaps)],
                    gistic_ranges[subjectHits(overlaps)]
                )
            )
            
            pcgr_widths <- width(pcgr_ranges[queryHits(overlaps)])
            gistic_widths <- width(gistic_ranges[subjectHits(overlaps)])
            
            ## Choose overlap method
            if (overlap_method == "pcgr_based") {
                overlap_fractions <- overlap_widths / pcgr_widths
            } else if (overlap_method == "gistic_based") {
                overlap_fractions <- overlap_widths / gistic_widths
            } else if (overlap_method == "reciprocal") {
                overlap_fractions <- pmin(
                    overlap_widths / pcgr_widths,
                    overlap_widths / gistic_widths
                )
            }
            
            ## Filter by minimum overlap
            valid_overlaps <- overlaps[overlap_fractions >= min_overlap]
            valid_overlap_fractions <- overlap_fractions[overlap_fractions >= min_overlap]
        } else {
            valid_overlaps <- overlaps
            valid_overlap_fractions <- rep(1, length(overlaps))
        }
        
        if (length(valid_overlaps) == 0) {
            warning("No overlaps found meeting the minimum overlap threshold for peak: ", i)
            failed_peaks <- c(failed_peaks, i)
            next
        }

        ## Filter GISTIC data for the correct type (Amp/Del)
        gistic_subset <- gistic_score_data[subjectHits(valid_overlaps), ] |>
            filter(Type == peak_type)
        
        if (nrow(gistic_subset) == 0) {
            warning("No GISTIC scores of type '", peak_type, "' found for peak: ", i)
            failed_peaks <- c(failed_peaks, i)
            next
        }
        
        ## Get corresponding PCGR data indices after filtering
        pcgr_subset_indices <- queryHits(valid_overlaps)[which(gistic_score_data[subjectHits(valid_overlaps), ]$Type == peak_type)]
        gistic_subset_indices <- subjectHits(valid_overlaps)[which(gistic_score_data[subjectHits(valid_overlaps), ]$Type == peak_type)]
        
        ## Create overlapped data and summarize
        peak_result <- data[pcgr_subset_indices, ] |>
            bind_cols(
                gistic_score_data[gistic_subset_indices, ] |>
                    select(Chromosome, Start, End, Type, G.score, X.log10.q.value., frequency)
            ) |>
            mutate(
                overlap_fraction = valid_overlap_fractions[which(gistic_score_data[subjectHits(valid_overlaps), ]$Type == peak_type)]
            ) |>
            group_by(cytoband, cn_status, cytoband_alteration, Chromosome, Type) |> 
            summarise(
                Start = as.integer(round(weighted.mean(Start, overlap_fraction, na.rm = TRUE))),
                End = as.integer(round(weighted.mean(End, overlap_fraction, na.rm = TRUE))),
                mean_gistic_score = if(all(is.na(G.score))) NA_real_ else weighted.mean(G.score, overlap_fraction, na.rm = TRUE),
                max_gistic_score = if(all(is.na(G.score))) NA_real_ else max(G.score, na.rm = TRUE),
                mean_overlap = mean(overlap_fraction, na.rm = TRUE),
                n_overlapping_intervals = n(),
                .groups = "drop"
            ) |>
            filter(!is.na(mean_gistic_score))  # Remove peaks with no valid GISTIC scores

        if (nrow(peak_result) > 0) {
            peaks_data[[i]] <- peak_result
            successful_peaks <- c(successful_peaks, i)
        } else {
            failed_peaks <- c(failed_peaks, i)
        }
    }

    ## Print summary
    message(paste(" - Successfully processed", length(successful_peaks), "peaks"))
    message(paste(" - Failed to process", length(failed_peaks), "peaks"))
    
    if (length(failed_peaks) > 0) {
        message(" - Failed peaks: ", paste(failed_peaks[1:min(5, length(failed_peaks))], collapse = ", "))
        if (length(failed_peaks) > 5) {
            message(" - ... and", length(failed_peaks) - 5, "more")
        }
    }

    if (length(peaks_data) == 0) {
        warning("No valid peak data found for any peaks")
        return(tibble())
    }

    result <- list_rbind(peaks_data)
    return(result)
}

AddPeaksPCGRCountData <- function(
    peaks_data,
    pcgr_cna_data,
    var,
    cn_column = "cn_status",
    event_type = NULL,
    group,
    is_somatic_matched = TRUE
) {
    ## Filter by event type
    if (!is.null(event_type)) {
        data <- pcgr_cna_data |> 
            filter(event_type %in% {{event_type}})
    } else {

        data <- pcgr_cna_data
    }

    ## We should use the total cohort size as the denominator for frequency 
    ## calculations as we are answer the question "What proportion of samples 
    ## have this alteration in my cohort?"
    
    if (is_somatic_matched) {

        clinical_info <- LoadClinicalInfo() |> 
            filter(Somatic.Status == "Matched")

    } else {

        clinical_info <- LoadClinicalInfo()

    }
    
    ## Total samples of each group in cohort
    n_sample <- nrow(clinical_info)

    n_sample_group <- clinical_info |> 
        filter(!!sym(var) %in% group) |> 
        nrow()

    count_data <- data |> 
        filter(sample_id %in% clinical_info$Sample.ID) |> 
        filter(cytoband_alteration %in% peaks_data$cytoband_alteration) |> 
        filter(!!sym(var) %in% group) |>
        group_by(cytoband_alteration, cytoband, !!sym(cn_column)) |>
        summarise(
            group_altered = n_distinct(sample_id),
            group_total = n_sample_group,
            freq = group_altered / group_total,
            .groups = "drop"
        ) |> 
        mutate(
            group = {{group}},
            peak_label = cytoband
            # peak_label = sapply(
            #     str_split(cytoband, ":"),
            #     function(x) x[2]
            # )
        ) |> 
        mutate(
            peak_label = paste0(
                peak_label, "(", group_altered, ")"
            )
        )

    ## Merge the count data
    peaks_data |> 
        left_join(
            count_data,
            by = c("cytoband", "cn_status", "cytoband_alteration")
        ) |> 
        filter(!is.na(group_altered))
}

LoadPCGRCytobandStatData <- function(
    is_somatic_matched = TRUE
) {

    dir <- "data/processed"

    if (is_somatic_matched) {
        
        filename <- "wes_pcgr_DFSP_cohort_somatic_matched_cytoband_stats"

    } else {
        
        filename <- "wes_pcgr_DFSP_cohort_all_cytoband_stats"
    }

    LoadData(
        filename = filename,
        dir = dir
    )
}

GetGisticCNData <- function(
    gistic_dir, 
    group = "all_tumors", 
    gistic_qval_thres = 0.5
) {

    gistic <- LoadGisticData(
        gistic_dir = gistic_dir,
        group = group,
        gistic_qval_thres = gistic_qval_thres
    )

    ## Get alternation matrix data
    cn_mat <- gistic@cnMatrix
    # num_mat = gistic@numericMatrix
    # dim(cn_mat)
    # cn_mat[1:5, 1:5]

    cn_cytoband_tbl <- cn_mat |> 
        as_tibble(rownames = "Unique_Name")

    cytoband_samples <- colnames(cn_cytoband_tbl)[-1]

    ## Get the gene-level data
    cn_gene_data <- gistic@data |> 
        as_tibble() |> 
        clean_names() |> 
        dplyr::rename(
            symbol = hugo_symbol,
            sample_id = tumor_sample_barcode,
            variant_class = variant_classification,
        )

    gene_samples <- cn_gene_data |> pull(sample_id) |> unique()

    missed_sample <- setdiff(cytoband_samples, gene_samples)

    cn_gene_tbl <- cn_gene_data |>
        group_by(symbol, sample_id) |>
        summarise(
            variant_class = paste(unique(variant_class), collapse = ","),
            .groups = "drop"
        ) |>
        pivot_wider(
            names_from = sample_id,
            values_from = variant_class,
            values_fill = ""
        ) |> 
        mutate(across(everything(), as.character))

    if (length(missed_sample > 0)) {

        cn_gene_tbl <- cn_gene_tbl |> 
            add_column(!!sym(missed_sample) := "")
    }

    list(
        cn_cytoband_tbl = cn_cytoband_tbl,
        cn_gene_tbl = cn_gene_tbl,
        cn_gene_data = cn_gene_data
    )
}

GetGisticGeneData <- function(
    gistic_dir, 
    group = "all_tumors", 
    gistic_qval_thres = 0.5
) {

    gistic <- LoadGisticData(
        gistic_dir = gistic_dir,
        group = group,
        gistic_qval_thres = gistic_qval_thres
    )

    gistic@gene.summary
    gistic@data

    ## Get alternation matrix data
    cn_mat <- gistic@cnMatrix
    # num_mat = gistic@numericMatrix
    # dim(cn_mat)
    # cn_mat[1:5, 1:5]

    cn_data <- cn_mat |> as_tibble(rownames = "Unique_Name")
    
    cn_data
}

GetGisticCytobandGroupStatRes <- function(
    cn_cytoband_tbl,
    group1_samples,
    group2_samples,
    group_comparison,
    alternative = "two.sided",
    p_adjust_method = "BH"
) {

    ## Get the matched group info
    group1_total <- length(group1_samples)
    group2_total <- length(group2_samples)

    # ## Load the Gistic2 data that with loose threshold using all samples as it 
    # ## contains more peaks data
    # gistic <- LoadGisticData(
    #     gistic_dir = gistic_dir,
    #     group = "all_tumors",
    #     gistic_qval_thres = 0.5
    # )

    # ## Cytoband summary data
    # cytoband_summary <- gistic@cytoband.summary |> as_tibble()

    # cn_data <- GetGisticCNData(
    #     gistic_dir = gistic_dir,
    #     group = "all_tumors",
    #     gistic_qval_thres = 0.5
    # ) |>
    #     pivot_longer(
    #         cols = -Unique_Name,
    #         names_to = "Sample.ID",
    #         values_to = "cn_status"
    #     )

    cn_data <- cn_cytoband_tbl |> 
        pivot_longer(
            cols = -Unique_Name,
            names_to = "Sample.ID",
            values_to = "cn_status"
        )

    ## Calculate the alteration frequencies
    group1_freq_data <- cn_data |> 
        filter(Sample.ID %in% group1_samples) |> 
        group_by(Unique_Name) |>
        summarise(
            group1_altered = sum(cn_status !=""),
            .groups = "drop"
        ) |> 
        mutate(group1_total = group1_total)

    group2_freq_data <- cn_data |> 
        filter(Sample.ID %in% group2_samples) |> 
        group_by(Unique_Name) |>
        summarise(
            group2_altered = sum(cn_status !=""),
            .groups = "drop"
        ) |> 
        mutate(group2_total = group2_total)

    freq_data <- group1_freq_data |> 
        left_join(group2_freq_data, by = "Unique_Name") |> 
        mutate(
            group1_freq = group1_altered / group1_total,
            group2_freq = group2_altered / group2_total,
            freq_diff = group2_freq - group1_freq
        )

    ## Do the fisher test
    all_cytobands <- unique(freq_data$Unique_Name)

    fisher_results <- map_dfr(
        all_cytobands,
        function(cytoband) {

            cur_data <- freq_data |> filter(Unique_Name == cytoband)

            group1_altered <- cur_data |> pull(group1_altered)

            group2_altered <- cur_data |> pull(group2_altered)

            group1_total <- cur_data |> pull(group1_total)

            group2_total <- cur_data |> pull(group2_total)

            ## Fisher's exact test
            contingency_table <- matrix(
                    c(
                        group1_altered, group1_total - group1_altered,
                        group2_altered, group2_total - group2_altered
                    ), 
                    byrow = TRUE,
                    nrow = 2
            )
            colnames(contingency_table) <- c("Altered", "Not_Altered")
            rownames(contingency_table) <- c("Group1", "Group2")

            fisher_test <- fisher.test(
                contingency_table,
                alternative = alternative
            )

            tibble(
                Unique_Name = cytoband,
                fisher_p_value = fisher_test$p.value,
                # fisher_odds_ratio = fisher_test$estimate
                # conf_int_low = fisher_test$conf.int[1],
                # conf_int_high = fisher_test$conf.int[2]
            )
        }
    ) |> 
    ## Adjust for multiple testing
    mutate(
        fisher_p_adj = p.adjust(
            fisher_p_value, 
            method = p_adjust_method
        )
    ) |> 
    arrange(fisher_p_value)
    
    ## Merge the data
    freq_data |> 
        left_join(fisher_results, by = "Unique_Name") |> 
        # left_join(cytoband_summary, by = "Unique_Name") |> 
        mutate(group_comparison = group_comparison, .after = Unique_Name)
}

GetGisticCNVGeneGroupStatRes <- function(
    cn_gene_tbl,
    group1_samples,
    group2_samples,
    group_comparison,
    alternative = "two.sided",
    p_adjust_method = "BH"
) {

    ## Get the matched group info
    group1_total <- length(group1_samples)
    group2_total <- length(group2_samples)

    # ## Load the Gistic2 data that with loose threshold using all samples as it 
    # ## contains more peaks data
    # gistic <- LoadGisticData(
    #     gistic_dir = gistic_dir,
    #     group = "all_tumors",
    #     gistic_qval_thres = 0.5
    # )

    # ## Cytoband summary data
    # cytoband_summary <- gistic@cytoband.summary |> as_tibble()

    # cn_data <- GetGisticCNData(
    #     gistic_dir = gistic_dir,
    #     group = "all_tumors",
    #     gistic_qval_thres = 0.5
    # ) |>
    #     pivot_longer(
    #         cols = -Unique_Name,
    #         names_to = "Sample.ID",
    #         values_to = "cn_status"
    #     )

    cn_data <- cn_gene_tbl |> 
        pivot_longer(
            cols = -symbol,
            names_to = "sample_id",
            values_to = "cn_status"
        )

    ## Calculate the alteration frequencies
    group1_freq_data <- cn_data |> 
        filter(sample_id %in% group1_samples) |> 
        group_by(symbol) |>
        summarise(
            group1_altered = sum(cn_status !=""),
            .groups = "drop"
        ) |> 
        mutate(group1_total = group1_total)

    group2_freq_data <- cn_data |> 
        filter(sample_id %in% group2_samples) |> 
        group_by(symbol) |>
        summarise(
            group2_altered = sum(cn_status !=""),
            .groups = "drop"
        ) |> 
        mutate(group2_total = group2_total)

    freq_data <- group1_freq_data |> 
        left_join(group2_freq_data, by = "symbol") |> 
        mutate(
            group1_freq = group1_altered / group1_total,
            group2_freq = group2_altered / group2_total,
            freq_diff = group2_freq - group1_freq
        )

    ## Do the fisher test
    all_symbols <- unique(freq_data$symbol)

    fisher_results <- map_dfr(
        all_symbols,
        function(gene) {

            cur_data <- freq_data |> filter(symbol == gene)

            group1_altered <- cur_data |> pull(group1_altered)

            group2_altered <- cur_data |> pull(group2_altered)

            group1_total <- cur_data |> pull(group1_total)

            group2_total <- cur_data |> pull(group2_total)

            ## Fisher's exact test
            contingency_table <- matrix(
                    c(
                        group1_altered, group1_total - group1_altered,
                        group2_altered, group2_total - group2_altered
                    ), 
                    byrow = TRUE,
                    nrow = 2
            )
            colnames(contingency_table) <- c("Altered", "Not_Altered")
            rownames(contingency_table) <- c("Group1", "Group2")

            fisher_test <- fisher.test(
                contingency_table,
                alternative = alternative
            )

            tibble(
                symbol = gene,
                fisher_p_value = fisher_test$p.value,
                # fisher_odds_ratio = fisher_test$estimate
                # conf_int_low = fisher_test$conf.int[1],
                # conf_int_high = fisher_test$conf.int[2]
            )
        }
    ) |> 
    ## Adjust for multiple testing
    mutate(
        fisher_p_adj = p.adjust(
            fisher_p_value, 
            method = p_adjust_method
        )
    ) |> 
    arrange(fisher_p_value)
    
    ## Merge the data
    freq_data |> 
        left_join(fisher_results, by = "symbol") |> 
        # left_join(cytoband_summary, by = "Unique_Name") |> 
        mutate(group_comparison = group_comparison, .after = symbol)
}

GetPCGRSNVIncelGroupStatRes <- function(
    snv_indel_tbl,
    group1_samples,
    group2_samples,
    group_comparison,
    alternative = "two.sided",
    p_adjust_method = "BH"
) {

    ## Get the matched group info
    group1_total <- length(group1_samples)
    group2_total <- length(group2_samples)

    ## Calculate the alteration frequencies
    group1_freq_data <- snv_indel_tbl |> 
        pivot_longer(
            cols = -symbol,
            names_to = "sample_id",
            values_to = "variant_class"
        ) |> 
        filter(sample_id %in% group1_samples) |> 
        group_by(symbol) |>
        summarise(
            group1_altered = sum(variant_class !=""),
            .groups = "drop"
        ) |> 
        mutate(group1_total = group1_total)

    group2_freq_data <- snv_indel_tbl |> 
        pivot_longer(
            cols = -symbol,
            names_to = "sample_id",
            values_to = "variant_class"
        ) |> 
        filter(sample_id %in% group2_samples) |> 
        group_by(symbol) |>
        summarise(
            group2_altered = sum(variant_class !=""),
            .groups = "drop"
        ) |> 
        mutate(group2_total = group2_total)

    freq_data <- group1_freq_data |> 
        left_join(group2_freq_data, by = "symbol") |> 
        mutate(
            group1_freq = group1_altered / group1_total,
            group2_freq = group2_altered / group2_total,
            freq_diff = group2_freq - group1_freq
        )

    ## Do the fisher test
    all_genes <- unique(freq_data$symbol)

    fisher_results <- map_dfr(
        all_genes,
        function(gene) {

            cur_data <- freq_data |> filter(symbol == gene)

            group1_altered <- cur_data |> pull(group1_altered)

            group2_altered <- cur_data |> pull(group2_altered)

            group1_total <- cur_data |> pull(group1_total)

            group2_total <- cur_data |> pull(group2_total)

            ## Fisher's exact test
            contingency_table <- matrix(
                    c(
                        group1_altered, group1_total - group1_altered,
                        group2_altered, group2_total - group2_altered
                    ), 
                    byrow = TRUE,
                    nrow = 2
            )
            colnames(contingency_table) <- c("Altered", "Not_Altered")
            rownames(contingency_table) <- c("Group1", "Group2")

            fisher_test <- fisher.test(
                contingency_table,
                alternative = alternative
            )

            tibble(
                symbol = gene,
                fisher_p_value = fisher_test$p.value,
                # fisher_odds_ratio = fisher_test$estimate
                # conf_int_low = fisher_test$conf.int[1],
                # conf_int_high = fisher_test$conf.int[2]
            )
        }
    ) |> 
    ## Adjust for multiple testing
    mutate(
        fisher_p_adj = p.adjust(
            fisher_p_value, 
            method = p_adjust_method
        )
    ) |> 
    arrange(fisher_p_value)
    
    ## Merge the data
    freq_data |> 
        left_join(fisher_results, by = "symbol") |> 
        mutate(group_comparison = group_comparison, .after = symbol)
}

GetGroupStatRes <- function(
    mat,
    group1_samples,
    group2_samples,
    variant_type = "symbol",
    group_comparison,
    alternative = "two.sided",
    p_adjust_method = "BH"
) {

    ## Get the matched group info
    group1_total <- length(group1_samples)
    group2_total <- length(group2_samples)

    stat_data <- as_tibble(mat, rownames = variant_type) |>
        pivot_longer(
            cols = -!!sym(variant_type),
            names_to = "sample_id",
            values_to = "variant_class"
        )

    # variants <- unique(stat_data$variant_class)
    # message(
    #     paste0(" - Variants: ", paste(variants, collapse = ", "))
    # )
    
    ## Calculate the alteration frequencies
    group1_freq_data <- stat_data |> 
        filter(sample_id %in% group1_samples) |> 
        group_by(!!sym(variant_type)) |>
        summarise(
            group1_altered = sum(variant_class !=""),
            .groups = "drop"
        ) |> 
        mutate(group1_total = group1_total)

    group2_freq_data <- stat_data |> 
        filter(sample_id %in% group2_samples) |> 
        group_by(!!sym(variant_type)) |>
        summarise(
            group2_altered = sum(variant_class !=""),
            .groups = "drop"
        ) |> 
        mutate(group2_total = group2_total)

    freq_data <- group1_freq_data |> 
        left_join(group2_freq_data, by = variant_type) |> 
        mutate(
            group1_freq = group1_altered / group1_total,
            group2_freq = group2_altered / group2_total,
            freq_diff = group2_freq - group1_freq
        )

    ## Do the fisher test
    all_variants <- unique(freq_data[[variant_type]])

    fisher_results <- map_dfr(
        all_variants,
        function(variant) {

            cur_data <- freq_data |> filter(!!sym(variant_type) == variant)

            group1_altered <- cur_data |> pull(group1_altered)

            group2_altered <- cur_data |> pull(group2_altered)

            group1_total <- cur_data |> pull(group1_total)

            group2_total <- cur_data |> pull(group2_total)

            ## Fisher's exact test
            contingency_table <- matrix(
                    c(
                        group1_altered, group1_total - group1_altered,
                        group2_altered, group2_total - group2_altered
                    ), 
                    byrow = TRUE,
                    nrow = 2
            )
            colnames(contingency_table) <- c("Altered", "Not_Altered")
            rownames(contingency_table) <- c("Group1", "Group2")

            fisher_test <- fisher.test(
                contingency_table,
                alternative = alternative
            )

            tibble(
                {{variant_type}} := variant,
                fisher_p_val = fisher_test$p.value,
                # fisher_odds_ratio = fisher_test$estimate
                # conf_int_low = fisher_test$conf.int[1],
                # conf_int_high = fisher_test$conf.int[2]
            )
        }
    ) |> 
    ## Adjust for multiple testing
    mutate(
        fisher_p_adj = p.adjust(
            fisher_p_val, 
            method = p_adjust_method
        )
    )
    
    ## Merge the data
    freq_data |> 
        left_join(fisher_results, by = variant_type) |> 
        mutate(
            group_comparison = group_comparison, 
            .after = {{variant_type}}
        ) |> 
        arrange(fisher_p_val)
}

LoadCBioPortalSarcomaData <- function(
    mutate_gene_freq_thres = 0.01,
    cna_gene_freq_thres = 0.1,
    structural_variants_freq_thres = 0.01

) {
    ## Public cohorts of sarcoma from cBioPortal:
    ##  - Adult Soft Tissue Sarcomas (TCGA, Cell 2017)
    ##  - Sarcoma (MSK, Nat Commun. 2022)
    ##  - Sarcoma (MSKCC/Broad, Nat Genet 2010)
    ##  - Sarcoma (TCGA, Firehose Legacy)
    ##  - Sarcoma (TCGA, PanCancer Atlas)
    ##  - Sarcoma (UCLA, Cell 2024)
    ##  - Soft Tissue Cancer (TCGA, GDC)
    ##  - Soft Tissue Sarcoma (MSK, 2025)
    ##  - Soft Tissue and Bone Sarcoma (MSK, Nat Commun 2022)

    ## "---------------------------------------------------------------------"
    ## Mutated gene data
    ## "---------------------------------------------------------------------"
    mutated_gene_file <- "data/public/cBioPortal_sarcoma_mutated_genes.txt"
    
    mutated_gene_data <- read_tsv(
        here(mutated_gene_file), 
        show_col_types = FALSE
    ) |>
        clean_names() |>
        select(gene, freq, is_cancer_gene_source_onco_kb) |>
        dplyr::rename(
            symbol = gene, 
            freq = freq,
            is_oncokb = is_cancer_gene_source_onco_kb
        ) |>
        mutate(freq = parse_number(freq)/100) |>
        arrange(desc(freq)) |>
        filter(freq >= mutate_gene_freq_thres)

    ## "---------------------------------------------------------------------"
    ## CNA gene data
    ## "---------------------------------------------------------------------"
    cna_gene_file <- "data/public/cBioPortal_sarcoma_CNA_genes.txt"

    cna_gene_data <- read_tsv(
            cna_gene_file, 
            show_col_types = FALSE
    ) |> 
        clean_names() |> 
        select(
            gene, cytoband, cna, freq, is_cancer_gene_source_onco_kb
        ) |> 
        dplyr::rename(
            symbol = gene,
            cytoband = cytoband,
            cn_status = cna,
            freq = freq,
            is_oncokb = is_cancer_gene_source_onco_kb
        ) |> 
        mutate(freq = parse_number(freq)/100) |> 
        arrange(desc(freq)) |> 
        filter(freq >= cna_gene_freq_thres)

    ## "---------------------------------------------------------------------"
    ## Structural variants
    ## "---------------------------------------------------------------------"
    structural_variants_file <- "data/public/cBioPortal_sarcoma_structural_variants.txt"

    structural_variants_data <- read_tsv(
        here(structural_variants_file), 
        show_col_types = FALSE
    ) |> 
        clean_names() |> 
        select(
            gene, number_structural_variant, freq, is_cancer_gene_source_onco_kb
        ) |> 
        dplyr::rename(
            symbol = gene,
            number_structural_variant = number_structural_variant,
            freq = freq,
            is_oncokb = is_cancer_gene_source_onco_kb
        ) |> 
        mutate(freq = parse_number(freq)/100) |> 
        arrange(desc(freq)) |> 
        filter(freq >= structural_variants_freq_thres)

    ## Return
    list(
        mutated_gene_data = mutated_gene_data,
        cna_gene_data = cna_gene_data,
        structural_variants_data = structural_variants_data
    )
}

LoadDFSPStudyData <- function() {
    
    ## There are two major DFSP study
    ## peng_2022 (https://doi.org/10.1111/bjd.20976)
    ## smith_2025 (https://www.sciencedirect.com/science/article/pii/S089339522500033X)
    
    list(
        ## Modern pathology
        smith_2025 = list(
            mutated_genes = c(
                "FANCC", "TERT", "CHEK2",
                "KMT2D", "BRIP1", "ERCC2", "KMT2A",
                "MGA", "MTOR", "PTEN", "RAD50",
                "RIT1", "TP53"
            ),
            cna_genes = c(
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
            mutated_genes = c(
                "MUC6", "MUC4", "KMT2C",
                "HERC2", "HLA-A", "PRSS3", "PDE4DIP",
                "PRSS1", "RBMXL1", "BRCA1", "CTBP2",
                "HLA-B", "HLA-DQA2", "LILRA6", "OR8U1",
                "PABPC1", "PCSK1", "ZXDB"
            ),
            cna_genes = c(
                "TERT", "CDKN2A", "CDKN2B",
                "AKT1", "NFKBIA", "BTBD7", "SPHK1", 
                "ITGB4", "COL1A1", "PDGFB", "PDGFD"
            )
        )
    )
}

GetPCGRVariantMat <- function(
    data,
    variant_type = "symbol",
    variant_class = "variant_class",
    is_somatic_matched = TRUE
) {

    ## Load clincial info
    clinical_info <- LoadClinicalInfo(is_somatic_matched = is_somatic_matched)

    all_samples <- unique(clinical_info$Sample.ID)

    ## Keep releveant data
    cur_data <- data |> filter(sample_id %in% all_samples)

    cur_samples <- unique(cur_data$sample_id)

    missed_sample <- setdiff(clinical_info$Sample.ID, cur_samples)

    ## Save a clearned data: row=gene name, column=sample_id
    variant_tbl <- cur_data |> 
        group_by(sample_id, !!sym(variant_type)) |>
        summarise(
            variant_class = paste(
                unique(!!sym(variant_class)), collapse = ","
            ),
            .groups = "drop"
        ) |> 
        pivot_wider(
            names_from = sample_id,
            values_from = variant_class,
            values_fill = ""
        )

    if (length(missed_sample) > 0) {
    
        message(
            paste0(" - No variants could be found for these samples: ")
        )
        message(
            paste0(" > ", missed_sample, collapse = "\n")
        )

        ## Add the missed sample
        for (i in missed_sample) {
            variant_tbl <- variant_tbl |> 
                add_column(!!i := "")
        }
    }

    variant_mat <- variant_tbl |> 
        column_to_rownames(var = variant_type) |>
        as.matrix()
    
    ## print out the dimension
    n_sample <- ncol(variant_mat)
    n_variant <- nrow(variant_mat)

    message(
        paste0(" - samples: ", n_sample, ", variants: ", n_variant)
    )

    variant_mat
}

GetVariantMat <- function(
    data,
    variant_type = "symbol",
    variant_class = "variant_class",
    is_somatic_matched = TRUE
) {

    ## Load clincial info
    clinical_info <- LoadClinicalInfo(is_somatic_matched = is_somatic_matched)

    all_samples <- unique(clinical_info$Sample.ID)

    ## Keep releveant data
    cur_data <- data |> filter(sample_id %in% all_samples)

    cur_samples <- unique(cur_data$sample_id)

    missed_sample <- setdiff(clinical_info$Sample.ID, cur_samples)

    ## Save a clearned data: row=gene name, column=sample_id
    variant_tbl <- cur_data |> 
        group_by(sample_id, !!sym(variant_type)) |>
        summarise(
            variant_class = paste(
                unique(!!sym(variant_class)), collapse = ","
            ),
            .groups = "drop"
        ) |> 
        pivot_wider(
            names_from = sample_id,
            values_from = variant_class,
            values_fill = ""
        )

    if (length(missed_sample) > 0) {
    
        message(
            paste0(" - No variants could be found for these samples: ")
        )
        message(
            paste0(" > ", missed_sample, collapse = "\n")
        )

        ## Add the missed sample
        for (i in missed_sample) {
            variant_tbl <- variant_tbl |> 
                add_column(!!i := "")
        }
    }

    variant_mat <- variant_tbl |> 
        column_to_rownames(var = variant_type) |>
        as.matrix()
    
    ## print out the dimension
    n_sample <- ncol(variant_mat)
    n_variant <- nrow(variant_mat)

    message(
        paste0(" - samples: ", n_sample, ", variants: ", n_variant)
    )

    variant_mat
}

GetVariantFreq <- function(
    mat,
    variant_type = "symbol"
) {

    as_tibble(mat, rownames = variant_type) |>
        pivot_longer(
            cols = -!!sym(variant_type),
            names_to = "sample_id",
            values_to = "variant_class"
        ) |>
        group_by(!!sym(variant_type)) |>
        summarise(
            n_altered = sum(variant_class != ""),
            n_total = n_distinct(sample_id),
            freq = n_altered / n_total
        ) |>
        arrange(desc(n_altered))
}