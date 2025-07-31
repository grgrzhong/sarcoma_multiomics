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
        library(ComplexHeatmap)
        library(reshape2)
    })
)
source(here::here("scripts/lib/study_lib.R"))

## Output directories
data_dir <- "data/processed"
plot_dir <- "figures/wes/pcgr"

## Capture size
capture_size <- 34 

## Clinical information
clinical_info <- LoadDFSPClinicalInfo()

## "==========================================================================="
## Collect PCGR CNV data --------------
## "==========================================================================="
## Merge all PCGR outputs
pcgr_dir <- "data/wes/PCGR"
pcgr_cna_raw <- CollectPCGRCNAData(dir = pcgr_dir)

## Define the cnv status based on total copy number
pcgr_cna_dat <- pcgr_cna_raw |> 
    mutate(
        cn_status = case_when(
            # Homozygous deletion (complete loss)
            total_cn == 0 ~ "homdel",

            # Heterozygous loss (one copy lost)
            total_cn == 1 ~ "hetloss",
            
            # Copy-neutral Loss of Heterozygosity (cnLOH)
            # Total CN = 2 but allelic imbalance (e.g., 2+0, not 1+1)
            total_cn == 2 & cn_minor == 0 ~ "cnloh",

            # Normal diploid (balanced)
            total_cn == 2 & cn_major == 1 & cn_minor == 1 ~ "neutral",

            # Low-level gain (3-4 copies)
            total_cn %in% 3:4 ~ "gain",

            # High-level amplification (more than 4 copies)
            total_cn > 4 ~ "amp",

            # Other cases (e.g., undefined)
            TRUE ~ "undefined"
        ),
        .after = total_cn
    ) |>  
    mutate(
        variant_class = case_when(
            total_cn == 0 ~ "del",
            total_cn >= 4 ~ "amp",
            TRUE ~ "undefined"
        )
    ) |> 
    relocate(variant_class, .after = cn_status)

## We should use the total cohort size as the denominator for frequency calculations
## as we are answer the question "What proportion of samples have this alteration in my cohort?"
total_cohort_size <- length(unique(pcgr_cna_dat$sample_id))

message(
    paste("Total samples with PCGR CNV data:", total_cohort_size)
)

## Analyze only the variants that are satify the thresholds
pcgr_cna_tbl <- pcgr_cna_dat |> 
    filter(total_cn == 0 | total_cn >= 4)

length(unique(pcgr_cna_tbl$sample_id))

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
