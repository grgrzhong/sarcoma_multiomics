#!/usr/bin/env Rscript

## "==========================================================================="
## General configurations ----
## "==========================================================================="
## Load required libraries and functions

suppressPackageStartupMessages(
    suppressWarnings({
        library(maftools)
    })
)

source(here::here("scripts/lib/study_lib.R"))

## Output directories
gistic_dir <- "data/wes/GISTIC2"
# group_name <- "all_tumors"

groups <- dir_ls(gistic_dir) |> path_file()

sample_groups <- LoadDFSPSampleGroups()
clinical_info <- LoadDFSPClinicalData()

for (group in groups) {
    ## Read the GISTIC2 data
    gistic_obj <- readGistic(
        gisticAllLesionsFile = here(
            gistic_dir, group, "all_lesions.conf_99.txt"
        ),
        gisticAmpGenesFile = here(gistic_dir, group, "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, group, "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, group, "scores.gistic"),
        verbose = TRUE
    )

    ## Plot parameters
    width <- 8
    height <- 4
    out_dir <- "figures/wes/gistic2/cnv_facets"

    ## Plot the chromosomal plot
    for (img in c("png", "pdf")) {
        file <- here(
            out_dir,
            paste0("wes_cnvfacets_gistic2_", group, "_chromplot.", img)
        )

        if (img == "png") {
            CairoPNG(
                file = file, width = width, height = height, res = 300,
                fonts = "Arial", units = "in"
            )
        } else {
            CairoPDF(
                file = file, width = width, height = height, fonts = "Arial"
            )
        }

        gisticChromPlot(
            gistic = gistic_obj,
            fdrCutOff = 0.25,
            txtSize = 0.6,
            cytobandTxtSize = 0.5,
            color = c("#D95F02", "#1B9E77"),
            markBands = "all",
            ref.build = "hg38",
            y_lims = c(-3, 3)
        )

        message(paste0("Saving plot: ", file))
        dev.off()
    }

    ## Bubble plot
    for (img in c("png", "pdf")) {
        file <- here(
            out_dir,
            paste0("wes_cnvfacets_gistic2_", group, "_bubbleplot.", img)
        )

        if (img == "png") {
            CairoPNG(
                file = file, width = width, height = height, res = 300,
                fonts = "Arial", units = "in"
            )
        } else {
            CairoPDF(
                file = file, width = width, height = height, fonts = "Arial"
            )
        }

        gisticBubblePlot(
            gistic = gistic_obj,
            color = c("#D95F02", "#1B9E77"),
            markBands = NULL,
            log_y = TRUE,
            fdrCutOff = 0.25,
            txtSize = 0.6
        )

        message(paste0("Saving plot: ", file))
        dev.off()
    }

    ## Oncoplot
    for (img in c("png", "pdf")) {
        file <- here(
            out_dir,
            paste0("wes_cnvfacets_gistic2_", group, "_oncoplot.", img)
        )

        if (img == "png") {
            CairoPNG(
                file = file, width = width, height = height, res = 300,
                fonts = "Arial", units = "in"
            )
        } else {
            CairoPDF(
                file = file, width = width, height = height, fonts = "Arial"
            )
        }

        gisticOncoPlot(
            gistic = gistic_obj,
            sortByAnnotation = FALSE,
            top = 10,
            gene_mar = 10,
            barcode_mar = 10,
            sepwd_genes = 0.5,
            bandsToIgnore = NULL,
            removeNonAltered = TRUE,
            colors = c(
                Amp = "#D95F02",
                Del = "#1B9E77"
            ),
            SampleNamefontSize = 0.6,
            fontSize = 0.8,
            legendFontSize = 0.7,
            annotationFontSize = 1.2,
            borderCol = "white",
            bgCol = "#CCCCCC"
        )

        message(paste0("Saving plot: ", file))
        dev.off()
    }
}

## "==========================================================================="
## Find enriched cytobands and genes ----
## "==========================================================================="
gistic_dir <- here("data/wes/GISTIC2/cnv_facets")

comparisons <- list(
    `Non-Meta_vs_Meta` = c("Non-Meta", "Meta"),
    `U-DFSP_vs_Pre-FST` = c("U-DFSP", "Pre-FST"),
    `U-DFSP_vs_Post-FST` = c("U-DFSP", "Post-FST"),
    `U-DFSP_vs_FS-DFSP` = c("U-DFSP", "FS-DFSP"),
    `Pre-FST_vs_Post-FST` = c("Pre-FST", "Post-FST"),
    `Pre-FST_vs_FS-DFSP` = c("Pre-FST", "FS-DFSP"),
    `Post-FST_vs_FS-DFSP` = c("Post-FST", "FS-DFSP")
)

stat_list <- list()

for (comparison in names(comparisons)) {

    group1 <- comparisons[[comparison]][1]
    group2 <- comparisons[[comparison]][2]

    message(paste0("Comparing: ", group1, " vs ", group2))

    stat_list[[comparison]] <- GetStatResGistic2(
        gistic_dir = gistic_dir,
        group1 = group1,
        group2 = group2,
        freq_thres = 0.2,
        qval_thres = 0.25
    )
}

for (name in names(stat_list)) {

    output <- stat_list[[name]]

    write_xlsx(
        output,
        path = here(
            "results/wes", paste0("wes_cnvfacets_gistic2_", name, ".xlsx")
        )
    )
}

test <- stat_list$`Non-Meta_vs_Meta`[["cytoband"]] |> 
    as_tibble()

test |> 
    filter(fisher_qvalue < 0.05) |> 
    filter(enriched_in_group2 | depleted_in_group2) |> 
    view()

stat_res <- GetStatResGistic2(
    gistic_dir = gistic_dir,
    group1 = group1,
    group2 = group2,
    freq_thres = 0.2,
    qval_thres = 0.25
)

# PlotGistic2Chrom <- 

# scores <- "/mnt/f/projects/sarcoma_multiomics/data/wes/GISTIC2/cnv_facets/FS-DFSP/scores.gistic"

# library(BSgenome.Hsapiens.UCSC.hg38)
# chrom_info <- tibble(
#     chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38),
#     chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)
# )
# chrom_info$chromNum <- 1:length(chrom_info$chromName)
# chrom_info <- chrom_info[1:22, ]

# ## Calculate cumulative chromosome lengths
# chrom_info$chromlengthCumsum <- cumsum(as.numeric(chrom_info$chromlength))

# ## Calculate chromosome start positions from 0
# chrom_info$chromStartPosFrom0 <- c(0, chrom_info$chromlengthCumsum[-nrow(chrom_info)])

# ## Calculate chromosome middle positions from 0
# tmp_middle <- diff(c(0, chrom_info$chromlengthCumsum)) / 2
# chrom_info$chromMiddlePosFrom0 <- chrom_info$chromStartPosFrom0 + tmp_middle

# scores <- read.table(
#     file = scores, 
#     header = TRUE, 
#     sep = "\t", 
#     stringsAsFactors = FALSE
# )
# ## Add chromosome ID and start/end positions
# chromID <- scores$Chromosome
# scores$StartPos <- scores$Start + chrom_info$chromStartPosFrom0[chromID]
# scores$EndPos <- scores$End + chrom_info$chromStartPosFrom0[chromID]

# range(scores$G.score)
# ## Adjust G-scores for deletions for visualization
# # scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1
# scores <- scores |> 
#     mutate(
#         G.score = case_when(
#             Type == "Del" ~ G.score * -1,
#             Type == "Amp" ~ G.score
#         )
#     )

# library(ggplot2)
# library(ggsci)

# ggplot(scores, aes(StartPos, G.score)) +
#     geom_area(aes(group = Type, fill = factor(Type, levels = c("Del", "Amp")))) +
#     scale_fill_lancet(guide = guide_legend(reverse = T)) +
#     geom_vline(data = chrom_info, mapping = aes(xintercept = chromlengthCumsum), linetype = 2) +
#     geom_text(data = chrom_info, aes(x = chromMiddlePosFrom0, y = 0.2, label = chromName)) +
#     scale_x_continuous(expand = c(0, -1000), limits = c(0, 2.9e9)) +
#     ylim(-0.3, 0.3) +
#     theme_minimal()

# chrom_info$ypos <- rep(c(0.2, 0.25), 11)

# ggplot(scores, aes(StartPos, G.score)) +
#     geom_area(aes(group = Type, fill = factor(Type, levels = c("Del", "Amp")))) +
#     scale_fill_lancet(guide = guide_legend(reverse = T), name = "Type") +
#     geom_vline(data = chrom_info, mapping = aes(xintercept = chromlengthCumsum), linetype = 2) +
#     geom_text(data = chrom_info, aes(x = chromMiddlePosFrom0, y = ypos, label = chromName)) +
#     scale_x_continuous(expand = c(0, -1000), limits = c(0, 2.9e9), name = NULL, labels = NULL) +
#     ylim(-0.3, 0.3) +
#     theme_minimal() +
#     theme(
#         legend.position = "top",
#         axis.text.y = element_text(color = "black", size = 14),
#         axis.title.y = element_text(color = "black", size = 16)
#     )
