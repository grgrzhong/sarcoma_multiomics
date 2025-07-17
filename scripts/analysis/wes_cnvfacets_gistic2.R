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
## Comparative Analysis for FST Transformation ----
## "==========================================================================="
## Run comparative analysis if both groups exist
if ("U-DFSP" %in% groups && "FS-DFSP" %in% groups) {
    message("Running FST transformation comparative analysis...")

    ## Load GISTIC objects
    untransformed_gistic <- readGistic(
        gisticAllLesionsFile = here(gistic_dir, "U-DFSP", "all_lesions.conf_99.txt"),
        gisticAmpGenesFile = here(gistic_dir, "U-DFSP", "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, "U-DFSP", "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, "U-DFSP", "scores.gistic"),
        verbose = FALSE
    )

    transformed_gistic <- readGistic(
        gisticAllLesionsFile = here(gistic_dir, "FS-DFSP", "all_lesions.conf_99.txt"),
        gisticAmpGenesFile = here(gistic_dir, "FS-DFSP", "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, "FS-DFSP", "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, "FS-DFSP", "scores.gistic"),
        verbose = FALSE
    )

    ## Perform comparative analysis
    transformation_cnas <- compare_gistic_groups(
        untransformed_gistic, transformed_gistic, "FST_transformation"
    )

    transformation_genes <- compare_gistic_genes(
        untransformed_gistic, transformed_gistic
    )

    ## Create comparison plots
    dir.create(here("figures/wes/gistic2/comparative"), recursive = TRUE, showWarnings = FALSE)

    ## Plot CNA frequency comparison
    for (img in c("png", "pdf")) {
        file <- here(
            "figures/wes/gistic2/comparative",
            paste0("fst_transformation_cna_comparison.", img)
        )

        if (img == "png") {
            CairoPNG(file = file, width = 10, height = 6, res = 300, units = "in")
        } else {
            CairoPDF(file = file, width = 10, height = 6)
        }

        plot_cna_comparison(transformation_cnas)
        dev.off()
    }

    ## Save gene-level results
    saveRDS(transformation_genes, here("results", "FST_transformation_genes.rds"))

    ## Print summary
    message("FST Transformation Analysis Complete:")
    message(paste(
        "- Total significant CNAs:",
        sum(transformation_cnas$transformation_enriched &
            transformation_cnas$significant_in_transformed)
    ))
    message(paste(
        "- Transformation-specific amplified genes:",
        length(transformation_genes$amp_specific)
    ))
    message(paste(
        "- Transformation-specific deleted genes:",
        length(transformation_genes$del_specific)
    ))
}

## Function to plot CNA comparison
plot_cna_comparison <- function(cna_data) {
    significant_cnas <- cna_data %>%
        filter(transformation_enriched == TRUE, significant_in_transformed == TRUE) %>%
        arrange(desc(frequency_difference))

    if (nrow(significant_cnas) > 0) {
        p <- ggplot(significant_cnas, aes(
            x = reorder(cytoband, frequency_difference),
            y = frequency_difference
        )) +
            geom_col(aes(fill = alteration_type)) +
            coord_flip() +
            scale_fill_manual(values = c("Amplification" = "#D95F02", "Deletion" = "#1B9E77")) +
            labs(
                title = "FST Transformation-Enriched CNAs",
                x = "Cytogenetic Band",
                y = "Frequency Difference (Transformed - Untransformed)",
                fill = "Alteration Type"
            ) +
            theme_minimal()

        print(p)
    }
}
