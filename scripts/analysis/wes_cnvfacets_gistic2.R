#!/usr/bin/env Rscript

suppressPackageStartupMessages(
    suppressWarnings({
        library(maftools)
        library(ComplexHeatmap)
        library(circlize)
    })
)

source(here::here("scripts/lib/study_lib.R"))
source(here::here("conf/study_conf.R"))

## "==========================================================================="
## Process EPIC CNV data for GISTIC2 input ----
## "==========================================================================="
epic_cnv_data <- read_tsv(
    here("data/epic/GISTIC2/all_tumors/all_tumors.tsv"),
    show_col_types = FALSE
)

sample_groups <- LoadSampleGroupInfo(
    is_somatic_matched = FALSE
)[["groups"]][1:11]

sample_groups$all_tumors <- NULL

## Check if there are any NA values for all columns
na_counts <- sapply(epic_cnv_data, function(x) sum(is.na(x)))
if (any(na_counts > 0)) {
    message("There are NA values in the EPIC CNV data:")
    print(na_counts[na_counts > 0])
} else {
    message("No NA values found in the EPIC CNV data.")
}

## Split the GISTIC2 input data by sample groups
for (sample_group in names(sample_groups)) {

    n_samples <- length(sample_groups[[sample_group]])
    
    message(
        paste0("Processing sample group: ", sample_group, " (", n_samples, ")")
    )

    group_data <- epic_cnv_data |> 
        filter(Sample %in% sample_groups[[sample_group]])

    out_dir <- here("data/epic/GISTIC2", paste0("/", sample_group))
    
    dir_create(out_dir)
    
    write.table(
        group_data,
        file = paste0(out_dir, "/", sample_group, ".tsv"),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
}

## "==========================================================================="
## Collect WES CNV Facet data  ----
## "==========================================================================="
cnv_facet_dir <- here("data/wes/CNV/cnv_facets")

tsv_files <- dir_ls(cnv_facet_dir, glob = "*.tsv", recurse = TRUE) |> 
    keep(~ str_detect(path_file(.x), "^[^.]+\\.tsv$"))

sample_names <- path_file(tsv_files) |> str_remove("\\.tsv$")
message(
    paste0(
        "Found ", length(sample_names), " sample with CNV Facets data."
    )
)

cnv_facet_data <- tibble(
    sample = sample_names,
    file = tsv_files
) |>
    mutate(
        data = map2(
            file,
            sample,
            \(f, s) {
                message(
                    paste0("Loading CNV Facets data = ", s)
                )
                dat <- read_tsv(
                    f, 
                    show_col_types = FALSE, 
                    col_types = cols(.default = "c")
                ) 
                dat
            },
            .progress = TRUE
        )
    ) |> 
    unnest(data)

## Save the combined CNV Facets data
SaveData(
    cnv_facet_data,
    dir = "data/processed",
    filename = "wes_cnv_facets_DFSP_cohort_raw"
)

## "==========================================================================="
## Export CNV facet purity, ploidy data ----
## "==========================================================================="
cnv_facet_data |> 
    select(sample, purity, ploidy, dipLogR, est_insert_size, emflags) |>
    distinct() |> 
    arrange(sample) |>
    write_xlsx(
        here("outputs", "wes_cnv_facet_DFSP_cohort_purity_ploidy.xlsx"),
    )

clinical_info <- LoadClinicalInfo()
matched_samples <- clinical_info |> 
    filter(Somatic.Status == "Matched") |>
    pull(Sample.ID)

cnv_facet_data |> 
    select(sample, purity, ploidy, dipLogR, est_insert_size, emflags) |>
    distinct() |> 
    arrange(dipLogR) |> 
    # filter(sample %in% matched_samples) |> 
    filter(is.na(purity) | is.na(ploidy)) 

cnv_facet_data |> pull(SVTYPE) |> unique() |> sort()

## "==========================================================================="
## Process WES CNV Facet data for GISTIC2 input ----
## "==========================================================================="
## Get the GISTIC2 input data
gistic2_data <- cnv_facet_data |> 
    filter(SVTYPE != "NEUTR") |> 
    select(sample, CHROM, POS, END, SVTYPE, TCN_EM, NUM_MARK, CNLR_MEDIAN, dipLogR) |> 
    mutate(
        CHROM = str_remove(CHROM, "chr"),
        POS = as.integer(POS),
        END = as.integer(END),
        NUM_MARK = as.integer(NUM_MARK),
        CNLR_MEDIAN = as.numeric(CNLR_MEDIAN),
        dipLogR = as.numeric(dipLogR)
    ) |> 
    # allele-specific copy number, 
    # Total copy number >=4 (amplified) and = 0 (deep deletion)
    # filter(TCN_EM == 0  | TCN_EM > 4) |>
    ## Remove the X, Y chromosomes
    filter(!CHROM %in% c("X", "Y")) |> 
    ## log(total copy number) (https://github.com/mskcc/facets/issues/167)
    mutate(Segment_Mean = CNLR_MEDIAN - dipLogR) |> 
    ## Filter out segments that are not estimateable
    ## Gistic2 run to problem with NA values
    filter(!is.na(Segment_Mean)) |> 
    dplyr::rename(
        Sample = sample,
        Chromosome = CHROM,
        Start = POS,
        End = END,
        Num_Probes = NUM_MARK
    ) |> 
    ## Make sure the the segments to be start < end
    mutate(
        Start_fixed = if_else(Start > End, End, Start),
        End_fixed = if_else(Start > End, Start, End),
        Chromosome = as.integer(Chromosome),
    ) |> 
    select(
        Sample, Chromosome, Start_fixed, End_fixed, Num_Probes, Segment_Mean
    ) |> 
    dplyr::rename(
        Start = Start_fixed,
        End = End_fixed
    ) |> 
    ## keep segment mean only 3 significant digits
    mutate(
        Segment_Mean = round(Segment_Mean, 3)
    )

## Check if there are any NA values for all columns
na_counts <- sapply(gistic2_data, function(x) sum(is.na(x)))
if (any(na_counts > 0)) {
    message("There are NA values in the GISTIC2 data:")
    print(na_counts[na_counts > 0])
} else {
    message("No NA values found in the GISTIC2 data.")
}

## Distribution of Segment Means
hist(
    gistic2_data$Segment_Mean,
    breaks = 100,
    main = "Distribution of Segment Means",
    xlab = "Segment Mean",
    col = "lightblue"
)
dev.off()

## Save the GISTIC2 input data
sample_groups <- LoadSampleGroupInfo(
    is_somatic_matched = FALSE
)[["groups"]][1:11]

write_tsv(
    gistic2_data,
    file = here("data/processed/wes_cnv_facets_DFSP_cohort_gistic2.tsv"),
    na = "NA",
    quote = "none"
)

## Split the GISTIC2 data by sample groups
somatic_matched <- LoadSampleGroupInfo(
    is_somatic_matched = TRUE
)[["groups"]][1:11]

somatic_unmatched <- LoadSampleGroupInfo(
    is_somatic_matched = FALSE
)[["groups"]][1:11]

sample_categories <- list(
    somatic_matched = somatic_matched,
    somatic_unmatched = somatic_unmatched
)

for (i in names(sample_categories)) {

    sample_groups <- sample_categories[[i]]

    for (sample_group in names(sample_groups)) {

        n_samples <- length(sample_groups[[sample_group]])
        
        message(
            paste0("Processing sample group: ", sample_group, " (", n_samples, ")")
        )

        group_data <- gistic2_data |> 
            filter(Sample %in% sample_groups[[sample_group]])

        out_dir <- here("data/wes/GISTIC2", paste0(i, "/", sample_group))
        
        dir_create(out_dir)
        
        write.table(
            group_data,
            file = paste0(out_dir, "/", sample_group, ".tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
    }
}

## Certain samples got not estimated CNV data
gistic2_data |> pull(Sample) |> unique() |> length()

## "==========================================================================="
## WES CNV Find significant changed regions in subgroups ----
## "==========================================================================="
## Do the analysis for both somatic matched (147) and all samples (161)
gistic_dirs <- list(
    ## Tumor samples with matched normal sample
    somatic_matched = here("data/wes/GISTIC2/somatic_matched"),
    ## All samples, including tumor samples without matched normal
    somatic_unmatched = here("data/wes/GISTIC2/somatic_unmatched")
)

## Run the fisher's exact test for each group comparison
group_stat_res_list <- list()

for (i in names(gistic_dirs)) {

    if (i == "somatic_unmatched") {
        
        group_comparisons <- LoadSampleGroupInfo(
            is_somatic_matched = FALSE
        )[["comparisons"]]

    } else if (i == "somatic_matched") {
        
        group_comparisons <- LoadSampleGroupInfo(
            is_somatic_matched = TRUE
        )[["comparisons"]]

    }

    gistic_dir <- gistic_dirs[[i]]
    
    message(paste0("Processing GISTIC2 data = ", gistic_dir))

    for (gc in names(group_comparisons)) {

        group1 <- str_split(gc, "_vs_")[[1]][1] |> 
            str_remove("\\(\\d+\\)")

        group2 <- str_split(gc, "_vs_")[[1]][2] |> 
            str_remove("\\(\\d+\\)")

        message(paste0(" - Processing ", gc))

        stat_res <- GetGistic2CytobandStatRes(
            gistic_dir = gistic_dir,
            group1 = group1,
            group2 = group2,
            freq_thres = 0.2,
            qval_thres = 0.25
        )

        ## Keep only for the significant altered regions based on gistic qvalue)
        group_stat_res_list[[i]][[gc]]  <- stat_res |> 
            dplyr::filter(sig_in_group1 | sig_in_group2)
    }

}

## Save the group_stat_res_list
SaveData(
    group_stat_res_list,
    dir = "data/processed",
    filename = "wes_cnvfacets_gistic2_group_stat_res_list"
)

# Combine all group_comparisons into one dataframe
sig_cytoband_list <- list()
for (i in names(group_stat_res_list)) {

    group_stat_re <- group_stat_res_list[[i]]
    
    group_comparisons <- names(group_stat_re)
    
    ## Focus on the main comparisons
    # group_comparisons <- c(
    #     group_comparisons[3],
    #     group_comparisons[4],
    #     group_comparisons[5],
    #     group_comparisons[1]
    # )

    ## Combine all the group comparisons into one dataframe
    sig_cytoband_list[[i]] <- map_dfr(
        group_comparisons,
        function(gc) {
            group_stat_re[[gc]] |>
                mutate(group_comparison = gc) |>
                # select(
                #     cytoband, alteration_type, group_comparison, 
                #     sig_in_group1, sig_in_group2,
                #     significantly_different, 
                #     enriched_in_group2, depleted_in_group2
                # ) |>
                mutate(
                    # Create combined cytoband_alteration identifier
                    cytoband_alteration = paste0(alteration_type, "_", cytoband),
                    # Create significance value
                    is_sig = case_when(
                        significantly_different ~ "Yes",
                        TRUE ~ "No"
                    )
                ) |>
                filter(is_sig == "Yes")
        }
    )
}

## Save the combined data
SaveData(
    sig_cytoband_list,
    dir = "data/processed",
    filename = "wes_cnvfacets_gistic2_sig_cytoband_list"
)

length(unique(sig_cytoband_list$somatic_matched$cytoband))
length(unique(sig_cytoband_list$somatic_unmatched$cytoband))

## "==========================================================================="
## Heatmap Plot the significant cytobands across group comparsions ----
## "==========================================================================="
sig_cytoband <- sig_cytoband_list[["somatic_matched"]]

cytoband_sort <- sig_cytoband |> 
    group_by(cytoband_alteration) |>
    summarise(
        n_sig = sum(is_sig == "Yes"),
        .groups = "drop"
    ) |>
    arrange(desc(n_sig)) |> 
    pull(cytoband_alteration)

group_sort <- sig_cytoband |> 
    pull(group_comparison) |>
    unique()

plot <- sig_cytoband |> 
    mutate(
        cytoband_alteration = factor(
            cytoband_alteration, 
            levels = rev(cytoband_sort)
        ),
        group_comparison = factor(group_comparison, levels = group_sort),
    ) |>
    ggplot(
        aes(
            x = group_comparison, 
            y = cytoband_alteration, 
            fill = factor(is_sig)
        )
    ) +
    geom_tile(color = "black", linewidth = 0.1) +
    scale_fill_manual(
        values = c("No" = "darkgray", "Yes" = "red"),
        labels = c(
            "No" = "No (q-value > 0.05)", "Yes" = "Yes (q-value < 0.05)"
        ),
        name = "Fisher's Exact Test\n q-value"
    ) +
    labs(
        x = NULL, y = NULL, 
        title = "GISTIC2 Group Comparison Significance"
    ) +
    # theme_minimal() +
    plot_theme() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank()
    )

filename <- "wes_cnvfacets_gistic2_group_comparison_significance_heatmap_all_groups"

SavePlot(
    plot = plot,
    width = 6, 
    height = 16, 
    only_png = FALSE, 
    dir = "figures/gistic2", 
    filename = filename
)

## "==========================================================================="
## Heatmap oncoplot single group ----
## "==========================================================================="


## "==========================================================================="
## Heatmap oncoplot combined group ----
## "==========================================================================="
## Entire cohort gistic2 data
all_tumour_gistic <- LoadGistic2Data(
    gistic_dir = "data/wes/GISTIC2/somatic_matched",
    group = "all_tumors"
)

all_tumour_gistic_tbl <- all_tumour_gistic@cytoband.summary |> 
    as_tibble()

all_tumour_gistic_cytobands <- all_tumour_gistic_tbl |> 
    pull(Cytoband) |> 
    unique()

sig_cytobands <- sig_cytoband |> pull(cytoband) |> unique()

setdiff(sig_cytobands, all_tumour_gistic_cytobands) |> sort()

show_cytobands <- c(sig_cytobands, all_tumour_gistic_cytobands) |> 
    unique()

## Load PCGR annotate data

pcgr_cna_raw <- LoadData(
    dir = "data/processed", 
    filename = "wes_pcgr_DFSP_cohort_cnv_raw_data"
)
names(pcgr_cna_raw)

pcgr_cna_raw |> 
    filter(symbol == "CDKN2A")  |> 
    select(1:11) |> 
    left_join(
        LoadClinicalInfo() |> 
            select(Sample.ID, FST.Group, Metastasis),
        by = c("sample_id" = "Sample.ID")
    ) |> 
    group_by(FST.Group) |>
    summarise(
        n_samples = n_distinct(sample_id),
        .groups = "drop"
    )

pcgr_cna_raw |> 
    pull(variant_class) |> 
    unique()

## Filter the significant cytobands in the PCGR CNA data
pcgr_sig_cytoband <- pcgr_cna_raw |> 
    select(sample_id:cytoband) |> 
    separate_wider_delim(
        cytoband, 
        delim = ":", 
        names = c("chromosome", "cytoband")
    ) |> 
    mutate(
        chr_num = str_extract(chromosome, "\\d+"),
        .after = chromosome
    ) |> 
    mutate(
        expanded_cytobands = map2_chr(
            chr_num,
            cytoband,
            function(chr, cb) {
                
                if (str_detect(cb, "-")) {
                    parts <- str_split(cb, "-")[[1]]
                    start_band <- parts[1]
                    end_band <- parts[2]

                    ## Create full cytoband chrosome
                    start_full <- paste0(chr, start_band)
                    end_full <- paste0(chr, end_band)

                    ## Return
                    paste(start_full, end_full, sep = ",")
                } else {
                    ## Create full cytoband chromosome
                    paste0(chr, cb)
                }
            },
            .progress = TRUE
        ),
        matched_cytoband = map_chr(
            expanded_cytobands,
            function(ec) {
                
                individual_cytobands <- str_split(ec, ",")[[1]]

                ## Check which cytoband is in the significant cytobands
                matches <- individual_cytobands[
                    individual_cytobands %in% show_cytobands
                ]

                if (length(matches) > 0) {
                    paste0(matches, collapse = ";")
                } else {
                    NA_character_
                }
            },
            .progress = TRUE
        ),
        is_sig_cytoband = !is.na(matched_cytoband)
    )

all_samples <- pcgr_cna_raw |> pull(sample_id) |> unique()

all_sig_cytobands <- pcgr_sig_cytoband |> 
    filter(is_sig_cytoband) |> 
    pull(matched_cytoband) |>
    str_split(";") |>
    unlist() |> 
    unique() |> 
    sort()

cytoband_freq <- pcgr_sig_cytoband |> 
    filter(is_sig_cytoband) |> 
    mutate(
        alteration = case_when(
            total_cn == 0 ~ "homdel",
            total_cn == 1 ~ "hetdel", 
            total_cn >= 3 ~ "gain",
            TRUE ~ "neutral"
        )
    ) |> 
    filter(alteration != "neutral") |> 
    separate_rows(matched_cytoband, sep = ";") |> 
    # Count unique samples with alterations per cytoband
    group_by(matched_cytoband) |> 
    summarise(
        n_altered_samples = n_distinct(sample_id),
        .groups = "drop"
    ) |> 
    mutate(
        alteration_frequency = n_altered_samples / length(all_samples)
    ) |> 
    arrange(desc(alteration_frequency))

oncoplot_matrix_data <- pcgr_sig_cytoband |> 
    filter(is_sig_cytoband) |> 
    mutate(
        alteration = case_when(
            total_cn == 0 ~ "homdel",
            total_cn == 1 ~ "hetdel", 
            total_cn >= 3 ~ "gain",
            TRUE ~ "neutral"
        )
    ) |> 
    filter(alteration != "neutral") |> 
    separate_rows(matched_cytoband, sep = ";") |> 
    # Handle potential duplicates by taking the most severe alteration
    group_by(sample_id, matched_cytoband) |> 
    summarise(
        alteration = case_when(
            "homdel" %in% alteration ~ "homdel",
            "gain" %in% alteration ~ "gain", 
            "hetdel" %in% alteration ~ "hetdel",
            TRUE ~ first(alteration)
        ),
        .groups = "drop"
    ) |>
    complete(
        sample_id = all_samples,
        matched_cytoband = all_sig_cytobands,
        fill = list(alteration = "")
    ) |> 
    pivot_wider(
        names_from = sample_id,
        values_from = alteration,
        values_fill = ""
    ) 

## Prepare the oncoplot matrix
oncoplot_matrix <- oncoplot_matrix_data |> 
    column_to_rownames("matched_cytoband") |> 
    as.matrix()

GenerateCytobandOncoplot(
    mat = oncoplot_matrix,
    sample_annotation = c("FST.Group", "Metastasis"),
    sample_sorted_by = "FST.Group",
    sample_sorted_level = c(
        "U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"
    ),
    width = 18,
    height = 12,
    dir = "figures/gistic2",
    filename = "wes_cnvfacets_gistic2_oncoplot_somatic_matched"
)

## "==========================================================================="
## Gistic oncoplot all tumors ----
## "==========================================================================="
gistic_dirs <- list(
    somatic_matched = here("data/wes/GISTIC2/somatic_matched"),
    somatic_unmatched = here("data/wes/GISTIC2/somatic_unmatched")
)

gistic_dir <- here("data/wes/GISTIC2/somatic_matched")

gistic_obj <- readGistic(
    gisticAllLesionsFile = here(gistic_dir, "all_tumors", "all_lesions.conf_99.txt"),
    gisticAmpGenesFile = here(gistic_dir, "all_tumors", "amp_genes.conf_99.txt"),
    gisticDelGenesFile = here(gistic_dir, "all_tumors", "del_genes.conf_99.txt"),
    gisticScoresFile = here(gistic_dir, "all_tumors", "scores.gistic"),
    verbose = FALSE
)

sig_cytobands <- sig_cytoband_data |> 
    pull(cytoband_alteration) |> 
    unique()

all_tumors_cytobands_data <- gistic_obj@data |> 
    as_tibble() |> 
    mutate(
        cytoband = sapply(
            Cytoband,
            function(x) str_split(x, ":")[[1]][2]
        )
    ) |> 
    mutate(
        cytoband_status = paste(Variant_Classification, cytoband, sep = "_")
    ) |> 
    filter(cytoband %in% all_cytobands$Cytoband) |> 
    pull(Cytoband) |> 
    unique() |> 
    sort()

table(sig_cytobands %in% all_tumors_cytobands_data$cytoband_status)

gistic_obj@gis.scores
gistic_obj@numericMatrix

clinical_data <- LoadClinicalInfo() |> 
    mutate(
        Tumor_Sample_Barcodes = Sample.ID
    ) |> 
    select(
        Tumor_Sample_Barcodes, 
        Specimen.Nature, 
        FST.Group, 
        Metastasis
    ) |> 
    as.data.frame()

gisticOncoPlot(
    gistic = gistic_obj,
    # bands = show_cytobands,
    clinicalData = clinical_data,
    clinicalFeatures = "FST.Group",
    showTumorSampleBarcodes = TRUE,
    sortByAnnotation = TRUE,
    # colors = c(
    #     Amp = "#D95F02",
    #     Del = "#1B9E77"
    # ),
    SampleNamefontSize = 0.6,
    fontSize = 0.8,
    legendFontSize = 0.7,
    annotationFontSize = 1.2,
    borderCol = "white",
    bgCol = "#CCCCCC"
)

## "==========================================================================="
## Gistic plot group comparsions ---------
## "==========================================================================="
plot_para <- list(
    plot_dir = "figures/gistic2",
    gistic_data = list(
        epic_cnv = list(
            gistic_dir = here("data/epic/GISTIC2"),
            plot_filename = "epic_cnv_gistic2"
        ),
        somatic_matched_samples = list(
            gistic_dir = here("data/wes/GISTIC2/somatic_matched"),
            plot_filename = "wes_cnvfacets_gistic2_somatic_matched_samples"
        ),
        all_samples = list(
            gistic_dir = here("data/wes/GISTIC2/somatic_unmatched"),
            plot_filename = "wes_cnvfacets_gistic2_all_samples"   
        )
    ),
    group_list = list(
        Metastatsis = c("Non-Meta", "Meta"),
        FST.subtype = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
    )
)

for (i in names(plot_para$gistic_data)) {

    if (i == "somatic_matched") {

        is_somatic_matched <- TRUE

    } else {
        
        is_somatic_matched <- FALSE
    }

    for (group in names(plot_para$group_list)) {
        
        message(
            paste0("\nPlotting GISTIC2 results for: ", group)
        )

        gistic_dir <- plot_para$gistic_data[[i]]$gistic_dir

        groups <- plot_para$group_list[[group]]

        filename <- paste0(
            plot_para$gistic_data[[i]]$plot_filename, "_", group
        )

        plot_dir <- plot_para$plot_dir

        PlotGisticGroupComparsion(
            gistic_dir = gistic_dir,
            groups = groups,
            is_somatic_matched = is_somatic_matched,
            fdrCutOff = 0.25,
            markBands = show_cytobands,
            color = c("#c82118", "#2c5496"),
            ref.build = "hg38",
            cytobandOffset = 0.05,
            txtSize = 0.8,
            cytobandTxtSize = 0.7,
            y_lims = NULL,
            width = 8,
            height = 8,
            plot_dir = plot_dir,
            filename = filename
        )
    }
}

## "==========================================================================="
## Explore the GISTIC2 results ----
## "==========================================================================="
clinical_info <- LoadClinicalInfo()

## Output directories
gistic_dir <- "data/wes/GISTIC2/somatic_matched"
# group_name <- "all_tumors"

sample_groups <- LoadSampleGroupInfo()
# groups <- dir_ls(gistic_dir) |> path_file()
groups <- names(sample_groups)[1:11]

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

    n_sample <- length(sample_groups[[group]])

    ## Plot parameters
    width <- 8
    height <- 4
    out_dir <- "figures/wes/gistic2/cnv_facets"

    dir_create(out_dir)

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
            color = c("#c82118", "#2c5496"),
            markBands = "all",
            ref.build = "hg38",
            y_lims = c(-2, 5)
        )
        title(
            main = paste0(group, " (", n_sample, ")"),
            family = "Arial"
        )
        
        message(paste0("Saving plot: ", file))
        dev.off()
    }
}
    # ## Bubble plot
    # for (img in c("png", "pdf")) {
    #     file <- here(
    #         out_dir,
    #         paste0("wes_cnvfacets_gistic2_", group, "_bubbleplot.", img)
    #     )

    #     if (img == "png") {
    #         CairoPNG(
    #             file = file, width = width, height = height, res = 300,
    #             fonts = "Arial", units = "in"
    #         )
    #     } else {
    #         CairoPDF(
    #             file = file, width = width, height = height, fonts = "Arial"
    #         )
    #     }

    #     gisticBubblePlot(
    #         gistic = gistic_obj,
    #         color = c("#D95F02", "#1B9E77"),
    #         markBands = NULL,
    #         log_y = TRUE,
    #         fdrCutOff = 0.25,
    #         txtSize = 0.6
    #     )

    #     message(paste0("Saving plot: ", file))
    #     dev.off()
    # }

    # ## Oncoplot
    # for (img in c("png", "pdf")) {
    #     file <- here(
    #         out_dir,
    #         paste0("wes_cnvfacets_gistic2_", group, "_oncoplot.", img)
    #     )

    #     if (img == "png") {
    #         CairoPNG(
    #             file = file, width = width, height = height, res = 300,
    #             fonts = "Arial", units = "in"
    #         )
    #     } else {
    #         CairoPDF(
    #             file = file, width = width, height = height, fonts = "Arial"
    #         )
    #     }

    #     gisticOncoPlot(
    #         gistic = gistic_obj,
    #         sortByAnnotation = FALSE,
    #         top = 10,
    #         gene_mar = 10,
    #         barcode_mar = 10,
    #         sepwd_genes = 0.5,
    #         bandsToIgnore = NULL,
    #         removeNonAltered = TRUE,
    #         colors = c(
    #             Amp = "#D95F02",
    #             Del = "#1B9E77"
    #         ),
    #         SampleNamefontSize = 0.6,
    #         fontSize = 0.8,
    #         legendFontSize = 0.7,
    #         annotationFontSize = 1.2,
    #         borderCol = "white",
    #         bgCol = "#CCCCCC"
    #     )

    #     message(paste0("Saving plot: ", file))
    #     dev.off()
    # }


## "==========================================================================="
## Non-Meta vs Meta chromosomal plot ---------
## "==========================================================================="
## Non-Meta vs Meta - Combined Plot with Aligned Axes
for (img in c("png", "pdf")) {
    file <- here(
        out_dir,
        paste0("wes_cnvfacets_gistic2_NonMeta_vs_Meta_combined.", img)
    )

    if (img == "png") {
        CairoPNG(
            file = file, width = 8, height = 8, res = 300,
            fonts = "Arial", units = "in"
        )
    } else {
        CairoPDF(
            file = file, width = 8, height = 8, fonts = "Arial"
        )
    }

    # Set up 2 rows, 1 column layout
    par(mfrow = c(2, 1), mar = c(3, 4, 3, 2))

    # Read GISTIC objects for both groups
    gistic_obj_nonmeta <- readGistic(
        gisticAllLesionsFile = here(gistic_dir, "Non-Meta", "all_lesions.conf_99.txt"),
        gisticAmpGenesFile = here(gistic_dir, "Non-Meta", "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, "Non-Meta", "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, "Non-Meta", "scores.gistic"),
        verbose = FALSE
    )
    
    gistic_obj_meta <- readGistic(
        gisticAllLesionsFile = here(gistic_dir, "Meta", "all_lesions.conf_99.txt"),
        gisticAmpGenesFile = here(gistic_dir, "Meta", "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, "Meta", "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, "Meta", "scores.gistic"),
        verbose = FALSE
    )

    n_sample_nonmeta <- length(sample_groups[["Non-Meta"]])

    n_sample_meta <- length(sample_groups[["Meta"]])

    # Plot 1: Non-Meta (top panel, no x-axis labels)
    gisticChromPlot(
        gistic = gistic_obj_nonmeta,
        fdrCutOff = 0.25,
        txtSize = 0.6,
        cytobandTxtSize = 0.5,
        color = c("#c82118", "#2c5496"),
        ref.build = "hg38",
        y_lims = c(-2, 5)  # Same y-limits for both plots
    )
    title(main = paste0("Non-Meta (n=", n_sample_nonmeta, ")"), family = "Arial", cex.main = 1.2)

    # Plot 2: Meta (bottom panel, with x-axis labels)
    par(mar = c(5, 4, 1, 2))  # More bottom margin for x-axis labels
    gisticChromPlot(
        gistic = gistic_obj_meta,
        fdrCutOff = 0.25,
        txtSize = 0.6,
        cytobandTxtSize = 0.5,
        color = c("#c82118", "#2c5496"),
        ref.build = "hg38",
        y_lims = c(-2, 5)  # Same y-limits for both plots
    )
    title(
        main = paste0("Meta (n=", n_sample_meta, ")"), 
        family = "Arial", cex.main = 1.2
    )

    # Reset par settings
    par(mfrow = c(1, 1))
    
    message(paste0("Saving combined plot: ", file))
    dev.off()
}


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
