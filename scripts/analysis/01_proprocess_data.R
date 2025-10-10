#!/usr/bin/env Rscript
##############################################################################
## Description: Preprocess EPIC, WES data
## Author:  Zhong Guorui
## Created: 2025-06-15
##############################################################################
## Load required functions and configurations
source(here::here("conf/study_conf.R"))

## "==========================================================================="
## Preprocess EPIC CNV data for GISTIC2 input ----
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
## Preprocess WES CNV FACETS segment data  ----
## "==========================================================================="
cnv_facet_dir <- here("data/WES/CNV/cnv_facets")

tsv_files <- dir_ls(cnv_facet_dir, glob = "*.tsv", recurse = TRUE) |> 
    keep(~ str_detect(path_file(.x), "^[^.]+\\.tsv$"))

sample_names <- path_file(tsv_files) |> str_remove("\\.tsv$")
message(
    paste0(
        "Found ", length(sample_names), " sample with CNV Facets data."
    )
)

## The following 14 samples are tumour samples without matched normal, 
## the outputs were generate using DFSP-336-N:
## DFSP-030-T-P1
## DFSP-030-T-P2
## DFSP-031-T
## DFSP-032-T
## DFSP-037-T
## DFSP-051-T-P1
## DFSP-051-T-P2
## DFSP-072-T
## DFSP-106-T
## DFSP-112-T
## DFSP-123-T
## DFSP-150-T
## DFSP-170-T
## DFSP-171-T

cnv_facet_raw <- tibble(
    sample = sample_names,
    file = tsv_files
) |>
    mutate(
        data = map2(
            file,
            sample,
            \(f, s) {
                message(
                    paste0(" *Loading CNV Facets data = ", s)
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
    cnv_facet_raw,
    dir = "data/processed",
    filename = "wes_cnv_facets_raw"
)

## Export CNV facet purity, ploidy data
clinical_info <- LoadClinicalInfo()
purity_ploidy_data <- cnv_facet_raw |> 
    select(sample, purity, ploidy, dipLogR, est_insert_size, emflags) |>
    distinct() |> 
    arrange(sample) |> 
    filter(sample %in% clinical_info$Sample.ID)

write_csv(
    purity_ploidy_data,
    here("outputs", "wes_cnv_facet_purity_ploidy.csv")
)

## Within the somatic samples (147), 5 samples have missing purity or ploidy values
## DFSP-042-T
## DFSP-060-T
## DFSP-331-T
## DFSP-341-T
## DFSP-354-T
missed_purity_ploidy <- purity_ploidy_data |> 
    filter(is.na(purity) | is.na(ploidy))

print(missed_purity_ploidy)

## "==========================================================================="
## Calculate WES CNV FACETS segment FGA  ----
## "==========================================================================="
fga_list <- list()

for (cnv_file in tsv_files) {
    
    sample_name <- path_file(cnv_file) |> str_remove("\\.tsv$")

    fga_list[[sample_name]] <- FacetCalculateFGA(cnv_file = cnv_file)
}

## Save data to a tibble
fga_tbl <- NULL

for (i in names(fga_list)) {
    fga_tbl <- rbind(
        fga_tbl,
        tibble(
            sample_id = i,
            fga = fga_list[[i]][["fga"]],
            altered_length = fga_list[[i]][["altered_length"]],
            total_length = fga_list[[i]][["total_length"]]
        )
    )
}
clinical_info <- LoadClinicalInfo()

SaveData(
    fga_tbl |> filter(sample_id %in% clinical_info$Sample.ID),
    filename = "wes_cnv_facet_fga",
    dir = "data/processed"
)

## "==========================================================================="
## Preprocess WES CNV FACETS CNV gene data  ----
## "==========================================================================="
## * bedtools annotated gene data
cnv_facet_dir <- here("data/WES/CNV/cnv_facets")

tsv_files <- dir_ls(cnv_facet_dir, glob = "*.annotated.tsv", recurse = TRUE)

sample_names <- path_file(tsv_files) |> 
    str_remove("\\.annotated.tsv$")

message(
    paste0(
        "Found ", length(sample_names), " sample with CNV Facets data."
    )
)

## The following 14 samples are tumour samples without somatic matched normal, 
## the facet outputs were generate using DFSP-336-N:
## DFSP-030-T-P1
## DFSP-030-T-P2
## DFSP-031-T
## DFSP-032-T
## DFSP-037-T
## DFSP-051-T-P1
## DFSP-051-T-P2
## DFSP-072-T
## DFSP-106-T
## DFSP-112-T
## DFSP-123-T
## DFSP-150-T
## DFSP-170-T
## DFSP-171-T

cnv_facet_gene_raw <- tibble(
    sample = sample_names,
    file = tsv_files
) |>
    mutate(
        data = map2(
            file,
            sample,
            \(f, s) {
                message(
                    paste0(" *Loading CNV Facets data = ", s)
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
    unnest(data) |> 
    rename(sample_id = sample)

unique(cnv_facet_gene_raw$Chromosome)
unique(cnv_facet_gene_raw$SVTYPE)
clinical_info <- LoadClinicalInfo()

## Keep samples with somatic matched normal only, 147 samples
cnv_facet_gene_tbl <- cnv_facet_gene_raw |> 
    filter(sample_id %in% clinical_info$Sample.ID) |>
    ## The not estimateable Lesser (minor) copy number segments were filtered out
    dplyr::filter(!(LCN_EM %in% c(".", NA, "NA"))) |>
    filter(!Chromosome %in% c("chrX", "chrY")) |> 
    ## Remove the neutral events
    filter(SVTYPE != "NEUTR") |> 
    mutate(
        Start_Position = as.integer(Start_Position),
        End_Position = as.integer(End_Position),
        TCN_EM = as.integer(TCN_EM),
        LCN_EM = as.integer(LCN_EM),
        ploidy = as.numeric(ploidy)
    ) |> 
    mutate(
        cn_major = as.integer(round(TCN_EM - LCN_EM)),
        cn_minor = as.integer(round(LCN_EM)),
        total_cn = as.integer(round(TCN_EM))
    ) |> 
    dplyr::select(
        sample_id, Chromosome, Start_Position, End_Position,
        cn_major, cn_minor, total_cn, purity, ploidy, SVTYPE, 
        Chr, Start, End, Gene
    ) |> 
    ## Define the CNV events
    mutate(
        cn_status = case_when(
            ## "Deep deletion"
            !is.na(ploidy) & total_cn == 0 ~ "HOMDEL",
            
            ## Deletion
            !is.na(ploidy) & (total_cn -ploidy) <= -1 ~ "DEL",
            
            ## Neutral
            !is.na(ploidy) & 
                ((total_cn -ploidy) > -1) &
                ((total_cn -ploidy) < 1) ~ "NEUTRAL",
            
            ## Amplification
            !is.na(ploidy) & 
                ((total_cn - ploidy) >= 1) &
                ((total_cn - ploidy) < 3) ~ "GAIN",

            ## Deep amplification
            !is.na(ploidy) & (total_cn - ploidy) >=3 ~ "AMP",
        ),
        .after = cn_minor
    ) |> 
    ## Keep only changed CNV events
    dplyr::filter(cn_status != "NEUTRAL") |> 
    clean_names() |> 
    rename(symbol = gene) |> 
    mutate(
        start = as.integer(start),
        end = as.integer(end)
    )

## Add the cytoband data
cnv_facet_gene_tbl <- MapGene2Cytoband(cnv_facet_gene_tbl)

cnv_facet_gene_tbl |> 
    filter(symbol == c("SPRY4", "KITLG"))

## Save the combined CNV bedtools annotated gene data
SaveData(
    cnv_facet_gene_tbl,
    dir = "data/processed",
    filename = "cnv_facet_gene_tbl"
)

## "==========================================================================="
## Prepare WES CNV FACETS segment data for GISTIC2 input ----
## "==========================================================================="
## Get the GISTIC2 input data
gistic2_data <- cnv_facet_raw |> 
    ## !!! Should use all segments
    # filter(SVTYPE != "NEUTR") |> 
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
    ## !!! Gistic2 run to problem with NA values
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
    file = here("outputs/wes_cnv_facets_gistic2_input_data.tsv"),
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

        out_dir <- here("data/WES/GISTIC2", paste0(i, "/", sample_group))
        
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

## "==========================================================================="
## Preprocess PCGR SNV/Indel data --------------
## "==========================================================================="
## Merge all SNV/Indel data from PCGR annotation
pcgr_snv_indel_raw <- CollectPCGRSNVINDELData(
    dir = "data/WES/PCGR", 
    sheet_name = "SOMATIC_SNV_INDEL"
)

## Total variants = 289599
total_variants <- nrow(pcgr_snv_indel_raw)
message(paste("Total number of variants = ", total_variants))

## Save the collected PCGR data
SaveData(
    pcgr_snv_indel_raw, 
    dir = "data/processed", 
    filename = "wes_pcgr_snv_indels_raw"
)

## Include all indels, actionable SNVs, and oncogenic SNVs, 
## and predictedly pathogenic SNVs
## total_variants = 197730, all samples
pcgr_snv_indel_tbl <- FilterPCGRSNVIndels(
    pcgr_snv_indel_raw,
    min_dp_tumor = NULL,
    min_vaf_tumor = NULL,
    min_dp_control = NULL,
    max_vaf_control = NULL,
    max_population_frequency = 0.001
)

## Total variants = 24177, all samples have variants
pcgr_snv_indel_tbl <- FilterPCGRSNVIndels(
    pcgr_snv_indel_raw,
    min_dp_tumor = 20,
    min_vaf_tumor = 0.03,
    min_dp_control = 10,
    max_vaf_control = 0.01,
    max_population_frequency = 0.001
)

## Total variants = 8003, 1 sample (DFSP-087-T) detect no variants
pcgr_snv_indel_tbl <- FilterPCGRSNVIndels(
    pcgr_snv_indel_raw,
    min_dp_tumor = 20,
    min_vaf_tumor = 0.05,
    min_dp_control = 10,
    max_vaf_control = 0.01,
    max_population_frequency = 0.001
)

## Save the data
SaveData(
    pcgr_snv_indel_tbl,
    dir = "data/processed",
    filename = "wes_pcgr_snv_indels_tbl"
)

## Calculate the TMB
maf_tbl <- ConvertPCGRToMaftools(pcgr_tbl = pcgr_snv_indel_tbl)

## Create a MAF object
maf_obj <- read.maf(
    maf = maf_tbl, 
    clinicalData = LoadClinicalInfo(is_somatic_matched = FALSE),
    # cnTable = cnv_data,
    verbose = TRUE
)

tcgaCompare(
    maf = maf_obj,
    cohortName = "FS-DFSP",
    logscale = TRUE, 
    capture_size = capture_size
)

tmb_data <- tmb(
    maf_obj, 
    captureSize = 34, 
    logScale = TRUE, 
    plotType = "classic", 
    verbose = TRUE
)

message(
    paste(
        "Median TMB = ", 
        round(median(tmb_data$total_perMB, na.rm = TRUE), 2)
    )
)

SaveData(
    tmb_data, 
    dir = "data/processed", 
    filename = "wes_pcgr_tmb_data"
)



## "==========================================================================="
## Preprocess PCGR CNV data --------------
## "==========================================================================="
pcgr_cnv_raw <- CollectPCGRCNVData()

## Save the raw PCGR CNV data
SaveData(
    pcgr_cnv_raw, 
    dir = "data/processed", 
    filename = "wes_pcgr_cnv_raw"
)

## 5 samples do not have facet estimated purity
## DFSP-042-T, DFSP-060-T, DFSP-331-T, DFSP-341-T, DFSP-354-T
cnv_facet_raw <- LoadData(
    dir = "data/processed",
    filename = "wes_cnv_facets_DFSP_cohort_raw"
)
cnv_facet_raw |> 
    filter(is.na(purity) | is.na(ploidy)) |> 
    select(sample, purity, ploidy) |> 
    distinct()

## Define the CNV events
## 10.1016/j.jtho.2021.02.023
pcgr_cnv_tbl <- pcgr_cnv_raw |> 
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
        chromosome = str_extract(var_id, "chr[0-9XY]+"), .after = "sample_id",
        start = as.integer(start),
        end = as.integer(end)
    ) |> 
    ## Exclude the sex chromosomes
    filter(!(grepl("chrx|chrX|chry|chrY", var_id))) |> 
    arrange(desc(facet_ploidy))

## One sample do not have non-neutral CNA: "DFSP-169-T"
length(unique(pcgr_cnv_tbl$sample_id))
setdiff(pcgr_cnv_raw$sample_id, pcgr_cnv_tbl$sample_id)

## Save the filtered Non-Neutral PCGR CNV data
SaveData(
    pcgr_cnv_tbl,
    dir = "data/processed",
    filename = "wes_pcgr_cnv_tbl"
)

## "==========================================================================="
## Preprocess ASCETS CNV arm data --------------
## "==========================================================================="
## https://github.com/beroukhim-lab/ascets
## Get the hg38 cytoband data
ascet_dir <- "scripts/modules/ascets"
# hg19_file <- here(ascet_dir, "cytoband_coordinates_hg19.txt")
# hg19_cytoband <- read_tsv(hg19_file, show_col_types = FALSE)
# data(cytobands.hg38)
# hg38_cytoband <- cytobands.hg38 |> 
#     as_tibble() |> 
#     dplyr::rename(arm = band) |> 
#     dplyr::select(-stain) |> 
#     dplyr::select(all_of(colnames(hg19_cytoband))) |> 
#     mutate(chrom = str_replace_all(chrom, "chr", ""))

# write_tsv(
#     hg38_cytoband,
#     file = here(ascet_dir, "cytoband_coordinates_hg38.txt"),
#     quote = "none"
# )

## Prepare the data for ASCETS analysis
required_input_cols <- c(
    "sample", "CHROM", "POS", "END", "NUM_MARK", "CNLR_MEDIAN"
)

clinical_info <- LoadClinicalInfo()
cnv_facet_raw <- LoadData(
    dir = "data/processed",
    filename = "wes_cnv_facets_raw"
)

## Convert to custom format
cnv_facet_data <- cnv_facet_raw |>
    dplyr::filter(sample %in% clinical_info$Sample.ID) |>
    dplyr::filter(SVTYPE != "NEUTR") |>
    ## The not estimateable Lesser (minor) copy number segments were filtered out
    dplyr::filter(!(LCN_EM %in% c(".", NA, "NA"))) |>
    dplyr::select(all_of(required_input_cols)) |>
    dplyr::rename(
        chrom = CHROM,
        segment_start = POS,
        segment_end = END,
        num_mark = NUM_MARK,
        log2ratio = CNLR_MEDIAN
    ) |>
    mutate(
        chrom = str_replace(chrom, "chr", ""),
        segment_start = as.integer(segment_start),
        segment_end = as.integer(segment_end)
    ) |>
    filter(!(chrom %in% c("X", "Y"))) |> 
    ## Make sure the the segments to be start < end
    mutate(
        segment_start_fixed = if_else(
            segment_start > segment_end, segment_end, segment_start
        ),
        segment_end_fixed = if_else(
            segment_start > segment_end, segment_start, segment_end
        )
    ) |>
    dplyr::select(
        sample,
        chrom,
        segment_start = segment_start_fixed, 
        segment_end = segment_end_fixed,
        num_mark, log2ratio
    ) |>
    ## Keep data types consistent
    mutate(
        chrom = as.character(chrom),
        segment_start = as.integer(segment_start),
        segment_end = as.integer(segment_end),
        num_mark = as.integer(num_mark),
        log2ratio = as.numeric(log2ratio)
    )

## Run ASCET arm level analysis
source(here(ascet_dir, "ascets_resources.R"))
# hg38_cytoband <- read_tsv(
#     file = here(ascet_dir, "cytoband_coordinates_hg38.txt"),
#     show_col_types = FALSE
# )

## arm-level
hg38_cytoband <- read_tsv(
    file = here(ascet_dir, "genomic_arm_coordinates_hg38.txt"),
    show_col_types = FALSE
)

ascets_input <- split(
    cnv_facet_data, 
    cnv_facet_data$sample
)

## Loop through the samples
ascets_arm_raw <- map(
    ascets_input,
    \(df) {
        message(paste0("\nProcessing ", unique(df$sample)))

        ascets_output <- ascets(
            cna = df,
            cytoband = hg38_cytoband,
            min_boc = 0.5,
            name = "DFSP-010-T-P1",
            keep_noisy = FALSE,
            threshold = 0.2,
            alteration_threshold = 0.7
        )
    },
    .progress = TRUE
)

print(ascets_arm_raw[[2]])
SaveData(
    ascets_arm_raw,
    dir = "data/processed/",
    filename = "wes_cnv_facets_ascets_arm_raw"
)

## Aneuploidy scores
ascets_aneuploidy_scores <- map_dfr(
    ascets_arm_raw,
    \(sample) {
        df <- sample[["aneu_scores"]]
    }
)

write_csv(
    ascets_aneuploidy_scores,
    file = here("outputs/wes_cnv_facets_ascets_aneuploidy_scores.csv")
)

## CNV arm matrix
ascets_arm_tbl <- map_dfr(
    ascets_arm_raw,
    \(sample) {
        df <- sample[["calls"]]
        df |> 
            pivot_longer(
                cols = -sample,
                names_to = "arm",
                values_to = "variant_class"
            ) |> 
            dplyr::rename(sample_id = sample) |> 
            dplyr::filter(
                !(variant_class %in% c(
                        ## insufficient data or coverage)
                        "NC",
                        ## low coverage, not enough data to make a reliable call
                        "LOWCOV",
                        ## no significant copy-number alteration
                        "NEUTRAL"
                    )
                )
            ) 
            ## Combine the chromosome and arm information
            # left_join(
            #     hg38_cytoband |> select(arm, chrom),
            #     by = "arm"
            # ) |> 
            # dplyr::mutate(
            #     arm = paste0("chr", chrom, ":", arm)
            # ) |> 
            # dplyr::select(-chrom)
    }
)

SaveData(
    ascets_arm_tbl,
    dir = "data/processed/",
    filename = "wes_cnv_facets_ascets_arm_tbl"
)

## "==========================================================================="
## Preprocess GISTIC2 CNV data --------------
## "==========================================================================="
gistic_dir <- here("data/WES/GISTIC2/somatic_matched")
gistic_data <- GetGisticCNData(
    gistic_dir = gistic_dir,
    group = "all_tumors",
    gistic_qval_thres = 0.5
)

SaveData(
    gistic_data,
    dir = "data/processed",
    filename = "wes_cnv_facets_gistic2_somatic_matched_all_tumors_qval0.5"
)

## "==========================================================================="
## Preprocess Pyclone input data --------------
## "==========================================================================="
## https://github.com/Roth-Lab/pyclone
## https://github.com/vanallenlab/PyCloneTSVGeneration_FACETS_or_TITAN

## Load the snv indel data
pcgr_snv_indel_raw <- LoadData(
    dir = "data/processed", 
    filename = "wes_pcgr_snv_indels_raw"
)

## Total variants = 24177, all samples have variants
pcgr_snv_indel_tbl <- FilterPCGRSNVIndels(
    pcgr_snv_indel_raw,
    min_dp_tumor = 20,
    min_vaf_tumor = 0.03,
    min_dp_control = 10,
    max_vaf_control = 0.01,
    max_population_frequency = 0.001
)

## Total variants = 8003, 1 sample (DFSP-087-T) detect no variants
pcgr_snv_indel_tbl <- FilterPCGRSNVIndels(
    pcgr_snv_indel_raw,
    min_dp_tumor = 20,
    min_vaf_tumor = 0.05,
    min_dp_control = 10,
    max_vaf_control = 0.01,
    max_population_frequency = 0.001
)

# pcgr_snv_indel_tbl <- LoadData(
#     dir = "data/processed",
#     filename = "wes_pcgr_snv_indels_tbl"
# )

clinical_info <- LoadClinicalInfo() |> 
    filter(FST.Group %in% c("Pre-FST", "Post-FST")) |> 
    filter(Number.of.samples.per.patient >=2)

snv_indel_data <- pcgr_snv_indel_tbl |> 
    filter(sample_id %in% clinical_info$Sample.ID)

pyclone_data <- GetPycloneInputData(
    snv_indel_data = snv_indel_data,
    vcf_dir = "data/WES/Mutect2",
    cnv_facet_dir = "data/WES/CNV/cnv_facets"
)

pyclone_dir <- "data/WES/PyClone_3"

## Save the Pyclone input data
for (i in names(pyclone_data)) {

    case_id <- str_extract(i, "DFSP-[0-9]+")

    # dir_delete(here("data/WES/PyClone", case_id))
    dir_create(here(pyclone_dir, case_id))

    message(
        paste0("Saving: ", i)
    )

    ## write the tsv
    write_tsv(
        pyclone_data[[i]],
        file = here(pyclone_dir, case_id, paste0(i, ".tsv")),
        na = "NA",
        quote = "none"
    )
    
}

rna_mat <- "data/RNA/merged.featureCounts.clean.txt"
rna_mat <- read_tsv(rna_mat, show_col_types = FALSE)
samples <- colnames(rna_mat)[-1]
length(samples)
clinical_info <- LoadClinicalInfo(is_somatic_matched = FALSE) |> 
    select(Sample.ID, FST.Group, Metastasis)

metadata <- tibble(Sample.ID = samples) |> 
    left_join(clinical_info, by = "Sample.ID")

write_csv(
    metadata,
    file = here("data/clinical/DFSP_RNAseq_metadata.csv")
)
