data_dir <- "data/processed"
figure_dir <- "figures/wes"

## "=========================================================================="
## Collect the annovar annotated variants
## "=========================================================================="
annovar_dir <- "data/wes/Mutect2"
annovar_tbl <- CollectAnnovarData(dir = annovar_dir)
filename <- "wes_annovar_DFSP_cohort_merged_tbl"
SaveData(annovar_tbl, dir = data_dir, filename = filename)

## Filter variants
annovar_tbl <- LoadData(dir = data_dir, filename = filename)

filtered_variants <- FilterAnnovarData(data = annovar_tbl)

## Find the enriched variants in FST ----
maf_obj <- read.maf(maf = filtered_variants)


SaveData <- function(obj, dir, filename) {
    
    ## Save data to qs file

    fs::dir_create(here(dir))

    qsave(obj, here(dir, paste0(filename, ".qs")))
}

LoadData <- function(dir, filename) {
    
    ## Load qs data
    qread(here(dir, paste0(filename, ".qs")))
    
}

loadCancerHotspot <- function(
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
    snv_hotspots <- loadCancerHotspot()[["snv_hotspots"]]

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