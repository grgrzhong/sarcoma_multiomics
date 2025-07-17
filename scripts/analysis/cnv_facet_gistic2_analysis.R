#!/usr/bin/env Rscript

source(here::here("scripts/lib/study_lib.R"))

gistic2_file <- here("data/wes/Processed/cnv_facets_DFSP_cohort_gistic2.tsv")
gistic2_data <- read_tsv(
    gistic2_file,
    show_col_types = FALSE
)

## Split sample into groups
epic_gistic_dir <- here("data/epic/GISTIC2")
group_paths <- dir_ls(epic_gistic_dir, glob = "*.seg", recurse = TRUE)
group_names <- path_file(group_paths) |> str_remove("\\.seg$")

## Load all GISTIC2 segment files
gistic_list <- map2(
    group_names,
    group_paths,
    function(name, path) {
        message(paste0(" - Processing ", name, " = ", path))
        read_tsv(
            path,
            col_types = cols(
                Sample = col_character(),
                Chromosome = col_character(),
                Start = col_integer(),
                End = col_integer(),
                Num_Probes = col_integer(),
                Segment_Mean = col_double()
            )
        )
    }
) |>
    setNames(group_names)

epic_sample_groups <- list()
for (name in names(gistic_list)) {

    gistic2_data <- gistic_list[[name]]

    samples <- unique(gistic2_data$Sample)

    message(paste0(" - Found ", length(samples), " samples in group: ", name))

    epic_sample_groups[[name]] <- samples
}

## Load DFSP clinical information
clinical_info <- LoadDFSPClinicalInfo()

wes_sample_groups <- list()
for (name in names(epic_sample_groups)) {

    samples <- epic_sample_groups[[name]]

    # Filter clinical info for these samples
    filtered_clinical_info <- clinical_info |> 
        filter(Sample.ID %in% samples)

    wes_sample_groups[[name]] <- filtered_clinical_info
}

## "=========================================================================="
## Generate the WES sample groups
## "=========================================================================="
## Patient levele analysi:
## Main: Choosen the sample with the highest grade: presence of FST
## Main.Met

wes_sample_groups <- list(

    ## All WES samples, EPIC=150, WES=161
    `all_tumors` = clinical_info,

    ## FS-DFSP group, EPIC=23, WES=23
    `FS-DFSP` = clinical_info |> 
        filter(FST.Group %in% c("FS-DFSP")),

    ## FS-DFSP Meta, EPIC=6, WES=11
    `FS-DFSP_Meta` = clinical_info |> 
        filter(FST.Group %in% c("FS-DFSP")) |> 
        filter(Metastasis == "Yes"),

    ## FS-DFSP Primary Rrecurrence, EPIC=17, WES=17
    `FS-DFSP_Pri_Rec` = clinical_info |> 
        filter(FST.Group %in% c("FS-DFSP")) |>
        filter(Specimen.Nature %in% c("Primary", "Recurrence")),
    
    ## Post-FST, EPIC=39, WES=41
    `Post-FST` = clinical_info |> 
        filter(FST.Group %in% c("Post-FST")),

    ## Pre-FST, EPIC=40, WES=45
    `Pre-FST` = clinical_info |> 
        filter(FST.Group %in% c("Pre-FST")),

    ## U-DFSP, EPIC=48, WES=52
    `U-DFSP` = clinical_info |> 
        filter(FST.Group %in% c("U-DFSP")),
    
    ## Patient level Metastasis
    Meta = clinical_info |> 
        filter(Main.Met == "Yes") |> 
        filter(Metastasis == "Yes"),
    
    `Non-Meta` = clinical_info |> 
        filter(Main.Met == "Yes") |> 
        filter(Metastasis == "No")
)

# Save the WES sample groups to an Excel file
write_xlsx(
    wes_sample_groups,
    here("data/clinical/DFSP_WES_sample_groups.xlsx")
)

sample_groups <- LoadDFSPSampleGroups()


facets_gistic_list <- map(
    sample_groups,
    \(group) {
        gistic2_data |> filter(Sample %in% group)
    }
)

## Save the GISTIC2 data for each sample group
for (name in names(facets_gistic_list)) {
    group_data <- facets_gistic_list[[name]]
    out_dir <- here("data/wes/GISTIC2", name)
    dir_create(out_dir)
    write.table(
        group_data,
        file = paste0(out_dir, "/", name, ".tsv"),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
}
