## Load required libraries and functions
source(here::here("scripts/lib/study_lib.R"))

pcgr_dir <- here("data/wes/pcgr")
pcgr_files <- dir_ls(pcgr_dir, glob = "*.pcgr.grch38.xlsx", recurse = TRUE)
sheet_names <- excel_sheets(pcgr_file[[1]])
message("Sheet names in the pcgr file: \n ", paste(sheet_names, collapse = "\n "))

# pcgr_file <- "/mnt/f/projects/250224_DFSP_Multiomics/data/wes/pcgr/DFSP-001-T/DFSP-001-T.pcgr.grch38.xlsx"
# pcgr_data_list <- map(
#     all_sheets,
#     ~ {
#         sheet <- .x
#         message("Processing sheet: ", sheet)
#         read_xlsx(pcgr_file, sheet = sheet)
#     }
# ) |>
#     setNames(all_sheets)

## Get sample names from PCGR files
sample_names <- pcgr_files |> 
    path_file() |> 
    str_remove(".pcgr.grch38.xlsx")

## Load all tbm data from PCGR files
tmb_data <- tibble(
    sample = sample_names,
    pcgr_file = pcgr_files
) |>
    mutate(
        pcgr_data = map2(
            pcgr_file,
            sample,
            \(f, s) {
                message("Processing sample = ", s)
                read_xlsx(f, sheet = "TMB")
            },
            .progress = TRUE
        )
    ) |> 
    unnest(pcgr_data)

group_info <- LoadDFSPGroupInfo()

plot_data <- tmb_data |> 
    select(SAMPLE_ID, TMB_UNIT, TMB_ESTIMATE, TMB_MEASURE, TMB_CSQ_REGEX) |> 
    left_join(group_info, by = c("SAMPLE_ID" = "Tumor_Sample_Barcode"))

