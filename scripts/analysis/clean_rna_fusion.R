
library(here)
library(tidyverse)
library(readxl)
library(writexl)

## QC metrics
qc_metrics_file <- "data/RNA/DFSP/Reports/qc_metrics_cohort_raw.xlsx"

qc_metrics_data <- read_excel(
    qc_metrics_file
)

pct_columns <- colnames(qc_metrics_data)
pct_columns <- grep("PCT", pct_columns, value = TRUE)

qc_metrics_tbl <- qc_metrics_data |> 
    mutate(
        across(
            all_of(pct_columns), 
            ~ as.numeric(str_remove(.x, "%"))
        )
    ) |> 
    ## is poor QC?
    mutate(
        is_Pass_QC = if_else(
            ## >= 30M reads
            TOTAL_READS >= 30e6 &
            
            PCT_PF_UQ_READS_ALIGNED >= 70 &

            ## mRNA bases
            (PCT_CODING_BASES + PCT_UTR_BASES) >= 75 &

            PCT_INTERGENIC_BASES <= 10 &
            PCT_INTRONIC_BASES <= 25 &
            PCT_RIBOSOMAL_BASES <= 10,

            "Yes",
            "No"
        ),
        .after = Sample
    )

## fusion results
fusion_report_file <- "data/RNA/DFSP/Reports/fusion_report_cohort_raw.xlsx"
fusion_report_data <- read_excel(
    fusion_report_file,
    sheet = "fusion_report"
)

colnames(fusion_report_data)
table(fusion_report_data$`arriba.reading-frame`)

pdgfd_pdgfb_info <- fusion_report_data |> 
    select(Fusion, Sample) |> 
    filter(grepl("PDGFB|PDGFD", Fusion))

known_drivers <- c(
    "COL1A1--PDGFB",
    "COL6A3--PDGFD",
    "SLC2A5--BTBD7",
    "TNC--PDGFD",
    "EMILIN2--PDGFD"
)

known_fusion_res <- fusion_report_data |> 
    filter(Fusion %in% known_drivers)

fusion_tbl <- filtered_fusion_res |> 
    bind_rows(known_fusion_res) |> 
    distinct()

filtered_fusion_res <- fusion_report_data |> 
    filter(Number_of_callers > 1) |> 
    bind_rows(known_fusion_res) |>
    distinct() |> 
    mutate(
        QC_Confidence = if_else(
            starfusion.junction_reads >=3 &
            starfusion.spanning_reads >=5 &
            (`arriba.reading-frame` %in% c("in-frame", "stop-codon")),
            "High",
            "Low"
        ),
        .after = Fusion
    ) |> 
    arrange(QC_Confidence) |> 
    filter(!grepl("ENSG|LINC", Fusion)) |> 
    filter(Fusion != "FGFR2--BICC1")

fusion_freq <- filtered_fusion_res |> 
    select(Fusion, Sample) |> 
    distinct() |> 
    group_by(Fusion) |> 
    tally() |> 
    arrange(desc(n))

write_xlsx(
    list(
        filtered_fusion_report = filtered_fusion_res,
        filtered_fusion_frequency = fusion_freq
    ),
    here("data/RNA/DFSP/Reports/fusion_report_cohort_filtered.xlsx")
)

qc_metrics_summary <- qc_metrics_tbl |> 
    left_join(
        pdgfd_pdgfb_info,
        by = "Sample"
    ) |>
    relocate(Fusion, .after = is_Pass_QC) |> 
    mutate(
        `is_PDGFD|PDGFB` = if_else(
            is.na(Fusion),
            "No",
            "Yes"
        ),
        .after = is_Pass_QC
    )

write_xlsx(
    qc_metrics_summary,
    here("data/RNA/DFSP/Reports/qc_metrics_cohort_summary.xlsx")
)
