
library(here)
library(tidyverse)
library(readxl)

## QC metrics
qc_metrics_file <- "data/RNA/DFSP/Reports/qc_metrics_cohort.xlsx"

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
fusion_report_file <- "data/RNA/DFSP/Reports/fusion_report_cohort.xlsx"
fusion_report_data <- read_excel(
    fusion_report_file,
    sheet = "fusion_report"
)

write_csv(
    fusion_report_data,
    here("data/RNA/DFSP/Reports/fusion_report_data.csv")
)

colnames(fusion_report_data)
table(fusion_report_data$`arriba.reading-frame`)

filtered_fusion_res <- fusion_report_data |> 
    filter(Number_of_callers > 1) |> 
    mutate(
        fusion_QC_confidence = if_else(
            starfusion.junction_reads >=3 &
            starfusion.spanning_reads >=5 &
            (`arriba.reading-frame` %in% c("in-frame", "stop-codon")),
            "High",
            "Low"
        ),
        .after = Fusion
    ) |> 
    arrange(fusion_QC_confidence)

known_fusion_res <- fusion_report_data |> 
    filter()
    select(starfusion.junction_reads, starfusion.spanning_reads)
