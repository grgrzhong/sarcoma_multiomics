##############################################################################
## Description: Configuration settings for the study
## Author:  Zhong Guorui
## Created: 2025-02-17
##############################################################################
## load customized functions
source(here::here("scripts/lib/study_lib.R"))

## Define colors
study_colors <- list(
    Alteration = c(
        "AMP" = "#d6604d",
        "GAIN" = "#f4a582",
        "HOMDEL" = "#053061",
        "DEL" = "#4393c3",
        "SNV" = "#4c9340",
        "substitution" = "#bf8125",
        "insertion" = "#b51e2b",
        "deletion" = "#66011b"
    ),
    cna_class = c(
        "MULTI" = "#b2182b",
        "AMP" = "#d6604d",
        "GAIN" = "#f4a582",
        "HOMDEL" = "#053061",
        "DEL" = "#4393c3" 
    ),
    gistic_cn_class = c(
        "Amp" = "#d6604d",
        "Del" = "#053061"
    ),
    FST.Group = c(
        "FS-DFSP" = "#E31A1C", 
        "Post-FST" = "#FF7F00", 
        "Pre-FST" = "#1F78B4", 
        "U-DFSP" = "#33A02C"
    ),
    FST.Type = c(
        FST = "#ea5b67",
        Classic = "darkgray"
    ),
    Metastasis = c(
        "No" = "#A6CEE3", 
        "Yes" = "#FB9A99"
    ),
    Specimen.Nature = c(
        "Metastasis" = "#CD5C5C",
        "Recurrence" = "#9370DB",
        "Residual" = "#DEB887",
        "Primary" = "gray"
    ),
    Meth.Subtype.Main.2Class = c(
        c(
            Meth1 = "#1f77b4", 
            Meth2 = "#ff7f0e"
        )
    ),
    Molecular.subtype = c(
        "PDGFB" = "white",
        "PDGFD" = "darkgray"
    ),
    Histology.subtype = c(
        "Classic" = "#8FA3B0",
        "FS" = "#C4A5C7",
        "Pigmented" = "#D4B896",
        "Myxoid" = "#B89B9B"
    ),
    HRD.cat = c(
        "High" = "#FFB6C1",
        "Low" = "darkgray"
    ),
    MSI.cat = c(
        "High" = "#ff7f0e",
        "Low" = "darkgray"
    )
)

FST_group_comparisons <- list(
    `U-DFSP_vs_Pre-FST` = list(
        group1 = "U-DFSP",
        group2 = "Pre-FST"
    ),
    `U-DFSP_vs_Post-FST` = list(
        group1 = "U-DFSP",
        group2 = "Post-FST"
    ),
    `U-DFSP_vs_FS-DFSP` = list(
        group1 = "U-DFSP",
        group2 = "FS-DFSP"
    ),
    `Pre-FST_vs_Post-FST` = list(
        group1 = "Pre-FST",
        group2 = "Post-FST"
    ),
    `Pre-FST_vs_FS-DFSP` = list(
        group1 = "Pre-FST",
        group2 = "FS-DFSP"
    ),
    `Post-FST_vs_FS-DFSP` = list(
        group1 = "Post-FST",
        group2 = "FS-DFSP"
    ),
    `U-DFSP+Pre-FST_vs_Post-FST+FS-DFSP` = list(
        group1 = c("U-DFSP", "Pre-FST"),
        group2 = c("Post-FST", "FS-DFSP")
    )
)

group_comparisons <- list(
    ## >>> To find genomic driver event for FST and metastasis
    ## >>> specifically to determine the “FST-selected genetic features”  and 
    ## >>> “metastases-selected genetic features” 
    ## "---------------------------------------------------------------------"
    ## To find early genetic event that may occur in classic region of cases 
    ## having FST
    ## "---------------------------------------------------------------------"
    `U-DFSP_vs_Pre-FST` = list(
        group1 = "U-DFSP",
        group2 = "Pre-FST"
    ),

    ## "---------------------------------------------------------------------"
    ## Unpaired primary FST (FS-DFSP) vs. Post-FST
    ## To find late genetic event in FST
    ## "---------------------------------------------------------------------"
    `Post-FST_vs_FS-DFSP` = list(
        group1 = "Post-FST",
        group2 = "FS-DFSP"
    ),

    ## "---------------------------------------------------------------------"
    ## Unpaired primary FST (FS-DFSP) and unpaired classic samples (U-DFSP)
    ## Compare the prevalence of different mutations and CNV in classic and FST regions, using Fisher’s exact tests
    ## This is **much more straight forward to implement, and can be used for the paired samples above for initial evaluation**.
    ## The list of genes can be compared with paired samples
    ## "---------------------------------------------------------------------"
    `U-DFSP_vs_FS-DFSP` = list(
        group1 = "U-DFSP",
        group2 = "FS-DFSP"
    ),

    ## "---------------------------------------------------------------------"
    ## Paired primary and metastatic samples
    ## "---------------------------------------------------------------------"
    Primary_vs_Metastasis = list(
        group1 = "Primary",
        group2 = "Metastasis"
    ),
    
    ## "---------------------------------------------------------------------"
    ## unpaired primary and metastatic samples in FS-DFSP
    ## "---------------------------------------------------------------------"
    Unpaired_Primary_vs_Metastasis = list(
        group1 = "Primary",
        group2 = "Metastasis"
    ),

    ## "---------------------------------------------------------------------"
    ## Paired classic and FST samples (Pre-FST and Post-FST)
    ## "---------------------------------------------------------------------"
    `Pre-FST_vs_Post-FST` = list(
        group1 = "Pre-FST",
        group2 = "Post-FST"
    )
)

FST.Group <- c(
    "U-DFSP",
    "Pre-FST",
    "Post-FST",
    "FS-DFSP"
)

HRD_genes <- c(
    "BRCA1", "BRCA2", "PALB2", "BRIP1", "RAD51C", "RAD51D", 
    "BARD1", "ATM", "CHEK2", "NBN", "RAD50", "MRE11A",
    "TP53", "PTEN", "RB1", "CDK12", "ARID1A"
)

MMR_genes <- c(
    "MLH1", "MSH2", "MSH6", "PMS2",
    "MSH3", "PMS1", "MLH3", "PCNA",
    "RFC1", "RPA1", "POLD1", "POLE"
)

## WES Capture size
capture_size <- 34 

## Parallel processing settings
mc_cores <- parallel::detectCores() - 10
mc_cores <- 16