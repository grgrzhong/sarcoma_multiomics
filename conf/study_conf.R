study_colors <- list(
    cna_class = c(
        "HOMDEL" = "#053061",
        "DEL" = "#4393c3", 
        "GAIN" = "#f4a582",
        "AMP" = "#d6604d",
        "MULTI" = "#b2182b"
    ),
    FST.Group = c(
        "FS-DFSP" = "#E31A1C", 
        "Post-FST" = "#FF7F00", 
        "Pre-FST" = "#1F78B4", 
        "U-DFSP" = "#33A02C"
    ),
    Metastasis = c(
        "No" = "#A6CEE3", 
        "Yes" = "#FB9A99"
    )
)

group_comparisons <- list(
    list(
        group1 = "U-DFSP",
        group2 = "Pre-FST",
        name = "U-DFSP_vs_Pre-FST"
    ),
    list(
        group1 = "Pre-FST",
        group2 = "Post-FST",
        name = "Pre-FST_vs_Post-FST"
    ),
    list(
        group1 = "Post-FST",
        group2 = "FS-DFSP",
        name = "Post-FST_vs_FS-DFSP"
    ),
    list(
        group1 = c("U-DFSP", "Pre-FST"),
        group2 = c("Post-FST", "FS-DFSP"),
        name = "U-DFSP+Pre-FST_vs_Post-FST+FS-DFSP"
    )
)
