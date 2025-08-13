study_colors <- list(
    cna_class = c(
        "MULTI" = "#b2182b",
        "AMP" = "#d6604d",
        "GAIN" = "#f4a582",
        "HOMDEL" = "#053061",
        "DEL" = "#4393c3" 
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
    HRD.cat = c(
        "High" = "#FFB6C1",
        "Low" = "darkgray"
    ),
    MSI.cat = c(
        "High" = "#ff7f0e",
        "Low" = "darkgray"
    )
)

group_comparisons <- list(
    list(
        group1 = "U-DFSP",
        group2 = "Pre-FST",
        name = "U-DFSP_vs_Pre-FST"
    ),
    list(
        group1 = "U-DFSP",
        group2 = "Post-FST",
        name = "U-DFSP_vs_Post-FST"
    ),
    list(
        group1 = "U-DFSP",
        group2 = "FS-DFSP",
        name = "U-DFSP_vs_FS-DFSP"
    ),
    list(
        group1 = "Pre-FST",
        group2 = "Post-FST",
        name = "Pre-FST_vs_Post-FST"
    ),
    list(
        group1 = "Pre-FST",
        group2 = "FS-DFSP",
        name = "Pre-FST_vs_FS-DFSP"
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

FST.Group <- c(
    "U-DFSP",
    "Pre-FST",
    "Post-FST",
    "FS-DFSP"
)

