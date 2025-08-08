
## Download data
library(TCGAbiolinks)
library(easyTCGA)
getsnvmaf("TCGA-COAD")

library(ComplexHeatmap)

## Read the Complexheatmap internal data
mat = read.table(
    system.file(
        "extdata", 
        package = "ComplexHeatmap", 
        "tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"
    ), 
    header = TRUE, 
    stringsAsFactors = FALSE, 
    sep = "\t"
)

dim(mat)
mat[1:5, 1:4]

mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]

mat = t(as.matrix(mat))
mat[1:5, 1:4]
dim(mat)

## Create oncoprint
col = c(
    "HOMDEL" = "blue", 
    "AMP" = "red", 
    "MUT" = "#008000"
)

alter_fun = list(
    background = alter_graphic("rect", fill = "#CCCCCC"),   
    HOMDEL = alter_graphic("rect", fill = col["HOMDEL"]),
    AMP = alter_graphic("rect", fill = col["AMP"]),
    MUT = alter_graphic("rect", height = 0.33, fill = col["MUT"])
)

heatmap_legend_param = list(
    title = "Alternations", 
    at = c("HOMDEL", "AMP", "MUT"), 
    labels = c("Deep deletion", "Amplification", "Mutation")
)

oncoPrint(
    mat,
    alter_fun = alter_fun,
    col = col, 
    heatmap_legend_param = heatmap_legend_param
)

load(
    here("output_snv/TCGA-COAD_maf.rdata")
)
snv


all <- read_tsv(
    here("data/wes/GISTIC2/somatic_matched/all_tumors/all_lesions.conf_99.txt")
)

colnames(all)[grepl("DFSP", colnames(all))]

amp <- read_tsv(
    here("data/wes/GISTIC2/somatic_matched/all_tumors/amp_genes.conf_99.txt")
)

del <- read_tsv(
    here("data/wes/GISTIC2/somatic_matched/all_tumors/del_genes.conf_99.txt")
)

cytobands <- map_chr(
    sig_cytobands,
    ~ str_split(.x, "_")[[1]][2]
)

table(cytobands %in% all$Descriptor)
view(all)

length(unique(all$Descriptor))
