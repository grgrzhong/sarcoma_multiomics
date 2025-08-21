
library(ComplexHeatmap)


mat <- read.table(
    system.file("extdata",
        package = "ComplexHeatmap",
        "tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"
    ),
    header = TRUE, stringsAsFactors = FALSE, sep = "\t"
)

mat[is.na(mat)] <- ""
rownames(mat) <- mat[, 1]
mat <- mat[, -1]
mat <- mat[, -ncol(mat)]
mat <- t(as.matrix(mat))
mat[1:3, 1:3]

col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000")

alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # big blue
    HOMDEL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["HOMDEL"], col = NA))
    },
    # big red
    AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["AMP"], col = NA))
    },
    # small green
    MUT = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
            gp = gpar(fill = col["MUT"], col = NA))
    }
)

# just for demonstration
alter_fun = list(
    background = alter_graphic("rect", fill = "#CCCCCC"),   
    HOMDEL = alter_graphic("rect", fill = col["HOMDEL"]),
    AMP = alter_graphic("rect", fill = col["AMP"]),
    MUT = alter_graphic("rect", height = 0.33, fill = col["MUT"])
)

column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling"
heatmap_legend_param = list(
    title = "Alternations", 
    at = c("HOMDEL", "AMP", "MUT"), 
    labels = c("Deep deletion", "Amplification", "Mutation")
)

oncoPrint(
    mat,
    alter_fun = alter_fun,
    col = col,
    column_title = column_title,
    heatmap_legend_param = heatmap_legend_param
)

## List of heatmap
set.seed(123)
mat1 = matrix(rnorm(80, 2), 8, 10)
mat1 = rbind(mat1, matrix(rnorm(40, -2), 4, 10))
rownames(mat1) = paste0("R", 1:12)
colnames(mat1) = paste0("C", 1:10)

mat2 = matrix(runif(60, max = 3, min = 1), 6, 10)
mat2 = rbind(mat2, matrix(runif(60, max = 2, min = 0), 6, 10))
rownames(mat2) = paste0("R", 1:12)
colnames(mat2) = paste0("C", 1:10)

le = sample(letters[1:3], 12, replace = TRUE)
names(le) = paste0("R", 1:12)

ind = sample(12, 12)
mat1 = mat1[ind, ]
mat2 = mat2[ind, ]
le = le[ind]

ht1 = Heatmap(mat1, name = "rnorm")
ht2 = Heatmap(mat2, name = "runif")
ht3 = Heatmap(le, name = "letters")

ht_list = ht1 + ht2 + ht3
class(ht_list)

# code only for demonstration
ht1 + ht_list
ht_list + ht1
ht_list + ht_list

col_rnorm = colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
col_runif = colorRamp2(c(0, 3), c("white", "orange"))
col_letters = c("a" = "pink", "b" = "purple", "c" = "blue")

ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm,
    row_title = "Heatmap 1", column_title = "Heatmap 1")
ht2 = Heatmap(mat2, name = "runif", col = col_runif,
    row_title = "Heatmap 2", column_title = "Heatmap 2")
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht_list = ht1 + ht2 + ht3

draw(
    ht_list, 
    row_title = "Three heatmaps, row title", 
    row_title_gp = gpar(col = "red"),
    column_title = "Three heatmaps, column title", 
    column_title_gp = gpar(fontsize = 16),
    ht_gap = unit(1, "cm")
)

ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, row_km = 2)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht2 + ht1 + ht3 # ht2 is the main heatmap and row_km in ht1 is ignored

ht_list = ht2 + ht1 + ht3
draw(ht_list, main_heatmap = "rnorm")

ht_list = ht2 + ht1 + ht3
draw(ht_list, main_heatmap = "rnorm", row_dend_side = "right", row_sub_title_side = "left")

ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)

ht_list = ht1 + ht2 + ht3
draw(ht_list, auto_adjust = FALSE)

ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, row_km = 2, cluster_rows = FALSE)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht_list = ht1 + ht2 + ht3
draw(ht_list, row_km = 1, row_split = le, cluster_rows = TRUE)

ha1 = HeatmapAnnotation(foo1 = 1:10, annotation_name_side = "left")
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, top_annotation = ha1)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht1 + ht2 + ht3

ha1 = HeatmapAnnotation(foo1 = 1:10, bar1 = anno_points(1:10), 
    annotation_name_side = "left")
ha2 = HeatmapAnnotation(bar2 = anno_barplot(1:10))
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, top_annotation = ha1)
ht2 = Heatmap(mat2, name = "runif", col = col_runif, top_annotation = ha2)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht_list = ht1 + ht2 + ht3
draw(ht_list, ht_gap = unit(c(6, 2), "mm"))
