library(purrr)
library(tximport)
library(DESeq2)
library(Rsubread)
library(openxlsx)
library(rtracklayer)

## Prepare GTF files with only protein-coding genes
# Load the GTF file
gtf <- import("/mnt/f/Reference/Gencode/gencode.hg38.v44/gencode.v44.primary_assembly.annotation.gtf")

# Filter for protein-coding genes
# Check if 'gene_biotype' or 'gene_type' exists in the GTF
if ("gene_biotype" %in% colnames(mcols(gtf))) {
    protein_coding_gtf <- gtf[gtf$gene_biotype == "protein_coding", ]
} else if ("gene_type" %in% colnames(mcols(gtf))) {
    protein_coding_gtf <- gtf[gtf$gene_type == "protein_coding", ]
} else {
    stop("No 'gene_biotype' or 'gene_type' field found in the GTF file.")
}

# Export the filtered GTF file
export(protein_coding_gtf, "/mnt/f/Reference/Gencode/gencode.hg38.v44/gencode.v44.protein_coding.annotation.gtf")


# Prepare expression matrix from featureCounts
basePath <- "/mnt/m/RNA-seq/DFSP/Output"
GTF_ANNO <- "D:/Reference/Gencode/gencode.v36.primary_assembly.annotation.gtf"
GTF_ANNO <- "/mnt/f/Reference/Gencode/gencode.hg38.v44/gencode.v44.protein_coding.annotation.gtf"
GTF_ANNO <- "/mnt/f/Reference/Gencode/gencode.hg38.v44/gencode.v44.primary_assembly.annotation.gtf"

features <- featureCounts(
    files = paste(basePath, list.files(basePath, pattern = "Aligned.sortedByCoord.out.bam$", recursive = TRUE), sep = "/"),
    annot.ext = GTF_ANNO, 
    isGTFAnnotationFile = TRUE, 
    GTF.attrType = "gene_name",
    useMetaFeatures = TRUE, 
    strandSpecific = 2, 
    isPairedEnd = TRUE, 
    countMultiMappingReads = FALSE,
    requireBothEndsMapped = TRUE, 
    nthreads = 16
)
# GTF.attrType: gene_id; transcript_id; exon_number; gene_name; gene_biotype; transcript_name; exon_id;




features <- featureCounts(
    files = paste(basePath, list.files(basePath, pattern = "Aligned.sortedByCoord.out.bam$", recursive = TRUE), sep = "/"),
    annot.ext = GTF_ANNO, isGTFAnnotationFile = TRUE, GTF.attrType = "exon_id",
    useMetaFeatures = TRUE, strandSpecific = 2, isPairedEnd = TRUE, countMultiMappingReads = FALSE,
    requireBothEndsMapped = TRUE, nthreads = 8
)

counts <- features$counts
colnames(counts) <- unlist(map(strsplit(colnames(counts), split = "/"), 5))

anno <- features$annotation

SampleSheet <- read.xlsx("Y:/RNA-seq/STUMP/SampleSheet_STUMP.xlsx")
SampleSheet$Age.Cat <- ifelse(SampleSheet$Age <= 40, "<=40", ">40")
SampleSheet$Size.Cat <- ifelse(SampleSheet$Size < 10, "<10", ">=10")

# Check case order is the same
table(SampleSheet$Sample_Name == colnames(counts))
SampleSheet <- SampleSheet[order(SampleSheet$Case.No), ]
rownames(SampleSheet) <- SampleSheet$Sample_Name # Assign rownames of SampleSheet for heatmap plot

## design = ~1 if no model
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = SampleSheet,
    design = ~Exp.Subtype
)


### pre-filtering
nrow(dds)
# keep only rows that have at least 10 reads total in 10% samples
keep <- rowSums(counts(dds) >= 10) > round(ncol(counts) * 0.1)
dds <- dds[keep, ]
# at least X samples with a count of 10 or more, where X can be chosen as the sample size of the smallest group of samples
# se_star <- se_star[rowSums(counts(se_star) >= 10) >= 4,]


# Set reference group
dds$Exp.Subtype <- relevel(dds$Exp.Subtype, ref = "Cluster1")
dds$Recurrence <- relevel(dds$Recurrence, ref = "No")

### Discover differential expressed genes
dds2 <- DESeq(dds)

###############################################################################
# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts <- log2(counts(dds2, normalized = TRUE) + 1)

fpkm <- fpkm(se_star2, robust = FALSE)

# tsne plot
library(M3C)
sds <- apply(norm_counts, 1, sd)
norm_counts_sort <- norm_counts[rev(order(sds))[1:1000], ]
tsne(norm_counts_sort,
    labels = SampleSheet$Exp.Subtype,
    perplex = 5,
    dotsize = 4,
    legendtextsize = 20
)



########################################
library(rtracklayer)

# Load the GTF file
gtf <- import(GTF_ANNO)

# Extract protein-coding gene IDs
protein_coding_genes <- unique(gtf$gene_name[gtf$gene_type == "protein_coding"])

# Filter the counts
filtered_counts <- counts$counts[rownames(counts$counts) %in% protein_coding_genes, ]

de_sig_filtered <- de_sig[rownames(de_sig) %in% protein_coding_genes & !grepl("^ENSG", rownames(de_sig)), ]
de_sig_filtered <- de_sig_filtered[order(de_sig_filtered$padj), ]



################################################################################

### Output results from dds object
resultsNames(dds2)
de <- results(
    object = dds2,
    name = "Exp.Subtype_Cluster2_vs_Cluster1"
)

# Shrinkage for better visualization of results
de_shrink <- lfcShrink(
    dds = dds2,
    coef = "Exp.Subtype_Cluster2_vs_Cluster1",
    type = "apeglm"
)

# Select significantly differentially expressed genes for heatmap
threshold <- (abs(de_shrink$log2FoldChange) > 1) & (de_shrink$padj < 0.1)

de_shrink$threshold <- threshold

de_sig <- data.frame(subset(de_shrink, threshold == TRUE))

de_sig_filtered <- de_sig[!grepl("^ENSG", rownames(de_sig)) & !grepl("-AS1$", rownames(de_sig)), ]

norm_sig <- norm_counts[rownames(de_sig), ]

# Select the top n number DEG for plotting if necessary
de_sig_top <- de_sig[head(order(de_sig$padj, decreasing = TRUE), 101), ]
norm_sig_top <- norm_counts[rownames(de_sig_top), ]

de_sig_top %in% de_sig

# Plot heatmap

anno <- c(
    "Recurrence", "Age.Cat", "Size.Cat",
    "Exp.Subtype"
)

library(pheatmap)
heatmap <- pheatmap(norm_sig_top,
    cluster_rows = T, show_rownames = T,
    labels_col = SampleSheet$Sample_Name, show_colnames = T,
    annotation = SampleSheet[, anno], cutree_cols = 2,
    border_color = NA, fontsize = 10, scale = "row",
    fontsize_row = 4, fontsize_col = 10, height = 20
)

# The resulting object contains the clustering dendrograms
column_dendrogram <- heatmap$tree_col

# You can cut the dendrogram to create clusters/subgroups
col_clusters <- cutree(column_dendrogram, k = 2)

SampleSheet$Exp.Subtype <- factor(col_clusters, labels = c("Cluster1", "Cluster2"))
## The most important one is scale="row", in which Z-scores are plotted,
## rather than the actual normalized count value. Z-scores are computed
## on a gene-by-gene basis by subtracting the mean and then dividing by
## the standard deviation. The Z-scores are computed after the clustering, so
## that it only affects the graphical aesthetics and the color visualization is improved.

summary(de_shrink, alpha = 0.1)

#################################################################################
### Generate distance matrix for heatmap and PCA plot
vsd <- vst(dds2)
rld <- rlog(se_star2, blind = TRUE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

# load libraries pheatmap to create the heatmap plot
library(pheatmap)

# calculate between-sample distance matrix
sampleDistMatrix <- as.matrix(dist(t(assay(rld))))
sampleDistMatrix <- assay(rld)[topVarGenes, ]
# Generate heatmap on distance matrix
pheatmap(sampleDistMatrix, labels_row = SampleSheet$Diagnosis, labels_col = SampleSheet$Diagnosis)

# Generate heatmap with col as subgroups and row as genes
vsd <- assay(vst(se_star2))
mads <- apply(vsd, 1, mad)


# Find the number of optimal classes
library(ConsensusClusterPlus)
results <- ConsensusClusterPlus(mad2k, # normalized expression matrix
    maxK = 6, # maximum k to make clusters
    reps = 50, # resampling
    pItem = 0.8, # item resampling
    pFeature = 1, # gene resampling
    title = "geneExp", # output title
    clusterAlg = "pam", # clustering algorithm
    distance = "spearman", # distance measure
    seed = 1262118388.71279, # random number
    plot = "pdf"
) # output file format

# Extract information about k classes
cck <- results[[k]]
cckClass <- data.frame(cck$consensusClass) # Class number of each sample

# Plot PCA graph
plotPCA(
    object = rld,
    intgroup = "LTK"
)



# Count plot
library(ggplot2)
library(ggsignif)
library(viridis)
topGene <- rownames(de_shrink)[which.min(de_shrink$padj)]
topGene2 <- rownames(de_shrink)[which.max(de_shrink$log2FoldChange)]

gene <- c("CCND1")
d <- plotCounts(se_star2, gene = gene, intgroup = "Subtype", returnData = TRUE)
p_value <- se_star

# Boxplot
ggplot(d, aes(x = Subtype, y = count, fill = Subtype)) +
    geom_boxplot(show.legend = TRUE) +
    labs(x = "Molecular subtype", y = "Normalized counts", title = paste(gene, "gene expression")) +
    geom_signif(
        comparisons = list(c("PDGFB", "PDGFD")),
        map_signif_level = TRUE, textsize = 6, step_increase = 0.1
    ) +
    theme_classic() +
    theme(
        legend.text = element_text(size = 20), # Adjust legend text size
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.key.size = unit(2.5, "lines"),
        axis.title.x = element_text(size = 18), # Adjust x-axis label size
        axis.title.y = element_text(size = 16), # Adjust y-axis label size
        axis.text.x = element_text(size = 16), # Adjust x-axis text size
        axis.text.y = element_text(size = 14)
    ) # Adjust y-axis text size

# Violin plot
ggplot(d, aes(x = Subtype, y = count, fill = Subtype)) +
    geom_violin(trim = FALSE, width = 1.2) +
    geom_boxplot(show.legend = TRUE, width = 0.1, fill = "white") +
    # scale_fill_viridis(discrete = TRUE) +
    labs(x = "DFSP molecular subtype", y = "Normalized counts", title = "PDGFB gene expression") +
    geom_signif(
        comparisons = list(c("PDGFB", "PDGFD")),
        map_signif_level = TRUE, textsize = 6
    ) +
    theme_classic() +
    theme(
        legend.text = element_text(size = 22), # Adjust legend text size
        legend.title = element_text(size = 24),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.key.size = unit(2.5, "lines"),
        axis.title.x = element_text(size = 16), # Adjust x-axis label size
        axis.title.y = element_text(size = 16), # Adjust y-axis label size
        axis.text.x = element_text(size = 14), # Adjust x-axis text size
        axis.text.y = element_text(size = 14)
    ) + # Adjust y-axis text size
    annotate("text", x = 1.5, y = 1000, label = paste("padj =", round(de_shrink$padj[rownames(de_shrink) == "PDGFB"], 3)), color = "black", size = 6)



topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 25)
mat <- assay(rld)[topVarGenes, ]
reference <- colnames(mat) %in% SampleSheet[SampleSheet$Subtype == "PDGFB", ]$Case.No
test <- colnames(mat) %in% SampleSheet[SampleSheet$Subtype == "PDGFD", ]$Case.No
mat <- mat - rowMeans(mat[, reference])
mat <- mat - rowMeans(mat)
pheatmap(mat,
    labels_row = rownames(se_star2)[topVarGenes], labels_col = SampleSheet$Subtype,
    cluster_cols = TRUE, cluster_rows = TRUE
)


# Volcano plot
library(EnhancedVolcano)
EnhancedVolcano(de_sig_filtered_reprocessed,
    lab = NA,
    x = "log2FoldChange",
    y = "padj",
    legendPosition = "right",
    drawConnectors = TRUE,
    labSize = 3,
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    legendIconSize = 2,
    legendLabSize = 12,
    xlim = c(-7, 7),
    ylim = c(0, 5),
    title = "RNA cluster 2 vs cluster 1"
)


# GSEA
library(clusterProfiler)
library(DOSE)
# Prepare input
# reading in data from deseq2 - differential expression output with log2fc

## The input for GSEA is a named factor of log2fc, sorted in descreasing value
gene_list <- de_sig_filtered[order(de_sig_filtered$log2FoldChange, decreasing = TRUE), ]

gene_list_input <- gene_list$log2FoldChange # extract the list to a vector
names(gene_list_input) <- rownames(gene_list) # name the vector
# sort the vector in decreasing order (required for clusterProfiler)
gene_list_input <- sort(gene_list_input, decreasing = TRUE)

# Load database
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

# check for keytypes of database
keytypes(org.Hs.eg.db)

# GO, gene ontology; BP, biological process; MF, molecular function; CC, cellular component.
gse <- gseGO(
    geneList = gene_list_input,
    ont = "ALL",
    keyType = "SYMBOL",
    minGSSize = 3,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    OrgDb = organism,
    pAdjustMethod = "BH"
)
gseaplot(gse,
    geneSetID = "PDGFB",
    by = "runningScore"
)

dotplot(gse)
goplot(gse)
barplot(gse)

gseaplot(gse, geneSetID = 2, title = gse$Description[2])
gseaplot2(gse, geneSetID = 2, title = gse$Description[2])

# https://www.spandidos-publications.com/10.3892/etm.2018.6884

# The top ranked GO terms according to gene count.
# ??Gene count?? is the number of genes enriched in a GO term.
# ??Gene ratio?? is the percentage of total DEGs in the given GO term
# (only input genes with at least one GO term annotation were included in the calculation).

eg <- bitr(de_symbols$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = FALSE)
de_symbols2 <- merge(de_symbols, eg$ENTREZID, by = eg$ENTREZID)

# GO classification
ggo <- groupGO(
    gene = eg$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    level = 4,
    readable = TRUE
)

head(ggo)

# GO Over-representation analysis
ego <- enrichGO(
    gene = eg$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)
head(ego)
goplot(ego)

dotplot(ego)

gene <- names(geneList)[abs(geneList) > 2]


gene_list
ggo <- groupGO(
    gene = gene_symbol,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "CC",
    level = 3,
    readable = TRUE
)

head(ggo)

## KEGG analysis
install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method", "auto")


kk2 <- gseKEGG(
    geneList = geneList,
    organism = "hsa",
    keyType = "kegg",
    minGSSize = 3,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    pAdjustMethod = "BH"
)

eg2np <- bitr_kegg(de_symbols$ENSG, fromType = "ncbi-geneid", toType = "kegg", organism = "hsa")


gseaplot(kk2,
    geneSetID = 1:3,
    title = kk2$Description[1]
)
library(enrichplot)
gseaplot2(kk2, geneSetID = 1:5)

#############################################

bam_list <- c(
    "D:/RNA-seq/Output/22BX18526-5/Aligned.sortedByCoord.out.bam",
    "D:/RNA-seq/Output/21B18984-1-1/Aligned.sortedByCoord.out.bam",
    "D:/RNA-seq/Output/22BX3762-1/Aligned.sortedByCoord.out.bam"
)
counts <- featureCounts(bam_list,
    annot.ext = "D:/Reference/Gencode/gencode.v36.primary_assembly.annotation.gtf",
    isGTFAnnotationFile = TRUE,
    useMetaFeatures = FALSE,
    strandSpecific = 2,
    isPairedEnd = TRUE,
    countMultiMappingReads = FALSE,
    requireBothEndsMapped = TRUE,
    nthreads = 8
)
