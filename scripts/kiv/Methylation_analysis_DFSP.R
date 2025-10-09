library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)
library(reshape2)
library(pheatmap)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

##Explanation of different files
# SampleSheet_DFSP.csv: Original sample sheet
# Beta_DFSP_all.csv: beta value of all DFSP
# phenoData_DFSP_all.csv: phenoData of all DFSP
# phenoData_DFSP_97pt: phenoData of 97 cases, with meth subgroups, CN group, Amp/Del peaks and IC %
# - Main = TRUE, Specimen.Nature NOT %in% c("Metastasis", "Recurrence (Post-imatinib)")
# - Samples excluded by DetP filtering
# phenoData_DFSP_96pt_Meth_immune: phenoData of 96 cases, with meth subgroups and IC %
# - Manually excluded "DFSP-356-T"
# beta_DFSP_97pt.csv: beta values of 97 cases
# beta_DFSP_96pt.csv: beta values of 96 cases

# DFSP_CNVseg_PDGFD.txt: seg file for PDGFD
# DFSP_CNVseg_PDGFBD.txt: seg file for PDGFB
# DFSP_CNVseg.txt: seg file for all DFSP

#phenoData <- read.csv("Y:/Methylation_profiling/DFSP/phenoData_DFSP_97pt_Meth_immune_CNV.csv", row.names = 1)
#phenoData <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_97pt_Meth_immune_CNV.csv", row.names = 1)



#Load these pre-calculated beta matrix and phenoData
#phenoData <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_97pts.csv", row.names = 1)
#phenoData <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_all.csv")
#beta <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/beta_DFSP_97pt.csv", row.names = 1)
#

beta_all <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/Beta_DFSP_all_tumour.csv", row.names = 1)
colnames(beta_all) <- gsub("X","",colnames(beta_all))
phenoData_all <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_allTumour_CNV_Meth.csv", row.names = 1) 

table(phenoData_all$IDAT==colnames(beta_all))
phenoData <- phenoData[match(colnames(beta),phenoData$IDAT),]

#Subset certain cases
phenoData <- phenoData_all[phenoData_all$FST.subgroup %in% c("FS-DFSP", "U-DFSP"),]
beta <- beta_all[,colnames(beta_all) %in% rownames (phenoData)]

table(phenoData$IDAT==colnames(beta))


RGSet_v2 <- read.metharray.exp(target = SampleSheet, verbose = TRUE, force = TRUE)

phenoData_v2 <- pData(RGSet_v2)
phenoData <- data.frame(phenoData_v2)

qcReport(RGSet_v2, sampNames=phenoData$Sample.ID, pdf="D:/Sarcoma/Methylation_profiling/DFSP/qcReport.pdf")

### V2
detP <- detectionP(RGSet_v2)

#Remove samples that have average detection p-value <0.05
#keep <- colMeans(detP) <0.05

#Remove samples that have probes detection p-value >0.01 in >10% of probes
remove <- colSums(detP > 0.01) > (nrow(RGSet_v2) * 0.11)
RGSet_v2_filtered <- RGSet_v2[,!(remove)]
print(paste("Removed", as.character(table(remove)[2]), "samples with detection p-value >0.01 in >10% probes"))
print(phenoData_v2$Sample.ID[phenoData_v2$IDAT %in% colnames(RGSet_v2)[remove]]) ##Removed "DFSP-072-T"
phenoData <- pData(RGSet_v2_filtered)
phenoData <- data.frame(phenoData)

#Between-array Functional normalization
MSet_v2 <- preprocessFunnorm(RGSet_v2_filtered)
any(duplicated(rownames(MSet_v2)))

#Remove probes with detection p-value >0.01 in >1% samples
detP <- detP[match(featureNames(MSet_v2),rownames(detP)),]
remove <- rowSums(detP > 0.01) > (ncol(MSet_v2) * 0.01)
MSet_filtered_v2 <- MSet_v2[!(remove),]
print(paste("Removed", as.character(table(remove)[2]), "probes with detection p-value >0.01 in >1% samples"))

#Remove probes with SNPs at CpG site
MSet_filtered_v2 <- dropLociWithSnps(MSet_filtered_v2)

#Remove probes on sex chromosome
keep <- !(featureNames(MSet_filtered_v2) %in% annEPICv2$Name[annEPICv2$chr %in% c("chrX","chrY")] )
MSet_filtered_v2 <- MSet_filtered_v2[keep,]
print(paste("Removed", as.character(table(keep)[1]), "probes on sex chromsome"))

#Remove cross-reactive probes
xReactiveProbes <- read.csv(file="D:/Sarcoma/Methylation_profiling/XProbes/48639-non-specific-probes-Illumina450k.csv",
                            stringsAsFactors=FALSE)

keep <- !(strsplit2(featureNames(MSet_filtered_v2), split="_")[,1] %in% xReactiveProbes$TargetID)
print(paste("Removed", as.character(table(keep)[1]), "cross-reacting probes"))
MSet_filtered_v2 <- MSet_filtered_v2[keep,]

beta <- getBeta(MSet_filtered_v2)
MVal <- getM(MSet_filtered_v2)

meth <- getMeth (MSet_filtered_v2)
unmeth <- getUnmeth (MSet_filtered_v2)
M <- log2((meth+100)/(unmeth+100))

##**start analysis here**##
beta_all <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/Beta_DFSP_all_tumour.csv", row.names = 1)
colnames(beta_all) <- gsub("X","",colnames(beta_all))
MVal_all <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/MVal_DFSP_all_tumour.csv", row.names = 1)
colnames(MVal_all) <- gsub("X","",colnames(MVal_all))
phenoData_all <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_allTumour_CNV_Meth.csv", row.names = 1)


#Subset certain cases for analysis
selection <- phenoData_all$Main == "Yes"
selection <- phenoData_all$FST.subgroup %in% c("Pre-FST", "Post-FST")
selection <- phenoData_all$FST.subgroup %in% c("U-DFSP", "FS-DFSP")
selection <- phenoData_all$FST.subgroup %in% c("U-DFSP", "Pre-FST")
selection <- phenoData_all$FST.subgroup %in% c("Post-FST", "FS-DFSP")
selection <- phenoData_all$Sample.ID %in% c("DFSP−047−T−P2","DFSP−047−T−P1",
                                            "DFSP−111−T−P2","DFSP−111−T−P1",
                                            "DFSP−145−T−P2","DFSP−145−T−P1",
                                            "DFSP−188−T−P2","DFSP−188−T−P1",
                                            "DFSP−191−T−P2","DFSP−191−T−P1",
                                            "DFSP−266−T−P2","DFSP−266−T−P1",	
                                            "DFSP−305−T−P2","DFSP−305−T−P1")

beta <- beta_all[,selection]
MVal <- MVal_all[,selection]
phenoData <- phenoData_all[selection,]

#Check the order of beta colnames and IDAT rownames of phenoData is the same 
table(colnames(beta) == phenoData$IDAT) #Expect all are TRUE


#Estimate purity
library(RFpurify)
MSet <- preprocessNoob(RGSet_v2_filtered)
rownames(MSet) <- strsplit2(rownames(MSet), split = "_")[,1]
absolute <- predict_purity(MSet,method="ABSOLUTE")
estimate <- predict_purity(MSet,method="ESTIMATE")


library(M3C)
mads <- apply(beta,1,mad)
beta_sort <- beta[rev(order(mads))[1:8000],]

tsne(mydata = beta_sort, 
     labels = phenoData$FST.subgroup,
     dotsize = 5,
     legendtextsize = 20,
     textlabelsize = 2,
     perplex = 10,
  #   controlscale = TRUE,
   #  scale = 3,
     colvec = c("skyblue","red"),
     text = phenoData$Sample.ID
)
#?Outlier = "DFSP-356-T" ,"DFSP-95-T", "DFSP-348-T"

library(Rtsne)
library(ggplot2)
phenoData$Meth.Subtype.2Class <- as.factor(sample(c("Group1", "Group2", "Group3"), 100, replace = TRUE))
df$factor_shape <- as.factor(sample(c("TypeA", "TypeB"), 100, replace = TRUE))
tsne <- Rtsne(beta_sort, dims = 2, perplexity = 30, verbose = TRUE)
tsne_data <- data.frame(X = tsne$Y[,1], Y = tsne$Y[,2])

plot(tsne_data, col=phenoData$Meth.Subtype.2Class)

ggplot(tsne_data, aes(x = X, y = Y, color = phenoData$Meth.Subtype.2Class, shape = phenoData$FST) +
  geom_point(alpha = 0.5, size = 3) +
  scale_color_manual(values = c("Meth1" = "red", "Meth2" = "blue")) +
  scale_shape_manual(values = c("Yes" = 16, "No" = 17)) +  # shapes: 16 is a filled circle, 17 is a filled triangle
  theme_minimal() +
  labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2")


ggplot()+
    geom_point(data=tsne_data,
               aes(x=X,y=Y,color=phenoData$FST))+
    theme_bw()
  

  
phenoData_paired <- phenoData[phenoData$Histology.Nature == "Paired classic" |
                                phenoData$Histology.Nature == "Paired FST",]
rownames(phenoData_paired) <- phenoData_paired$IDAT
beta_paired <- beta[, colnames(beta) %in% phenoData_paired$IDAT]
mads <- apply(beta_paired,1,mad)
beta_sort <- beta_paired[rev(order(mads))[1:10000],]

phenoData$Histology.Nature.Group[phenoData$Histology.Nature == "Paired classic"] <- "Paired classic" 
phenoData$Histology.Nature.Group[phenoData$Histology.Nature == "Paired FST"] <- "Paired FST" 
phenoData$Histology.Nature.Group[!(phenoData$Histology.Nature == "Paired classic" |
                                     phenoData$Histology.Nature == "Paired FST")] <- "Others" 
phenoData$Histology.Nature.Group[phenoData$Specimen.Class == "Normal"] <- "Normal" 


##Plot heatmap
library(pheatmap)

sampleMatrix <- as.matrix(beta_sort)

anno <- c("Local.recurrence", "Metastasis","Sex","Molecular.subtype", "Specimen.Nature", "Histology.subtype", "Meth.Subtype.Main.2Class")

anno <- c("Sex", "Site.Category.Broad", "Molecular.subtype",
           "Histology.Subtype", "FST.subgroup", "Meth.Subtype.2Class"
                    )

annotation_colors <- list(
  FST.subgroup = c("FS-DFSP" = "tomato", "Post-FST" = "seagreen", "U-DFSP" = "dodgerblue", "Pre-FST" = "orchid"),
  Metastasis = c("No" = "steelblue", "Yes" = "gold"),
  Local.recurrence = c("No" = "mediumseagreen", "Yes" = "slateblue"),
  Molecular.subtype = c("PDGFB" = "firebrick","PDGFD" = "orchid"),
  Site.Category.Broad = c(
    "Trunk" = "magenta", "Extremity" = "coral", "H&N" = "navy"),
  Histology.Subtype = c(
    "Classic" = "#377EB8",
    "FS" = "#E41A1C",
    "Myxoid" = "#4DAF4A",
    "Pigmented" = "#984EA3"
  )
)

heatmap <- pheatmap(sampleMatrix, cluster_rows = T, show_rownames=F,
                    labels_col = phenoData$Sample.ID, show_colnames = T,
                    annotation = phenoData[,anno], cutree_cols = 2,
                    clustering_method = "ward.D2", 
                    annotation_colors = annotation_colors,
                    border_color=NA, fontsize = 10,
                    fontsize_row = 10, fontsize_col = 10, height=30)

# The resulting object contains the clustering dendrograms
column_dendrogram <- heatmap$tree_col

# You can cut the dendrogram to create clusters/subgroups
col_clusters <- cutree(column_dendrogram, k = 2)  # adjust 'k' to the number of clusters/subgroups you want

# You can then add this information back to your original data
phenoData$Meth.subgroup <- factor(col_clusters, labels = c("Meth1", "Meth2", "Meth3", "Meth4"))  # adjust labels as needed
phenoData$Meth.subgroup.WardD <- factor(col_clusters, labels = c("Meth1", "Meth2", "Meth3", "Meth4"))
phenoData$Meth.subgroup.WardD2 <- factor(col_clusters, labels = c("Meth1", "Meth2", "Meth3", "Meth4"))
phenoData$Meth.subgroup.probesel <- factor(col_clusters, labels = c("Meth1", "Meth2", "Meth3", "Meth4"))
phenoData$Meth.Subtype.Main.2Class <- factor(col_clusters, labels = c("Meth1", "Meth2"))
phenoData$Meth.Subtype.Main.3Class <- factor(col_clusters, labels = c("Meth1", "Meth2", "Meth3"))

write.csv(phenoData, "D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_97pts.csv", row.names = TRUE)
write.csv(beta_sort, "D:/Sarcoma/Methylation_profiling/DFSP/beta_probeselection_top8000.csv", row.names = TRUE)

library(ConsensusClusterPlus)
setwd("F:/Documents/Research/DFSP/Data analysis")
results = ConsensusClusterPlus(beta_sort,maxK=10,reps=1000,pItem=0.8,pFeature=1, 
                               title = "DFSP Consensus Cluster 97 pts",
                               clusterAlg="hc",distance="euclidean",
                               seed=1262118388.71279,plot= "pdf")

results

icl = calcICL(results,title= "DFSP Consensus Cluster",plot="pdf")


##**Survival analysis**##
library(survival)
library(survminer)

cox_ph_model <- coxph(Surv(time_to_event, event_status) ~ x1_continuous + x2_categorical, data = your_data)

# Convert survival months into numeric variables
phenoData$OS.time <- as.numeric(phenoData$OS.time)
phenoData$RFS.time <- as.numeric(phenoData$RFS.time)

# Produce survival object
survival_obj_OS <- Surv(phenoData$OS.time, phenoData$OS.status == "1:DECEASED")
survival_obj_RFS <- Surv(phenoData$RFS.time, phenoData$RFS.status== "1:Recurrence")

#Subset to FST
phenoData_FST <- phenoData[phenoData$FST == "Yes",]
survival_obj <- Surv(phenoData_FST$OS.time, phenoData_FST$OS.status == "1:DECEASED")
survival_obj <- Surv(phenoData_FST$RFS.time, phenoData_FST$RFS.status== "1:Recurrence")

km_curves_OS <- survfit(survival_obj_OS ~ Meth.Subtype.Main.2Class, data = phenoData_new_subset)
km_curves_RFS <- survfit(survival_obj_RFS ~ Meth.Subtype.Main.2Class, data = phenoData_new_subset)

break.levels <- c("CN1", "CN2", "CN3", "CN4")
break.levels <- c("PDGFB","PDGFD")
break.levels <- c("Meth1","Meth2")
break.levels <- c("No FST","FST")

surv_plot <- ggsurvplot(km_curves_RFS, data = phenoData, conf.int = FALSE, pval = TRUE,
                        pval.method = TRUE, risk.table = TRUE,
                        xlab = "Time (months)", ylab = "OS probability",
                        legend.labs = break.levels, 
                        legend.title = "CN Subgroup",
                        color = "CN.subgroup4", # set color to grouping variable
                        palette = c("tomato","seagreen","dodgerblue", "violet"), # specify custom color palette
                        font.legend = 15, font.x = 15, font.y = 15, pval.size = 5)

surv_plot <- ggsurvplot(km_curves_OS, data = phenoData_new_subset, conf.int = FALSE, pval = TRUE,
                        pval.method = TRUE, risk.table = TRUE,
                        xlab = "Time (months)", ylab = "OS probability",
                        legend.labs = break.levels, 
                        legend.title = "Meth subgroup",
                       font.legend = 15, font.x = 15, font.y = 15, pval.size = 5)

surv_plot <- ggsurvplot(km_curves_RFS, data = phenoData_new_subset, conf.int = FALSE, pval = TRUE,
                        pval.method = TRUE, risk.table = TRUE,
                        xlab = "Time (months)", ylab = "RFS probability",
                        legend.labs = break.levels, 
                        legend.title = "FST",
                        color = "CN.subgroup", # set color to grouping variable
                        palette = c("tomato","seagreen","dodgerblue"), # specify custom color palette
                        font.legend = 15, font.x = 15, font.y = 15, pval.size = 5)
surv_plot

phenoData$CN.bin[phenoData$CN.subgroup=="CN1"] <- "CN1"
phenoData$CN.bin[phenoData$CN.subgroup %in% c("CN2","CN3")] <- "CN2+3"

phenoData$Amp17q22q[phenoData$AP_7.17q21.32=="Amp"|
                      phenoData$AP_8.17q22 == "Amp"|
                      phenoData$AP_9.17q22=="Amp"|
                      phenoData$AP_10.17q23.3=="Amp"|
                      phenoData$AP_11.17q24.3=="Amp"|
                      phenoData$AP_12.17q25.3=="Amp"|
                      phenoData$AP_13.17q25.3=="Amp"|
                      phenoData$AP_14.22q11.22=="Amp"|
                      phenoData$AP_15.22q12.3=="Amp"|
                      phenoData$AP_16.22q13.1=="Amp"] <- "Amp"
phenoData$Amp17q22q[!(phenoData$AP_7.17q21.32=="Amp"|
                      phenoData$AP_8.17q22 == "Amp"|
                      phenoData$AP_9.17q22=="Amp"|
                      phenoData$AP_10.17q23.3=="Amp"|
                      phenoData$AP_11.17q24.3=="Amp"|
                      phenoData$AP_12.17q25.3=="Amp"|
                      phenoData$AP_13.17q25.3=="Amp"|
                      phenoData$AP_14.22q11.22=="Amp"|
                      phenoData$AP_15.22q12.3=="Amp"|
                      phenoData$AP_16.22q13.1=="Amp")] <- "Non-Amp"


##** limma differential methylation probes analysis **##
# For EPIC_V2, please run these lines to remove replicate and non-cg probes:
# remove non-cg probes
MVal_filtered <- MVal[grepl("^cg", rownames(MVal)),]
# select replicate probes with best sensitivity
MVal_filtered <- DMRcate::rmPosReps(MVal_filtered, filter.strategy="sensitivity")

design <- model.matrix(~FST.subgroup, data=phenoData)

lfit1 <- lmFit(MVal_filtered, design = design)
lfit2 <- eBayes(lfit1) # Stage 1 analysis

summary(decideTests(lfit2, lfc = 0.3))

##**GSEA on CpG using missMethyl**##
top <- topTable(lfit2, number = Inf)
topCpGs<-topTable(lfit2,number=10000)
sigCpGs <- rownames(topCpGs)

gst <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(top), collection="GO",
              array.type = "EPIC_V2",plot.bias=TRUE, 
              genomic.features = c("TSS200", "TSS1500","1stExon"),
              sig.genes = TRUE)


topGSA(gst, n=20)


#** Perform RUV adjustment and fit **#
## The M-Value must be generated:  
#meth <- getMeth (MSet_filtered_v2)
#unmeth <- getUnmeth (MSet_filtered_v2)
#M <- log2((meth+100)/(unmeth+100))
#MVal <- getM(MSet_filtered_v2) ==> will get error
head(top)
ctl3 <- rownames(MVal_filtered) %in% rownames(top[top$adj.P.Val > 0.5,])
table(ctl3)
grp <- factor(phenoData$FST.subgroup, labels=c("U-DFSP","FS-DFSP"))
rfit5 <- RUVfit(Y = MVal_filtered, X = grp, ctl = ctl3) # Stage 2 analysis
rfit6 <- RUVadj(Y = MVal_filtered, fit = rfit5)


AnnoEPICv2 <- data.frame(annEPICv2)
AnnoEPICv2Sub <- AnnoEPICv2[match(rownames(beta_subset), annEPICv2$Name),
                          c(1:4,12:19,28:ncol(AnnoEPICv2))]

# Display the top n DMPs
DMPs <- topTable(lfit2, number = nrow(beta_subset), coef=ncol(design), sort.by="p", genelist = AnnoEPICv2Sub)
DMPs_subset <- topTable(lfit2, number = nrow(beta), coef=ncol(design), sort.by="p", genelist = AnnoEPICv2Sub,
                 p.value = 0.05)

DMG <- c("RFFL","DMTN","GNA12","AQP1","MIR145","DDR2","PLEC","APPL2","TNS1","PDLIM1","NHERF1","ZEB2")

DMPs_subset_genelist <- strsplit2(DMPs_subset$UCSC_RefGene_Name,split=";")[,1]
keep <- DMPs_subset_genelist %in% DMG
DMPs_subset_DMG_genelist <- DMPs_subset_genelist[keep]
DMPs_subset_DMG <- DMPs_subset[keep,]
names(DMPs_subset_DMG_genelist) <- rownames(DMPs_subset[keep,])
DMPs_subset_probe <- rownames(DMPs_subset[keep,])

DMPs_genelist <- strsplit2(DMPs$UCSC_RefGene_Name,split=";")[,1]
keep <- DMPs_genelist %in% DMG
DMPs_DMG_genelist <- DMPs_genelist[keep]
DMPs_DMG <- DMPs[keep,]
DMPs_DMG_df <- data.frame(probes = names(DMPs_DMG_genelist), gene = DMPs_DMG_genelist)
names(DMPs_DMG_genelist) <- rownames(DMPs[keep,])
genelist_final <- DMPs_DMG_genelist[!duplicated(DMPs_DMG_genelist)]

DMPs_subset_probe <- rownames(DMPs_subset[keep,])





head(DMPs)
write.csv(DMPs_subset, "D:/Sarcoma/Methylation_profiling/DFSP/DMPs_PDGFB_vs_PDGFD_excludeMetastasis.csv", 
          row.names = TRUE)

par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(beta, cpg=cpg, pheno=phenoData$Molecular.subtype, ylab = "Beta values")
})

library(EnhancedVolcano)

keyvals <- ifelse(
  DMPs$logFC < 0 & DMPs$adj.P.Val <0.05, 'blue',
  ifelse(DMPs$logFC > 0 & DMPs$adj.P.Val <0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Hypermethylated'
names(keyvals)[keyvals == 'black'] <- 'N.S.'
names(keyvals)[keyvals == 'blue'] <- 'Hypomethylated'

EnhancedVolcano(DMPs,
                lab = NA,
            #    selectLab = rownames(DMPs)[which(names(keyvals) %in% c('Hypermethylated', 'Hypomethylated'))],
             #   colCustom = keyvals,
            #    colAlpha = 1,
                x = 'logFC',
                y = 'adj.P.Val',
                legendPosition = 'right',
                drawConnectors = FALSE,
                labSize = 3,
                pCutoff = 0.05,
                FCcutoff = 0.3,
                pointSize = 1.0,
                legendIconSize = 2,
                legendLabSize = 12,
                max.overlaps = Inf,
                xlim = c(-1, 1),
                ylim = c(0,35),
            #    legendLabels = NA,
                title = 'PDGFB DFSP versus PDGFD DFSP')






##**RUV analysis**##
# For EPIC_V2, please run these lines to remove replicate and non-cg probes:
# remove non-cg probes
MVal <- MVal[grepl("^cg", rownames(MVal)),]
# select replicate probes with best sensitivity
MVal <- DMRcate::rmPosReps(MVal, filter.strategy="sensitivity")

# setup the factor of interest
grp <- factor(phenoData$FST.subgroup, labels=c("U-DFSP","FS-DFSP"), levels = c("U-DFSP","FS-DFSP"))
# extract Illumina negative control data
INCs <- getINCs(RGSet_v2)
head(INCs)
INCs_all <- INCs
INCs <- INCs_all[,selection]

# add negative control data to M-values
MVal_INC <- rbind(MVal,INCs)
# create vector marking negative controls in data matrix
ctl1 <- rownames(MVal_INC) %in% rownames(INCs)
table(ctl1)

rfit1 <- RUVfit(Y = MVal_INC, X = grp, ctl = ctl1) # Stage 1 analysis
rfit2 <- RUVadj(Y = MVal_INC, fit = rfit1)

top1 <- topRUV(rfit2, num=Inf, p.BH = 1)
head(top1)

ctl2 <- rownames(M) %in% rownames(top1[top1$p.BH_X1.PDGFD > 0.5,])
table(ctl2)

# Perform RUV adjustment and fit
rfit3 <- RUVfit(Y = M, X = grp, ctl = ctl2) # Stage 2 analysis
rfit4 <- RUVadj(Y = M, fit = rfit3)

# Look at table of top results
top <- topRUV(rfit4, number = Inf, sort.by = "p", p.BH = 1)
top_cpg <- strsplit2(rownames(top), split = "_")[,1]
table(top$p.BH_X1.PDGFD < 0.05)

beta <- getBeta(MSet_filtered_v2)
# make sure that order of beta values matches orer after analysis
beta <- beta[match(rownames(top),rownames(beta)),]
beta_PDGFB <- rowMeans(beta[,grp=="PDGFB"])
beta_PDGFD <- rowMeans(beta[,grp=="PDGFD"])
Delta_beta <- beta_PDGFD - beta_PDGFB
sigDM <- top$p.BH_X1.PDGFD < 0.05 & abs(Delta_beta) > 0.3
table(sigDM)
sigCpGs <- strsplit2(names(sigDM[sigDM]), split = "_")[,1]

##Else if SigDM is very large, select the top k CpG sites
topCpGs <- topRUV(rfit4, number=10000)
sigCpGs <- strsplit2(rownames(topCpGs), split = "_")[,1]

###***Results generated from limma***###
top <- topTable(lfit2, number = Inf)
table(top$adj.P.Val < 0.05)
# setup the factor of interest
grp <- factor(phenoData$Molecular.subtype, labels=c("PDGFB","PDGFD"), levels = c("PDGFB","PDGFD"))

beta <- beta[match(rownames(top),rownames(beta)),]
beta_PDGFB <- rowMeans(beta[,grp=="PDGFB"])
beta_PDGFD <- rowMeans(beta[,grp=="PDGFD"])
Delta_beta <- beta_PDGFD - beta_PDGFB
sigDM <- top$adj.P.Val < 0.05 & abs(Delta_beta) > 0.3
table(sigDM)
sigCpGs <- strsplit2(names(sigDM[sigDM]), split = "_")[,1]

top_cpg <- strsplit2(rownames(top), split = "_")[,1]
#topCpGs <- topTable(lfit2, number=10000)

# Check number of genes that significant CpGs are annotated to
check <- getMappedEntrezIDs(sig.cpg = sigCpGs)
length(check$sig.eg)

gst <- gometh(sig.cpg=sigCpGs, all.cpg=top_cpg, collection="GO",
              plot.bias=TRUE, sig.genes = TRUE)
topGSA(gst, n=10)
gst.order <- topGSA(gst, n=Inf)

write.csv(gst.order, "D:/Sarcoma/Methylation_profiling/DFSP/DMP_GSEA_GO_PDGFBvsPDGFD.csv")

gst.go.prom <- gometh(sig.cpg=sigCpGs, all.cpg=top_cpg, 
                        collection="GO", sig.genes = TRUE,
                      genomic.features = c("TSS200", "TSS1500","1stExon"))
topGSA(gst.go.prom, n=10)
gst.kegg <- gometh(sig.cpg=sigCpGs, all.cpg=top_cpg, collection="KEGG",
                        sig.genes = TRUE)
topGSA(gst.kegg, n = 10)

gst.kegg.prom <- gometh(sig.cpg=sigCpGs, all.cpg=top_cpg, collection="KEGG",
                   sig.genes = TRUE,genomic.features = c("TSS200", "TSS1500","1stExon"))

topGSA(gst.kegg.prom, n=10)

# Add a new column to the dataframe for -log of "P.DE"
gst.order$log <- -log(gst.order$P.DE)

# Create the plot
ggplot(gst.order[1:6,], aes(x=log, y=reorder(TERM, log))) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=TERM), hjust=1, color="white", nudge_x = -0.5, size = 5) +
  labs(x="-log P", y="") +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())


##**DMR analysis**##
# DMR analysis
library("DMRcate")

plot(density(beta[,1]), col="forestgreen", xlab="Beta value", ylim=c(0, 6),
     main="DFSP", lwd=2)
invisible(sapply(2:5, function (x) lines(density(beta[,x]), col="forestgreen", lwd=2)))
invisible(sapply(6:10, function (x) lines(density(beta[,x]), col="orange", lwd=2)))
legend("topleft", c("PDGFB", "PDGFD"), text.col=c("forestgreen", "orange"))

MVal <- minfi::logit2(beta)

MVal.noSNPs <- rmSNPandCH(MVal, rmcrosshyb = FALSE)
#MVal.noSNPs.repmean <- rmPosReps(MVal.noSNPs, filter.strategy="mean")
design <- model.matrix(~FST.subgroup, data = phenoData)
myannotation <- cpg.annotate("array", object=MVal.noSNPs, what = "M",
                             arraytype = "EPICv2", epicv2Filter = "sensitivity",
                             epicv2Remap = TRUE, analysis.type="differential", 
                             design=design, coef=2,fdr=0.05)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, betacutoff = 0.2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg38")

results.df <- data.frame(results.ranges)

write.csv(results.df, "D:/Sarcoma/Methylation_profiling/DFSP/DMR_FS-DFSPvsU-DFSP.csv")

# Plot one of the top DMRs
groups <- c("FS-DFSP"="forestgreen", "U-DFSP"="orange") 
cols <- groups[as.character(phenoData$FST.subgroup)] 
DMR_plot <- DMR.plot(ranges=results.ranges, dmr=3, CpGs=myannotation, what = "Beta", 
         arraytype = "EPICv2", phen.col=cols, genome = "hg38")
pdf(DMR_plot, "D:Sarcoma/Methylation_profiling/DFSP/DMR.pdf")

#Results were ranked by Fisher's multiple comparison statistic and filtered for
#those DMRs with both FDR and Stouffer scores less than 0.001. DMRs were then annotated to the nearest
#gene transcriptional start sites, based on ENSEMBL genome
#annotations. Gene Ontology (GO) analysis of differentially
#methylated gene regions was performed using
#the gometh function in the missMethyl package

cpg <- strsplit2(rownames(MVal.noSNPs), split = "_")[,1]
cpg <- rownames(MVal.noSNPs)
gst.region <- goregion(results.ranges, all.cpg=cpg,
                       collection = "GO", array.type="EPIC_V2", plot.bias = TRUE)
gst.region.order <- topGSA(gst.region, n=Inf)

write.csv(gst.region.order, "D:/Sarcoma/Methylation_profiling/DFSP/GSEA_GO_PDGFBvsPDGFD.csv")

##**Perform GSEA on DMR##**##
gst.region.prom <- goregion(results.ranges, all.cpg=rownames(MVal.noSNPs),
                       genomic.features = c("TSS200", "TSS1500", "1stExon"),  #Limit to promoter region
                       collection = "GO", array.type="EPIC_V2", plot.bias = TRUE,
                       sig.genes = TRUE)
gst.region.prom.order <- topGSA(gst.region.prom, n=Inf)
gst.region.prom.order2 <- topGSA(gst.region.prom, n=Inf)
gst.region.prom.order3 <- topGSA(gst.region.prom, n=Inf)

write.csv(gst.region.prom.order2[1:20,], "D:/Sarcoma/Methylation_profiling/DFSP/GSEA_GO_Prom_Pre-FSTvsPost-FST.csv")


gst.region.kegg <- goregion(results.ranges, all.cpg=rownames(MVal.noSNPs),
                       collection = "KEGG", array.type="EPIC_V2", plot.bias = TRUE)
topGSA(gst.region.kegg, n=10)



##metylglm
DMP <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/myDMP_ChAMP.csv")
DMP_subset <- DMP[DMP$deltaBeta >0.3,]
CPG.PVAL <- DMP_subset$P.Value
names(CPG.PVAL) <- DMP_subset$Name
cpg.gene <- data.frame(CpG = rownames(annEPICv2.df), 
                       Gene = strsplit2(annEPICv2.df$UCSC_RefGene_Name, split = ";")[,1])
FullAnnot = prepareAnnot(cpg.gene) 
test <- methylglm(cpg.pval = CPG.PVAL, FullAnnot = FullAnnot, minsize = 200, maxsize = 500,
                  GS.type = "GO")



##**Conumee for CNV analysis**##
library(conumee2)
data(exclude_regions)
data(detail_regions)  # example detail regions for hg19 (not compatible with hg38)
data(detail_regions.hg38) # detail regions for hg38 (only for EPICv2)

anno_overlap <- CNV.create_anno(bin_minprobes = 15, bin_minsize = 50000,
                                bin_maxsize = 5000000, array_type = "overlap", chrXY = FALSE,
                                exclude_regions = exclude_regions, detail_regions = detail_regions)
anno_450k <- CNV.create_anno(bin_minprobes = 15, bin_minsize = 50000,
                             bin_maxsize = 5000000, array_type = "450k", chrXY = FALSE,
                             exclude_regions = exclude_regions, detail_regions = detail_regions)

anno_EPIC <- CNV.create_anno(bin_minprobes = 15, bin_minsize = 50000,
                             bin_maxsize = 5000000, array_type = "EPIC", chrXY = FALSE,
                             exclude_regions = exclude_regions, detail_regions = detail_regions)

anno_EPICv2 <- CNV.create_anno(bin_minprobes = 15, bin_minsize = 50000,
                               bin_maxsize = 5000000, array_type = "EPICv2", chrXY = FALSE,
                               exclude_regions = exclude_regions, detail_regions = detail_regions)

#For EPIC probes, run this command
anno_EPIC@probes <- anno_EPIC@probes[names(anno_EPIC@probes) %in% names(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19))]

anno3@probes <- anno3@probes[names(anno3@probes) %in% names(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19))]

query_name <- c("DFSP-067-T", "DFSP-037-T")
control_names <- c("DFSP-166-N", "DFSP-332-N", "DFSP-336-N", "DFSP-337-N", 
                   "DFSP-341-N", "DFSP-348-N", "DFSP-354-N", "DFSP-356-N")
Metastasis <- c("DFSP-037-T", "DFSP-141-T-M1", "DFSP-141-T-M2", 
                "DFSP-272-T", "DFSP-289-T-M1", "DFSP-319-T-M1")
# Convert RGSet to MSet
MSet_CNV <- preprocessIllumina(RGSet_v2_filtered)
#rownames(MSet) <- strsplit2(rownames(MSet), split = "_")[,1]

# Combine intensity values
minfi.data <- CNV.load(MSet_CNV, names = phenoData$Sample.ID)

minfi.controls <- phenoData$Sample.ID %in% control_names

minfi.tumour <- phenoData$Sample.ID %in% query_name
minfi.tumour <- !(phenoData$Sample.ID %in% control_names)

minfi.metastasis <- phenoData$Sample.ID %in% Metastasis

# Perform CNV analysis
x <- CNV.fit(minfi.data[minfi.tumour], minfi.data[minfi.controls], anno_EPICv2)
x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)

x <- CNV.focal(x, sig_cgenes = TRUE)

CNV.summaryplot(x[phenoData_all_tumour$Sample.ID[phenoData_all_tumour$FST.subgroup== "DFSP"]])
CNV.summaryplot(x[phenoData_all_tumour$Sample.ID[phenoData_all_tumour$FST.subgroup== "Pre-FST"]])
CNV.summaryplot(x[phenoData_all_tumour$Sample.ID[phenoData_all_tumour$FST.subgroup== "Post-FST"]])
CNV.summaryplot(x[phenoData_all_tumour$Sample.ID[phenoData_all_tumour$FST.subgroup== "FS-DFSP"]])

segments <- CNV.write(x, what = "segments") 
detail <- CNV.write(x, what = "detail") 
gistic <- CNV.write(x, what = "gistic")
focal <- CNV.write(x, what = "focal")
#what can choose from "probes", "bins", "detail" , "segments", "gistic" (for downstream processing) or "focal" (results from CNV.focal) 


x.df <- data.frame(x)
CNV.genomeplot(x)
CNV.genomeplot(x, chr = 'chr9')
CNV.detailplot(x["DFSP-362-T-P2"], name = 'CCNE1')
CNV.detailplot_wrap(x)

#Extract the segmentation file
seg.summary <- x@seg$summary
seg.summary.df <- do.call(rbind, seg.summary)
colnames(seg.summary.df)[1] <- "Sample.ID"
seg.summary.df.pheno <- merge(seg.summary.df, phenoData, by = "Sample.ID", all.x = TRUE)
write.table(seg.summary.df.pheno, "D:/Sarcoma/Methylation_profiling/DFSP/DFSP_CNVseg_alltumour.txt", sep = "\t",row.names = FALSE )

#Separate seg files according to molecular subtype
seg.summary.df.PDGFD <- seg.summary.df.pheno[seg.summary.df.pheno$Molecular.subtype=="PDGFD",]
seg.summary.df.PDGFB <- seg.summary.df.pheno[seg.summary.df.pheno$Molecular.subtype=="PDGFB" &
                                               seg.summary.df.pheno$Main=="Yes" &
                                               !(seg.summary.df.pheno$Specimen.Nature %in% c("Metastasis", "Recurrence (Post-imatinib)")),]
seg.summary.df.PDGFB.classic <- seg.summary.df.pheno[seg.summary.df.pheno$Molecular.subtype=="PDGFB" &
                                               seg.summary.df.pheno$Main=="Yes" &
                                               !(seg.summary.df.pheno$Specimen.Nature %in% c("Metastasis", "Recurrence (Post-imatinib)")) &
                                                 seg.summary.df.pheno$Histology.Subtype == "Classic",]
write.table(seg.summary.df.PDGFD, "D:/Sarcoma/Methylation_profiling/DFSP/DFSP_CNVseg_PDGFD.txt", sep = "\t",row.names = FALSE)
write.table(seg.summary.df.PDGFB, "D:/Sarcoma/Methylation_profiling/DFSP/DFSP_CNVseg_PDGFB.txt", sep = "\t",row.names = FALSE)
write.table(seg.summary.df.PDGFB.classic, "D:/Sarcoma/Methylation_profiling/DFSP/DFSP_CNVseg_PDGFB_classic.txt", sep = "\t",row.names = FALSE)


# Extract copy number of the detailed region
detailed.data <- x@detail$ratio
detailed.df <- as.data.frame(do.call(rbind, detailed.data))
write.table(detailed.df, "D:/Sarcoma/Methylation_profiling/DFSP/DFSP_CNV_detailedRegions.txt", sep = "\t",row.names = FALSE)


#Read dataframe containing CNV segmentation info with clinical data
seg.summary.df.pheno <- read.delim("D:/Sarcoma/Methylation_profiling/DFSP/DFSP_CNVseg_alltumour.txt", 
sep = "\t")
Sample.list <- unique (seg.summary.df.pheno$Sample.ID)

#Calculate the fraction of genome altered by CNV in each case
#CNV is defined as seg.mean >0.3, divided by total number of segments
seg.summary.df <- read.delim("D:/Sarcoma/Methylation_profiling/DFSP/DFSP_CNVseg_alltumour.txt", 
                                   sep = "\t")

result <- seg.summary.df |>
  group_by(Sample.ID) |>
  summarize(FGA.CNV = round((sum(num.mark[abs(seg.mean) > 0.3], na.rm = TRUE)/ sum(num.mark, na.rm = TRUE)) * 100,2)
            )



phenoData_subset <- phenoData[phenoData$Histology.Nature!="Skin and subcutis",]
phenoData_subset <- phenoData[phenoData$Main=="Yes" & phenoData$Molecularly.confirmed =="Yes" & phenoData$Histology.Subtype =="Classic" &
                                !(phenoData$Specimen.Nature %in% c("Metastasis", "Recurrence (Post-imatinib)")),]

phenoData_subset$`CDKN2A/B` <- detailed.df$`CDKN2A/B`
phenoData_subset$CDKN2A.status[phenoData_subset$`CDKN2A`< -0.5] <- "Deletion"
phenoData_subset$CDKN2A.status[phenoData_subset$`CDKN2A/B`>= -0.5] <- "No deletion"


library(ggpubr)


##Plots for multiple comparisons:
comparisons = list(c("Pre-FST", "Post-FST"), 
                   c("Post-FST", "FS-DFSP"),
                   c("Pre-FST", "FS-DFSP"),
                   c("U-DFSP", "FS-DFSP"))

p <- ggboxplot(phenoData_new, x = "FST.subgroup", y = "FGA.CNV", fill = "FST.subgroup",
                     palette = "npg", add = "jitter") + 
  stat_compare_means(vjust = 2.8, label.x = 1.2) +
  stat_compare_means(comparisons = comparisons,
                   label = "p.signif",
                   symnum.args = list(
                     cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                     symbols = c("****", "***", "**", "*", "ns")
                   ))

legend_data <- data.frame(
  Symbol = c("****", "***", "**", "*", "ns"),
  Meaning = c("p < 0.0001", "p < 0.001", "p < 0.01", "p < 0.05", "p ≥ 0.05")
)

p + geom_point(data = legend_data, aes(x = Inf, y = Inf, shape = Symbol), alpha = 0) +
  scale_shape_manual(name = "", 
                     values = c("****" = 16, "***" = 16, "**" = 16, "*" = 16, "ns" = 16),
                     labels = c("**** p < 0.0001", "***  p < 0.001", "**   p < 0.01", "*    p < 0.05", "ns   p ≥ 0.05")) +
  theme(legend.position = "right",
        legend.text = element_text(margin = margin(r = 20)))

#Plots for 2 groups only
ggboxplot(phenoData_new, x = "FST.status", y = "FGA.CNV", fill = "FST.status",
               palette = "npg", add = "jitter") + 
  stat_compare_means(vjust = 2.8, label.x = 1.2)


##dpylr way to create new columns
phenoData_new <- phenoData_new |>
  mutate(FST.status = case_when(
    FST.subgroup %in% c("U-DFSP", "Pre-FST") ~ "Non-FST",
    FST.subgroup %in% c("Post-FST", "FS-DFSP") ~ "FST"
  ))
  

  
##**Immune cell deconvolution**##
library(EpiDISH)
beta_mod <- beta
rownames(beta_mod) <- strsplit2(rownames(beta_mod), split = "_")[,1]
out.l <- epidish(beta.m = beta_mod, ref.m = centEpiFibIC.m, method = "RPC") 
out.l2 <- epidish(beta.m = beta_mod, ref.m = centDHSbloodDMC.m, method = "RPC") 

out.l$estF
dim(out.l$ref)
dim(out.l$dataREF)
boxplot(out.l$estF)

frac.m <- hepidish(beta.m = beta, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
frac.m

pheno.meth <- phenoData$Meth.subgroup
celldmc.o <- CellDMC(beta, pheno.meth, out.l$estF)
which(colSums(model.matrix(~ frac.m + pheno.meth:frac.m)[, -1]) == 0)

library(ggplot2)
library(viridis)
library(rstatix)
library(ggpubr)
table(rownames(phenoData) == rownames(out.l$estF))
phenoData.immune <- cbind(phenoData, out.l$estF, out.l2$estF)
phenoData.immune2 <- cbind(phenoData, frac.m)
phenoData.subset <- phenoData.immune[phenoData.immune$Meth.subgroup == "Meth4",]

comparisons <- list(c("Meth1", "Meth2"), c("Meth1", "Meth3"), c("Meth1", "Meth4"), 
                    c("Meth2", "Meth3"), c("Meth2", "Meth4"), c("Meth3", "Meth4"))
ggplot(phenoData, 
        aes (x = Meth.subgroup.probesel, y = IC, fill = Meth.subgroup.probesel)) + 
 # geom_violin (trim = FALSE, width = 1) + 
  geom_boxplot (width = 0.5)+ 
  scale_fill_viridis(discrete = TRUE)  + 
  labs(x="Methylation subtype", y = "IC") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  #stat_compare_means(comparisons = comparisons, vjust = 0, label.x = 1.5, size = 6) +
  stat_compare_means(label.y = 0.5, size = 6) +
  ggtitle("Immune cell percentage")


##**Clustering analysis on CNV**##
#Load CN values from GISTIC output (all_lesion.config_95.txt)
CN_raw <- read.delim("Y:/Methylation_profiling/DFSP/CNV/gistic2/all_lesions.conf_99.txt")
#Preprocess the datasheet
CN <- CN_raw[CN_raw$Amplitude.Threshold == "Actual Copy Change Given",]
del <- c("Descriptor", "Wide.Peak.Limits", "Peak.Limits", "Region.Limits",
         "Broad.or.Focal", "Amplitude.Threshold", "Unique.Name", "X",
         "q.values", "Residual.q.values.after.removing.segments.shared.with.higher.peaks")
CN$Unique.Name <- strsplit2(CN$Unique.Name, split = " - ")[,1]
CN$Unique.Name <- paste0(substr(CN$Unique.Name, 1, 3), " ",
                         strsplit2(CN$Unique.Name, split = " ")[,2], " ", 
                         strsplit2(CN$Unique.Name, split = "Peak ")[,2], " ",
                         CN$Descriptor)
rownames(CN) <- CN$Unique.Name
CN[,del] <- NULL
colnames(CN) <- gsub("\\.", "-", colnames(CN))

#Convert the sample sheet to matrix
sampleMatrix <- as.matrix(CN)

#Load samplesheet
phenoData <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_allTumour_CNV_Meth.csv", row.names = 1)
rownames(phenoData) <- phenoData$Sample.ID

table(rownames(phenoData) == colnames(sampleMatrix))

phenoData<- phenoData[phenoData$FST.subgroup %in% c("Pre-FST", "Post-FST"),]
sampleMatrix <- sampleMatrix[, colnames(sampleMatrix) %in% phenoData$Sample.ID]

phenoData <- phenoData[match(colnames(sampleMatrix),rownames(phenoData)),]


anno <- c("Sex", "Age.cat", "Site.Category.Broad", "Local.recurrence", 
          "Metastasis", "Molecular.subtype", "Histology.subtype", "Ploidy.cat",
          "MSI.cat", "FST.subgroup"
)



annotation_colors <- list(
  CN.subgroup = c("CN1" = "tomato", "CN2" = "seagreen", "CN3" = "dodgerblue"),
  FST = c("No" = "firebrick", "Yes" = "orchid"),
  Metastasis = c("No" = "steelblue", "Yes" = "gold"),
  Local.recurrence = c("No" = "mediumseagreen", "Yes" = "slateblue"),
  Molecular.subtype = c("PDGFB" = "firebrick","PDGFD" = "orchid"),
  Site.Category.Broad = c(
    "Trunk" = "magenta", "Extremity" = "coral", "H&N" = "navy"),
  Histology.subtype = c(
    "Classic" = "#377EB8",
    "FS" = "#E41A1C",
    "Myxoid" = "#4DAF4A",
    "Pigmented" = "#984EA3")
  )


library(RColorBrewer)
color_palette1 <- brewer.pal(4, "Set1")
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-4, 4, length.out = 100) 

heatmap <- pheatmap(sampleMatrix, cluster_rows = T, show_rownames=T,
                    labels_col = phenoData$Sample.ID, show_colnames = T,
                    annotation = phenoData[,anno], 
                    annotation_colors = annotation_colors, #breaks = breaks,
                    cutree_cols = 4, color = my_palette,
                    clustering_method = "ward.D", treeheight_row = 0,
                    border_color=NA, fontsize = 10, scale="row",
                    fontsize_row = 10, fontsize_col = 5, height=30)

#one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#Default "complete"
# The resulting object contains the clustering dendrograms
column_dendrogram <- heatmap$tree_col

# You can cut the dendrogram to create clusters/subgroups
col_clusters <- cutree(column_dendrogram, k = 4)  # adjust 'k' to the number of clusters/subgroups you want

# You can then add this information back to your original data
phenoData$CN.subgroup4 <- factor(col_clusters, labels = c("CN2", "CN3", "CN1", "CN4"))  # adjust labels as needed
phenoData$CN.subgroup4 <- relevel(phenoData$CN.subgroup4, ref = "CN1")

write.csv(phenoData, "D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_97pt_Meth_immune_CNV.csv", row.names = TRUE)


##Plot Sample level CNV
CN <- read.delim("D:/Sarcoma/Methylation_profiling/DFSP/DFSP_CNVseg.txt")
CN_PDGFD <- CN[CN$Molecular.subtype == "PDGFD",]
CN_PDGFD <- CN_PDGFD[,c("Sample.ID", "chrom","loc.start","loc.end","seg.mean")]
colnames(CN_PDGFD) <- c("sample","chromosome","start","end","segmean")
CN_PDGFB <- CN[CN$Main=="Yes" & CN$Molecular.subtype =="PDGFB",]
CN_PDGFB <- CN_PDGFB[,c("Sample.ID", "chrom","loc.start","loc.end","seg.mean")]
colnames(CN_PDGFB) <- c("sample","chromosome","start","end","segmean")

plot <- cnFreq(CN_PDGFD, genome = "hg38", CN_low_cutoff = -0.3, CN_high_cutoff = 0.3, 
       plotLayer = ggplot2::facet_grid(.~chromosome, scales="free"))

plot <- cnFreq(CN_PDGFB, genome = "hg38", CN_low_cutoff = -0.3, CN_high_cutoff = 0.3, 
       plotLayer = ggplot2::facet_grid(.~chromosome, scales="free"))

plot$data$chromosome <- factor(plot$data$chromosome, 
                               levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))

print(plot)

##**iCluster analysis**##
library(iClusterPlus)

Meth_iCluster <- t(beta_sort)
Meth_iCluster <- Meth_iCluster[match(phenoData$IDAT, rownames(Meth_iCluster)), ]
CN_iCluster <- t(CN)

fit.single=iClusterPlus(dt1=Meth_iCluster,dt2=CN_iCluster,
                        type=c("gaussian","gaussian"),
                        lambda=c(0.61,0.90),K=2,maxiter=10)


for(k in 1:5){
  cv.fit = tune.iClusterPlus(dt1=Meth_iCluster, dt2=CN_iCluster, cpus = 1,
                             type=c("gaussian","gaussian"),K=k,n.lambda=NULL,
                             scale.lambda=c(1,1),maxiter=20)
  save(cv.fit, file=paste("D:/Sarcoma/Methylation_profiling/DFSP/iCluster/cv.fit.k",k,".Rdata",sep=""))
  }

output=alist()
files=grep("D:/Sarcoma/Methylation_profiling/DFSP/iCluster/cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
  }
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)
minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

clusters=getClusters(output)
rownames(clusters)=rownames(gbm.exp)
colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
#write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
k=2
best.cluster=clusters[,k]
best.fit=output[[k]]$fit[[which.min(BIC[,k])]]
plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     + ylab="%Explained Variation")