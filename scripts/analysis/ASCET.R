source("G:/Sarcoma/ascets-1.1.3/ascets_resources.R")
cna <- read.delim("G:/Sarcoma/GLI1/CNV-GATK/GLI-004-T.called.igv.seg")
cna$Call<- NULL
colnames(cna) <- c("sample", "chrom", "segment_start", "segment_end", "num_mark", "log2ratio")
cna$chrom <- gsub("chr", "", cna$chrom)

cytoband <- read.delim ("G:/Sarcoma/ascets-1.1.3/genomic_arm_coordinates_hg38.txt")
noise <- read.delim("G:/Sarcoma/ascets-1.1.3/sample_data/input/lcr_example.txt")

ascets_output <- ascets(cna, 
                        cytoband, 
                        min_boc = 0.5, 
                        name = "GLI-004", 
                        keep_noisy = FALSE, 
                        threshold = 0.2, alteration_threshold = 0.7)

write_outputs_to_file(ascets_output, location = "G:/Sarcoma/ASCETS/")
