
gene_set = strsplit(readLines(here("data/Test/c5.go.bp.v2025.1.Hs.symbols.gmt")), "\t")
names(gene_set) = sapply(gene_set, "[", 1)
gene_set = lapply(gene_set, "[", -(1:2))
