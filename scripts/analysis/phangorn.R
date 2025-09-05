# Install packages if needed
# install.packages(c("phangorn", "ape"))

library(phangorn)
library(ape)

# Sample presence/absence matrix (replace with your data, e.g., read.csv("your_matrix.csv", row.names=1))
# Rows: samples (taxa), Columns: mutations (characters), Values: 0/1
data_matrix <- matrix(c(
  0,0,0,0,  # Normal
  1,1,0,0,  # PT1
  1,1,1,0,  # PT2
  1,0,1,1,  # LN1
  1,0,1,1   # DM1
), nrow=5, byrow=TRUE)
rownames(data_matrix) <- c("Normal", "PT1", "PT2", "LN1", "DM1")
colnames(data_matrix) <- c("TP53", "EGFR", "PTEN", "Tp-gain")

# Convert to phyDat (binary user-defined levels)
phy_data <- as.phyDat(data_matrix, type="USER", levels=c(0,1))

# Reconstruct tree with parsimony ratchet
tree_ratchet <- pratchet(phy_data, trace=0, minit=100)
tree_ratchet <- acctran(tree_ratchet, phy_data)  # Assign branch lengths
tree_ratchet <- di2multi(tree_ratchet)  # Handle multifurcations
tree_ratchet <- unique(tree_ratchet)  # Ensure unique trees if multiple

# Bootstrapping (1000 replicates)
bs_fun <- function(x) pratchet(x, trace=0, minit=100)
bs_trees <- bootstrap.phyDat(phy_data, bs_fun, bs=1000)

# Root the tree at the midpoint
rooted_tree <- midpoint(tree_ratchet)

# Plot with bootstrap support (basic plot)
plotBS(midpoint(tree_ratchet), bs_trees, type="phylogram", main="Parsimony Tree with Bootstrap")

# Parsimony score
parsimony(tree_ratchet, phy_data)