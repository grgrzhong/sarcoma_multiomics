# Example: Testing clonality enrichment for driver vs passenger mutations

# Create contingency table
# Rows: Mutation type (Driver, Passenger)
# Columns: Clonality status (Clonal, Subclonal)

clonality_data <- matrix(
    c(
        # Clonal  Subclonal
        45, 15, # Driver mutations
        120, 200 # Passenger mutations
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        c("Driver", "Passenger"),
        c("Clonal", "Subclonal")
    )
)

print(clonality_data)

# Perform Fisher's exact test
fisher_result <- fisher.test(clonality_data)
print(fisher_result)

# Extract key results
cat("P-value:", fisher_result$p.value, "\n")
cat("Odds Ratio:", fisher_result$estimate, "\n")
cat("95% CI:", fisher_result$conf.int[1], "-", fisher_result$conf.int[2], "\n")

# Chi-square test for the same data
chisq_result <- chisq.test(clonality_data)
print(chisq_result)

# Check expected frequencies (should be >5 for valid chi-square)
print("Expected frequencies:")
print(chisq_result$expected)

# Extract results
cat("Chi-square statistic:", chisq_result$statistic, "\n")
cat("P-value:", chisq_result$p.value, "\n")
cat("Degrees of freedom:", chisq_result$parameter, "\n")

# Example: Testing multiple candidate driver genes
library(dplyr)

# Create sample data for multiple genes
set.seed(123)
gene_data <- data.frame(
  gene = rep(c("TP53", "CTNNB1", "TERT", "ARID1A", "RB1"), each = 100),
  mutation_type = "driver",
  clonality = sample(c("clonal", "subclonal"), 500, replace = TRUE, prob = c(0.7, 0.3))
)

# Add passenger mutations as control
passenger_data <- data.frame(
  gene = "passengers",
  mutation_type = "passenger", 
  clonality = sample(c("clonal", "subclonal"), 1000, replace = TRUE, prob = c(0.4, 0.6))
)

# Combine data
all_data <- rbind(gene_data, passenger_data)

# Test each driver gene vs passengers
driver_genes <- c("TP53", "CTNNB1", "TERT", "ARID1A", "RB1")
results <- data.frame()

for (gene in driver_genes) {
    # Create 2x2 table for current gene vs passengers
    gene_subset <- subset(all_data, gene == gene | gene == "passengers")

    contingency <- table(gene_subset$mutation_type, gene_subset$clonality)

    # Fisher's exact test
    fisher_test <- fisher.test(contingency)

    # Store results
    results <- rbind(results, data.frame(
        gene = gene,
        p_value = fisher_test$p.value,
        odds_ratio = fisher_test$estimate,
        ci_lower = fisher_test$conf.int[1],
        ci_upper = fisher_test$conf.int[2]
    ))
}

# Apply multiple testing correction
results$p_adjusted <- p.adjust(results$p_value, method = "BH")

# Display results
print(results)

# Identify significantly enriched genes (p < 0.05)
significant_genes <- subset(results, p_adjusted < 0.05)
cat("Genes with significant clonality enrichment:\n")
print(significant_genes$gene)

# Visualize clonality patterns
library(ggplot2)

# Create summary data for plotting
summary_data <- all_data |>
  group_by(gene, clonality) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(gene) |>
  mutate(proportion = count / sum(count))

# Plot clonal vs subclonal proportions
ggplot(summary_data, aes(x = gene, y = proportion, fill = clonality)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Clonality Distribution by Gene",
       x = "Gene", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
