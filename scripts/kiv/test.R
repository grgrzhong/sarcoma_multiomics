# For categorical variables:
# - Chi-square test (χ2 test): For larger sample, expected frequenty >5; can be used with two or more categories
# - Fisher's exact test: For smaller sample, expected frequency <5; optimized for 2 x 2 table
   
# For continuous variables:
# If 2 groups only (e.g. metastasis versus no metastasis)
# - T-test - independent samples t-test; for normally distributed data
# - non-parametric Wilcoxon or Mann-Whitney U test: if the data is not normally distributed
# If more than 2 groups:
# - One-way ANOVA: for normally distributed
# - non-parametric Kruskal-Wallis test: if the data is not normally distributed

# To determine if your data is normally distributed
# Shapiro-Wilk test: The null hypothesis is that the data is normally distributed. 
# If the p-value is >0.05, fail to reject the null hypothesis --> normally distributed.
# If p<0.05 --> NOT normally distributed


library(colorspace)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(rstatix)
library(viridis)

##**Load the samplesheet**##

DFSP <- read.csv(
    "/mnt/m/WES/DFSP/clinical/DFSP_WES_clinical.csv",
    header = TRUE
)

DFSP <- DFSP[DFSP$Main == "Yes",]

as_tibble(DFSP)

###***Change the comparison group name as needed***###
outcome <- "Metastasis" 

###***Set categorical factors***###
prognostic_factors_cat <- c("Sex", "Molecular.subtype", "Site.Category.Broad",
                            "Local.recurrence", 
                            "MI.Score",
                            "FNCLCC.Grade", "Contour", "Depth", "FST", 
                            "Histology.subtype",  "Age.cat",
                            "T.Stage", "HRD.cat", "MSI.cat", "Cellularity.cat",
                            "Ploidy.cat"
                            )
#"Metastasis", "CDKN2A.status",

# Function to perform Fisher's Exact Test and return the p-value, contingency table, and proportion table
get_fisher_exact_info <- function(factor_name, data) {
  tbl <- table(data[[factor_name]], data[[outcome]])    
  prop_tbl <- prop.table(tbl, margin = 2)
  fisher_result <- fisher.test(tbl)
  return(list(p_value = fisher_result$p.value, contingency_table = tbl, proportion_table = prop_tbl))
}

# Loop through each factor and output the p-value, contingency table, and proportion table
for (factor_name in prognostic_factors_cat) {
  fisher_exact_info <- get_fisher_exact_info(factor_name, DFSP)
  factor<- data.frame(fisher_exact_info$contingency_table)
  colnames(factor) <- c(factor_name,outcome,"Freq")
  print(factor)
  bar_chart <- ggplot(factor, aes(fill=!!sym(factor_name), y=Freq, x=!!sym(outcome))) + 
    geom_bar(position="fill", stat="identity") +
    theme_classic() +   
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) + 
    annotate("text", x = 1.5, y = 1.05, label = paste("p = ", round(fisher_exact_info$p_value,5)), color = "black", size = 6) + 
    labs(x=outcome, y = "Percent (%)")
  print(bar_chart)
  cat(factor_name, "p-value:", fisher_exact_info$p_value, "\n")
  cat("Contingency table:\n")
  print(fisher_exact_info$contingency_table)
  cat("Proportion table:\n")
  print(fisher_exact_info$proportion_table)
  cat("\n")
}


###***Set continuous factors***###
prognostic_factors_cont <- c("Age", "Size", "Mitotic.count", 
                             "HRDscore", "MSIsensor.Score", 
                             "Sequenza.Cellularity", "Sequenza.Ploidy")

# Dynamically create color mappings based on unique outcome groups
unique_groups <- unique(DFSP[[outcome]])
fill_colors <- setNames(rainbow(length(unique_groups)), unique_groups)  # Use rainbow() or any other palette
point_colors <- setNames(darken(fill_colors, amount = 0.5), unique_groups)  # Darken fill colors for points


get_summary_stats <- function(var_name, data, outcome) {
  # Create the formula for aggregation
  formula_str <- paste(var_name, "~", outcome)
  formula_obj <- as.formula(formula_str)
  
  # Calculate median, mean, and range
  median_df <- aggregate(formula_obj, data, FUN = median)
  mean_df <- aggregate(formula_obj, data, FUN = mean)
  range_df <- aggregate(formula_obj, data, FUN = function(x) {
    range_values <- range(x)
    data.frame(Range.Lower = range_values[1], Range.Upper = range_values[2])
  })
  
  # Rename columns for clarity
  colnames(median_df)[2] <- "Median"
  colnames(mean_df)[2] <- "Mean"
  
  # Flatten the range_df (it might be a nested data frame)
  range_df <- do.call(data.frame, range_df)
  
  # Merge the data frames
  merged_df <- median_df %>%
    left_join(mean_df, by = outcome) %>%
    left_join(range_df, by = outcome)
  
  return(merged_df)
}


# Function to perform Wilcoxon test
get_wilcox_pvalue <- function(var_name, data, outcome) {
  # Create the formula for the Wilcoxon test
  formula_str <- paste(var_name, "~", outcome)
  formula_obj <- as.formula(formula_str)
  
  # Perform the Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(formula_obj, data = data, exact = FALSE)
  
  # Extract the p-value
  p_value <- wilcox_test$p.value
  
  return(p_value)
}

# Loop through each prognostic factor and create a boxplot
for (factor in prognostic_factors_cont) {
  p <- ggplot(DFSP, aes(x = !!sym(outcome), y = .data[[factor]], fill = !!sym(outcome))) +  # Add fill for color
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Add boxplot with transparency
    geom_jitter(width = 0.2, alpha = 0.5, aes(color = !!sym(outcome))) +  # Add jittered points
    stat_compare_means(method = "wilcox.test", label = "p.format", 
                       label.x = 1.2, size = 5) +  # Add p-value
    labs(title = "",
         x = outcome,
         y = factor) +
    theme_minimal() +
    scale_fill_manual(values = fill_colors) +  # Custom fill colors
    scale_color_manual(values = point_colors)  + # Custom point colors
  theme(
    axis.title = element_text(size = 14),  # Enlarge axis title text
    axis.text = element_text(size = 12),   # Enlarge axis tick labels
    legend.title = element_text(size = 14),  # Enlarge legend title text
    legend.text = element_text(size = 12)    # Enlarge legend labels
  )
  print(p)
  
  summary_stats <- get_summary_stats(factor, DFSP, outcome)
  p_value <- get_wilcox_pvalue(factor, DFSP, outcome)
  print(factor)
  print(paste("p-value:",p_value))
  print(summary_stats)

  
}



##**Perform univariate and multivariate analysis using GLM**##
DFSP$Metastasis <- as.factor(DFSP$Metastasis)
logistic_regression <- glm(Metastasis~Sex, data = DFSP, family = binomial)

#Output p value of each variable
summary(logistic_regression)

#Output odds ratio (OR)
exp(coef(logistic_regression))  

#Output CI of OR
exp(confint(logistic_regression))

DFSP_mol_confirmed$Tumour.size..cm.


###***Input the prognostic factors to be analyzed***###
prognostic_factors_cat <- c("Sex", "Site.Category.Broad", 
                             "Contour","Depth", "FST", "Molecular.subtype", 
                            "FNCLCC.Grade", "Age.cat" ,"T.Staging","MSI.cat",
                            "HRD.cat"
                             )
prognostic_factors_cont <- c("Age", "Size", "HRDscore", "MSIsensor.Score")
prognostic_factors <- c(prognostic_factors_cat, prognostic_factors_cont)


###Perform univariate analysis for each prognostic factor
univariate_analysis <- function(factor) {
  formula <- as.formula(paste("Metastasis ~", factor))
  model <- glm(formula, data = DFSP, family = "binomial")
  
  # Calculate the odds ratio
  odds_ratio <- exp(coef(model)[2])
  
  # Calculate the 95% confidence interval
  ci <- exp(confint(model)[2, ])
  
  # Calculate the p-value
  p_value <- summary(model)$coefficients[2, 4]
  
  return(c(odds_ratio, ci, p_value))
}

# Perform univariate analysis for each prognostic factor
results <- t(sapply(prognostic_factors, univariate_analysis))

# Create a data frame with the results
results_df <- data.frame(
  PrognosticFactor = prognostic_factors,
  OddsRatio = results[, 1],
  LowerCI = results[, 2],
  UpperCI = results[, 3],
  PValue = results[, 4]
)

# Print the results data frame
print(results_df)

###Print the count and proportion table for each categorical prognostic factor
# # Function to generate a counting table and proportions table for a given prognostic factor
generate_table <- function(factor) {
  freq_table <- table(DFSP[[factor]], DFSP$Metastasis)
  prop_table <- prop.table(table(DFSP[[factor]], DFSP$Metastasis), margin = 1)*100
  return(list(freq_table = freq_table, prop_table = prop_table))
}

# Generate counting tables for each categorical prognostic factor
tables_list <- lapply(prognostic_factors_cat, generate_table)

# Print the tables
for (i in seq_along(tables_list)) {
  cat("\nTables for", prognostic_factors[i], ":\n")
  cat("\nFrequency table:\n")
  print(tables_list[[i]]$freq_table)
  cat("\nProportions table:\n")
  print(tables_list[[i]]$prop_table)
}


###Print the median, mean and range of each continuous prognostic factor in different groups
# Split the data by adverse event group
grouped_data <- split(DFSP_mol_confirmed, DFSP_mol_confirmed$Metastasis)

# Function to calculate the median, mean, and range of a given prognostic factor
calculate_stats <- function(values) {
  median_value <- median(values, na.rm = TRUE)
  mean_value <- mean(values, na.rm = TRUE)
  range_value <- range(na.omit(values))
  return(list(median = median_value, mean = mean_value, range = range_value))
}

# Calculate the statistics for each prognostic factor in each adverse event group
stats_list <- lapply(grouped_data, function(group) {
  lapply(prognostic_factors_cont, function(factor) {
    calculate_stats(group[[factor]])
  })
})

# Print the statistics
for (group in names(grouped_data)) {
  cat("\nMetastasis:", group, "\n")
  for (i in seq_along(prognostic_factors_cont)) {
    cat("\nStatistics for", prognostic_factors_cont[i], ":\n")
    cat("Median:", stats_list[[group]][[i]]$median, "\n")
    cat("Mean:", stats_list[[group]][[i]]$mean, "\n")
    cat("Range:", paste(stats_list[[group]][[i]]$range, collapse = " to "), "\n")
  }
}


###Multivariate analysis

model <- glm(Metastasis ~ FST + 
               MSI.cat + Age + Depth + Site.Category.Broad, data = DFSP, family = 'binomial')

# Output a summary of the model
summary(model)

# Extract the p-values for each predictor
p_values <- summary(model)$coefficients[, "Pr(>|z|)"]

# Print the p-values
cat("P-values for each predictor:\n")
print(p_values)


#General formula for Cox proportional hazard ratio
library(survival)
cox_ph_model <- coxph(Surv(time_to_event, event_status) ~ x1_continuous + x2_categorical, data = your_data)

DFSP_mol_confirmed_FST <- DFSP_mol_confirmed %>% filter (FST == "Yes")

# Convert survival months into numeric variables
DFSP$OS.time <- as.numeric(DFSP$OS.time)
DFSP$RFS.time <- as.numeric(DFSP$RFS.time)

# Produce survival object
survival_obj_OS <- Surv(DFSP$OS.time, DFSP$OS.status == "1:DECEASED")
survival_obj_RFS <- Surv(DFSP$RFS.time, DFSP$RFS.status== "1:Recurrence")

# Survival object from FST subset
survival_obj <- Surv(DFSP_mol_confirmed_FST$OS.time, DFSP_mol_confirmed_FST$OS.status == "1:DECEASED")
survival_obj <- Surv(DFSP_mol_confirmed_FST$RFS.time, DFSP_mol_confirmed_FST$RFS.status== "1:Recurrence")


cox_ph_model <- coxph(survival_obj ~ FST + risk, data = DFSP_mol_confirmed)
summary(cox_ph_model)
ggforest(cox_ph_model, data = DFSP_mol_confirmed, refLabel = "Reference")
# exp(coef) = hazard ratio
# p-values (Pr(>|z|)) indicate the statistical significance of the predictor variables
# The log-likelihood, likelihood ratio test, Wald test, and Score (log-rank) test are global test statistics 
# that assess the overall significance of the model. 
# <0.05 indicates that the model as a whole is statistically significant.

#Calculate CI of HR
hazard_ratio_conf_int <- exp(confint(cox_ph_model))

km_curves_OS <- survfit(survival_obj_OS ~ MSI.cat, data = DFSP)
km_curves_RFS <- survfit(survival_obj_RFS ~ MSI.cat, data = DFSP)

#FST subset
km_curves <- survfit(survival_obj_OS ~ CDKN2A.status, data = DFSP_mol_confirmed_FST) 

library(survminer)
surv_plot_OS <- ggsurvplot(km_curves_OS, data = DFSP, conf.int = TRUE, pval = TRUE,
                        pval.method = TRUE, risk.table = TRUE)
surv_plot_RFS <- ggsurvplot(km_curves_RFS, data = DFSP, conf.int = TRUE, pval = TRUE,
                           pval.method = TRUE, risk.table = TRUE)

#FST subset
surv_plot <- ggsurvplot(km_curves, data = DFSP_mol_confirmed_FST, conf.int = TRUE, pval = TRUE,
                        pval.method = TRUE, risk.table = TRUE)
surv_plot_OS
surv_plot_RFS
