# Load necessary libraries
library(dplyr)

# Load and filter the data
data <- micro_analysis_with_MRI_200225_final %>%
  filter(Micro_YN == 1, Diagnosis_number %in% c(1, 3, 4))

# Convert Diagnosis_number to factor
data$Diagnosis_number <- as.factor(data$Diagnosis_number)

# List of categorical variables
categorical_vars <- c("Antidepressant_YN", "Benzodiazepines_YN", "Cholinesterase_YN", "Sex")

# Initialize results dataframe
chisq_results <- data.frame(
  Variable = character(),
  ChiSq = numeric(),
  df = numeric(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each categorical variable
for (var in categorical_vars) {
  if (!var %in% names(data)) {
    warning(paste("Skipping:", var, "- not found in dataset"))
    next
  }
  
  # Create contingency table
  tbl <- table(data[[var]], data$Diagnosis_number)
  
  # Print counts
  cat("\n========================================\n")
  cat("Variable:", var, "\n")
  print(tbl)
  
  # Check if table is valid
  if (any(rowSums(tbl) == 0) || any(colSums(tbl) == 0)) {
    warning(paste("Skipping:", var, "- table has empty rows or columns"))
    next
  }
  
  # Perform chi-squared test
  test <- chisq.test(tbl)
  
  # Append to results
  chisq_results <- rbind(chisq_results, data.frame(
    Variable = var,
    ChiSq = test$statistic,
    df = test$parameter,
    PValue = test$p.value
  ))
}

# Format p-values
chisq_results$PValue <- formatC(chisq_results$PValue, format = "e", digits = 5)

# Print chi-squared results summary
cat("\n========================================\n")
cat("Chi-squared Test Results Summary:\n")
print(chisq_results)

# Optional: Save results
# write.csv(chisq_results, "C:/Users/jacka/Downloads/chisq_results_with_counts.csv", row.names = FALSE)
