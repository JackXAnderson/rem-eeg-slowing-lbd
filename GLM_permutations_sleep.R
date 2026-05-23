# Load necessary libraries
library(dplyr)
library(boot)        # For permutation testing
library(multcomp)    # (Optional: Was used before for glht, not needed now)

# Load and clean data
data <- micro_analysis_with_MRI_200225_final

# Filter for Micro_YN == 1 and groups 1, 3, 4
data <- data %>%
  filter(Micro_YN == 1, Diagnosis_number %in% c(1, 3, 4))

# Convert to factors after filtering
data$Sex <- as.factor(data$Sex)
data$Diagnosis_number <- as.factor(data$Diagnosis_number)

# Define the list of sleep variables
sleep_variables <- c(
  "Total_time_in_bed_(min)", "Sleep_latency_(min)", "Total_sleep_time_(min)", 
  "REM_latency_(min)", "Wake_after_sleep_onset_(min)", "_NREM_sleep_(min)", 
  "Sleep_efficiency", "REM_sleep_(min)", "Total_AHI_(events/hr)", 
  "Total_RDI_(events/hr)_", "RERA_index_", "RDI_in_NREM", "RDI_in_REM", 
  "Minimum_SpO2_during_sleep_(%", "Total_arousal_index_(arousals/hr)", 
  "PLM_arousal_index_(arousals/hr)_", "WASO", "N1_Duration", "N1_%_TST", 
  "N2_Duration", "N2_%_TST", "N3_Duration", "N3_%_TST", 
  "REM_Duration", "REM_%_TST"
)

# Function to compute F-statistic for GLM
compute_f_statistic <- function(data, measure) {
  if (!measure %in% names(data) || all(is.na(data[[measure]]))) {
    return(NA)
  }
  
  data[[measure]] <- as.numeric(data[[measure]])
  formula <- as.formula(paste0("`", measure, "` ~ Age_at_study_date + Sex + Diagnosis_number"))
  model <- glm(formula, data = data, family = gaussian())
  anova_results <- anova(model)
  
  if ("Diagnosis_number" %in% rownames(anova_results)) {
    f_statistic <- anova_results["Diagnosis_number", "F"]
  } else {
    f_statistic <- NA
  }
  
  return(f_statistic)
}

# Function for residual-based pairwise permutation post hoc
residual_posthoc_test <- function(residuals, groups, g1, g2, nperm = 10000) {
  idx <- which(groups %in% c(g1, g2))
  if (length(unique(groups[idx])) < 2 || length(idx) < 6) return(NULL)
  
  group <- droplevels(factor(groups[idx], levels = c(g1, g2)))
  obs_diff <- abs(mean(residuals[idx][group == g1], na.rm = TRUE) - mean(residuals[idx][group == g2], na.rm = TRUE))
  
  perm_diffs <- replicate(nperm, {
    shuffled <- sample(group)
    abs(mean(residuals[idx][shuffled == g1], na.rm = TRUE) - mean(residuals[idx][shuffled == g2], na.rm = TRUE))
  })
  
  p_val <- mean(perm_diffs >= obs_diff, na.rm = TRUE)
  data.frame(Comparison = paste(g1, "vs", g2), Difference = obs_diff, P_value = p_val)
}

# Initialize result dataframes
glm_results_sleep <- data.frame(
  Variable = character(),
  FStatistic = numeric(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Post hoc results in wide format
posthoc_list <- list()

# Loop through each sleep variable
for (measure in sleep_variables) {
  if (!measure %in% names(data) || all(is.na(data[[measure]]))) {
    warning(paste("Skipping variable:", measure, "- not found or all values are NA"))
    next
  }
  
  data[[measure]] <- as.numeric(data[[measure]])
  observed_f_statistic <- compute_f_statistic(data, measure)
  if (is.na(observed_f_statistic)) next
  
  # Permutation test for GLM
  n_permutations <- 1000
  permuted_f_statistics <- numeric(n_permutations)
  for (i in 1:n_permutations) {
    permuted_data <- data
    permuted_data$Diagnosis_number <- sample(permuted_data$Diagnosis_number)
    permuted_f_statistics[i] <- compute_f_statistic(permuted_data, measure)
  }
  p_value <- mean(permuted_f_statistics >= observed_f_statistic, na.rm = TRUE)
  
  # Save GLM result
  glm_results_sleep <- rbind(glm_results_sleep, data.frame(
    Variable = measure,
    FStatistic = observed_f_statistic,
    PValue = p_value
  ))
  
  # Residual-based post hoc test
  covariate_formula <- as.formula(paste0("`", measure, "` ~ Age_at_study_date + Sex"))
  lm_model <- lm(covariate_formula, data = data)
  residuals <- resid(lm_model)
  groups <- droplevels(data$Diagnosis_number)
  
  posthoc_row <- data.frame(
    Variable = measure,
    `1_vs_3` = NA,
    `1_vs_4` = NA,
    `3_vs_4` = NA,
    stringsAsFactors = FALSE
  )
  
  pairwise_comparisons <- list(c("1", "3"), c("1", "4"), c("3", "4"))
  for (pair in pairwise_comparisons) {
    g1 <- pair[1]; g2 <- pair[2]
    if (g1 %in% levels(groups) && g2 %in% levels(groups)) {
      res <- residual_posthoc_test(residuals, groups, g1, g2)
      if (!is.null(res)) {
        col_name <- paste0(g1, "_vs_", g2)
        posthoc_row[[col_name]] <- res$P_value
      }
    }
  }
  
  posthoc_list[[measure]] <- posthoc_row
}

# Combine all post hoc results into one wide dataframe
posthoc_results_sleep <- do.call(rbind, posthoc_list)

# Format for display
glm_results_sleep$PValue <- formatC(glm_results_sleep$PValue, format = "e", digits = 5)
posthoc_results_sleep[, -1] <- lapply(posthoc_results_sleep[, -1], function(x) {
  x_numeric <- suppressWarnings(as.numeric(x))
  formatC(x_numeric, format = "e", digits = 5)
})

# Print results
print(glm_results_sleep)
print(posthoc_results_sleep)

# Optional: save results
# write.csv(glm_results_sleep, "C:/Users/jacka/Downloads/glm_results_sleep.csv", row.names = FALSE)
# write.csv(posthoc_results_sleep, "C:/Users/jacka/Downloads/posthoc_results_sleep.csv", row.names = FALSE)
