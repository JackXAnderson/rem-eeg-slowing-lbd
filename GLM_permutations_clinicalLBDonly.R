# Load necessary libraries
library(dplyr)
library(boot)

# Load and filter data
data <- micro_analysis_with_MRI_200225_final %>%
  filter(Micro_YN == 1, Diagnosis_number %in% c(3, 4))

# Convert Diagnosis_number to factor
data$Diagnosis_number <- as.factor(data$Diagnosis_number)

# Define list of variables
demographic_variables <- c(
  "Age_at_study_date", "Sex", "RBD_PSG", "ESS", "DDE", "CAF_Severity", "Weight", "Height", "BMI", 
  "Evening_SBP", "Evening_DBP", "Morning_SBP", "Morning_DBP", "ESS_macro", "SCOPA_S_Nocturl_", "SCOPA_S_Daytime", "MOCA_total_macro",
  "MMSE_result_macro",
  "MDUPDRS_Section_I_Total_macro", "MDUPDRS_Section_II_Total_macro", 
  "MDUPDRS_Section_III_Total_macro", "MDUPDRS_Section_IV_Total_macro", "REMSBD_Total", 
  "Years_of_education", "RT", "MMSE_result_macro", 
  "Verbal_fluency_letters", "Verbal_fluency_animals", 
  "Digit_Span___Forward", "Digit_Span___Backward", "Longest_Digit_Span___Forward", 
  "Longest_Digit_Span___Backwards", "Digit_span_ASS_macro", "ASS_1", "ASS_2",
  "ASS_%", 
  "Trails_A_Z_score_macro", 
  "Trails_B_Z_score_macro", "Stroop__1_raw", "Stroop_1_ASS_macro", "Stroop_2_raw", 
  "Stroop_2_ASS_macro", "Stroop_3_raw", "Stroop_3_ASS_macro", "Stroop_4_raw", 
  "Stroop_4_ASS_macro", 
  "RAVLT_A1_5_Total_ASS", "RAVLT_B1", "RAVLT_A6", "RAVLT_A7", "RAVLT_A6_A5",  "RAVLT_Delayed_recall_ASS", 
  "Clock_Drawing___Total", "Boston_ming_Spon_Correct", "Grand_total",
  "Total_Sniffin_Sticks_macro", "Contrast_Discrimition",
  "Total_Colour_macro"
)

# Function to compute F-statistic
compute_f_statistic <- function(data, measure) {
  data[[measure]] <- as.numeric(data[[measure]])
  if (length(unique(data$Diagnosis_number)) < 2) return(NA)
  formula <- as.formula(paste0("`", measure, "` ~ Diagnosis_number"))
  model <- glm(formula, data = data, family = gaussian())
  anova_results <- anova(model)
  if ("Diagnosis_number" %in% rownames(anova_results)) {
    return(anova_results["Diagnosis_number", "F"])
  } else {
    return(NA)
  }
}

# Initialize result dataframe
glm_results_demo <- data.frame(
  Variable = character(),
  FStatistic = numeric(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Loop through variables
for (measure in demographic_variables) {
  if (!measure %in% names(data) || all(is.na(data[[measure]]))) {
    warning(paste("Skipping:", measure, "- not found or all NA"))
    next
  }
  
  temp <- data[!is.na(data[[measure]]), ]
  if (length(unique(temp$Diagnosis_number)) < 2) {
    warning(paste("Skipping:", measure, "- only one Diagnosis group with data"))
    next
  }
  
  temp[[measure]] <- suppressWarnings(as.numeric(temp[[measure]]))
  if (all(is.na(temp[[measure]]))) {
    warning(paste("Skipping:", measure, "- not numeric or all NA after coercion"))
    next
  }
  
  # F-statistic
  f_stat <- compute_f_statistic(temp, measure)
  if (is.na(f_stat)) next
  
  # Permutation test
  n_permutations <- 1000
  permuted_f_statistics <- numeric(n_permutations)
  for (i in 1:n_permutations) {
    permuted_temp <- temp
    permuted_temp$Diagnosis_number <- sample(permuted_temp$Diagnosis_number)
    permuted_f_statistics[i] <- compute_f_statistic(permuted_temp, measure)
  }
  p_value <- mean(permuted_f_statistics >= f_stat, na.rm = TRUE)
  
  # Save result
  glm_results_demo <- rbind(glm_results_demo, data.frame(
    Variable = measure,
    FStatistic = f_stat,
    PValue = p_value
  ))
}

# Format results
glm_results_demo$PValue <- formatC(glm_results_demo$PValue, format = "e", digits = 5)

# Print results
print(glm_results_demo)

# Optional: Save to CSV
# write.csv(glm_results_demo, "C:/Users/jacka/Downloads/glm_results_group3_vs_4.csv", row.names = FALSE)
