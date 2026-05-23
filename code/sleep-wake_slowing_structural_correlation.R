library(dplyr)
library(broom)

# Filter for groups 3 and 4
group34_data <- micro_analysis_with_MRI_200225_final %>%
  filter(Diagnosis_number %in% c(3, 4))

# Define occipital volume variables
occipital_vars <- c(
  "Total_bankssts_volume_Adjusted", "Total_caudalanteriorcingulate_volume_Adjusted", 
  "Total_caudalmiddlefrontal_volume_Adjusted", "Total_cuneus_volume_Adjusted", 
  "Total_entorhinal_volume_Adjusted", "Total_fusiform_volume_Adjusted", 
  "Total_inferiorparietal_volume_Adjusted", "Total_inferiortemporal_volume_Adjusted", 
  "Total_isthmuscingulate_volume_Adjusted", "Total_lateraloccipital_volume_Adjusted", 
  "Total_lateralorbitofrontal_volume_Adjusted", "Total_lingual_volume_Adjusted", 
  "Total_medialorbitofrontal_volume_Adjusted", "Total_middletemporal_volume_Adjusted", 
  "Total_parahippocampal_volume_Adjusted", "Total_paracentral_volume_Adjusted", 
  "Total_parsopercularis_volume_Adjusted", "Total_parsorbitalis_volume_Adjusted", 
  "Total_parstriangularis_volume_Adjusted", "Total_pericalcarine_volume_Adjusted", 
  "Total_postcentral_volume_Adjusted", "Total_posteriorcingulate_volume_Adjusted", 
  "Total_precentral_volume_Adjusted", "Total_precuneus_volume_Adjusted", 
  "Total_rostralanteriorcingulate_volume_Adjusted", "Total_rostralmiddlefrontal_volume_Adjusted", 
  "Total_superiorfrontal_volume_Adjusted", "Total_superiorparietal_volume_Adjusted", 
  "Total_superiortemporal_volume_Adjusted", "Total_supramarginal_volume_Adjusted", 
  "Total_frontalpole_volume_Adjusted", "Total_temporalpole_volume_Adjusted", 
  "Total_transversetemporal_volume_Adjusted", "Total_insula_volume_Adjusted"
)

# Define sleep variables
sleep_vars <- c(
  # Central Absolute NREM
  "Central_Absolute_NREM_Delta", "Central_Absolute_NREM_Theta", "Central_Absolute_NREM_Alpha",
  "Central_Absolute_NREM_Sigma", "Central_Absolute_NREM_Beta", "Central_Absolute_NREM_SFratio1", "Central_Absolute_NREM_SFratio2",
  
  # Frontal Absolute NREM
  "Frontal_Absolute_NREM_Delta", "Frontal_Absolute_NREM_Theta", "Frontal_Absolute_NREM_Alpha",
  "Frontal_Absolute_NREM_Sigma", "Frontal_Absolute_NREM_Beta", "Frontal_Absolute_NREM_SFratio1", "Frontal_Absolute_NREM_SFratio2",
  
  # Occipital Absolute NREM
  "Occipital_Absolute_NREM_Delta", "Occipital_Absolute_NREM_Theta", "Occipital_Absolute_NREM_Alpha",
  "Occipital_Absolute_NREM_Sigma", "Occipital_Absolute_NREM_Beta", "Occipital_Absolute_NREM_SFratio1", "Occipital_Absolute_NREM_SFratio2",
  
  # Central Absolute REM
  "Central_Absolute_REM_Delta", "Central_Absolute_REM_Theta", "Central_Absolute_REM_Alpha",
  "Central_Absolute_REM_Sigma", "Central_Absolute_REM_Beta", "Central_Absolute_REM_SFratio1", "Central_Absolute_REM_SFratio2",
  
  # Frontal Absolute REM
  "Frontal_Absolute_REM_Delta", "Frontal_Absolute_REM_Theta", "Frontal_Absolute_REM_Alpha",
  "Frontal_Absolute_REM_Sigma", "Frontal_Absolute_REM_Beta", "Frontal_Absolute_REM_SFratio1", "Frontal_Absolute_REM_SFratio2",
  
  # Occipital Absolute REM
  "Occipital_Absolute_REM_Delta", "Occipital_Absolute_REM_Theta", "Occipital_Absolute_REM_Alpha",
  "Occipital_Absolute_REM_Sigma", "Occipital_Absolute_REM_Beta", "Occipital_Absolute_REM_SFratio1", "Occipital_Absolute_REM_SFratio2"
  
)

# Initialize results dataframe
results <- data.frame()

cat("Starting model loop (occipital x sleep)...\n")

# Loop through sleep-occipital combinations
for (sleep_var in sleep_vars) {
  for (occipital_var in occipital_vars) {
    
    cat("Running model for:", sleep_var, "~", occipital_var, "\n")
    
    # Prepare data
    valid_data <- group34_data %>%
      filter(
        !is.na(.data[[sleep_var]]),
        !is.na(.data[[occipital_var]]),
        !is.na(Age_at_study_date),
        !is.na(Sex)
      ) %>%
      mutate(
        sleep = as.numeric(.data[[sleep_var]]),
        occipital = as.numeric(.data[[occipital_var]])
      ) %>%
      filter(!is.na(sleep), !is.na(occipital))
    
    cat("Valid rows:", nrow(valid_data), "\n")
    
    # Run model if enough data
    if (nrow(valid_data) >= 10) {
      model <- lm(sleep ~ occipital + Age_at_study_date + Sex, data = valid_data)
      model_summary <- summary(model)
      
      beta <- model_summary$coefficients["occipital", "Estimate"]
      p_value <- model_summary$coefficients["occipital", "Pr(>|t|)"]
      r_squared <- model_summary$r.squared
      
      results <- rbind(results, data.frame(
        Sleep_Variable = sleep_var,
        Occipital_Variable = occipital_var,
        Beta = beta,
        P_Value = p_value,
        R_Squared = r_squared
      ))
    }
  }
}

# Final output
if (nrow(results) > 0) {
  results <- results %>% arrange(P_Value)
  print(results)
  write.csv(results, "model_results_occipital_sleep.csv", row.names = FALSE)
  cat("Results saved to 'model_results_occipital_sleep.csv'\n")
} else {
  cat("No valid models were run due to insufficient data.\n")
}
