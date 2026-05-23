# Load libraries
library(dplyr)
library(broom)
library(ggplot2)

# Filter for groups 3 and 4
group34_data <- micro_analysis_with_MRI_200225_final %>%
  filter(Diagnosis_number %in% c(3, 4))

# Updated functional connectivity variables (your provided list)
func_variables <- c(
  "VISnbm",
  "VANnbm",
  "TEMPnbm",
  "DMNnbm",
  "FPNnbm",
  "LIMnbm",
  "DANnbm",
  "SMNnbm"
  
)

# Sleep slowing variables
sleep_vars <- c(
  "Central_Absolute_REM_SFratio1",
  "Frontal_Absolute_REM_SFratio1",
  "Occipital_Absolute_REM_SFratio1"
)

# Create plots folder
if (!dir.exists("plots")) dir.create("plots")

# Storage for results
results <- data.frame()

cat("Starting analysis...\n")

# Loop over sleep and functional variables
for (sleep_var in sleep_vars) {
  for (func_var in func_variables) {
    
    cat("\nAnalyzing combination:", sleep_var, "and", func_var, "\n")
    
    # Filter valid data
    valid_data <- group34_data %>%
      filter(!is.na(.data[[sleep_var]]), !is.na(.data[[func_var]]),
             !is.na(Age_at_study_date), !is.na(Sex)) %>%
      mutate(
        sleep = as.numeric(.data[[sleep_var]]),
        func = as.numeric(.data[[func_var]])
      ) %>%
      filter(!is.na(sleep), !is.na(func))
    
    cat("  Initial n after NA filtering:", nrow(valid_data), "\n")
    
    if (nrow(valid_data) >= 10) {
      model <- lm(sleep ~ func + Age_at_study_date + Sex, data = valid_data)
      model_summary <- summary(model)
      
      beta <- model_summary$coefficients["func", "Estimate"]
      p_value <- model_summary$coefficients["func", "Pr(>|t|)"]
      r_squared <- model_summary$r.squared
      
      # Save results
      results <- rbind(results, data.frame(
        Sleep_Variable = sleep_var,
        Functional_Connectivity_Variable = func_var,
        Beta = beta,
        P_Value = p_value,
        R_Squared = r_squared
      ))
      
      # Plot
      plot_title <- paste0(sleep_var, " vs ", func_var)
      p <- ggplot(valid_data, aes(x = func, y = sleep)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(
          title = plot_title,
          subtitle = paste0("Beta = ", round(beta, 3),
                            ", p = ", signif(p_value, 3),
                            ", R² = ", round(r_squared, 3)),
          x = func_var,
          y = sleep_var
        ) +
        theme_minimal(base_size = 14)
      
      print(p)
      
      plot_filename <- paste0("plots/", sleep_var, "_", func_var, "_plot.png")
      ggsave(plot_filename, p, width = 7, height = 5)
    } else {
      cat("  Skipped due to insufficient data (<10 participants).\n")
    }
  }
}

# Output results
if (nrow(results) > 0) {
  results <- results %>% arrange(P_Value)
  print(results)
  write.csv(results, "func_sleep_model_results.csv", row.names = FALSE)
  cat("Results saved to 'func_sleep_model_results.csv'.\n")
} else {
  cat("No valid results generated.\n")
}
