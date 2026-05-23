library(dplyr)
library(broom)
library(ggplot2)

# ── Data ──────────────────────────────────────────────────────────────────────
data <- micro_analysis_with_MRI_200225_final

# ── REM slowing ratio outcomes ────────────────────────────────────────────────
slowing_vars <- c(
  "Occipital_Absolute_REM_SFratio1"
)

# ── Disease burden predictors ─────────────────────────────────────────────────
# CAF_Severity: DLB only (group 4) — cognitive fluctuations only assessed in DLB
# All others: PD + DLB (groups 3 & 4)
burden_vars <- c(
  "ESS_macro",                      # Epworth Sleepiness Scale
  "MDUPDRS_Section_I_Total_macro",  # MDS-UPDRS Part I (non-motor)
  "MDUPDRS_Section_III_Total_macro",# MDS-UPDRS Part III (motor)
  "REMSBD_Total",                   # RBD Screening Questionnaire
  "CAF_Severity"                    # Cognitive fluctuations — DLB only
)

# Groups for each predictor
burden_groups <- list(
  ESS_macro                       = c(3, 4),
  MDUPDRS_Section_I_Total_macro   = c(3, 4),
  MDUPDRS_Section_III_Total_macro = c(3, 4),
  REMSBD_Total                    = c(3, 4),
  CAF_Severity                    = c(4)       # DLB only
)

# ── Helper: nice label for plots ──────────────────────────────────────────────
var_labels <- c(
  ESS_macro                        = "Epworth Sleepiness Scale",
  MDUPDRS_Section_I_Total_macro    = "MDS-UPDRS Part I (Non-Motor)",
  MDUPDRS_Section_III_Total_macro  = "MDS-UPDRS Part III (Motor)",
  REMSBD_Total                     = "RBD Screening Questionnaire",
  CAF_Severity                     = "CAF Severity (DLB only)",
  Central_Absolute_REM_SFratio1    = "Central REM Slowing Ratio",
  Frontal_Absolute_REM_SFratio1    = "Frontal REM Slowing Ratio",
  Occipital_Absolute_REM_SFratio1  = "Occipital REM Slowing Ratio"
)

# ── Storage ───────────────────────────────────────────────────────────────────
results <- data.frame()

# ── Optional: save all plots to a single PDF ──────────────────────────────────
# Uncomment the pdf() line below and the dev.off() line at the very end
# to save all 15 plots to one file instead of printing to the RStudio plots pane.
# pdf("disease_burden_correlations.pdf", width = 8, height = 6)

# ── Main loop ─────────────────────────────────────────────────────────────────
for (burden_var in burden_vars) {
  
  # Select appropriate groups for this predictor
  groups_to_use <- burden_groups[[burden_var]]
  working_data  <- data %>%
    filter(Diagnosis_number %in% groups_to_use) %>%
    mutate(Group = factor(Diagnosis_number,
                          levels = c(3, 4),
                          labels = c("PD", "DLB")))
  
  for (slowing_var in slowing_vars) {
    
    # Complete cases + 3SD outlier removal
    valid_data <- working_data %>%
      filter(
        !is.na(.data[[slowing_var]]),
        !is.na(.data[[burden_var]]),
        !is.na(Age_at_study_date),
        !is.na(Sex)
      ) %>%
      mutate(
        outcome   = as.numeric(.data[[slowing_var]]),
        predictor = as.numeric(.data[[burden_var]]),
        out_z     = as.numeric(scale(outcome)),
        pred_z    = as.numeric(scale(predictor))
      ) %>%
      filter(abs(out_z) <= 3, abs(pred_z) <= 3)
    
    n <- nrow(valid_data)
    
    if (n < 6) {
      cat("Skipping", burden_var, "vs", slowing_var,
          "— insufficient n after filtering (n =", n, ")\n")
      next
    }
    
    # ── Linear model: slowing ratio ~ predictor + age + sex ──────────────────
    model         <- lm(outcome ~ predictor + Age_at_study_date + Sex,
                        data = valid_data)
    model_summary <- summary(model)
    
    beta      <- model_summary$coefficients["predictor", "Estimate"]
    se        <- model_summary$coefficients["predictor", "Std. Error"]
    p_value   <- model_summary$coefficients["predictor", "Pr(>|t|)"]
    r_squared <- model_summary$r.squared
    
    # ── Partial correlation (predictor ~ slowing after removing age+sex) ──────
    resid_outcome   <- resid(lm(outcome   ~ Age_at_study_date + Sex, data = valid_data))
    resid_predictor <- resid(lm(predictor ~ Age_at_study_date + Sex, data = valid_data))
    partial_r       <- cor(resid_outcome, resid_predictor, use = "complete.obs")
    
    # ── Save results ──────────────────────────────────────────────────────────
    results <- rbind(results, data.frame(
      Burden_Variable  = var_labels[burden_var],
      Slowing_Variable = var_labels[slowing_var],
      Groups           = ifelse(length(groups_to_use) > 1, "PD + DLB", "DLB only"),
      n                = n,
      Beta             = round(beta,      4),
      SE               = round(se,        4),
      Partial_r        = round(partial_r, 3),
      P_Value          = round(p_value,   4),
      R_Squared        = round(r_squared, 3),
      stringsAsFactors = FALSE
    ))
    
    # ── Scatter plot ──────────────────────────────────────────────────────────
    # Colour by group only when both PD and DLB present
    if (length(groups_to_use) > 1) {
      p <- ggplot(valid_data, aes(x = predictor, y = outcome, color = Group)) +
        scale_color_manual(values = c("PD" = "#2196F3", "DLB" = "#E53935"))
    } else {
      p <- ggplot(valid_data, aes(x = predictor, y = outcome)) +
        scale_color_manual(values = c("DLB" = "#E53935"))
    }
    
    p <- p +
      geom_point(size = 3, alpha = 0.85) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      labs(
        title    = paste0(var_labels[slowing_var], "\nvs ", var_labels[burden_var]),
        subtitle = paste0("β = ", round(beta, 3),
                          ", partial r = ", round(partial_r, 3),
                          ", p = ", signif(p_value, 3),
                          ", n = ", n),
        x        = var_labels[burden_var],
        y        = var_labels[slowing_var],
        color    = "Group"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    
    print(p)
  }
}

# ── Summary table ─────────────────────────────────────────────────────────────
if (nrow(results) > 0) {
  
  # Sort by p-value
  results <- results %>% arrange(P_Value)
  
  cat("\n=====================================================================\n")
  cat(" DISEASE BURDEN ~ REM SLOWING RATIO — SUMMARY\n")
  cat(" Covariate-adjusted (age + sex) | 3SD outliers removed\n")
  cat(" Raw p-values (no correction)\n")
  cat("=====================================================================\n")
  print(results, row.names = FALSE)
  
  # Flag nominally significant findings
  sig <- results %>% filter(P_Value < 0.05)
  if (nrow(sig) > 0) {
    cat("\n--- Nominally significant results (p < 0.05) ---\n")
    print(sig[, c("Burden_Variable","Slowing_Variable","Groups","n",
                  "Beta","Partial_r","P_Value")],
          row.names = FALSE)
  } else {
    cat("\nNo nominally significant results (p < 0.05).\n")
  }
  
} else {
  cat("No valid models completed.\n")
}

# Uncomment if you opened a pdf() device above:
# dev.off()
