library(dplyr)

data <- micro_analysis_with_MRI_200225_final

# ── Variable name for AHI — adjust if column name differs in your dataframe ───
# Using backticks to handle special characters in the column name
ahi_var <- "Total_AHI_(events/hr)"

# ── Helper: run ESS ~ AHI + Age + Sex regression and print clean summary ──────
run_ess_ahi <- function(df, group_label) {
  
  valid <- df %>%
    filter(
      !is.na(ESS_macro),
      !is.na(.data[[ahi_var]]),
      !is.na(Age_at_study_date),
      !is.na(Sex)
    ) %>%
    mutate(
      ESS = as.numeric(ESS_macro),
      AHI = as.numeric(.data[[ahi_var]]),
      ESS_z = as.numeric(scale(ESS)),
      AHI_z = as.numeric(scale(AHI))
    ) %>%
    filter(abs(ESS_z) <= 3, abs(AHI_z) <= 3)  # 3SD outlier removal
  
  n <- nrow(valid)
  cat("\n─────────────────────────────────────────────\n")
  cat(" Group:", group_label, "| n =", n, "\n")
  cat("─────────────────────────────────────────────\n")
  
  if (n < 6) {
    cat(" Insufficient n — skipping.\n")
    return(invisible(NULL))
  }
  
  # Full model: ESS ~ AHI + Age + Sex
  model_full <- lm(ESS ~ AHI + Age_at_study_date + Sex, data = valid)
  s_full     <- summary(model_full)
  
  # Null model (age + sex only) for comparison of R²
  model_null <- lm(ESS ~ Age_at_study_date + Sex, data = valid)
  s_null     <- summary(model_null)
  
  # Partial r for AHI
  resid_ess <- resid(lm(ESS ~ Age_at_study_date + Sex, data = valid))
  resid_ahi <- resid(lm(AHI ~ Age_at_study_date + Sex, data = valid))
  partial_r  <- cor(resid_ess, resid_ahi, use = "complete.obs")
  
  # Extract AHI coefficient
  beta    <- s_full$coefficients["AHI", "Estimate"]
  se      <- s_full$coefficients["AHI", "Std. Error"]
  p_value <- s_full$coefficients["AHI", "Pr(>|t|)"]
  r2_full <- s_full$r.squared
  r2_null <- s_null$r.squared
  delta_r2 <- r2_full - r2_null   # variance in ESS uniquely explained by AHI
  
  cat(sprintf(" AHI beta:       %+.3f (SE = %.3f)\n", beta, se))
  cat(sprintf(" AHI partial r:  %.3f\n", partial_r))
  cat(sprintf(" AHI p-value:    %.4f\n", p_value))
  cat(sprintf(" R² full model:  %.3f\n", r2_full))
  cat(sprintf(" R² age+sex only:%.3f\n", r2_null))
  cat(sprintf(" ΔR² for AHI:    %.3f  (%.1f%% of ESS variance)\n",
              delta_r2, delta_r2 * 100))
  cat("\n Full model summary:\n")
  print(s_full$coefficients, digits = 3)
  
  # Return results invisibly for any downstream use
  invisible(list(
    group     = group_label,
    n         = n,
    beta      = beta,
    se        = se,
    partial_r = partial_r,
    p_value   = p_value,
    r2_full   = r2_full,
    r2_null   = r2_null,
    delta_r2  = delta_r2
  ))
}

# ── Run in PD + DLB combined ──────────────────────────────────────────────────
pd_dlb <- data %>% filter(Diagnosis_number %in% c(3, 4))
res_combined <- run_ess_ahi(pd_dlb, "PD + DLB (groups 3 & 4)")

# ── Run in DLB only ───────────────────────────────────────────────────────────
dlb_only <- data %>% filter(Diagnosis_number == 4)
res_dlb <- run_ess_ahi(dlb_only, "DLB only (group 4)")

# ── Run in PD only (for comparison) ──────────────────────────────────────────
pd_only <- data %>% filter(Diagnosis_number == 3)
res_pd <- run_ess_ahi(pd_only, "PD only (group 3)")

cat("\n=====================================================\n")
cat(" INTERPRETATION GUIDE\n")
cat("=====================================================\n")
cat(" If AHI p > 0.05 and ΔR² is small (e.g. <5%):\n")
cat("   → Somnolence is NOT well explained by apnea severity\n")
cat("   → Supports neurodegeneration-driven somnolence argument\n\n")
cat(" If AHI p < 0.05 and ΔR² is substantial (e.g. >10%):\n")
cat("   → Apnea explains meaningful ESS variance\n")
cat("   → Consider ESS residuals for downstream correlations\n")
cat("=====================================================\n")