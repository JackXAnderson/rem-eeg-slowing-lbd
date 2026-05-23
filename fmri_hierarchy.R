library(dplyr)

data <- micro_analysis_with_MRI_200225_final

# ── Variables ──────────────────────────────────────────────────────────────────
outcome   <- "Occipital_Absolute_REM_SFratio1"
ppn_var   <- "totalTHAppn"
nbm_vars  <- c("VISnbm", "VANnbm", "FPNnbm", "DMNnbm")  # your four FDR-sig networks

# ── Filter to MRI subsample (PD + DLB only, non-missing MRI) ──────────────────
mri_data <- data %>%
  filter(
    Diagnosis_number %in% c(3, 4),
    !is.na(.data[[outcome]]),
    !is.na(.data[[ppn_var]]),
    !is.na(Age_at_study_date),
    !is.na(Sex)
  ) %>%
  mutate(
    outcome_z = as.numeric(scale(.data[[outcome]])),
    ppn_z     = as.numeric(scale(.data[[ppn_var]]))
  ) %>%
  filter(abs(outcome_z) <= 3, abs(ppn_z) <= 3)  # 3SD outlier removal

cat("══════════════════════════════════════════════════════════\n")
cat(" MRI SUBSAMPLE\n")
cat(sprintf(" Total n = %d\n", nrow(mri_data)))
cat(" Group breakdown:\n")
print(table(mri_data$Diagnosis_number))
cat("══════════════════════════════════════════════════════════\n\n")

# ── Helper: print a clean model comparison block ──────────────────────────────
print_model_block <- function(m0, m1, m2, labels) {
  
  s0 <- summary(m0)
  s1 <- summary(m1)
  s2 <- summary(m2)
  
  cat(sprintf("  %-30s R² = %.3f\n", labels[1], s0$r.squared))
  cat(sprintf("  %-30s R² = %.3f  ΔR² = %.3f  (%.1f%%)\n",
              labels[2], s1$r.squared,
              s1$r.squared - s0$r.squared,
              (s1$r.squared - s0$r.squared) * 100))
  cat(sprintf("  %-30s R² = %.3f  ΔR² = %.3f  (%.1f%%)\n\n",
              labels[3], s2$r.squared,
              s2$r.squared - s1$r.squared,
              (s2$r.squared - s1$r.squared) * 100))
  
  # Coefficient table for final model
  cat("  Final model coefficients:\n")
  coef_tbl <- s2$coefficients
  for (row in rownames(coef_tbl)) {
    cat(sprintf("    %-35s  β = %+.4f  p = %.4f\n",
                row,
                coef_tbl[row, "Estimate"],
                coef_tbl[row, "Pr(>|t|)"]))
  }
  cat("\n")
}

# ── Main loop: four NBM networks ──────────────────────────────────────────────
results_summary <- data.frame()

for (nbm in nbm_vars) {
  
  # Complete cases for this NBM variable
  df <- mri_data %>%
    filter(!is.na(.data[[nbm]])) %>%
    mutate(nbm_z = as.numeric(scale(.data[[nbm]]))) %>%
    filter(abs(nbm_z) <= 3)
  
  n <- nrow(df)
  
  cat("══════════════════════════════════════════════════════════\n")
  cat(sprintf(" NBM NETWORK: %s  |  n = %d\n", nbm, n))
  cat("══════════════════════════════════════════════════════════\n\n")
  
  if (n < 10) {
    cat("  Insufficient n — skipping.\n\n")
    next
  }
  
  # ── Base model: age + sex ──────────────────────────────────────────────────
  m_base <- lm(as.formula(paste(outcome, "~ Age_at_study_date + Sex")), data = df)
  
  # ── ORDER A: PPN first, then NBM ──────────────────────────────────────────
  cat(" ORDER A: Age+Sex → +PPN → +NBM\n")
  cat(" (Does NBM explain variance OVER AND ABOVE PPN?)\n\n")
  
  m_A1 <- lm(as.formula(paste(outcome, "~ Age_at_study_date + Sex +", ppn_var)), data = df)
  m_A2 <- lm(as.formula(paste(outcome, "~ Age_at_study_date + Sex +", ppn_var, "+", nbm)), data = df)
  
  print_model_block(m_base, m_A1, m_A2, c(
    "Step 1: Age + Sex",
    paste0("Step 2: + PPN (", ppn_var, ")"),
    paste0("Step 3: + NBM (", nbm, ")")
  ))
  
  # ── ORDER B: NBM first, then PPN ──────────────────────────────────────────
  cat(" ORDER B: Age+Sex → +NBM → +PPN\n")
  cat(" (Does PPN explain variance OVER AND ABOVE NBM?)\n\n")
  
  m_B1 <- lm(as.formula(paste(outcome, "~ Age_at_study_date + Sex +", nbm)), data = df)
  m_B2 <- lm(as.formula(paste(outcome, "~ Age_at_study_date + Sex +", nbm, "+", ppn_var)), data = df)
  
  print_model_block(m_base, m_B1, m_B2, c(
    "Step 1: Age + Sex",
    paste0("Step 2: + NBM (", nbm, ")"),
    paste0("Step 3: + PPN (", ppn_var, ")")
  ))
  
  # ── Partial correlations ───────────────────────────────────────────────────
  resid_out <- resid(m_base)
  resid_ppn <- resid(lm(as.formula(paste(ppn_var, "~ Age_at_study_date + Sex")), data = df))
  resid_nbm <- resid(lm(as.formula(paste(nbm,     "~ Age_at_study_date + Sex")), data = df))
  
  pr_ppn <- cor(resid_out, resid_ppn, use = "complete.obs")
  pr_nbm <- cor(resid_out, resid_nbm, use = "complete.obs")
  
  cat(sprintf(" Partial r (PPN ~ slowing, adj age+sex): %.3f\n", pr_ppn))
  cat(sprintf(" Partial r (NBM ~ slowing, adj age+sex): %.3f\n\n", pr_nbm))
  
  # ── Save summary row ───────────────────────────────────────────────────────
  s_A2 <- summary(m_A2)
  s_B2 <- summary(m_B2)
  
  ppn_beta_joint <- s_A2$coefficients[ppn_var, "Estimate"]
  ppn_p_joint    <- s_A2$coefficients[ppn_var, "Pr(>|t|)"]
  nbm_beta_joint <- s_A2$coefficients[nbm,     "Estimate"]
  nbm_p_joint    <- s_A2$coefficients[nbm,     "Pr(>|t|)"]
  
  deltaR2_ppn_over_nbm <- s_B2$r.squared - summary(m_B1)$r.squared
  deltaR2_nbm_over_ppn <- s_A2$r.squared - summary(m_A1)$r.squared
  
  results_summary <- rbind(results_summary, data.frame(
    NBM_Network        = nbm,
    n                  = n,
    Partial_r_PPN      = round(pr_ppn, 3),
    Partial_r_NBM      = round(pr_nbm, 3),
    PPN_Beta_joint     = round(ppn_beta_joint, 4),
    PPN_p_joint        = round(ppn_p_joint,    4),
    NBM_Beta_joint     = round(nbm_beta_joint, 4),
    NBM_p_joint        = round(nbm_p_joint,    4),
    DeltaR2_PPN_overNBM = round(deltaR2_ppn_over_nbm, 3),
    DeltaR2_NBM_overPPN = round(deltaR2_nbm_over_ppn, 3),
    R2_full            = round(s_A2$r.squared, 3),
    stringsAsFactors   = FALSE
  ))
}

# ── Final summary table ────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════\n")
cat(" SUMMARY: Joint model results across NBM networks\n")
cat(" Outcome: Occipital REM EEG Slowing Ratio\n")
cat(" Both PPN and NBM in same model, adjusted for age + sex\n")
cat("══════════════════════════════════════════════════════════\n")
print(results_summary, row.names = FALSE)