library(dplyr)

data <- micro_analysis_with_MRI_200225_final

# ── DLB only, with ChEI status ─────────────────────────────────────────────────
dlb <- data %>%
  filter(Diagnosis_number == 4) %>%
  mutate(
    ChEI       = factor(Cholinesterase_YN, levels = c(0, 1),
                        labels = c("No ChEI", "ChEI")),
    occipital  = as.numeric(Occipital_Absolute_REM_SFratio1),
    mmse       = as.numeric(MMSE_result_macro),
    updrs3     = as.numeric(MDUPDRS_Section_III_Total_macro),
    rbdsq      = as.numeric(REMSBD_Total),
    age        = as.numeric(Age_at_study_date),
    sex        = as.factor(Sex)
  )

cat("══════════════════════════════════════════════\n")
cat(" DLB SAMPLE — ChEI STATUS\n")
cat("══════════════════════════════════════════════\n")
print(table(dlb$ChEI, useNA = "ifany"))
cat("\n")

# ── Descriptives by ChEI group ─────────────────────────────────────────────────
cat("── Group descriptives ─────────────────────────────────────\n")
desc <- dlb %>%
  group_by(ChEI) %>%
  summarise(
    n             = n(),
    Age_mean      = round(mean(age,       na.rm = TRUE), 1),
    Age_sd        = round(sd(age,         na.rm = TRUE), 1),
    MMSE_mean     = round(mean(mmse,      na.rm = TRUE), 1),
    MMSE_sd       = round(sd(mmse,        na.rm = TRUE), 1),
    UPDRS3_mean   = round(mean(updrs3,    na.rm = TRUE), 1),
    UPDRS3_sd     = round(sd(updrs3,      na.rm = TRUE), 1),
    RBDSQ_mean    = round(mean(rbdsq,     na.rm = TRUE), 1),
    RBDSQ_sd      = round(sd(rbdsq,       na.rm = TRUE), 1),
    Occipital_mean= round(mean(occipital, na.rm = TRUE), 3),
    Occipital_sd  = round(sd(occipital,   na.rm = TRUE), 3),
    .groups = "drop"
  )
print(as.data.frame(desc), row.names = FALSE)
cat("\n")

# ── Helper: clean regression output ───────────────────────────────────────────
run_model <- function(df, outcome_var, label) {
  
  df_clean <- df %>%
    filter(
      !is.na(.data[[outcome_var]]),
      !is.na(ChEI),
      !is.na(mmse),
      !is.na(updrs3),
      !is.na(age),
      !is.na(sex)
    ) %>%
    mutate(y = as.numeric(.data[[outcome_var]]),
           y_z = as.numeric(scale(y))) %>%
    filter(abs(y_z) <= 3)
  
  n_chei    <- sum(df_clean$ChEI == "ChEI",    na.rm = TRUE)
  n_nochei  <- sum(df_clean$ChEI == "No ChEI", na.rm = TRUE)
  
  cat(sprintf("── %s ──────────────────────────────────────\n", label))
  cat(sprintf("   n (ChEI): %d  |  n (No ChEI): %d\n\n", n_chei, n_nochei))
  
  if (n_chei < 3 || n_nochei < 3) {
    cat("   Insufficient n in one group — skipping.\n\n")
    return(invisible(NULL))
  }
  
  # Unadjusted
  m_unadj <- lm(y ~ ChEI, data = df_clean)
  s_unadj  <- summary(m_unadj)
  
  # Adjusted for age, sex, MMSE, UPDRS III
  m_adj   <- lm(y ~ ChEI + age + sex + mmse + updrs3, data = df_clean)
  s_adj   <- summary(m_adj)
  
  # Group means (raw)
  means <- df_clean %>%
    group_by(ChEI) %>%
    summarise(mean = round(mean(y, na.rm = TRUE), 3),
              sd   = round(sd(y,   na.rm = TRUE), 3),
              .groups = "drop")
  
  cat("   Group means (raw):\n")
  print(as.data.frame(means), row.names = FALSE)
  cat("\n")
  
  cat("   Unadjusted model (ChEI effect):\n")
  coef_u <- s_unadj$coefficients
  cat(sprintf("     β = %+.4f  SE = %.4f  p = %.4f\n\n",
              coef_u["ChEIChEI", "Estimate"],
              coef_u["ChEIChEI", "Std. Error"],
              coef_u["ChEIChEI", "Pr(>|t|)"]))
  
  cat("   Adjusted model (ChEI + age + sex + MMSE + UPDRS III):\n")
  coef_a <- s_adj$coefficients
  if ("ChEIChEI" %in% rownames(coef_a)) {
    cat(sprintf("     ChEI:    β = %+.4f  SE = %.4f  p = %.4f\n",
                coef_a["ChEIChEI",  "Estimate"],
                coef_a["ChEIChEI",  "Std. Error"],
                coef_a["ChEIChEI",  "Pr(>|t|)"]))
  }
  if ("mmse" %in% rownames(coef_a)) {
    cat(sprintf("     MMSE:    β = %+.4f  SE = %.4f  p = %.4f\n",
                coef_a["mmse",    "Estimate"],
                coef_a["mmse",    "Std. Error"],
                coef_a["mmse",    "Pr(>|t|)"]))
  }
  if ("updrs3" %in% rownames(coef_a)) {
    cat(sprintf("     UPDRS3:  β = %+.4f  SE = %.4f  p = %.4f\n",
                coef_a["updrs3",  "Estimate"],
                coef_a["updrs3",  "Std. Error"],
                coef_a["updrs3",  "Pr(>|t|)"]))
  }
  cat(sprintf("     R² (full model): %.3f\n\n", s_adj$r.squared))
  
  invisible(list(m_unadj = m_unadj, m_adj = m_adj))
}

# ══════════════════════════════════════════════════════════
# ANALYSIS 1: Occipital REM Slowing Ratio
# ══════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════════════\n")
cat(" ANALYSIS 1: Occipital REM Slowing Ratio\n")
cat("══════════════════════════════════════════════\n\n")

run_model(dlb, "occipital", "Occipital REM SFratio1")

# ══════════════════════════════════════════════════════════
# ANALYSIS 2: RBDSQ Total
# ══════════════════════════════════════════════════════════
cat("══════════════════════════════════════════════\n")
cat(" ANALYSIS 2: RBD Screening Questionnaire\n")
cat("══════════════════════════════════════════════\n\n")

run_model(dlb, "rbdsq", "RBDSQ Total")

# ══════════════════════════════════════════════════════════
# ANALYSIS 3: MMSE as outcome (does ChEI associate with
# cognition independently — sanity check)
# ══════════════════════════════════════════════════════════
cat("══════════════════════════════════════════════\n")
cat(" ANALYSIS 3: MMSE (sanity check)\n")
cat(" Expected: ChEI group more cognitively impaired\n")
cat("══════════════════════════════════════════════\n\n")

dlb_mmse <- dlb %>%
  filter(!is.na(mmse), !is.na(ChEI), !is.na(age), !is.na(sex)) %>%
  mutate(y = mmse)

m_mmse <- lm(y ~ ChEI + age + sex, data = dlb_mmse)
s_mmse <- summary(m_mmse)
coef_m <- s_mmse$coefficients

cat(sprintf("   n (ChEI): %d  |  n (No ChEI): %d\n\n",
            sum(dlb_mmse$ChEI == "ChEI"),
            sum(dlb_mmse$ChEI == "No ChEI")))
cat("   MMSE means by group:\n")
dlb_mmse %>% group_by(ChEI) %>%
  summarise(mean = round(mean(y, na.rm=TRUE), 1),
            sd   = round(sd(y,   na.rm=TRUE), 1),
            .groups="drop") %>%
  as.data.frame() %>% print(row.names=FALSE)
cat("\n")
if ("ChEIChEI" %in% rownames(coef_m)) {
  cat(sprintf("   ChEI effect on MMSE (adj age+sex):\n"))
  cat(sprintf("     β = %+.3f  p = %.4f\n\n",
              coef_m["ChEIChEI", "Estimate"],
              coef_m["ChEIChEI", "Pr(>|t|)"]))
}

cat("══════════════════════════════════════════════\n")
cat(" NOTE ON INTERPRETATION\n")
cat("══════════════════════════════════════════════\n")
cat(" If ChEI group has MORE slowing: consistent with\n")
cat("   ChEI prescribed to more impaired patients\n")
cat("   (confounding by indication)\n\n")
cat(" If ChEI group has LESS slowing: consistent with\n")
cat("   cholinergic restoration improving REM EEG\n")
cat("   (supports reviewer hypothesis)\n\n")
cat(" The adjusted model attempts to control for\n")
cat("   disease severity (MMSE + UPDRS III) to\n")
cat("   isolate the ChEI effect\n")
cat("══════════════════════════════════════════════\n")