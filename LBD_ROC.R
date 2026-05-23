library(dplyr)
library(pROC)

data <- micro_analysis_with_MRI_200225_final

# ── Prepare data ───────────────────────────────────────────────────────────────
# Keep HC (1), PD (3), DLB (4) with complete cases
clf_data <- data %>%
  filter(
    Diagnosis_number %in% c(1, 3, 4),
    !is.na(Occipital_Absolute_REM_SFratio1),
    !is.na(Age_at_study_date),
    !is.na(Sex)
  ) %>%
  mutate(
    ratio     = as.numeric(Occipital_Absolute_REM_SFratio1),
    ratio_z   = as.numeric(scale(ratio))
  ) %>%
  filter(abs(ratio_z) <= 3) %>%   # 3SD outlier removal
  mutate(
    Group = factor(Diagnosis_number,
                   levels = c(1, 3, 4),
                   labels = c("HC", "PD", "DLB"))
  )

cat("══════════════════════════════════════════════\n")
cat(" SAMPLE\n")
cat("══════════════════════════════════════════════\n")
print(table(clf_data$Group))
cat("\n")

# ── Helper: age+sex adjusted residuals for fair ROC ───────────────────────────
# We use age+sex adjusted residuals so the ROC reflects the
# biomarker's discriminative value independent of demographic differences
get_residuals <- function(df) {
  model <- lm(ratio ~ Age_at_study_date + Sex, data = df)
  resid(model)
}

# ── Helper: LOOCV AUC ─────────────────────────────────────────────────────────
loocv_auc <- function(df, predictor_col, response_col) {
  n        <- nrow(df)
  pred_out <- numeric(n)
  
  for (i in 1:n) {
    train <- df[-i, ]
    test  <- df[i,  ]
    
    # Fit logistic regression on training fold
    formula <- as.formula(paste(response_col, "~", predictor_col))
    m <- glm(formula, data = train, family = binomial())
    
    # Predict on left-out case
    pred_out[i] <- predict(m, newdata = test, type = "response")
  }
  
  roc_obj <- roc(df[[response_col]], pred_out,
                 quiet = TRUE, direction = "<")
  list(auc = as.numeric(auc(roc_obj)), roc = roc_obj, preds = pred_out)
}

# ── Helper: print clean ROC summary ──────────────────────────────────────────
print_roc <- function(roc_obj, label) {
  ci  <- ci.auc(roc_obj, conf.level = 0.95)
  cat(sprintf("  AUC = %.3f  (95%% CI: %.3f\u2013%.3f)\n",
              as.numeric(auc(roc_obj)), ci[1], ci[3]))
  
  # Best threshold by Youden index
  coords_best <- coords(roc_obj, "best", ret = c("threshold","sensitivity","specificity"),
                        best.method = "youden", transpose = FALSE)
  cat(sprintf("  Best threshold (Youden): %.3f\n",      coords_best$threshold))
  cat(sprintf("  Sensitivity: %.1f%%\n", coords_best$sensitivity * 100))
  cat(sprintf("  Specificity: %.1f%%\n", coords_best$specificity * 100))
  cat("\n")
}

# ══════════════════════════════════════════════════════════
# COMPARISON: DLB vs PD
# ══════════════════════════════════════════════════════════
cat("══════════════════════════════════════════════\n")
cat(" COMPARISON: DLB vs PD\n")
cat("══════════════════════════════════════════════\n\n")

dlb_pd <- clf_data %>%
  filter(Group %in% c("PD", "DLB")) %>%
  droplevels() %>%
  mutate(
    resid_ratio = get_residuals(.),
    binary      = as.integer(Group == "DLB")
  )

cat(sprintf(" n: PD=%d  DLB=%d\n\n",
            sum(dlb_pd$Group == "PD"),
            sum(dlb_pd$Group == "DLB")))

# Raw ratio ROC
roc_dlb_pd_raw <- roc(dlb_pd$binary, dlb_pd$ratio,
                      quiet = TRUE, direction = "<")
cat(" RAW ratio (unadjusted):\n")
print_roc(roc_dlb_pd_raw, "DLB vs PD raw")

# Age+sex adjusted residuals ROC
roc_dlb_pd_adj <- roc(dlb_pd$binary, dlb_pd$resid_ratio,
                      quiet = TRUE, direction = "<")
cat(" Age+sex ADJUSTED residuals:\n")
print_roc(roc_dlb_pd_adj, "DLB vs PD adjusted")

# LOOCV
cat(" LOOCV (age+sex adjusted, logistic regression):\n")
loocv_pd <- loocv_auc(dlb_pd, "resid_ratio", "binary")
cat(sprintf("  LOOCV AUC = %.3f\n\n", loocv_pd$auc))

# ══════════════════════════════════════════════════════════
# ROC PLOT
# ══════════════════════════════════════════════════════════
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

plot(roc_dlb_pd_raw,
     main = "DLB vs PD — Raw ratio",
     col  = "#C0392B", lwd = 2,
     print.auc = TRUE, print.auc.y = 0.15,
     auc.polygon = TRUE, auc.polygon.col = "#FADBD8")
abline(a = 0, b = 1, lty = 2, col = "grey60")

plot(roc_dlb_pd_adj,
     main = "DLB vs PD — Age+sex adjusted",
     col  = "#1F4E79", lwd = 2,
     print.auc = TRUE, print.auc.y = 0.15,
     auc.polygon = TRUE, auc.polygon.col = "#D6E4F0")
abline(a = 0, b = 1, lty = 2, col = "grey60")

# ══════════════════════════════════════════════════════════
# CI for adjusted AUC
# ══════════════════════════════════════════════════════════
ci_adj <- ci.auc(roc_dlb_pd_adj, conf.level = 0.95)
coords_best <- coords(roc_dlb_pd_adj, "best",
                      ret = c("threshold","sensitivity","specificity"),
                      best.method = "youden", transpose = FALSE)

# ══════════════════════════════════════════════════════════
# CLEAN SUMMARY
# ══════════════════════════════════════════════════════════
cat("══════════════════════════════════════════════\n")
cat(" SUMMARY: DLB vs PD Classification\n")
cat(" Occipital REM EEG Slowing Ratio\n")
cat("══════════════════════════════════════════════\n")
cat(sprintf(" n PD:                %d\n",   sum(dlb_pd$binary == 0)))
cat(sprintf(" n DLB:               %d\n",   sum(dlb_pd$binary == 1)))
cat(sprintf(" AUC (raw):           %.3f\n", as.numeric(auc(roc_dlb_pd_raw))))
cat(sprintf(" AUC (age+sex adj):   %.3f  (95%% CI: %.3f-%.3f)\n",
            as.numeric(auc(roc_dlb_pd_adj)), ci_adj[1], ci_adj[3]))
cat(sprintf(" AUC (LOOCV):         %.3f\n", loocv_pd$auc))
cat(sprintf(" Best threshold:      %.3f\n", coords_best$threshold))
cat(sprintf(" Sensitivity:         %.1f%%\n", coords_best$sensitivity * 100))
cat(sprintf(" Specificity:         %.1f%%\n", coords_best$specificity * 100))
cat("\nNote: AUC_adj uses age+sex residuals. LOOCV = leave-one-out cross-validation.\n")