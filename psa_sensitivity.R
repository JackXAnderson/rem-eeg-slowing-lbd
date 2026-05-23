library(dplyr)

# Load and clean data
data <- micro_analysis_with_MRI_200225_final

# Filter for Micro_YN == 1 and groups 1, 3, 4
data <- data %>%
  filter(Micro_YN == 1, Diagnosis_number %in% c(1, 3, 4))

# ── Zero-imputation for medication/dose variables ──────────────────────────────
# Controls will not be prescribed dopaminergic agents, ChEIs, or antidepressants,
# so NA in these columns genuinely means "not taking" (= 0 / No).
# We impute before converting to factor so 0 maps cleanly to the reference level.

data$DDE <- ifelse(is.na(data$DDE), 0, data$DDE)

data$Antidepressant_YN <- ifelse(is.na(data$Antidepressant_YN), 0,
                                 data$Antidepressant_YN)
data$Cholinesterase_YN <- ifelse(is.na(data$Cholinesterase_YN), 0,
                                 data$Cholinesterase_YN)

# RDI_in_REM: only impute if missing — true missingness (no REM recorded)
# is handled by complete-case filtering per variable below, so we leave it as-is.
# If you know a missing RDI_in_REM also means 0 events, uncomment the next line:
# data$RDI_in_REM <- ifelse(is.na(data$RDI_in_REM), 0, data$RDI_in_REM)

# Convert categorical variables to factors
data$Diagnosis_number  <- as.factor(data$Diagnosis_number)
data$Sex               <- as.factor(data$Sex)
data$Antidepressant_YN <- as.factor(data$Antidepressant_YN)
data$Cholinesterase_YN <- as.factor(data$Cholinesterase_YN)

# ── Variable list ──────────────────────────────────────────────────────────────
demographic_variables <- c(
  "Central_Absolute_REM_SFratio1",
  "Frontal_Absolute_REM_SFratio1",
  "Occipital_Absolute_REM_SFratio1"
)

# ── Covariates ─────────────────────────────────────────────────────────────────
# Original: Age + Sex
# New sensitivity: + RDI_in_REM + Antidepressant_YN + Cholinesterase_YN + DDE
base_covariates        <- c("Age_at_study_date", "Sex")
sensitivity_covariates <- c("Age_at_study_date", "Sex",
                            "RDI_in_REM", "Antidepressant_YN",
                            "Cholinesterase_YN", "DDE")

# Covariate terms used in formula strings
base_cov_str        <- paste(base_covariates,        collapse = " + ")
sensitivity_cov_str <- paste(sensitivity_covariates, collapse = " + ")

# ── Helper: F-statistic for Diagnosis term via permutation ────────────────────
compute_f_stat <- function(df, measure, cov_str) {
  formula_str <- paste0("`", measure, "` ~ Diagnosis_number + ", cov_str)
  model <- glm(as.formula(formula_str), data = df, family = gaussian())
  aov   <- anova(model)
  if ("Diagnosis_number" %in% rownames(aov)) aov["Diagnosis_number", "F"] else NA
}

permutation_p <- function(df, measure, cov_str, n_perm = 1000) {
  obs_f <- compute_f_stat(df, measure, cov_str)
  if (is.na(obs_f)) return(NA)
  perm_f <- replicate(n_perm, {
    df_perm <- df
    df_perm$Diagnosis_number <- sample(df_perm$Diagnosis_number)
    compute_f_stat(df_perm, measure, cov_str)
  })
  mean(perm_f >= obs_f, na.rm = TRUE)
}

# ── Helper: covariate beta + p from full sensitivity model ───────────────────
extract_covariate_effects <- function(df, measure, cov_str) {
  formula_str <- paste0("`", measure, "` ~ Diagnosis_number + ", cov_str)
  model <- lm(as.formula(formula_str), data = df)
  coef_tbl <- summary(model)$coefficients
  # Grab rows for our extra covariates only
  extra_covs <- c("RDI_in_REM", "Antidepressant_YN1", "Cholinesterase_YN1", "DDE")
  out <- list()
  for (cov in extra_covs) {
    if (cov %in% rownames(coef_tbl)) {
      out[[cov]] <- c(Beta = coef_tbl[cov, "Estimate"],
                      P    = coef_tbl[cov, "Pr(>|t|)"])
    } else {
      out[[cov]] <- c(Beta = NA, P = NA)
    }
  }
  out
}

# ── Helper: residual post-hoc permutation ─────────────────────────────────────
residual_posthoc <- function(df, measure, cov_str, g1, g2, n_perm = 10000) {
  resid_formula <- paste0("`", measure, "` ~ ", cov_str)
  resid_model   <- lm(as.formula(resid_formula), data = df)
  residuals     <- resid(resid_model)
  groups        <- droplevels(df$Diagnosis_number)
  
  idx <- which(groups %in% c(g1, g2))
  if (length(unique(groups[idx])) < 2 || length(idx) < 6) return(NA)
  grp <- droplevels(factor(groups[idx], levels = c(g1, g2)))
  
  obs_diff <- abs(mean(residuals[idx][grp == g1], na.rm = TRUE) -
                    mean(residuals[idx][grp == g2], na.rm = TRUE))
  perm_diffs <- replicate(n_perm, {
    shuf <- sample(grp)
    abs(mean(residuals[idx][shuf == g1], na.rm = TRUE) -
          mean(residuals[idx][shuf == g2], na.rm = TRUE))
  })
  mean(perm_diffs >= obs_diff, na.rm = TRUE)
}

# No FDR correction — raw permutation p-values reported throughout

# ── Main loop ─────────────────────────────────────────────────────────────────
glm_results     <- data.frame()
posthoc_results <- data.frame()
covariate_effects_list <- list()

pairs <- list(c("1","3"), c("1","4"), c("3","4"))

for (measure in demographic_variables) {
  
  # ── Data prep ──
  if (!measure %in% names(data) || all(is.na(data[[measure]]))) {
    warning(paste("Skipping:", measure, "- missing")); next
  }
  data[[measure]] <- as.numeric(data[[measure]])
  
  # Complete cases: medication/dose NAs already imputed to 0 above.
  # Only drop rows where the EEG measure itself or RDI_in_REM is missing.
  temp_sens <- data %>%
    filter(!is.na(.data[[measure]]),
           !is.na(RDI_in_REM)) %>%
    droplevels()
  
  if (length(unique(temp_sens$Diagnosis_number)) < 2) {
    warning(paste("Skipping:", measure, "- fewer than 2 groups")); next
  }
  
  # ── Original model p (base covariates, same sample as sensitivity for fairness) ──
  p_orig <- permutation_p(temp_sens, measure, base_cov_str)
  
  # ── Sensitivity model p ──
  p_sens <- permutation_p(temp_sens, measure, sensitivity_cov_str)
  
  # ── GLM results ──
  glm_results <- rbind(glm_results, data.frame(
    Variable    = measure,
    P_Original  = p_orig,
    P_Sensitive = p_sens,
    stringsAsFactors = FALSE
  ))
  
  # ── Covariate effects from sensitivity model ──
  cov_fx <- extract_covariate_effects(temp_sens, measure, sensitivity_cov_str)
  covariate_effects_list[[measure]] <- data.frame(
    Variable          = measure,
    RDI_in_REM_Beta   = cov_fx[["RDI_in_REM"]][["Beta"]],
    RDI_in_REM_P      = cov_fx[["RDI_in_REM"]][["P"]],
    Antidep_Beta      = cov_fx[["Antidepressant_YN1"]][["Beta"]],
    Antidep_P         = cov_fx[["Antidepressant_YN1"]][["P"]],
    ChEI_Beta         = cov_fx[["Cholinesterase_YN1"]][["Beta"]],
    ChEI_P            = cov_fx[["Cholinesterase_YN1"]][["P"]],
    DDE_Beta          = cov_fx[["DDE"]][["Beta"]],
    DDE_P             = cov_fx[["DDE"]][["P"]],
    stringsAsFactors  = FALSE
  )
  
  # ── Post-hoc (sensitivity model residuals) ──
  ph_row <- data.frame(Variable = measure,
                       `1_vs_3` = NA, `1_vs_4` = NA, `3_vs_4` = NA,
                       stringsAsFactors = FALSE)
  for (pair in pairs) {
    col <- paste0(pair[1], "_vs_", pair[2])
    ph_row[[col]] <- residual_posthoc(temp_sens, measure, sensitivity_cov_str,
                                      pair[1], pair[2])
  }
  posthoc_results <- rbind(posthoc_results, ph_row)
}

# ── Format numeric columns ────────────────────────────────────────────────────
fmt <- function(x) formatC(suppressWarnings(as.numeric(x)), format = "f", digits = 4)

covariate_effects <- do.call(rbind, covariate_effects_list)

glm_results$P_Original  <- fmt(glm_results$P_Original)
glm_results$P_Sensitive <- fmt(glm_results$P_Sensitive)

posthoc_results[, c("1_vs_3","1_vs_4","3_vs_4")] <-
  lapply(posthoc_results[, c("1_vs_3","1_vs_4","3_vs_4")], fmt)

for (col in c("RDI_in_REM_Beta","RDI_in_REM_P","Antidep_Beta","Antidep_P",
              "ChEI_Beta","ChEI_P","DDE_Beta","DDE_P")) {
  covariate_effects[[col]] <- fmt(covariate_effects[[col]])
}

# ── Print results ─────────────────────────────────────────────────────────────

cat("\n========================================================\n")
cat(" TABLE 1: Group effect — Original vs Sensitivity model\n")
cat(" REM bands only | Groups: 1=HC  3=PD  4=DLB\n")
cat(" Sensitivity covariates: Age, Sex, RDI_in_REM,\n")
cat("   Antidepressant_YN, Cholinesterase_YN, DDE\n")
cat(" Raw permutation p-values (no FDR correction)\n")
cat("========================================================\n")
print(glm_results, row.names = FALSE)

cat("\n========================================================\n")
cat(" TABLE 2: Post-hoc pairwise p-values (sensitivity model)\n")
cat(" Residual permutation | Raw p-values\n")
cat(" 1=HC  3=PD  4=DLB\n")
cat("========================================================\n")
print(posthoc_results, row.names = FALSE)

cat("\n========================================================\n")
cat(" TABLE 3: Covariate effects within sensitivity model\n")
cat(" Beta = unstandardised | Raw p-values\n")
cat("========================================================\n")
print(covariate_effects, row.names = FALSE)