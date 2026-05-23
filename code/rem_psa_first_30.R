library(dplyr)
library(tidyr)
library(openxlsx)

select <- dplyr::select
set.seed(123)

psd <- read.table(
  "/Users/jackanderson/Documents/luna/harmonized_edfs/final/out/rem_psd_bands.txt",
  header = TRUE, stringsAsFactors = FALSE
)

demo <- read.table(
  "/Users/jackanderson/Documents/luna/masterfile/master_test.txt",
  header = TRUE, stringsAsFactors = FALSE
)

cat("Bands available:\n")
print(unique(psd$B))

slowing <- psd %>%
  filter(B %in% c("DELTA", "THETA", "ALPHA", "SIGMA", "BETA")) %>%
  group_by(ID, CH, B) %>%
  summarise(PSD = mean(PSD, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = B, values_from = PSD) %>%
  mutate(SLOWING_RATIO = (DELTA + THETA) / (ALPHA + SIGMA + BETA)) %>%
  dplyr::select(ID, CH, SLOWING_RATIO)

psd_demo <- merge(psd, demo[, c("ID", "Diagnosis_number", "Age_at_study_date")],
                  by = "ID", all.x = TRUE)
psd_demo <- merge(psd_demo, slowing, by = c("ID", "CH"), all.x = TRUE)

psd_demo <- psd_demo %>%
  filter(Diagnosis_number %in% c(1, 3, 4)) %>%
  filter(!ID %in% c("PD_00567_DH", "PD_00672_IR"))

psd_demo$Diagnosis_number <- factor(psd_demo$Diagnosis_number,
                                    levels = c(1, 3, 4), labels = c("HC", "PD", "DLB"))
psd_demo$CH <- factor(psd_demo$CH, levels = c("F3", "C3", "O1"))

compute_f <- function(data, outcome = "PSD") {
  formula <- as.formula(paste(outcome, "~ Age_at_study_date + Diagnosis_number"))
  model <- lm(formula, data = data)
  anova_res <- anova(model)
  if ("Diagnosis_number" %in% rownames(anova_res))
    return(anova_res["Diagnosis_number", "F value"])
  return(NA)
}

residual_posthoc <- function(residuals, groups, g1, g2, nperm = 5000) {
  idx <- which(groups %in% c(g1, g2))
  group <- droplevels(groups[idx])
  obs <- abs(mean(residuals[idx][group == g1]) - mean(residuals[idx][group == g2]))
  perm <- replicate(nperm, {
    shuffled <- sample(group)
    abs(mean(residuals[idx][shuffled == g1]) - mean(residuals[idx][shuffled == g2]))
  })
  mean(perm >= obs)
}

results   <- data.frame()
channels  <- levels(psd_demo$CH)
bands     <- unique(psd_demo$B)

# ── Band-level PSD loop ───────────────────────────────────────────────────────
for (ch in channels) {
  for (b in bands) {
    
    df <- psd_demo %>% filter(CH == ch, B == b)
    if (nrow(df) < 10) next
    
    # ── MEAN AND SD per group ──
    stats <- df %>%
      group_by(Diagnosis_number) %>%
      summarise(m = mean(PSD, na.rm = TRUE),
                s = sd(PSD,   na.rm = TRUE),
                .groups = "drop")
    
    means <- stats %>% dplyr::select(Diagnosis_number, m) %>%
      pivot_wider(names_from = Diagnosis_number, values_from = m)
    sds   <- stats %>% dplyr::select(Diagnosis_number, s) %>%
      pivot_wider(names_from = Diagnosis_number, values_from = s)
    
    obsF  <- compute_f(df, outcome = "PSD")
    nperm <- 2000
    permF <- numeric(nperm)
    for (i in 1:nperm) {
      tmp <- df; tmp$Diagnosis_number <- sample(tmp$Diagnosis_number)
      permF[i] <- compute_f(tmp, outcome = "PSD")
    }
    p_perm <- mean(permF >= obsF, na.rm = TRUE)
    
    lm_res    <- lm(PSD ~ Age_at_study_date, data = df)
    residuals <- resid(lm_res)
    groups    <- df$Diagnosis_number
    
    results <- rbind(results, data.frame(
      Channel   = ch,
      Band      = b,
      HC_Mean   = means$HC,   HC_SD   = sds$HC,
      PD_Mean   = means$PD,   PD_SD   = sds$PD,
      DLB_Mean  = means$DLB,  DLB_SD  = sds$DLB,
      F_stat    = obsF,
      P_perm    = p_perm,
      HC_vs_PD  = residual_posthoc(residuals, groups, "HC", "PD"),
      HC_vs_DLB = residual_posthoc(residuals, groups, "HC", "DLB"),
      PD_vs_DLB = residual_posthoc(residuals, groups, "PD", "DLB")
    ))
  }
}

# ── Slowing ratio loop ────────────────────────────────────────────────────────
slowing_results <- data.frame()

for (ch in channels) {
  
  df <- psd_demo %>%
    filter(CH == ch) %>%
    distinct(ID, CH, Diagnosis_number, Age_at_study_date, SLOWING_RATIO) %>%
    filter(!is.na(SLOWING_RATIO))
  
  if (nrow(df) < 10) next
  
  # ── MEAN AND SD per group ──
  stats <- df %>%
    group_by(Diagnosis_number) %>%
    summarise(m = mean(SLOWING_RATIO, na.rm = TRUE),
              s = sd(SLOWING_RATIO,   na.rm = TRUE),
              .groups = "drop")
  
  means <- stats %>% dplyr::select(Diagnosis_number, m) %>%
    pivot_wider(names_from = Diagnosis_number, values_from = m)
  sds   <- stats %>% dplyr::select(Diagnosis_number, s) %>%
    pivot_wider(names_from = Diagnosis_number, values_from = s)
  
  obsF <- tryCatch({
    model <- lm(SLOWING_RATIO ~ Age_at_study_date + Diagnosis_number, data = df)
    anova(model)["Diagnosis_number", "F value"]
  }, error = function(e) NA)
  
  nperm <- 2000
  permF <- numeric(nperm)
  for (i in 1:nperm) {
    tmp <- df; tmp$Diagnosis_number <- sample(tmp$Diagnosis_number)
    permF[i] <- tryCatch({
      model <- lm(SLOWING_RATIO ~ Age_at_study_date + Diagnosis_number, data = tmp)
      anova(model)["Diagnosis_number", "F value"]
    }, error = function(e) NA)
  }
  p_perm <- mean(permF >= obsF, na.rm = TRUE)
  
  lm_res    <- lm(SLOWING_RATIO ~ Age_at_study_date, data = df)
  residuals <- resid(lm_res)
  groups    <- df$Diagnosis_number
  
  slowing_results <- rbind(slowing_results, data.frame(
    Channel   = ch,
    Band      = "SLOWING_RATIO",
    HC_Mean   = means$HC,   HC_SD   = sds$HC,
    PD_Mean   = means$PD,   PD_SD   = sds$PD,
    DLB_Mean  = means$DLB,  DLB_SD  = sds$DLB,
    F_stat    = obsF,
    P_perm    = p_perm,
    HC_vs_PD  = residual_posthoc(residuals, groups, "HC", "PD"),
    HC_vs_DLB = residual_posthoc(residuals, groups, "HC", "DLB"),
    PD_vs_DLB = residual_posthoc(residuals, groups, "PD", "DLB")
  ))
}

# ── Combine and save ──────────────────────────────────────────────────────────
results_combined <- bind_rows(results, slowing_results) %>%
  arrange(Channel, Band)

print(results_combined)

write.csv(results_combined,
          "/Users/jackanderson/Documents/luna/harmonized_edfs/final/out/REM_PSD_perm_results.csv",
          row.names = FALSE)
write.xlsx(results_combined,
           "/Users/jackanderson/Documents/luna/harmonized_edfs/final/out/REM_PSD_perm_results.xlsx")

cat("\nDONE\n")
