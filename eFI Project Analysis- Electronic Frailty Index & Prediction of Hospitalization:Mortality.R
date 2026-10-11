# ******************************************************************************
# eFI Project Analysis
# Electronic Frailty Index & Prediction of Hospitalization/Mortality
# Xena Eriksen
# June - July 2026
# Reference: Pajewski et al. (2019)
# ******************************************************************************


# ******************************************************************************
# TABLE OF CONTENTS ----
# ******************************************************************************
# 0. Setup & Packages
# 1. Data Import  
# 2. Data Cleaning & Type Conversion
# 3. Missingness / Selection Bias Check
# 4. Cohort Summary
# 5. Table: Baseline Characteristics
# 6. Table: Clinical Comorbidity Profile
# 7. Table: Stratified by Frailty Category
# 8. Regression Analysis Setup
# 9. Firth's Logistic Regression
# 10. Table: Full Results (Plain)
# 11. Table: AUC Summary (No Color)
# 12. Color Palette for Heat Maps
# 13. Table: Color-Coded AUC Heat Map
# 14. Table: Color-Coded Full Table
# 15. Table: Odds Ratios (with AUC)
# 16. Table: Odds Ratios (with AUC) — Color-Coded
# 17. Table: Full Tables Color-Coded with Significance Markers (Asterisks/Italic/Bold) — Color-Coded, Sorted by Outcome
# 18. Risk Gradient Visualization + Final Summary
# ******************************************************************************
 

# **********************************************
# --- 0. SETUP & PACKAGES ----
# **********************************************


## -- 0.2 install packages----
# install.packages(c("tidyverse", "readxl", "rms", "pROC", "gt", "logistf"))

## -- 0.3 load libraries ----
library(tidyverse)
library(readxl)
library(rms)
library(pROC)
library(gt)
library(logistf)


# **********************************************
# --- 1.0 IMPORT DATA ----
# **********************************************

## -- 1.1 Load data ----
efidata <- read_excel("~/Downloads/combined efi ATS Clean.xlsx")
View(efidata)

## -- 1.2 Quick look----
glimpse(efidata)
head(efidata)
summary(efidata)


# **********************************************
# --- 2.0 DATA CLEANING & TYPE CONVERSION ----
# **********************************************

## -- 2.1 fix outcome variables ----
# Excel imported these as "Yes"/"No" strings — convert to 1/0
outcome_vars <- c(
  "HospitalizedWithin6Months",
  "HospitalizedWithin12Months",
  "MortalityWithin6Months",
  "MortalityWithin12Months"
)

efidata <- efidata %>%
  mutate(across(all_of(outcome_vars), ~ ifelse(.x == "Yes", 1, 0)))

## -- 2.2 fix EFI score variables ----
# character -> numeric
efidata$efi_score <- as.numeric(efidata$efi_score)
efidata$efi_percentile <- as.numeric(efidata$efi_percentile)

## -- 2.3 fix other variables ----
# character -> numeric
char_scores_to_fix <- c(
  "smokingscore", "obesescore", "uwscore", "sbpscore", "dbpscore",
  "healthstatusscore", "adlscore", "iadlscore", "awvhearingscore",
  "cogscore", "awvstressscore", "awvpainscore", "polypharmacyscore",
  "smoking_tob_use_c"
)

for (var in char_scores_to_fix) {
  efidata[[var]] <- as.numeric(efidata[[var]])
}

## -- 2.4 fix chair score ----
# logical -> numeric
efidata$chairscore <- as.numeric(efidata$chairscore)

## -- 2.5 create completeness flag ----
# This lets us compare included vs excluded patients
efidata <- efidata %>%
  mutate(efi_complete = !is.na(efi_score))


# **********************************************
# --- 3.0 MISSINGNESS / SELECTION BIAS CHECK ----
# **********************************************
# 6,651 patients had no calculable eFI — why? is this random?

## -- 3.1 explore age distributions ----
age_comparison <- efidata %>%
  group_by(efi_complete) %>%
  summarise(
    n = n(),
    age_mean = round(mean(age, na.rm = TRUE), 2),
    age_sd = round(sd(age, na.rm = TRUE), 2),
    age_min = min(age, na.rm = TRUE),
    age_max = max(age, na.rm = TRUE)
  )
print(age_comparison)

## -- 3.2 T-test for age difference ----
t_test_result <- t.test(age ~ efi_complete, data = efidata)
print(t_test_result)


# **********************************************
# --- 4.0 COHORT SUMMARY ----
# **********************************************

## -- 4.1 filter to complete cases for analysis ----
efi_complete_data <- efidata %>% filter(!is.na(efi_score))

## -- 4.2 summary statistics for complete cases ----
efi_summary <- efi_complete_data %>%
  summarise(
    Total_n = n(),
    age_mean = round(mean(age, na.rm = TRUE), 2),
    age_sd = round(sd(age, na.rm = TRUE), 2),
    efi_score_mean = round(mean(efi_score, na.rm = TRUE), 3),
    efi_score_sd = round(sd(efi_score, na.rm = TRUE), 3),
    efi_percentile_mean = round(mean(efi_percentile, na.rm = TRUE), 2),
    efi_percentile_median = round(median(efi_percentile, na.rm = TRUE), 2),
    hosp_6mo_n = sum(HospitalizedWithin6Months, na.rm = TRUE),
    hosp_6mo_pct = round(mean(HospitalizedWithin6Months, na.rm = TRUE) * 100, 2),
    hosp_12mo_n = sum(HospitalizedWithin12Months, na.rm = TRUE),
    hosp_12mo_pct = round(mean(HospitalizedWithin12Months, na.rm = TRUE) * 100, 2),
    mort_6mo_n = sum(MortalityWithin6Months, na.rm = TRUE),
    mort_6mo_pct = round(mean(MortalityWithin6Months, na.rm = TRUE) * 100, 2),
    mort_12mo_n = sum(MortalityWithin12Months, na.rm = TRUE),
    mort_12mo_pct = round(mean(MortalityWithin12Months, na.rm = TRUE) * 100, 2)
  )
print(efi_summary)

## -- 4.3 breakdown by frailty category ----
efi_by_category <- efi_complete_data %>%
  group_by(efi_category) %>%
  summarise(
    n = n(),
    pct = round(n / sum(n) * 100, 2),
    age_mean = round(mean(age, na.rm = TRUE), 1),
    efi_score_mean = round(mean(efi_score, na.rm = TRUE), 2),
    hosp_6mo_rate = round(mean(HospitalizedWithin6Months, na.rm = TRUE) * 100, 2),
    hosp_12mo_rate = round(mean(HospitalizedWithin12Months, na.rm = TRUE) * 100, 2),
    mort_6mo_rate = round(mean(MortalityWithin6Months, na.rm = TRUE) * 100, 2),
    mort_12mo_rate = round(mean(MortalityWithin12Months, na.rm = TRUE) * 100, 2)
  ) %>%
  arrange(desc(n))
print(efi_by_category)


# **********************************************
# --- 5.0 TABLE: BASELINE CHARACTERISTICS ----
# **********************************************

## -- 5.1 format helpers - continuous and binary ----
### - 5.1.1 continuous variable formatter ----
fmt_cont <- function(x) {
  paste0(round(mean(x, na.rm = TRUE), 1), " ± ", round(sd(x, na.rm = TRUE), 1),
         ", Median: ", round(median(x, na.rm = TRUE), 1),
         " (", min(x, na.rm = TRUE), "–", max(x, na.rm = TRUE), ")")
}

### - 5.1.2 binary variable formatter ----
fmt_bin <- function(x) {
  n_pos <- sum(x == 1, na.rm = TRUE)
  paste0(n_pos, " (", round(n_pos / sum(!is.na(x)) * 100, 1), "%)")
}

## -- 5.2 variable definitions ----
cont_vars <- list(
  "Age (years)" = "age",
  "EFI Score" = "efi_score",
  "EFI Percentile" = "efi_percentile"
)

bin_outcomes <- list(
  "Hospitalized Within 6 Months" = "HospitalizedWithin6Months",
  "Hospitalized Within 12 Months" = "HospitalizedWithin12Months",
  "Mortality Within 6 Months" = "MortalityWithin6Months",
  "Mortality Within 12 Months" = "MortalityWithin12Months"
)

## -- 5.3 build table ----
efi_table <- bind_rows(
  imap_dfr(cont_vars, ~ tibble(Variable = .y,
                               n = sum(!is.na(efi_complete_data[[.x]])),
                               Summary = fmt_cont(efi_complete_data[[.x]]))),
  efi_complete_data %>%
    count(efi_category) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    transmute(Variable = paste0("  ", efi_category), n = n,
              Summary = paste0(n, " (", pct, "%)")),
  imap_dfr(bin_outcomes, ~ tibble(Variable = .y,
                                  n = sum(efi_complete_data[[.x]] == 1, na.rm = TRUE),
                                  Summary = fmt_bin(efi_complete_data[[.x]])))
)

### -- EFI table ----
cat("\n========== TABLE: EFI Table -- Core Variables ==========\n")
print(efi_table)

## -- 5.4 create and save the gt_tbl objects ----
efi_table_gt <- efi_table %>%
  gt() %>%
  tab_header(title = "Baseline Characteristics") %>%
  tab_options(table.font.size = 11)


# **********************************************
# --- 6.0 TABLE: CLINICAL COMORBIDITY PROFILE ----
# **********************************************

## -- 6.1 clinical score definitions ----
# All 36 eFI components — report n (%) present
clinical_scores <- c(
  "anemiascore", "rascore", "afscore", "strokescore", "renalscore",
  "diabetesscore", "dizzyscore", "dyspneascore", "fallsscore",
  "fracturescore", "hearingscore", "hfscore", "valvescore", "htnscore",
  "hyposcore", "amiscore", "cadscore", "melanomascore", "bloodcancerscore",
  "cancerscore", "demscore", "osteoscore", "pkscore", "pepticscore",
  "pvdscore", "suscore", "thyroidscore", "incontinencescore",
  "urinaryscore", "visionscore", "wlscore", "depressscore", "liverscore",
  "cpainscore", "deliriumscore"
)

## -- 6.2 display names for table ----
display_names <- c(
  "Anemia", "Rheumatoid Arthritis", "Atrial Fibrillation", "Stroke",
  "Renal Disease", "Diabetes", "Dizziness", "Dyspnea", "Falls",
  "Fracture", "Hearing Loss", "Heart Failure", "Valvular Disease",
  "Hypertension", "Hypothyroidism", "Acute MI", "Coronary Artery Disease",
  "Melanoma", "Blood Cancer", "Cancer (Other)", "Dementia", "Osteoporosis",
  "Parkinson's", "Peptic Ulcer", "Peripheral Vascular Disease", "Syncope",
  "Thyroid Disease", "Incontinence", "Urinary Issues", "Vision Loss",
  "Weight Loss", "Depression", "Liver Disease", "Chronic Pain", "Delirium"
)

## -- 6.3 build clinical comorbidity table ----
clinical_table <- map2_dfr(clinical_scores, display_names, function(score, label) {
  if (!(score %in% names(efidata))) return(NULL)
  col_data <- efidata[[score]]
  n_total <- sum(!is.na(col_data))
  unique_vals <- unique(col_data[!is.na(col_data)])
  
  if (all(unique_vals %in% c(0, 1))) {
    n_present <- sum(col_data == 1, na.rm = TRUE)
    tibble(Variable = label, n = n_total, Present_n = n_present,
           Summary = paste0(n_present, " (", round(n_present / n_total * 100, 1), "%)"))
  } else {
    tibble(Variable = label, n = n_total, Present_n = NA_integer_,
           Summary = paste0(round(mean(col_data, na.rm = TRUE), 3), " ± ",
                            round(sd(col_data, na.rm = TRUE), 3)))
  }
})

### -- clinical comorbidity table ----
cat("\n========== TABLE: Clinical Scores ==========\n")
print(clinical_table)

## -- 5.4 create and save the gt_tbl objects for clinical comorbidity table----
clinical_table_gt <- clinical_table %>%
  select(Variable, n, Summary) %>%
  gt() %>%
  tab_header(title = "Clinical Comorbidity Profile") %>%
  tab_options(table.font.size = 11)



# **********************************************
# --- 7.0 TABLE: STRATIFIED BY FRAILTY CATEGORY ----
# **********************************************

## --- 7.1 calc P-value ----
calc_p <- function(df, var_name, is_cont = FALSE) {
  ### -- 7.1.1 continuous variable: Kruskal-Wallis test ----
  if (is_cont) {
    kruskal.test(df[[var_name]] ~ df$efi_category)$p.value
  } 
  ### -- 7.1.2 binary variable: Chi-square test ----
  else {
    suppressWarnings(chisq.test(table(df$efi_category, df[[var_name]]))$p.value)
  }
}

### -- 7.1.3 format P-value ----
fmt_p <- function(p) {
  if (is.na(p)) return(NA)
  if (p < 0.001) return("<0.001")
  unname(formatC(round(p, 3), format = "f", digits = 3))
}

## -- 7.2 Row builder ----
get_row_stats <- function(data, var, is_cont = FALSE) {
  ### -- 7.2.1 continuous variables ----  
  if (is_cont) {
    fit_m  <- mean(data[[var]][data$efi_category == "Fit"], na.rm = TRUE)
    pre_m  <- mean(data[[var]][data$efi_category == "Pre-Frail"], na.rm = TRUE)
    fra_m  <- mean(data[[var]][data$efi_category == "Frail"], na.rm = TRUE)
    fit_sd <- sd(data[[var]][data$efi_category == "Fit"], na.rm = TRUE)
    pre_sd <- sd(data[[var]][data$efi_category == "Pre-Frail"], na.rm = TRUE)
    fra_sd <- sd(data[[var]][data$efi_category == "Frail"], na.rm = TRUE)
    
    p_val <- calc_p(data, var, is_cont = TRUE) 
    
    tibble(
      Type = "Continuous",
      Fit = paste0(round(fit_m, 1), " ± ", round(fit_sd, 1)),
      Pre_Frail = paste0(round(pre_m, 1), " ± ", round(pre_sd, 1)),
      Frail = paste0(round(fra_m, 1), " ± ", round(fra_sd, 1)),
      P_Value = fmt_p(p_val)
    )
  } 
  ### -- 7.2.2 binary variables ----
  else {
    dat <- data %>% filter(!is.na(.data[[var]]))
    
    n_f  <- sum(dat[[var]] == 1 & dat$efi_category == "Fit")
    n_pr <- sum(dat[[var]] == 1 & dat$efi_category == "Pre-Frail")
    n_fr <- sum(dat[[var]] == 1 & dat$efi_category == "Frail")
    
    p_f  <- round(n_f  / sum(dat$efi_category == "Fit") * 100, 1)
    p_pr <- round(n_pr / sum(dat$efi_category == "Pre-Frail") * 100, 1)
    p_fr <- round(n_fr / sum(dat$efi_category == "Frail") * 100, 1)
    
    p_val <- calc_p(data, var, is_cont = FALSE)
    
    tibble(
      Type = "Binary",
      Fit = paste0(n_f, " (", p_f, "%)"),
      Pre_Frail = paste0(n_pr, " (", p_pr, "%)"),
      Frail = paste0(n_fr, " (", p_fr, "%)"),
      P_Value = fmt_p(p_val)
    )
  }
}

## -- 7.3 Build stratified table ----
table_stratified <- bind_rows(
  get_row_stats(efi_complete_data, "age", is_cont = TRUE) %>% 
    mutate(Variable = "Age (years)"),
  get_row_stats(efi_complete_data, "HospitalizedWithin6Months")  %>% 
    mutate(Variable = "Hospitalized < 6mo"),
  get_row_stats(efi_complete_data, "HospitalizedWithin12Months") %>% 
    mutate(Variable = "Hospitalized < 12mo"),
  get_row_stats(efi_complete_data, "MortalityWithin6Months")     %>% 
    mutate(Variable = "Mortality < 6mo"),
  get_row_stats(efi_complete_data, "MortalityWithin12Months")    %>% 
    mutate(Variable = "Mortality < 12mo")
) %>%
  select(Variable, Fit, Pre_Frail, Frail, P_Value) %>%
  gt() %>%
  cols_label(
    Fit = "Fit (n=3,554)",
    Pre_Frail = "Pre-Frail (n=8,768)",
    Frail = "Frail (n=5,956)",
    P_Value = "P-Value"
  ) %>%
  tab_header(title = "Table: Outcomes Stratified by Frailty Status") %>%
  tab_source_note(md("Values are Mean ± SD or n (%). P-values from Kruskal-Wallis (continuous) or Chi-Square (binary).")) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(rows = 1)) %>%
  cols_align("left", columns = Variable) %>%
  cols_align("center", columns = -Variable) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )

### -- stratified table ----
cat("\n========== TABLE: Stratified (Table_Stratified.html) ==========\n")
gtsave(table_stratified, "Table_Stratified.html")


# **********************************************
# --- 8.0 REGRESSION ANALYSIS SETUP ----
# **********************************************

## -- 8.1 data prep -- shorthand outcome names ----
dat <- efidata %>%
  filter(!is.na(efi_score)) %>%
  mutate(
    hosp_6mo  = as.numeric(HospitalizedWithin6Months),
    hosp_12mo = as.numeric(HospitalizedWithin12Months),
    mort_6mo  = as.numeric(MortalityWithin6Months),
    mort_12mo = as.numeric(MortalityWithin12Months),
    efi_score = as.numeric(efi_score)
  ) %>%
  select(where(~ !all(is.na(.))))

## -- 8.2 remove columns that are all NA ----

## -- 8.3 distribution summaries for predictor variables for regression modeling ----
ddist <- datadist(dat)
options(datadist = "ddist")


# **********************************************
# --- 9.0 FIRTH'S LOGISTIC REGRESSION ----
# **********************************************

## -- 9.1 Variable list ----
variables <- list(
  efi_score         = "Overall EFI Score",
  hfscore           = "Heart Failure",
  diabetesscore     = "Diabetes",
  dyspneascore      = "Dyspnea",
  deliriumscore     = "Delirium",
  afscore           = "Atrial Fibrillation",
  renalscore        = "Renal Disease",
  anemiascore       = "Anemia",
  suscore           = "Syncope",
  fallsscore        = "Falls",
  fracturescore     = "Fracture",
  htnscore          = "Hypertension",
  hyposcore         = "Hypothyroidism",
  urinaryscore      = "Urinary Issues",
  cadscore          = "CAD",
  depressscore      = "Depression",
  cpainscore        = "Chronic Pain",
  rascore           = "Rheumatoid Arthritis",
  dizzyscore        = "Dizziness",
  valvescore        = "Valvular Disease",
  cancerscore       = "Cancer (Other)",
  osteoscore        = "Osteoporosis",
  wlscore           = "Weight Loss",
  strokescore       = "Stroke",
  pepticscore       = "Peptic Ulcer",
  pvdscore          = "Peripheral Vascular Disease",
  incontinencescore = "Incontinence",
  amiscore          = "Acute MI",
  thyroidscore      = "Thyroid Disease", 
  visionscore       = "Vision Loss",
  hearingscore      = "Hearing Loss", 
  demscore          = "Dementia",
  melanomascore     = "Melanoma", 
  bloodcancerscore  = "Blood Cancer",
  pkscore           = "Parkinson's Disease", 
  liverscore        = "Liver Disease"
)

## -- 9.2 outcome variables ----
outcomes <- c("hosp_6mo", "hosp_12mo", "mort_6mo", "mort_12mo")

## -- 9.3 Helper: run one model, get AUC + OR + P ----
run_model <- function(data, outcome, predictor) {
  ## -- 9.4 Restricted cubic spline for the continuous eFI score, plain binary for everything else ----
  if (predictor == "efi_score") {
    form <- as.formula(paste(outcome, "~ rcs(efi_score)"))
  } else {
    form <- as.formula(paste(outcome, "~", predictor))
  }
  
  ## -- 9.5 Standard GLM for AUC (predicted probabilities) ----
  fit <- tryCatch(
    glm(form, family = "binomial", data = data),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(AUC = NA, OR = NA, P = NA))
  }
  
  pred <- predict(fit, type = "response")
  roc_obj <- roc(data[[outcome]], pred, quiet = TRUE)
  auc_val <- round(auc(roc_obj), 3)
  
  ## -- 9.6 Firth's penalized logistic regression for OR and P ----
  # For eFI score, use linear term (single OR) since RCS produces multiple coefficients
  # AUC still uses RCS above — this only affects the OR/P reporting
  if (predictor == "efi_score") {
    firth_form <- as.formula(paste(outcome, "~ efi_score"))
  } else {
    firth_form <- form
  }
  
  fit_firth <- tryCatch(
    logistf(firth_form, data = data),
    error = function(e) NULL
  )
  if (is.null(fit_firth)) {
    return(list(AUC = unname(formatC(auc_val, format = "f", digits = 3)),
                OR = "NE", P = NA))
  }
  
  ## -- 9.7 OR, CI, P value directly from the logistf object ----
  or_raw <- unname(exp(fit_firth$coefficients[2]))
  ci_lo  <- unname(exp(fit_firth$ci.lower[2]))
  ci_hi  <- unname(exp(fit_firth$ci.upper[2]))
  p_val  <- unname(fit_firth$prob[2])
  
  or_str <- if (is.finite(or_raw) && is.finite(ci_lo) && is.finite(ci_hi)) {
    paste0(round(or_raw, 2), " (", round(ci_lo, 2), "-", round(ci_hi, 2), ")")
  } else "NE"
  
  p_str <- if (is.na(p_val)) NA else if (p_val < 0.001) "<0.001" else 
    unname(formatC(round(p_val, 3), format = "f", digits = 3))
  
  list(AUC = unname(formatC(auc_val, format = "f", digits = 3)),
       OR  = or_str,
       P   = p_str)
}

### -- 9.7 Loop through all 36 variables × 4 outcomes ----
results <- tibble()

for (v in seq_along(variables)) {
  var_code  <- names(variables)[v]
  var_label <- variables[[v]]
  
  ##### -- 9.7.1 Skip if the column doesn't exist or is constant ----
  if (!(var_code %in% names(dat))) next
  if (length(unique(na.omit(dat[[var_code]]))) < 2) next
  
  hosp6  <- run_model(dat, "hosp_6mo",  var_code)
  hosp12 <- run_model(dat, "hosp_12mo", var_code)
  mort6  <- run_model(dat, "mort_6mo",  var_code)
  mort12 <- run_model(dat, "mort_12mo", var_code)
  
  new_row <- tibble(
    Variable       = var_label,
    hosp_6mo_AUC   = hosp6$AUC,
    hosp_6mo_OR    = hosp6$OR,
    hosp_6mo_P     = hosp6$P,
    hosp_12mo_AUC  = hosp12$AUC,
    hosp_12mo_OR   = hosp12$OR,
    hosp_12mo_P    = hosp12$P,
    mort_6mo_AUC   = mort6$AUC,
    mort_6mo_OR    = mort6$OR,
    mort_6mo_P     = mort6$P,
    mort_12mo_AUC  = mort12$AUC,
    mort_12mo_OR   = mort12$OR,
    mort_12mo_P    = mort12$P
  )
  
  results <- bind_rows(results, new_row)
}


### -- 9.8 save csv full results all outcomes----
cat("\n========== CSV: Full Results All Outcomes ) ==========\n")
write_csv(results, "Full_Results_All_Outcomes.csv")




# **********************************************
# --- 10.0 TABLE 3: FULL RESULTS (PLAIN) ----
# **********************************************


## -- 10.1 Build gt table ----
gt_pub <- results %>%
  gt() %>%
  tab_spanner("Hospitalized < 6 Months",  c(hosp_6mo_AUC, hosp_6mo_OR, hosp_6mo_P)) %>%
  tab_spanner("Hospitalized < 12 Months", c(hosp_12mo_AUC, hosp_12mo_OR, hosp_12mo_P)) %>%
  tab_spanner("Mortality < 6 Months",     c(mort_6mo_AUC, mort_6mo_OR, mort_6mo_P)) %>%
  tab_spanner("Mortality < 12 Months",    c(mort_12mo_AUC, mort_12mo_OR, mort_12mo_P)) %>%
  cols_label(
    Variable      = "EFI Variable",
    hosp_6mo_AUC  = "AUC",  hosp_6mo_OR  = "OR (95% CI)",  hosp_6mo_P  = "P",
    hosp_12mo_AUC = "AUC",  hosp_12mo_OR = "OR (95% CI)",  hosp_12mo_P = "P",
    mort_6mo_AUC  = "AUC",  mort_6mo_OR  = "OR (95% CI)",  mort_6mo_P  = "P",
    mort_12mo_AUC = "AUC",  mort_12mo_OR = "OR (95% CI)",  mort_12mo_P = "P"
  ) %>%
  tab_header(
    title    = md("**Table: EFI Variables and Adverse Outcomes**"),
    subtitle = md(paste0("Complete Case Analysis (N = ", format(nrow(dat), big.mark = ","), ")"))
  ) %>%
  tab_source_note(md("EFI Score modeled with RCS (5 knots). Other components modeled with unadjusted Firth's penalized logistic regression.")) %>%
  tab_source_note(md("AUC = Area Under ROC Curve; OR = Odds Ratio; CI = Confidence Interval.")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = Variable) %>%
  cols_align("center", columns = -Variable) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )

### -- 10.2 save table: all binary outcomes ----
cat("\n========== TABLE: All Binary Outcomes (Table_All_Binary_Outcomes.html) ==========\n")
gtsave(gt_pub, "Table_All_Binary_Outcomes.html")


# **********************************************
# --- 11.0 TABLE 4: AUC SUMMARY (NO COLOR) ----
# **********************************************

## -- 11.1 prepare data ----
auc_table <- results %>%
  select(Variable, hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC)

## -- 11.2 build gt table ----
gt_auc <- auc_table %>%
  gt() %>%
  cols_label(
    Variable      = "EFI Component",
    hosp_6mo_AUC  = "Hosp < 6 Mo",
    hosp_12mo_AUC = "Hosp < 12 Mo",
    mort_6mo_AUC  = "Mort < 6 Mo",
    mort_12mo_AUC = "Mort < 12 Mo"
  ) %>%
  tab_header(
    title    = md("**Table: AUC Summary: EFI Variables and Adverse Outcomes**"),
    subtitle = md(paste0("Complete Case Analysis (N = ", format(nrow(dat), big.mark = ","), ")"))
  ) %>%
  tab_source_note(md("AUC = Area Under the ROC Curve. EFI Score modeled with RCS (5 knots).")) %>%
  tab_source_note(md("All other components modeled with unadjusted Firth's penalized logistic regression.")) %>%
  tab_style(
    style = list(cell_fill(color = "#E8E8F0"), cell_text(weight = "bold")),
    locations = cells_body(rows = 1)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = Variable) %>%
  cols_align("center", columns = -Variable) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )


### -- 11.3 save table: AUC summary ----
cat("\n========== TABLE: AUC summary (Table_AUC_Summary.html) ==========\n")
gtsave(gt_auc, "Table_AUC_Summary.html")



# **********************************************
# --- 12.0 COLOR PALETTE FOR HEAT MAPS ----
# **********************************************

# Red (~0.50) -> yellow (~0.57) -> green (~0.76+)
# color_pal <- scales::col_numeric(
#   palette = c("#FFB3B3", "#FFF5BA", "#B3FFB3"),
#   domain = c(0.50, 0.76)
# )

## 12.1 WFU color palette for heat maps ----
# White (~0.50) -> beige (~0.57) -> gold (~0.76+)
color_pal <- scales::col_numeric(
  palette = c("white", "#CFBF9C", "#9E7E38"),
  domain = c(0.50, 0.76)
)

## -- 12.2 Style function----
style_wfsom <- function(gt_table) {
  gt_table %>%
    tab_options(table.background.color = "white") %>%
    tab_style(
      style = list(cell_fill(color = wfu_gold), cell_text(color = "white", weight = "bold")),
      locations = cells_column_labels()
    ) %>%
    tab_style(
      style = cell_text(color = wfu_black),
      locations = cells_body()
    ) %>%
    tab_options(
      table.border.top.width = px(2),
      table.border.top.color = wfu_gold,
      table.border.bottom.width = px(2),
      table.border.bottom.color = wfu_gold,
      column_labels.border.bottom.width = px(1.5),
      column_labels.border.bottom.color = wfu_gold
    )
}


# **********************************************
# --- 13.0 TABLE 4A: COLOR-CODED AUC HEAT MAP ----
# **********************************************

## -- 13.1 prepare data ---
auc_color_data <- results %>%
  select(Variable, hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC) %>%
  mutate(across(ends_with("_AUC"), as.numeric))

## -- 13.2 build gt table ----
gt_auc_color <- auc_color_data %>%
  gt() %>%
  data_color(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    fn = color_pal,
    autocolor_text = FALSE
  ) %>%
  fmt_number(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    decimals = 3
  ) %>%
  cols_label(
    Variable      = "EFI Component",
    hosp_6mo_AUC  = "Hosp < 6 Mo",
    hosp_12mo_AUC = "Hosp < 12 Mo",
    mort_6mo_AUC  = "Mort < 6 Mo",
    mort_12mo_AUC = "Mort < 12 Mo"
  ) %>%
  tab_header(
    title    = md("**Table: AUC Heat Map: EFI Variables and Adverse Outcomes**"),
    subtitle = md(paste0("Complete Case Analysis (N = ", format(nrow(dat), big.mark = ","), ")"))
  ) %>%
  tab_source_note(md("Color scale: White = poor discrimination (AUC ~0.50), Gold = good discrimination (AUC ~0.65+).")) %>%
  tab_source_note(md("EFI Score modeled with RCS (5 knots). Other components modeled with unadjusted Firth's penalized logistic regression.")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = Variable) %>%
  cols_align("center", columns = -Variable) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )

### -- 13.3 save table: AUC summary ----
cat("\n========== TABLE: AUC summary - color coded (Table_AUC_Summary_Color_Coded.html) ==========\n")
gtsave(gt_auc_color, "Table_AUC_Color_Coded.html")


# **********************************************
# --- 14.0 TABLE 3A: COLOR-CODED FULL TABLE ----
# **********************************************


## -- 14.1 Prepare data ----
full_color_data <- results %>%
  mutate(across(ends_with("_AUC"), as.numeric))

## -- 14.2 Build gt table ----
gt_full_color <- full_color_data %>%
  gt() %>%
  style_wfsom() %>%
  tab_spanner("Hospitalized < 6 Months",   c(hosp_6mo_AUC, hosp_6mo_OR, hosp_6mo_P)) %>%
  tab_spanner("Hospitalized < 12 Months", c(hosp_12mo_AUC, hosp_12mo_OR, hosp_12mo_P)) %>%
  tab_spanner("Mortality < 6 Months",     c(mort_6mo_AUC, mort_6mo_OR, mort_6mo_P)) %>%
  tab_spanner("Mortality < 12 Months",    c(mort_12mo_AUC, mort_12mo_OR, mort_12mo_P)) %>%
  data_color(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    fn = color_pal,
    autocolor_text = FALSE
  ) %>%
  fmt_number(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    decimals = 3
  ) %>%
  cols_label(
    Variable      = "EFI Component",
    hosp_6mo_AUC  = "AUC",  hosp_6mo_OR  = "OR (95% CI)",  hosp_6mo_P  = "P",
    hosp_12mo_AUC = "AUC",  hosp_12mo_OR = "OR (95% CI)",  hosp_12mo_P = "P",
    mort_6mo_AUC  = "AUC",  mort_6mo_OR  = "OR (95% CI)",  mort_6mo_P  = "P",
    mort_12mo_AUC = "AUC",  mort_12mo_OR = "OR (95% CI)",  mort_12mo_P = "P"
  ) %>%
  tab_header(
    title    = md("**Table: EFI Variables and Adverse Outcomes**"),
    subtitle = md(paste0("Complete Case Analysis (N = ", format(nrow(dat), big.mark = ","), ")"))
  ) %>%
  tab_source_note(md("AUC cells color-coded: White = poor discrimination (~0.50), Gold = good discrimination (~0.65+).")) %>%
  tab_source_note(md("EFI Score modeled with RCS (5 knots). Other components modeled with unadjusted Firth's penalized logistic regression.")) %>%
  tab_source_note(md("AUC = Area Under ROC Curve; OR = Odds Ratio; CI = Confidence Interval.")) %>%
  tab_style(
    #bold first row (headings)
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)
  ) %>%
  tab_style(
    #bold second row (Overall EFI Score)
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = Variable) %>%
  cols_align("center", columns = -Variable) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )


### -- 14.3 save table: full color coded----
cat("\n========== TABLE: full summary - color coded (Table_Full_Color_Coded.html) ==========\n")
gtsave(gt_full_color, "Table_Full_Color_Coded.html")


# **********************************************
# --- 15.0 TABLE 5: EFI CATEGORY (3-LEVEL) AS PREDICTOR ----
# **********************************************

## -- 15.1 set up efi_category as a factor with Fit as reference ----
dat <- dat %>%
  mutate(efi_category = factor(efi_category,
                               levels = c("Fit", "Pre-Frail", "Frail")))

## -- 15.2 category model helper ----
run_category_model <- function(data, outcome) {
  
## -- 15.3 standard GLM for AUC ----
  fit <- glm(
    as.formula(paste(outcome, "~ efi_category")),
    family = "binomial", data = data
  )
  
  pred <- predict(fit, type = "response")
  roc_obj <- roc(data[[outcome]], pred, quiet = TRUE)
  auc_val <- round(auc(roc_obj), 3)
  
## -- 15.4 Firth's for ORs and P-values ----
  fit_firth <- logistf(
    as.formula(paste(outcome, "~ efi_category")),
    data = data
  )
  
## -- 15.5 Strip names to prevent "variable names limited to 10000 bytes" error ----
  or_prefrail <- unname(exp(fit_firth$coefficients[2]))
  ci_lo_pf    <- unname(exp(fit_firth$ci.lower[2]))
  ci_hi_pf    <- unname(exp(fit_firth$ci.upper[2]))
  p_prefrail  <- unname(fit_firth$prob[2])
  
  or_frail <- unname(exp(fit_firth$coefficients[3]))
  ci_lo_f  <- unname(exp(fit_firth$ci.lower[3]))
  ci_hi_f  <- unname(exp(fit_firth$ci.upper[3]))
  p_frail  <- unname(fit_firth$prob[3])
  
  fmt_or <- function(or, lo, hi) {
    paste0(round(or, 2), " (", round(lo, 2), "-", round(hi, 2), ")")
  }
  
  list(
    AUC       = formatC(auc_val, format = "f", digits = 3),
    OR_PF     = fmt_or(or_prefrail, ci_lo_pf, ci_hi_pf),
    P_PF      = fmt_p(p_prefrail),
    OR_Frail  = fmt_or(or_frail, ci_lo_f, ci_hi_f),
    P_Frail   = fmt_p(p_frail)
  )
}

## -- 15.6 Run all 4 outcomes ----
res_hosp6  <- run_category_model(dat, "hosp_6mo")
res_hosp12 <- run_category_model(dat, "hosp_12mo")
res_mort6  <- run_category_model(dat, "mort_6mo")
res_mort12 <- run_category_model(dat, "mort_12mo")

## -- 15.7 Build results table ----
category_results <- tibble(
  Comparison = c("AUC (Overall)",
                 "Pre-Frail vs Fit: OR (95% CI)",
                 "Pre-Frail vs Fit: P",
                 "Frail vs Fit: OR (95% CI)",
                 "Frail vs Fit: P"),
  
  `Hosp < 6 Mo`  = c(res_hosp6$AUC, res_hosp6$OR_PF, res_hosp6$P_PF,
                     res_hosp6$OR_Frail, res_hosp6$P_Frail),
  
  `Hosp < 12 Mo` = c(res_hosp12$AUC, res_hosp12$OR_PF, res_hosp12$P_PF,
                     res_hosp12$OR_Frail, res_hosp12$P_Frail),
  
  `Mort < 6 Mo`  = c(res_mort6$AUC, res_mort6$OR_PF, res_mort6$P_PF,
                     res_mort6$OR_Frail, res_mort6$P_Frail),
  
  `Mort < 12 Mo` = c(res_mort12$AUC, res_mort12$OR_PF, res_mort12$P_PF,
                     res_mort12$OR_Frail, res_mort12$P_Frail)
)

### -- 15.8 print category results ----
cat("\n===== EFI Category (3-Level) Results =====\n\n")
print(category_results, n = Inf)


### -- 15.9 save csv full results all outcomes----
cat("\n========== CSV: EFI category results ==========\n")
write_csv(category_results, "EFI_Category_Results.csv")

## --- 15.10 Build gt table ----
gt_cat <- category_results %>%
  gt() %>%
  cols_label(Comparison = "EFI Category Model") %>%
  tab_header(
    title    = md("**Table: EFI Frailty Category as Predictor of Adverse Outcomes**"),
    subtitle = md(paste0("Complete Case Analysis (N = ", format(nrow(dat), big.mark = ","), ")"))
  ) %>%
  tab_source_note(md("Reference group: Fit (eFI <= 0.10). Categories: Pre-Frail (0.10-0.21), Frail (>0.21).")) %>%
  tab_source_note(md("AUC from standard logistic regression. ORs and P-values from Firth's penalized logistic regression.")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = Comparison) %>%
  cols_align("center", columns = -Comparison) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )


### -- 15.11 save table EFI category results all outcomes----
cat("\n========== TABLE: EFI category results (Table_EFI_Category.html) ==========\n")
gtsave(gt_cat, "Table_EFI_Category.html")



# **********************************************
# --- 16.0 TABLE: RCS PREDICTIONS (CONTINUOUS EFI) ----
# **********************************************

options(datadist = "ddist")

## -- 16.1 percentile definitions ----
percentiles <- c(10, 25, 50, 75, 90, 95)
names(percentiles) <- c("10th", "25th", "50th", "75th", "90th", "95th")

## -- 16.2 reference point (10th percentile ≈ 0.05) ----
ref_eFi <- quantile(dat$efi_score, probs = 0.10, na.rm = TRUE)

## -- 16.3 function for RCS predictions at percentiles ----
get_rcs_predictions <- function(data, outcome) {
  
  fit <- glm(
    as.formula(paste(outcome, "~ rcs(efi_score, 5)")),
    data = data
  )
  
  results <- tibble()
  
  for (p in names(percentiles)) {
    efi_val <- quantile(data$efi_score, probs = percentiles[p] / 100, na.rm = TRUE)
    
    pred_prob <- predict(fit, newdata = data.frame(efi_score = efi_val), type = "response")
    
## 16.4 OR relative to reference (10th percentile) ----
    prob_ref <- predict(fit, newdata = data.frame(efi_score = ref_eFi), type = "response")
    log_odds_ref  <- log(prob_ref / (1 - prob_ref))
    log_odds_curr <- log(pred_prob / (1 - pred_prob))
    or_vs_ref <- exp(log_odds_curr - log_odds_ref)
    
    new_row <- tibble(
      Percentile = p,
      eFI_Value  = round(efi_val, 3),
      Pred_Prob  = as.numeric(pred_prob),
      OR_vs_Ref  = round(as.numeric(or_vs_ref), 2)
    )
    
    results <- bind_rows(results, new_row)
  }
  
  return(results)
}

## -- 16.5 run all 4 predictions ----
cat("Generating RCS predictions...\n\n")

pred_hosp6  <- get_rcs_predictions(dat, "hosp_6mo")
pred_hosp12 <- get_rcs_predictions(dat, "hosp_12mo")
pred_mort6  <- get_rcs_predictions(dat, "mort_6mo")
pred_mort12 <- get_rcs_predictions(dat, "mort_12mo")

## --- 16.6 Combine into wide table ----
rcs_final <- pred_hosp6 %>%
  select(Percentile, eFI_Value) %>%
  mutate(
    Hosp_6mo_Prob  = pred_hosp6$Pred_Prob * 100,
    Hosp_6mo_OR    = pred_hosp6$OR_vs_Ref,
    Hosp_12mo_Prob = pred_hosp12$Pred_Prob * 100,
    Hosp_12mo_OR   = pred_hosp12$OR_vs_Ref,
    Mort_6mo_Prob  = pred_mort6$Pred_Prob * 100,
    Mort_6mo_OR    = pred_mort6$OR_vs_Ref,
    Mort_12mo_Prob = pred_mort12$Pred_Prob * 100,
    Mort_12mo_OR   = pred_mort12$OR_vs_Ref
  )


## -- 16.7 save csv RCS prediction summary ----
cat("\n===== RCS Prediction Summary =====\n\n")
print(rcs_final, n = Inf)
write_csv(rcs_final, "RCS_Prediction_Table.csv")
cat("\nSaved as RCS_Prediction_Table.csv\n")

## -- 16.8 build full gt table ----
gt_rcs <- rcs_final %>%
  mutate(
    Hosp_6mo_Prob  = paste0(format(Hosp_6mo_Prob,  nsmall = 1), "%"),
    Hosp_12mo_Prob = paste0(format(Hosp_12mo_Prob, nsmall = 1), "%"),
    Mort_6mo_Prob  = paste0(format(Mort_6mo_Prob,  nsmall = 1), "%"),
    Mort_12mo_Prob = paste0(format(Mort_12mo_Prob, nsmall = 1), "%"),
    Hosp_6mo_OR   = format(Hosp_6mo_OR,  nsmall = 2),
    Hosp_12mo_OR  = format(Hosp_12mo_OR, nsmall = 2),
    Mort_6mo_OR   = format(Mort_6mo_OR,  nsmall = 2),
    Mort_12mo_OR  = format(Mort_12mo_OR, nsmall = 2)
  ) %>%
  gt() %>%
  tab_spanner("Hospitalization",  c(Hosp_6mo_Prob, Hosp_6mo_OR, Hosp_12mo_Prob, Hosp_12mo_OR)) %>%
  tab_spanner("Mortality",        c(Mort_6mo_Prob, Mort_6mo_OR, Mort_12mo_Prob, Mort_12mo_OR)) %>%
  cols_label(
    Percentile    = "eFI Percentile",
    eFI_Value     = "eFI Score",
    Hosp_6mo_Prob  = "6-mo Prob",  Hosp_6mo_OR   = "OR vs 10th",
    Hosp_12mo_Prob = "12-mo Prob", Hosp_12mo_OR  = "OR vs 10th",
    Mort_6mo_Prob  = "6-mo Prob",  Mort_6mo_OR   = "OR vs 10th",
    Mort_12mo_Prob = "12-mo Prob", Mort_12mo_OR  = "OR vs 10th"
  ) %>%
  tab_header(
    title    = md("**Table: RCS Predictions: Continuous eFI Score and Adverse Outcomes**"),
    subtitle = md(paste0("Reference = 10th percentile (eFI ≈ ", round(ref_eFi, 3), ")"))
  ) %>%
  tab_source_note(md("Modeled with restricted cubic splines (5 knots) using rms package.")) %>%
  tab_source_note(md("Probabilities are predicted risks; ORs are relative to the 10th percentile eFI score.")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = c(Percentile, eFI_Value)) %>%
  cols_align("center", columns = -c(Percentile, eFI_Value)) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )

## -- 16.9 save table RCS predictions  ----
cat("\n========== TABLE: RCS Predictions (Table_RCS_Predictions.html) ==========\n")
gtsave(gt_rcs, "Table_RCS_Predictions.html")

## -- 16.10 Build simplified gt table - only probabilities ----
gt_rcs_simple <- rcs_final %>%
  select(Percentile, eFI_Value, Hosp_6mo_Prob, Hosp_12mo_Prob, Mort_6mo_Prob, Mort_12mo_Prob) %>%
  mutate(
    Hosp_6mo_Prob  = paste0(format(Hosp_6mo_Prob,  nsmall = 1), "%"),
    Hosp_12mo_Prob = paste0(format(Hosp_12mo_Prob, nsmall = 1), "%"),
    Mort_6mo_Prob  = paste0(format(Mort_6mo_Prob,  nsmall = 1), "%"),
    Mort_12mo_Prob = paste0(format(Mort_12mo_Prob, nsmall = 1), "%")
  ) %>%
  gt() %>%
  tab_spanner("Hospitalization",  c(Hosp_6mo_Prob, Hosp_12mo_Prob)) %>%
  tab_spanner("Mortality",        c(Mort_6mo_Prob, Mort_12mo_Prob)) %>%
  cols_label(
    Percentile    = "eFI Percentile",
    eFI_Value     = "eFI Score",
    Hosp_6mo_Prob  = "6-mo Risk",
    Hosp_12mo_Prob = "12-mo Risk",
    Mort_6mo_Prob  = "6-mo Risk",
    Mort_12mo_Prob = "12-mo Risk"
  ) %>%
  tab_header(
    title    = md("**Table: Predicted Risks Across eFI Distribution (RCS)**"),
    subtitle = md(paste0("N = ", format(nrow(dat), big.mark = ",")))
  ) %>%
  tab_source_note(md("Restricted cubic splines (5 knots). Predicted probabilities from logistic regression models.")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = c(Percentile, eFI_Value)) %>%
  cols_align("center", columns = -c(Percentile, eFI_Value)) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )

## -- 16.9 save table RCS predictions  ----
cat("\n========== TABLE: RCS Predictions Simplified (Table_RCS_Predictions_Simple.html) ==========\n")
gtsave(gt_rcs_simple, "Table_RCS_Predictions_Simple.html")

# **********************************************
# --- 16.0 TABLE: RCS PREDICTIONS (CONTINUOUS EFI) WITH FIRTH ----
# **********************************************

# WFU Gold Color Palette
wfu_gold_primary <- "#9E7E38"
wfu_gold_dark <- "#8C6D2C"
wfu_gold_light <- "#F5F0E1"
wfu_black <- "#000000"
wfu_cool_gray <- "#A7A8B3"

# %% -- 16.1 percentile definitions ----
percentiles <- c(10, 25, 50, 75, 90, 95)
names(percentiles) <- c("10th", "25th", "50th", "75th", "90th", "95th")

# %% -- 16.2 reference point (10th percentile ≈ 0.05) ----
ref_efi <- quantile(dat$efi_score, probs = 0.10, na.rm = TRUE)

# %% -- 16.3 function for Firth + RCS predictions at percentiles ----
get_firth_rcs_predictions <- function(data, outcome) {
  
  ## -- 16.3.1 Fit with Firth's penalized logistic regression ----
  fit <- logistf(
    as.formula(paste(outcome, "~ rcs(efi_score, 5)")),
    data = data
  )
  
  results <- tibble()
  
  for (p in names(percentiles)) {
    efi_val <- quantile(data$efi_score, probs = percentiles[p] / 100, na.rm = TRUE)
    
    ## -- 16.3.2 Manual prediction from Firth model ----
    lp <- predict(fit, newdata = data.frame(efi_score = efi_val), type = "link")
    pred_prob <- plogis(lp)  # Convert log-odds to probability
    
    ## -- 16.3.3 OR relative to reference (10th percentile) ----
    lp_ref  <- predict(fit, newdata = data.frame(efi_score = ref_efi), type = "link")
    or_vs_ref <- exp(lp - lp_ref)
    
    new_row <- tibble(
      Percentile = p,
      eFI_Value  = round(efi_val, 3),
      Pred_Prob  = as.numeric(pred_prob),
      OR_vs_Ref  = round(as.numeric(or_vs_ref), 2)
    )
    
    results <- bind_rows(results, new_row)
  }
  
  return(results)
}

# %% -- 16.4 run all 4 predictions ----
cat("Generating Firth + RCS predictions...\n\n")

pred_hosp6  <- get_firth_rcs_predictions(dat, "hosp_6mo")
pred_hosp12 <- get_firth_rcs_predictions(dat, "hosp_12mo")
pred_mort6  <- get_firth_rcs_predictions(dat, "mort_6mo")
pred_mort12 <- get_firth_rcs_predictions(dat, "mort_12mo")

# %% -- 16.5 combine into wide table ----
rcs_final_firth <- pred_hosp6 %>%
  select(Percentile, eFI_Value) %>%
  mutate(
    Hosp_6mo_Prob  = pred_hosp6$Pred_Prob * 100,
    Hosp_6mo_OR    = pred_hosp6$OR_vs_Ref,
    Hosp_12mo_Prob = pred_hosp12$Pred_Prob * 100,
    Hosp_12mo_OR   = pred_hosp12$OR_vs_Ref,
    Mort_6mo_Prob  = pred_mort6$Pred_Prob * 100,
    Mort_6mo_OR    = pred_mort6$OR_vs_Ref,
    Mort_12mo_Prob = pred_mort12$Pred_Prob * 100,
    Mort_12mo_OR   = pred_mort12$OR_vs_Ref
  )

# %% -- 16.6 save csv ----
cat("\n===== Firth + RCS Prediction Summary =====\n\n")
print(rcs_final_firth, n = Inf)
write_csv(rcs_final_firth, "RCS_Prediction_Table_Firth.csv")
cat("\nSaved as RCS_Prediction_Table_Firth.csv\n")


# %% -- 16.7 build gt HTML table with WFU Gold styling ----
gt_rcs_firth <- rcs_final_firth %>%
  mutate(
    Hosp_6mo_Prob  = paste0(format(Hosp_6mo_Prob,  nsmall = 1), "%"),
    Hosp_12mo_Prob = paste0(format(Hosp_12mo_Prob, nsmall = 1), "%"),
    Mort_6mo_Prob  = paste0(format(Mort_6mo_Prob,  nsmall = 1), "%"),
    Mort_12mo_Prob = paste0(format(Mort_12mo_Prob, nsmall = 1), "%"),
    Hosp_6mo_OR    = format(Hosp_6mo_OR,  nsmall = 2),
    Hosp_12mo_OR   = format(Hosp_12mo_OR, nsmall = 2),
    Mort_6mo_OR    = format(Mort_6mo_OR,  nsmall = 2),
    Mort_12mo_OR   = format(Mort_12mo_OR, nsmall = 2)
  ) %>%
  gt() %>%
  style_wfsom() %>%
  tab_spanner("Hospitalization",  c(Hosp_6mo_Prob, Hosp_6mo_OR, Hosp_12mo_Prob, Hosp_12mo_OR)) %>%
  tab_spanner("Mortality",        c(Mort_6mo_Prob, Mort_6mo_OR, Mort_12mo_Prob, Mort_12mo_Prob)) %>%
  cols_label(
    Percentile    = "eFI Percentile",
    eFI_Value     = "eFI Score",
    Hosp_6mo_Prob  = "6-mo Prob",  Hosp_6mo_OR   = "OR vs 10th",
    Hosp_12mo_Prob = "12-mo Prob", Hosp_12mo_OR  = "OR vs 10th",
    Mort_6mo_Prob  = "6-mo Prob",  Mort_6mo_OR   = "OR vs 10th",
    Mort_12mo_Prob = "12-mo Prob", Mort_12mo_OR  = "OR vs 10th"
  ) %>%
  tab_header(
    title    = md("**RCS Predictions with Firth Penalization: Continuous eFI Score and Adverse Outcomes**")
   ) %>%
  tab_source_note(md("Modeled with restricted cubic splines (5 knots) using Firth's penalized logistic regression.")) %>%
  tab_source_note(md("Probabilities are predicted risks; ORs are relative to the 10th percentile eFI score.")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = c(Percentile, eFI_Value)) %>%
  cols_align("center", columns = -c(Percentile, eFI_Value)) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = wfu_gold_light
  )

# %% -- 16.8 save HTML ----
cat("\n========== TABLE: Firth RCS Predictions (HTML) ==========\n")
gtsave(gt_rcs_firth, "Table_RCS_Predictions_Firth.html")
cat("Saved as Table_RCS_Predictions_Firth.html\n")

# %% -- 16.9 create PNG visualization with WFU Gold ----
# Reshape data for plotting
plot_data_long <- rcs_final_firth %>%
  pivot_longer(
    cols = starts_with("Hosp") | starts_with("Mort"),
    names_to = "Outcome_Time",
    values_to = "Value"
  ) %>%
  mutate(
    Outcome = case_when(
      str_detect(Outcome_Time, "^Hosp_6") ~ "Hosp 6-mo",
      str_detect(Outcome_Time, "^Hosp_12") ~ "Hosp 12-mo",
      str_detect(Outcome_Time, "^Mort_6") ~ "Mort 6-mo",
      str_detect(Outcome_Time, "^Mort_12") ~ "Mort 12-mo"
    ),
    Metric = if_else(str_detect(Outcome_Time, "Prob"), "Risk (%)", "OR"),
    Value = case_when(
      str_detect(Outcome_Time, "Prob") ~ Value,
      TRUE ~ Value * 10  # Scale OR for visibility on same plot
    )
  )

png_plot <- ggplot(rcs_final_firth) +
  # Main risk lines in WFU Gold variants
  geom_line(aes(x = eFI_Value, y = Hosp_6mo_Prob, color = "Hospitalization 6-mo"), 
            linewidth = 1.5, alpha = 0.9) +
  geom_line(aes(x = eFI_Value, y = Hosp_12mo_Prob, color = "Hospitalization 12-mo"), 
            linewidth = 1.5, alpha = 0.9) +
  geom_line(aes(x = eFI_Value, y = Mort_6mo_Prob, color = "Mortality 6-mo"), 
            linewidth = 1.5, alpha = 0.9, linetype = "dashed") +
  geom_line(aes(x = eFI_Value, y = Mort_12mo_Prob, color = "Mortality 12-mo"), 
            linewidth = 1.5, alpha = 0.9, linetype = "dashed") +
  
  # WFU Gold color scale for lines
  scale_color_manual(
    values = c("Hospitalization 6-mo" = wfu_gold_primary,
               "Hospitalization 12-mo" = wfu_gold_dark,
               "Mortality 6-mo" = "#A88545",  # Medium gold variant
               "Mortality 12-mo" = "#7D632B"), # Dark gold variant
    name = ""
  ) +
  
  # Shaded area highlighting pre-frail range in light gold
  annotate("rect", xmin = 0.10, xmax = 0.21, ymin = -Inf, ymax = Inf,
           fill = wfu_gold_primary, alpha = 0.15, inherit.aes = FALSE) +
  annotate("text", x = 0.15, y = max(rcs_final_firth$Hosp_12mo_Prob) - 2,
           label = "Pre-Frail", size = 4, fontface = "italic", color = wfu_black, alpha = 0.7) +
  
  annotate("rect", xmin = 0.21, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = wfu_gold_primary, alpha = 0.20, inherit.aes = FALSE) +
  annotate("text", x = 0.30, y = max(rcs_final_firth$Hosp_12mo_Prob) - 4,
           label = "Frail", size = 4, fontface = "italic", color = wfu_black, alpha = 0.7) +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(color = wfu_black),
    axis.title = element_text(face = "bold", color = wfu_black),
    panel.grid.major = element_line(color = "#D9D9D9"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = wfu_cool_gray, linewidth = 0.5)
  ) +
  
  labs(
    title = "Predicted Hospitalization & Mortality Risk by eFI Score",
    subtitle = paste0("Firth Penalized Logistic Regression with RCS (5 knots) | Ref = 10th percentile (eFI ≈ ", round(ref_efi, 3), ")"),
    x = "eFI Score",
    y = "Predicted Probability (%)"
  )

cat("\n========== PLOT: Firth RCS Predictions (PNG) ==========\n")
ggsave("RCS_Predictions_Firth.png", plot = png_plot, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved as RCS_Predictions_Firth.png\n")

# %% -- 16.10 PNG for OR visualization with WFU Gold ----
or_plot <- ggplot(rcs_final_firth %>% pivot_longer(cols = ends_with("_OR"), names_to = "Outcome", values_to = "OR")) +
  geom_line(aes(x = eFI_Value, y = OR, color = Outcome), linewidth = 1.5) +
  geom_hline(yintercept = 1, linetype = "dotted", color = wfu_cool_gray) +
  scale_color_manual(
    values = c("Hosp_6mo_OR" = wfu_gold_primary,
               "Hosp_12mo_OR" = wfu_gold_dark,
               "Mort_6mo_OR" = "#A88545",
               "Mort_12mo_OR" = "#7D632B"),
    name = ""
  ) +
  labs(
    title = "Odds Ratios Relative to 10th Percentile (Firth + RCS)",
    x = "eFI Score",
    y = "Odds Ratio (vs 10th percentile)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(color = wfu_black),
    axis.title = element_text(face = "bold", color = wfu_black),
    panel.grid.major = element_line(color = "#D9D9D9"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = wfu_cool_gray, linewidth = 0.5)
  )

ggsave("RCS_Predictions_Firth_OR.png", plot = or_plot, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved as RCS_Predictions_Firth_OR.png\n")

# %% -- 16.11 Final summary ----
cat("\n========== RCS Firth Outputs COMPLETE ==========\n")
cat("Output files:\n")
cat("  Table_RCS_Predictions_Firth.html\n")
cat("  RCS_Predictions_Firth.png\n")
cat("  RCS_Predictions_Firth_OR.png\n")
cat("  RCS_Prediction_Table_Firth.csv\n")
# **********************************************
# --- 17.0 COLOR-CODED FULL TABLES WITH SIGNIFICANCE MARKERS ----
# **********************************************

## -- 17.1 Prepare data ----
# Convert AUC to numeric for coloring, parse P-values for significance testing
parse_p_num <- function(p) {
  if (is.na(p)) return(1)
  if (grepl("^<", p)) return(0)
  as.numeric(p)
}

full_color_data <- results %>%
  mutate(across(ends_with("_AUC"), as.numeric)) %>%
  mutate(
    hosp_6mo_P_num  = sapply(hosp_6mo_P, parse_p_num),
    hosp_12mo_P_num = sapply(hosp_12mo_P, parse_p_num),
    mort_6mo_P_num  = sapply(mort_6mo_P, parse_p_num),
    mort_12mo_P_num = sapply(mort_12mo_P, parse_p_num)
  )

## -- 17.2 Helper: add significance asterisks to OR strings ----
add_stars <- function(or_str, p_str) {
  if (is.na(p_str) || is.na(or_str)) return(or_str)
  p_num <- parse_p_num(p_str)
  if (p_num >= 0.05) return(or_str)
  if (p_num < 0.001) return(paste0(or_str, "***"))
  if (p_num < 0.01)  return(paste0(or_str, "**"))
  paste0(or_str, "*")
}



## -- 17.3 TABLE A: ASTERISKS ----


full_stars <- full_color_data %>%
  mutate(
    hosp_6mo_OR  = mapply(add_stars, hosp_6mo_OR, hosp_6mo_P),
    hosp_12mo_OR = mapply(add_stars, hosp_12mo_OR, hosp_12mo_P),
    mort_6mo_OR  = mapply(add_stars, mort_6mo_OR, mort_6mo_P),
    mort_12mo_OR = mapply(add_stars, mort_12mo_OR, mort_12mo_P)
  ) %>%
  select(-ends_with("_num"))

gt_full_stars <- full_stars %>%
  gt() %>%
  style_wfsom() %>%
  tab_spanner("Hospitalized < 6 Months",   c(hosp_6mo_AUC, hosp_6mo_OR, hosp_6mo_P)) %>%
  tab_spanner("Hospitalized < 12 Months", c(hosp_12mo_AUC, hosp_12mo_OR, hosp_12mo_P)) %>%
  tab_spanner("Mortality < 6 Months",     c(mort_6mo_AUC, mort_6mo_OR, mort_6mo_P)) %>%
  tab_spanner("Mortality < 12 Months",    c(mort_12mo_AUC, mort_12mo_OR, mort_12mo_P)) %>%
  data_color(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    fn = color_pal,
    autocolor_text = FALSE
  ) %>%
  fmt_number(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    decimals = 3
  ) %>%
  cols_label(
    Variable      = "EFI Component",
    hosp_6mo_AUC  = "AUC",  hosp_6mo_OR  = "OR (95% CI)*",  hosp_6mo_P  = "P",
    hosp_12mo_AUC = "AUC",  hosp_12mo_OR = "OR (95% CI)*",  hosp_12mo_P = "P",
    mort_6mo_AUC  = "AUC",  mort_6mo_OR  = "OR (95% CI)*",  mort_6mo_P  = "P",
    mort_12mo_AUC = "AUC",  mort_12mo_OR = "OR (95% CI)*",  mort_12mo_P = "P"
  ) %>%
  tab_header(
    title    = md("**Table: EFI Variables and Adverse Outcomes (Asterisks)**"),
    subtitle = md(paste0("Complete Case Analysis (N = ", format(nrow(dat), big.mark = ","), ")"))
  ) %>%
  tab_source_note(md("*p < 0.05, **p < 0.01, ***p < 0.001.")) %>%
  tab_source_note(md("AUC cells color-coded: White = poor discrimination (~0.50), Gold = good discrimination (~0.65+).")) %>%
  tab_source_note(md("EFI Score modeled with RCS (5 knots). Other components modeled with unadjusted Firth's penalized logistic regression.")) %>%
  tab_source_note(md("AUC = Area Under ROC Curve; OR = Odds Ratio; CI = Confidence Interval.")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = Variable) %>%
  cols_align("center", columns = -Variable) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )



## -- 17.4 TABLE B: ITALICIZED ----


full_plain <- full_color_data %>%
  select(-ends_with("_num"))

gt_full_italic <- full_plain %>%
  gt() %>%
  style_wfsom() %>%
  tab_spanner("Hospitalized < 6 Months",   c(hosp_6mo_AUC, hosp_6mo_OR, hosp_6mo_P)) %>%
  tab_spanner("Hospitalized < 12 Months", c(hosp_12mo_AUC, hosp_12mo_OR, hosp_12mo_P)) %>%
  tab_spanner("Mortality < 6 Months",     c(mort_6mo_AUC, mort_6mo_OR, mort_6mo_P)) %>%
  tab_spanner("Mortality < 12 Months",    c(mort_12mo_AUC, mort_12mo_OR, mort_12mo_P)) %>%
  data_color(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    fn = color_pal,
    autocolor_text = FALSE
  ) %>%
  fmt_number(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    decimals = 3
  ) %>%
  cols_label(
    Variable      = "EFI Component",
    hosp_6mo_AUC  = "AUC",  hosp_6mo_OR  = "OR (95% CI)",  hosp_6mo_P  = "P",
    hosp_12mo_AUC = "AUC",  hosp_12mo_OR = "OR (95% CI)",  hosp_12mo_P = "P",
    mort_6mo_AUC  = "AUC",  mort_6mo_OR  = "OR (95% CI)",  mort_6mo_P  = "P",
    mort_12mo_AUC = "AUC",  mort_12mo_OR = "OR (95% CI)",  mort_12mo_P = "P"
  ) %>%
  tab_header(
    title    = md("**Table: EFI Variables and Adverse Outcomes (Italicized)**"),
    subtitle = md(paste0("Complete Case Analysis (N = ", format(nrow(dat), big.mark = ","), ")"))
  ) %>%
  tab_source_note(md("Italic text indicates p < 0.05 (statistically significant).")) %>%
  tab_source_note(md("AUC cells color-coded: White = poor discrimination (~0.50), Gold = good discrimination (~0.65+).")) %>%
  tab_source_note(md("EFI Score modeled with RCS (5 knots). Other components modeled with unadjusted Firth's penalized logistic regression.")) %>%
  tab_source_note(md("AUC = Area Under ROC Curve; OR = Odds Ratio; CI = Confidence Interval.")) %>%
  # Italicize hosp_6mo significant cells
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = c(hosp_6mo_OR, hosp_6mo_P),
      rows = which(full_color_data$hosp_6mo_P_num < 0.05)
    )
  ) %>%
  # Italicize hosp_12mo significant cells
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = c(hosp_12mo_OR, hosp_12mo_P),
      rows = which(full_color_data$hosp_12mo_P_num < 0.05)
    )
  ) %>%
  # Italicize mort_6mo significant cells
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = c(mort_6mo_OR, mort_6mo_P),
      rows = which(full_color_data$mort_6mo_P_num < 0.05)
    )
  ) %>%
  # Italicize mort_12mo significant cells
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = c(mort_12mo_OR, mort_12mo_P),
      rows = which(full_color_data$mort_12mo_P_num < 0.05)
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = Variable) %>%
  cols_align("center", columns = -Variable) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )



## -- 17.5 TABLE C: BOLDED ----


gt_full_bold <- full_plain %>%
  gt() %>%
  style_wfsom() %>%
  tab_spanner("Hospitalized < 6 Months",   c(hosp_6mo_AUC, hosp_6mo_OR, hosp_6mo_P)) %>%
  tab_spanner("Hospitalized < 12 Months", c(hosp_12mo_AUC, hosp_12mo_OR, hosp_12mo_P)) %>%
  tab_spanner("Mortality < 6 Months",     c(mort_6mo_AUC, mort_6mo_OR, mort_6mo_P)) %>%
  tab_spanner("Mortality < 12 Months",    c(mort_12mo_AUC, mort_12mo_OR, mort_12mo_P)) %>%
  data_color(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    fn = color_pal,
    autocolor_text = FALSE
  ) %>%
  fmt_number(
    columns = c(hosp_6mo_AUC, hosp_12mo_AUC, mort_6mo_AUC, mort_12mo_AUC),
    decimals = 3
  ) %>%
  cols_label(
    Variable      = "EFI Component",
    hosp_6mo_AUC  = "AUC",  hosp_6mo_OR  = "OR (95% CI)",  hosp_6mo_P  = "P",
    hosp_12mo_AUC = "AUC",  hosp_12mo_OR = "OR (95% CI)",  hosp_12mo_P = "P",
    mort_6mo_AUC  = "AUC",  mort_6mo_OR  = "OR (95% CI)",  mort_6mo_P  = "P",
    mort_12mo_AUC = "AUC",  mort_12mo_OR = "OR (95% CI)",  mort_12mo_P = "P"
  ) %>%
  tab_header(
    title    = md("**EFI Variables and Adverse Outcomes**"),
    subtitle = md(paste0("Complete Case Analysis (N = ", format(nrow(dat), big.mark = ","), ")"))
  ) %>%
  tab_source_note(md("Bold text indicates p < 0.05 (statistically significant).")) %>%
  tab_source_note(md("AUC cells color-coded: White = poor discrimination (~0.50), Gold = good discrimination (~0.65+).")) %>%
  tab_source_note(md("EFI Score modeled with RCS (5 knots). Other components modeled with unadjusted Firth's penalized logistic regression.")) %>%
  tab_source_note(md("AUC = Area Under ROC Curve; OR = Odds Ratio; CI = Confidence Interval.")) %>%
  # Bold hosp_6mo significant cells
  tab_style(
    style = cell_text(weight = "bold"), 
    locations = cells_body(
      columns = c(hosp_6mo_AUC, hosp_6mo_OR, hosp_6mo_P),
      rows = which(full_color_data$hosp_6mo_P_num < 0.05)
    )
  ) %>%
  # Bold hosp_12mo significant cells
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(hosp_12mo_AUC, hosp_12mo_OR, hosp_12mo_P),
      rows = which(full_color_data$hosp_12mo_P_num < 0.05)
    )
  ) %>%
  # Bold mort_6mo significant cells
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(mort_6mo_AUC, mort_6mo_OR, mort_6mo_P),
      rows = which(full_color_data$mort_6mo_P_num < 0.05)
    )
  ) %>%
  # Bold mort_12mo significant cells
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(mort_12mo_AUC, mort_12mo_OR, mort_12mo_P),
      rows = which(full_color_data$mort_12mo_P_num < 0.05)
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = 1)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  cols_align("left", columns = Variable) %>%
  cols_align("center", columns = -Variable) %>%
  tab_options(
    data_row.padding = px(6),
    table.font.size = 12,
    heading.title.font.size = 16,
    source_notes.font.size = 10,
    column_labels.background.color = "#F0F0F5"
  )



## -- 17.6 Save all three outputs ----


gtsave(gt_full_stars,  "Table_Full_Color_Asterisks.html")
cat("\nAsterisks table saved as Table_Full_Color_Asterisks.html\n")

gtsave(gt_full_italic, "Table_Full_Color_Italic.html")
cat("\nItalicized table saved as Table_Full_Color_Italic.html\n")

gtsave(gt_full_bold,   "Table_Full_Color_Bold.html")
cat("\nBolded table saved as Table_Full_Color_Bold.html\n")





# **********************************************
# --- 18.0 VISUALIZE RISK GRADIENT ----
# **********************************************

## -- 18.1 Build plot ----
risk_plot <-ggplot(efi_by_category, aes(x = reorder(efi_category, hosp_6mo_rate),
                            y = hosp_6mo_rate, fill = efi_category)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(round(hosp_6mo_rate, 1), "%")), vjust = -0.5, size = 5) +
  labs(title = "Hospitalization Rate Increases with Frailty",
       x = "EFI Category",
       y = "Hospitalized Within 6 Months (%)") +
  theme_minimal() +
  scale_fill_manual(values = c("Fit" = "#4DAF4A", "Pre-Frail" = "#FF7F00", "Frail" = "#E41A1C")) +
  theme(legend.position = "none")

## -- 18.2 Save plot: frality risk gradient ----
cat("\n===== PLOT: Frality risk gradient =====\n")
ggsave("frailty_risk_gradient.png", width = 6, height = 4, dpi = 300)


## -- 18.3 Final summary ----
cat("\n========== ALL TABLES COMPLETE ==========\n")
cat("Output files:\n")
cat("  Table1_Stratified.html\n")
cat("  Table_All_Binary_Outcomes.html\n")
cat("  Table_AUC_Summary.html\n")
cat("  Table_AUC_Color_Coded.html\n")
cat("  Table_Full_Color_Coded.html\n")
cat("  Table_EFI_Category.html\n")
cat("  Table_RCS_Predictions.html\n")
cat("  Table_RCS_Predictions_Simple.html\n")
cat("  frailty_risk_gradient.png\n")
cat("  Full_Results_All_Outcomes.csv\n")
cat("  EFI_Category_Results.csv\n")
cat("  RCS_Prediction_Table.csv\n")
cat("  Table_Full_Color_Asterisks.html\n")
cat(". Table_Full_Color_Italic.html\n")
cat(". Table_Full_Color_Bold.html\n")

## -- 18.4 Export HTML tables as PNG ----
library(webshot2)

# List of your HTML output files
html_files <- c(
  # "Table1_Stratified.html",
  # "Table_All_Binary_Outcomes.html",
  # "Table_AUC_Summary.html",
  # "Table_AUC_Color_Coded.html",
  # "Table_Full_Color_Coded.html",
  # "Table_EFI_Category.html",
  # "Table_RCS_Predictions_Firth.html",
  # "Table_RCS_Predictions.html",
  # "Table_RCS_Predictions_Simple.html",
  # "Table_Full_Color_Asterisks.html",
  # "Table_Full_Color_Italic.html",
  "Table_Full_Color_Bold.html"
)

# Loop through and convert each
for (html_file in html_files) {
  if (file.exists(html_file)) {
    png_file <- sub("\\.html$", ".png", html_file)
    webshot2::webshot(html_file, png_file, delay = 2)
    cat("Saved:", png_file, "\n")
  } else {
    cat("WARNING: File not found -", html_file, "\n")
  }
}

cat("\n========== PNG EXPORT COMPLETE ==========\n")


# --- 19.0 SAVE WORKSPACE FOR DASHBOARD ----
save.image("efi_dashboard_data.RData")
cat("Workspace saved for dashboard!\n")











library(ggplot2)
library(scales)
# 1. Define your data
df <- data.frame(
  Category = factor(c("Fit", "Pre-Frail", "Frail"), 
                    levels = c("Fit", "Pre-Frail", "Frail")),
  Rate = c(8.2, 14.1, 24.3),
  lower = c(6.1, 11.8, 21.5), 
  upper = c(10.3, 16.4, 27.1)  
)

# 2. Define color
wf_gold_primary <- "#9E7E38"

# 3. Build the plot
chart <- ggplot(df, aes(x = Category, y = Rate, fill = Category)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = rep(wf_gold_primary, 3)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "black", linewidth = 0.6) +
  geom_text(aes(label = paste0(Rate, "%"), y = upper + 1.5), 
            vjust = 0, 
            size = 5, 
            fontface = "bold", 
            color = "black") +
  scale_y_continuous(
    limits = c(0, 35),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(x = "", y = "6-Month Hospitalization Rate (%)") +
  theme_minimal(base_size = 12) +
  theme(
    # Explicitly remove all borders on backgrounds
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    
    # Remove panel border and all axis lines/ticks
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    # Keep only the faint horizontal grid lines
    panel.grid.major.y = element_line(color = "#D9D9D9", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    # Text styling
    axis.text = element_text(color = "#000000", size = 11),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.margin = margin(15, 15, 15, 15)
  )
# 4. Save the chart
ggsave("frailty_risk_by_eFI.png", 
       plot = chart, 
       width = 5, 
       height = 4, 
       dpi = 300, 
       bg = "white")

print("Chart saved successfully!")