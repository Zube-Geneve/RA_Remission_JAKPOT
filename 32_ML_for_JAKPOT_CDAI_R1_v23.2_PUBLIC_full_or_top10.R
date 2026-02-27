########## S T A R T   C O D E ##########


#Libraries
library(tidyr)
library(tidyverse)
library(dplyr)
library(caret)
library(xgboost)
library(data.table)
library(mice)
library(MatchIt)
library(tableone)
library(twang)
library(survival)
library(cmprsk)
library(survminer)
library(car)
library(glmnet)
library(RANN)
library(pROC)
library(corrplot)
library(cluster)
library(fpc)
library(smotefamily)

# Set the seed for random number generation 
set.seed(7201)


# -------------------------------------------------------------------------
# PUBLIC REPOSITORY SAFETY + REPRODUCIBILITY SETTINGS
# -------------------------------------------------------------------------
# This script is designed to be shareable without sharing any registry data.
# By default it will NOT export patient-level datasets or trained model binaries.
# To run internally, set the environment variables below or edit paths locally.

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required. Please install it with install.packages('here').")
  }
})

PROJECT_DIR <- here::here()

# Helper: fail fast if a required env var is missing
need_env <- function(name) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) {
    stop(sprintf("Missing environment variable '%s'. See README for setup.", name))
  }
  normalizePath(val, winslash = "/", mustWork = FALSE)
}

# Private inputs: MUST be set locally, never committed
DATA_FILE     <- need_env("JAKPOT_RDATA_FILE")      # e.g. /secure/path/clean_data.RData
# PREV_DATA_FILE <- Sys.getenv("JAKPOT_PREV_RDATA_FILE", unset = "")  # optional (unused in public script)
HAQ_SCQM_FILE  <- Sys.getenv("JAKPOT_HAQ_SCQM_FILE", unset = "")    # optional

# Outputs: default to a temp folder to reduce accidental git commits
OUT_DIR   <- Sys.getenv("JAKPOT_OUT_DIR", unset = file.path(tempdir(), "jakpot_outputs"))
MODEL_DIR <- Sys.getenv("JAKPOT_MODEL_DIR", unset = file.path(tempdir(), "jakpot_models"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# Model specification
# ------------------------------------------------------------
# Choose which model to run:
# - "full": full feature set (as per primary analysis)
# - "top10": simplified model using the 10-variable set reported in the manuscript
MODEL_TYPE <- Sys.getenv("JAKPOT_MODEL_TYPE", unset = "full")
MODEL_TYPE <- tolower(MODEL_TYPE)
if (!MODEL_TYPE %in% c("full", "top10")) {
  stop("MODEL_TYPE must be either 'full' or 'top10'.")
}

# Top-10 base variables (pre-encoding). The script will automatically include any
# dummy-encoded columns derived from these variables (caret::dummyVars naming).
TOP10_BASE_VARS <- c(
  "PGA",
  "TJC28",
  "Prev_btsDMARD",
  "HAQ0",
  "PhGA",
  "Age",
  "CDAI0",
  "CRP",
  "Pain",
  "AntiCCPevernever"
)

# Safety flags (default OFF for public sharing)
EXPORT_PATIENT_LEVEL   <- FALSE
EXPORT_DERIVED_OUTPUTS <- FALSE
SAVE_MODEL_ARTIFACTS   <- FALSE
WRITE_LOGS             <- TRUE


fp_out   <- function(...) file.path(OUT_DIR, ...)
fp_model <- function(...) file.path(MODEL_DIR, ...)

# -------------------------------------------------------------------------

# Working directory: use an RStudio Project and relative paths via `here::here()`.
RUN_TUNING <- TRUE

# load the current dataset
load(DATA_FILE)

jakpot.data_eff1 <- jakpot.data_eff
jakpot.data_eff1_baseline <- jakpot.data_effbaseline
jakpot.data_eff1_ae <- jakpot.data_ae


# Load the corrected HAQ-DI dataset (optional, registry-specific)
# If the file is not provided, the pipeline continues using the original HAQ0 values.
haq_scqm <- NULL
if (nzchar(HAQ_SCQM_FILE) && file.exists(HAQ_SCQM_FILE)) {
  haq_scqm <- readRDS(HAQ_SCQM_FILE) %>%
    rename(ID_original = ID) %>%
    mutate(Visit_date = as.Date(Visit_date))

  # Join SCQM corrected data (only if columns exist)
  if (all(c("ID_original", "Visit_date") %in% names(jakpot.data_eff1_baseline)) &&
      all(c("ID_original", "Visit_date", "HAQ_corrected") %in% names(haq_scqm))) {

    jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>%
      left_join(haq_scqm %>% select(ID_original, Visit_date, HAQ_corrected),
                by = c("ID_original", "Visit_date"))

    # Store original HAQ0 then replace where corrected values exist
    if ("HAQ0" %in% names(jakpot.data_eff1_baseline)) {
      jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>%
        mutate(HAQ0_original = HAQ0,
               HAQ0 = ifelse(!is.na(HAQ_corrected), HAQ_corrected, HAQ0))
    }

    # Count replacements (informative only)
    if (all(c("Code", "HAQ_corrected", "HAQ0_original", "HAQ0") %in% names(jakpot.data_eff1_baseline))) {
      num_replaced <- sum(
        jakpot.data_eff1_baseline$Code == "CH" &
          !is.na(jakpot.data_eff1_baseline$HAQ_corrected) &
          jakpot.data_eff1_baseline$HAQ0_original != jakpot.data_eff1_baseline$HAQ0,
        na.rm = TRUE
      )
      message("SCQM HAQ correction applied, rows replaced: ", num_replaced)
    }

    # Cleanup helper columns
    if ("HAQ_corrected" %in% names(jakpot.data_eff1_baseline)) {
      jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>% select(-HAQ_corrected)
    }
    if ("HAQ0_original" %in% names(jakpot.data_eff1_baseline)) {
      jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>% select(-HAQ0_original)
    }
  } else {
    warning("HAQ correction file provided but required join columns are missing, skipping correction.")
  }
} else {
  message("No HAQ correction file provided, skipping SCQM HAQ correction.")
}

# View the result
num_replaced


# 
duplicated_pairs <- haq_scqm %>%
  filter(!is.na(pair_id)) %>%
  distinct(ID_original, pair_id)

duplicated_pairs %>% 
  count(pair_id) %>%
  filter(n > 1)  # Expect 9 pairs with n = 2

# keep only the lowest ID per pair_id
ids_to_keep <- duplicated_pairs %>%
  group_by(pair_id) %>%
  slice_min(ID_original) %>%
  ungroup()

ids_to_drop <- setdiff(duplicated_pairs$ID_original, ids_to_keep$ID_original)

jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>%
  filter(!ID_original %in% ids_to_drop)

# Drop HAQ_corrected
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>%
  select(-HAQ0_original)

# CREATE THE CDAI OUTCOME
# Add a row number column to identify the first visit per ID and treatment course
jakpot.data_eff1 <- jakpot.data_eff1 %>%
  arrange(ID, ttt_course, Visit_date) %>%
  group_by(ID, ttt_course) %>%
  mutate(visit_order = row_number()) %>%
  ungroup()

# Filter out the first visit for each ID and treatment course combination
followup_visits <- jakpot.data_eff1 %>%
  filter(visit_order > 1)

# Filter visits with valid CDAI within the follow-up period (6 to 24 months after treatment start)
followup_visits_filtered <- followup_visits %>%
  filter(Visit_date >= ttt_startDate + months(6), 
         Visit_date <= ttt_startDate + years(2))

# Create a dataset identifying whether remission occurred during the filtered period
remission_cdai <- followup_visits_filtered %>%
  filter(!is.na(CDAI)) %>% # Ensure CDAI is not missing
  group_by(ID, ttt_course) %>%
  summarise(
    has_remission_cdai = ifelse(any(CDAI <= 2.8), 1, 0), # Check if any valid CDAI value <= 2.8
    .groups = 'drop'
  )


# Identify early-window remission (<6 months)
early_visits <- followup_visits %>%
  filter(Visit_date < ttt_startDate + months(6), 
         Visit_date > ttt_startDate,   # exclude baseline
         !is.na(CDAI)) %>%
  group_by(ID, ttt_course) %>%
  summarise(
    early_remission = as.integer(any(CDAI <= 2.8)),
    .groups = "drop"
  )

# Identify late-window remission (6–24 months, you already did this)
late_visits <- followup_visits %>%
  filter(Visit_date >= ttt_startDate + months(6), 
         Visit_date <= ttt_startDate + years(2),
         !is.na(CDAI)) %>%
  group_by(ID, ttt_course) %>%
  summarise(
    late_remission = as.integer(any(CDAI <= 2.8)),
    .groups = "drop"
  )

# Combine both
remission_groups <- full_join(early_visits, late_visits, by = c("ID","ttt_course")) %>%
  mutate(
    early_remission = ifelse(is.na(early_remission), 0, early_remission),
    late_remission  = ifelse(is.na(late_remission), 0, late_remission),
    remission_group = case_when(
      early_remission == 1 & late_remission == 0 ~ "Early only",
      early_remission == 0 & late_remission == 1 ~ "Late only",
      early_remission == 1 & late_remission == 1 ~ "Early + Late",
      early_remission == 0 & late_remission == 0 ~ "Never"
    )
  )

# Count patients in each group
table(remission_groups$remission_group)


# Add a column to indicate if a patient has at least one valid CDAI during the follow-up period
followup_flag <- followup_visits_filtered %>%
  filter(!is.na(CDAI)) %>% # Ensure follow-up visits with valid CDAI
  group_by(ID, ttt_course) %>%
  summarise(has_followup = 1, .groups = 'drop')

# Ensure only participants with at least one valid CDAI during follow-up visits in the 6-24 month window are included
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>%
  left_join(followup_flag, by = c("ID", "ttt_course")) %>%
  left_join(remission_cdai, by = c("ID", "ttt_course")) %>%
  mutate(
    CDAI24_10 = ifelse(
      has_followup == 1, 
      has_remission_cdai, # Use remission status directly (1 for remission, 0 for non-remission)
      NA # No follow-up visits with valid CDAI
    )
  ) %>%
  # Exclude rows where CDAI24_10 is NA (no follow-up data)
  filter(!is.na(CDAI24_10))


# Check the unique country codes in the datasets
#unique(jakpot.data_eff1$Code)
unique(jakpot.data_eff1_baseline$Code)
#unique(jakpot.data_eff1_ae$Code)

#Sort the data for each ID is in the order of their visits
jakpot.data_eff1_baseline <- arrange(jakpot.data_eff1_baseline,ID,Visit_date)

#remove ttt_index because it is confusing. Denis used this for his calculation et ca sert rien.
jakpot.data_eff1_baseline = subset(jakpot.data_eff1_baseline, select = -c(ttt_index) )

#Calculate new treatment_duration_calculation because the data provided assumed the last visit date as 2023-11-10 for all countries which is not correct for some 
jakpot.data_eff1_baseline$extraction_date_calc <- jakpot.data_eff1_baseline$extraction_date

# Setting italian extraction date a bit later, because we received the data late
jakpot.data_eff1_baseline[Code %in% "IT", extraction_date_calc := as.Date("2023-12-01")]

# Focusing only on data before 2021 for NO, because we have AE only until dec 2020
jakpot.data_eff1_baseline[Code %in% "NO", extraction_date_calc := as.Date("2021-01-01")]
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline[!(Code %in% "NO") | (Code %in% "NO" & ttt_startDate <= extraction_date_calc), ]


# Calculate treatment_duration and cap it at 24 months 
jakpot.data_eff1_baseline[, treatment_duration_calc := pmin(
  ifelse(
    is.na(ttt_stopDate) & is.na(lostFUDate),
    (extraction_date_calc - ttt_startDate) / 365.25,  # Calculate based on extraction_date_calc
    ifelse(
      is.na(ttt_stopDate) & !is.na(lostFUDate),
      (lostFUDate - ttt_startDate) / 365.25,  # Calculate based on lostFUDate
      (ttt_stopDate - ttt_startDate) / 365.25  # Calculate based on ttt_stopDate
    )
  ),
  2  # Cap the duration at 2 year
), by = ID_ttt]



# Recalculate BMI
jakpot.data_eff1_baseline$BMI_recalc <- jakpot.data_eff1_baseline$Weight / ((jakpot.data_eff1_baseline$Height / 100) ^ 2)

jakpot.data_eff1_baseline$BMI <- NULL


#
# Filter for disease == "ra"
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>%
  filter(disease == "ra")


# Create the new variable CDAI_response
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>%
  mutate(CDAI_response = as.factor(CDAI24_10))




#Create afterJaki variable that represents the data took place after JAKI was introduced to that country
jakpot.data_eff1_baseline[, afterJaki := ifelse((Code == "AU"& year>=2017)|
                                                  (Code == "CA"& year>=2014)|
                                                  (Code == "CH"& year>=2013)|
                                                  (Code == "CZ"& year>=2018)|
                                                  (Code == "DE"& year>=2017)|
                                                  (Code == "FI"& year>=2017)|
                                                  (Code == "GR"& year>=2018)|
                                                  (Code == "IL"& year>=2014)|
                                                  (Code == "IT"& year>=2014)|
                                                  (Code == "NO"& year>=2017)|
                                                  (Code == "PT"& year>=2017)|
                                                  (Code == "RO"& year>=2017)|
                                                  (Code == "RU"& year>=2013)|
                                                  (Code == "SI"& year>=2017)|
                                                  (Code == "SP"& year>=2017)|
                                                  (Code == "TR"& year>=2014)|
                                                  (Code == "UK"& year>=2017)
                                                
                                                ,1,0)]

#We are not excluding the treatments before JAKI commenced so comment out the below
#jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline[jakpot.data_eff1_baseline$afterJaki==1,]


#Exclude countries that do not have CDAI info
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline[!jakpot.data_eff1_baseline$Code %in% c("DE", "SP", "TR", "UK"), ]


# Create the treatment_group column
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline %>%
  mutate(
    treatment_TNFi = generic_name %in% c("adalimumab", "certolizumab", 
                                         "etanercept", "golimumab", "infliximab", "atnf"),
    treatment_JAKi = generic_name %in% c("baricitinib", "filgotinib", "tofacitinib", "upadacitinib"),
    treatment_IL6  = (generic_name %in% c("sarilumab ", "tocilizumab")),
    # treatment_IL17  = (generic_name %in% c("bimekizumab ", "ixekizumab", "secukinumab")),
    # treatment_IL23  = (generic_name %in% c("guselkumab ", "risankizumab", "tildrakizumab", "ustekinumab")),
    treatment_ABA  = (generic_name %in% c("abatacept")),
    treatment_RITX  = (generic_name %in% c("rituximab")),
    treatment_group = case_when(
      treatment_TNFi ~ "TNFi",
      treatment_JAKi ~ "JAKi",
      treatment_IL6 ~ "IL6",
      treatment_ABA ~ "ABA",
      treatment_RITX ~ "RITX"
    )
  )



# Select only "TNFi", "JAKi", "IL6", "ABA", "RITX")
# We do not include ILL7 and IL23 due to they are not RA drugs
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline[treatment_group %in% 
                                                         c("TNFi", "JAKi", "IL6", "ABA", "RITX")]


#
#Write them down :)
if (EXPORT_DERIVED_OUTPUTS) write.csv(jakpot.data_eff1, file=fp_out("jakpot_data_eff1.csv"))
if (EXPORT_DERIVED_OUTPUTS) write.csv(jakpot.data_eff1_baseline, file=fp_out("jakpot_data_eff1_baseline.csv"))
if (EXPORT_DERIVED_OUTPUTS) write.csv(remission_cdai, file=fp_out("remission_cdai.csv"))


combined_dataset <- arrange(jakpot.data_eff1_baseline,ID,Visit_date)


### REMOVE the rows with missing outcomes
combined_dataset_original_count <- nrow(combined_dataset)

# Remove rows with missing values or NaNs in the outcome
combined_dataset <- combined_dataset[complete.cases(combined_dataset$CDAI_response), ]



combined_dataset_cleaned_count <- nrow(combined_dataset)

# Calculate the number of removed rows
combined_dataset_removed_count <- combined_dataset_original_count - combined_dataset_cleaned_count

# Print the information
cat("Number of rows before removing missing values:", combined_dataset_original_count, "\n")
cat("Number of rows removed:", combined_dataset_removed_count , "\n")
cat("Number of rows after removing missing values:", combined_dataset_cleaned_count, "\n")

##

combined_dataset_df <- combined_dataset

# Remove date variables as I do not need them

remove_date_columns <- function(data) {
  date_vars <- sapply(data, function(x) inherits(x, "Date") || inherits(x, "POSIXct"))
  non_date_columns <- names(data)[!date_vars]
  return(data[, ..non_date_columns, with = FALSE])
}

# Apply the function
combined_dataset <- remove_date_columns(combined_dataset)



# 1. Categorization of continuous variables with 10 or fewer unique values
# Create a data frame to store variables and their unique value counts
continuous_to_categorical <- data.frame(
  Variable = character(),
  UniqueValues = integer(),
  stringsAsFactors = FALSE
)

# Loop through the data frame and identify continuous variables with 10 or fewer unique values
for(colname in names(combined_dataset)) {
  if(is.numeric(combined_dataset[[colname]])) {
    unique_values <- n_distinct(combined_dataset[[colname]])
    if(unique_values <= 10) {
      continuous_to_categorical <- rbind(continuous_to_categorical, data.frame(
        Variable = colname,
        UniqueValues = unique_values
      ))
    }
  }
}

# Print the variables that will be converted along with the unique value counts
print(continuous_to_categorical)

# Write the information to a CSV file
if (EXPORT_DERIVED_OUTPUTS) write.csv(continuous_to_categorical, fp_out("continuous_to_categorical.csv"), row.names = FALSE)
# Convert these numeric variables with 10 or fewer unique values to factors
combined_dataset <- combined_dataset %>%
  mutate_at(vars(continuous_to_categorical$Variable), as.factor)


# Identify categorical and continuous variables based on the number of unique values
categorical_vars <- character()
continuous_vars <- character()

for (var in names(combined_dataset)) {
  unique_values <- unique(combined_dataset[[var]])
  if (length(unique_values) <= 10) {
    categorical_vars <- c(categorical_vars, var)
  } else {
    continuous_vars <- c(continuous_vars, var)
  }
}

# Convert character columns to factors
combined_dataset_conv <- lapply(combined_dataset, function(x) {
  if (is.character(x)) {
    as.factor(x)
  } else {
    x
  }
})

# Convert the list back to a data frame (may not be necessary)
combined_dataset_conv <- as.data.frame(combined_dataset_conv)

# Check the structure to confirm the conversion
str(combined_dataset_conv)




# Write down the dataset
if (EXPORT_DERIVED_OUTPUTS) write.csv(combined_dataset_conv, file=fp_out("combined_dataset_conv.csv"))
#Remove the useless variables
# Define the columns to exclude explicitly
exclude_columns <- c("region_id", "center_id", "Visit_date",
                     "brand_name", "ttt_startDate", "ttt_stopDate", "ttt_course",
                     "Date_birth", "Date_symptoms", "Date_diagnosis",
                     "ID_original", "ID_ttt", "extraction_date_calc", "bmi", 
                     "last_visit", "generic_name", "treatment_JAKi", "treatment_TNFi",
                     "treatment_IL6", "treatment_IL17", "treatment_IL23", "treatment_ABA", "treatment_RITX",
                     "treatment_OMA", "ttt", "ID", "bts_index", "N", "dift", "tmp",
                     "stopreasonremission", 
                     "Stop_remission", "Stop_pregnancy",
                     "stopreasonpregnancy", "Stop_ineffectiveness",
                     "stopreasonineffectiveness", "Stop_safety", 
                     "Stop_other", "stopreason", "treatment_duration", "Education",
                     "totalpatienYear", "stopany", "visitnumber", "year", "difft", "Prev_btsDMARD1", 
                     "Prev_btsDMARD3", "comorb_na", "ConcDMARD_tot", "disease", "comorb" ,
                     "Comorb_tobacco", "Comorb_tobacco_date", "lostFUDate", 
                     "cDAPSA", "skin_VAS_max", "DAPSA",  "DAPSA28",
                     "DAS28", "DAS28CRP4v", "DAS28ESR4v", "DAS28CRP3v", "DAS28ESR3v",
                     "CDAI", "HAQ", "ESR", "RF", "seropositive", "AntiCCP",
                     "DAS28ESR", "DAS28CRP", "concomitantCsDMARD", 
                     "peripheral_arthritis", "enthesitis", "dactylitis", "axial_involvement", "uveitis",
                     "colitis", "dapsa", "dapsa28", "cdapsa", "skin_vas", "skin_vas_max", "pasi", "sj66", "tj68",
                     "skin_VAS",
                     "CDAI24_10", "has_remission_cdai", "has_followup",
                     "treatment_duration_calc", "afterJaki", "concomitantCsDMARD_tot"
)

# Identify columns that start with "EQ5"
EQ5_columns <- grep("^EQ5", names(combined_dataset_conv), value = TRUE)

# Identify columns that start with "previous"
previous_columns <- grep("^previous", names(combined_dataset_conv), value = TRUE)

# Identify columns that start with "SF36_"
SF36_columns <- grep("^SF36_", names(combined_dataset_conv), value = TRUE)

# Identify columns that start with "reasonStop"
reasonStop_columns <- grep("^reasonStop", names(combined_dataset_conv), value = TRUE)

# Identify columns that start with "time"
time_columns <- grep("^time", names(combined_dataset_conv), value = TRUE)

# Identify columns that end with "_n"
n_columns <- grep("_n$", names(combined_dataset_conv), value = TRUE)

# Identify columns that end with "_12"
f12_columns <- grep("12$", names(combined_dataset_conv), value = TRUE)


# Combine both lists of columns to exclude
all_exclude_columns <- c(exclude_columns, EQ5_columns, previous_columns, SF36_columns, reasonStop_columns, time_columns, n_columns, f12_columns)

# Exclude the specified columns from filtered_data
combined_dataset_conv <- combined_dataset_conv[, !(names(combined_dataset_conv) %in% all_exclude_columns)]

original_combined_dataset_before_encoding <- combined_dataset_conv

#ONE_HOT ENCODING
# STEP 1: Specify variables to one-hot encode 
vars_to_encode <- c(
  "treatment_group", "AntiCCPevernever", "Comorb_diabetes", "Comorb_hiv", "Comorb_hyperlipidemia", 
  "Comorb_hypertension", "Comorb_inf", "Comorb_lung", "Comorb_malignancy", 
  "Comorb_myocardial", "Comorb_neuropsy", "Comorb_other_cv", "Comorb_stroke", 
  "Comorb_tobacco_ever_never", "GC", "GC_route", "HCQ", "LEF", "MTX", 
  "MTX_delivery", "Other_DMARD", "PASI", "RFevernever", "SSZ", "Sex", 
  "antiaggregant", "anticoagulant", "currentNSAIDS", "currentNSAIDScox", 
  "delivery", "educational_status", "family_history", "pastNSAIDS", "pastNSAIDScox", 
  "tsDMARD_schedule"
)

# STEP 2: Replace NAs with "missing" and convert to factor
for (col in vars_to_encode) {
  if (col %in% names(combined_dataset_conv)) {
    combined_dataset_conv[[col]] <- as.character(combined_dataset_conv[[col]])
    combined_dataset_conv[[col]][is.na(combined_dataset_conv[[col]])] <- "missing"
    combined_dataset_conv[[col]] <- as.factor(combined_dataset_conv[[col]])
  }
}

# STEP 3: Preserve outcome and grouping variables
preserve_vars <- c("CDAI_response", "Code")
preserve_data <- combined_dataset_conv[, preserve_vars]

# STEP 4: Create dummy variables
# Filter out variables with only 1 unique level
vars_to_encode_original <- vars_to_encode
vars_to_encode <- vars_to_encode[sapply(combined_dataset_conv[vars_to_encode], function(x) length(unique(x)) > 1)]
excluded_due_to_low_variation <- setdiff(vars_to_encode_original, vars_to_encode)
cat("Variables excluded from dummyVars due to single level:\n")
print(excluded_due_to_low_variation)

dummies <- dummyVars(~ ., data = combined_dataset_conv[, vars_to_encode], fullRank = FALSE)
encoded_data <- as.data.frame(predict(dummies, newdata = combined_dataset_conv[, vars_to_encode]))

# STEP 5: Drop original encoded columns, then combine everything
combined_dataset_conv <- cbind(
  combined_dataset_conv[, !(names(combined_dataset_conv) %in% vars_to_encode)],
  encoded_data,
  preserve_data  # Restore outcome and grouping variables
)

combined_dataset_conv <- combined_dataset_conv[, !grepl("^CDAI_response\\.", names(combined_dataset_conv))]
combined_dataset_conv <- combined_dataset_conv[, !grepl("^Code\\.", names(combined_dataset_conv))]

# Save
if (EXPORT_DERIVED_OUTPUTS) write.csv(combined_dataset_conv, file=fp_out("combined_dataset_conv_apres_exc.csv"), row.names = FALSE)
#Check if there are still missing predictors
na_counts <- sapply(combined_dataset_conv, function(x) sum(is.na(x)))
na_counts[na_counts > 0]


# generalize predictor variables
outcomeName <- c('CDAI_response')
predictorsNames <- names(combined_dataset_conv)[!(names(combined_dataset_conv) %in% outcomeName)]


# 3. Final Preprocessing

# Convert all NaN values to NA in the data frame
combined_dataset_conv[] <- lapply(combined_dataset_conv, function(x) {
  x[is.nan(x)] <- NA
  return(x)
})

# Remove any unintended one-hot encoded outcome columns
combined_dataset_conv <- combined_dataset_conv[, !grepl("^CDAI_response\\.", names(combined_dataset_conv))]



#Reserve Suisse data for external testing
# Subset data for external testing where the country code is 'CH'
external_data <- combined_dataset_conv[combined_dataset_conv$Code == "CH", ]

external_data_factor <- original_combined_dataset_before_encoding [original_combined_dataset_before_encoding $Code == "CH", ]

# Subset data for internal use where the country code is not 'CH'
internal_data <- combined_dataset_conv[combined_dataset_conv$Code != "CH", ]

internal_data_factor <- original_combined_dataset_before_encoding[original_combined_dataset_before_encoding$Code != "CH", ]



# Drop the 'Code' column from external_data
external_data <- external_data[, -which(names(external_data) == "Code")]



# Drop the 'Code' column from internal_data
internal_data <- internal_data[, -which(names(internal_data) == "Code")]

# ------------------------------------------------------------
# Feature set selection (applies to the MODELING dataset only)
# ------------------------------------------------------------
# Expand base variable names to include dummy-encoded columns, e.g. "AntiCCPevernever.X1"
expand_dummy_features <- function(all_names, base_vars) {
  out <- character(0)
  for (v in base_vars) {
    # exact match (numeric variables)
    if (v %in% all_names) out <- c(out, v)
    # dummyVars typically prefixes derived columns with "var." (sometimes "var_")
    out <- c(out, grep(paste0("^", gsub("([\\.^$|()\[\]{}*+?])", "\\\\1", v), "(\\.|_)"), all_names, value = TRUE))
  }
  unique(out)
}

# Determine full feature set (after preprocessing, before splitting)
target_column <- "CDAI_response"
full_feature_set <- setdiff(names(internal_data), target_column)

if (MODEL_TYPE == "top10") {
  selected_features <- expand_dummy_features(names(internal_data), TOP10_BASE_VARS)

  # Ensure we keep at least something sensible
  if (length(selected_features) == 0) {
    stop("MODEL_TYPE='top10' but none of TOP10_BASE_VARS were found in the preprocessed dataset. Check variable naming.")
  }

  # Keep only selected predictors + outcome in internal/external datasets
  keep_cols <- c(selected_features, target_column)
  internal_data <- internal_data[, intersect(keep_cols, names(internal_data))]
  external_data <- external_data[, intersect(keep_cols, names(external_data))]

  cat(sprintf("MODEL_TYPE=top10. Using %d predictors (expanded from %d base variables).\n",
              length(setdiff(names(internal_data), target_column)),
              length(TOP10_BASE_VARS)))
} else {
  selected_features <- full_feature_set
  cat(sprintf("MODEL_TYPE=full. Using %d predictors.\n", length(selected_features)))
}




# Print dimensions to verify the number of rows and columns
cat("External Data Dimensions: ", dim(external_data), "\n")
cat("Internal Data Dimensions: ", dim(internal_data), "\n")
if (EXPORT_DERIVED_OUTPUTS) write.csv(external_data, file=fp_out("external_data.csv"))
if (EXPORT_DERIVED_OUTPUTS) write.csv(internal_data, file=fp_out("internal_data.csv"))
# Split internal data into training and validation sets
set.seed(7202)

# Drop unused levels from the factor
internal_data$CDAI_response <- droplevels(internal_data$CDAI_response)


# Create a 80-20% balanced split for the outcome variable which is the country code.
splitIndex <- createDataPartition(internal_data$CDAI_response, p = .80, list = FALSE, times = 1)

# Create train and test datasets 
train_data <- internal_data[ splitIndex,]

train_data_factor <- internal_data_factor[ splitIndex,]

test_data  <- internal_data[-splitIndex,]

test_data_factor  <- internal_data_factor[-splitIndex,]
if (EXPORT_DERIVED_OUTPUTS) write.csv(train_data, file=fp_out("train_data.csv"))
if (EXPORT_DERIVED_OUTPUTS) write.csv(test_data, file=fp_out("test_data.csv"))
# ---------------------------------------------
# Percentage and number of outcome
# ---------------------------------------------
# Percentage and number of outcome in training dataset
train_data_outcome_percent <- prop.table(table(train_data$CDAI_response)) * 100
print(table(train_data$CDAI_response))
print(train_data_outcome_percent)

# Percentage and number of outcome in internal test dataset
test_data_outcome_percent <- prop.table(table(test_data$CDAI_response)) * 100
print(table(test_data$CDAI_response))
print(test_data_outcome_percent)


# Percentage and number of outcome in external test dataset
external_data_outcome_percent <- prop.table(table(external_data$CDAI_response)) * 100
print(table(external_data$CDAI_response))
print(external_data_outcome_percent)


# ---------------------------------------------
# Percentage and number of baseline remission
# ---------------------------------------------

# Create baseline remission indicator
train_data$baseline_remission <- ifelse(train_data$CDAI0 <= 2.8, 1, 0)
test_data$baseline_remission  <- ifelse(test_data$CDAI0 <= 2.8, 1, 0)
external_data$baseline_remission <- ifelse(external_data$CDAI0 <= 2.8, 1, 0)

# Training dataset
train_baseline_rem_percent <- prop.table(table(train_data$baseline_remission)) * 100
print(table(train_data$baseline_remission))
print(train_baseline_rem_percent)

# Internal test dataset
test_baseline_rem_percent <- prop.table(table(test_data$baseline_remission)) * 100
print(table(test_data$baseline_remission))
print(test_baseline_rem_percent)

# External validation dataset
external_baseline_rem_percent <- prop.table(table(external_data$baseline_remission)) * 100
print(table(external_data$baseline_remission))
print(external_baseline_rem_percent)


# ============================
# Missingness Reporting
# ============================

datasets <- list(train_data = train_data_factor, test_data = test_data_factor, external_data = external_data_factor)

for (dataset_name in names(datasets)) {
  data <- datasets[[dataset_name]]
  
  # Create empty report
  na_report <- data.frame(
    Variable = character(),
    Missing_Count = integer(),
    Missing_Percent = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var in names(data)) {
    n_miss <- sum(is.na(data[[var]]))
    if (n_miss > 0) {
      pct_miss <- round(100 * n_miss / nrow(data), 1)
      na_report <- rbind(na_report, data.frame(
        Variable = var,
        Missing_Count = n_miss,
        Missing_Percent = pct_miss
      ))
    }
  }
  
  # Save to file
if (EXPORT_DERIVED_OUTPUTS) write.csv(na_report, file = paste0(fp_out("missingness_"), dataset_name, ".csv"), row.names = FALSE)
}





####
#### Create Table One
####
library(dplyr)
library(tableone)
# Add a 'Dataset' column to identify the source in the combined data
train_data_factor$Dataset <- "Train"
test_data_factor$Dataset <- "Test"
external_data_factor$Dataset <- "External"

# Combine the data
combined_data_T1 <- rbind(train_data_factor, test_data_factor, external_data_factor)

#selected_columns <- c("Dataset", "Prev_btsDMARD", "bDMARD_currentDose", "bDMARD_currentInt", "tsDMARD_dose",
#                      "tsDMARD_schedule", "delivery", "Age", "Sex", "Height", "Weight", "BMI", "Disease_duration",
#                      "SJC28", "TJC28", "CRP", "PhGA", "PGA", "Pt_health", "Pain", "Fatigue", "MTX", "MTX_dose", 
#                      "MTX_delivery", "LEF", "LEF_dose", "SSZ", "SSZ_dose", "HCQ", "HCQ_dose", "GC", "GC_dose", 
#                      "GC_route", "Other_DMARD", "antiaggregant", "anticoagulant", "currentNSAIDS", "currentNSAIDScox",
#                      "pastNSAIDS", "pastNSAIDScox", "Comorb_tobacco_ever_never", "Comorb_myocardial", "Comorb_stroke",
#                      "Comorb_other_cv", "family_history", "Comorb_lung", "Comorb_hypertension", "Comorb_hyperlipidemia",
#                      "Comorb_hiv", "Comorb_inf", "Comorb_malignancy", "Comorb_diabetes", "Comorb_neuropsy", "PASI", 
#                      "SJ66", "TJ68", "Code", "RFevernever", "AntiCCPevernever", "educational_status", "CDAI0", "DAS280",
#                      "ESR0", "HAQ0", "CDAI_response", "treatment_group")
#combined_selected_T1 <- combined_data_T1[selected_columns]

# Create Table 1
factorVars <- c(
  "treatment_group", "AntiCCPevernever", "Comorb_diabetes", "Comorb_hiv", "Comorb_hyperlipidemia", 
  "Comorb_hypertension", "Comorb_inf", "Comorb_lung", "Comorb_malignancy", 
  "Comorb_myocardial", "Comorb_neuropsy", "Comorb_other_cv", "Comorb_stroke", 
  "Comorb_tobacco_ever_never", "GC", "GC_route", "HCQ", "LEF", "MTX", 
  "MTX_delivery", "Other_DMARD", "PASI", "RFevernever", "SSZ", "Sex", 
  "antiaggregant", "anticoagulant", "currentNSAIDS", "currentNSAIDScox", 
  "delivery", "educational_status", "family_history", "pastNSAIDS", "pastNSAIDScox", 
  "tsDMARD_schedule"
)

table1 <- CreateTableOne(
  strata = "Dataset", 
  data = combined_data_T1,
  factorVars  = factorVars
)

print(table1, printToggle = TRUE, noSpaces = TRUE)
# Convert TableOne object to a matrix format
table_matrix <- as.matrix(print(table1, printToggle = FALSE, noSpaces = TRUE, showAllLevels = TRUE))

# Capture the output of the table print
captured_output <- capture.output(
  print(table1, printToggle = TRUE, noSpaces = TRUE)
)

# Write the captured output to a text file
writeLines(captured_output, fp_out("table1.txt"))
####


# 4. Model training

# Open a file connection for writing the outputs
outputs_path <- OUT_DIR
tool_path <- PROJECT_DIR
# Create a filename with the current date
file_name <- paste("Analysis_output", Sys.Date(), "txt", sep = "_")

# Combine the path and the filename
full_path <- file.path(outputs_path, file_name)

# Open the file connection for appending
file_conn <- file(full_path, "a")

cat(paste("Start of the analysis for training", "\n"), file = file_conn)

#set seed for the CDAI_response outcome
set.seed(7203)

# Convert outcome to a factor and ensure valid R variable names for levels
# It is already factor but in case if I just repeat this code block for another outcome, and also
# make.names() ensures that the levels are syntactically valid names in R. 
train_data$CDAI_response <- factor(train_data$CDAI_response)
levels(train_data$CDAI_response) <- make.names(levels(train_data$CDAI_response))


target_column <- 'CDAI_response' 

# Convert the target column from factor to numeric 0 and 1
train_data[[target_column]] <- as.numeric(factor(train_data[[target_column]], levels = c("X0", "X1"))) - 1

# Ensure the target column is binary
if (any(!train_data[[target_column]] %in% c(0, 1))) {
  stop("Target column contains non-binary values")
}


# Convert logical columns to numeric (there are none but it is just a good practice)
logical_cols <- sapply(train_data, is.logical)
train_data[, logical_cols] <- lapply(train_data[, logical_cols], as.numeric)


# Convert factor columns to numeric 

# Start of the code for convert factor
factor_cols <- sapply(train_data, is.factor)

# Conversion function
convert_factor_to_numeric <- function(x) {
  # Check if 'x' is a factor
  if (is.factor(x)) {
    # Get unique levels from the factor
    unique_levels <- levels(x)
    
    # Create numeric labels corresponding to the index of each level
    labels <- seq(0, length(unique_levels) - 1)
    
    # Convert factor to numeric using these indices as labels
    return(as.numeric(factor(x, levels = unique_levels, labels = labels)))
  } else {
    # If not a factor, simply return x (this handles numeric and character columns that are not factors)
    return(x)
  }
}

# Apply conversion to all factor columns
train_data[factor_cols] <- lapply(train_data[factor_cols], convert_factor_to_numeric)

# End of the code for convert factor


# Convert character columns to numeric (there is none but it i just a good practice
char_cols <- sapply(train_data, is.character)
train_data[, char_cols] <- lapply(train_data[, char_cols], function(x) as.numeric(as.factor(x)))

# Ensure no non-numeric columns are left before proceeding
if (any(sapply(train_data[, -which(names(train_data) == target_column)], is.logical))) {
  stop("Non-numeric columns still exist in the dataset.")
}


# Identify low variance features in the training predictors (never drop the outcome)
predictor_cols <- setdiff(names(train_data), target_column)
low_variance_metrics <- nearZeroVar(train_data[, predictor_cols, drop = FALSE], saveMetrics = TRUE)
low_variance_features <- rownames(low_variance_metrics[low_variance_metrics$nzv, , drop = FALSE])

# Remove low variance features from training and test/external predictors
if (length(low_variance_features) > 0) {
  train_data <- train_data[, !names(train_data) %in% low_variance_features, drop = FALSE]
  test_data  <- test_data[,  !names(test_data)  %in% low_variance_features, drop = FALSE]
  if (exists("external_data")) {
    external_data <- external_data[, !names(external_data) %in% low_variance_features, drop = FALSE]
  }
}



#################
# TRAINING METHOD - DIRECT USE OF XGBOOST

cat(paste("Starting training analysis with direct use of xgboost", "\n\n"), file = file_conn)


# Create the matrix excluding the target column
train_matrix <- as.matrix(train_data[, -which(names(train_data) == target_column)])

# Create DMatrix
dtrain <- xgb.DMatrix(data = train_matrix, label = train_data[[target_column]])


# Extreme Gradient Boosting (XGB) Model: Imbalanced data
params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eta = 0.01,
  max_depth = 6,
  subsample = 0.7,
  colsample_bytree = 0.7
)

nrounds <- 100  # number of boosting rounds

cat(paste("Initial parameters for training: ", params, nrounds, "\n\n"), file = file_conn)


model <- xgb.train(params = params, data = dtrain, nrounds = nrounds)


# Make predictions
preds <- predict(model, dtrain)

# Convert predictions to binary labels
pred_labels <- ifelse(preds > 0.5, 1, 0)

# Specify that the positive class is 1
confusion_matrix <- confusionMatrix(as.factor(pred_labels), as.factor(train_data$CDAI_response), positive = "1")
print(confusion_matrix)
cat(paste("Training - Confusion Matrix: ", confusion_matrix, "\n\n"), file = file_conn)


# Compute AUC
roc_curve <- roc(train_data$CDAI_response, preds)
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))
cat(paste("Training - AUC: ", auc_value, "\n\n"), file = file_conn)


# SHAP analysis
library(SHAPforxgboost)


# Extract predictor data from train_data
X_train_df <- as.data.frame(train_data[, -which(names(train_data) == target_column)])

# Convert to numeric matrix
X_train_mat <- data.matrix(X_train_df)

# Check structure (optional sanity check)
str(X_train_mat)


# Compute SHAP values from the final model
shap_values <- shap.prep(xgb_model = model, X_train = X_train_mat)

# Save SHAP summary to file
if (EXPORT_DERIVED_OUTPUTS) write.csv(shap_values, file = fp_out("shap_values_training_feature.csv"), row.names = FALSE)
# Plot SHAP summary to PDF
pdf(fp_out("SHAP_Summary_Plot_training_feature.pdf"))
shap.plot.summary(shap_values)

#SHAP summary Plot for top 10
# 1. Compute mean absolute SHAP value for ranking
mean_shap <- shap_values %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value), na.rm = TRUE)) %>%
  arrange(desc(mean_abs_shap))

# 2. Get top 10 variable names
top_vars <- mean_shap$variable[1:10]

# 3. Subset the SHAP object *and* redefine factor levels
shap_top10 <- shap_values %>%
  filter(variable %in% top_vars) %>%
  mutate(variable = factor(variable, levels = top_vars))  # reverse order for plot

# 4. Plot only these 10
pdf(fp_out("SHAP_Summary_Plot_training_feature_top10.pdf"), width = 8, height = 6)
shap.plot.summary(shap_top10)

dev.off()

cat(paste("Training - Shap summary plot saved at: ", "\n\n"), file = file_conn)


#Save the model
if (SAVE_MODEL_ARTIFACTS) xgb.save(model, fp_model("xgb_model_full_no_finetune.bin"))
cat(paste("End of Training. ", "\n\n"), file = file_conn)


# END TRAINING METHOD - DIRECT USE OF XGBOOST
#################

#################
# FINE-TUNING HYPERPARAMETERS


#set seed for fine-tuning
set.seed(7204)

cat(paste("Starting Training fine-tuning : ", "\n\n"), file = file_conn)

if (RUN_TUNING) {
#BEGIN of tuning 
## Define a more comprehensive grid of hyperparameters
  tuning_grid <- expand.grid(
    eta = c(0.01, 0.05, 0.1),
    max_depth = c(3, 6, 9),
    gamma = c(0, 0.1, 0.2),
    colsample_bytree = c(0.5, 0.8),
    min_child_weight = c(1, 3, 5),
    subsample = c(0.6, 0.8, 1.0),
    nrounds = c(50, 100, 150)
  )
#
## Use xgb.cv for cross-validation with each parameter set
  cv_results <- lapply(1:nrow(tuning_grid), function(i) {
    params <- list(
      booster = "gbtree",
      objective = "binary:logistic",
      eta = tuning_grid[i, "eta"],
      max_depth = tuning_grid[i, "max_depth"],
      gamma = tuning_grid[i, "gamma"],
      colsample_bytree = tuning_grid[i, "colsample_bytree"],
      min_child_weight = tuning_grid[i, "min_child_weight"],
      subsample = tuning_grid[i, "subsample"]
    )
  
    xgb.cv(
      params = params,
      data = dtrain,
      nrounds = tuning_grid[i, "nrounds"],
      nfold = 5,
      metrics = "auc",
      showsd = TRUE,
      stratified = TRUE,
      print_every_n = 50,
      early_stopping_rounds = 10,
      maximize = TRUE
    )
  })
#
### Review the results to find the best performing set of parameters
  best_model_index <- which.max(sapply(cv_results, function(x) max(x$evaluation_log$test_auc_mean)))
  best_params <- tuning_grid[best_model_index, ]
  print(best_params)
##
##
##
## Update the parameter list with the best parameters
  params <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eta = best_params$eta,
    max_depth = best_params$max_depth,
    gamma = best_params$gamma,
    colsample_bytree = best_params$colsample_bytree,
    min_child_weight = best_params$min_child_weight,
    subsample = best_params$subsample
  )
#
### Define the number of boosting rounds to the best or a specified value
  best_nrounds = best_params$nrounds  # Using the number of rounds from the tuning results
#
#
#END of tuning 

} else {
  
## hard-coded params and best_nrounds
  best_eta <- 0.05
  best_max_depth <- 9
  best_gamma <- 0
  best_colsample_bytree <- 0.5
  best_min_child_weight <- 5
  best_subsample <- 1
  best_nrounds <- 150

# Update the parameter list with the best parameters
  params <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eta = best_eta,
    max_depth = best_max_depth,
    gamma = best_gamma,
    colsample_bytree = best_colsample_bytree,
    min_child_weight = best_min_child_weight,
    subsample = best_subsample
  )
}
#
cat(paste("Training fine-tuned - Best parameters : ", params, best_nrounds, "\n\n"), file = file_conn)


# Train the final model with the best number of rounds
final_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds
)


# Save the model for later use
if (SAVE_MODEL_ARTIFACTS) xgb.save(final_model, fp_model("final_xgb_model.bin"))
# Make predictions for tuned_model
preds <- predict(final_model, dtrain)

# Convert predictions to binary labels
pred_labels <- ifelse(preds > 0.5, 1, 0)

# Specify that the positive class is 1
confusion_matrix <- confusionMatrix(as.factor(pred_labels), as.factor(train_data$CDAI_response), positive = "1")
print(confusion_matrix)
cat(paste("Training fine-tuned - Confusion Matrix: ", confusion_matrix, "\n\n"), file = file_conn)


# Compute AUC
roc_curve <- roc(train_data$CDAI_response, preds)
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))
cat(paste("Training fine-tuned - AUC: ", auc_value, "\n\n"), file = file_conn)


# Find the optimal threshold
optimal_threshold <- coords(roc_curve, "best", ret = "threshold")$threshold
print(paste("Optimal Threshold:", optimal_threshold))
cat(paste("Training fine-tuned - Optimal Threshold: ", optimal_threshold, "\n\n"), file = file_conn)

# Convert predictions to binary labels using the optimal threshold
pred_labels_optimal <- ifelse(preds > optimal_threshold, 1, 0)


# Recalculate confusion matrix with the optimal threshold
confusion_matrix_train_optimal <- confusionMatrix(as.factor(pred_labels_optimal), as.factor(train_data$CDAI_response), positive = "1")
print(confusion_matrix_train_optimal)
cat(paste("Training fine-tuned with Optimal Threshold - Confusion Matrix: ", confusion_matrix_train_optimal, "\n\n"), file = file_conn)
f1_score <- confusion_matrix_train_optimal$byClass["F1"]
print(paste("F1 score for Positive Class (class '1') on Training Data with Optimal Threshold:", f1_score))
cat(paste("F1 score for Positive Class (class '1') on Training Data with Optimal Threshold: ", f1_score, "\n\n"), file = file_conn)

#Now calculate the metrics for negative class
cm_table <- confusion_matrix_train_optimal$table
TN <- cm_table[1, 1]  # True Negatives: correct predictions of class "0"
FP <- cm_table[1, 2]  # False Positives: actual "0" predicted as "1"
FN <- cm_table[2, 1]  # False Negatives: actual "1" predicted as "0"
TP <- cm_table[2, 2]  # True Positives: correct predictions of class "1"
NPV <- TN / (TN + FN)
Specificity <- TN / (TN + FP)
F1_neg <- 2 * (NPV * Specificity) / (NPV + Specificity)
print(paste("F1 Score for Negative Class (class '0') on Training Data with Optimal Threshold):", F1_neg))
cat(paste("F1 score for Negative Class (class '0') on Training Data with Optimal Threshold: ", F1_neg, "\n\n"), file = file_conn)


# Compute AUC with the optimal threshold on training data
auc_value_train_optimal <- auc(roc_curve)
print(paste("AUC on Training Data with Optimal Threshold:", auc_value_train_optimal))
cat(paste("Training fine-tuned with Optimal Threshold - AUC: ", auc_value_train_optimal, "\n\n"), file = file_conn)

# Calculate the confidence intervals for the ROC curve
ci_tuned <- ci(auc_value_train_optimal, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
print(paste("95% Confidence Intervals for AUC on Training Data with optimum threshold:"))
print(paste("Lower bound:", ci_tuned[1]))
print(paste("Upper bound:", ci_tuned[3]))
cat(paste("95% Confidence Intervals for AUC on Training Data with optimum threshold: ", "\n"), file = file_conn)
cat(paste("Lower bound: ", ci_tuned[1], "\n"), file = file_conn)
cat(paste("Upper bound: ", ci_tuned[3], "\n"), file = file_conn)


## ROC plots start
# Let's draw the ROC curve

# Define the full path for the PDF file
pdf_file_xgb_train_opt <- file.path(outputs_path, "roc_curve_plot_xgb_orig_TRAIN_opt.pdf")

# Open the PDF device
pdf(pdf_file_xgb_train_opt)


# Create a data frame with sensitivity and specificity values
roc_data <- data.frame(Sensitivity = roc_curve$sensitivities,
                       Specificity = 1 - roc_curve$specificities)


# Plot the ROC curve
plot(roc_data$Specificity, roc_data$Sensitivity, type = "l",
     main = "ROC Curve", xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)", xlim = c(0, 1), ylim = c(0, 1))

# Add a diagonal line representing random classifier (AUC = 0.500)
abline(0, 1, col = "gray", lty = 2)


# Retrieve the TPR and FPR at the optimal threshold
selected_tpr <- roc_curve$sensitivities[roc_curve$thresholds == optimal_threshold]
selected_fpr <- 1 - roc_curve$specificities[roc_curve$thresholds == optimal_threshold]


# Print the selected TPR and FPR
cat("Selected TPR: ", selected_tpr, "\n", file = file_conn)
cat("Selected FPR: ", selected_fpr, "\n", file = file_conn)



# Calculate ROC curve values for optimal threshold
roc_optimal <- roc(train_data[, outcomeName], pred_labels_optimal)



# Create a data frame with sensitivity and specificity values for the optimal roc
roc_data_optimal <- data.frame(Sensitivity = roc_optimal$sensitivities,
                               Specificity = 1 - roc_optimal$specificities)


##
# Find the closest index for the selected FPR and TPR in the optimal ROC data
closest_spec_idx <- which.min(abs(roc_data_optimal$Specificity - selected_fpr))
closest_sens_idx <- which.min(abs(roc_data_optimal$Sensitivity - selected_tpr))

# Since the indices might differ, ensure using only one index based on a primary criterion, e.g., specificity
selected_spec <- roc_data_optimal$Specificity[closest_spec_idx]
selected_sens <- roc_data_optimal$Sensitivity[closest_spec_idx]

# Add point for the selected TPR and FPR on the ROC curve
points(selected_spec, selected_sens, pch = 20, col = "red")

if (length(selected_spec) > 0 && length(selected_sens) > 0) {
  points(selected_spec, selected_sens, pch = 20, col = "red")
} else {
  cat("No close matches found for selected TPR and FPR on the ROC curve.\n")
}

##

# Add point for the selected TPR and FPR on the ROC curve
#points(roc_data_optimal$Specificity[roc_data_optimal$Specificity == selected_fpr],
#       roc_data_optimal$Sensitivity[roc_data_optimal$Sensitivity == selected_tpr],
#       pch = 20, col = "red")


formatted_optimal_threshold <- sprintf("%.3f", optimal_threshold)
formatted_tpr <- sprintf("%.3f", selected_tpr)
formatted_fpr <- sprintf("%.3f", selected_fpr)

text(selected_fpr + 0.10, selected_tpr - 0.10, 
     labels = paste("Threshold =", formatted_optimal_threshold, 
                    "\nTPR =", formatted_tpr, 
                    "\nFPR =", formatted_fpr),
     cex = 1.2, col = "blue", adj = c(0, 0))


cat("roc_data_optimal$Sensitivity: ", roc_data_optimal$Sensitivity, "\n\n", file = file_conn)
cat("roc_data_optimal$Specificity: ", roc_data_optimal$Specificity, "\n\n", file = file_conn)

# Add legend with AUC and 95% CI and F1-score for positive class
# Get the rounded value of F1 score and auc_value
f1_value_str <- sprintf("%.3f", f1_score)
F1_neg_str <- sprintf("%.3f", F1_neg)
auc_value <- round(auc_value_train_optimal, 3)
ci <- ci_tuned 

# Construct the legend text to include AUC, its 95% CI, and the F1 score
legend_text <- c(
  paste("AUC: ", auc_value, " (", round(ci[1], 3), " to ", round(ci[3], 3), ")", sep = ""),
  paste("F1 (positive class): ", f1_value_str),
  paste("F1 (negative class): ", F1_neg_str)
)

# Add the legend to the plot
legend("bottomright", legend = legend_text,
       lty = c(0, 0), col = "black", bty = "n", cex = 1.2)





# Close the PDF device and save the plot
dev.off()

## ROC plots end


# SHAP analysis
library(SHAPforxgboost)


# Extract predictor data from train_data
X_train_df <- as.data.frame(train_data[, -which(names(train_data) == target_column)])

# Convert to numeric matrix
X_train_mat <- data.matrix(X_train_df)

# Check structure (optional sanity check)
str(X_train_mat)

# Compute SHAP values from the final model
shap_values <- shap.prep(xgb_model = final_model, X_train = X_train_mat)

# Save SHAP summary to file
if (EXPORT_DERIVED_OUTPUTS) write.csv(shap_values, fp_out("shap_values_finetuned_optimal_threshold.csv"), row.names = FALSE)
# Plot SHAP summary to PDF
pdf(fp_out("SHAP_Summary_Plot_finetuned_optimal_threshold.pdf"))
shap.plot.summary(shap_values)


#SHAP summary Plot for top 10
# 1. Compute mean absolute SHAP value for ranking
mean_shap <- shap_values %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value), na.rm = TRUE)) %>%
  arrange(desc(mean_abs_shap))

# 2. Get top 10 variable names
top_vars <- mean_shap$variable[1:10]

# 3. Subset the SHAP object *and* redefine factor levels
shap_top10 <- shap_values %>%
  filter(variable %in% top_vars) %>%
  mutate(variable = factor(variable, levels = top_vars))  # reverse order for plot

# 4. Plot only these 10
pdf(fp_out("SHAP_Summary_Plot_finetuned_optimal_threshold_top10.pdf"), width = 8, height = 6)
shap.plot.summary(shap_top10)


dev.off()

cat(paste("Training fine-tuned with Optimal Threshold - SHAP summary plot saved: ",  "\n\n"), file = file_conn)

cat(paste("End of Training fine-tuning. ", "\n\n"), file = file_conn)


# END - FINE-TUNING HYPERPARAMETERS
#################



#################
#  EVALUATE USING TEST DATASET 

#set seed for the CDAI_response outcome
set.seed(7205)

cat(paste("Starting evaluation using Test dataset : ", "\n\n"), file = file_conn)


# Convert outcome to a factor and ensure valid R variable names for levels
# It is already factor but in case if I just repeat this code block for another outcome, and also
# make.names() ensures that the levels are syntactically valid names in R. 
test_data$CDAI_response <- factor(test_data$CDAI_response)
levels(test_data$CDAI_response) <- make.names(levels(test_data$CDAI_response))

target_column <- 'CDAI_response' 

# Convert the target column from factor to numeric 0 and 1
test_data[[target_column]] <- as.numeric(factor(test_data[[target_column]], levels = c("X0", "X1"))) - 1

# Ensure the target column is binary
if (any(!test_data[[target_column]] %in% c(0, 1))) {
  stop("Target column contains non-binary values")
}


# Convert logical columns to numeric
logical_cols <- sapply(test_data, is.logical)
test_data[, logical_cols] <- lapply(test_data[, logical_cols], as.numeric)


# Convert factor columns to numeric
factor_cols <- sapply(test_data, is.factor)


# Apply conversion to all factor columns
test_data[factor_cols] <- lapply(test_data[factor_cols], convert_factor_to_numeric)


# Convert character columns to numeric if there are any
char_cols <- sapply(test_data, is.character)
test_data[, char_cols] <- lapply(test_data[, char_cols], function(x) as.numeric(as.factor(x)))

# Ensure no non-numeric columns are left before proceeding
if (any(sapply(test_data[, -which(names(test_data) == target_column)], is.logical))) {
  stop("Non-numeric columns still exist in the dataset.")
}

# Remove low variance features from both training and test data
test_data <- test_data[, !names(test_data) %in% low_variance_features]

# Create the matrix excluding the target column
test_matrix <- as.matrix(test_data[, -which(names(test_data) == target_column)])

# Create DMatrix
dtest <- xgb.DMatrix(data = test_matrix, label = test_data[[target_column]])

# Make predictions
preds <- predict(final_model, dtest)

# Convert predictions to binary labels
pred_labels <- ifelse(preds > 0.5, 1, 0)

# Specify that the positive class is 1
confusion_matrix <- confusionMatrix(as.factor(pred_labels), as.factor(test_data$CDAI_response), positive = "1")
print(confusion_matrix)
cat(paste("Test No Optimal Threshold - Confusion Matrix: ", confusion_matrix, "\n\n"), file = file_conn)


# Compute AUC
roc_curve <- roc(test_data$CDAI_response, preds)
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))
cat(paste("Test No Optimal Threshold - AUC: ", auc_value, "\n\n"), file = file_conn)


#OPTIMAL threshold
# Convert predictions to binary labels using the optimal threshold
pred_labels_optimal <- ifelse(preds > optimal_threshold, 1, 0)

# Specify that the positive class is 1
confusion_matrix_test_optimal <- confusionMatrix(as.factor(pred_labels_optimal), as.factor(test_data$CDAI_response), positive = "1")
print(confusion_matrix_test_optimal)
cat(paste("Test with Optimal Threshold - Confusion Matrix: ", confusion_matrix_test_optimal, "\n\n"), file = file_conn)

f1_score <- confusion_matrix_test_optimal$byClass["F1"]
print(paste("F1 score for Positive Class (class '1') on Test Data with Optimal Threshold:", f1_score))
cat(paste("F1 score for Positive Class (class '1') on Test Data with Optimal Threshold: ", f1_score, "\n\n"), file = file_conn)

#Now calculate the metrics for negative class
cm_table <- confusion_matrix_test_optimal$table
TN <- cm_table[1, 1]  # True Negatives: correct predictions of class "0"
FP <- cm_table[1, 2]  # False Positives: actual "0" predicted as "1"
FN <- cm_table[2, 1]  # False Negatives: actual "1" predicted as "0"
TP <- cm_table[2, 2]  # True Positives: correct predictions of class "1"
NPV <- TN / (TN + FN)
Specificity <- TN / (TN + FP)
F1_neg <- 2 * (NPV * Specificity) / (NPV + Specificity)
print(paste("F1 Score for Negative Class (class '0') on Test Data with Optimal Threshold):", F1_neg))
cat(paste("F1 score for Negative Class (class '0') on Test Data with Optimal Threshold: ", F1_neg, "\n\n"), file = file_conn)



# Compute AUC
roc_curve <- roc(test_data$CDAI_response, preds)
auc_value_test_optimal <- auc(roc_curve)
print(paste("AUC on Test Data with Optimal Threshold::", auc_value_test_optimal))
cat(paste("Test with Optimal Threshold - AUC: ", auc_value_test_optimal, "\n\n"), file = file_conn)

# Calculate the confidence intervals for the ROC curve
ci_tuned <- ci(auc_value_test_optimal, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
print(paste("95% Confidence Intervals for AUC on Test Data with optimum threshold:"))
print(paste("Lower bound:", ci_tuned[1]))
print(paste("Upper bound:", ci_tuned[3]))
cat(paste("95% Confidence Intervals for AUC on Test Data with optimum threshold: ", "\n"), file = file_conn)
cat(paste("Lower bound: ", ci_tuned[1], "\n"), file = file_conn)
cat(paste("Upper bound: ", ci_tuned[3], "\n"), file = file_conn)


## ROC plots start
# Let's draw the ROC curve

# Define the full path for the PDF file
pdf_file_xgb_test_opt <- file.path(outputs_path, "roc_curve_plot_xgb_orig_TEST_opt.pdf")

# Open the PDF device
pdf(pdf_file_xgb_test_opt)


# Create a data frame with sensitivity and specificity values
roc_data <- data.frame(Sensitivity = roc_curve$sensitivities,
                       Specificity = 1 - roc_curve$specificities)


# Plot the ROC curve
plot(roc_data$Specificity, roc_data$Sensitivity, type = "l",
     main = "ROC Curve", xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)", xlim = c(0, 1), ylim = c(0, 1))

# Add a diagonal line representing random classifier (AUC = 0.500)
abline(0, 1, col = "gray", lty = 2)


# Find index of the closest threshold
closest_idx <- which.min(abs(roc_curve$thresholds - optimal_threshold))

# Retrieve the TPR and FPR at the closest threshold
selected_tpr <- roc_curve$sensitivities[closest_idx]
selected_fpr <- 1 - roc_curve$specificities[closest_idx]


# Print the selected TPR and FPR
cat("Selected TPR: ", selected_tpr, "\n", file = file_conn)
cat("Selected FPR: ", selected_fpr, "\n", file = file_conn)



# Calculate ROC curve values for optimal threshold
roc_optimal <- roc(test_data[, outcomeName], pred_labels_optimal)

# Create a data frame with sensitivity and specificity values for the optimal roc
roc_data_optimal <- data.frame(Sensitivity = roc_optimal$sensitivities,
                               Specificity = 1 - roc_optimal$specificities)



##
# Find the closest index for the selected FPR and TPR in the optimal ROC data
closest_spec_idx <- which.min(abs(roc_data_optimal$Specificity - selected_fpr))
closest_sens_idx <- which.min(abs(roc_data_optimal$Sensitivity - selected_tpr))

# Since the indices might differ, ensure using only one index based on a primary criterion, e.g., specificity
selected_spec <- roc_data_optimal$Specificity[closest_spec_idx]
selected_sens <- roc_data_optimal$Sensitivity[closest_spec_idx]

# Add point for the selected TPR and FPR on the ROC curve
points(selected_spec, selected_sens, pch = 20, col = "red")

if (length(selected_spec) > 0 && length(selected_sens) > 0) {
  points(selected_spec, selected_sens, pch = 20, col = "red")
} else {
  cat("No close matches found for selected TPR and FPR on the ROC curve.\n")
}

##

formatted_optimal_threshold <- sprintf("%.3f", optimal_threshold)
formatted_tpr <- sprintf("%.3f", selected_tpr)
formatted_fpr <- sprintf("%.3f", selected_fpr)

text(selected_fpr + 0.10, selected_tpr - 0.10, 
     labels = paste("Threshold =", formatted_optimal_threshold, 
                    "\nTPR =", formatted_tpr, 
                    "\nFPR =", formatted_fpr),
     cex = 1.2, col = "blue", adj = c(0, 0))


cat("roc_data_optimal$Sensitivity: ", roc_data_optimal$Sensitivity, "\n\n", file = file_conn)
cat("roc_data_optimal$Specificity: ", roc_data_optimal$Specificity, "\n\n", file = file_conn)

# Add legend with AUC and 95% CI and F1-score for positive class
# Get the rounded value of F1 score and auc_value
f1_value_str <- sprintf("%.3f", f1_score)
F1_neg_str <- sprintf("%.3f", F1_neg)
auc_value <- round(auc_value_test_optimal, 3)
ci <- ci_tuned 

# Construct the legend text to include AUC, its 95% CI, and the F1 score
legend_text <- c(
  paste("AUC: ", auc_value, " (", round(ci[1], 3), " to ", round(ci[3], 3), ")", sep = ""),
  paste("F1 (positive class): ", f1_value_str),
  paste("F1 (negative class): ", F1_neg_str)
)

# Add the legend to the plot
legend("bottomright", legend = legend_text,
       lty = c(0, 0), col = "black", bty = "n", cex = 1.2)





# Close the PDF device and save the plot
dev.off()

## ROC plots end

cat(paste("End of evaluation using Test dataset. ", "\n\n"), file = file_conn)


#  END - EVALUATE USING TEST DATASET 
#################



#################
#  EVALUATE USING EXTERNAL DATASET (SUISSE)

#set seed for the CDAI_response outcome
set.seed(7206)

cat(paste("Starting evaluation using External dataset : ", "\n\n"), file = file_conn)


# Convert outcome to a factor and ensure valid R variable names for levels
# It is already factor but in case if I just repeat this code block for another outcome, and also
# make.names() ensures that the levels are syntactically valid names in R. 
external_data$CDAI_response <- factor(external_data$CDAI_response)
levels(external_data$CDAI_response) <- make.names(levels(external_data$CDAI_response))

target_column <- 'CDAI_response' 

# Convert the target column from factor to numeric 0 and 1
external_data[[target_column]] <- as.numeric(factor(external_data[[target_column]], levels = c("X0", "X1"))) - 1

# Ensure the target column is binary
if (any(!external_data[[target_column]] %in% c(0, 1))) {
  stop("Target column contains non-binary values")
}


# Convert logical columns to numeric
logical_cols <- sapply(external_data, is.logical)
external_data[, logical_cols] <- lapply(external_data[, logical_cols], as.numeric)


# Convert factor columns to numeric
factor_cols <- sapply(external_data, is.factor)


# Apply conversion to all factor columns
external_data[factor_cols] <- lapply(external_data[factor_cols], convert_factor_to_numeric)


# Convert character columns to numeric if there are any
char_cols <- sapply(external_data, is.character)
external_data[, char_cols] <- lapply(external_data[, char_cols], function(x) as.numeric(as.factor(x)))

# Ensure no non-numeric columns are left before proceeding
if (any(sapply(external_data[, -which(names(external_data) == target_column)], is.logical))) {
  stop("Non-numeric columns still exist in the dataset.")
}

# Remove low variance features from both training and test data
external_data <- external_data[, !names(external_data) %in% low_variance_features]

# Create the matrix excluding the target column
external_matrix <- as.matrix(external_data[, -which(names(external_data) == target_column)])

# Create DMatrix
dexternal <- xgb.DMatrix(data = external_matrix, label = external_data[[target_column]])

# Make predictions
preds <- predict(final_model, dexternal)

# Convert predictions to binary labels
pred_labels <- ifelse(preds > 0.5, 1, 0)

# Specify that the positive class is 1
confusion_matrix <- confusionMatrix(as.factor(pred_labels), as.factor(external_data$CDAI_response), positive = "1")
print(confusion_matrix)
cat(paste("External test no Optimal Threshold - Confusion Matrix: ", confusion_matrix, "\n\n"), file = file_conn)


# Compute AUC
roc_curve <- roc(external_data$CDAI_response, preds)
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))
cat(paste("External test no Optimal Threshold - AUC: ", auc_value, "\n\n"), file = file_conn)


#OPTIMAL Threshold
# Convert predictions to binary labels using the optimal threshold
pred_labels_optimal <- ifelse(preds > optimal_threshold, 1, 0)

# Specify that the positive class is 1
confusion_matrix_external_optimal <- confusionMatrix(as.factor(pred_labels_optimal), as.factor(external_data$CDAI_response), positive = "1")
print(confusion_matrix_external_optimal)
cat(paste("External test with Optimal Threshold - Confusion Matrix: ", confusion_matrix_external_optimal, "\n\n"), file = file_conn)

f1_score <- confusion_matrix_external_optimal$byClass["F1"]
print(paste("F1 score for Positive Class (class '1') on External Data with Optimal Threshold:", f1_score))
cat(paste("F1 score for Positive Class (class '1') on External Data with Optimal Threshold: ", f1_score, "\n\n"), file = file_conn)

#Now calculate the metrics for negative class
cm_table <- confusion_matrix_external_optimal$table
TN <- cm_table[1, 1]  # True Negatives: correct predictions of class "0"
FP <- cm_table[1, 2]  # False Positives: actual "0" predicted as "1"
FN <- cm_table[2, 1]  # False Negatives: actual "1" predicted as "0"
TP <- cm_table[2, 2]  # True Positives: correct predictions of class "1"
NPV <- TN / (TN + FN)
Specificity <- TN / (TN + FP)
F1_neg <- 2 * (NPV * Specificity) / (NPV + Specificity)
print(paste("F1 Score for Negative Class (class '0') on External Data with Optimal Threshold):", F1_neg))
cat(paste("F1 score for Negative Class (class '0') on External Data with Optimal Threshold: ", F1_neg, "\n\n"), file = file_conn)


# Compute AUC
roc_curve <- roc(external_data$CDAI_response, preds)
auc_value_external_optimal <- auc(roc_curve)
print(paste("AUC on External Data with Optimal Threshold:", auc_value_external_optimal))
cat(paste("External test with Optimal Threshold - AUC: ", auc_value_external_optimal, "\n\n"), file = file_conn)

# Calculate the confidence intervals for the ROC curve
ci_tuned <- ci(auc_value_external_optimal, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
print(paste("95% Confidence Intervals for AUC on External Data with optimum threshold:"))
print(paste("Lower bound:", ci_tuned[1]))
print(paste("Upper bound:", ci_tuned[3]))
cat(paste("95% Confidence Intervals for AUC on External Data with optimum threshold: ", "\n"), file = file_conn)
cat(paste("Lower bound: ", ci_tuned[1], "\n"), file = file_conn)
cat(paste("Upper bound: ", ci_tuned[3], "\n"), file = file_conn)


## ROC plots start
# Let's draw the ROC curve

# Define the full path for the PDF file
pdf_file_xgb_external_opt <- file.path(outputs_path, "roc_curve_plot_xgb_orig_EXTERNAL_opt.pdf")

# Open the PDF device
pdf(pdf_file_xgb_external_opt)


# Create a data frame with sensitivity and specificity values
roc_data <- data.frame(Sensitivity = roc_curve$sensitivities,
                       Specificity = 1 - roc_curve$specificities)


# Plot the ROC curve
plot(roc_data$Specificity, roc_data$Sensitivity, type = "l",
     main = "ROC Curve", xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)", xlim = c(0, 1), ylim = c(0, 1))

# Add a diagonal line representing random classifier (AUC = 0.500)
abline(0, 1, col = "gray", lty = 2)


# Find index of the closest threshold
closest_idx <- which.min(abs(roc_curve$thresholds - optimal_threshold))

# Retrieve the TPR and FPR at the closest threshold
selected_tpr <- roc_curve$sensitivities[closest_idx]
selected_fpr <- 1 - roc_curve$specificities[closest_idx]


# Print the selected TPR and FPR
cat("Selected TPR: ", selected_tpr, "\n", file = file_conn)
cat("Selected FPR: ", selected_fpr, "\n", file = file_conn)



# Calculate ROC curve values for optimal threshold
roc_optimal <- roc(external_data[, outcomeName], pred_labels_optimal)

# Create a data frame with sensitivity and specificity values for the optimal roc
roc_data_optimal <- data.frame(Sensitivity = roc_optimal$sensitivities,
                               Specificity = 1 - roc_optimal$specificities)


##
# Find the closest index for the selected FPR and TPR in the optimal ROC data
closest_spec_idx <- which.min(abs(roc_data_optimal$Specificity - selected_fpr))
closest_sens_idx <- which.min(abs(roc_data_optimal$Sensitivity - selected_tpr))

# Since the indices might differ, ensure using only one index based on a primary criterion, e.g., specificity
selected_spec <- roc_data_optimal$Specificity[closest_spec_idx]
selected_sens <- roc_data_optimal$Sensitivity[closest_spec_idx]

# Add point for the selected TPR and FPR on the ROC curve
points(selected_spec, selected_sens, pch = 20, col = "red")

if (length(selected_spec) > 0 && length(selected_sens) > 0) {
  points(selected_spec, selected_sens, pch = 20, col = "red")
} else {
  cat("No close matches found for selected TPR and FPR on the ROC curve.\n")
}

##


formatted_optimal_threshold <- sprintf("%.3f", optimal_threshold)
formatted_tpr <- sprintf("%.3f", selected_tpr)
formatted_fpr <- sprintf("%.3f", selected_fpr)

text(selected_fpr + 0.10, selected_tpr - 0.10, 
     labels = paste("Threshold =", formatted_optimal_threshold, 
                    "\nTPR =", formatted_tpr, 
                    "\nFPR =", formatted_fpr),
     cex = 1.2, col = "blue", adj = c(0, 0))


cat("roc_data_optimal$Sensitivity: ", roc_data_optimal$Sensitivity, "\n\n", file = file_conn)
cat("roc_data_optimal$Specificity: ", roc_data_optimal$Specificity, "\n\n", file = file_conn)

# Add legend with AUC and 95% CI and F1-score for positive class
# Get the rounded value of F1 score and auc_value
f1_value_str <- sprintf("%.3f", f1_score)
F1_neg_str <- sprintf("%.3f", F1_neg)
auc_value <- round(auc_value_external_optimal, 3)
ci <- ci_tuned 

# Construct the legend text to include AUC, its 95% CI, and the F1 score
legend_text <- c(
  paste("AUC: ", auc_value, " (", round(ci[1], 3), " to ", round(ci[3], 3), ")", sep = ""),
  paste("F1 (positive class): ", f1_value_str),
  paste("F1 (negative class): ", F1_neg_str)
)

# Add the legend to the plot
legend("bottomright", legend = legend_text,
       lty = c(0, 0), col = "black", bty = "n", cex = 1.2)





# Close the PDF device and save the plot
dev.off()

## ROC plots end


cat(paste("End of  evaluation using External ", "\n\n"), file = file_conn)


#  END - EVALUATE USING EXTERNAL DATASET (LA SUISSE)
#################


##############################################
# R1 - BEGIN Sensitivity analysis: baseline CDAI > 10 (EVALUATION ONLY, no retraining)
##############################################

cat("===== Sensitivity analysis: baseline CDAI > 10 (evaluation only) =====\n",
    file = file_conn, append = TRUE)

# sanity check: CDAI0 must exist in the factor datasets (pre-encoding)
if (!"CDAI0" %in% names(test_data_factor)) {
  stop("CDAI0 not found in test_data_factor. Check the baseline CDAI variable name.")
}
if (!"CDAI0" %in% names(external_data_factor)) {
  stop("CDAI0 not found in external_data_factor. Check the baseline CDAI variable name.")
}

stopifnot(nrow(test_data) == nrow(test_data_factor))
stopifnot(nrow(external_data) == nrow(external_data_factor))


# Indices on the factor (human-readable) datasets
idx_test_cdai10 <- !is.na(test_data_factor$CDAI0) & test_data_factor$CDAI0 > 10
idx_ext_cdai10  <- !is.na(external_data_factor$CDAI0) & external_data_factor$CDAI0 > 10


cat(paste0("N test with baseline CDAI0 > 10: ", sum(idx_test_cdai10), "\n"),
    file = file_conn, append = TRUE)
cat(paste0("N external with baseline CDAI0 > 10: ", sum(idx_ext_cdai10), "\n"),
    file = file_conn, append = TRUE)

# Subset the encoded datasets (must match row order of factor datasets, which they do in the pipeline)
test_cdai10     <- test_data[idx_test_cdai10, , drop = FALSE]
external_cdai10 <- external_data[idx_ext_cdai10, , drop = FALSE]

# Predict with the SAME final_model
dtest_cdai10 <- xgb.DMatrix(
  data  = as.matrix(test_cdai10[, -which(names(test_cdai10) == "CDAI_response")]),
  label = test_cdai10$CDAI_response
)
dext_cdai10 <- xgb.DMatrix(
  data  = as.matrix(external_cdai10[, -which(names(external_cdai10) == "CDAI_response")]),
  label = external_cdai10$CDAI_response
)


pred_test_cdai10 <- predict(final_model, dtest_cdai10)
pred_ext_cdai10  <- predict(final_model, dext_cdai10)

# AUC
roc_test_cdai10 <- pROC::roc(test_cdai10$CDAI_response, pred_test_cdai10)
roc_ext_cdai10  <- pROC::roc(external_cdai10$CDAI_response, pred_ext_cdai10)

auc_test_cdai10 <- pROC::auc(roc_test_cdai10)
auc_ext_cdai10  <- pROC::auc(roc_ext_cdai10)


cat(paste0("AUC test (CDAI0 > 10): ", auc_test_cdai10, "\n"),
    file = file_conn, append = TRUE)
cat(paste0("AUC external (CDAI0 > 10): ", auc_ext_cdai10, "\n"),
    file = file_conn, append = TRUE)

# Thresholded metrics using the SAME optimal_threshold from training
pred_test_lbl <- ifelse(pred_test_cdai10 > optimal_threshold, 1, 0)
pred_ext_lbl  <- ifelse(pred_ext_cdai10  > optimal_threshold, 1, 0)

cm_test_cdai10 <- caret::confusionMatrix(
  as.factor(pred_test_lbl),
  as.factor(test_cdai10$CDAI_response),
  positive = "1"
)
cm_ext_cdai10 <- caret::confusionMatrix(
  as.factor(pred_ext_lbl),
  as.factor(external_cdai10$CDAI_response),
  positive = "1"
)

cat("Confusion matrix test (CDAI0 > 10):\n",
    paste(capture.output(cm_test_cdai10), collapse = "\n"),
    "\n", file = file_conn, append = TRUE)

cat("Confusion matrix external (CDAI0 > 10):\n",
    paste(capture.output(cm_ext_cdai10), collapse = "\n"),
    "\n", file = file_conn, append = TRUE)


## SHAP for CDAI > 10 subset (same final model)

# sanity checks
stopifnot(nrow(train_data) == nrow(train_data_factor))

idx_train_cdai10 <- !is.na(train_data_factor$CDAI0) & train_data_factor$CDAI0 > 10

cat(paste0("N train with baseline CDAI0 > 10: ",
           sum(idx_train_cdai10), "\n"),
    file = file_conn, append = TRUE)

train_cdai10 <- train_data[idx_train_cdai10, , drop = FALSE]

#ensure same feature space
train_cdai10 <- train_cdai10[, !names(train_cdai10) %in% low_variance_features, drop = FALSE]

#compute SHAP using final_model
library(SHAPforxgboost)

X_train_cdai10 <- as.matrix(
  train_cdai10[, setdiff(names(train_cdai10), target_column)]
)

shap_cdai10 <- shap.prep(
  xgb_model = final_model,
  X_train   = X_train_cdai10
)

#save + plot (top 10)
if (EXPORT_DERIVED_OUTPUTS) write.csv(
  shap_cdai10,
  file = fp_out("SHAP_values_train_CDAIgt10.csv"),
  row.names = FALSE
)

pdf(fp_out("SHAP_Summary_Plot_train_CDAIgt10.pdf"), width = 8, height = 6)
shap.plot.summary(shap_cdai10)
dev.off()

# Top 10 bar plot (mean |SHAP|)
mean_shap_cdai10 <- shap_cdai10 |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_abs_shap = mean(abs(value), na.rm = TRUE)) |>
  dplyr::arrange(desc(mean_abs_shap))

if (EXPORT_DERIVED_OUTPUTS) write.csv(
  mean_shap_cdai10,
  file = fp_out("SHAP_meanAbs_train_CDAIgt10.csv"),
  row.names = FALSE
)



################################################
# END - Sensitivity analysis: baseline CDAI > 10
################################################

###############################################################
# R1 - Subgroup analysis by baseline b/ts-DMARD therapy class
#      (MoA / treatment_group)
#      Evaluation ONLY, no retraining
#      - Applies the final trained XGBoost model to each subgroup
#      - Uses a single fixed threshold derived from TRAIN ROC
#      - Outputs subgroup performance tables + subgroup SHAP top10
###############################################################

# --- safety checks
if (!("treatment_group" %in% names(train_data_factor))) stop("train_data_factor missing treatment_group")
if (!("treatment_group" %in% names(test_data_factor)))  stop("test_data_factor missing treatment_group")
if (!("treatment_group" %in% names(external_data_factor))) stop("external_data_factor missing treatment_group")

# Ensure outcome is 0/1 in these datasets (should already be, but keep it explicit)
TARGET_COLUMN <- target_column  # should be "CDAI_response"

# helper: compute subgroup metrics at a fixed threshold
subgroup_metrics_moa <- function(y, p, group, thr) {
  df <- data.frame(y = y, p = p, group = group) |>
    dplyr::filter(!is.na(group))
  
  df |>
    dplyr::group_by(group) |>
    dplyr::do({
      d <- .
      auc_val <- if (length(unique(d$y)) == 2) {
        as.numeric(pROC::auc(pROC::roc(d$y, d$p, quiet = TRUE)))
      } else {
        NA_real_
      }
      
      pred <- ifelse(d$p > thr, 1, 0)
      
      cm <- caret::confusionMatrix(
        factor(pred, levels = c(0, 1)),
        factor(d$y, levels = c(0, 1)),
        positive = "1"
      )
      
      data.frame(
        n                 = nrow(d),
        auc               = auc_val,
        sensitivity       = unname(cm$byClass["Sensitivity"]),
        specificity       = unname(cm$byClass["Specificity"]),
        ppv               = unname(cm$byClass["Pos Pred Value"]),
        npv               = unname(cm$byClass["Neg Pred Value"]),
        balanced_accuracy = unname(cm$byClass["Balanced Accuracy"])
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::rename(treatment_group = group)
}

# --- build aligned feature matrix (strictly the same feature space as the final model)
# Note: low_variance_features were removed earlier, keep consistent here.
train_moa <- train_data[, !names(train_data) %in% low_variance_features, drop = FALSE]
test_moa  <- test_data[,  !names(test_data)  %in% low_variance_features, drop = FALSE]
ext_moa   <- external_data[, !names(external_data) %in% low_variance_features, drop = FALSE]

model_features_moa <- setdiff(names(train_moa), TARGET_COLUMN)

# sanity: check identical feature sets
stopifnot(identical(model_features_moa, setdiff(names(test_moa), TARGET_COLUMN)))
stopifnot(identical(model_features_moa, setdiff(names(ext_moa),  TARGET_COLUMN)))

dtrain_moa <- xgb.DMatrix(as.matrix(train_moa[, model_features_moa, drop = FALSE]), label = train_moa[[TARGET_COLUMN]])
dtest_moa  <- xgb.DMatrix(as.matrix(test_moa[,  model_features_moa, drop = FALSE]), label = test_moa[[TARGET_COLUMN]])
dext_moa   <- xgb.DMatrix(as.matrix(ext_moa[,   model_features_moa, drop = FALSE]), label = ext_moa[[TARGET_COLUMN]])

pred_train_moa <- predict(final_model, dtrain_moa)
pred_test_moa  <- predict(final_model, dtest_moa)
pred_ext_moa   <- predict(final_model, dext_moa)

# --- threshold from TRAIN ROC (consistent with Methods)
roc_train_moa <- pROC::roc(train_moa[[TARGET_COLUMN]], pred_train_moa, quiet = TRUE)
optimal_threshold_moa <- pROC::coords(roc_train_moa, "best", ret = "threshold")$threshold

cat(paste0("

Subgroup analysis by treatment_group (MoA)
",
           "Threshold from train ROC (best/Youden): ", round(optimal_threshold_moa, 6), "
"),
    file = file_conn, append = TRUE)

# --- subgroup metrics: internal TEST + external
res_test_moa <- subgroup_metrics_moa(
  y     = test_moa[[TARGET_COLUMN]],
  p     = pred_test_moa,
  group = test_data_factor$treatment_group,
  thr   = optimal_threshold_moa
)

res_ext_moa <- subgroup_metrics_moa(
  y     = ext_moa[[TARGET_COLUMN]],
  p     = pred_ext_moa,
  group = external_data_factor$treatment_group,
  thr   = optimal_threshold_moa
)

# save results
if (EXPORT_DERIVED_OUTPUTS) write.csv(res_test_moa, fp_out("subgroup_metrics_TEST_by_treatment_group.csv"), row.names = FALSE)
if (EXPORT_DERIVED_OUTPUTS) write.csv(res_ext_moa,  fp_out("subgroup_metrics_EXTERNAL_by_treatment_group.csv"), row.names = FALSE)
cat("Wrote subgroup metrics tables:
",
    " - Outputs/subgroup_metrics_TEST_by_treatment_group.csv
",
    " - Outputs/subgroup_metrics_EXTERNAL_by_treatment_group.csv
",
    file = file_conn, append = TRUE)

# also print to log
cat("
Internal TEST by treatment_group:
",
    paste(capture.output(res_test_moa), collapse = "
"),
    "

External by treatment_group:
",
    paste(capture.output(res_ext_moa), collapse = "
"),
    "
", file = file_conn, append = TRUE)

# --- SHAP top10 by group (using already-saved training SHAP file)
# Requires Outputs/shap_values_finetuned_optimal_threshold.csv (created earlier in this script)
if (!file.exists(fp_out("shap_values_finetuned_optimal_threshold.csv"))) {
  warning("SHAP file not found: Outputs/shap_values_finetuned_optimal_threshold.csv. Skipping subgroup SHAP.")
} else {
  
  shap_values_moa <- readr::read_csv(fp_out("shap_values_finetuned_optimal_threshold.csv"), show_col_types = FALSE)
  
  if (!all(c("ID", "variable", "value") %in% names(shap_values_moa))) {
    stop("SHAP file missing required columns: ID, variable, value")
  }
  
  id_map_moa <- tibble::tibble(
    ID = 1:nrow(train_data_factor),
    treatment_group = train_data_factor$treatment_group
  )
  
  shap_top10_by_group_moa <- shap_values_moa |>
    dplyr::left_join(id_map_moa, by = "ID") |>
    dplyr::filter(!is.na(treatment_group)) |>
    # 1. REMOVE treatment group dummy variables
    dplyr::filter(!stringr::str_detect(variable, "^treatment_group\\.")) |>
    
    # 2. Compute mean absolute SHAP per variable per group
    dplyr::group_by(treatment_group, variable) |>
    dplyr::summarise(mean_abs_shap = mean(abs(value), na.rm = TRUE), .groups = "drop") |>
    
    # 3. Re-rank variables WITHIN each treatment group
    dplyr::group_by(treatment_group) |>
    dplyr::arrange(desc(mean_abs_shap)) |>
    
    # 4. Take the top 10 AFTER exclusion
    dplyr::slice_head(n = 10) |>
    dplyr::ungroup()
  
  if (EXPORT_DERIVED_OUTPUTS) write.csv(shap_top10_by_group_moa,
            fp_out("SHAP_top10_by_treatment_group_finetuned.csv"),
            row.names = FALSE)
  
  cat("
Wrote subgroup SHAP table:
",
      " - Outputs/SHAP_top10_by_treatment_group_finetuned.csv
",
      file = file_conn, append = TRUE)
}

###############################################################
# END - Subgroup analysis by baseline therapy class (MoA)
###############################################################


close(file_conn)