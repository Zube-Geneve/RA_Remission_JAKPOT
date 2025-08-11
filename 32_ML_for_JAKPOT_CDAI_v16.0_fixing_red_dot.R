#Version history
#
# v1.2  Initial creation from discont_TNF_v1.2
# v1.3  Updates based on feedback from the presentation on 20.June.2024
# v1.4  Make treatment_group based on ttt
# v1.4.2 Change treatment_group
# v1.4.3 Missing values imputed as single value
# v1.4.4 Reserve Suisse for external test
# V1.4.5 Not to impute and let XGBOOST to handle the missing values
# V2.0   Tidy up with reduced model etc.
# v3.0   Further tuning hyperparameters
# v4.0 Using already selected optimal hyperparameters
# v5.0 Add Turkey
# v6.0 Fixing the issue about including baseline CDAI when check the CDAI follow-up. 
# v7.0 Exclude DE, SP, TR and UK because they do not have CDAI values due to not to have PGA and PhGA (CDAI=TJC+SJC+PGA+PhGA)
#      Remove adding previous year dataset (i.e., TR, IL etc. as they have no CDAI)
#  NOTE: This code is now using no remission (CDAI > 2.8 ) as positive class. Also the fine tuning needs to be re-done.
# v8.0 Find the new tuning parameters
#
# v9.0 Applying the changes for ROC plots to test and external data
#
# v10.0 Selecting only the patients with at least one follow-up visit between 6 to 24 months after treatment start
#
# v11.0 Selecting only the patients with at least one CDAI info during follow-up visit between 6 to 24 months after treatment start
#
# v12.0 Same as v11.0 but for the selecting the patients I tried a different logic which is
#       - Participants are included if they have at least one valid CDAI measurement during the 6–24 month follow-up period
#       - If any CDAI value ≤ 2.8 during the valid follow-up visits, the patient is classified as achieving remission 
#       - If all available CDAI values are > 2.8, the patient is classified as not achieving remission
#       - If no valid CDAI is recorded during the follow-up period (6–24 months), those participants are excluded from the final dataset 
#       - 
# v13.0 Applying SHAP as well as using the dummies.
# v14.0 Treating treatment_group as one-hot encode. Also fixing the Table 1 to report factor variables not dummy.
# v15.0 Aligning remission as positive class in all dataset analyses
# v16.0 Recalculate BMI if Height and Weight are present
#
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

# set working directory
setwd("C:/Users/zubey/OneDrive/Documents/Projects/32. ML for JAKPOT/CDAI")


# load the previous data set for the missing countries in the new dataset
#load('~/PostDoc/JAKPOT/2023-07-12/clean_data_2023-07-12.RData')

#write.csv(jakpot.data_eff, file="Outputs/jakpot.data_eff_prev.csv")

# Filter records for "TR", "IL", "RU" from the previous datasets
#additional_countries_data <- filter(jakpot.data_eff, Code %in% c("TR", "IL"))
#additional_countries_data_baseline <- filter(jakpot.data_effbaseline, Code %in% c("TR", "IL"))
#additional_countries_data_ae <- filter(jakpot.data_ae, Code %in% c("TR", "IL"))


# load the current dataset
load('~/PostDoc/JAKPOT/clean_data_2023-12-12.RData')

# Combine the new data set with the filtered records from the previous dataset
#jakpot.data_eff1 <- rbind(jakpot.data_eff, additional_countries_data)
#jakpot.data_eff1_baseline <- rbind(jakpot.data_effbaseline, additional_countries_data_baseline)
#jakpot.data_eff1_ae <- rbind(jakpot.data_eff_ae, additional_countries_data_ae)
#Above RBIND is giving error: " Item 2 has 255 columns, inconsistent with item 1 which has 256 columns".

# Get the column names of both datasets
#columns_dataset1 <- colnames(jakpot.data_eff)
#columns_dataset2 <- colnames(additional_countries_data)

# Find columns that are in dataset1 but not in dataset2
#missing_in_dataset2 <- setdiff(columns_dataset1, columns_dataset2)

# Find columns that are in dataset2 but not in dataset1
#missing_in_dataset1 <- setdiff(columns_dataset2, columns_dataset1)

# Print the missing columns
#print("Columns missing in dataset2:")
#print(missing_in_dataset2)

#print("Columns missing in dataset1:")
#print(missing_in_dataset1)

# I think we just combine them anyway 
# Combine datasets with possibly different column structures
#jakpot.data_eff1 <- dplyr::bind_rows(jakpot.data_eff, additional_countries_data)
#jakpot.data_eff1_baseline <- dplyr::bind_rows(jakpot.data_effbaseline, additional_countries_data_baseline)
#jakpot.data_eff1_ae <- dplyr::bind_rows(jakpot.data_ae, additional_countries_data_ae)




jakpot.data_eff1 <- jakpot.data_eff
jakpot.data_eff1_baseline <- jakpot.data_effbaseline
jakpot.data_eff1_ae <- jakpot.data_ae

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


#Calculate new treatment_duration_calculation because Denis assumed the last visit date as 2023-11-10 for all countries which is not correct for some 
jakpot.data_eff1_baseline$extraction_date_calc <- jakpot.data_eff1_baseline$extraction_date

# Setting italian extraction date a bit later, because we received the data late
jakpot.data_eff1_baseline[Code %in% "IT", extraction_date_calc := as.Date("2023-12-01")]

# Setting UK and DE extraction date a bit earlier because we don't have the most recent data
jakpot.data_eff1_baseline[Code %in% "DE", extraction_date_calc := as.Date("2023-07-01")]
jakpot.data_eff1_baseline[Code %in% "UK", extraction_date_calc := as.Date("2023-10-14")]

# Setting TR, RU, and FI extraction date much earlier because we don't have the most recent data (for RU the latest stop date was after the visit date so I used that. for others last visit date)
jakpot.data_eff1_baseline[Code %in% "TR", extraction_date_calc := as.Date("2020-11-20")]
jakpot.data_eff1_baseline[Code %in% "RU", extraction_date_calc := as.Date("2021-02-02")]
jakpot.data_eff1_baseline[Code %in% "FI", extraction_date_calc := as.Date("2022-12-30")]

# Focusing only on data before 2021 for NO, because we have AE only until dec 2020
jakpot.data_eff1_baseline[Code %in% "NO", extraction_date_calc := as.Date("2021-01-01")]
jakpot.data_eff1_baseline <- jakpot.data_eff1_baseline[!(Code %in% "NO") | (Code %in% "NO" & ttt_startDate <= extraction_date_calc), ]


# Calculate treatment_duration
#jakpot.data_eff1_baseline[, treatment_duration_calc := ifelse(
#  is.na(ttt_stopDate) & is.na(lostFUDate),
#  (extraction_date_calc - ttt_startDate) / 365.25,
#  ifelse(
#    is.na(ttt_stopDate) & !is.na(lostFUDate),
#    (lostFUDate - ttt_startDate) / 365.25,
#    (ttt_stopDate - ttt_startDate) / 365.25
#  )),
#  by = ID_ttt
#]

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
        # treatment_IL17 ~ "IL17",
        # treatment_IL23 ~ "IL23",
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

write.csv(jakpot.data_eff1, file="Outputs/jakpot_data_eff1.csv")
write.csv(jakpot.data_eff1_baseline, file="Outputs/jakpot_data_eff1_baseline.csv")
write.csv(remission_cdai, file="Outputs/remission_cdai.csv")
#write.csv(jakpot.data_eff1_ae, file="Outputs/jakpot_data_eff1_ae.csv")

combined_dataset <- arrange(jakpot.data_eff1_baseline,ID,Visit_date)


### REMOVE the rows with missing outcomes
combined_dataset_original_count <- nrow(combined_dataset)

# Remove rows with missing values or NaNs in the outcome
# Note that in fact there is no missing outcome because we have made the worst case assumption i.e., if missing there is no outcome
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
write.csv(continuous_to_categorical, "Outputs/continuous_to_categorical.csv", row.names = FALSE)

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
write.csv(combined_dataset_conv, file="Outputs/combined_dataset_conv.csv")

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
write.csv(combined_dataset_conv, file="Outputs/combined_dataset_conv_apres_exc.csv", row.names = FALSE)


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


# Print dimensions to verify the number of rows and columns
cat("External Data Dimensions: ", dim(external_data), "\n")
cat("Internal Data Dimensions: ", dim(internal_data), "\n")


write.csv(external_data, file="Outputs/external_data.csv")
write.csv(internal_data, file="Outputs/internal_data.csv")


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


write.csv(train_data, file="Outputs/train_data.csv")
write.csv(test_data, file="Outputs/test_data.csv")

# Percentage and number of outcome in training dataset
train_data_outcome_percent <- prop.table(table(train_data$CDAI_response)) * 100
print(table(train_data$CDAI_response))
print(train_data_outcome_percent)

# Percentage and number of outcome in test dataset
test_data_outcome_percent <- prop.table(table(test_data$CDAI_response)) * 100
print(table(test_data$CDAI_response))
print(test_data_outcome_percent)


# Percentage and number of outcome in test dataset
external_data_outcome_percent <- prop.table(table(external_data$CDAI_response)) * 100
print(table(external_data$CDAI_response))
print(external_data_outcome_percent)


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
  write.csv(na_report, file = paste0("Outputs/missingness_", dataset_name, ".csv"), row.names = FALSE)
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
writeLines(captured_output, "Outputs/table1.txt")
####


# 4. Model training

# Open a file connection for writing the outputs
outputs_path <- "C:/Users/zubey/OneDrive/Documents/Projects/32. ML for JAKPOT/CDAI/Outputs/"
tool_path <- "C:/Users/zubey/OneDrive/Documents/Projects/32. ML for JAKPOT/CDAI/Tool/"

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


# Identify low variance features in the training data
low_variance_features <- nearZeroVar(train_data, saveMetrics = TRUE)
low_variance_features <- rownames(low_variance_features[low_variance_features$nzv, ])

# Remove low variance features from both training and test data
train_data <- train_data[, !names(train_data) %in% low_variance_features]


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




##Importance matrix
#importance_matrix <- xgb.importance(model = model)
#print(importance_matrix)
#cat(paste("Training - Importance Matrix: ", importance_matrix, "\n\n"), file = file_conn)
#
## Plot feature importance
#xgb.plot.importance(importance_matrix)
#
## Plot feature importance
#importance_plot <- xgb.plot.importance(importance_matrix)

## Specify the path and file name for the plot image
#plot_file_name <- paste("Training_Feature_Importance_Plot", Sys.Date(), "png", sep = ".")
#plot_full_path <- file.path(outputs_path, plot_file_name)

## Open a graphics device
#png(filename = plot_full_path, width = 10, height = 8, units = "in", res = 300)

## Plot feature importance
#xgb.plot.importance(importance_matrix)

## Close the device
##dev.off()


# Write a reference to the plot file in the text file
#cat(paste("Training - Importance Matrix plot saved at: ", plot_full_path, "\n\n"), file = file_conn)



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
write.csv(shap_values, file = "Outputs/shap_values_training_feature.csv", row.names = FALSE)

# Plot SHAP summary to PDF
pdf("Outputs/SHAP_Summary_Plot_training_feature.pdf")
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
pdf("Outputs/SHAP_Summary_Plot_training_feature_top110.pdf", width = 8, height = 6)
shap.plot.summary(shap_top10)

dev.off()

cat(paste("Training - Shap summary plot saved at: ", "\n\n"), file = file_conn)


#Save the model
xgb.save(model, "Models/xgb_model_full_no_finetune.bin")

cat(paste("End of Training. ", "\n\n"), file = file_conn)


# END TRAINING METHOD - DIRECT USE OF XGBOOST
#################

#################
# FINE-TUNING HYPERPARAMETERS


#set seed for fine-tuning
set.seed(7204)

cat(paste("Starting Training fine-tuning : ", "\n\n"), file = file_conn)


#START of the comment out - Because when I do rerun I don't want to spent 24 hours for hypertuning
##
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
#END of the comment out 


#This is the best parameters
#> print(best_params)
#eta       max_depth gamma colsample_bytree min_child_weight subsample nrounds
# 0.05         9       0       0.5                5            1        150

# Hard-coded best parameters
#best_eta <- 0.05
#best_max_depth <- 9
#best_gamma <- 0
#best_colsample_bytree <- 0.5
#best_min_child_weight <- 5
#best_subsample <- 1
#best_nrounds <- 150

# Update the parameter list with the best parameters
#params <- list(
#  booster = "gbtree",
#  objective = "binary:logistic",
#  eta = best_eta,
#  max_depth = best_max_depth,
#  gamma = best_gamma,
#  colsample_bytree = best_colsample_bytree,
#  min_child_weight = best_min_child_weight,
#  subsample = best_subsample
#)
#
cat(paste("Training fine-tuned - Best parameters : ", params, best_nrounds, "\n\n"), file = file_conn)


# Train the final model with the best number of rounds
final_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds
)


# Save the model for later use
xgb.save(final_model, "Models/final_xgb_model.bin")


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


#### Importance matrix 
#importance_matrix <- xgb.importance(model = final_model)
#print(importance_matrix)
#cat(paste("Training fine-tuned with Optimal Threshold - Importance Matrix: ", importance_matrix, "\n\n"), file = file_conn)
##
#
#
## Plot feature importance
#xgb.plot.importance(importance_matrix)
#
#
# Plot feature importance
#importance_plot <- xgb.plot.importance(importance_matrix)

## Specify the path and file name for the plot image
#plot_file_name <- paste("Training_Feature_Importance_Plot_finetuned_optimal_threshold", Sys.Date(), "png", sep = ".")
#plot_full_path <- file.path(outputs_path, plot_file_name)


# Open a graphics device
#png(filename = plot_full_path, width = 10, height = 8, units = "in", res = 300)

## Plot feature importance
#xgb.plot.importance(importance_matrix)

# Close the device
#dev.off()

# Write a reference to the plot file in the text file
#cat(paste("Training fine-tuned with Optimal Threshold - Importance Matrix plot saved at: ", plot_full_path, "\n\n"), file = file_conn)


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
write.csv(shap_values, "Outputs/shap_values_finetuned_optimal_threshold.csv", row.names = FALSE)

# Plot SHAP summary to PDF
pdf("Outputs/SHAP_Summary_Plot_finetuned_optimal_threshold.pdf")
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
pdf("Outputs/SHAP_Summary_Plot_finetuned_optimal_threshold_top10.pdf", width = 8, height = 6)
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


# Add point for the selected TPR and FPR on the ROC curve
#points(roc_data_optimal$Specificity[roc_data_optimal$Specificity == selected_fpr],
#       roc_data_optimal$Sensitivity[roc_data_optimal$Sensitivity == selected_tpr],
#       pch = 20, col = "red")

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
#  EVALUATE USING EXTERNAL DATASET (LA SUISSE)

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


close(file_conn)
