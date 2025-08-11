# Load libraries
library(caret)
library(pROC)
library(mice)
library(dplyr)

# Set working directory
setwd("~/Projects/32. ML for JAKPOT/CDAI_2.8_Short")

# Load datasets
train_data <- read.csv("Outputs - Fine_tuning/train_data.csv")
test_data <- read.csv("Outputs - Fine_tuning/test_data.csv")
external_data <- read.csv("Outputs - Fine_tuning/external_data.csv")

# Define predictors and outcome
predictors <- c("PGA",
                "TJC28",
                "Prev_btsDMARD",
                "HAQ0",
                "PhGA",
                "Age",
                "CDAI0",
                "CRP",
                "Pain",
                "AntiCCPevernever")
external_predictors <- setdiff(predictors, "Pain")
outcome <- "CDAI_response"
all_vars <- c(predictors, outcome)

# Subset and format
train_data <- train_data[, all_vars]
test_data <- test_data[, all_vars]
external_data <- external_data[, c(external_predictors, outcome)]

train_data$CDAI_response <- as.factor(train_data$CDAI_response)
test_data$CDAI_response <- as.factor(test_data$CDAI_response)
external_data$CDAI_response <- as.factor(external_data$CDAI_response)

# === Step 1: Multiple Imputation on Training Data ===
imp_train <- mice(train_data, m = 10, method = 'pmm', seed = 42)

# === Step 2: Initialize storage for AUC results ===
auc_train_list <- numeric(10)
auc_test_list <- numeric(10)
auc_external_list <- numeric(10)

# === Step 3: Loop over imputations ===
for (i in 1:10) {
  train_i <- complete(imp_train, i)
  
  # Full model (with Pain) for internal test
  model_full <- glm(CDAI_response ~ ., data = train_i, family = binomial)
  
  # Predict + AUC on training set
  probs_train <- predict(model_full, newdata = train_i, type = "response")
  auc_train_list[i] <- auc(roc(train_i$CDAI_response, probs_train))
  
  # Remove "Pain" for external model training
  train_i_no_pain <- train_i[, !(names(train_i) %in% "Pain")]
  model_no_pain <- glm(CDAI_response ~ ., data = train_i_no_pain, family = binomial)
  
  # Impute test and external using same logic
  test_i <- complete(mice(test_data, m = 1, method = 'pmm', seed = 100 + i), 1)
  external_i <- complete(mice(external_data, m = 1, method = 'pmm', seed = 200 + i), 1)
  
  # Predict + AUC on test set (Pain included)
  probs_test <- predict(model_full, newdata = test_i, type = "response")
  auc_test_list[i] <- auc(roc(test_i$CDAI_response, probs_test))
  
  # Predict + AUC on external set (Pain excluded)
  probs_external <- predict(model_no_pain, newdata = external_i, type = "response")
  auc_external_list[i] <- auc(roc(external_i$CDAI_response, probs_external))
}


# === Step 4: Pool and report AUC results ===
mean_auc_train <- mean(auc_train_list)
sd_auc_train <- sd(auc_train_list)
cat("\nPooled Logistic Regression Performance (Training Set):\n")
cat(sprintf("AUC (mean ± SD): %.3f ± %.3f\n", mean_auc_train, sd_auc_train))

mean_auc_test <- mean(auc_test_list)
sd_auc_test <- sd(auc_test_list)
mean_auc_external <- mean(auc_external_list)
sd_auc_external <- sd(auc_external_list)

cat("\nPooled Logistic Regression Performance (Test Set):\n")
cat(sprintf("AUC (mean ± SD): %.3f ± %.3f\n", mean_auc_test, sd_auc_test))

cat("\nPooled Logistic Regression Performance (External Validation):\n")
cat(sprintf("AUC (mean ± SD): %.3f ± %.3f\n", mean_auc_external, sd_auc_external))
