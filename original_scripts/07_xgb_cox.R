# 07_xgb_cox.R
# XGBoost model with Cox proportional hazards objective:
# training, C-index evaluation, time-dependent AUC calculation,
# and feature importance extraction.
#
# Inputs:  outputs/metabric_cleaned.csv
# Outputs: outputs/XGB_COX_Cindex.csv,
#          outputs/XGB_COX_AUC.csv,
#          outputs/XGB_COX_feature_importance.csv

# Setup: Working Directory

proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"

if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)

setwd(proj_dir); if (!dir.exists("outputs")) dir.create("outputs")
cat("Working in:", getwd(), "\n")


# Load Libraries

library(xgboost)    # Gradient boosting framework
library(survival)   # Survival analysis utilities
library(timeROC)    # Time-dependent ROC and AUC computation
library(Matrix)     # Sparse matrix operations
library(readr)      # Reading and writing CSV files

# Load and preprocess data
df0 <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = FALSE)

# Variables to be used
vars <- c("event","time","age_at_diagnosis","tumor_size","lymph_nodes_examined_positive",
          "er_status","pr_status","her2_status")

# Retain complete cases for variables of interest
df <- df0[complete.cases(df0[, vars]), vars]

# One-hot encode categorical variables and prepare features
X <- model.matrix(~ age_at_diagnosis + tumor_size + lymph_nodes_examined_positive +
                    er_status + pr_status + her2_status, data=df)[, -1, drop=FALSE]

# Prepare survival response variable for Cox objective
# Positive times with event =1, negative times for censored observations
y <- ifelse(df$event==1, df$time, -df$time)

# Create DMatrix for xgboost training
dtrain <- xgb.DMatrix(data = X, label = y)

# Define XGBoost parameters for Cox proportional hazards model
params <- list(objective="survival:cox", eval_metric="cox-nloglik",
               eta=0.05, max_depth=3, subsample=0.8, colsample_bytree=0.8,
               tree_method="hist")

# Train model
set.seed(42)
xgb_cox <- xgb.train(params=params, data=dtrain, nrounds=800, verbose=0)

# Predict risk scores for training data
risk <- as.numeric(predict(xgb_cox, dtrain))

# Calculate and save Concordance Index (C-index)
cobj <- survival::concordance(Surv(df$time, df$event) ~ risk, data=df)
cindex_val <- cobj$concordance; cindex_se <- sqrt(cobj$var)
write_csv(tibble(Model="XGB-COX", C_index=cindex_val, SE=cindex_se),
          "outputs/XGB_COX_Cindex.csv")
cat(sprintf("XGB-COX C-index = %.3f (SE=%.3f)\n", cindex_val, cindex_se))

# Compute time-dependent AUCs at specified months
eval_times <- c(12,36,60)
roc_xgb <- timeROC(T=df$time, delta=df$event, marker=risk, cause=1, times=eval_times, iid=FALSE)
auc_tbl <- tibble(Time_in_Months = roc_xgb$times, AUC = as.numeric(roc_xgb$AUC))
write_csv(auc_tbl, "outputs/XGB_COX_AUC.csv"); print(auc_tbl)

# Extract and save feature importance from XGBoost model
feat_names <- colnames(X)
imp <- xgb.importance(model=xgb_cox, feature_names=feat_names)
write_csv(imp, "outputs/XGB_COX_feature_importance.csv"); print(head(imp, 10))

cat("XGB-COX results saved in outputs/\n")
