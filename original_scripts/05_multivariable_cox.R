# 05_multivariable_cox.R
# Multivariable Cox Proportional Hazards Model for breast cancer survival
# Includes model fitting, proportional hazards test, and model evaluation
#
# Inputs:  outputs/metabric_cleaned.csv
# Outputs: outputs/cox_multivariable.csv,
#          outputs/cox_global_test.txt,
#          outputs/CPH_Cindex.csv,
#          outputs/CPH_AUC.csv,
#          outputs/CPH_survival_rates.csv

# Setup: Working Directory
proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"

if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)
setwd(proj_dir)

if (!dir.exists("outputs")) dir.create("outputs")

cat("Working in:", getwd(), "\n")

# Load Libraries
library(survival)   # Survival modeling and analysis
library(survminer)  # Visualization for survival analysis
library(readr)      # CSV read/write utilities
library(dplyr)      # Data handling
library(timeROC)    # Time-dependent ROC and AUC

# Load cleaned dataset
df <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = TRUE)

# Convert survival columns to numeric for analysis
df$event <- as.numeric(df$event)
df$time  <- as.numeric(df$time)


# Step 1: Prepare dataset for multivariable Cox model

# Select relevant variables and remove rows with missing data
df_model <- na.omit(df[, c("time","event","age_at_diagnosis",
                           "tumor_size","lymph_nodes_examined_positive",
                           "er_status","pr_status","her2_status")])

cat("Analysis dataset created with", nrow(df_model), "patients\n")


# Step 2: Fit Multivariable Cox Proportional Hazards Model

fit <- coxph(Surv(time, event) ~ age_at_diagnosis + tumor_size +
               lymph_nodes_examined_positive + er_status +
               pr_status + her2_status, data = df_model)

summary_fit <- summary(fit)

# Extract and format coefficient results
out <- data.frame(
  Variable = rownames(summary_fit$coefficients),
  B        = summary_fit$coefficients[, "coef"],
  SE       = summary_fit$coefficients[, "se(coef)"],
  z        = summary_fit$coefficients[, "z"],
  p        = summary_fit$coefficients[, "Pr(>|z|)"],
  HR       = summary_fit$coefficients[, "exp(coef)"],
  HR_LCL   = summary_fit$conf.int[, "lower .95"],
  HR_UCL   = summary_fit$conf.int[, "upper .95"]
)

print(out)
write_csv(out, "outputs/cox_multivariable.csv")
cat("Saved: outputs/cox_multivariable.csv\n")

# Optional: Global Proportional Hazards Test

global_test <- cox.zph(fit)
print(global_test)
sink("outputs/cox_global_test.txt"); print(global_test); sink()

# Step 3: Model Evaluation (C-index, AUCs, survival rates)


# Calculate and save Concordance Index (C-index)
df_model$risk <- predict(fit, type="risk")
cobj <- concordance(Surv(time, event) ~ risk, data=df_model)
cindex_val <- cobj$concordance
cindex_se  <- sqrt(cobj$var)
cindex_lcl <- cindex_val - 1.96*cindex_se
cindex_ucl <- cindex_val + 1.96*cindex_se

write_csv(data.frame(Model="CPH", C_index=cindex_val, SE=cindex_se,
                     LCL=cindex_lcl, UCL=cindex_ucl),
          "outputs/CPH_Cindex.csv")
cat(sprintf("CPH C-index = %.3f (SE=%.3f) [%.3f, %.3f]\n",
            cindex_val, cindex_se, cindex_lcl, cindex_ucl))

# Calculate and save time-dependent AUCs at 12, 36, 60 months
eval_times <- c(12, 36, 60)
roc_cph <- timeROC(T=df_model$time,
                   delta=df_model$event,
                   marker=df_model$risk,
                   cause=1,
                   times=eval_times,
                   iid=FALSE)

auc_tbl <- data.frame(Time_in_Months=roc_cph$times,
                      AUC=as.numeric(roc_cph$AUC))
print(auc_tbl)
write_csv(auc_tbl, "outputs/CPH_AUC.csv")

# Estimate and save survival probabilities at 5 and 10 years
fit_km <- survfit(fit, newdata=df_model)
tgrid <- fit_km$time
idx5  <- which.min(abs(tgrid - 60))
idx10 <- which.min(abs(tgrid - 120))
surv5  <- summary(fit_km, times=60)$surv
surv10 <- summary(fit_km, times=120)$surv
sr <- data.frame(Model="CPH",
                 Surv5y=mean(surv5, na.rm=TRUE),
                 Surv10y=mean(surv10, na.rm=TRUE))
print(sr)
write_csv(sr, "outputs/CPH_survival_rates.csv")

cat("Multivariable CPH results saved in outputs/\n")
