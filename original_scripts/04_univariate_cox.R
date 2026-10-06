# 04_univariate_cox.R
# Univariate Cox Proportional Hazards models for continuous variables:
# Age at diagnosis, Tumor size, and Positive lymph nodes count.
#
# Inputs:  outputs/metabric_cleaned.csv
# Outputs: outputs/cox_univariate.csv

# Setup: Working Directory
proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"

# Create project directory (and parent directories) if missing
if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)

setwd(proj_dir)

# Ensure output folder exists
if (!dir.exists("outputs")) dir.create("outputs")
cat("Working in:", getwd(), "\n")

# Load Libraries

library(survival)      # Survival analysis functions
library(dplyr)         # Data manipulation
library(readr)         # Reading/writing CSV

# Load cleaned dataset

df <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = TRUE)

run_uni <- function(var) {
  # Create formula for Cox PH model
  f <- as.formula(paste0("Surv(time, event) ~ ", var))
  
  # Fit Cox model
  fit <- coxph(f, data = df)
  
  # Extract summary
  s <- summary(fit)
 
  # Return tidy tibble with estimates and statistics
   tibble(
    Variable = var,
    B = unname(s$coefficients[,"coef"]),
    SE = unname(s$coefficients[,"se(coef)"]),
    z = unname(s$coefficients[,"z"]),
    p = unname(s$coefficients[,"Pr(>|z|)"]),
    HR = unname(s$coefficients[,"exp(coef)"]),
    HR_LCL = unname(s$conf.int[,"lower .95"]),
    HR_UCL = unname(s$conf.int[,"upper .95"])
  )
}


# Variables to analyze
vars <- c("age_at_diagnosis","tumor_size","lymph_nodes_examined_positive")

# Run univariate Cox models and combine results

res <- bind_rows(lapply(vars, run_uni))

# Save results and print summary

write_csv(res, "outputs/cox_univariate.csv")
print(res)
cat("Saved: outputs/cox_univariate.csv\n")
