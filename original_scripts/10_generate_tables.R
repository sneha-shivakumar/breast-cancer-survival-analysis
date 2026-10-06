# 10_generate_tables_and_figures.R
# Produce polished CSV tables mirroring Ozgur-style tables
#
# Inputs : outputs/* from previous scripts
# Outputs: outputs/logrank_summary.csv          # Corresponds to paper's Table 1
#          outputs/univariate_cox_results.csv    # Corresponds to paper's Table 2
#          outputs/multivariable_cox_results.csv # Corresponds to paper's Table 3
#          outputs/model_comparison_summary.csv  # Corresponds to paper's Table 

# Setup: Working Directory

proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"
if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)
setwd(proj_dir); if (!dir.exists("outputs")) dir.create("outputs")
cat("Working in:", getwd(), "\n")

# Load libraries

library(readr); library(dplyr); library(stringr); library(survival)

# Log-rank summary table (paper Table 1)
df <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = TRUE)
LR <- tibble(
  Variable = c("ER status","PR status","HER2 status"),
  Chisq = c(
    survdiff(Surv(time, event) ~ er_status, data=df)$chisq,
    survdiff(Surv(time, event) ~ pr_status, data=df)$chisq,
    survdiff(Surv(time, event) ~ her2_status, data=df)$chisq
  ),
  df = 1L,
  p = pchisq(Chisq, df=1, lower.tail = FALSE)
)
write_csv(LR, "outputs/logrank_summary.csv")

# Univariate Cox results (paper Table 2)
file.copy("outputs/cox_univariate.csv", "outputs/univariate_cox_results.csv", overwrite = TRUE)

# Multivariable Cox results (paper Table 3)
file.copy("outputs/cox_multivariable.csv", "outputs/multivariable_cox_results.csv", overwrite = TRUE)

# Model comparison summary (paper Table 5)
file.copy("outputs/Model_Comparison_Summary.csv", "outputs/model_comparison_summary.csv", overwrite = TRUE)

cat("Final tables saved in outputs/: logrank_summary.csv, univariate_cox_results.csv, multivariable_cox_results.csv, model_comparison_summary.csv\n")
