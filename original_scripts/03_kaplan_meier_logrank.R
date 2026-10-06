# 03_kaplan_meier_logrank.R
# Major Project - Survival Analysis in Breast Cancer
#
# This script generates Kaplan–Meier survival curves for
# ER, PR, and HER2 receptor status groups, performs log-rank tests,
# and saves the results.
#
# Inputs:  outputs/metabric_cleaned.csv
# Outputs: outputs/KM_ER.png
#          outputs/KM_PR.png
#          outputs/KM_HER2.png
#          outputs/logrank_results.txt
#          outputs/Marker_Distributions.csv



# Setup: Working Directory

proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"
if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)
setwd(proj_dir)

if (!dir.exists("outputs")) dir.create("outputs")

cat("Working in:", getwd(), "\n")


# Step 1: Load libraries

library(survival)     # Survival analysis functions
library(survminer)    # Survival plotting utilities
library(readr)        # CSV file reading/writing

# Step 2: Load cleaned dataset

df <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = TRUE)

# Ensure survival data columns are numeric for compatibility with survival analysis
df$event <- as.numeric(df$event)
df$time  <- as.numeric(df$time)


# Step 3: Kaplan–Meier curves


# Generate KM curves by ER status
fit_er <- survfit(Surv(time, event) ~ er_status, data = df)
p_er <- ggsurvplot(fit_er, data = df, pval = TRUE, conf.int = TRUE,
                   risk.table = TRUE, legend.labs = c("ER Negative", "ER Positive"),
                   title = "Kaplan–Meier Survival by ER Status")
ggsave("outputs/KM_ER.png", p_er$plot, width = 6, height = 5)

# Generate KM curves by PR status
fit_pr <- survfit(Surv(time, event) ~ pr_status, data = df)
p_pr <- ggsurvplot(fit_pr, data = df, pval = TRUE, conf.int = TRUE,
                   risk.table = TRUE, legend.labs = c("PR Negative", "PR Positive"),
                   title = "Kaplan–Meier Survival by PR Status")
ggsave("outputs/KM_PR.png", p_pr$plot, width = 6, height = 5)

# Generate KM curves by HER2 status
fit_her2 <- survfit(Surv(time, event) ~ her2_status, data = df)
p_her2 <- ggsurvplot(fit_her2, data = df, pval = TRUE, conf.int = TRUE,
                     risk.table = TRUE, legend.labs = c("HER2 Negative", "HER2 Positive"),
                     title = "Kaplan–Meier Survival by HER2 Status")
ggsave("outputs/KM_HER2.png", p_her2$plot, width = 6, height = 5)


# Step 4: Log-Rank tests

# Save log-rank test results to a text file
sink("outputs/logrank_results.txt")
cat("Log-Rank Test Results\n")
cat("=====================\n\n")

cat("ER Status:\n")
print(survdiff(Surv(time, event) ~ er_status, data = df))
cat("\n--------------------------------\n")

cat("PR Status:\n")
print(survdiff(Surv(time, event) ~ pr_status, data = df))
cat("\n--------------------------------\n")

cat("HER2 Status:\n")
print(survdiff(Surv(time, event) ~ her2_status, data = df))
cat("\n--------------------------------\n")

sink()

cat("Kaplan–Meier plots and Log-Rank results saved in 'outputs/'\n")


# Step 5: Marker distributions

proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"
setwd(proj_dir)

# Load cleaned dataset
df <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = TRUE)

# Function to compute counts and percentages for negative and positive markers
get_distribution <- function(var) {
  tab <- table(df[[var]])
  total <- sum(tab)
  tibble(
    Marker = var,
    Negative = paste0(tab[1], " (", round(100*tab[1]/total,1), "%)"),
    Positive = paste0(tab[2], " (", round(100*tab[2]/total,1), "%)"),
    Total = total
  )
}

# Compute distributions for ER, PR, and HER2
dist_tbl <- bind_rows(
  get_distribution("er_status"),
  get_distribution("pr_status"),
  get_distribution("her2_status")
)

print(dist_tbl)

# Print and save marker distributions table
write_csv(dist_tbl, "outputs/Marker_Distributions.csv")
cat("Marker distributions saved to outputs/Marker_Distributions.csv\n")
