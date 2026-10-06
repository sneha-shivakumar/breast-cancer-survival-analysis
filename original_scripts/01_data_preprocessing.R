# 01_data_preprocessing.R
# Major Project - Survival Analysis in Breast Cancer
#
# This script loads the raw METABRIC dataset, cleans variables,
# creates survival columns (event, time), and saves a cleaned dataset
# for downstream survival analysis.
# It also ensures an "outputs" folder exists to store results.


# Setup: Working Directory

# Set working directory to the project folder.
# Modify the path below to your actual project location.
setwd("~/Documents/Major_Project_Code")

# Confirm the working directory is set correctly
cat("Current working directory:", getwd(), "\n")


# Step 1: Create outputs folder

# Create an 'outputs' folder to save result files if it doesn't already exist.
if (!dir.exists("outputs")) {
  dir.create("outputs")
  cat("Created 'outputs' folder\n")
} else {
  cat("'outputs' folder already exists\n")
}


# Step 2: Load libraries

# Load necessary libraries for data manipulation and CSV handling
library(dplyr)
library(readr)


# Step 3: Load dataset

# Read raw METABRIC dataset CSV file into a dataframe
# stringsAsFactors = FALSE ensures character columns are not converted to factors
metabric <- read.csv("METABRIC_RNA_Mutation.csv", stringsAsFactors = FALSE)


# Step 4: Create survival outcome columns

# 'event' column: 1 if death occurred, else 0 (censored)
# 'time' column: survival time in months
metabric$event <- as.numeric(metabric$overall_survival)
metabric$time  <- as.numeric(metabric$overall_survival_months)


# Step 5: Clean receptor status variables

# Function to standardise receptor status values by:
# - Lowercasing text
# - Trimming whitespace
# - Fixing typos (e.g. "positve" -> "positive")
# - Mapping various synonymous inputs to 'positive' or 'negative'
# - Replacing unrecognised values with NA
norm_status <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("positve", "positive", x)   # Fix common typo
  x[x %in% c("pos","pos.","positive","1")] <- "positive"
  x[x %in% c("neg","neg.","negative","0")] <- "negative"
  x[!(x %in% c("positive","negative"))] <- NA
  x
}

# Apply normalisation function to ER, PR, and HER2 receptor status columns
metabric$er_status   <- norm_status(metabric$er_status)
metabric$pr_status   <- norm_status(metabric$pr_status)
metabric$her2_status <- norm_status(metabric$her2_status)


# Step 6: Save cleaned dataset

# Export the cleaned dataframe to CSV for use in subsequent analysis scripts
write_csv(metabric, "outputs/metabric_cleaned.csv")

cat("Data preprocessing complete. Cleaned file saved to outputs/metabric_cleaned.csv\n")
