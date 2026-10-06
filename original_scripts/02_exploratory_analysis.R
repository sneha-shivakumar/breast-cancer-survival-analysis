# 02_exploratory_analysis.R
# Major Project - Survival Analysis in Breast Cancer
#
# This script performs exploratory data analysis (EDA) on the
# METABRIC dataset, including summary statistics, correlation
# analysis, and distribution plots.
#
# Inputs:  outputs/metabric_cleaned.csv
# Outputs: outputs/summary_statistics.csv
#          outputs/EDA_age_distribution.png
#          outputs/EDA_tumor_size_distribution.png
#          outputs/EDA_lymph_nodes_distribution.png
#          outputs/correlation_matrix.png



# Setup: Working Directory
# Define project directory path and set it as the working directory
proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"

# Create project directory if missing (creates parent dirs recursively if needed)
if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)

# Set working directory to project folder
setwd(proj_dir)

# Ensure 'outputs' folder exists for saving results
if (!dir.exists("outputs")) dir.create("outputs")

cat("Working in:", getwd(), "\n")


# Step 1: Load libraries

# Load required packages for plotting, correlation visualisation and data handling
library(ggplot2)    # plotting
library(corrplot)   # correlation matrix visualisation
library(readr)      # reading/writing CSV files
library(dplyr)      # data manipulation


# Step 2: Load cleaned dataset

# Load the cleaned METABRIC dataset from previous preprocessing step
df <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = FALSE)

# Select key numeric variables to focus on for exploratory analysis
num_vars <- c("age_at_diagnosis", "tumor_size", "lymph_nodes_examined_positive")

# Subset dataframe to only the numeric variables of interest
df_num <- df[, num_vars]


# Step 3: Summary statistics

# Generate summary statistics of numeric variables for overview
summary_stats <- summary(df_num)

# Print summary statistics to console for inspection
print(summary_stats)

# Save summary statistics text output to CSV file in 'outputs' folder
capture.output(summary_stats, file = "outputs/summary_statistics.csv")


# Step 4: Distribution plots

# Plot histograms showing distributions of key numeric variables

# Age at diagnosis distribution plot
p_age <- ggplot(df, aes(x = age_at_diagnosis)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Age Distribution", x = "Age at Diagnosis", y = "Count") +
  theme_minimal(base_size = 12)
ggsave("outputs/EDA_age_distribution.png", p_age, width = 6, height = 4)

# Tumor size distribution plot
p_tumor <- ggplot(df, aes(x = tumor_size)) +
  geom_histogram(binwidth = 5, fill = "salmon", color = "black") +
  labs(title = "Tumor Size Distribution", x = "Tumor Size (mm)", y = "Count") +
  theme_minimal(base_size = 12)
ggsave("outputs/EDA_tumor_size_distribution.png", p_tumor, width = 6, height = 4)

# Lymph nodes examined positive distribution plot
p_nodes <- ggplot(df, aes(x = lymph_nodes_examined_positive)) +
  geom_histogram(binwidth = 1, fill = "lightgreen", color = "black") +
  labs(title = "Lymph Nodes Positive Distribution", x = "Count", y = "Frequency") +
  theme_minimal(base_size = 12)
ggsave("outputs/EDA_lymph_nodes_distribution.png", p_nodes, width = 6, height = 4)


# Step 5: Correlation matrix

# Calculate correlation matrix for numeric variables using complete cases only
cor_matrix <- cor(df_num, use = "complete.obs")

# Save correlation plot as PNG to outputs folder
png("outputs/correlation_matrix.png", width = 800, height = 600)

# Generate colorful upper-triangle correlation plot with coefficient labels
corrplot(cor_matrix, method = "color", type = "upper", addCoef.col = "black",
         tl.col = "black", tl.srt = 45, title = "Correlation Matrix", mar=c(0,0,1,0))
dev.off()

cat("Exploratory analysis complete. Results saved in 'outputs/'\n")
