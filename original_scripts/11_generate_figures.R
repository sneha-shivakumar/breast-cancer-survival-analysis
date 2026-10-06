# 11_generate_figures.R
# Generate and save publication-ready figures
#
# Inputs : outputs/* from previous scripts
# Outputs: outputs/figures/*.png


# Setup: Working Directory

proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"
if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)
setwd(proj_dir)

# Create output folders if missing (main and figures)
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("outputs/figures")) dir.create("outputs/figures")

cat("Working in:", getwd(), "\n")

# Load libraries

library(ggplot2)    # Data visualization
library(corrplot)   # Correlation matrix visualization
library(readr)      # CSV reading/writing
library(dplyr)      # Data manipulation
library(survival)   # Survival analysis
library(survminer)  # Survival curve plotting


# 1. Distribution plots (EDA)

df <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = TRUE)

# Age distribution
p_age <- ggplot(df, aes(x = age_at_diagnosis)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title="Age Distribution", x="Age at Diagnosis", y="Count") +
  theme_minimal(base_size=12)
ggsave("outputs/figures/Fig1_Age_Distribution.png", p_age, width=6, height=4)

# Tumor size
p_tumor <- ggplot(df, aes(x = tumor_size)) +
  geom_histogram(binwidth = 5, fill = "salmon", color = "black") +
  labs(title="Tumor Size Distribution", x="Tumor Size (mm)", y="Count") +
  theme_minimal(base_size=12)
ggsave("outputs/figures/Fig2_Tumor_Size_Distribution.png", p_tumor, width=6, height=4)

# Lymph nodes examined positive
p_nodes <- ggplot(df, aes(x = lymph_nodes_examined_positive)) +
  geom_histogram(binwidth = 1, fill = "lightgreen", color = "black") +
  labs(title="Lymph Nodes Positive Distribution", x="Count", y="Frequency") +
  theme_minimal(base_size=12)
ggsave("outputs/figures/Fig3_Lymph_Nodes_Distribution.png", p_nodes, width=6, height=4)


# 2. Correlation matrix

num_vars <- c("age_at_diagnosis", "tumor_size", "lymph_nodes_examined_positive")
cor_matrix <- cor(df[, num_vars], use="complete.obs")
png("outputs/figures/Fig4_Correlation_Matrix.png", width=800, height=600)
corrplot(cor_matrix, method="color", type="upper", addCoef.col="black",
         tl.col="black", tl.srt=45, title="Correlation Matrix", mar=c(0,0,1,0))
dev.off()


# 3. Kaplan–Meier curves

df$event <- as.numeric(df$event)
df$time  <- as.numeric(df$time)

# By ER status

fit_er <- survfit(Surv(time, event) ~ er_status, data=df)
p_er <- ggsurvplot(fit_er, data=df, pval=TRUE, conf.int=TRUE,
                   risk.table=TRUE, legend.labs=c("ER Negative","ER Positive"),
                   title="Kaplan–Meier Survival by ER Status")
ggsave("outputs/figures/Fig5_KM_ER.png", p_er$plot, width=6, height=5)

# By PR status

fit_pr <- survfit(Surv(time, event) ~ pr_status, data=df)
p_pr <- ggsurvplot(fit_pr, data=df, pval=TRUE, conf.int=TRUE,
                   risk.table=TRUE, legend.labs=c("PR Negative","PR Positive"),
                   title="Kaplan–Meier Survival by PR Status")
ggsave("outputs/figures/Fig6_KM_PR.png", p_pr$plot, width=6, height=5)

## By HER2 status

fit_her2 <- survfit(Surv(time, event) ~ her2_status, data=df)
p_her2 <- ggsurvplot(fit_her2, data=df, pval=TRUE, conf.int=TRUE,
                     risk.table=TRUE, legend.labs=c("HER2 Negative","HER2 Positive"),
                     title="Kaplan–Meier Survival by HER2 Status")
ggsave("outputs/figures/Fig7_KM_HER2.png", p_her2$plot, width=6, height=5)


# 7) Save overall Kaplan–Meier curve as PNG
p <- ggsurvplot(
  fit_all, data = df, conf.int = TRUE, risk.table = TRUE,
  ggtheme = theme_minimal(base_size = 12),
  title = "Overall Kaplan–Meier Survival (All Patients)",
  xlab = "Time (months)", ylab = "Survival probability"
)
ggsave("outputs/figures/Fig8_Overall_KM_Survival.png", p$plot, width = 7, height = 5, dpi = 300)



# 4. Variable importance plots (RSF, XGB if available)

if (file.exists("outputs/RSF_variable_importance.csv")) {
  vip_rsf <- read_csv("outputs/RSF_variable_importance.csv")
  p_vip_rsf <- ggplot(vip_rsf, aes(x=reorder(variable, importance), y=importance)) +
    geom_col(fill="steelblue") + coord_flip() +
    labs(title="RSF Variable Importance", x=NULL, y="Permutation Importance") +
    theme_minimal(base_size=12)
  ggsave("outputs/figures/Fig9_RSF_VariableImportance.png", p_vip_rsf, width=7, height=5)
}

if (file.exists("outputs/XGB_COX_feature_importance.csv")) {
  vip_xgb <- read_csv("outputs/XGB_COX_feature_importance.csv")
  p_vip_xgb <- ggplot(vip_xgb[1:10,], aes(x=reorder(Feature, Gain), y=Gain)) +
    geom_col(fill="darkorange") + coord_flip() +
    labs(title="XGB-COX Top 10 Feature Importance", x=NULL, y="Gain") +
    theme_minimal(base_size=12)
  ggsave("outputs/figures/Fig10_XGB_FeatureImportance.png", p_vip_xgb, width=7, height=5)
}

cat("All figures saved in outputs/figures/\n")
