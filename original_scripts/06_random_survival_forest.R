# 06_random_survival_forest.R
# Random Survival Forest (RSF) modeling, evaluation, and variable importance
#
# Inputs:  outputs/metabric_cleaned.csv
# Outputs: outputs/RSF_Cindex.csv,
#          outputs/RSF_AUC.csv,
#          outputs/RSF_variable_importance.csv,
#          outputs/RSF_variable_importance.png,
#          outputs/RSF_ROC_overlay.png,
#          outputs/RSF_survival_rates.csv (5 and 10 years)


# Setup: Working Directory

proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"

if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)

setwd(proj_dir); if (!dir.exists("outputs")) dir.create("outputs")
cat("Working in:", getwd(), "\n")

# Load Libraries

library(survival)          # For survival analysis functions
library(randomForestSRC)   # Random Survival Forest modeling
library(timeROC)           # Time-dependent ROC and AUC metrics
library(ggplot2)           # Plotting utilities
library(readr)             # Reading/writing CSV
library(dplyr)             # Data manipulation

# Load and preprocess data
df0 <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = FALSE)

# Select required variables and keep only complete cases
vars <- c("event","time","age_at_diagnosis","tumor_size",
          "lymph_nodes_examined_positive","er_status","pr_status","her2_status")
df <- df0[complete.cases(df0[, vars]), vars]

# Convert receptor markers to factor with consistent levels
df$er_status   <- factor(df$er_status,   levels = c("negative","positive"))
df$pr_status   <- factor(df$pr_status,   levels = c("negative","positive"))
df$her2_status <- factor(df$her2_status, levels = c("negative","positive"))


# Fit Random Survival Forest model

set.seed(42)  # For reproducibility

rsf_model <- rfsrc(Surv(time, event) ~ ., data=df, ntree=1000, importance=TRUE)

# Calculate and save C-index (Harrell’s concordance index)

# Use out-of-bag predicted mortality for evaluation
pred_oob <- rsf_model$predicted.oob
cobj <- survival::concordance(Surv(time, event) ~ pred_oob, data=df, reverse=TRUE)

cindex_val <- cobj$concordance
cindex_se <- sqrt(cobj$var)

write_csv(tibble(Model="RSF", C_index=cindex_val, SE=cindex_se),
          "outputs/RSF_Cindex.csv")
cat(sprintf("RSF C-index = %.3f (SE=%.3f)\n", cindex_val, cindex_se))

# Calculate and save time-dependent AUC at 12, 36, 60 months

eval_times <- c(12,36,60)

roc_rsf <- timeROC(T=df$time, delta=df$event, marker=pred_oob,
                   cause=1, times=eval_times, iid=FALSE)

auc_tbl <- tibble(Time_in_Months = roc_rsf$times, AUC = as.numeric(roc_rsf$AUC))

write_csv(auc_tbl, "outputs/RSF_AUC.csv")
print(auc_tbl)

# Plot ROC curves overlay for evaluation times

png("outputs/RSF_ROC_overlay.png", width=1800, height=1200, res=220)

plot(roc_rsf, time=eval_times[1], col=1, title=FALSE)

for (i in 2:length(eval_times)) plot(roc_rsf, time=eval_times[i], col=i, add=TRUE)
legend("bottomright", legend=paste(eval_times, "mo"), col=seq_along(eval_times), lwd=2)

title(main="Time-dependent ROC (RSF)")
dev.off()

# Variable importance calculation and plotting

vip <- sort(rsf_model$importance, decreasing=TRUE)
vip_df <- tibble(variable=names(vip), importance=as.numeric(vip))
write_csv(vip_df, "outputs/RSF_variable_importance.csv")
p_vip <- ggplot(vip_df, aes(x=reorder(variable, importance), y=importance)) +
  geom_col() + coord_flip() +
  labs(title="RSF Variable Importance", x=NULL, y="Permutation importance") +
  theme_minimal(base_size=12)
ggsave("outputs/RSF_variable_importance.png", p_vip, width=7, height=5, dpi=300)

# Compute and save survival probabilities at 5 and 10 years

rsf_pred <- predict(rsf_model, newdata=df)

tgrid <- rsf_pred$time.interest

idx5  <- which.min(abs(tgrid - 60))
idx10 <- which.min(abs(tgrid - 120))

surv5  <- rsf_pred$survival[, idx5]
surv10 <- rsf_pred$survival[, idx10]

sr <- tibble(Model="RSF", Surv5y=mean(surv5, na.rm=TRUE), Surv10y=mean(surv10, na.rm=TRUE))
write_csv(sr, "outputs/RSF_survival_rates.csv"); print(sr)

cat("RSF results saved in outputs/\n")
