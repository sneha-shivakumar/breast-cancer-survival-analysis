# 09_model_comparison.R
# Consolidate model results from previous scripts: C-index,
# time-dependent AUC, and survival rates into comparison table
#
# Inputs:  outputs/* files from prior modeling scripts
# Outputs: outputs/Model_Comparison_Summary.csv

# Setup: Working Directory

proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"

if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)

setwd(proj_dir); if (!dir.exists("outputs")) dir.create("outputs")
cat("Working in:", getwd(), "\n")

# Load libraries

library(readr)   # Reading and writing CSV files
library(dplyr)   # Data manipulation
library(tidyr)   # Tools for reshaping and tidying data (e.g., pivoting,
                 #separating, combining columns) to create tidy datasets


# Load model C-index and SE

c_cph <- tryCatch(read_csv("outputs/CPH_Cindex.csv"), error=function(e) NULL)
# Optional if computed
c_rsf <- read_csv("outputs/RSF_Cindex.csv")
c_xgb <- read_csv("outputs/XGB_COX_Cindex.csv")

# DeepSurv C-index from text file (optional)
ds_cidx <- tryCatch({
  val <- readLines("outputs/DeepSurv_cindex.txt")
  as.numeric(sub(".*: ", "", val[1]))
}, error=function(e) NA_real_)

# Load model AUCs
a_rsf <- read_csv("outputs/RSF_AUC.csv")  %>% mutate(Model="RSF")
a_xgb <- read_csv("outputs/XGB_COX_AUC.csv") %>% mutate(Model="XGB-COX")

# Load survival rates (5y/10y)

s_cph <- tryCatch(read_csv("outputs/CPH_survival_rates.csv"), error=function(e) NULL)
s_rsf <- read_csv("outputs/RSF_survival_rates.csv")

# XGB and DeepSurv can be NA/blank if not computed (or)
# Placeholders for models without survival rates computed
s_xgb <- tibble(Model="XGB-COX", Surv5y=NA_real_, Surv10y=NA_real_)
s_ds  <- tibble(Model="DeepSurv", Surv5y=NA_real_, Surv10y=NA_real_)


# Helper function to format C-index with confidence interval

fmt_ci <- function(cidx, se) {
  if (is.na(cidx) || is.na(se)) return(NA_character_)
  lcl <- cidx - 1.96*se; ucl <- cidx + 1.96*se
  sprintf("%.3f ± %.3f (%.3f–%.3f)", cidx, se, lcl, ucl)
}

# Build comparison table rows

row_rsf <- tibble(
  Model = "RSF",
  `C-index ± SE (95% CI)` = fmt_ci(c_rsf$C_index[1], c_rsf$SE[1]),
  C_index = c_rsf$C_index[1], SE = c_rsf$SE[1],
  AUC12 = a_rsf$AUC[a_rsf$Time_in_Months==12],
  AUC36 = a_rsf$AUC[a_rsf$Time_in_Months==36],
  AUC60 = a_rsf$AUC[a_rsf$Time_in_Months==60]
) %>% left_join(s_rsf, by=character()) %>%
  rename(`5y_Survival`=Surv5y, `10y_Survival`=Surv10y)

row_xgb <- tibble(
  Model = "XGB-COX",
  `C-index ± SE (95% CI)` = fmt_ci(c_xgb$C_index[1], c_xgb$SE[1]),
  C_index = c_xgb$C_index[1], SE = c_xgb$SE[1],
  AUC12 = a_xgb$AUC[a_xgb$Time_in_Months==12],
  AUC36 = a_xgb$AUC[a_xgb$Time_in_Months==36],
  AUC60 = a_xgb$AUC[a_xgb$Time_in_Months==60],
  `5y_Survival` = NA_real_,
  `10y_Survival` = NA_real_
)

row_ds <- tibble(
  Model = "DeepSurv",
  `C-index ± SE (95% CI)` = ifelse(is.na(ds_cidx), NA_character_, sprintf("%.3f (unstable)", ds_cidx)),
  C_index = ds_cidx, SE = NA_real_,
  AUC12 = NA_real_, AUC36 = NA_real_, AUC60 = NA_real_,
  `5y_Survival` = NA_real_, `10y_Survival` = NA_real_
)

# Optionally add CPH row if available

row_cph <- tryCatch({
  ci <- read_csv("outputs/CPH_Cindex.csv")  # if you create it
  sr <- read_csv("outputs/CPH_survival_rates.csv")
  tibble(
    Model = "CPH",
    `C-index ± SE (95% CI)` = fmt_ci(ci$C_index[1], ci$SE[1]),
    C_index = ci$C_index[1], SE = ci$SE[1],
    AUC12 = NA_real_, AUC36 = NA_real_, AUC60 = NA_real_
  ) %>% bind_cols(sr %>% select(Surv5y, Surv10y)) %>%
    rename(`5y_Survival`=Surv5y, `10y_Survival`=Surv10y)
}, error=function(e) NULL)


# Combine rows into final table and save
model_comparison_summary <- bind_rows(list(row_cph, row_rsf, row_xgb, row_ds)) %>% filter(!is.na(Model))

write_csv(model_comparison_summary, "outputs/Model_Comparison_Summary.csv")
print(model_comparison_summary)
cat("Saved: outputs/Model_Comparison_Summary.csv\n")
