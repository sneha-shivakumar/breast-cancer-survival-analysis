# 12_survival_overall_and_fixed_times.R (FINAL)
# Overall and time-fixed survival analysis for METABRIC cohort
#
# - Computes KM survival at 5 & 10 years, distribution beyond 10 years,
#   and overall median survival.
# - Produces summary tables and publication-ready figure.
#
# Outputs:
#   outputs/Survival_Summary_Table.csv
#   outputs/Followup_Distribution.csv
#   outputs/Median_Survival.txt
#   outputs/Overall_KM_Survival.png

# Setup: Working Directory

proj_dir <- "/Users/snehakamanahallishivakumar/Documents/Major_Project_Code"
setwd(proj_dir)
if (!dir.exists("outputs")) dir.create("outputs")

# Load Libraries

library(survival)    # Survival analysis functions
library(survminer)   # KM curves and plotting
library(dplyr)       # Data manipulation
library(readr)       # Reading/writing CSV files


# 1) Load cleaned dataset
df <- read.csv("outputs/metabric_cleaned.csv", stringsAsFactors = TRUE)
df$time  <- as.numeric(df$time)
df$event <- as.numeric(df$event)

# 2) Overall Kaplan–Meier fit (all patients)
fit_all <- survfit(Surv(time, event) ~ 1, data = df)

# 3) Survival probability at fixed times (5 years, 10 years)
times <- c(60, 120) # 5 yrs = 60 mo, 10 yrs = 120 mo
s_fix <- summary(fit_all, times = times)

fixed_tbl <- tibble(
  Time_label   = c("5-year (≤60m)", "10-year (≤120m)"),
  Time_months  = s_fix$time,
  Survival     = s_fix$surv,
  Lower_95CI   = s_fix$lower,
  Upper_95CI   = s_fix$upper,
  N_at_risk    = s_fix$n.risk
)

# 4) Median overall survival
med_tab  <- summary(fit_all)$table
med_line <- sprintf(
  "Median survival: %.1f months (95%% CI: %.1f–%.1f)",
  as.numeric(med_tab["median"]),
  as.numeric(med_tab["0.95LCL"]),
  as.numeric(med_tab["0.95UCL"])
)
writeLines(med_line, "outputs/Median_Survival.txt")

# 5) Follow-up time distribution bands (for >10 years, reporting count/percent)
cuts <- cut(
  df$time,
  breaks = c(-Inf, 60, 120, Inf),
  labels = c("0–5 years (≤60m)", "5–10 years (60–120m)", ">10 years (>120m)"),
  right = TRUE
)
followup_tbl <- as.data.frame(table(cuts))
colnames(followup_tbl) <- c("Followup_Band", "Count")
followup_tbl <- followup_tbl %>% mutate(Percent = round(100 * Count / sum(Count), 1))
write_csv(followup_tbl, "outputs/Followup_Distribution.csv")

# Row for >10 years: distribution only, not a KM point
gt10 <- followup_tbl %>% filter(Followup_Band == ">10 years (>120m)")
gt10_row <- tibble(
  Time_label   = ">10 years (>120m)",
  Time_months  = NA_real_,
  Survival     = NA_real_,                 # reporting proportion only
  Lower_95CI   = NA_real_,
  Upper_95CI   = NA_real_,
  N_at_risk    = as.numeric(gt10$Count)
)

# 6) Build a neat summary table: 5y, 10y, >10y, median survival
median_row <- tibble(
  Time_label   = "Overall (Median survival)",
  Time_months  = as.numeric(med_tab["median"]),
  Survival     = 0.50,  # by definition at the median
  Lower_95CI   = as.numeric(med_tab["0.95LCL"]),
  Upper_95CI   = as.numeric(med_tab["0.95UCL"]),
  N_at_risk    = NA_real_
)

# Combine all rows into reporting table
summary_table <- bind_rows(fixed_tbl, gt10_row, median_row) %>%
  mutate(
    `Survival Probability` = ifelse(is.na(Survival), "",
                                    sprintf("%.3f", Survival)),
    `95% CI` = case_when(
      Time_label == "Overall (Median survival)" ~ sprintf("%.1f–%.1f months", Lower_95CI, Upper_95CI),
      is.na(Survival)                           ~ "",
      TRUE                                      ~ sprintf("%.3f–%.3f", Lower_95CI, Upper_95CI)
    )
  ) %>%
  select(
    `Time Point` = Time_label,
    `Time (months)` = Time_months,
    `Survival Probability`,
    `95% CI`,
    `n at risk` = N_at_risk
  )

write_csv(summary_table, "outputs/Survival_Summary_Table.csv")

# 7) Save overall Kaplan–Meier curve as PNG
p <- ggsurvplot(
  fit_all, data = df, conf.int = TRUE, risk.table = TRUE,
  ggtheme = theme_minimal(base_size = 12),
  title = "Overall Kaplan–Meier Survival (All Patients)",
  xlab = "Time (months)", ylab = "Survival probability"
)
ggsave("outputs/Overall_KM_Survival.png", p$plot, width = 7, height = 5, dpi = 300)

cat("Saved:\n",
    " - outputs/Survival_Summary_Table.csv\n",
    " - outputs/Followup_Distribution.csv\n",
    " - outputs/Median_Survival.txt\n",
    " - outputs/Overall_KM_Survival.png\n")
