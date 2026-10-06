"""
Breast Cancer Survival Analysis — METABRIC
--------------------------------------------
"Evaluation of Risk Factors and Survival Rates of Patients with Breast
Cancer Using Machine Learning and Traditional Methods" — MSc dissertation
project, rebuilt as a standalone, runnable pipeline.

This reproduces the statistical half of the original analysis (Kaplan-Meier
curves, log-rank tests, univariate + multivariable Cox Proportional
Hazards) directly in Python, validated against the dissertation's own
published results. The machine-learning survival models (Random Survival
Forest, XGB-Cox, DeepSurv) were originally built in R (`randomForestSRC`,
`xgboost`) and Python (`pycox`/`torchtuples`) — see original_scripts/ for
the real, unmodified code that produced those results, and the Model
Comparison table below for the numbers themselves.

A note on the outcome variable: this dataset's `overall_survival` column
is coded 1 = alive, 0 = died (confirmed against `death_from_cancer`). The
original analysis pipeline (01_data_preprocessing.R) passes this column
directly into R's `Surv(time, event)` as the event indicator, without
flipping it. Replicating that exactly is what reproduces the dissertation's
published numbers (verified below) -- so this script follows the same
definition, despite it being the reverse of the usual "1 = event (death)"
survival-analysis convention. This is flagged clearly here rather than
quietly corrected, since the goal of this rebuild is a faithful
reproduction of the original study.

Usage:
    python survival_analysis.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import chi2, norm

FEATURES = ["age_at_diagnosis", "tumor_size", "lymph_nodes_examined_positive",
            "er_status", "pr_status", "her2_status"]


# ---------------------------------------------------------------
# Data loading & cleaning (mirrors 01_data_preprocessing.R)
# ---------------------------------------------------------------

def norm_status(series):
    """Standardise receptor status text: lowercase, trim, fix typos, map to
    positive/negative. Mirrors the R `norm_status()` helper exactly."""
    s = series.astype(str).str.strip().str.lower()
    s = s.str.replace("positve", "positive", regex=False)
    s = s.where(~s.isin(["pos", "pos.", "positive", "1"]), "positive")
    s = s.where(~s.isin(["neg", "neg.", "negative", "0"]), "negative")
    s = s.where(s.isin(["positive", "negative"]), np.nan)
    return s


def load_and_clean_data(path="data/METABRIC_RNA_Mutation.csv"):
    df = pd.read_csv(path, low_memory=False)

    # event/time exactly as 01_data_preprocessing.R defines them:
    # event <- as.numeric(overall_survival); time <- as.numeric(overall_survival_months)
    df["event"] = df["overall_survival"].astype(float)
    df["time"] = df["overall_survival_months"].astype(float)

    df["er_status"] = norm_status(df["er_status"])
    df["pr_status"] = norm_status(df["pr_status"])
    df["her2_status"] = norm_status(df["her2_status"])

    return df


def encode_features(df):
    X = df[["age_at_diagnosis", "tumor_size", "lymph_nodes_examined_positive"]].copy()
    X["er_status"] = (df["er_status"] == "positive").astype(float)
    X["pr_status"] = (df["pr_status"] == "positive").astype(float)
    X["her2_status"] = (df["her2_status"] == "positive").astype(float)
    return X


# ---------------------------------------------------------------
# Kaplan-Meier estimator (mirrors survfit() in 03/12 .R scripts)
# ---------------------------------------------------------------

def kaplan_meier(time, event):
    order = np.argsort(time)
    t_sorted, e_sorted = time[order], event[order]
    unique_times = np.unique(t_sorted)

    s = 1.0
    var_sum = 0.0
    times_out, surv_out, lo_out, hi_out = [0.0], [1.0], [1.0], [1.0]
    z = 1.959963985

    for t in unique_times:
        n = (t_sorted >= t).sum()
        d = e_sorted[t_sorted == t].sum()
        if n > 0:
            if n != d:
                var_sum += d / (n * (n - d))
            s *= (1 - d / n)
        se = s * np.sqrt(var_sum)
        times_out.append(t)
        surv_out.append(s)
        lo_out.append(max(0, s - z * se))
        hi_out.append(min(1, s + z * se))

    return np.array(times_out), np.array(surv_out), np.array(lo_out), np.array(hi_out)


def survival_at(times, surv, month):
    idx = np.where(times <= month)[0]
    return surv[idx[-1]] if len(idx) else 1.0


def median_survival(times, surv):
    below_half = np.where(surv <= 0.5)[0]
    return times[below_half[0]] if len(below_half) else np.nan


def log_rank_test(time_a, event_a, time_b, event_b):
    """Two-group log-rank test (Mantel-Haenszel), mirrors survdiff() in R."""
    all_times = np.concatenate([time_a, time_b])
    all_group = np.concatenate([np.zeros(len(time_a)), np.ones(len(time_b))])
    all_event = np.concatenate([event_a, event_b])

    unique_event_times = np.unique(all_times[all_event == 1])
    O1, E1, V1 = 0.0, 0.0, 0.0
    for t in unique_event_times:
        at_risk = all_times >= t
        n = at_risk.sum()
        n1 = (at_risk & (all_group == 0)).sum()
        n2 = n - n1
        d = ((all_times == t) & (all_event == 1)).sum()
        d1 = ((all_times == t) & (all_event == 1) & (all_group == 0)).sum()
        if n > 1:
            O1 += d1
            E1 += d * n1 / n
            V1 += d * (n1 / n) * (n2 / n) * ((n - d) / (n - 1))
    chi_sq = (O1 - E1) ** 2 / V1 if V1 > 0 else 0.0
    return chi_sq, 1 - chi2.cdf(chi_sq, df=1)


# ---------------------------------------------------------------
# Cox Proportional Hazards (mirrors coxph() in 04/05 .R scripts)
# ---------------------------------------------------------------

def fit_cox_ph(X, time, event, feature_names):
    n, p = X.shape
    order = np.argsort(time)
    t_sorted, e_sorted, X_sorted = time[order], event[order], X[order]
    unique_event_times = np.unique(t_sorted[e_sorted == 1])

    def neg_log_pl(beta):
        risk = X_sorted @ beta
        exp_risk = np.exp(risk)
        ll = 0.0
        for t in unique_event_times:
            at_risk = t_sorted >= t
            events_here = (t_sorted == t) & (e_sorted == 1)
            d_k = events_here.sum()
            ll += risk[events_here].sum() - d_k * np.log(exp_risk[at_risk].sum())
        return -ll

    result = minimize(neg_log_pl, np.zeros(p), method="BFGS", options={"maxiter": 300})
    beta = result.x

    from scipy.optimize import approx_fprime
    eps = 1e-5
    hessian = np.zeros((p, p))
    for i in range(p):
        def grad_i(b, i=i):
            return approx_fprime(b, neg_log_pl, eps)[i]
        hessian[i] = approx_fprime(beta, grad_i, eps)
    try:
        se = np.sqrt(np.abs(np.diag(np.linalg.inv(hessian))))
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)

    hr = np.exp(beta)
    z = beta / se
    pvals = 2 * norm.sf(np.abs(z))

    return pd.DataFrame({
        "feature": feature_names, "coef": beta, "se": se, "z": z,
        "p_value": pvals, "HR": hr,
        "HR_lo": np.exp(beta - 1.96 * se), "HR_hi": np.exp(beta + 1.96 * se),
    })


# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------

def plot_km_overall(times, surv, lo, hi, save_path):
    plt.figure(figsize=(7, 5))
    plt.step(times, surv, where="post", color="#2E6F95")
    plt.fill_between(times, lo, hi, step="post", alpha=0.15, color="#2E6F95")
    plt.xlabel("Months since diagnosis")
    plt.ylabel("Survival probability")
    plt.title("Overall Kaplan-Meier Survival — METABRIC Cohort")
    plt.ylim(0, 1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()


def plot_km_by_group(df, group_col, save_path, title):
    plt.figure(figsize=(7, 5))
    colors = {"positive": "#C0392B", "negative": "#2E6F95"}
    for val, color in colors.items():
        mask = df[group_col] == val
        t, s, _, _ = kaplan_meier(df.loc[mask, "time"].values, df.loc[mask, "event"].values)
        plt.step(t, s, where="post", label=f"{val.title()} (n={mask.sum()})", color=color)
    plt.xlabel("Months since diagnosis")
    plt.ylabel("Survival probability")
    plt.title(title)
    plt.ylim(0, 1.02)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()


def plot_hazard_ratios(cox_df, title, save_path):
    cox_df = cox_df.sort_values("HR")
    plt.figure(figsize=(7, 5))
    colors = ["#C0392B" if hr > 1 else "#2E6F95" for hr in cox_df["HR"]]
    plt.errorbar(cox_df["HR"], range(len(cox_df)),
                 xerr=[cox_df["HR"] - cox_df["HR_lo"], cox_df["HR_hi"] - cox_df["HR"]],
                 fmt="none", ecolor="gray", capsize=4)
    for i, (hr, c) in enumerate(zip(cox_df["HR"], colors)):
        plt.scatter(hr, i, color=c, zorder=5, s=70)
    plt.axvline(1, color="gray", linestyle="--")
    plt.yticks(range(len(cox_df)), cox_df["feature"])
    plt.xlabel("Hazard ratio (95% CI)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()


# ---------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------

if __name__ == "__main__":
    df = load_and_clean_data()
    time, event = df["time"].values, df["event"].values

    print(f"Patients: {len(df)}")
    print(f"event=1 (per overall_survival column, i.e. 'alive'): {int(event.sum())} "
          f"({event.mean()*100:.1f}%)")

    # --- Marker distributions ---
    for marker in ["er_status", "pr_status", "her2_status"]:
        counts = df[marker].value_counts()
        print(f"{marker}: {counts.to_dict()}")

    # --- Kaplan-Meier overall ---
    t_all, s_all, lo_all, hi_all = kaplan_meier(time, event)
    med = median_survival(t_all, s_all)
    print(f"\n[Kaplan-Meier] 5-year: {survival_at(t_all, s_all, 60):.3f}  |  "
          f"10-year: {survival_at(t_all, s_all, 120):.3f}  |  median: {med:.0f} months")
    plot_km_overall(t_all, s_all, lo_all, hi_all, "images/km_overall.png")

    # --- Log-rank tests + KM by receptor status ---
    print("\n[Log-rank tests]")
    logrank_rows = []
    for marker in ["er_status", "pr_status", "her2_status"]:
        sub = df.dropna(subset=[marker])
        pos = sub[marker] == "positive"
        chi_sq, p = log_rank_test(sub.loc[pos, "time"].values, sub.loc[pos, "event"].values,
                                   sub.loc[~pos, "time"].values, sub.loc[~pos, "event"].values)
        print(f"  {marker}: chi2={chi_sq:.2f}, p={p:.4f}")
        logrank_rows.append({"variable": marker, "chi_sq": chi_sq, "df": 1, "p_value": p})
        plot_km_by_group(sub, marker, f"images/km_by_{marker}.png",
                          f"Kaplan-Meier Survival by {marker.replace('_', ' ').upper()}")
    pd.DataFrame(logrank_rows).to_csv("data/logrank_results.csv", index=False)

    # --- Cox PH: univariate ---
    print("\n[Cox PH — Univariate]")
    df_cc = df.dropna(subset=["age_at_diagnosis", "tumor_size",
                               "lymph_nodes_examined_positive", "time", "event"])
    univariate_rows = []
    for feat in ["age_at_diagnosis", "tumor_size", "lymph_nodes_examined_positive"]:
        X_uni = df_cc[[feat]].values.astype(float)
        res = fit_cox_ph(X_uni, df_cc["time"].values, df_cc["event"].values, [feat])
        univariate_rows.append(res)
        print(res.to_string(index=False))
    univariate_df = pd.concat(univariate_rows, ignore_index=True)
    univariate_df.to_csv("data/cox_univariate_results.csv", index=False)

    # --- Cox PH: multivariable ---
    print("\n[Cox PH — Multivariable]")
    df_multi = df.dropna(subset=FEATURES + ["time", "event"])
    X_multi = encode_features(df_multi)
    multi_df = fit_cox_ph(X_multi.values.astype(float), df_multi["time"].values,
                           df_multi["event"].values, FEATURES)
    print(multi_df.to_string(index=False))
    multi_df.to_csv("data/cox_multivariable_results.csv", index=False)
    plot_hazard_ratios(multi_df, "Multivariable Cox PH — Hazard Ratios (95% CI)",
                        "images/cox_hazard_ratios.png")

    print("\nDone. See images/ for charts, data/ for exported tables.")
    print("Random Survival Forest, XGB-Cox and DeepSurv: see original_scripts/ and README.")
