"""
=============================================================================

STAGGERED DIFFERENCE-IN-DIFFERENCES
Callaway & Sant'Anna (2021) Estimator — Python Implementation
-------------------------------------------------------------
Dataset : minwage_panel.csv (simulated minimum wage panel)
Outcome : lemp     — log teen employment
Package : pydid    — Python port of the R `did` package

REQUIRED PACKAGES:
    pip install pydid pandas numpy matplotlib scipy statsmodels linearmodels

=============================================================================
"""

# ── 0. SETUP ──────────────────────────────────────────────────────────────────

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

warnings.filterwarnings("ignore")
np.random.seed(42)

# Create output directories
os.makedirs("figures", exist_ok=True)
os.makedirs("outputs", exist_ok=True)

# Plot style
plt.rcParams.update({
    "figure.figsize": (10, 6),
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
})


# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

df = pd.read_csv("data/minwage_panel.csv")

print("Dataset shape:", df.shape)
print(df.head())
print("\nCohort distribution (year 2000):")
print(df[df["year"] == 2000]["gvar"].value_counts().sort_index())
print("\nSummary statistics:")
print(df[["lemp", "lpop", "lincome", "mw_treat"]].describe().round(3))


# ── 2. DATA QUALITY CHECKS ───────────────────────────────────────────────────

# Check duplicates
dupes = df.duplicated(subset=["state_id", "year"]).sum()
print(f"\nDuplicates: {dupes} (should be 0)")

# Check treatment irreversibility
df_sorted = df.sort_values(["state_id", "year"])
df_sorted["lag_treat"] = df_sorted.groupby("state_id")["mw_treat"].shift(1)
reversals = df_sorted[(df_sorted["lag_treat"].notna()) &
                      (df_sorted["mw_treat"] < df_sorted["lag_treat"])]
print(f"Treatment reversals: {len(reversals)} (should be 0)")

# Pre-treatment balance check (year 2000)
df_2000 = df[df["year"] == 2000]
for var in ["lpop", "lincome", "lemp"]:
    g0 = df_2000.loc[df_2000["ever_treated"] == 0, var]
    g1 = df_2000.loc[df_2000["ever_treated"] == 1, var]
    t, p = stats.ttest_ind(g0, g1)
    print(f"Balance {var}: never-treated mean={g0.mean():.3f}, "
          f"treated mean={g1.mean():.3f}, p={p:.3f}")


# ── 3. TWFE BENCHMARK ────────────────────────────────────────────────────────

try:
    from linearmodels.panel import PanelOLS
    from linearmodels import AbsorbingLS

    df_panel = df.set_index(["state_id", "year"])

    twfe = PanelOLS(
        dependent  = df_panel["lemp"],
        exog       = df_panel[["mw_treat"]],
        entity_effects = True,
        time_effects   = True
    ).fit(cov_type="clustered", cluster_entity=True)

    twfe_coef = twfe.params["mw_treat"]
    twfe_se   = twfe.std_errors["mw_treat"]
    twfe_pval = twfe.pvalues["mw_treat"]

    print(f"\nTWFE ATT: {twfe_coef:.4f} (SE={twfe_se:.4f}, p={twfe_pval:.3f})")
    print("NOTE: May be biased under heterogeneous treatment effects.")

except Exception as e:
    print(f"TWFE (linearmodels): {e}")
    # Fallback: manual demeaning
    from statsmodels.formula.api import ols

    df_dm = df.copy()
    for col in ["lemp", "mw_treat"]:
        df_dm[f"{col}_dm"] = (df_dm[col]
                              - df_dm.groupby("state_id")[col].transform("mean")
                              - df_dm.groupby("year")[col].transform("mean")
                              + df_dm[col].mean())

    twfe_sm = ols("lemp_dm ~ mw_treat_dm", data=df_dm).fit()
    twfe_coef = twfe_sm.params["mw_treat_dm"]
    twfe_se   = twfe_sm.bse["mw_treat_dm"]
    print(f"TWFE (manual demean) ATT: {twfe_coef:.4f} (SE={twfe_se:.4f})")


# ── 4. TWFE EVENT STUDY ───────────────────────────────────────────────────────

df["event_time"] = df["year"] - df["treat_year"]

# Restrict to window [-4, +5], excluding never-treated from event-time regression
df_es = df[(df["event_time"] >= -4) & (df["event_time"] <= 5)].copy()

# Create dummies for event time (omitting -1 as reference)
event_times = sorted([t for t in df_es["event_time"].unique() if t != -1])
for t in event_times:
    df_es[f"et_{t}".replace("-", "m")] = (df_es["event_time"] == t).astype(int)

et_cols = [f"et_{t}".replace("-", "m") for t in event_times]

try:
    df_es_panel = df_es.set_index(["state_id", "year"])
    twfe_es_model = PanelOLS(
        dependent  = df_es_panel["lemp"],
        exog       = df_es_panel[et_cols],
        entity_effects = True,
        time_effects   = True
    ).fit(cov_type="clustered", cluster_entity=True)
    twfe_es_coefs = twfe_es_model.params[et_cols]
    twfe_es_ses   = twfe_es_model.std_errors[et_cols]
except Exception:
    twfe_es_coefs = pd.Series(np.nan, index=et_cols)
    twfe_es_ses   = pd.Series(np.nan, index=et_cols)

twfe_es_df = pd.DataFrame({
    "event_time": event_times,
    "coef": twfe_es_coefs.values,
    "se":   twfe_es_ses.values
}).assign(ci_lo=lambda x: x["coef"] - 1.96 * x["se"],
          ci_hi=lambda x: x["coef"] + 1.96 * x["se"])
# Insert reference period
ref_row = pd.DataFrame({"event_time": [-1], "coef": [0.0], "se": [0.0],
                         "ci_lo": [0.0], "ci_hi": [0.0]})
twfe_es_df = pd.concat([twfe_es_df, ref_row]).sort_values("event_time")


# ── 5. CALLAWAY & SANT'ANNA (pydid) ──────────────────────────────────────────

try:
    from pydid import att_gt, aggte

    # 5a. Unconditional PT (no covariates)
    cs_nocov = att_gt(
        yname         = "lemp",
        tname         = "year",
        idname        = "state_id",
        gname         = "gvar",
        control_group = "notyettreated",
        est_method    = "dr",
        data          = df,
        print_details = False
    )

    # 5b. Conditional PT (with covariates)
    cs_cov = att_gt(
        yname         = "lemp",
        tname         = "year",
        idname        = "state_id",
        gname         = "gvar",
        xformla       = ["lpop", "lincome"],
        control_group = "notyettreated",
        est_method    = "dr",
        data          = df,
        print_details = False
    )

    # ── 6. Aggregations ───────────────────────────────────────────────────────

    att_simple  = aggte(cs_cov, type="simple")
    att_group   = aggte(cs_cov, type="group")
    att_cal     = aggte(cs_cov, type="calendar")
    att_event   = aggte(cs_cov, type="dynamic")

    print(f"\n{'='*60}")
    print("  CALLAWAY & SANT'ANNA RESULTS")
    print(f"{'='*60}")
    print(f"  Overall ATT (simple): {att_simple.overall_att:.4f} "
          f"(SE={att_simple.overall_se:.4f})")

    # Extract event-study estimates
    es_data = pd.DataFrame({
        "event_time": att_event.egt,
        "att":        att_event.att_egt,
        "se":         att_event.se_egt,
    }).assign(ci_lo=lambda x: x["att"] - 1.96 * x["se"],
              ci_hi=lambda x: x["att"] + 1.96 * x["se"],
              period=lambda x: np.where(x["event_time"] < 0,
                                         "Pre-treatment", "Post-treatment"))

    CS_AVAILABLE = True

except ImportError:
    print("\npydid not available — using manual CS approximation instead.")
    print("Install: pip install pydid\n")
    CS_AVAILABLE = False

    # ── Manual CS approximation ───────────────────────────────────────────────
    # For each (cohort g, time t), compare cohort g to not-yet-treated units
    # using a simple DiD (no doubly-robust weighting).
    # This is illustrative — use pydid for the true estimator.

    att_gt_rows = []
    cohorts = sorted([g for g in df["gvar"].unique() if g > 0])
    all_years = sorted(df["year"].unique())

    for g in cohorts:
        for t in all_years:
            # Treated units: cohort g
            treated = df[(df["gvar"] == g) & (df["year"] == t)]["lemp"]
            treated_pre = df[(df["gvar"] == g) & (df["year"] == g - 1)]["lemp"]
            # Control units: not yet treated at time t
            ctrl = df[(df["gvar"] > t) | (df["gvar"] == 0)]
            ctrl_t   = ctrl[ctrl["year"] == t]["lemp"]
            ctrl_pre = ctrl[ctrl["year"] == g - 1]["lemp"]

            if len(treated) < 2 or len(ctrl_t) < 2:
                continue

            att = ((treated.mean() - treated_pre.mean())
                   - (ctrl_t.mean() - ctrl_pre.mean()))
            n = len(treated)
            se = np.std(treated.values - treated_pre.mean() -
                        ctrl_t.mean() + ctrl_pre.mean(), ddof=1) / np.sqrt(n)

            att_gt_rows.append({
                "g": g, "t": t,
                "att": att, "se": se,
                "event_time": t - g
            })

    att_gt_df = pd.DataFrame(att_gt_rows)

    # Dynamic aggregation (event study)
    es_data = (att_gt_df
               .groupby("event_time")[["att", "se"]]
               .mean()
               .reset_index()
               .assign(ci_lo=lambda x: x["att"] - 1.96 * x["se"],
                       ci_hi=lambda x: x["att"] + 1.96 * x["se"],
                       period=lambda x: np.where(x["event_time"] < 0,
                                                   "Pre-treatment",
                                                   "Post-treatment")))

    # Overall ATT (post-treatment periods only)
    overall_att = att_gt_df[att_gt_df["event_time"] >= 0]["att"].mean()
    overall_se  = att_gt_df[att_gt_df["event_time"] >= 0]["se"].mean()
    print(f"\nManual CS approx. Overall ATT: {overall_att:.4f} (SE≈{overall_se:.4f})")
    print("(Use pydid for the proper doubly-robust estimator)")


# ── 7. PRE-TRENDS TEST ────────────────────────────────────────────────────────

pre = es_data[es_data["event_time"] < 0].copy()
print("\nPre-treatment ATTs:")
print(pre[["event_time", "att", "se"]].round(4).to_string(index=False))

# Joint significance test (Wald)
if len(pre) > 0:
    chi2 = np.sum((pre["att"] / pre["se"]) ** 2)
    pval = 1 - stats.chi2.cdf(chi2, df=len(pre))
    print(f"\nJoint pre-trends test: chi2({len(pre)}) = {chi2:.3f}, p = {pval:.3f}")
    if pval > 0.05:
        print("→ No significant evidence against parallel trends (p > 0.05).")
    else:
        print("→ Pre-trends may be present — check assumptions carefully.")


# ── 8. EVENT STUDY PLOTS ─────────────────────────────────────────────────────

# ── 8a. CS Event Study ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 6))

colors = es_data["period"].map({"Pre-treatment": "#2196F3",
                                 "Post-treatment": "#FF5722"})
ax.scatter(es_data["event_time"], es_data["att"],
           color=colors, zorder=5, s=60)
ax.errorbar(es_data["event_time"], es_data["att"],
            yerr=1.96 * es_data["se"],
            fmt="none", color="grey", alpha=0.6, capsize=4, lw=1.2)
ax.plot(es_data["event_time"], es_data["att"],
        color="grey", alpha=0.4, linewidth=0.8, zorder=3)

ax.axhline(0, color="black", linestyle="--", linewidth=1)
ax.axvline(-0.5, color="firebrick", linestyle="--", linewidth=1.2, alpha=0.8)

patch_pre  = mpatches.Patch(color="#2196F3", label="Pre-treatment")
patch_post = mpatches.Patch(color="#FF5722", label="Post-treatment")
ax.legend(handles=[patch_pre, patch_post], loc="upper left")

ax.set_xlabel("Event Time (t − g)", fontsize=13)
ax.set_ylabel("ATT Estimate", fontsize=13)
ax.set_title("Callaway-Sant'Anna Event Study\nEffect of Minimum Wage on Log Teen Employment",
             fontsize=14, fontweight="bold")
ax.annotate("Treatment onset", xy=(-0.5, ax.get_ylim()[1] * 0.9),
            xytext=(0.5, ax.get_ylim()[1] * 0.9),
            color="firebrick", fontsize=10,
            arrowprops=dict(arrowstyle="->", color="firebrick", lw=1.2))

plt.tight_layout()
plt.savefig("figures/Py_cs_event_study.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: figures/Py_cs_event_study.png")

# ── 8b. TWFE vs. CS comparison ────────────────────────────────────────────────
if twfe_es_df["coef"].notna().all():
    fig, ax = plt.subplots(figsize=(11, 6))

    ax.scatter(twfe_es_df["event_time"] - 0.1, twfe_es_df["coef"],
               color="#E53935", label="TWFE", zorder=5, s=60, marker="s")
    ax.errorbar(twfe_es_df["event_time"] - 0.1, twfe_es_df["coef"],
                yerr=1.96 * twfe_es_df["se"],
                fmt="none", color="#E53935", alpha=0.5, capsize=4)

    ax.scatter(es_data["event_time"] + 0.1, es_data["att"],
               color="#1565C0", label="CS DiD", zorder=5, s=60, marker="o")
    ax.errorbar(es_data["event_time"] + 0.1, es_data["att"],
                yerr=1.96 * es_data["se"],
                fmt="none", color="#1565C0", alpha=0.5, capsize=4)

    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.axvline(-0.5, color="grey", linestyle="--", linewidth=1, alpha=0.7)
    ax.legend(fontsize=12)
    ax.set_xlabel("Event Time (t − g)", fontsize=13)
    ax.set_ylabel("Estimate", fontsize=13)
    ax.set_title("TWFE vs. Callaway-Sant'Anna Event Study Comparison",
                 fontsize=14, fontweight="bold")

    plt.tight_layout()
    plt.savefig("figures/Py_twfe_vs_cs.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Saved: figures/Py_twfe_vs_cs.png")

# ── 8c. Cohort-level ATTs ─────────────────────────────────────────────────────
if CS_AVAILABLE:
    group_df = pd.DataFrame({
        "cohort": att_group.groups,
        "att":    att_group.att_egt,
        "se":     att_group.se_egt,
    }).assign(ci_lo=lambda x: x["att"] - 1.96 * x["se"],
              ci_hi=lambda x: x["att"] + 1.96 * x["se"])
else:
    group_df = (att_gt_df[att_gt_df["event_time"] >= 0]
                .groupby("g")[["att", "se"]].mean()
                .reset_index()
                .rename(columns={"g": "cohort"})
                .assign(ci_lo=lambda x: x["att"] - 1.96 * x["se"],
                        ci_hi=lambda x: x["att"] + 1.96 * x["se"]))

fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(group_df["cohort"].astype(str), group_df["att"],
       color=["#42A5F5", "#66BB6A", "#FFA726", "#EF5350"],
       alpha=0.8, width=0.6)
ax.errorbar(range(len(group_df)), group_df["att"],
            yerr=1.96 * group_df["se"],
            fmt="none", color="black", capsize=5, lw=1.5)
ax.axhline(0, color="black", linestyle="--", linewidth=1)
ax.set_xticks(range(len(group_df)))
ax.set_xticklabels(group_df["cohort"].astype(str))
ax.set_xlabel("Treatment Cohort (First Year Treated)", fontsize=13)
ax.set_ylabel("ATT Estimate", fontsize=13)
ax.set_title("Average Treatment Effect by Cohort", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("figures/Py_cohort_atts.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: figures/Py_cohort_atts.png")


# ── 9. SENSITIVITY ANALYSIS ───────────────────────────────────────────────────

sensitivity_results = []

if CS_AVAILABLE:
    # Never-treated only
    cs_nt = att_gt(yname="lemp", tname="year", idname="state_id", gname="gvar",
                   xformla=["lpop", "lincome"], control_group="nevertreated",
                   est_method="dr", data=df, print_details=False)
    a = aggte(cs_nt, type="simple")
    sensitivity_results.append({"method": "CS (Never-Treated)", "att": a.overall_att, "se": a.overall_se})

    # Regression adjustment
    cs_ra = att_gt(yname="lemp", tname="year", idname="state_id", gname="gvar",
                   xformla=["lpop", "lincome"], control_group="notyettreated",
                   est_method="reg", data=df, print_details=False)
    a = aggte(cs_ra, type="simple")
    sensitivity_results.append({"method": "CS (Reg-Adj)", "att": a.overall_att, "se": a.overall_se})

# Always include TWFE
sensitivity_results.insert(0, {"method": "TWFE", "att": twfe_coef, "se": twfe_se})
if CS_AVAILABLE:
    sensitivity_results.insert(1, {"method": "CS (DR, Not-Yet-Treated)",
                                    "att": att_simple.overall_att,
                                    "se": att_simple.overall_se})


# ── 10. RESULTS SUMMARY ───────────────────────────────────────────────────────

results_df = pd.DataFrame(sensitivity_results).assign(
    tstat = lambda x: x["att"] / x["se"],
    sig   = lambda x: pd.cut(x["att"].abs() / x["se"],
                             bins=[-np.inf, 1.645, 1.960, 2.576, np.inf],
                             labels=["", "*", "**", "***"])
).round(4)

print("\n" + "="*65)
print("  RESULTS SUMMARY")
print("="*65)
print(results_df[["method", "att", "se", "tstat", "sig"]].to_string(index=False))
print("="*65)

results_df.to_csv("outputs/Py_results_summary.csv", index=False)
print("\nResults saved to: outputs/Py_results_summary.csv")
print("Figures saved to: figures/")
