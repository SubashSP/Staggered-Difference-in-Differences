# =============================================================================
#
#  STAGGERED DIFFERENCE-IN-DIFFERENCES
#  Callaway & Sant'Anna (2021) Estimator — R Implementation
#  ---------------------------------------------------------
#  Dataset : minwage_panel.csv (simulated minimum wage panel)
#  Outcome : lemp     — log teen employment
#  Package : did      — by Callaway & Sant'Anna
#
#  REQUIRED PACKAGES:
#    install.packages(c("did", "tidyverse", "ggplot2", "fixest",
#                       "bacondecomp", "modelsummary", "patchwork"))
#
# =============================================================================

# ── 0. SETUP ──────────────────────────────────────────────────────────────────

# Install packages if needed (run once)
# install.packages(c("did", "tidyverse", "ggplot2", "fixest",
#                    "bacondecomp", "modelsummary", "patchwork"))

library(did)           # Callaway & Sant'Anna estimator
library(tidyverse)     # data manipulation
library(ggplot2)       # plotting
library(fixest)        # fast TWFE with feols()
library(modelsummary)  # regression tables
library(patchwork)     # combine plots (optional)

# Create output directories if they don't exist
dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

set.seed(42)


# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

df <- read_csv("data/minwage_panel.csv", show_col_types = FALSE)

cat("Dataset dimensions:", nrow(df), "rows ×", ncol(df), "cols\n")
glimpse(df)

# Check cohort distribution
df |>
  filter(year == 2000) |>
  count(gvar) |>
  print()


# ── 2. DATA QUALITY CHECKS ───────────────────────────────────────────────────

# Check for duplicates
dupes <- df |> group_by(state_id, year) |> filter(n() > 1)
cat("Duplicates:", nrow(dupes), "\n")

# Verify treatment irreversibility
reversals <- df |>
  arrange(state_id, year) |>
  group_by(state_id) |>
  mutate(lag_treat = lag(mw_treat, default = 0),
         reversal  = mw_treat < lag_treat) |>
  filter(reversal)
cat("Treatment reversals:", nrow(reversals), "(should be 0)\n")

# Pre-treatment balance: 2000 baseline
balance <- df |>
  filter(year == 2000) |>
  group_by(ever_treated) |>
  summarise(across(c(lemp, lpop, lincome), mean), n = n())
print(balance)


# ── 3. NAIVE TWFE BENCHMARK ──────────────────────────────────────────────────

# Standard TWFE using fixest::feols (handles FE and clustering efficiently)
twfe_model <- feols(lemp ~ mw_treat | state_id + year,
                    data    = df,
                    cluster = ~state_id)
summary(twfe_model)

cat("\nTWFE ATT:", coef(twfe_model)["mw_treat"], "\n")
cat("NOTE: This estimate may be biased under heterogeneous treatment effects.\n")


# ── 4. TWFE EVENT STUDY (BENCHMARK) ──────────────────────────────────────────

# Create event-time variable
df <- df |>
  mutate(event_time = if_else(!is.na(treat_year), year - treat_year, NA_real_))

# TWFE event study: restrict to event window [-4, +5], omit l=-1 as reference
twfe_es <- feols(
  lemp ~ i(event_time, ref = -1) | state_id + year,
  data    = df |> filter(event_time >= -4 & event_time <= 5 | is.na(event_time)),
  cluster = ~state_id
)

# Plot TWFE event study
iplot(twfe_es,
      main    = "TWFE Event Study (Benchmark)",
      xlab    = "Event Time (years since MW increase)",
      ylab    = "Coefficient estimate",
      col     = "firebrick")
abline(h = 0, lty = 2, col = "black")
abline(v = -0.5, lty = 2, col = "grey50")


# ── 5. GOODMAN-BACON DECOMPOSITION ───────────────────────────────────────────

# Decompose TWFE into weighted 2×2 DiD comparisons
library(bacondecomp)

bacon_out <- bacon(lemp ~ mw_treat,
                   data            = df |> filter(!is.na(gvar)),
                   id_var          = "state_id",
                   time_var        = "year")

cat("\nGoodman-Bacon Decomposition:\n")
print(bacon_out |> group_by(type) |>
        summarise(avg_estimate = weighted.mean(estimate, wt),
                  total_weight = sum(wt)))

# Plot
ggplot(bacon_out, aes(x = wt, y = estimate, colour = type)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title   = "Goodman-Bacon Decomposition",
       subtitle = "Each point is a 2×2 DiD comparison",
       x        = "Weight",
       y        = "2×2 DiD Estimate",
       colour   = "Comparison Type") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("figures/R_bacon_decomposition.png", width = 9, height = 6, dpi = 150)


# ── 6. CALLAWAY & SANT'ANNA — MAIN ESTIMATION ────────────────────────────────

# ── 6a. Unconditional Parallel Trends (no covariates) ────────────────────────
cs_nocov <- att_gt(
  yname         = "lemp",
  tname         = "year",
  idname        = "state_id",
  gname         = "gvar",           # 0 for never-treated
  control_group = "notyettreated",  # use not-yet-treated as controls
  est_method    = "dr",             # doubly-robust (RECOMMENDED)
  data          = df,
  print_details = FALSE
)

summary(cs_nocov)

# ── 6b. Conditional Parallel Trends (with covariates) ────────────────────────
cs_cov <- att_gt(
  yname         = "lemp",
  tname         = "year",
  idname        = "state_id",
  gname         = "gvar",
  xformla       = ~ lpop + lincome,  # pre-treatment covariates
  control_group = "notyettreated",
  est_method    = "dr",
  data          = df,
  print_details = FALSE
)

summary(cs_cov)


# ── 7. AGGREGATIONS ───────────────────────────────────────────────────────────

# 7a. Simple (overall) ATT
att_simple <- aggte(cs_cov, type = "simple")
summary(att_simple)
cat("\nOverall ATT:", att_simple$overall.att,
    "| SE:", att_simple$overall.se, "\n")

# 7b. Group/Cohort ATTs
att_group <- aggte(cs_cov, type = "group")
summary(att_group)

# 7c. Calendar-Time ATTs
att_cal <- aggte(cs_cov, type = "calendar")
summary(att_cal)

# 7d. Dynamic / Event-Study ATTs (MOST INFORMATIVE)
att_event <- aggte(cs_cov, type = "dynamic")
summary(att_event)


# ── 8. TESTING PRE-TRENDS ────────────────────────────────────────────────────

# Extract pre-treatment ATTs from the event-study aggregation
pre_periods <- att_event$egt[att_event$egt < 0]
pre_atts    <- att_event$att.egt[att_event$egt < 0]
pre_ses     <- att_event$se.egt[att_event$egt < 0]

cat("\nPre-treatment ATTs:\n")
data.frame(event_time = pre_periods,
           att        = round(pre_atts, 4),
           se         = round(pre_ses, 4),
           tstat      = round(pre_atts / pre_ses, 3)) |>
  print()

cat("If |t-stat| < 1.96 for all pre-periods → no evidence against parallel trends\n")


# ── 9. EVENT STUDY PLOTS ─────────────────────────────────────────────────────

# ── 9a. CS event study using did::ggdid ──────────────────────────────────────
p_cs <- ggdid(att_event) +
  labs(title    = "Callaway-Sant'Anna Event Study",
       subtitle = "Effect of Minimum Wage on Log Teen Employment",
       x        = "Years Relative to Treatment Onset",
       y        = "ATT Estimate") +
  theme_bw(base_size = 13)

ggsave("figures/R_cs_event_study.png", p_cs, width = 10, height = 6, dpi = 150)
print(p_cs)

# ── 9b. Group-level ATTs ─────────────────────────────────────────────────────
p_group <- ggdid(att_group) +
  labs(title = "ATT by Treatment Cohort",
       x     = "Cohort (First Treatment Year)",
       y     = "ATT Estimate") +
  theme_bw(base_size = 13)

ggsave("figures/R_cohort_atts.png", p_group, width = 8, height = 5, dpi = 150)
print(p_group)

# ── 9c. Manual ggplot2 event study (publication quality) ─────────────────────
es_df <- data.frame(
  event_time = att_event$egt,
  att        = att_event$att.egt,
  se         = att_event$se.egt
) |>
  mutate(
    ci_lo  = att - 1.96 * se,
    ci_hi  = att + 1.96 * se,
    period = if_else(event_time < 0, "Pre-treatment", "Post-treatment")
  )

p_manual <- ggplot(es_df, aes(x = event_time, y = att, colour = period)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_line(colour = "grey60", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
  scale_colour_manual(values = c("Pre-treatment" = "steelblue",
                                  "Post-treatment" = "darkorange")) +
  labs(title    = "Callaway-Sant'Anna Event Study (Manual ggplot2)",
       subtitle = "95% confidence intervals. Red dashed line = treatment onset.",
       x        = "Event Time (t − g)",
       y        = "ATT(g,t) Estimate",
       colour   = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave("figures/R_event_study_manual.png", p_manual, width = 10, height = 6, dpi = 150)
print(p_manual)


# ── 10. SENSITIVITY ANALYSIS ─────────────────────────────────────────────────

# 10a. Never-treated only as controls
cs_nevertreated <- att_gt(
  yname         = "lemp",
  tname         = "year",
  idname        = "state_id",
  gname         = "gvar",
  xformla       = ~ lpop + lincome,
  control_group = "nevertreated",   # more conservative
  est_method    = "dr",
  data          = df,
  print_details = FALSE
)
att_never_simple <- aggte(cs_nevertreated, type = "simple")
cat("\nCS (never-treated controls) overall ATT:", att_never_simple$overall.att, "\n")

# 10b. Regression adjustment only
cs_ra <- att_gt(
  yname         = "lemp",
  tname         = "year",
  idname        = "state_id",
  gname         = "gvar",
  xformla       = ~ lpop + lincome,
  control_group = "notyettreated",
  est_method    = "reg",            # regression adjustment only
  data          = df,
  print_details = FALSE
)
att_ra_simple <- aggte(cs_ra, type = "simple")
cat("CS (regression adjustment) overall ATT:", att_ra_simple$overall.att, "\n")

# 10c. IPW only
cs_ipw <- att_gt(
  yname         = "lemp",
  tname         = "year",
  idname        = "state_id",
  gname         = "gvar",
  xformla       = ~ lpop + lincome,
  control_group = "notyettreated",
  est_method    = "ipw",            # IPW only
  data          = df,
  print_details = FALSE
)
att_ipw_simple <- aggte(cs_ipw, type = "simple")
cat("CS (IPW) overall ATT:", att_ipw_simple$overall.att, "\n")


# ── 11. SUMMARY TABLE ────────────────────────────────────────────────────────

results_summary <- tibble(
  Method           = c("TWFE", "CS (No Cov, DR)", "CS (Cov, DR)",
                        "CS (Never-Treated)", "CS (Reg-Adj)", "CS (IPW)"),
  ATT              = c(coef(twfe_model)["mw_treat"],
                        aggte(cs_nocov, type = "simple")$overall.att,
                        att_simple$overall.att,
                        att_never_simple$overall.att,
                        att_ra_simple$overall.att,
                        att_ipw_simple$overall.att),
  SE               = c(se(twfe_model)["mw_treat"],
                        aggte(cs_nocov, type = "simple")$overall.se,
                        att_simple$overall.se,
                        att_never_simple$overall.se,
                        att_ra_simple$overall.se,
                        att_ipw_simple$overall.se)
) |>
  mutate(
    t_stat  = round(ATT / SE, 3),
    ATT     = round(ATT, 4),
    SE      = round(SE, 4),
    sig     = case_when(abs(t_stat) > 2.576 ~ "***",
                        abs(t_stat) > 1.960 ~ "**",
                        abs(t_stat) > 1.645 ~ "*",
                        TRUE ~ "")
  )

print(results_summary)
write_csv(results_summary, "outputs/R_results_summary.csv")

cat("\n============================================================\n")
cat("  CALLAWAY & SANT'ANNA — KEY RESULTS SUMMARY\n")
cat("============================================================\n")
cat(sprintf("  Overall ATT (DR, not-yet-treated): %.4f (SE = %.4f)\n",
            att_simple$overall.att, att_simple$overall.se))
cat("  Figures saved to: figures/\n")
cat("  Results saved to: outputs/R_results_summary.csv\n")
cat("============================================================\n")
