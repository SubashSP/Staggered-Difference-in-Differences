# Staggered Difference-in-Differences
## Theory, Assumptions & Implementation
### Focusing on the Callaway & Sant'Anna (2021) Estimator

---

## Table of Contents

1. [Introduction & Motivation](#1-introduction--motivation)
2. [Review: Classical DiD](#2-review-classical-difference-in-differences)
3. [The Staggered Treatment Problem](#3-the-staggered-treatment-problem)
4. [Key Concepts in Staggered DiD](#4-key-concepts-in-staggered-did)
5. [Core Assumptions](#5-core-assumptions)
6. [Callaway & Sant'Anna (2021) Estimator](#6-callaway--santanna-2021-estimator)
7. [Installing Required Stata Packages](#7-installing-required-stata-packages)
8. [Data Setup and Requirements](#8-data-setup-and-requirements)
9. [Implementing csdid in Stata: Full Walkthrough](#9-implementing-csdid-in-stata-full-walkthrough)
10. [Interpretation of Results](#10-interpretation-of-results)
11. [Testing Assumptions & Pre-trends](#11-testing-assumptions--pre-trends)
12. [Aggregation Strategies](#12-aggregation-strategies)
13. [Event Study Plots](#13-event-study-plots)
14. [Comparison with Other Methods](#14-comparison-with-other-methods)
15. [Common Pitfalls & Troubleshooting](#15-common-pitfalls--troubleshooting)
16. [Full Empirical Example: Minimum Wage & Employment](#16-full-empirical-example-minimum-wage--employment)
17. [Summary & Further Reading](#17-summary--further-reading)

---

## 1. Introduction & Motivation

Difference-in-Differences (DiD) is one of the most widely used quasi-experimental methods in applied economics and social sciences. In its classical form, it compares the change in outcomes over time between a treatment group and a control group. The identifying assumption is that, absent treatment, both groups would have followed **parallel trends**.

However, modern policy and program evaluations rarely involve a single, simultaneous rollout. Instead, treatment is often adopted at different times by different units — a setting known as **staggered treatment adoption** or staggered DiD.

> ⚠️ **The Problem with Standard TWFE in Staggered Settings**
>
> Recent econometric research (Goodman-Bacon 2021; Callaway & Sant'Anna 2021; Sun & Abraham 2021; de Chaisemartin & D'Haultfoeuille 2020) has shown that the standard Two-Way Fixed Effects (TWFE) estimator — the workhorse DiD regression — can produce **severely biased estimates** when treatment effects are heterogeneous across time and groups.
>
> The TWFE estimator uses already-treated units as controls when estimating treatment effects for later-treated units. This leads to **negative weights** on some comparisons, and can even reverse the sign of the true treatment effect.

These problems were largely unknown or overlooked for decades. The solution requires using estimators explicitly designed for staggered settings — of which the Callaway & Sant'Anna (2021) estimator, implemented in Stata as `csdid`, is the most comprehensive and widely adopted.

---

## 2. Review: Classical Difference-in-Differences

### 2.1 The 2×2 Case

The canonical DiD setup involves two groups (treated and control) and two time periods (before and after). The DiD estimator is:

$$\delta_{\text{DiD}} = (\bar{Y}^1_1 - \bar{Y}^1_0) - (\bar{Y}^0_1 - \bar{Y}^0_0)$$

Where $\bar{Y}^g_t$ = mean outcome for group $g$ at time $t$; superscript 1 = treated, 0 = control; subscript 1 = post, 0 = pre.

### 2.2 TWFE Regression

In practice, DiD is estimated via a Two-Way Fixed Effects regression:

```stata
* Standard TWFE DiD regression
reghdfe Y D_it, absorb(i t) vce(cluster i)

* Where:
* Y    = outcome variable
* D_it = binary treatment indicator (1 if unit i treated at time t)
* i    = unit fixed effects
* t    = time fixed effects
```

This absorbs unit-specific time-invariant heterogeneity (unit FE) and common time shocks (time FE). The coefficient on `D_it` is the DiD estimate of the Average Treatment Effect on the Treated (ATT).

---

## 3. The Staggered Treatment Problem

### 3.1 What is Staggered Treatment?

Staggered treatment occurs when different units adopt treatment at different calendar times. A common example: U.S. states adopting a policy (minimum wage increases, marijuana legalization, Medicaid expansion) in different years.

| **Unit/State** | **2015** | **2016** | **2017** | **2018** | **Cohort** |
|----------------|----------|----------|----------|----------|------------|
| State A        | 0        | 0        | **1**    | **1**    | 2017       |
| State B        | 0        | **1**    | **1**    | **1**    | 2016       |
| State C        | 0        | 0        | 0        | **1**    | 2018       |
| State D        | 0        | 0        | 0        | 0        | Never      |

In this setting, TWFE uses State B (already treated in 2016) as a control when estimating the treatment effect for State A (first treated in 2017). If treatment effects grow over time, this contamination leads to bias.

### 3.2 The Goodman-Bacon Decomposition

Goodman-Bacon (2021) showed that the TWFE estimator is a **weighted average of all possible 2×2 DiD comparisons**. Critically, some comparisons use already-treated units as controls — and when treatment effects are heterogeneous, the weights can be negative.

```stata
* Install and run Bacon Decomposition
ssc install bacondecomp

bacondecomp Y D_it, ddetail
* Decomposes TWFE into 2×2 DiD comparisons and shows the weight on each
```

---

## 4. Key Concepts in Staggered DiD

### 4.1 Cohorts (Treatment Groups)

A **cohort** $g$ is defined as the group of units that first receive treatment in calendar time period $g$. Units that are never treated form the "never-treated" group.

**Notation**: $G_i$ = first period unit $i$ is treated ($G_i = \infty$ for never-treated units)

### 4.2 Group-Time Average Treatment Effects — ATT(g,t)

The key object of interest is the **Group-Time ATT**:

$$\text{ATT}(g, t) = \mathbb{E}[Y_t(g) - Y_t(0) \mid G = g]$$

The average treatment effect at calendar time $t$, for units whose treatment first began at time $g$ (cohort $g$), compared to their counterfactual potential outcome $Y_t(0)$ under no treatment.

### 4.3 Event Time vs. Calendar Time

| Concept | Definition |
|---------|-----------|
| **Calendar time** $t$ | Actual year/period |
| **Event time** $\ell = t - g$ | Time relative to when the unit was first treated |
| $\ell = 0$ | Period treatment begins |
| $\ell = -1$ | Pre-treatment period just before adoption |

### 4.4 Clean Controls

A key innovation in Callaway & Sant'Anna is using only **"clean" control units** — units that have not yet been treated or are never treated — when forming comparisons. This avoids the contamination problem of TWFE.

---

## 5. Core Assumptions

### 5.1 Summary of Assumptions

| Assumption | Description |
|-----------|-------------|
| **Staggered Adoption (Irreversibility)** | Once treated, units remain treated. No de-adoption. If unit $i$ is treated at $g$, then $D_{it} = 1$ for all $t \geq g$. |
| **Parallel Trends (Conditional)** | Absent treatment, the average outcome for cohort $g$ would have followed the same trend as the comparison group, possibly after conditioning on covariates $X$. |
| **No Anticipation** | Units do not change behavior before treatment begins. Formally: $\text{ATT}(g,t) = 0$ for all $t < g$. |
| **Overlap (Positivity)** | For each cohort $g$ and time $t$, there is a positive probability of being in the comparison group. |
| **Sampling Independence** | Units are sampled independently. Cluster SEs at the unit level to address within-unit serial correlation. |

### 5.2 Parallel Trends in Detail

| **Unconditional PT** | **Conditional PT** |
|----------------------|-------------------|
| $\mathbb{E}[Y_t(0) \mid G=g] - \mathbb{E}[Y_{g-1}(0) \mid G=g] = \mathbb{E}[Y_t(0) \mid G=\infty] - \mathbb{E}[Y_{g-1}(0) \mid G=\infty]$ | Parallel trends holds after conditioning on pre-treatment covariates $X_i$ |
| Does not require conditioning on covariates. **Stronger assumption.** | Weaker, more credible assumption. |
| `csdid Y, ivar(id) tvar(t) gvar(g)` | `csdid Y X1 X2, ivar(id) tvar(t) gvar(g)` |

### 5.3 No Anticipation

Units do not change behavior in anticipation of future treatment. In practice, some anticipation may be acceptable if researchers re-define the treatment date. Significant pre-trends at $\ell = -2$ or earlier may indicate violations of no-anticipation or parallel trends.

---

## 6. Callaway & Sant'Anna (2021) Estimator

### 6.1 Overview

Callaway & Sant'Anna (2021), published in the *Journal of Econometrics*, propose a two-step approach:
1. **Estimate** group-time ATTs for every $(g, t)$ pair
2. **Aggregate** these ATTs into meaningful summary statistics

### 6.2 Step 1: Estimating ATT(g,t)

For each cohort $g$ and time period $t$, compare cohort $g$ to a comparison group. Four estimation methods are available:

| Method | Stata Option | Description |
|--------|-------------|-------------|
| **Regression Adjustment (RA)** | `method(reg)` | OLS regression on outcome differences. Parametric. |
| **Inverse Probability Weighting (IPW)** | `method(ipw)` | Propensity score reweighting. Semiparametric. |
| **Doubly Robust (DR)** | `method(dripw)` | Combines RA + IPW. Consistent if **either** model is correct. **DEFAULT & RECOMMENDED.** |
| **Doubly Robust (DRDID)** | `method(drimp)` | Alternative DR using improved propensity score estimation. |

### 6.3 Step 2: Aggregation

The ATT(g,t) estimates are aggregated into four summary parameters:

- **Simple Average ATT** — single weighted average across all (g,t) pairs
- **Group/Cohort ATT** — average ATT for each cohort, averaged over post-treatment periods
- **Calendar-Time ATT** — ATT at each calendar time, averaged over treated cohorts
- **Dynamic/Event-Study ATT** — ATT at each event time $\ell = t - g$ — **most useful for visualizing treatment dynamics**

---

## 7. Installing Required Stata Packages

```stata
* 1. Main package: csdid (Callaway-Sant'Anna DiD)
ssc install csdid

* 2. DRDID (doubly-robust DiD building block)
ssc install drdid

* 3. For event-study plots
ssc install coefplot

* 4. For Bacon decomposition diagnostic
ssc install bacondecomp

* 5. For parallel trends testing
ssc install pretrends

* 6. For comparison with Sun-Abraham estimator
ssc install eventstudyinteract

* 7. For TWFE benchmark
ssc install reghdfe
ssc install ftools

* Verify installation
which csdid
which drdid
```

> 💡 **Version Note**: `csdid` requires **Stata 16 or higher** due to its use of Stata frames. For older versions, use `csdid_old` or the `didimputation` package.

---

## 8. Data Setup and Requirements

### 8.1 Required Data Structure

`csdid` requires a balanced or unbalanced panel in **long format**:

| Variable | Requirement | Description |
|----------|------------|-------------|
| **Unit ID** (`ivar`) | Numeric integer | Unique unit identifier |
| **Time** (`tvar`) | Numeric, consecutive integers | Calendar time period |
| **Cohort** (`gvar`) | Numeric or 0/missing | First period treated; **0 or missing = never-treated** |
| **Outcome** Y | Numeric | The outcome variable |
| **Covariates** X | Numeric | Pre-treatment or time-invariant covariates |

### 8.2 Creating the Cohort Variable

```stata
* Method 1: From an existing treat_year variable
gen gvar = treat_year
replace gvar = 0 if treat_year == .       // never-treated = 0

* Method 2: Derive from a 0/1 treatment indicator
bysort id (year): gen first_treat = year if treated == 1
bysort id: egen gvar = min(first_treat)
replace gvar = 0 if gvar == .

* Verify
tab gvar                                  // check distribution
list id year treated gvar in 1/20         // spot check
```

### 8.3 Panel Structure Checks

```stata
xtset id year                             // declare panel
xtdescribe                                // check for gaps
duplicates report id year                 // check for duplicates
tab gvar, missing                         // check cohort variable

* Verify irreversibility
bysort id (year): gen reversal = (treated < L.treated)
tab reversal                              // should be all zeros
drop reversal
```

---

## 9. Implementing csdid in Stata: Full Walkthrough

### 9.1 Basic Syntax

```
csdid depvar [indepvars], ivar(varname) tvar(varname) gvar(varname) [options]
```

- `depvar` — outcome variable
- `indepvars` — (optional) covariates for conditional parallel trends
- `ivar()` — unit identifier
- `tvar()` — time variable
- `gvar()` — cohort variable (first treatment period; 0 = never treated)

### 9.2 Unconditional Parallel Trends

```stata
csdid Y, ivar(id) tvar(year) gvar(gvar)
estimates store csdid_base
```

### 9.3 Conditional Parallel Trends (With Covariates)

```stata
* Doubly-robust (DEFAULT, RECOMMENDED)
csdid Y covX1 covX2, ivar(id) tvar(year) gvar(gvar) method(dripw)

* Regression adjustment only
csdid Y covX1 covX2, ivar(id) tvar(year) gvar(gvar) method(reg)

* IPW only
csdid Y covX1 covX2, ivar(id) tvar(year) gvar(gvar) method(ipw)
```

### 9.4 Control Group Options

```stata
* Default: never-treated + not-yet-treated as controls
csdid Y, ivar(id) tvar(year) gvar(gvar) notyet

* Only never-treated (more conservative)
csdid Y, ivar(id) tvar(year) gvar(gvar) notyetnever
```

> When there are few never-treated units, `notyet` is preferred (more data, weaker assumption).

### 9.5 Inference Options

```stata
* Clustered SEs at unit level (RECOMMENDED for panel data)
csdid Y, ivar(id) tvar(year) gvar(gvar) cluster(id)

* Wild bootstrap (better for small samples, N < 40 clusters)
csdid Y, ivar(id) tvar(year) gvar(gvar) wboot reps(999) seed(12345)

* Multiplier bootstrap (asymptotic, faster for large datasets)
csdid Y, ivar(id) tvar(year) gvar(gvar) reps(999)
```

### 9.6 Full Recommended Specification

```stata
csdid Y covX1 covX2,            ///
    ivar(id)        ///          unit identifier
    tvar(year)      ///          time variable
    gvar(gvar)      ///          cohort (first treatment year; 0 = never treated)
    method(dripw)   ///          doubly-robust (RECOMMENDED)
    notyet          ///          use not-yet-treated as additional controls
    cluster(id)     ///          cluster SEs at unit level
    asinr            ///          store results for csdid_plot

estimates store csdid_main

* Aggregations
estat simple      // overall ATT
estat group       // ATT by cohort
estat calendar    // ATT by calendar time
estat event       // ATT by event time (MOST COMMON)
```

---

## 10. Interpretation of Results

### 10.1 Reading the ATT(g,t) Table

```stata
* Example ATT(g,t) output (stylized):
*
* ATT(g,t)       | Coeff   Std.Err.   t    P>|t|   [95% Conf. Interval]
* ───────────────+──────────────────────────────────────────────────────
* g=2016, t=2015 | -0.021   0.018   -1.17  0.243   -0.056   0.014
* g=2016, t=2016 |  0.087   0.025    3.48  0.001    0.038   0.136
* g=2016, t=2017 |  0.112   0.031    3.61  0.000    0.051   0.173
* g=2017, t=2015 |  0.011   0.022    0.50  0.617   -0.032   0.054
* g=2017, t=2016 | -0.018   0.019   -0.95  0.343   -0.055   0.019
* g=2017, t=2017 |  0.064   0.028    2.29  0.022    0.009   0.119
*
* g=2016, t=2015: PRE-TREATMENT for 2016 cohort → should be ~0
* g=2016, t=2016: TREATMENT ONSET → should be significant if there is an effect
* g=2016, t=2017: POST-TREATMENT → magnitude and direction of the effect
```

### 10.2 Summary Aggregations

```stata
estimates restore csdid_main

* 1. Overall average ATT
estat simple
* → Single weighted average ATT across all groups and periods

* 2. Cohort-specific ATTs
estat group
* → ATT for each cohort; useful to detect heterogeneous responses

* 3. Calendar-time ATTs
estat calendar
* → ATT at each calendar year, averaged over treated cohorts

* 4. Dynamic/Event-study ATTs (MOST RECOMMENDED)
estat event
* → l < 0: PRE-TRENDS (should be ~0 if parallel trends holds)
* → l = 0: effect at treatment onset
* → l > 0: dynamic post-treatment effects
```

---

## 11. Testing Assumptions & Pre-trends

### 11.1 Pre-trends Test

```stata
csdid Y covX1 covX2, ivar(id) tvar(year) gvar(gvar) method(dripw)
estat event, post

* Test individual pre-period coefficients
test _b[Tm2] = 0
test _b[Tm3] = 0

* Joint test of all pre-treatment periods
test _b[Tm3] = 0 & _b[Tm2] = 0
* p-value > 0.05 = no evidence against parallel trends

* Using pretrends package
pretrends, numpre(3)
```

### 11.2 Checking Treatment Irreversibility

```stata
sort id year
by id: gen ever_treated   = max(treated)
by id: gen once_untreated = (treated < L.treated)
tab once_untreated        // should be all zeros
```

### 11.3 Bacon Decomposition

```stata
bacondecomp Y D_it, ddetail
* If large share of weight is from 'Already Treated vs Timing'
* → strong motivation to use csdid instead of TWFE

* Plot the decomposition
bacondecomp Y D_it, ddetail grcolor(red%40 blue%40 green%40)
```

---

## 12. Aggregation Strategies

### 12.1 csdid_plot

```stata
csdid Y, ivar(id) tvar(year) gvar(gvar) asinr

csdid_plot,                                          ///
    title("Event Study: Effect of Policy on Y")      ///
    xtitle("Years Relative to Treatment Onset")      ///
    ytitle("ATT Estimate")                           ///
    xline(0, lpattern(dash) lcolor(red))             ///
    yline(0, lpattern(dot))

graph export "figures/event_study.png", replace width(1200)
```

### 12.2 Publication-Quality Event Study with coefplot

```stata
csdid Y covX1, ivar(id) tvar(year) gvar(gvar) method(dripw)
estat event, post

coefplot,                                            ///
    vertical                                         ///
    yline(0, lwidth(thin) lpattern(dot))             ///
    xline(3.5, lwidth(medium) lpattern(dash) lcolor(red)) ///
    title("Callaway-Sant'Anna Event Study")          ///
    xtitle("Event Time (t − g)")                     ///
    ytitle("ATT(g,t) Estimate")                      ///
    msymbol(circle) mcolor(navy)                     ///
    ciopts(lcolor(navy%60))
```

### 12.3 Collecting All Aggregations

```stata
* Run and store each aggregation
csdid Y, ivar(id) tvar(year) gvar(gvar)

estat simple, post;    estimates store att_simple
estat group, post;     estimates store att_group
estat calendar, post;  estimates store att_calendar
estat event, post;     estimates store att_event

* Present with esttab
esttab att_simple att_group,                         ///
    cells(b(fmt(3)) se(par fmt(3)))                  ///
    title("Callaway-Sant'Anna ATT Estimates")
```

---

## 13. Event Study Plots

### 13.1 Standard Event Study

```stata
csdid Y covX1 covX2, ivar(id) tvar(year) gvar(gvar) method(dripw) cluster(id)
estat event

* Output labels:
* Pre_avg / Tm3, Tm2, Tm1 = ATTs at l=-3,-2,-1 (pre-trends; should be ~0)
* Post_avg / T0, T1, T2   = ATTs at l=0,1,2,... (post-treatment dynamics)
```

### 13.2 Comparing TWFE vs. CS Event Study

```stata
* TWFE event study
reghdfe Y ib(-1).event_time_var, absorb(id year) cluster(id)
estimates store twfe_es

* CS event study
csdid Y, ivar(id) tvar(year) gvar(gvar)
estat event, post
estimates store cs_es

* Plot both on same graph
coefplot (twfe_es, label(TWFE)   mcolor(red)  ciopts(lcolor(red%40)))  ///
         (cs_es,   label(CS DiD) mcolor(blue) ciopts(lcolor(blue%40))), ///
    vertical yline(0) xline(3.5, lpattern(dash))                         ///
    title("TWFE vs. Callaway-Sant'Anna Event Study")
```

---

## 14. Comparison with Other Methods

| Estimator | Stata Command | Handles Het. TE? | Clean Controls? | Best For |
|-----------|--------------|:-----------------:|:---------------:|----------|
| TWFE DiD | `reghdfe` | ✗ No | ✗ No | Simple settings, homogeneous TE |
| Goodman-Bacon | `bacondecomp` | ✔ Diagnostic | Partial | Diagnosing TWFE bias |
| **Callaway & Sant'Anna** | **`csdid`** | **✔ Yes** | **✔ Yes** | **General staggered DiD (RECOMMENDED)** |
| Sun & Abraham | `eventstudyinteract` | ✔ Yes | ✔ Yes | Event-study focused |
| de Chaisemartin & D'H | `did_multiplegt` | ✔ Yes | ✔ Yes | Allows treatment reversals |
| Borusyak et al. | `didimputation` | ✔ Yes | ✔ Yes | Imputation-based approach |

---

## 15. Common Pitfalls & Troubleshooting

> ⚠️ **Common Mistakes & How to Fix Them**

- **`gvar` not set to 0 for never-treated:** Must be 0 or missing, not a large number like 9999.
  ```stata
  replace gvar = 0 if never_treated == 1
  ```

- **Non-consecutive time variable:** `csdid` requires `tvar` to be consecutive integers.
  ```stata
  egen year_seq = group(year), label   // recode to 1, 2, 3, ...
  ```

- **Time-varying covariates:** Include only pre-treatment or time-invariant covariates. Post-treatment covariates cause post-treatment bias.

- **Small sample with bootstrap:** Wild bootstrap (`wboot`) is recommended when N < 40 clusters.

- **No never-treated units:** Use `notyet` option to include not-yet-treated as controls. Without any controls, CS cannot run.

- **Wrong Stata version:** `csdid` requires Stata 16+. Use `csdid_old` or `didimputation` for older versions.

---

## 16. Full Empirical Example: Minimum Wage & Employment

### 16.1 Data Loading & Preparation

```stata
use "minwage_panel.dta", clear

* Variables:
*   countyid   - FIPS county code (unit)
*   year       - calendar year (2000–2019)
*   lemp       - log(teen employment)  [OUTCOME]
*   mw_treat   - =1 if state raised MW in or before this year
*   treat_year - year state first raised MW above federal level
*   lpop       - log(county population)  [COVARIATE]
*   lincome    - log(median household income)  [COVARIATE]

xtset countyid year
xtdescribe
```

### 16.2 Construct Cohort Variable

```stata
gen gvar = treat_year
replace gvar = 0 if treat_year == .       // never-treated = 0

* Verify irreversibility
sort countyid year
by countyid: gen reversal = (mw_treat < L.mw_treat)
assert reversal == 0
drop reversal

tab gvar, missing
```

### 16.3 TWFE Benchmark

```stata
reghdfe lemp mw_treat, absorb(countyid year) cluster(countyid)
estimates store twfe
```

### 16.4 Callaway-Sant'Anna

```stata
* No covariates
csdid lemp, ivar(countyid) tvar(year) gvar(gvar)                ///
    method(dripw) notyet cluster(countyid)
estimates store cs_nocov

* With covariates
csdid lemp lpop lincome, ivar(countyid) tvar(year) gvar(gvar)   ///
    method(dripw) notyet cluster(countyid)
estimates store cs_cov

* Event study plot
csdid lemp lpop lincome, ivar(countyid) tvar(year) gvar(gvar)   ///
    method(dripw) notyet cluster(countyid) asinr

csdid_plot, title("CS Event Study: Minimum Wage on Log Employment")  ///
    xtitle("Years Since MW Increase") ytitle("ATT (log employment)")

graph export "cs_event_study.png", replace
```

### 16.5 Extract and Interpret Results

```stata
estimates restore cs_cov

estat simple    // overall ATT: e.g. -0.032 → ~3.2% decrease in teen employment
estat group     // cohort-specific ATTs
estat event     // dynamic ATTs: check Tm2, Tm1 ≈ 0 (no pre-trends)
```

### 16.6 Sensitivity Checks

```stata
* Never-treated controls only
csdid lemp lpop lincome, ivar(countyid) tvar(year) gvar(gvar) ///
    method(dripw) notyetnever cluster(countyid)
estat simple

* Regression adjustment
csdid lemp lpop lincome, ivar(countyid) tvar(year) gvar(gvar) ///
    method(reg) notyet cluster(countyid)
estat simple

* Wild bootstrap
csdid lemp lpop lincome, ivar(countyid) tvar(year) gvar(gvar) ///
    method(dripw) notyet wboot reps(999) seed(42)
estat simple
```

---

## 17. Summary & Further Reading

### 17.1 Key Takeaways

1. **Standard TWFE DiD is biased** when treatment effects are heterogeneous across groups and time in staggered adoption settings.

2. The **Callaway & Sant'Anna (2021)** estimator avoids this bias by estimating ATT(g,t) parameters using only clean comparison units.

3. **Key assumptions**: staggered adoption (irreversibility), parallel trends (conditional or unconditional), no anticipation, and overlap.

4. **Always validate assumptions**: test for pre-trends, check treatment irreversibility, and run Bacon decomposition to quantify TWFE contamination.

5. Use **doubly-robust estimation** (`method(dripw)`) as the default — it is consistent if either the outcome or propensity score model is correct.

6. Report **dynamic/event-study ATTs** via `estat event` to visualize both pre-trends validity and post-treatment dynamics.

7. Consider **sensitivity analysis**: vary the control group (`notyet` vs. `notyetnever`), estimation method, and inference approach.

### 17.2 Recommended Further Reading

| Paper | Key Contribution |
|-------|-----------------|
| **Callaway & Sant'Anna (2021)**, *J. of Econometrics* | CS estimator, group-time ATTs, aggregation strategies |
| **Goodman-Bacon (2021)**, *J. of Econometrics* | Decomposition of TWFE, negative weights problem |
| **Sun & Abraham (2021)**, *J. of Econometrics* | Interaction-weighted estimator for event studies |
| **de Chaisemartin & D'Haultfoeuille (2020)**, *AER* | DiD with treatment reversals, local ATT |
| **Roth et al. (2023)**, *J. of Econometrics* | Comprehensive survey of modern DiD methods |
| **Baker et al. (2022)**, *J. of Finance* | Application guide for staggered DiD in finance |
| **Borusyak, Jaravel & Spiess (2024)** | Imputation-based DiD estimator (`didimputation`) |

### 17.3 Quick Reference: csdid Cheat Sheet

```stata
* ================================================================
*  CSDID QUICK REFERENCE CHEAT SHEET
* ================================================================

* INSTALLATION
ssc install csdid; ssc install drdid; ssc install coefplot

* DATA: long format panel — id, t, gvar (0 = never-treated), Y
*       tvar must be consecutive integers

* BASIC
csdid Y, ivar(id) tvar(t) gvar(gvar)

* WITH COVARIATES (conditional PT)
csdid Y X1 X2, ivar(id) tvar(t) gvar(gvar) method(dripw)

* CONTROL GROUP:
*   notyet      = never-treated + not-yet-treated (DEFAULT)
*   notyetnever = only never-treated (more conservative)

* INFERENCE:
*   cluster(id) = clustered SEs (RECOMMENDED)
*   wboot       = wild bootstrap (small samples)

* AGGREGATIONS (run after csdid)
estat simple     // overall ATT
estat group      // ATT by cohort
estat calendar   // ATT by calendar time
estat event      // ATT by event time → EVENT STUDY

* PLOTTING
csdid Y, ivar(id) tvar(t) gvar(gvar) asinr
csdid_plot       // automatic event study plot

* ================================================================
```

---

*End of Lecture Notes — Staggered DiD & Callaway-Sant'Anna Estimator*
